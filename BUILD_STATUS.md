# De novo optimization redesign — build status

Replaces the retired map-back sweep (archived in `_retired_sweep/`) with a
weight-free rank-aggregation selector over multiple non-circular signals.

## Design (agreed)
- Three branches off `repair.out`: production (intact, untouched), assembly
  (intact subset -> candidate references), concordance (6 high-depth individuals
  split 50/50, used ONLY in stage 2). Assembly subset and concordance individuals
  may OVERLAP (default).
- Signals: 1a contig/length inflection, 1b internal redundancy, 2 biological
  locus-count proximity (optional, input-gated), 5b NB-mixture for cutoff1,
  4 bcftools SNP recovery (expensive), concordance (expensive).
- Two-stage: cheap signals -> provisional rank -> gap-based top-N (clamped
  [3,10]) -> bcftools on survivors -> weight-free rank aggregation -> cutoffs.
- Soft constraints (locus anchor, redundancy) enter as rankings, never gates.
- `--cutoff1/--cutoff2/--cluster_similarity` still pin/override.

## DESIGN INVARIANT (Position B — chosen): full-sample candidates
- Candidate references are built from ALL samples (NOT a subset). Each candidate
  is a real, production-grade assembly. The WINNER *is* the production reference.
- finalize (chunk 6) does NOT re-assemble — it SELECTS the winning candidate's
  existing denovo_reference.fa by meta.id and PUBLISHES it to a clean stable path.
  finalize becomes the single publishDir for the production reference; candidate
  publishDirs stay FOR NOW (dev inspection) and get removed once finalize owns it.
- Rationale: cutoff2 (= n individuals) means the same thing during optimization
  and in the final reference only if both use the same sample count. A subset
  would optimize cutoff2 on the wrong N. Position B removes that mismatch and the
  chunk-6 re-assembly step entirely.
- Compute: every candidate runs on all samples, but the k2/s0.80 pathological
  corner is already pruned (cutoff2_floor=3, no s0.80), so the surviving grid is
  the fast region.

## map_reads cardinality bug — only 1 sample mapped (FIXED)
Symptom: downstream map_reads/samtools_stats ran 1 of 1 (not 45). This is the
old "n=1 samples" report bug, now root-caused.
Cause: map_reads(repair.out, genome_indexed) has TWO queue-channel inputs. The
genome (genome_indexed = prepare_genome_local.out.genome) emits ONCE. A process
with two queue inputs consumes them in lockstep and STOPS when the shorter
(1-emission genome) is exhausted -> only the first read sample maps.
FIX: genome must be a VALUE channel so it broadcasts to all 45 read emissions:
  map_reads( repair.out, genome_indexed.first() )
Applied to BOTH map_reads calls (denovo branch line ~353 AND reference-genome
branch line ~295 — same latent bug in both).
NOTE: this also fixes the separate "Mapped Reads n=1 / 1.3% retention" report
issue we shelved earlier — same root cause.

## Full-sample candidate contig counts (Position B, real data 10-Jun)
Grid behaves sensibly. cutoff1 (c) dominates; cutoff2 (k3 vs k4) barely matters
at 45 individuals; similarity nudges up monotonically:
  c2 ~88-93k, c3 ~71-75k, c4 ~59-62k, c5 ~48-51k contigs.
Still ~2x size spread across grid -> size-bias concern remains live -> validates
using non-circular signals (NB-mixture/redundancy/anchor/concordance), NOT map rate.

## CHUNK 3a/Position B — finalize filename-collision fix
Position B made every candidate's reference literally `denovo_reference.fa`.
The old finalize collected ALL candidate fastas into one process -> Nextflow
"input file name collision: multiple input files denovo_reference.fa".
FIX (also the correct chunk-6 selection logic): finalize now SELECTS only the
winning candidate by joining best_id against the candidates channel:
  winning_candidate = candidates.map{meta,fa->[id,meta,fa]}
      .combine(best_id_ch).filter{id==best}.map{->[meta,fa]}
  -> single (meta,fa) passed to finalize (no collect, no collision).
finalize now:
  - takes tuple(meta, path(winner, stageAs:'winner_input.fa'))  [stageAs avoids
    cp-onto-itself since the real file is already named denovo_reference.fa]
  - publishDir ${outdir}/denovo_assembly  [single publish point for production ref]
  - tag finalize:${meta.id}
Run before this fix: grid = 8 filter -> 24 candidates (8 cutoff x 3 sims), all
real full-sample assemblies built OK; chain ran through aggregate; only finalize
errored on the collision. So Position B assembly + scoring topology all work.

## PARAM rename (Position B)
- optimize_sample_pct -> snp_sample_pct: now the fraction of samples used for the
  STAGE-2 SNP-recovery signal (chunk 5), NOT an assembly subset.
- optimize_seed: retained — now governs the pseudo-rep split + stage-2 SNP-sample
  selection (no longer an assembly subset seed).

## (superseded) earlier invariant: subset vs final reference
- CANDIDATE references (the ~18 grid points) are built from the SUBSET
  (optimize_sample_pct) — cheap, throwaway, used ONLY for scoring.
- The FINAL/PRODUCTION reference is re-assembled at the WINNING (c1,c2,sim) on
  ALL samples at full depth. This is the job of finalize (chunk 6, currently a
  STUB). Until chunk 6, downstream mapping uses the stubbed reference — ignore.
- Subsets/pseudo-reps NEVER become the production genome. Confirmed matches the
  original plan (subset to optimize, all samples for the real reference).

## Build order & status
- [x] CHUNK 0  Strip/archive old sweep files (_retired_sweep/)
- [x] CHUNK 1  Params + config labels  -> CHUNK1_params.txt  (APPLIED by user)
- [x] CHUNK 2  Three-branch split + pseudo-rep subworkflow skeleton (stubs)
        - [x] split_pseudo_rep.nf  (tested: mate-pairing, conservation, deterministic)
        - [x] optimize_denovo.nf  orchestrator skeleton: 3-branch topology + STUB processes
              (stub_build_candidate, stub_cheap_signals, stub_provisional_rank,
               stub_stage2, stub_aggregate, stub_finalize_reference)
        - [x] denovo_assembly.nf rewired to branch on params.do_optimize
- [~] CHUNK 3  Split into 3a (subset+grid+real assembly) and 3b (cheap signals)
        - [x] CHUNK 3a: seeded subset, real floor..knee grid, per-(c1,c2) filter
              crossed with sims -> real candidate assemblies. Cheap signals STILL
              STUBBED. New modules: filter_unique_seqs_candidate.nf,
              assemble_rainbow_candidate.nf (full N50/stats block + PE-spacer count).
        - [x] CHUNK 3b: real cheap signals + provisional rank (weight-free aggregation)
              NEW modules: compute_cheap_signals.nf (per-candidate: contig stats=1a
                inputs, self-cluster redundancy=1b, n_contigs for anchor=2),
                fit_nb_mixture.nf (global NB-mixture cutoff1=5b),
                provisional_rank.nf (runs rank R script).
              NEW r_scripts: fit_nb_mixture.R (2-comp NB EM on coverage_freq ->
                posterior-crossover cutoff1; falls back to knee), provisional_rank.R
                (per-signal RANKS -> mean aggregate, weight-free; anchor dropped if
                no enzyme params; NB enters as |c1 - nb_cutoff1| proximity rank).
              EDITED: assembly_diagnostics.nf now emits coverage_freq.txt (for NB).
              optimize_denovo.nf: stubs removed, real stage-1 wired, expected-loci
                anchor computed from enzyme/genome/size-select params (cut-site model).
              main.nf: + params.optimize_redundancy_identity = 0.98
              Validated (python port, no R here): NB EM recovers 2 components +
                sensible crossover; rank aggregation on REAL contig counts picks an
                interior candidate (not size-extreme). 5b participates as ranking
                signal (option a). NB still STUBBED downstream: stage2/aggregate/
                finalize unchanged (chunks 4-6).
- [ ] CHUNK 4  Gap handoff (top-N selection by provisional-rank gap)
- [ ] CHUNK 5  Stage-2 bcftools (SNP recovery + concordance) + depth-based pseudo-rep selection
- [ ] CHUNK 6  Aggregation R script (weight-free rank aggregation) + full-depth finalize assembly
- [ ] CHUNK 7  Rewire generate_report.nf + main.nf integration (optimize_summary/plot, cutoff plots)

## CHUNK 2 install set (drop into repo, then run topology test)
  modules/split_pseudo_rep.nf
  workflows/optimize_denovo.nf
  workflows/denovo_assembly.nf   (overwrite existing)
  main.nf                        (overwrite existing — emit-fix applied, see below)

## main.nf emit-fix (applied at chunk-2 topology test)
  Live main.nf still had sweep-era report wiring out of sync with the project
  snapshot. Fixed minimally to clear "No such property: sweep_plot":
    - line 343: denovo_assembly.out.sweep_plot    -> .optimize_plot
    - line 344: denovo_assembly.out.sweep_summary -> .optimize_summary
    - params.do_sweep -> params.do_optimize (everywhere)
  NOTE: internal var names (sweep_plot_ch, final_sweep_plot, placeholder slots
  11/12, "No sweep performed" strings) were LEFT as-is — cosmetic, repoint in
  chunk 7 when the report section is properly rewired. They feed positional args
  to generate_report() and work correctly despite the stale names.

## >>> CHECKPOINT PASSED — chunk 2 verified, project resynced <<<
  Topology test (09-Jun-2026): ALL processes ran in order, pipeline completed
  (157 succeeded / 592 cached, no errors). Verified:
    - 3 branches wired: extract_uniq(45) + diagnostics; stub candidates(2);
      split_pseudo_rep(6); stub_stage2(2) -> aggregate -> finalize -> map_reads(45).
    - two-stage handoff routed (cheap -> survivors -> stage2).
    - finalize-as-reference fed production mapping (val/path fix held, no stall).
    - ISOLATION CONFIRMED: cleaned_reads/ has NO _a/_b pseudo-rep halves.
  Project memory resynced; Claude now reads true on-disk state.
  RESUME: build CHUNK 3.

## (historical) CHECKPOINT (user resyncing project memory)
  After installing the chunk-2 files + fixed main.nf, user will re-upload the
  full repo to project memory so the snapshot matches disk (stops Claude reading
  stale files). RESUME at: run the chunk-2 topology test, verify pseudo-rep
  isolation, then build CHUNK 3.

## CHUNK 3a — TWO problems found in the long 3a run (diagnosed from logs)

### Problem 2: cache miss on resume (FIXED — non-deterministic subset)
Root cause: the seeded subset shuffled cleaned_reads.toList() WITHOUT sorting
first. cleaned_reads emission order is not stable across runs, so seed+unsorted
gave a DIFFERENT 12-sample subset each run (only 3/12 overlap between runs).
Different sample set -> different collected *.uniq.seqs -> different filter hash
-> every filter_unique_seqs_candidate re-ran -> every assembly re-ran (cache bust).
FIX: sort rows by sample_id (a[0] <=> b[0]) BEFORE Collections.shuffle. Applied to
BOTH the assembly subset and the concordance .take(N) selection. Now seed+sorted
input = identical subset every run = cache-stable.
(extract_unique_seqs itself cached fine per-sample; the bust was downstream.)

### Problem 1 RESOLUTION (params changed)
- Split single cutoff_floor into independent params:
    params.cutoff1_floor = 2  (per-individual coverage; permissive is fine)
    params.cutoff2_floor = 3  (n individuals; 2 = junk + slow CD-HIT, skipped)
- params.optimize_cluster_similarity = [0.85, 0.90, 0.95]  (dropped 0.80)
- assemble_rainbow passes all 64 cores to CD-HIT already (-T) — resources OK.
- Net: removes all four multi-hour k2 jobs and the slow s0.80 word-size case.
  Grid now c1{2..knee} x c2{3..knee} x sim{0.85,0.90,0.95}.
- Files changed: main.nf (params), optimize_denovo.nf (two floors),
  denovo_assembly.nf (log line).

### Problem 1: assemblies take hours (NOT a bug — CD-HIT cost on worst corner)
Per-candidate durations (prev run): the 4 slowest were ALL *_k2_s0.80:
  c2_k2_s0.80 = 11h43m,  c3_k2_s0.80 = 6h30m,  c4_k2_s0.80 = 4h50m,
  c5_k2_s0.80 = 3h55m.  Every other candidate < ~53m, most < 10m.
Cause: CD-HIT-EST runtime ~quadratic in #sequences AND much worse at low
similarity (word size -n drops: 0.80->n5 vs 0.95->n10, weaker prefilter). So
low cutoff2 (k2 = 50-90k contigs) x low sim (0.80) explodes. Same c2_k2 input
runs 24m at s0.95 but 11h+ at s0.80.
RECOMMENDED MITIGATION (params, no code change):
  - cutoff_floor = 3   -> drops all k2 candidates (kills all 4 multi-hour jobs;
    also removes the biologically junky 2-individual-locus corner = the bias
    region we already flagged).
  - optimize_cluster_similarity = [0.85, 0.90, 0.95]  -> drop 0.80 (was never
    competitive in the old sweep AND is the slow-word-size case).
  Together: grid 32 -> ~18 candidates, no multi-hour jobs, compute concentrated
  on the sensible region. Also check denovo_assembly label cpus (CD-HIT -T).
  DECISION PENDING with user.

## CHUNK 3a — combine-spreading fix (applied after first 3a test)
First 3a run: subset worked (extract_unique_seqs 12 of 12 = 25% of 45 ✓),
diagnostics + split_pseudo_rep ran. FAILED at cutoff_combos with
"Invalid method invocation call with arguments: [2,3,4,5,2,3]" — c1_vals/c2_vals
lists were SPREAD by .combine() instead of staying as single slots (same class of
bug as the old sweep). FIX: wrap each with .map{ [it] }:
   cutoff_combos = c1_vals.map{ [it] }.combine( c2_vals.map{ [it] } ).flatMap{...}
Audited all other .combine() calls — the rest already wrap lists or combine
scalars/tuples correctly.

## CHUNK 3a — install + test
Install:
  modules/filter_unique_seqs_candidate.nf   (new)
  modules/assemble_rainbow_candidate.nf      (new)
  workflows/optimize_denovo.nf               (overwrite — real subset/grid/assembly)

Run (same command as chunk 2 topology test). With knees ~c1=5,c2=4 and floor=2,
expect the grid to expand to:
  - filter_unique_seqs_candidate: 12 tasks (one per c1,c2 in floor..knee)
  - assemble_rainbow_candidate:   48 tasks (12 x 4 similarities)
  - candidate ids like c2_k2_s0.8 .. c5_k4_s0.95 (trailing zero dropped; unique)
Each candidate publishes a REAL denovo_reference.fa under
  illumina_qc/.../denovo_assembly/optimize/candidates/<id>/
with real contig counts + full assembly_stats.txt (N50, length dist, PE-spacer loci).

Cheap signals/stage-2/finalize STILL STUBBED, so:
  - stub_finalize_reference looks for candidate_<id>.fa but real candidates are
    named denovo_reference.fa -> finalize falls back to its 1-contig stub.
    => downstream mapping still runs on a stub reference. EXPECTED in 3a; chunk 6
       rewrites finalize to re-assemble the winner at full depth.
  - Verify the 12/48 task fan-out and that candidate references have sane,
    cutoff-dependent contig counts (k2 >> k4), like the old sweep did.

Known 3a id detail: Groovy renders 0.80->"0.8", 0.90->"0.9" in ids. Internally
consistent (join matches), just cosmetic. Will tidy display in report (chunk 7).

## DECISION: go 1-PASS (SNP-call all candidates); cheap signals become ranking inputs
84-candidate run diagnosis: cheap signals keep proving unreliable to GATE on.
Redundancy died, inflection/NB collapse onto c1, and CV/Gini are SIZE-CONFOUNDED
and ANTI-CORRELATED (corr n_contigs~CV=+0.80, ~Gini=-0.91, CV~Gini=-0.56) -> they
cancel. A signal too unreliable to select on is too unreliable to gate on. So:
- ABANDON 2-pass gating. SNP recovery (the non-circular signal that motivated the
  whole redesign) runs on ALL candidates. Cheap signals become ranking INPUTS, not
  a gate. CHUNK 4 (gap-based top-N) is DELETED from the plan — no longer needed.
- stage2_min/max_candidates params become irrelevant.
- provisional_rank + stage2 + aggregate collapse into ONE rank_and_select.

GRID INSIGHT (84-run): cutoff2 was NOT flat — k3/k4 were similar but k5/k6 collapse
contig count hard (c5_k3=50289 -> c5_k5=9402). Wide grid was the right call; winner
moved to interior high-stringency (c5_k5_s0.90, ~9400 contigs).

STACKS validation (Paris et al 2017 r80 rule): STACKS optimizes on # POLYMORPHIC
loci shared across >=80% of individuals (r80) — deliberately size-robust (junk/
paralog loci fail the 80% share test). Adopting: SNP signal = concordance (primary,
pseudo-rep genotype agreement) + r80 shared-polymorphic-loci + SNPs-per-locus
(secondary). r80 threshold = param. Compute all, prune after inspection.

SNP module structure: 2 bcftools passes per candidate — concordance from pseudo-rep
halves; r80 + SNP-density from snp_sample_pct subset.

## CHUNK 5a — coverage signal upgraded to per-base depth (DONE; build sequenced 5a then 5b)
The CV/Gini size-confound was an idxstats artifact (read-count dominated by # junk
contigs). FIX: compute_coverage_cv.nf now uses samtools depth -a -> per-contig MEAN
depth -> MEDIAN-NORMALIZED statistics. Emits THREE (prune later):
  frac_hi (PRIMARY) = fraction of contigs with mean depth > cv_pileup_mult*median
                      (direct paralog-pileup count; param cv_pileup_mult=2)
  gini, cv (secondary, reported for comparison)
Validation (synthetic, python): median-normalization makes ALL three size-robust
(small vs big match at equal paralog level); frac>2x separates clean(2%) vs
collapsed(15%) ~7.5x (best), gini 2.3x, cv 1.5x, MAD too dull (dropped).
awk = computation only (two awks: depth->per-contig-mean, then stats), no multi-
line strings (anchor-bug lesson).
provisional_rank.R: reads new 7-col cv (id c1 c2 sim frac_hi gini cv); agg_rank now
= mean(rank_nb, rank_frac_hi, rank_anchor) — frac_hi is the uniformity signal in the
aggregate; gini/cv reported but NOT aggregated (they cancelled). main.nf +cv_pileup_mult=2.
INSTALL 5a: modules/compute_coverage_cv.nf, r_scripts/provisional_rank.R,
  workflows/optimize_denovo.nf (pileup_mult arg), main.nf (cv_pileup_mult).
RUN: still 2-pass for now (5a is isolated). compute_coverage_cv x84 RE-RUN (per-base
depth). Check cv rows: does frac_hi separate candidates AND is it size-robust (not
just tracking n_contigs)? Then 5b = SNP module + 1-pass collapse.

## write_anchor_calc — awk newline-escaping bug FIXED

## write_anchor_calc — awk newline-escaping bug FIXED
First run of the refactored module failed: awk printf strings contained \n that
became REAL newlines through the Groovy-triple-quote -> bash -> awk layers ->
"unterminated string" (a string literal can't span lines in awk).
FIX: awk does ARITHMETIC ONLY now — prints `key=value` lines (no newline escapes
anywhere in awk). The human-readable report is built by a plain bash heredoc
(<<EOF) where newlines are literal. Shell loads awk's values via
  while IFS='=' read k v; do eval "CV_$k=$v"; done < calc_vars.txt
and the heredoc references them as \${CV_...} (Groovy-escaped so bash expands at
runtime). Display numbers formatted in awk via sprintf (%.0f sites, %.3e/%.5f
rates) — safe, no newlines. Validated end-to-end in bash: 652 loci, value file
correct, report renders clean. Reinstall modules/write_anchor_calc.nf and -resume.

## write_anchor_calc — calculation moved INTO the process (refactor)

## write_anchor_calc — calculation moved INTO the process (refactor)
The anchor arithmetic was in Groovy (optimize_denovo.nf) with the process a dumb
echo of a pre-baked string. Moved the WHOLE calculation into write_anchor_calc.nf
(awk) so the model is self-contained, testable, single-source-of-truth.
- write_anchor_calc.nf now takes raw vals (enabled flag, genome_size, two site
  lengths, insert min/max, enzyme_pair) and DOES the math in awk:
  rare/common from max/min site len, n_rare=G/4^rare, common_rate=1/4^common,
  p_one=exp(-rate*imin)-exp(-rate*imax), p_win=1-(1-p_one)^2, loci=round(n_rare*p_win).
  Emits expected_loci.value (int or "NA") + expected_loci_calculation.txt (published).
  awk validated vs python: 652 loci for 1.3Gb/SbfI8/EcoRI6/insert150-221 (exact).
- optimize_denovo.nf: Groovy now ONLY sets anchor_enabled boolean + passes raw
  params; reads expected_loci from write_anchor_calc.out.expected_loci.map{f->text}.
  provisional_rank now gets expected_loci_ch (value channel from file) not a Groovy
  string; rank R drops anchor on "NA" as before.
Install: workflows/optimize_denovo.nf + modules/write_anchor_calc.nf (this round).

## GRID CEILING fix + expected-loci output (after first 3c run)

## GRID CEILING fix + expected-loci output (after first 3c run)
First 3c run (real ddRAD params: 1.3Gb, SbfI-8/EcoRI-6, insert 150-221) results:
- NB plot: clean 2-component fit, error mu=1.3 / locus mu=9.0, crossover=5.
  Locus component visually tiny (most UNIQUE SEQS are errors; loci hold the READS).
- CV/Gini DID discriminate (cv ~4.7-5.8, gini ~0.79-0.82) — 3c goal met, signals
  vary across c2/sim. Top candidate c4_k3_s0.90 (interior, not size-extreme). Good.

PROBLEM Jason caught: NB crossover=5 == grid c1 ceiling (knee capped at 5 in
assembly_diagnostics.R). So rank_nb is degenerate — it can only ever favor c1=5
(the highest available = the NB value), can't bracket an optimum that might be
higher. rank_nb is also coarse: just 4 c1-levels -> 4 tied ranks (3.5/9.5/15.5/
21.5). Circular: NB asked to discriminate within a range truncated at the NB value.

FIX: explicit grid CEILINGS decoupled from the auto-detect knee.
  NEW params: cutoff1_ceiling=8, cutoff2_ceiling=6 (main.nf). Grid is now
  floor..ceiling (pure param range); diagnostics knee still computed for the
  report + as NB fallback, but NO LONGER bounds the optimization grid.
  optimize_denovo.nf: c1_vals/c2_vals = Channel.value(floor..ceiling), no longer
  read cutoff*_value for the grid.
  Grid: c1{2..8}x c2{3..6} x sim{.85,.90,.95} = 7x4x3 = 84 candidates (was 24).
  Jason chose widest bracketing (84) over tighter options. Assembly is cached
  after first run so the cost is paid once.

EXPECTED-LOCI OUTPUT (Jason asked where the anchor calc goes): NEW module
write_anchor_calc.nf writes expected_loci_calculation.txt to
denovo_assembly/optimize/ — full arithmetic (rare/common sites, per-side window
prob, P_window, expected loci) or a DISABLED note listing missing params. Wired
in optimize_denovo.nf; was previously only in the nextflow log.

INSTALL (this round): main.nf, workflows/optimize_denovo.nf,
modules/write_anchor_calc.nf (new), modules/fit_nb_mixture.nf + r_scripts/
fit_nb_mixture.R (NB plot, from prior turn). Run with:
  --genome_size_est 1.3e9 --size_select_min 150 --size_select_max 221
  --enzyme1_site_len 8 --enzyme2_site_len 6
Expect 84 candidate assemblies (fresh; grid changed). Then re-inspect
provisional_rank.tsv: with c1 now 2..8, does rank_nb bracket an INTERIOR c1
(i.e. is the best c1 < 8)? And does anchor (~650 expected) add info vs rank_nb?

## Anchor model FIX + NB plot (pre-3c-test)

## Anchor model FIX + NB plot (pre-3c-test)
Jason's real ddRAD params: genome ~1.3 Gb, SbfI (8-cutter) + EcoRI (6-cutter),
insert ~150-221 bp (NOT the 280-350 bp fragment-trace window, which includes
adapters — use the INSERT).
BUG in old anchor model: 2*cuts1*cuts2/G * (win/mean_frag) underflowed to 0 for
realistic genomes and wasn't rare-cutter anchored.
FIX (rare-cutter-anchored exponential-spacing model): a ddRAD locus = a RARE-
cutter site (larger recognition seq = SbfI 8) whose NEAREST COMMON-cutter site
(EcoRI 6) falls in the insert window. Common spacing ~Exp(1/4^common_len);
P_one_side = exp(-rate*min) - exp(-rate*max); P = 1-(1-p)^2;
expected_loci = n_rare_sites * P. For Jason's params -> ~650 loci (validated in
python). Model auto-assigns larger site_len as rare; param order doesn't matter.
NOTE: candidates have 48k-93k contigs (~100x expected) — anchor will pull toward
aggressively-pruned (high-c1) candidates, so it may partly reinforce NB; watch
whether it adds independent info.
Jason's run should add: --genome_size_est 1.3e9 --size_select_min 150
  --size_select_max 221 --enzyme1_site_len 8 --enzyme2_site_len 6

NB PLOT added: fit_nb_mixture.R now emits nb_mixture_plot.png (observed coverage
histogram + both fitted NB component densities scaled to counts + red crossover
line). Placeholder plot on fallback/failure. Process emits+publishes .plot to
denovo_assembly/optimize/. optimize_rank env already has ggplot2 (via r-tidyverse)
— no config change. (Report wiring = chunk 7.)

## CHUNK 3c — cheap-signal rework (diagnosed from first 3b run)

## CHUNK 3c — cheap-signal rework (diagnosed from first 3b run)
PROBLEM (first 3b run): cheap stage did not discriminate. provisional_rank.tsv
showed redundancy=0 for ALL 24 candidates (self-cluster at 0.98 on already-
CD-HIT'd refs finds nothing -> DEAD signal), rank_anchor=NA (no enzyme params),
and rank_inflect + rank_nb are BOTH essentially functions of c1 pulling opposite
directions -> averaging cancels -> near-tied arbitrary winner. Three of four
signals were c1-functions; only redundancy was meant to be orthogonal and it was
dead.
NB-mixture itself worked great: low mu=1.27 (errors), high mu=8.96 (loci),
crossover cutoff1=5 (genuine 2-component fit; matched knee by coincidence).

FIX (chunk 3c):
- DROP redundancy (1b): inherently ~0 for CD-HIT'd refs.
- DROP paired-locus fraction (option 4): all contigs paired by construction
  (N-spacer joins fwd/rev before assembly) -> constant.
- DROP length-distribution (option 3): ddRAD contigs are read-length-governed
  (~2x read_len) -> flat across grid.
- DEMOTE inflection (1a) to REPORT-ONLY (kept as dist_inflect column, EXCLUDED
  from agg_rank; duplicates NB on c1 axis).
- KEEP NB cutoff1 (5b) as the c1-axis authority.
- ADD coverage-uniformity (option 2): map a seeded cv_sample_n subset (default 4)
  to ALL candidates in STAGE 1; compute CV and Gini of per-contig depth. This is
  the QUALITY signal that discriminates along c2/similarity (where NB is blind).
  Paralog collapse -> fat high-depth tail -> high CV/Gini.
  Depth proxy: samtools idxstats read-count per contig, normalized reads-per-bp
  (so contig-length variation isn't mistaken for depth non-uniformity).
  NOTE: idxstats is a proxy; upgrade to per-base `samtools depth` if CV is noisy.
- This MOVES the stage line: STAGE 1 now includes light CV mapping on all
  candidates; STAGE 2 = heavy bcftools on survivors only.

Ranking signals now ORTHOGONAL: NB(c1) + CV + Gini (c2/sim quality) + anchor(opt).
awk CV+Gini logic validated vs python: uniform->0/0, [1,1,1,1,100]->CV=1.9038/
Gini=0.7615 (exact match).

NEW: modules/compute_coverage_cv.nf (label map_reads), param cv_sample_n=4.
REVISED: compute_cheap_signals.nf (redundancy removed, 7-col output),
  provisional_rank.R (NB+CV+Gini+anchor; inflection report-only; reads cheap+cv,
  joins on id), provisional_rank.nf (+cv_rows input), optimize_denovo.nf (CV
  subset+wiring), main.nf (cv_sample_n replaces optimize_redundancy_identity).
Validated: rank aggregation balance + column contracts (cheap 7, cv 6) + cv_in
channel arity (each candidate paired with flat 8-file subset; glob finds 4 R1).

INSTALL (chunk 3c):
  modules/compute_coverage_cv.nf   (new)
  modules/compute_cheap_signals.nf (overwrite — redundancy removed)
  modules/provisional_rank.nf      (overwrite — +cv_rows)
  r_scripts/provisional_rank.R     (overwrite — new signal set)
  workflows/optimize_denovo.nf     (overwrite — CV wiring)
  main.nf                          (overwrite — cv_sample_n)
RUN: cached thru assemble_rainbow (24). NEW: compute_coverage_cv x24 (light map
of 4-sample subset to each candidate), compute_cheap_signals x24 (re-run, now
7-col), provisional_rank x1. Inspect denovo_assembly/optimize/provisional_rank.tsv
-> CV/Gini columns should VARY across c2/sim (the whole point); agg_rank should
pick interior, not size-extreme.

## CHUNK 3b — outputs now published (added after first 3b run)

## CHUNK 3b — outputs now published (added after first 3b run)
First 3b run succeeded but outputs were invisible (no publishDir — they sat in
work/ only). Added publishDir to ${outdir}/denovo_assembly/optimize/:
  - provisional_rank.nf  -> provisional_rank.tsv
  - fit_nb_mixture.nf    -> nb_mixture_fit.txt
(compute_cheap_signals per-candidate rows stay unpublished — intermediate.)
To inspect the FIRST run without rerunning, pull from work dirs:
  cat $(find work/05/9696eb* -name provisional_rank.tsv)
  cat $(find work/b1/ec738f* -name nb_mixture_fit.txt)
Reinstall provisional_rank.nf + fit_nb_mixture.nf and -resume to publish them
(these two processes will re-run since their definition changed; cheap.

## CHUNK 3b — install + test
Install:
  modules/compute_cheap_signals.nf   (new)
  modules/fit_nb_mixture.nf          (new)
  modules/provisional_rank.nf        (new)
  r_scripts/fit_nb_mixture.R         (new)
  r_scripts/provisional_rank.R       (new)
  modules/assembly_diagnostics.nf    (overwrite — emits coverage_freq.txt)
  workflows/optimize_denovo.nf       (overwrite — real stage-1)
  main.nf                            (overwrite — + optimize_redundancy_identity)

Config: ensure withLabel:optimize_rank exists (r-base + tidyverse) — used by
fit_nb_mixture and provisional_rank. compute_cheap_signals reuses denovo_assembly
label (needs cd-hit-est).

Run: same command. Cached up through assemble_rainbow_candidate (24). NEW work:
  - compute_cheap_signals: 24 tasks (self-cluster each candidate at 0.98)
  - fit_nb_mixture: 1 (global NB on coverage_freq)
  - provisional_rank: 1 (rank R script) -> survivors.txt
Then stage2/aggregate/finalize STILL STUBBED -> downstream mapping on stub winner.
Check: provisional_rank.tsv has per-signal ranks + agg_rank; survivors plausible
(interior candidates, not pure size-extremes). nb_mixture_fit.txt shows 2-comp fit.

To test the anchor (signal 2), add e.g.:
  --genome_size_est 1.2e9 --size_select_min 300 --size_select_max 500
  (and optionally --enzyme1_site_len 6 --enzyme2_site_len 4)

## Topology test — what to expect (STUB run)

## Topology test — what to expect (STUB run)
Command (after install; upstream cached so this is cheap):
  nextflow run gcl_illumina_qc/main.nf -profile slurm -resume \
    --reads "data/fq_raw/guttata*.{1,2}.fq.gz" --assembly_mode denovo \
    --sequencing_type ddrad --decontam_conffile configs/decontam.conf \
    --outdir "illumina_qc/guttata_denovo" \
    --do_optimize --n_pseudo_reps 6 --optimize_sample_pct 25

Expected STUB behavior (NOT real results — just topology validation):
  - extract_unique_seqs runs on all samples (subset stubbed to "use all")
  - assembly_diagnostics runs (real knee/curves)
  - 2 stub candidates built (hardcoded grid c2_k2_s0.90, c3_k3_s0.90)
  - stub_cheap_signals x2 -> stub_provisional_rank -> survivors = up to 3 ids
  - split_pseudo_rep runs on first 6 samples (real split, deterministic)
  - stub_stage2 on survivors, stub_aggregate -> best_id, stub_finalize_reference
  - downstream mapping/report run on the STUB reference (1 tiny contig) — will look
    nonsensical (near-zero mapping) — that's FINE, we're testing wiring not results.

CRITICAL CHECK during/after the run:
  - Confirm pseudo-rep halves (`*_a.r1.fq.gz`, `*_b.*`) appear ONLY in work/ dirs,
    NEVER in illumina_qc/guttata_denovo/cleaned_reads/ (production output).
    Run: ls illumina_qc/guttata_denovo/cleaned_reads/ | grep -E '_(a|b)\.r[12]' 
    -> should return NOTHING.

## Resolved design review notes (chunk 2)
- best_id passed to finalize as a FILE (not val) to avoid val/path cardinality mix.
- candidates queue channel forked to two consumers (survivor join + all-fastas
  collect) — Nextflow auto-forks queue channels to multiple consumers; OK.
- cleaned_reads (workflow take: input) fed to extract_unique_seqs, .take(N) for
  concordance, etc. — each operator forks; safe to reuse.
- split_pseudo_rep seed passed as plain value -> auto-wrapped as value channel input.

## Open notes
- split_pseudo_rep needs bash for `<(...)`; Nextflow uses /bin/bash by default. OK.
- The NB-mixture (5b) uses the glmmTMB/NB stack already in r_analysis env.
- Locus-count model: cut-site freq = 1/4^site_len per enzyme; expected fragments
  in [size_select_min,max] from genome_size_est. Refine once real inputs arrive.
- STUB reference is tiny -> chunk 2 test will show broken mapping %; ignore until
  chunk 6 wires the real finalize assembly.

---

## CHUNK 5b — 1-PASS SNP-signal rework (BUILT, not yet installed/run)

### Why (the coverage saga conclusion)
Every coverage-uniformity statistic tried (idxstats CV/Gini, per-base frac_hi/Gini/CV)
stayed correlated with assembly SIZE on real data — latest run: corr(n_contigs,
frac_hi) = -0.80. Coverage uniformity reads pruning-aggressiveness, not an independent
quality axis, so it can't select OR gate. DECISION: drop coverage entirely (no mapping
for it) and go 1-PASS — run SNP recovery on ALL candidates. The 2-stage gate
(provisional_rank -> survivors -> stage2) is removed; cheap signals are ranking inputs.

### Final 1-pass signal set (per candidate, all 84)
1. NB cutoff1 proximity   |c1 - nb_cutoff1|        (c1 axis)
2. anchor proximity       |n_contigs - expected|   (size/biological target; ~652)
3. concordance (PRIMARY)  depth-weighted genotype agreement between pseudo-rep halves
4. r80_loci               STACKS r80: # contigs polymorphic & genotyped in >=80% of
                          SNP-subset samples (locus-level; size-robust — junk loci fail
                          the share test). threshold = params.r80_threshold (0.8).
5. snps_per_locus         n_snps / n_contigs (density, secondary)
   inflection             REPORT-ONLY (excluded from aggregate; duplicates NB)
Aggregation: weight-free MEAN of available signal ranks; winner = min agg rank.

### r80 = STACKS rule (confirmed via web: Paris et al 2017, catchenlab)
LOCUS-level, not site-level: a polymorphic locus (contig with >=1 SNP) genotyped in
>=80% of individuals. Chosen because loci shared across 80% are likely real; low-cov/
junk/paralog loci fail the share test (the size-robustness we need).

### Files (in /mnt/user-data/outputs)
NEW  modules/compute_snp_signals.nf   (label stage2_call) — runs on ALL candidates;
       index + PASS A (map each pseudo-rep half separately, bcftools mpileup -a FORMAT/DP
       | call -mv, query GT+DP, awk depth-weighted concordance) + PASS B (map snp subset,
       joint call, awk r80 locus count + snps_per_locus). Emits 9-col:
       id c1 c2 sim concordance r80_loci snps_per_locus n_snps n_called_contigs.
       publishDir optimize/snp_signals/.
NEW  r_scripts/rank_and_select.R + modules/rank_and_select.nf (label optimize_rank) —
       single weight-free aggregation. Reads cheap_all.tsv(7) + snp_all.tsv(9 ciidddddi).
       Ranks nb/anchor (asc on |.|) + conc/r80/snpdens (desc). agg=rowMeans(available).
       Emits final_rank.tsv + best_id.value + optimize_plot.png (top-15 faceted).
       publishDir optimize/.
EDIT workflows/optimize_denovo.nf — removed compute_coverage_cv + provisional_rank
       includes/wiring + stub_stage2/stub_aggregate processes + survivors/splitText/gate.
       New flow: cheap_signals -> fit_nb_mixture -> write_anchor_calc -> split_pseudo_rep
       -> pseudo_fastqs bag + snp_subset_reads bag (snp_sample_pct, sorted-before-take,
       min 2) -> snp_in = candidates.combine(pseudo).combine(snp) -> compute_snp_signals
       -> rank_and_select -> finalize selects winner by best_id join.
       emit optimize_summary=final_rank, optimize_plot=plot (denovo_assembly.nf consumes
       these names — confirmed unchanged).
EDIT main.nf — removed cv_sample_n/cv_pileup_mult; +r80_threshold=0.8; stage2_min/max
       params removed (1-pass = no gate).
DELETE/stop-using modules/compute_coverage_cv.nf, modules/provisional_rank.nf,
       r_scripts/provisional_rank.R, modules/compute_coverage_cv (+ any 5a per-base
       coverage module).

### Validated before shipping
- concordance join awk: hand-calc 60/74 = 0.8108 -> awk gives 0.8108 ✓
- concordance formula (python): clean ref 0.978 vs collapsed 0.756 ✓ (separates)
- r80 count awk: ctg1(0.8 PASS)+ctg3(1.0 PASS), ctg2(0.2 FAIL) -> 2 ✓
- all files brace/paren/bracket balanced ✓
- column contract: SNP emit 9 == R read 9 (ciidddddi) ✓
- stage2_call conda has bwa-mem2+samtools+bcftools (cpus=8) ✓

### To INSTALL then RUN
Install the 5 files above; delete/stop-using the coverage+provisional modules.
Run (full biological params):
  nextflow run gcl_illumina_qc/main.nf -profile slurm -resume \
    --reads "data/fq_raw/guttata*.{1,2}.fq.gz" --assembly_mode denovo \
    --sequencing_type ddrad --decontam_conffile configs/decontam.conf \
    --outdir "illumina_qc/guttata_denovo" --do_optimize --n_pseudo_reps 6 \
    --genome_size_est 1.3e9 --size_select_min 150 --size_select_max 221 \
    --enzyme1_site_len 8 --enzyme2_site_len 6
Cached through assemble_rainbow (84). NEW heavy work: compute_snp_signals x84 (2 bcftools
passes each — the accepted compute jump), rank_and_select x1.

### Inspect after run
- final_rank.tsv: do concordance/r80/snpdens DISCRIMINATE and are they size-orthogonal
  (corr each vs n_contigs ~ 0, unlike the -0.80 coverage failures)? Does winner make
  biological sense vs ~652 anchor? Are concordance & r80 redundant (may prune one)?

---

## PIVOT — grid model + r80-vs-n_contigs ELBOW selector (BUILT, not yet installed/run)

### Why (the conclusion from the 1-pass SNP run)
The 84-candidate 1-pass run showed EVERY signal is monotone in assembly size:
corr(n_contigs, concordance)=-0.92, r80=+0.83, snps_per_locus=-0.94; the three SNP
signals are mutually collinear (|r|0.95-0.98). So genotype signals are NOT size-
orthogonal either — the premise behind the 1-pass redesign was wrong on real data.
ROOT CAUSE: the cutoff grid (c1,c2) is a pure STRINGENCY dial — no interior optimum
exists to find; every quality stat just reads pruning-aggressiveness = size.

WHAT DOES have a turning point: r80 vs n_contigs. r80 rises steeply then PLATEAUS
(real data: marginal r80/1000-contigs decays 583->146->93->28->18->~0, slightly neg
at the largest assemblies). STACKS optimizes exactly this (Paris 2017): maximize
broadly-shared polymorphic loci, pick the plateau. The dDocent analogue = the ELBOW
of r80(n_contigs). Validated: Kneedle on the r80 upper-envelope -> elbow at ~27,498
contigs (matches the visual ~25k read off Jason's plot); winner among <=elbow by
max-r80 then fewest-contigs = c8_k3_s0.95 (27,498 contigs, r80=3403). Interior, not
an edge. (~42x the 652 anchor, expected for a raw dDocent ref; downstream filters.)

Param research (dDocent docs, ddocent.com): initial CD-HIT similarity is DELIBERATELY
loose (0.8) and a minor knob — confirmed on data (sweeping it moved outputs ~4% vs
~25% for one c1 step). The real STACKS-M-analogue tradeoff (collapse paralogs vs split
alleles) lives in Rainbow div -f / merge -r and the FINAL CD-HIT similarity. So the
sweep PIVOTS off cutoffs onto those.

### New design
GRID MODEL: 6 axes, each a SCALAR (fixed) or LIST (swept); grid = Cartesian product.
  Defaults: cutoff1=null(->NB crossover), cutoff2=3, cluster_similarity(init)=0.8 [all
  pinned]; final_similarity=[0.90,0.95], div_f=[0.1,0.2,0.5], merge_r=2 [swept defaults].
  div_K=10 fixed. max_grid_candidates=64 guardrail (warns).
SELECTION: r80-vs-n_contigs ELBOW. Kneedle on the r80 upper-envelope (param-free), or
  marginal-slope override via r80_elbow_min_slope (default 30, used only if >0). Winner
  = among candidates at/below the elbow contig count, MAX r80, ties -> FEWEST contigs.
  concordance/anchor/NB/snps_per_locus REPORTED but do NOT drive (concordance is
  size-monotone). Parameter-agnostic: operates on the OUTCOME (size<->r80), so it works
  no matter which axes were swept.
SNP sampling: snp_sample_pct default 75% (was 25) — r80 stability scales with samples.

### Files (in /mnt/user-data/outputs)
EDIT main.nf — replaced floor/ceiling cutoff params with scalar-or-list axis model
  (cutoff1=null, cutoff2=3, cluster_similarity=0.8, final_similarity=[.90,.95],
  div_f=[.1,.2,.5], merge_r=2, div_K=10); +max_grid_candidates=64, +r80_elbow_min_slope=30,
  snp_sample_pct 25->75. Removed duplicate old cutoff block.
EDIT denovo_assembly.nf — single-assembly (non-optimize) path now coerces any list
  axis to its first element (firstOf helper) so scalar-or-list params don't break it.
  Updated optimize log line.
EDIT optimize_denovo.nf — MAJOR: fit_nb_mixture moved BEFORE grid construction (c1=null
  resolves to NB value via nb_c1_ch.flatMap). asList() normalizes each axis; full
  Cartesian product into grid_metas carrying [c1,c2,isim,divf,mr,fsim,id_cutoff,id];
  grid-size warn. Filter runs ONCE per distinct (c1,c2) (cutoff_pairs.unique{id_cutoff}),
  re-expanded via combine(by:0) [correct vs join: left has dup keys, right has 1/key].
  assemble_rainbow_candidate(assemble_in, params.div_K). rank_and_select now takes
  min_slope arg. finalize echo uses new meta fields.
EDIT assemble_rainbow_candidate.nf — reads ALL assembly params from meta (isim/divf/mr/
  fsim); only div_K is a global arg. id now encodes all 6 axes. Stats heredoc updated.
EDIT compute_cheap_signals.nf + compute_snp_signals.nf — 4th row col (sim) now = meta.fsim
  (was meta.sim, which no longer exists). Row schemas otherwise unchanged (cheap 7, snp 9).
REWRITE r_scripts/rank_and_select.R — ELBOW selector. Reads cheap_all(7)+snp_all(9).
  Builds r80 upper-envelope vs n_contigs, Kneedle elbow (or min_slope override), winner=
  fewest-contigs at plateau. Emits final_rank.tsv + best_id.value + optimize_plot.png
  (elbow curve, winner gold-circled) + optimize_params_plot.png (r80 faceted vs each
  SWEPT param — shows which knob moves r80). r80-absent fallback = fewest contigs.
EDIT rank_and_select.nf — +min_slope input, +optimize_params_plot.png output.

### Validated (python ports; no R/bioconda in sandbox)
- elbow Kneedle-on-envelope = 27,498 contigs on real data (matches visual ~25k) ✓
- winner selection (max-r80<=elbow, tie fewest-contigs) = c8_k3_s0.95, interior ✓
- all 8 edited files brace/paren/bracket balanced ✓
- no stale params (floor/ceiling/cv_*/stage2_* gone); all new params defined once ✓
- combine(by:0) correct for dup-key-left/single-key-right (per NF docs) ✓
- filter module meta fields (c1,c2,id_cutoff) subset of cutoff_pairs meta ✓
- column contracts: cheap 7-col, snp 9-col unchanged; R readers match ✓

### To INSTALL then RUN
Overwrite: main.nf, workflows/denovo_assembly.nf, workflows/optimize_denovo.nf,
  modules/assemble_rainbow_candidate.nf, modules/compute_cheap_signals.nf,
  modules/compute_snp_signals.nf, modules/rank_and_select.nf, r_scripts/rank_and_select.R
Run (NB-pinned c1, c2=3, sweep div_f x final_sim => ~12 candidates default):
  nextflow run gcl_illumina_qc/main.nf -profile slurm -resume \
    --reads "data/fq_raw/guttata*.{1,2}.fq.gz" --assembly_mode denovo \
    --sequencing_type ddrad --decontam_conffile configs/decontam.conf \
    --outdir "illumina_qc/guttata_denovo" --do_optimize --n_pseudo_reps 6 \
    --genome_size_est 1.3e9 --size_select_min 150 --size_select_max 221 \
    --enzyme1_site_len 8 --enzyme2_site_len 6
NOTE: grid changed (assembly params now in meta.id) => assemble_rainbow_candidate work
  dirs will NOT cache-hit the old 84-candidate run; expect fresh assemblies. Default
  grid is only ~12 candidates though (cheaper than the 84-run).

### Inspect after run
- optimize_plot.png: does r80-vs-contigs show the plateau + is the gold winner at the
  elbow (not an edge)? optimize_params_plot.png: does div_f actually move r80, or is it
  inert like initial-sim was? (If div_f is flat, the real lever is elsewhere.)
- Is the winner's size defensible for downstream genotyping (enough shared loci, not
  pure junk)? With only ~12 candidates the elbow may be coarse — widen grid if so.

---

## PIVOT FIXES (round 2) — after the 6-candidate pinned-grid run

Run with pinned c1/c2 produced only 6 candidates over a 4% size range (46.1k-48.0k
contigs). Three problems surfaced; all fixed.

1. ELBOW NEEDS A SIZE-SPANNING GRID. With c1/c2 pinned, the r80-vs-contigs "curve"
   was a short monotone rising line with NO interior elbow -> Kneedle picked an
   arbitrary endpoint (the winner flipped between top/bottom depending on tie-break
   noise). FIX: default now SWEEPS the size-driving axes again — cutoff1=[2..8],
   cutoff2=[3,4] — plus div_f=[.1,.2,.5], final_similarity=[.90,.95]. = 84 candidates
   (max_grid_candidates raised 64->100). The cutoff sweep GENERATES the curve; the
   r80 elbow is still the SELECTOR. Any axis set to a SCALAR is pinned (not swept) —
   pin-any-axis preserved. (cluster_similarity=0.8 and merge_r=2 pinned by default.)
   Confirmed: div_f DOES move r80 (0.5->0.1 raised r80 2789->2967 / 2829->3173), but
   along the size axis — which is fine, the elbow is size-aware.

2. SWEPT PARAMS WEREN'T COLUMNS -> params_plot showed only final_sim facet. FIX:
   compute_cheap_signals + compute_snp_signals now emit ALL SIX axes as columns.
   New schemas: cheap = id c1 c2 isim divf mr fsim n_contigs total_len mean_len (10,
   "ciiddidiii"); snp = id c1 c2 isim divf mr fsim concordance r80_loci snps_per_locus
   n_snps n_called_contigs (12, "ciiddidddddi"). final_rank.tsv now writes c1,c2,isim,
   divf,mr,fsim. rank_and_select.R: params_plot facets over EVERY varying axis (colored
   by n_contigs to expose size-confounding); elbow plot encodes the swept-param combo
   via color (prefers divf/fsim) + shape (when grid <=40 combos), single color above.

3. NB PLOT AXIS BROKEN (image 1: "1 0 0 0" labels, histogram a thin band). Cause:
   scale_y_continuous(comma_format + log10_trans) on counts that hit 0 baseline ->
   log10(0)=-Inf, degenerate breaks. FIX in fit_nb_mixture.R: floor n and component
   densities at 1, scale_y_log10 with trans_breaks/trans_format(math_format(10^.x)) +
   annotation_logticks, coord ylim=c(1,ymax). NOTE: even fixed, the two NB components
   barely separate on this data (locus comp is 17% mass, shallow shoulder) — the
   cutoff1=5 VALUE is sound (matches knee) but the PICTURE is inherently muddy here.

### Validated (python ports; no R in sandbox)
- new col_types lengths match (cheap 10, snp 12); type assignment correct ✓
- aesthetic-axis detection: divf+fsim flagged varying -> 2 facets (not 1) ✓
- elbow on 6-pt pinned grid is UNSTABLE (endpoint flip) = confirms why wide grid needed ✓
- all touched files balanced ✓; default grid 84 < guardrail 100 ✓

### Files changed this round
main.nf (cutoff1/cutoff2 back to wide swept defaults; max_grid 100),
compute_cheap_signals.nf + compute_snp_signals.nf (6-axis columns),
r_scripts/rank_and_select.R (wide readers, full faceting, param-combo aesthetics),
r_scripts/fit_nb_mixture.R (log-axis fix).

### RERUN
Same command as before (genome_size_est/size_select/enzyme params). Default now 84
candidates spanning ~5k-90k contigs so the r80 elbow is interior & stable (~27k from
the earlier full run). Inspect optimize_params_plot.png: which axis (c1? c2? div_f?
fsim?) actually moves r80 vs just tracking n_contigs color.
