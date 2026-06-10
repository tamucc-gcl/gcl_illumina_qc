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
