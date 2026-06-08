// workflows/assembly_sweep.nf
// Joint sweep over cutoff1 x cutoff2 x cluster_similarity.
// - cutoff grids run floor..knee (knees from assembly_diagnostics) unless a
//   cutoff is pinned via params; cluster_similarity grid from params.
// - each (cutoff1,cutoff2) is filtered; each filtered set assembled at each sim;
//   a percentage subset of samples is mapped back to score every candidate;
//   the winning reference is selected without re-assembly.
//
// Replaces cluster_similarity_sweep.nf.

nextflow.enable.dsl = 2

include { filter_unique_seqs_sweep } from '../modules/filter_unique_seqs_sweep.nf'
include { assemble_rainbow_sweep }   from '../modules/assemble_rainbow_sweep.nf'
include { sweep_map_score }          from '../modules/sweep_map_score.nf'
include { score_assemblies }         from '../modules/score_assemblies.nf'

workflow assembly_sweep {
    take:
        all_uniq_seqs    // collected list of *.uniq.seqs (single emission)
        ceil1_value      // file: cutoff1 knee (grid ceiling)
        ceil2_value      // file: cutoff2 knee (grid ceiling)
        cleaned_reads    // tuple(sid, r1, r2) for all samples

    main:
        floor = params.cutoff_sweep_floor

        // cutoff value channels: pinned param -> single value; else floor..knee
        c1_vals = (params.cutoff1 != null)
            ? Channel.of(params.cutoff1 as int)
            : ceil1_value.map{ it.text.trim() as int }
                         .flatMap{ k -> def hi = (k < floor) ? floor : k; (floor..hi).toList() }
        c2_vals = (params.cutoff2 != null)
            ? Channel.of(params.cutoff2 as int)
            : ceil2_value.map{ it.text.trim() as int }
                         .flatMap{ k -> def hi = (k < floor) ? floor : k; (floor..hi).toList() }

        // cutoff combos + the shared uniq files -> (meta[c1,c2], files)
        filter_in = c1_vals.combine(c2_vals)
            .combine(all_uniq_seqs)
            .map{ c1, c2, files ->
                tuple([c1: c1, c2: c2, id_cutoff: "c${c1}_k${c2}"], files)
            }

        filter_unique_seqs_sweep( filter_in )

        // cross each filtered set with the cluster_similarity grid -> full candidate meta
        sims = Channel.fromList(params.sweep_cluster_similarity)
        asm_in = filter_unique_seqs_sweep.out.filtered
            .combine(sims)
            .map{ meta, fasta, tuniq, sim ->
                def m = meta + [sim: sim, id: "c${meta.c1}_k${meta.c2}_s${sim}"]
                tuple(m, fasta, tuniq)
            }

        assemble_rainbow_sweep(
            asm_in, params.div_f, params.div_K, params.merge_r, params.final_similarity
        )

        // percentage subset of samples (seeded), 100 = all
        frac = (params.sweep_sample_pct as double) / 100.0
        reads_subset = cleaned_reads.toList().flatMap { lst ->
            def sh = new ArrayList(lst)
            Collections.shuffle(sh, new Random(params.sweep_seed))
            def n = Math.max(1, (int) Math.ceil(lst.size() * frac))
            sh.take(n)
        }
        reads_flat = reads_subset.map{ sid, r1, r2 -> [r1, r2] }.flatten().collect()

        // attach the read list to every candidate reference -> (meta, ref, [reads])
        score_in = assemble_rainbow_sweep.out.reference.combine(reads_flat.map{ [it] })
        sweep_map_score( score_in )

        score_assemblies( sweep_map_score.out.score.collect() )
        best_id = score_assemblies.out.best_id.map{ f -> f.text.trim() }

        // select winning reference / stats from the already-built candidates
        winner_ref = assemble_rainbow_sweep.out.reference
            .combine(best_id)
            .filter{ meta, ref, bid -> meta.id == bid }
            .map{ meta, ref, bid -> ref }

        winner_stats = assemble_rainbow_sweep.out.stats
            .combine(best_id)
            .filter{ meta, st, bid -> meta.id == bid }
            .map{ meta, st, bid -> st }

    emit:
        reference      = winner_ref
        assembly_stats = winner_stats
        sweep_summary  = score_assemblies.out.summary
        sweep_plot     = score_assemblies.out.plot
        best_id        = best_id
}
