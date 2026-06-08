// workflows/cluster_similarity_sweep.nf
// Sweep over a grid of cluster_similarity values: assemble a candidate reference
// for each, map a subset of cleaned reads back to score it, and return the
// winning reference (no re-assembly of the winner).

nextflow.enable.dsl = 2

include { assemble_rainbow_sweep } from '../modules/assemble_rainbow_sweep.nf'
include { sweep_map_score }        from '../modules/sweep_map_score.nf'
include { score_assemblies }       from '../modules/score_assemblies.nf'

workflow cluster_similarity_sweep {
    take:
        filtered_fasta   // single-emission: uniq.fasta
        totaluniqseq     // single-emission: totaluniqseq
        cleaned_reads    // tuple(sid, r1, r2) for all samples

    main:
        // Candidate grid (value channel broadcasts the fixed inputs to each)
        sims = Channel.fromList(params.sweep_cluster_similarity)

        // Pair each candidate sim with the (shared) filtered fasta + totaluniqseq
        asm_in = sims
            .combine(filtered_fasta)
            .combine(totaluniqseq)          // -> (sim, fasta, totaluniqseq)

        assemble_rainbow_sweep(
            asm_in,
            params.div_f,
            params.div_K,
            params.merge_r,
            params.final_similarity
        )

        // Subset of samples used only for scoring the candidates
        reads_subset = cleaned_reads.randomSample(params.sweep_n_samples, 42)
        reads_flat   = reads_subset
            .map{ sid, r1, r2 -> [r1, r2] }
            .flatten()
            .collect()                       // single emission: [all subset files]

        // Attach the read list to every candidate reference -> (sim, ref, [reads])
        score_in = assemble_rainbow_sweep.out.reference.combine(reads_flat.map{ [it] })

        sweep_map_score(score_in)

        // Rank candidates and pick the winner
        score_assemblies( sweep_map_score.out.score.collect() )

        best_sim = score_assemblies.out.best_sim.map{ f -> f.text.trim() }

        // Select the winning reference from the already-built candidates
        winner_ref = assemble_rainbow_sweep.out.reference
            .combine(best_sim)
            .filter{ sim, ref, best -> sim.toString() == best.toString() }
            .map{ sim, ref, best -> ref }

        winner_stats = assemble_rainbow_sweep.out.stats
            .combine(best_sim)
            .filter{ sim, st, best -> sim.toString() == best.toString() }
            .map{ sim, st, best -> st }

    emit:
        reference      = winner_ref
        assembly_stats = winner_stats
        sweep_summary  = score_assemblies.out.summary
        sweep_plot     = score_assemblies.out.plot
        best_sim       = best_sim
}
