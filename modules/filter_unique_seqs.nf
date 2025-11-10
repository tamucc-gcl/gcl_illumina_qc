// modules/filter_unique_seqs.nf
process filter_unique_seqs {
    label 'denovo_assembly'
    tag "filtering"
    
    //publishDir "${params.outdir}/denovo_assembly/filtered", mode: 'copy'
    
    input:
        path(uniq_seq_files)
        val(cutoff1)  // Minimum count per individual (default: 4)
        val(cutoff2)  // Minimum number of individuals (default: 4)
    
    output:
        path("uniq.fasta"), emit: filtered_fasta
        path("totaluniqseq"), emit: totaluniqseq
        path("filter_stats.txt"), emit: stats
    
    script:
    """
    echo "Starting sequence filtering with cutoff1=${cutoff1}, cutoff2=${cutoff2}"
    
    # Combine all unique sequences
    cat *.uniq.seqs > uniq.seqs
    echo "Total unique sequences across all samples: \$(wc -l < uniq.seqs)"
    
    # Keep only sequences present at least cutoff1 times per individual
    echo "Filtering sequences with minimum ${cutoff1} reads per individual..."
    
    # Fix: Use single quotes around the entire mawk command to prevent shell expansion
    parallel --no-notice -j ${task.cpus} 'mawk -v x=${cutoff1} "\\\$1 >= x" {} | cut -f2' ::: *.uniq.seqs | \\
        perl -e 'while (<>) {chomp; \$z{\$_}++;} while((\$k,\$v) = each(%z)) {print "\$v\\t\$k\\n";}' > uniqCperindv
    
    echo "Sequences passing per-individual filter: \$(wc -l < uniqCperindv)"
    
    # Second filtering based on how many individuals have the sequence
    echo "Filtering sequences present in at least ${cutoff2} individuals..."
    mawk -v x=${cutoff2} '\$1 >= x' uniqCperindv > uniq.k.${cutoff2}.c.${cutoff1}.seqs
    
    echo "Sequences passing both filters: \$(wc -l < uniq.k.${cutoff2}.c.${cutoff1}.seqs)"
    
    # Turn back into fasta
    cut -f2 uniq.k.${cutoff2}.c.${cutoff1}.seqs > totaluniqseq
    mawk '{c = c + 1; print ">Contig_" c "\\n" \$1}' totaluniqseq > uniq.fasta
    
    # Generate filtering statistics
    cat > filter_stats.txt <<-EOF
	Filtering Statistics
	====================
	Cutoff 1 (min reads per individual): ${cutoff1}
	Cutoff 2 (min individuals): ${cutoff2}
	
	Total unique sequences: \$(wc -l < uniq.seqs)
	Sequences after per-individual filter: \$(wc -l < uniqCperindv)
	Sequences after individual count filter: \$(wc -l < uniq.k.${cutoff2}.c.${cutoff1}.seqs)
	Final sequences for assembly: \$(grep -c '^>' uniq.fasta)
	EOF
    
    cat filter_stats.txt
    """
}