process generate_report {
    label 'generate_report'
    tag "final_report"
    
    publishDir "${params.outdir}", mode: params.publish_dir_mode
    
    input:
        path qc_plot
        path read_summary
        path multiqc_reports
        val genome_source
        path initial_histogram
        path mapped_histogram
        path mapping_summary
        path insert_size_violin 
        path assembly_stats
        path filter_stats
        path species_blast_tsv
        path species_raw_pie
        path species_summary_pie
        path species_top_hits
        
    output:
        path "qc_pipeline_report.md"
        path "qc_pipeline_report.html"
    
    script:
    """
    # Check if mapping was performed
    if [ "${mapping_summary}" == "NO_MAPPING" ] || [[ "${mapping_summary}" == *"NO_MAPPING"* ]]; then
        echo "No mapping performed" > mapping_summary.txt
        export MAPPING_PERFORMED="false"
    else
        export MAPPING_PERFORMED="true"
    fi
    
    # Check if insert size violin exists and mapping was performed
    if [ "\$MAPPING_PERFORMED" == "true" ] && [ -f "${insert_size_violin}" ] && [[ "${insert_size_violin}" != *"NO_"* ]]; then
        export INSERT_SIZE_VIOLIN_PERFORMED="true"
    else
        export INSERT_SIZE_VIOLIN_PERFORMED="false"
    fi

    # Check if assembly stats exist
    if [[ "${assembly_stats}" == *"no_assembly"* ]] || [[ "${assembly_stats}" == *"NO_ASSEMBLY"* ]]; then
        export ASSEMBLY_PERFORMED="false"
    else
        export ASSEMBLY_PERFORMED="true"
    fi
    
    # Check if filter stats exist
    if [[ "${filter_stats}" == *"no_filter"* ]] || [[ "${filter_stats}" == *"NO_FILTER"* ]]; then
        export FILTER_PERFORMED="false"
    else
        export FILTER_PERFORMED="true"
    fi
    
    # Check if species ID was performed
    if [[ "${species_blast_tsv}" == *"NO_SPECIES"* ]] || [[ "${species_blast_tsv}" == *"no_species"* ]]; then
        export SPECIES_ID_PERFORMED="false"
    else
        export SPECIES_ID_PERFORMED="true"
    fi

    cat <<'PYEOF' > generate_report.py
#!/usr/bin/env python3

import re
import statistics
import os
import subprocess
import math
from pathlib import Path

# Parse genome source
genome_source = "${genome_source}"

# Check flags from environment
mapping_performed = os.environ.get("MAPPING_PERFORMED", "false") == "true"
insert_size_violin_performed = os.environ.get("INSERT_SIZE_VIOLIN_PERFORMED", "false") == "true"
assembly_performed = os.environ.get("ASSEMBLY_PERFORMED", "false") == "true"
filter_performed = os.environ.get("FILTER_PERFORMED", "false") == "true"
species_id_performed = os.environ.get("SPECIES_ID_PERFORMED", "false") == "true"

# Initialize variables
species_name = ""
reference_line = ""
assembly_section = ""

if genome_source.startswith("denovo:"):
    species_name = ""  # No species name for de novo
    
    # Read assembly statistics if available
    assembly_info = []
    
    if filter_performed:
        try:
            # Read filter stats
            with open("${filter_stats}", 'r') as f:
                filter_content = f.read()
                
                # Extract filtering parameters
                cutoff1_match = re.search(r'Cutoff 1.*?: (\\d+)', filter_content)
                cutoff2_match = re.search(r'Cutoff 2.*?: (\\d+)', filter_content)
                total_uniq_match = re.search(r'Total unique sequences: (\\d+)', filter_content)
                after_indiv_match = re.search(r'Sequences after per-individual filter: (\\d+)', filter_content)
                after_count_match = re.search(r'Sequences after individual count filter: (\\d+)', filter_content)
                final_seqs_match = re.search(r'Final sequences for assembly: (\\d+)', filter_content)
                
                if cutoff1_match and cutoff2_match:
                    assembly_info.append(f"**Filtering Parameters:**")
                    assembly_info.append(f"- Minimum reads per individual (cutoff1): {cutoff1_match.group(1)}")
                    assembly_info.append(f"- Minimum individuals required (cutoff2): {cutoff2_match.group(1)}")
                    assembly_info.append("")
                
                if total_uniq_match:
                    assembly_info.append(f"**Filtering Statistics:**")
                    if total_uniq_match:
                        assembly_info.append(f"- Total unique sequences across all samples: {int(total_uniq_match.group(1)):,}")
                    if after_indiv_match:
                        assembly_info.append(f"- Sequences after per-individual filter: {int(after_indiv_match.group(1)):,}")
                    if after_count_match:
                        assembly_info.append(f"- Sequences after individual count filter: {int(after_count_match.group(1)):,}")
                    if final_seqs_match:
                        assembly_info.append(f"- Final sequences for assembly: {int(final_seqs_match.group(1)):,}")
                    assembly_info.append("")
                    
        except Exception as e:
            print(f"Could not read filter stats: {e}")
        
        try:
            # Read assembly stats
            with open("${assembly_stats}", 'r') as f:
                content = f.read()
                
                # Extract assembly parameters
                cluster_sim_match = re.search(r'Initial clustering: ([\\d.]+)', content)
                div_f_match = re.search(r'Rainbow div -f: ([\\d.]+)', content)
                div_K_match = re.search(r'Rainbow div -K: (\\d+)', content)
                merge_r_match = re.search(r'Rainbow merge -r: (\\d+)', content)
                final_cluster_match = re.search(r'Final clustering: ([\\d.]+)', content)
                
                if cluster_sim_match:
                    assembly_info.append(f"**Rainbow Assembly Parameters:**")
                    assembly_info.append(f"- Initial CD-HIT clustering similarity: {cluster_sim_match.group(1)}")
                    if div_f_match:
                        assembly_info.append(f"- Rainbow div -f parameter: {div_f_match.group(1)}")
                    if div_K_match:
                        assembly_info.append(f"- Rainbow div -K parameter: {div_K_match.group(1)}")
                    if merge_r_match:
                        assembly_info.append(f"- Rainbow merge -r parameter: {merge_r_match.group(1)}")
                    if final_cluster_match:
                        assembly_info.append(f"- Final CD-HIT clustering similarity: {final_cluster_match.group(1)}")
                    assembly_info.append("")
                
                # Extract assembly metrics
                input_seqs_match = re.search(r'Input sequences: (\\d+)', content)
                final_contigs_match = re.search(r'Final reference contigs: (\\d+)', content)
                total_bases_match = re.search(r'Total bases: (\\d+)', content)
                n50_match = re.search(r'N50: (\\d+)', content)
                min_contig_match = re.search(r'Min contig: (\\d+)', content)
                max_contig_match = re.search(r'Max contig: (\\d+)', content)
                mean_contig_match = re.search(r'Mean contig: (\\d+)', content)
                
                if final_contigs_match:
                    assembly_info.append(f"**Assembly Metrics:**")
                    if input_seqs_match:
                        assembly_info.append(f"- Input sequences to Rainbow: {int(input_seqs_match.group(1)):,}")
                    if final_contigs_match:
                        assembly_info.append(f"- Final reference contigs: {int(final_contigs_match.group(1)):,}")
                    if total_bases_match:
                        total_bases = int(total_bases_match.group(1))
                        if total_bases >= 1e9:
                            assembly_info.append(f"- Total assembly size: {total_bases/1e9:.2f} Gbp")
                        elif total_bases >= 1e6:
                            assembly_info.append(f"- Total assembly size: {total_bases/1e6:.2f} Mbp")
                        else:
                            assembly_info.append(f"- Total assembly size: {total_bases/1e3:.2f} Kbp")
                    if n50_match:
                        assembly_info.append(f"- N50: {int(n50_match.group(1)):,} bp")
                    if mean_contig_match:
                        assembly_info.append(f"- Mean contig length: {int(mean_contig_match.group(1)):,} bp")
                    if min_contig_match and max_contig_match:
                        assembly_info.append(f"- Contig length range: {int(min_contig_match.group(1)):,} - {int(max_contig_match.group(1)):,} bp")
                    assembly_info.append("")
                
                # Extract contig size distribution
                dist_10kb = re.search(r'>10kb: (\\d+)', content)
                dist_5_10kb = re.search(r'5-10kb: (\\d+)', content)
                dist_1_5kb = re.search(r'1-5kb: (\\d+)', content)
                dist_500_1kb = re.search(r'500bp-1kb: (\\d+)', content)
                dist_sub500 = re.search(r'<500bp: (\\d+)', content)
                
                if any([dist_10kb, dist_5_10kb, dist_1_5kb, dist_500_1kb, dist_sub500]):
                    assembly_info.append(f"**Contig Size Distribution:**")
                    if dist_10kb:
                        assembly_info.append(f"- >10kb: {int(dist_10kb.group(1)):,} contigs")
                    if dist_5_10kb:
                        assembly_info.append(f"- 5-10kb: {int(dist_5_10kb.group(1)):,} contigs")
                    if dist_1_5kb:
                        assembly_info.append(f"- 1-5kb: {int(dist_1_5kb.group(1)):,} contigs")
                    if dist_500_1kb:
                        assembly_info.append(f"- 500bp-1kb: {int(dist_500_1kb.group(1)):,} contigs")
                    if dist_sub500:
                        assembly_info.append(f"- <500bp: {int(dist_sub500.group(1)):,} contigs")
                    
        except Exception as e:
            print(f"Could not read assembly stats: {e}")
    
    if assembly_info:
        assembly_section = "\\n".join(assembly_info)
        reference_line = f"## De Novo Assembly\\n\\n{assembly_section}"
    else:
        reference_line = "## De Novo Assembly\\n\\nReference genome assembled de novo from cleaned reads using Rainbow pipeline optimized for ddRAD data."
    
elif genome_source.startswith("accession:"):
    accession = genome_source.replace("accession:", "")
    
    # Try to fetch species name from NCBI using datasets CLI if available
    try:
        result = subprocess.run(
            ['datasets', 'summary', 'genome', 'accession', accession, '--as-json-lines'],
            capture_output=True, text=True, timeout=30
        )
        
        if result.returncode == 0:
            import json
            data = json.loads(result.stdout)
            if 'reports' in data and len(data['reports']) > 0:
                report = data['reports'][0]
                if 'organism' in report:
                    if 'organism_name' in report['organism']:
                        species_name = report['organism']['organism_name']
                    elif 'sci_name' in report['organism']:
                        species_name = report['organism']['sci_name']
            print(f"Found species: {species_name}")
        else:
            print(f"Could not fetch species info: {result.stderr}")
            
    except Exception as e:
        print(f"Could not fetch species name from NCBI: {e}")
        species_name = f"Species for {accession}"
    
    reference_line = f"Reference genome used: [{accession}](https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/)"
    
elif genome_source.startswith("local:"):
    genome_path = genome_source.replace("local:", "")
    species_name = ""
    reference_line = f"Reference genome used: Local file - `{genome_path}`"
elif genome_source.startswith("none:"):
    species_name = ""
    reference_line = "No reference genome used - cleaned reads output only"
else:
    species_name = "Unknown"
    reference_line = "Reference genome used: Unknown source"

# ----------------------------------------------------------------
# Read the read count summary to get statistics
# ----------------------------------------------------------------
initial_stats = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0, "total": 0}
cleaned_stats = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0, "total": 0}
final_stats   = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0, "map_pct": 0, "total": 0}
pp_stats      = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0}

# Variables for calculating total base pairs
total_bases_raw = 0
total_bases_cleaned = 0
total_bases_mapped = 0
read_length_estimate = 150  # Default estimate

try:
    with open("${read_summary}", 'r') as f:
        lines = f.readlines()
        
        if len(lines) > 0:
            header = lines[0].strip().split('\\t')
            
            col_indices = {}
            for i, col in enumerate(header):
                col_lower = col.lower()
                if col_lower in ['raw', 'map', 'pp', 'repr', 'sample_id']:
                    col_indices[col_lower] = i
            
            print(f"Found columns: {col_indices}")
            
            raw_reads = []
            cleaned_reads_list = []
            final_reads = []
            properly_paired = []
            
            for line in lines[1:]:
                parts = line.strip().split('\\t')
                
                # Get raw reads
                if 'raw' in col_indices and col_indices['raw'] < len(parts):
                    try:
                        val = parts[col_indices['raw']]
                        if val and val != 'NA' and val != '':
                            raw_val = int(float(val))
                            raw_reads.append(raw_val)
                            initial_stats["total"] += raw_val
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing raw value: {e}")
                
                # Get cleaned reads (repr = after repair, the final QC step)
                if 'repr' in col_indices and col_indices['repr'] < len(parts):
                    try:
                        val = parts[col_indices['repr']]
                        if val and val != 'NA' and val != '':
                            repr_val = int(float(val))
                            cleaned_reads_list.append(repr_val)
                            cleaned_stats["total"] += repr_val
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing repr value: {e}")
                
                # Get mapped reads
                if 'map' in col_indices and col_indices['map'] < len(parts):
                    try:
                        val = parts[col_indices['map']]
                        if val and val != 'NA' and val != '':
                            map_val = int(float(val))
                            final_reads.append(map_val)
                            final_stats["total"] += map_val
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing map value: {e}")
                
                # Get properly paired if available
                if 'pp' in col_indices and col_indices['pp'] < len(parts):
                    try:
                        val = parts[col_indices['pp']]
                        if val and val != 'NA' and val != '':
                            properly_paired.append(int(float(val)))
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing pp value: {e}")
            
            print(f"Raw reads: {len(raw_reads)} samples, total: {initial_stats['total']}")
            print(f"Cleaned reads: {len(cleaned_reads_list)} samples, total: {cleaned_stats['total']}")
            print(f"Mapped reads: {len(final_reads)} samples, total: {final_stats['total']}")
            print(f"Properly paired: {len(properly_paired)} samples")
            
            if raw_reads:
                initial_stats["mean"] = statistics.mean(raw_reads)
                initial_stats["sd"] = statistics.stdev(raw_reads) if len(raw_reads) > 1 else 0
                initial_stats["min"] = min(raw_reads)
                initial_stats["max"] = max(raw_reads)
                initial_stats["n"] = len(raw_reads)
            
            if cleaned_reads_list:
                cleaned_stats["mean"] = statistics.mean(cleaned_reads_list)
                cleaned_stats["sd"] = statistics.stdev(cleaned_reads_list) if len(cleaned_reads_list) > 1 else 0
                cleaned_stats["min"] = min(cleaned_reads_list)
                cleaned_stats["max"] = max(cleaned_reads_list)
                cleaned_stats["n"] = len(cleaned_reads_list)
                
            if final_reads:
                final_stats["mean"] = statistics.mean(final_reads)
                final_stats["sd"] = statistics.stdev(final_reads) if len(final_reads) > 1 else 0
                final_stats["min"] = min(final_reads)
                final_stats["max"] = max(final_reads)
                final_stats["n"] = len(final_reads)
                try:
                    with open("${mapping_summary}", 'r') as f:
                        map_lines = f.readlines()
                        if len(map_lines) > 1:
                            total_reads_all = 0
                            mapped_reads_all = 0
                            for line in map_lines[1:]:
                                parts = line.strip().split('\\t')
                                if len(parts) >= 5:
                                    try:
                                        total_reads_all += int(parts[1])
                                        mapped_reads_all += int(parts[2])
                                    except (ValueError, IndexError):
                                        continue
                            
                            if total_reads_all > 0:
                                overall_map_rate = (mapped_reads_all / total_reads_all) * 100
                                final_stats["map_pct"] = overall_map_rate
                            else:
                                final_stats["map_pct"] = 0
                except Exception as e:
                    print(f"Could not calculate mapping rate from mapping_summary: {e}")
                    final_stats["map_pct"] = 0
            
            if properly_paired:
                pp_stats["mean"] = statistics.mean(properly_paired)
                pp_stats["sd"] = statistics.stdev(properly_paired) if len(properly_paired) > 1 else 0
                pp_stats["min"] = min(properly_paired)
                pp_stats["max"] = max(properly_paired)
                pp_stats["n"] = len(properly_paired)
                
except Exception as e:
    print(f"Could not read statistics from read summary: {e}")

# Try to get actual read length from MultiQC general stats
try:
    for mqc_file in os.listdir('.'):
        if mqc_file.startswith('multiqc_raw_fastqc') and mqc_file.endswith('_general_stats.txt'):
            with open(mqc_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:
                    header = lines[0].strip().split('\\t')
                    for i, col in enumerate(header):
                        if 'sequence_length' in col.lower() or 'avg_sequence_length' in col.lower():
                            lengths = []
                            for line in lines[1:]:
                                parts = line.strip().split('\\t')
                                if i < len(parts):
                                    try:
                                        if '-' in parts[i]:
                                            range_parts = parts[i].split('-')
                                            avg_len = (int(range_parts[0]) + int(range_parts[1])) / 2
                                            lengths.append(avg_len)
                                        else:
                                            lengths.append(float(parts[i]))
                                    except:
                                        continue
                            if lengths:
                                read_length_estimate = int(statistics.mean(lengths))
                                print(f"Found actual read length from FastQC: {read_length_estimate}")
                            break
                    break
except Exception as e:
    print(f"Could not extract read length from FastQC data, using default {read_length_estimate}: {e}")

# Calculate total base pairs
total_bases_raw = initial_stats["total"] * read_length_estimate * 2
total_bases_cleaned = cleaned_stats["total"] * read_length_estimate * 2
total_bases_mapped = final_stats["total"] * read_length_estimate * 2

def format_bases(n_bases):
    if n_bases >= 1e9:
        return f"{n_bases/1e9:.2f} Gbp"
    elif n_bases >= 1e6:
        return f"{n_bases/1e6:.2f} Mbp"
    elif n_bases >= 1e3:
        return f"{n_bases/1e3:.2f} Kbp"
    else:
        return f"{n_bases:.0f} bp"

def fmt_num(n):
    return f"{n:,.0f}"

# Find MultiQC reports and sort them
multiqc_files = [f for f in os.listdir('.') if f.startswith('multiqc_') and f.endswith('.html')]

stage_order = ['raw_fastqc', 'fastp_trim_3', 'clumpify', 'fastp_trim_5', 'fastq_screen', 'repair', 'mapping']

sorted_multiqc = []
for stage in stage_order:
    for f in multiqc_files:
        if stage in f:
            sorted_multiqc.append(f)
            break

stage_names = {
    'raw_fastqc': 'Raw FastQC',
    'fastp_trim_3': "3' Trimming (Fastp)",
    'clumpify': 'Deduplication (Clumpify)',
    'fastp_trim_5': "5' Trimming (Fastp)",
    'fastq_screen': 'Contamination Screening',
    'repair': 'Read Repair',
    'mapping': 'Mapping to Reference',
    'mapping_denovo': 'Mapping to De Novo Assembly'
}

multiqc_links = []
for mqc_file in sorted_multiqc:
    for stage, name in stage_names.items():
        if stage in mqc_file:
            multiqc_links.append(f"- [{name}](${params.outdir}/qc/multiqc_reports/{mqc_file})")
            break

# ----------------------------------------------------------------
# Process species identification results if available
# ----------------------------------------------------------------
species_id_section = ""
if species_id_performed:
    try:
        import csv
        from collections import Counter
        
        with open("${species_top_hits}", 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        
        if rows:
            species_counts = {}
            taxonomic_levels = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
            
            for col in taxonomic_levels:
                if col in rows[0]:
                    values = [row[col] for row in rows if row.get(col) and row[col].strip()]
                    if values:
                        counter = Counter(values)
                        most_common = counter.most_common(1)
                        if most_common:
                            species_counts[col] = most_common[0][0]
            
            species_id_section = "## Species Identification (BLAST Analysis)\\n\\n"
            
            if species_counts:
                species_id_section += "**Taxonomic Classification (Most Likely):**\\n"
                for level, taxon in species_counts.items():
                    if taxon and taxon != '':
                        species_id_section += f"- {level.capitalize()}: {taxon}\\n"
                species_id_section += "\\n"
            
            n_samples = len(rows)
            species_id_section += f"**Sample Summary:**\\n"
            species_id_section += f"- Number of samples analyzed: {n_samples}\\n"
            
            if 'species' in rows[0]:
                species_values = [row.get('species', '') for row in rows if row.get('species', '').strip()]
                if species_values:
                    species_counter = Counter(species_values)
                    most_common_species = species_counter.most_common(1)
                    if most_common_species:
                        top_species = most_common_species[0][0]
                        n_agree = most_common_species[0][1]
                        species_id_section += f"- Samples with consensus species: {n_agree}/{n_samples} ({n_agree/n_samples*100:.1f}%)\\n"
            
            species_id_section += "\\n"
            
            species_id_section += "### BLAST Hit Distributions\\n\\n"
            species_id_section += "#### Raw BLAST Results by Taxonomic Level\\n"
            species_id_section += "![Raw BLAST Pie Charts](${params.outdir}/species_id/blast_raw_pie.png)\\n\\n"
            
            species_id_section += "#### Summarized Species Identification\\n"
            species_id_section += "![Summary BLAST Pie Charts](${params.outdir}/species_id/blast_summary_pie.png)\\n\\n"
            
            species_id_section += "### Detailed Results\\n\\n"
            species_id_section += "- [Combined BLAST results (TSV)](${params.outdir}/species_id/blast_results.tsv)\\n"
            species_id_section += "- [Top BLAST hits per sample (CSV)](${params.outdir}/species_id/top_blast_hits.csv)\\n"
            species_id_section += "- [Posterior probabilities by taxonomic level](${params.outdir}/species_id/blast_posteriors/)\\n"
            
    except Exception as e:
        print(f"Could not process species identification results: {e}")
        species_id_section = "## Species Identification\\n\\nSpecies identification was run but results could not be processed.\\n"
else:
    print("Species identification was not performed or no results available")


# ----------------------------------------------------------------
# Build the cleaned reads section (always present)
# ----------------------------------------------------------------
cleaned_retention_pct = ""
if initial_stats["total"] > 0 and cleaned_stats["total"] > 0:
    cleaned_retention_pct = f"- Retention from raw: {(cleaned_stats['total']/initial_stats['total']*100):.1f}% of initial read pairs"

cleaned_section = (
    "## Final Cleaned Reads\\n"
    "\\n"
    f"**Summary Statistics (n={cleaned_stats['n']} samples):**\\n"
    f"- Mean reads per sample: {fmt_num(cleaned_stats['mean'])}\\n"
    f"- Standard deviation: {fmt_num(cleaned_stats['sd'])}\\n"
    f"- Min reads: {fmt_num(cleaned_stats['min'])}\\n"
    f"- Max reads: {fmt_num(cleaned_stats['max'])}\\n"
    "\\n"
    "**Total Cleaned Output:**\\n"
    f"- Total cleaned read pairs: {fmt_num(cleaned_stats['total'])}\\n"
    f"- Total cleaned bases: {format_bases(total_bases_cleaned)} (assuming {read_length_estimate}bp reads)\\n"
    f"- Total individual reads: {fmt_num(cleaned_stats['total'] * 2)} (paired-end)\\n"
    f"{cleaned_retention_pct}"
)


# ----------------------------------------------------------------
# Build the mapping section (only if mapping was performed)
# ----------------------------------------------------------------
mapping_section = ""
if mapping_performed and final_stats["n"] > 0:
    mapped_retention_pct = ""
    if initial_stats["total"] > 0 and final_stats["total"] > 0:
        mapped_retention_pct = f"- Retention from raw: {(final_stats['total']/initial_stats['total']*100):.1f}% of initial read pairs"

    insert_size_line = ""
    if insert_size_violin_performed:
        insert_size_line = "\\n![Insert Size Distribution](${params.outdir}/qc/insert_size_violin.png)\\n"

    mapping_section = (
        "## Mapped Reads\\n"
        "\\n"
        "![Mapped Read Distribution](${params.outdir}/qc/mapped_reads_histogram.png)\\n"
        f"{insert_size_line}"
        "\\n"
        f"**Summary Statistics (n={final_stats['n']} samples):**\\n"
        f"- Mean reads per sample: {fmt_num(final_stats['mean'])}\\n"
        f"- Standard deviation: {fmt_num(final_stats['sd'])}\\n"
        f"- Min reads: {fmt_num(final_stats['min'])}\\n"
        f"- Max reads: {fmt_num(final_stats['max'])}\\n"
        f"- Overall mapping rate: {final_stats['map_pct']:.2f}%\\n"
        "\\n"
        "**Total Mapped Output:**\\n"
        f"- Total read pairs mapped: {fmt_num(final_stats['total'])}\\n"
        f"- Total bases mapped: {format_bases(total_bases_mapped)} (assuming {read_length_estimate}bp reads)\\n"
        f"{mapped_retention_pct}\\n"
        "\\n"
        "### Final Mapping Statistics\\n"
        "See [${mapping_summary}](${params.outdir}/qc/${mapping_summary}) for detailed mapping statistics per sample."
    )


# ----------------------------------------------------------------
# Generate the markdown report
# ----------------------------------------------------------------
markdown_content = f\'\'\'# GCL Illumina QC Pipeline Report

{"## " + species_name if species_name else ""}

{reference_line}

## Initial Sequencing

![Initial Read Distribution](${params.outdir}/qc/initial_reads_histogram.png)

**Summary Statistics (n={initial_stats["n"]} samples):**
- Mean reads per sample: {fmt_num(initial_stats["mean"])}
- Standard deviation: {fmt_num(initial_stats["sd"])}
- Min reads: {fmt_num(initial_stats["min"])}
- Max reads: {fmt_num(initial_stats["max"])}

**Total Sequencing Output:**
- Total read pairs sequenced: {fmt_num(initial_stats["total"])}
- Total bases sequenced: {format_bases(total_bases_raw)} (assuming {read_length_estimate}bp reads)
- Total individual reads: {fmt_num(initial_stats["total"] * 2)} (paired-end)

## Sequencing QC

### Read Retention Through QC Pipeline
![QC Summary Plot](${params.outdir}/qc/qc_summary_plot.png)

### MultiQC Reports (download locally to view)
{chr(10).join(multiqc_links) if multiqc_links else "No MultiQC reports found"}

### Stage-by-Stage Comparison
See [stage_comparison.txt](${params.outdir}/qc/stage_comparison.txt) for detailed retention rates between stages.

{cleaned_section}

{mapping_section}

{species_id_section}

---
*QC Report generated on: {subprocess.check_output(['date']).decode().strip()}*
\'\'\'

# Write markdown file
with open("qc_pipeline_report.md", 'w') as f:
    f.write(markdown_content)

print("Markdown report generated successfully!")
print(f"\\nReport Summary:")
print(f"  Total read pairs sequenced: {fmt_num(initial_stats['total'])}")
print(f"  Total bases sequenced: {format_bases(total_bases_raw)}")
print(f"  Total cleaned read pairs: {fmt_num(cleaned_stats['total'])}")
print(f"  Total cleaned bases: {format_bases(total_bases_cleaned)}")
if mapping_performed:
    print(f"  Total read pairs mapped: {fmt_num(final_stats['total'])}")
    print(f"  Total bases mapped: {format_bases(total_bases_mapped)}")
else:
    print(f"  Mapping: Not performed")
if species_id_performed:
    print(f"  Species identification: Completed")

# Convert to HTML using pandoc if available
try:
    subprocess.run(['pandoc', 'qc_pipeline_report.md', '-o', 'qc_pipeline_report.html', '--standalone', '--metadata', 'title=GCL QC Pipeline Report'], check=True)
    print("HTML report generated successfully!")
except:
    html_content = f\'\'\'<!DOCTYPE html>
<html>
<head>
    <title>GCL QC Pipeline Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; max-width: 1200px; margin: 0 auto; padding: 20px; }}
        h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
        h2 {{ color: #34495e; margin-top: 30px; }}
        h3 {{ color: #7f8c8d; }}
        img {{ max-width: 100%; height: auto; border: 1px solid #ddd; padding: 5px; margin: 10px 0; }}
        code {{ background-color: #f4f4f4; padding: 2px 5px; border-radius: 3px; }}
        pre {{ background-color: #f4f4f4; padding: 10px; border-radius: 5px; overflow-x: auto; }}
    </style>
</head>
<body>
<pre>{markdown_content}</pre>
</body>
</html>\'\'\'
    with open("qc_pipeline_report.html", 'w') as f:
        f.write(html_content)
    print("Basic HTML report generated (pandoc not available)")
PYEOF

    python3 generate_report.py
    """
}