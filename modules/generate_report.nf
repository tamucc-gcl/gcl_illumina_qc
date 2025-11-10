process generate_report {
    label 'generate_report'
    tag "final_report"
    
    publishDir "${params.outdir}", mode: 'copy'
    
    input:
        path qc_plot
        path read_summary
        path multiqc_reports
        val genome_source
        path initial_histogram
        path mapped_histogram
        path mapping_summary
        
    output:
        path "qc_pipeline_report.md"
        path "qc_pipeline_report.html"
    
    script:
    """
    # Check if mapping was performed
    if [ "${mapping_summary}" == "NO_MAPPING" ]; then
        echo "No mapping performed" > mapping_summary.txt
        MAPPING_PERFORMED="false"
    else
        MAPPING_PERFORMED="true"
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

# Initialize variables
species_name = ""
reference_line = ""

if genome_source.startswith("accession:"):
    accession = genome_source.replace("accession:", "")
    
    # Try to fetch species name from NCBI using datasets CLI if available
    try:
        # Try using datasets CLI to get species info
        result = subprocess.run(
            ['datasets', 'summary', 'genome', 'accession', accession, '--as-json-lines'],
            capture_output=True, text=True, timeout=30
        )
        
        if result.returncode == 0:
            import json
            data = json.loads(result.stdout)
            # Navigate the JSON structure to find species name
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
    
    # Set reference line with NCBI link
    reference_line = f"Reference genome used: [{accession}](https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/)"
    
elif genome_source.startswith("local:"):
    genome_path = genome_source.replace("local:", "")
    species_name = ""  # Leave blank for local genomes
    reference_line = f"Reference genome used: Local file - `{genome_path}`"
else:
    species_name = "Unknown"
    reference_line = "Reference genome used: Unknown source"

# Read the read count summary to get statistics
initial_stats = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0, "total": 0}
final_stats = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0, "map_pct": 0, "total": 0}
pp_stats = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0}

# Variables for calculating total base pairs
total_bases_raw = 0
total_bases_mapped = 0
read_length_estimate = 150  # Default estimate, will try to get actual from FastQC

try:
    with open("${read_summary}", 'r') as f:
        lines = f.readlines()
        
        if len(lines) > 0:
            # Read header to get stage names (wide format)
            header = lines[0].strip().split('\\t')
            
            # Find column indices for stages we care about
            col_indices = {}
            for i, col in enumerate(header):
                col_lower = col.lower()
                if col_lower in ['raw', 'map', 'pp', 'repr', 'sample_id']:
                    col_indices[col_lower] = i
            
            print(f"Found columns: {col_indices}")
            
            # Process each sample (each row after header)
            raw_reads = []
            final_reads = []
            properly_paired = []
            
            for line in lines[1:]:  # Skip header
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
                
                # Get mapped reads as final
                if 'map' in col_indices and col_indices['map'] < len(parts):
                    try:
                        val = parts[col_indices['map']]
                        if val and val != 'NA' and val != '':
                            map_val = int(float(val))
                            final_reads.append(map_val)
                            final_stats["total"] += map_val
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing map value: {e}")
                elif 'repr' in col_indices and col_indices['repr'] < len(parts):
                    try:
                        val = parts[col_indices['repr']]
                        if val and val != 'NA' and val != '':
                            repr_val = int(float(val))
                            final_reads.append(repr_val)
                            final_stats["total"] += repr_val
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing repr value: {e}")
                
                # Get properly paired if available
                if 'pp' in col_indices and col_indices['pp'] < len(parts):
                    try:
                        val = parts[col_indices['pp']]
                        if val and val != 'NA' and val != '':
                            properly_paired.append(int(float(val)))
                    except (ValueError, IndexError) as e:
                        print(f"Error parsing pp value: {e}")
            
            print(f"Raw reads: {len(raw_reads)} samples, total: {initial_stats['total']}")
            print(f"Final reads: {len(final_reads)} samples, total: {final_stats['total']}")
            print(f"Properly paired: {len(properly_paired)} samples")
            
            if raw_reads:
                initial_stats["mean"] = statistics.mean(raw_reads)
                initial_stats["sd"] = statistics.stdev(raw_reads) if len(raw_reads) > 1 else 0
                initial_stats["min"] = min(raw_reads)
                initial_stats["max"] = max(raw_reads)
                initial_stats["n"] = len(raw_reads)
                
            if final_reads:
                final_stats["mean"] = statistics.mean(final_reads)
                final_stats["sd"] = statistics.stdev(final_reads) if len(final_reads) > 1 else 0
                final_stats["min"] = min(final_reads)
                final_stats["max"] = max(final_reads)
                final_stats["n"] = len(final_reads)
                try:
                    with open("${mapping_summary}", 'r') as f:
                        map_lines = f.readlines()
                        if len(map_lines) > 1:  # Skip header
                            total_reads_all = 0
                            mapped_reads_all = 0
                            for line in map_lines[1:]:
                                parts = line.strip().split('\\t')
                                if len(parts) >= 5:
                                    try:
                                        total_reads_all += int(parts[1])  # Total_Reads column
                                        mapped_reads_all += int(parts[2])  # Mapped_Reads column
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
            
            # Calculate properly paired statistics if available
            if properly_paired:
                pp_stats["mean"] = statistics.mean(properly_paired)
                pp_stats["sd"] = statistics.stdev(properly_paired) if len(properly_paired) > 1 else 0
                pp_stats["min"] = min(properly_paired)
                pp_stats["max"] = max(properly_paired)
                pp_stats["n"] = len(properly_paired)
                
except Exception as e:
    print(f"Could not read statistics from read summary: {e}")

# Try to get actual read length from MultiQC general stats (FastQC data)
# Look for the raw_fastqc MultiQC report's general stats
try:
    for mqc_file in os.listdir('.'):
        if mqc_file.startswith('multiqc_raw_fastqc') and mqc_file.endswith('_general_stats.txt'):
            with open(mqc_file, 'r') as f:
                lines = f.readlines()
                if len(lines) > 1:  # Has header and data
                    header = lines[0].strip().split('\\t')
                    # Look for sequence length column (FastQC reports this)
                    for i, col in enumerate(header):
                        if 'sequence_length' in col.lower() or 'avg_sequence_length' in col.lower():
                            # Get average from all samples
                            lengths = []
                            for line in lines[1:]:
                                parts = line.strip().split('\\t')
                                if i < len(parts):
                                    try:
                                        # Handle ranges like "150-151" by taking the average
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
# Note: We multiply by 2 for paired-end reads (R1 + R2)
total_bases_raw = initial_stats["total"] * read_length_estimate * 2
total_bases_mapped = final_stats["total"] * read_length_estimate * 2

# Convert to appropriate units (Gbp, Mbp, etc.)
def format_bases(n_bases):
    if n_bases >= 1e9:
        return f"{n_bases/1e9:.2f} Gbp"
    elif n_bases >= 1e6:
        return f"{n_bases/1e6:.2f} Mbp"
    elif n_bases >= 1e3:
        return f"{n_bases/1e3:.2f} Kbp"
    else:
        return f"{n_bases:.0f} bp"

# Find MultiQC reports and sort them
multiqc_files = [f for f in os.listdir('.') if f.startswith('multiqc_') and f.endswith('.html')]

# Define the expected order
stage_order = ['raw_fastqc', 'fastp_trim_3', 'clumpify', 'fastp_trim_5', 'fastq_screen', 'repair', 'mapping']

# Sort MultiQC files by stage order
sorted_multiqc = []
for stage in stage_order:
    for f in multiqc_files:
        if stage in f:
            sorted_multiqc.append(f)
            break

# Create MultiQC links section
multiqc_links = []
stage_names = {
    'raw_fastqc': 'Raw FastQC',
    'fastp_trim_3': "3' Trimming (Fastp)",
    'clumpify': 'Deduplication (Clumpify)',
    'fastp_trim_5': "5' Trimming (Fastp)",
    'fastq_screen': 'Contamination Screening',
    'repair': 'Read Repair',
    'mapping': 'Mapping to Reference'
}

for mqc_file in sorted_multiqc:
    # Extract stage from filename
    for stage, name in stage_names.items():
        if stage in mqc_file:
            multiqc_links.append(f"- [{name}](${params.outdir}/multiqc_reports/{mqc_file})")
            break

# Format numbers with thousands separator
def fmt_num(n):
    return f"{n:,.0f}"

# Generate properly paired section
pp_section = ""
if pp_stats["n"] > 0:
    pp_section = "### Properly Paired Reads\\n"
    pp_section += f"**Summary Statistics (n={pp_stats['n']} samples):**\\n"
    pp_section += f"- Mean reads: {fmt_num(pp_stats['mean'])}\\n"
    pp_section += f"- Standard deviation: {fmt_num(pp_stats['sd'])}\\n"
    pp_section += f"- Min reads: {fmt_num(pp_stats['min'])}\\n"
    pp_section += f"- Max reads: {fmt_num(pp_stats['max'])}"

# Generate the markdown report
markdown_content = f'''# GCL Illumina QC Pipeline Report

{"## " + species_name if species_name else ""}

{reference_line}

## Initial Sequencing

![Initial Read Distribution](${params.outdir}/qc_analysis/initial_reads_histogram.png)

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
![QC Summary Plot](${params.outdir}/qc_analysis/qc_summary_plot.png)

### MultiQC Reports (download locally to view)
{chr(10).join(multiqc_links) if multiqc_links else "No MultiQC reports found"}

### Stage-by-Stage Comparison
See [stage_comparison.txt](${params.outdir}/qc_analysis/stage_comparison.txt) for detailed retention rates between stages.

## Post QC

### Mapped Reads

![Mapped Read Distribution](${params.outdir}/qc_analysis/mapped_reads_histogram.png)

**Summary Statistics (n={final_stats["n"]} samples):**
- Mean reads per sample: {fmt_num(final_stats["mean"])}
- Standard deviation: {fmt_num(final_stats["sd"])}
- Min reads: {fmt_num(final_stats["min"])}
- Max reads: {fmt_num(final_stats["max"])}
- Overall mapping rate: {final_stats["map_pct"]:.2f}%

**Total Mapped Output:**
- Total read pairs mapped: {fmt_num(final_stats["total"])}
- Total bases mapped: {format_bases(total_bases_mapped)} (assuming {read_length_estimate}bp reads)
- Retention from raw: {(final_stats["total"]/initial_stats["total"]*100):.1f}% of initial read pairs

### Final Mapping Statistics
See [${mapping_summary}](${params.outdir}/qc_analysis/${mapping_summary}) for detailed mapping statistics per sample.

---
*QC Report generated on: {subprocess.check_output(['date']).decode().strip()}*
'''

# Write markdown file
with open("qc_pipeline_report.md", 'w') as f:
    f.write(markdown_content)

print("Markdown report generated successfully!")
print(f"\\nReport Summary:")
print(f"  Total read pairs sequenced: {fmt_num(initial_stats['total'])}")
print(f"  Total bases sequenced: {format_bases(total_bases_raw)}")
print(f"  Total read pairs mapped: {fmt_num(final_stats['total'])}")
print(f"  Total bases mapped: {format_bases(total_bases_mapped)}")

# Convert to HTML using pandoc if available
try:
    subprocess.run(['pandoc', 'qc_pipeline_report.md', '-o', 'qc_pipeline_report.html', '--standalone', '--metadata', 'title=GCL QC Pipeline Report'], check=True)
    print("HTML report generated successfully!")
except:
    # If pandoc not available, create a simple HTML wrapper
    html_content = f'''<!DOCTYPE html>
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
</html>'''
    with open("qc_pipeline_report.html", 'w') as f:
        f.write(html_content)
    print("Basic HTML report generated (pandoc not available)")
PYEOF

    python3 generate_report.py
    """
}