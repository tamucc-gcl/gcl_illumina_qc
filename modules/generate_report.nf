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
        path soft_clip_violin      // NEW
        path aln_score_violin      // NEW
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

    # Check if soft clipping / alignment score plots exist
    if [ "\$MAPPING_PERFORMED" == "true" ] && [ -f "${soft_clip_violin}" ] && [[ "${soft_clip_violin}" != *"no_"* ]]; then
        export SOFT_CLIP_PERFORMED="true"
    else
        export SOFT_CLIP_PERFORMED="false"
    fi

    if [ "\$MAPPING_PERFORMED" == "true" ] && [ -f "${aln_score_violin}" ] && [[ "${aln_score_violin}" != *"no_"* ]]; then
        export ALN_SCORE_PERFORMED="true"
    else
        export ALN_SCORE_PERFORMED="false"
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
mapping_performed          = os.environ.get("MAPPING_PERFORMED", "false") == "true"
insert_size_violin_performed = os.environ.get("INSERT_SIZE_VIOLIN_PERFORMED", "false") == "true"
soft_clip_performed        = os.environ.get("SOFT_CLIP_PERFORMED", "false") == "true"
aln_score_performed        = os.environ.get("ALN_SCORE_PERFORMED", "false") == "true"
assembly_performed         = os.environ.get("ASSEMBLY_PERFORMED", "false") == "true"
filter_performed           = os.environ.get("FILTER_PERFORMED", "false") == "true"
species_id_performed       = os.environ.get("SPECIES_ID_PERFORMED", "false") == "true"

# ----------------------------------------------------------------
# Genome / reference section  (unchanged from your current version)
# ----------------------------------------------------------------
species_name = ""
reference_line = ""
assembly_section = ""

if genome_source.startswith("denovo:"):
    species_name = ""
    assembly_info = []

    if filter_performed:
        try:
            with open("${filter_stats}", 'r') as f:
                filter_content = f.read()
            cutoff1_match   = re.search(r'Cutoff 1.*?: (\\d+)', filter_content)
            cutoff2_match   = re.search(r'Cutoff 2.*?: (\\d+)', filter_content)
            total_uniq_match   = re.search(r'Total unique sequences: (\\d+)', filter_content)
            after_indiv_match  = re.search(r'Sequences after per-individual filter: (\\d+)', filter_content)
            after_count_match  = re.search(r'Sequences after individual count filter: (\\d+)', filter_content)
            final_seqs_match   = re.search(r'Final sequences for assembly: (\\d+)', filter_content)
            if cutoff1_match and cutoff2_match:
                assembly_info.append("**Filtering Parameters:**")
                assembly_info.append(f"- Minimum reads per individual (cutoff1): {cutoff1_match.group(1)}")
                assembly_info.append(f"- Minimum individuals required (cutoff2): {cutoff2_match.group(1)}")
                assembly_info.append("")
            if total_uniq_match:
                assembly_info.append("**Filtering Statistics:**")
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
            with open("${assembly_stats}", 'r') as f:
                content = f.read()
            cluster_sim_match  = re.search(r'Initial clustering: ([\\d.]+)', content)
            div_f_match        = re.search(r'Rainbow div -f: ([\\d.]+)', content)
            div_K_match        = re.search(r'Rainbow div -K: (\\d+)', content)
            merge_r_match      = re.search(r'Rainbow merge -r: (\\d+)', content)
            final_cluster_match= re.search(r'Final clustering: ([\\d.]+)', content)
            if cluster_sim_match:
                assembly_info.append("**Rainbow Assembly Parameters:**")
                assembly_info.append(f"- Initial CD-HIT clustering similarity: {cluster_sim_match.group(1)}")
                if div_f_match:   assembly_info.append(f"- Rainbow div -f parameter: {div_f_match.group(1)}")
                if div_K_match:   assembly_info.append(f"- Rainbow div -K parameter: {div_K_match.group(1)}")
                if merge_r_match: assembly_info.append(f"- Rainbow merge -r parameter: {merge_r_match.group(1)}")
                if final_cluster_match: assembly_info.append(f"- Final CD-HIT clustering similarity: {final_cluster_match.group(1)}")
                assembly_info.append("")
            input_seqs_match    = re.search(r'Input sequences: (\\d+)', content)
            final_contigs_match = re.search(r'Final reference contigs: (\\d+)', content)
            total_bases_match   = re.search(r'Total bases: (\\d+)', content)
            n50_match           = re.search(r'N50: (\\d+)', content)
            min_contig_match    = re.search(r'Min contig: (\\d+)', content)
            max_contig_match    = re.search(r'Max contig: (\\d+)', content)
            mean_contig_match   = re.search(r'Mean contig: (\\d+)', content)
            if final_contigs_match:
                assembly_info.append("**Assembly Metrics:**")
                if input_seqs_match:    assembly_info.append(f"- Input sequences to Rainbow: {int(input_seqs_match.group(1)):,}")
                assembly_info.append(f"- Final reference contigs: {int(final_contigs_match.group(1)):,}")
                if total_bases_match:
                    tb = int(total_bases_match.group(1))
                    size_str = f"{tb/1e9:.2f} Gbp" if tb >= 1e9 else f"{tb/1e6:.2f} Mbp" if tb >= 1e6 else f"{tb/1e3:.2f} Kbp"
                    assembly_info.append(f"- Total assembly size: {size_str}")
                if n50_match:           assembly_info.append(f"- N50: {int(n50_match.group(1)):,} bp")
                if mean_contig_match:   assembly_info.append(f"- Mean contig length: {int(mean_contig_match.group(1)):,} bp")
                if min_contig_match and max_contig_match:
                    assembly_info.append(f"- Contig length range: {int(min_contig_match.group(1)):,} - {int(max_contig_match.group(1)):,} bp")
                assembly_info.append("")
            for tag, label in [('10kb', '>10kb'), ('5-10kb', '5-10kb'), ('1-5kb', '1-5kb'),
                                ('500bp-1kb', '500bp-1kb'), ('500bp', '<500bp')]:
                m = re.search(rf'{re.escape(label)}: (\\d+)', content)
                if not m:
                    m = re.search(rf'>?{re.escape(tag)}: (\\d+)', content)
            dist_matches = {}
            for pat, lbl in [(r'>10kb: (\\d+)', '>10kb'), (r'5-10kb: (\\d+)', '5-10kb'),
                             (r'1-5kb: (\\d+)', '1-5kb'), (r'500bp-1kb: (\\d+)', '500bp-1kb'),
                             (r'<500bp: (\\d+)', '<500bp')]:
                m = re.search(pat, content)
                if m: dist_matches[lbl] = int(m.group(1))
            if dist_matches:
                assembly_info.append("**Contig Size Distribution:**")
                for lbl, cnt in dist_matches.items():
                    assembly_info.append(f"- {lbl}: {cnt:,} contigs")
        except Exception as e:
            print(f"Could not read assembly stats: {e}")

    if assembly_info:
        assembly_section = "\\n".join(assembly_info)
        reference_line = f"## De Novo Assembly\\n\\n{assembly_section}"
    else:
        reference_line = "## De Novo Assembly\\n\\nReference genome assembled de novo from cleaned reads using Rainbow pipeline optimized for ddRAD data."

elif genome_source.startswith("accession:"):
    accession = genome_source.replace("accession:", "")
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
                    species_name = report['organism'].get('organism_name') or report['organism'].get('sci_name', '')
    except Exception as e:
        print(f"Could not fetch species name from NCBI: {e}")
        species_name = f"Species for {accession}"
    reference_line = f"Reference genome used: [{accession}](https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/)"

elif genome_source.startswith("local:"):
    genome_path = genome_source.replace("local:", "")
    species_name = ""
    reference_line = f"Reference genome used: Local file - '{genome_path}'"
elif genome_source.startswith("none:"):
    species_name = ""
    reference_line = "No reference genome used - cleaned reads output only"
else:
    species_name = "Unknown"
    reference_line = "Reference genome used: Unknown source"

# ----------------------------------------------------------------
# Read count statistics  (unchanged)
# ----------------------------------------------------------------
initial_stats = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0, "total": 0}
cleaned_stats = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0, "total": 0}
final_stats   = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0, "map_pct": 0, "total": 0}
pp_stats      = {"mean": 0, "sd": 0, "min": 0, "max": 0, "n": 0}
total_bases_raw = 0
total_bases_cleaned = 0
total_bases_mapped = 0
read_length_estimate = 150

try:
    with open("${read_summary}", 'r') as f:
        lines = f.readlines()
    if lines:
        header = lines[0].strip().split('\\t')
        col_indices = {col.lower(): i for i, col in enumerate(header) if col.lower() in ['raw', 'map', 'pp', 'repr', 'sample_id']}
        raw_reads, cleaned_reads_list, final_reads, properly_paired = [], [], [], []
        for line in lines[1:]:
            parts = line.strip().split('\\t')
            for key, lst, stat in [('raw', raw_reads, initial_stats), ('repr', cleaned_reads_list, cleaned_stats),
                                   ('map', final_reads, final_stats)]:
                if key in col_indices and col_indices[key] < len(parts):
                    val = parts[col_indices[key]]
                    if val and val not in ('NA', ''):
                        try:
                            v = int(float(val)); lst.append(v); stat['total'] += v
                        except: pass
            if 'pp' in col_indices and col_indices['pp'] < len(parts):
                val = parts[col_indices['pp']]
                if val and val not in ('NA', ''):
                    try: properly_paired.append(int(float(val)))
                    except: pass
        for lst, stat in [(raw_reads, initial_stats), (cleaned_reads_list, cleaned_stats),
                          (final_reads, final_stats), (properly_paired, pp_stats)]:
            if lst:
                stat['mean'] = statistics.mean(lst)
                stat['sd']   = statistics.stdev(lst) if len(lst) > 1 else 0
                stat['min']  = min(lst); stat['max'] = max(lst); stat['n'] = len(lst)
        if final_reads:
            try:
                with open("${mapping_summary}", 'r') as f:
                    mlines = f.readlines()
                tr, mr = 0, 0
                for ml in mlines[1:]:
                    p = ml.strip().split('\\t')
                    if len(p) >= 3:
                        try: tr += int(p[1]); mr += int(p[2])
                        except: pass
                final_stats['map_pct'] = (mr / tr * 100) if tr > 0 else 0
            except: pass
except Exception as e:
    print(f"Could not read statistics: {e}")

try:
    for mqf in os.listdir('.'):
        if mqf.startswith('multiqc_raw_fastqc') and mqf.endswith('_general_stats.txt'):
            with open(mqf) as f:
                hdr = f.readline().strip().split('\\t')
                for i, col in enumerate(hdr):
                    if 'sequence_length' in col.lower():
                        lengths = []
                        for ln in f:
                            p = ln.strip().split('\\t')
                            if i < len(p):
                                try:
                                    s = p[i]
                                    lengths.append((int(s.split('-')[0]) + int(s.split('-')[1]))/2 if '-' in s else float(s))
                                except: pass
                        if lengths: read_length_estimate = int(statistics.mean(lengths))
                        break
            break
except: pass

total_bases_raw     = initial_stats["total"] * read_length_estimate * 2
total_bases_cleaned = cleaned_stats["total"] * read_length_estimate * 2
total_bases_mapped  = final_stats["total"]   * read_length_estimate * 2

def format_bases(n):
    if n >= 1e9: return f"{n/1e9:.2f} Gbp"
    if n >= 1e6: return f"{n/1e6:.2f} Mbp"
    if n >= 1e3: return f"{n/1e3:.2f} Kbp"
    return f"{n:.0f} bp"

def fmt_num(n): return f"{n:,.0f}"

# MultiQC links
multiqc_files = sorted([f for f in os.listdir('.') if f.startswith('multiqc_') and f.endswith('.html')])
stage_order = ['raw_fastqc', 'fastp_trim_3', 'clumpify', 'fastp_trim_5', 'fastq_screen', 'repair', 'mapping']
stage_names = {
    'raw_fastqc':   'Raw FastQC',
    'fastp_trim_3': "3' Trimming (Fastp)",
    'clumpify':     'Deduplication (Clumpify)',
    'fastp_trim_5': "5' Trimming (Fastp)",
    'fastq_screen': 'Contamination Screening',
    'repair':       'Read Repair',
    'mapping':      'Mapping to Reference',
    'mapping_denovo': 'Mapping to De Novo Assembly'
}
sorted_multiqc = []
for stage in stage_order:
    for f in multiqc_files:
        if stage in f:
            sorted_multiqc.append(f); break
multiqc_links = []
for mqf in sorted_multiqc:
    for stage, name in stage_names.items():
        if stage in mqf:
            multiqc_links.append(f"- [{name}](${params.outdir}/qc/multiqc_reports/{mqf})")
            break

# ----------------------------------------------------------------
# Species ID section  (unchanged)
# ----------------------------------------------------------------
species_id_section = ""
if species_id_performed:
    try:
        import csv
        from collections import Counter
        with open("${species_top_hits}", 'r') as f:
            rows = list(csv.DictReader(f))
        if rows:
            species_counts = {}
            for col in ['kingdom','phylum','class','order','family','genus','species']:
                if col in rows[0]:
                    vals = [r[col] for r in rows if r.get(col,'').strip()]
                    if vals:
                        top = Counter(vals).most_common(1)
                        if top: species_counts[col] = top[0][0]
            species_id_section  = "## Species Identification (BLAST Analysis)\\n\\n"
            if species_counts:
                species_id_section += "**Taxonomic Classification (Most Likely):**\\n"
                for level, taxon in species_counts.items():
                    if taxon: species_id_section += f"- {level.capitalize()}: {taxon}\\n"
                species_id_section += "\\n"
            n_samples = len(rows)
            species_id_section += f"**Sample Summary:**\\n- Number of samples analyzed: {n_samples}\\n"
            if 'species' in rows[0]:
                spp = [r.get('species','') for r in rows if r.get('species','').strip()]
                if spp:
                    top_spp = Counter(spp).most_common(1)
                    if top_spp:
                        n_agree = top_spp[0][1]
                        species_id_section += f"- Samples with consensus species: {n_agree}/{n_samples} ({n_agree/n_samples*100:.1f}%)\\n"
            species_id_section += (
                "\\n### BLAST Hit Distributions\\n\\n"
                "#### Raw BLAST Results by Taxonomic Level\\n"
                f"![Raw BLAST Pie Charts](${params.outdir}/species_id/blast_raw_pie.png)\\n\\n"
                "#### Summarized Species Identification\\n"
                f"![Summary BLAST Pie Charts](${params.outdir}/species_id/blast_summary_pie.png)\\n\\n"
                "### Detailed Results\\n\\n"
                f"- [Combined BLAST results (TSV)](${params.outdir}/species_id/blast_results.tsv)\\n"
                f"- [Top BLAST hits per sample (CSV)](${params.outdir}/species_id/top_blast_hits.csv)\\n"
                f"- [Posterior probabilities by taxonomic level](${params.outdir}/species_id/blast_posteriors/)\\n"
            )
    except Exception as e:
        print(f"Could not process species ID results: {e}")
        species_id_section = "## Species Identification\\n\\nResults could not be processed.\\n"

# ----------------------------------------------------------------
# Cleaned reads section
# ----------------------------------------------------------------
cleaned_retention_pct = ""
if initial_stats["total"] > 0 and cleaned_stats["total"] > 0:
    cleaned_retention_pct = f"- Retention from raw: {cleaned_stats['total']/initial_stats['total']*100:.1f}% of initial read pairs"

cleaned_section = (
    "## Final Cleaned Reads\\n\\n"
    f"**Summary Statistics (n={cleaned_stats['n']} samples):**\\n"
    f"- Mean reads per sample: {fmt_num(cleaned_stats['mean'])}\\n"
    f"- Standard deviation: {fmt_num(cleaned_stats['sd'])}\\n"
    f"- Min reads: {fmt_num(cleaned_stats['min'])}\\n"
    f"- Max reads: {fmt_num(cleaned_stats['max'])}\\n\\n"
    "**Total Cleaned Output:**\\n"
    f"- Total cleaned read pairs: {fmt_num(cleaned_stats['total'])}\\n"
    f"- Total cleaned bases: {format_bases(total_bases_cleaned)} (assuming {read_length_estimate}bp reads)\\n"
    f"- Total individual reads: {fmt_num(cleaned_stats['total'] * 2)} (paired-end)\\n"
    f"{cleaned_retention_pct}"
)

# ----------------------------------------------------------------
# Mapping section — now includes soft clipping and alignment score
# ----------------------------------------------------------------
mapping_section = ""
if mapping_performed and final_stats["n"] > 0:
    mapped_retention_pct = ""
    if initial_stats["total"] > 0 and final_stats["total"] > 0:
        mapped_retention_pct = f"- Retention from raw: {final_stats['total']/initial_stats['total']*100:.1f}% of initial read pairs"

    insert_size_line = ""
    if insert_size_violin_performed:
        insert_size_line = f"\\n![Insert Size Distribution](${params.outdir}/qc/insert_size_violin.png)\\n"

    soft_clip_line = ""
    if soft_clip_performed:
        soft_clip_line = f"\\n![Soft Clipping Distribution](${params.outdir}/qc/soft_clipping_violin.png)\\n"

    aln_score_line = ""
    if aln_score_performed:
        aln_score_line = f"\\n![Alignment Score Distribution](${params.outdir}/qc/alignment_score_violin.png)\\n"

    mapping_section = (
        "## Mapped Reads\\n\\n"
        f"![Mapped Read Distribution](${params.outdir}/qc/mapped_reads_histogram.png)\\n"
        f"{insert_size_line}"
        f"{soft_clip_line}"
        f"{aln_score_line}"
        "\\n"
        f"**Summary Statistics (n={final_stats['n']} samples):**\\n"
        f"- Mean reads per sample: {fmt_num(final_stats['mean'])}\\n"
        f"- Standard deviation: {fmt_num(final_stats['sd'])}\\n"
        f"- Min reads: {fmt_num(final_stats['min'])}\\n"
        f"- Max reads: {fmt_num(final_stats['max'])}\\n"
        f"- Overall mapping rate: {final_stats['map_pct']:.2f}%\\n\\n"
        "**Total Mapped Output:**\\n"
        f"- Total read pairs mapped: {fmt_num(final_stats['total'])}\\n"
        f"- Total bases mapped: {format_bases(total_bases_mapped)} (assuming {read_length_estimate}bp reads)\\n"
        f"{mapped_retention_pct}\\n\\n"
        "### Final Mapping Statistics\\n"
        f"See [${mapping_summary}](${params.outdir}/qc/${mapping_summary}) for detailed mapping statistics per sample."
    )

# ----------------------------------------------------------------
# Build and write the report
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
*QC Report generated on: {subprocess.check_output(["date"]).decode().strip()}*
\'\'\'

with open("qc_pipeline_report.md", 'w') as f:
    f.write(markdown_content)

print("Markdown report generated successfully!")

try:
    subprocess.run(
        ['pandoc', 'qc_pipeline_report.md', '-o', 'qc_pipeline_report.html',
         '--standalone', '--metadata', 'title=GCL QC Pipeline Report'],
        check=True
    )
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