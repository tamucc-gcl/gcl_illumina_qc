# GCL Illumina QC Pipeline [![DOI](https://zenodo.org/badge/1019725764.svg)](https://doi.org/10.5281/zenodo.19666434)

<!---![](pipeline_dag.png)--->

This pipeline processes demultiplexed Illumina paired-end reads through a series of quality control steps and outputs a set of high-quality reads aligned to a specified reference genome. Each step is modular and follows the [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) standard.

---

## 🔧 QC Steps

The following QC tools are integrated as modules in the pipeline:

1. **[fastp (3' trim)](https://github.com/OpenGene/fastp)** – Removes low-quality bases and trims Illumina adapters from the 3' end.
2. **[Clumpify](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify/)** – De-duplicates reads and optimizes file structure for compression and alignment.
3. **[fastp (5' trim)](https://github.com/OpenGene/fastp)** – Trims polyG/polyX sequences and further filters reads at the 5' end.
4. **[FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)** – Screens reads for contamination by mapping to multiple reference genomes.
5. **[BBMap Repair](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/repair/)** – Ensures read pairing integrity after filtering.
6. **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)** – Quality check after each step to assess sequence quality, duplication levels, and adapter content.
7. **[MultiQC](https://multiqc.info/)** – Aggregates results from FastQC and other tools into a single report.
8. **[NCBI Datasets CLI](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/)** – Downloads reference genomes by accession number.
9. **[BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)** – Aligns reads to the reference genome.
10. **[SAMtools](http://www.htslib.org/)** – Summarizes mapping statistics (`samtools stats` and `flagstat`).

---

## 📁 Expected Directory Structure

Clone this repository into a project folder structured like this:

```
📁 project-root/
├── 📁 gcl_illumina_qc/               # ✅ Must be present (cloned repo)
│   ├── 📄 main.nf
│   ├── 📄 nextflow.config
│   ├── 📄 README.md
│   ├── 📄 run_qc.sbatch
│   └── 📁 modules/
│       ├── clumpify.nf
│       ├── fastp_trim_3.nf
│       ├── fastp_trim_5.nf
│       ├── fastq_screen.nf
│       ├── fetch_genome.nf
│       ├── index_genome.nf
│       ├── map_reads.nf
│       ├── multiqc.nf
│       ├── repair.nf
│       ├── fastqc.nf
│       ├── samtools_stats.nf
│       ├── samtools_summary.nf
│       └── analyze_read_stats.nf
├── 📁 data/
│   ├── 📁 fq_raw/                    # ✅ Must be present (input FASTQ files)
│   ├── 📁 fq_fp1_clmp_fp2_scrn_rpr/         # 🚀 Created by pipeline
│   ├── 📁 bam/                    # 🚀 Created by pipeline
│   └── ... (other subfolders by step)       # 🚀 Created by pipeline
├── 📁 genome/                        # 🚀 Created by pipeline (genome + index)
├── 📁 logs/                          # 🚀 Created by pipeline (SLURM + Nextflow logs)
└── 📁 results/
    └── 📁 multiqc/                   # 🚀 Created by pipeline (MultiQC reports)
```

---

## 🚀 Running the Pipeline

The pipeline accepts reference genomes in two ways:
- **Local genome file** – Use `--genome` to specify a path to a FASTA file on your system
- **NCBI download** – Use `--accession` to automatically download from NCBI

### Option 1: Local/Interactive Run

#### Using a local genome file:
```bash
# Run pipeline with local genome
nextflow run gcl_illumina_qc/main.nf \
    -profile local \
    -resume \
    --reads "data/fq_raw/*.{1,2}.fq.gz" \
    --genome "/path/to/reference/genome.fasta" \
    --sequencing_type "whole_genome" \
    --decontam_conffile "configs/contam_db.conf" \
    --outdir "results"
```

#### Using NCBI accession:
```bash
# Run pipeline with NCBI genome download
nextflow run gcl_illumina_qc/main.nf \
    -profile local \
    -resume \
    --reads "data/fq_raw/*.{1,2}.fq.gz" \
    --accession "GCA_042920385.1" \
    --sequencing_type "whole_genome" \
    --decontam_conffile "configs/contam_db.conf" \
    --outdir "results"
```

#### Optional: Visualize DAG
```bash
nextflow run gcl_illumina_qc/main.nf -with-dag flowchart.dot
dot -Tpng flowchart.dot -o flowchart.png
```

### Option 2: SLURM Submission

The `run_qc.sbatch` script automatically detects whether you're providing a local genome file or an NCBI accession:

#### Using a local genome file:
```bash
# Automatically detects file and uses --genome parameter
sbatch run_qc.sbatch /path/to/genome.fasta ddrad configs/contam_db.conf
sbatch run_qc.sbatch ./references/my_species.fa whole_genome configs/contam_db.conf
```

#### Using NCBI accession:
```bash
# Automatically detects accession pattern and uses --accession parameter
sbatch run_qc.sbatch GCA_042920385.1 whole_genome configs/contam_db.conf
sbatch run_qc.sbatch GCF_000001405.40 ddrad configs/contam_db.conf
```

Edit `gcl_illumina_qc/run_qc.sbatch` to customize resources used by orchestra conductor job (e.g., CPUs, memory, partition).

Edit `gcl_illumina_qc/nextflow.config` to customize resources used by each QC stage (e.g., CPUs, memory, partition).

### 📝 Parameter Reference

| Parameter | Description | Required | Example |
|-----------|-------------|----------|---------|
| `--reads` | Input paired-end FASTQ files | ✅ Yes | `"data/fq_raw/*.{1,2}.fq.gz"` |
| `--genome` | Path to local genome FASTA file | ⚠️ One required | `/path/to/genome.fa` |
| `--accession` | NCBI assembly accession | ⚠️ One required | `GCA_042920385.1` |
| `--sequencing_type` | Library type | 🔧 Optional | `ddrad` or `whole_genome` |
| `--decontam_conffile` | FastQ Screen config | ✅ Yes | `configs/contam_db.conf` |
| `--outdir` | Output directory | 🔧 Optional | `results` |

**Note:** Either `--genome` OR `--accession` must be specified, but not both.

### 🧬 Supported Genome Formats

When using local genome files (`--genome`), the following formats are supported:
- `.fa` – FASTA format
- `.fasta` – FASTA format  
- `.fna` – FASTA nucleotide format
- `.fa.gz` – Gzipped FASTA
- `.fasta.gz` – Gzipped FASTA
- `.fna.gz` – Gzipped FASTA nucleotide

---

## 📊 Post-Pipeline QC Analysis

The final step of the pipeline summarizes alignment and quality metrics across all samples. This includes:

- **Read Retention Statistics** – Compiled at each QC step (via `analyze_read_stats.nf`).
- **Mapping Metrics** – Generated using `samtools flagstat` and `samtools stats`.
- **MultiQC** – Integrates all FastQC, fastp, and mapping metrics into one HTML report (`results/multiqc/`).

These reports help evaluate the effectiveness of each filtering and trimming step, ensuring high-quality downstream analyses.

---

## 💡 Notes

- Each module runs independently with `publishDir` to maintain organized outputs.
- The pipeline can be extended with additional modules (e.g., variant calling or quantification).
- Uses Conda environments (no containers required by default).
- If job is interupted prior to finishing the full pipeline it can be resumed from where it was interupted with the same command that started the run
- After satisfactory completion the `./work` directory can be deleted to free-up disk-space. The sub-folders in `./data/fq_fp1*` can also be removed as intermediate files.

---

## 📬 Contact

For questions or feedback, contact Jason Selwyn — [jason.selwyn@tamucc.edu](mailto:jason.selwyn@tamucc.edu)
