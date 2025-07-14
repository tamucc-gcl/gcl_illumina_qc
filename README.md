# gcl_illumina_qc

## Folder Structure
<pre> 📁 <b>project-root/</b> ├── 📁 <b>gcl_illumina_qc/</b> # All Nextflow workflow logic │ ├── 📄 main.nf # Entry point for pipeline │ ├── 📄 nextflow.config # SLURM and container config │ ├── 📄 README.md # Usage and documentation │ └── 📁 modules/ # Individual DSL2 module processes │ ├── clumpify.nf │ ├── fastp_trim_3.nf │ ├── fastp_trim_5.nf │ ├── fastq_screen.nf │ ├── fetch_genome.nf │ ├── index_genome.nf │ ├── map_reads.nf │ ├── multiqc.nf │ └── repair.nf │ ├── 📁 data/ # Input FASTQ files and QC outputs ├── 📁 genome/ # Downloaded and indexed reference ├── 📁 logs/ # SLURM/Nextflow job logs └── 📁 results/ # Final results (e.g. MultiQC reports) └── 📁 multiqc/ # Per-step and consolidated MultiQC outputs </pre>

## To Run 
1. Clone repo into your directory
2. Run code
    - `nextflow run gcl_illumina_qc/main.nf -profile standard -resume`