# gcl_illumina_qc

## Folder Structure
project-root/
├── gcl_illumina_qc/              # All Nextflow logic goes here
│   ├── main.nf                   # Entry script
│   ├── nextflow.config           # Cluster + container config
│   ├── README.md                 # Optional usage docs
│   └── modules/                  # All DSL2 process scripts
│       ├── clumpify.nf
│       ├── fastp_trim_3.nf
│       ├── fastp_trim_5.nf
│       ├── fastq_screen.nf
│       ├── fetch_genome.nf
│       ├── index_genome.nf
│       ├── map_reads.nf
│       ├── multiqc.nf
│       └── repair.nf
│
├── data/                         # Input FASTQs and intermediate files
│
├── genome/                       # Genome files and indexes
│
├── logs/                         # Logs from SLURM and Nextflow
│
├── results/                      # Final outputs
│   └── multiqc/                  # Consolidated MultiQC reports


## To Run 
1. Clone repo into your directory
2. Run code
    - `nextflow run gcl_illumina_qc/main.nf -profile standard -resume`