# gcl_illumina_qc
![](pipeline_dag.png)
## Folder Structure
```
📁 project-root/
├── 📁 gcl_illumina_qc/               # All Nextflow workflow logic
│   ├── 📄 main.nf                    # Entry point for pipeline
│   ├── 📄 nextflow.config            # SLURM and container config
│   ├── 📄 README.md                  # Usage and documentation
│   └── 📁 modules/                   # Individual DSL2 module processes
│       ├── clumpify.nf
│       ├── fastp_trim_3.nf
│       ├── fastp_trim_5.nf
│       ├── fastq_screen.nf
│       ├── fetch_genome.nf
│       ├── index_genome.nf
│       ├── map_reads.nf
│       ├── multiqc.nf
│       └── repair.nf
├── 📁 data/                          # Input FASTQ files and QC outputs
├── 📁 genome/                        # Downloaded and indexed reference
├── 📁 logs/                          # SLURM/Nextflow job logs
└── 📁 results/                       # Final results (e.g. MultiQC reports)
    └── 📁 multiqc/                   # Per-step and consolidated MultiQC outputs
```

## To Run

1. Clone repo into your directory.
2. Run the pipeline from the project root:

```bash
module load miniconda3
source activate nextflow
nextflow run gcl_illumina_qc/main.nf -resume -with-dag flowchart.dot
dot -Tpng flowchart.dot -o flowchart.png

nextflow run gcl_illumina_qc/main.nf \
    -profile standard \
    -resume \
    --reads "data/fq_raw/*.{1,2}.fq.gz" \
    --accession "GCA_042920385.1" \
	--decontam_conffile "configs/contam_db.conf" \
    --outdir "results"
```
3. Or Run in SLURM `sbatch run_qc.sbatch`
