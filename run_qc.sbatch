#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=QC_Handler
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5

module load miniconda3
source activate nextflow

nextflow run gcl_illumina_qc/main.nf \
    -profile standard \
    -resume \
    --reads "data/fq_raw/*.{1,2}.fq.gz" \
    --accession "GCA_042920385.1" \
	  --decontam_conffile "configs/contam_db.conf" \
    --outdir "results"
