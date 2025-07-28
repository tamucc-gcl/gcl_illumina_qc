#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=QC_Handler
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5

#module load miniconda3
#source activate nextflow
module load nextflow/23.10.1
module load miniconda3

accession_number=${1}
# accession_number=GCA_042920385.1

nextflow run gcl_illumina_qc/main.nf \
    -profile slurm \
    -resume \
    --reads "data/fq_raw/*.{1,2}.fq.gz" \
    --accession ${accession_number} \
    --decontam_conffile "configs/contam_db.conf" \
    --outdir "results"
