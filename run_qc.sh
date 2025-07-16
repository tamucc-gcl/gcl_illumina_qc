#!/bin/bash

module load miniconda3
source activate nextflow

nextflow run gcl_illumina_qc/main.nf \
    -profile standard \
    -resume \
    --reads "data/fq_raw/*.{1,2}.fq.gz" \
    --accession "GCA_042920385.1" \
	--decontam_conffile "configs/contam_db.conf" \
    --outdir "results"