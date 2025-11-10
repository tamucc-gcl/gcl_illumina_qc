#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=QC_Handler
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5

# Load required modules
module load nextflow/23.10.1
module load miniconda3

# Parse arguments
genome_or_accession=${1}
sequencing_type=${2}  # ddrad or whole_genome
contamination_file=${3}

# Instructions
# Example usage for run_qc.sbatch

echo "==========================================="
echo "QC Pipeline - Usage Examples"
echo "==========================================="
echo ""
echo "The script automatically detects whether you're using a local genome file,"
echo "an NCBI accession number, de novo assembly mode, or no-genome mode."
echo ""
echo "SYNTAX:"
echo "  sbatch run_qc.sbatch <genome_or_accession|denovo|none> <sequencing_type> <contamination_config>"
echo ""
echo "-------------------------------------------"
echo "EXAMPLE 1: Using a local genome file"
echo "-------------------------------------------"
echo "sbatch run_qc.sbatch ./genomes/my_species.fasta ddrad configs/contam_db.conf"
echo ""
echo "sbatch run_qc.sbatch /shared/references/genome_v2.fa whole_genome configs/contam_db.conf"
echo ""
echo "sbatch run_qc.sbatch ~/project/ref.fa.gz ddrad configs/contam_db.conf"
echo ""
echo "-------------------------------------------"
echo "EXAMPLE 2: Using NCBI accession"
echo "-------------------------------------------"
echo "sbatch run_qc.sbatch GCA_042920385.1 whole_genome configs/contam_db.conf"
echo ""
echo "sbatch run_qc.sbatch GCF_000001405.40 ddrad configs/contam_db.conf"
echo ""
echo "-------------------------------------------"
echo "EXAMPLE 3: De novo assembly mode"
echo "-------------------------------------------"
echo "sbatch run_qc.sbatch denovo ddrad configs/contam_db.conf"
echo ""
echo "sbatch run_qc.sbatch denovo whole_genome configs/contam_db.conf"
echo ""
echo "-------------------------------------------"
echo "EXAMPLE 4: No genome mode (output cleaned reads only)"
echo "-------------------------------------------"
echo "sbatch run_qc.sbatch none ddrad configs/contam_db.conf"
echo ""
echo "sbatch run_qc.sbatch none whole_genome configs/contam_db.conf"
echo ""
echo "-------------------------------------------"
echo "PARAMETERS:"
echo "-------------------------------------------"
echo "  Arg 1: genome_or_accession"
echo "    - Local file path (e.g., ./genome.fa, /path/to/genome.fasta)"
echo "    - NCBI accession (e.g., GCA_042920385.1, GCF_000001405.40)"
echo "    - 'denovo' for de novo assembly from cleaned reads"
echo "    - 'none' for no genome (output cleaned reads only)"
echo ""
echo "  Arg 2: sequencing_type"
echo "    - 'ddrad' for ddRAD sequencing"
echo "    - 'whole_genome' for whole genome sequencing"
echo ""
echo "  Arg 3: contamination_config"
echo "    - Path to FastQ Screen configuration file"
echo "    - e.g., configs/contam_db.conf"
echo ""
echo "-------------------------------------------"
echo "SUPPORTED GENOME FILE FORMATS:"
echo "-------------------------------------------"
echo "  - .fa      (FASTA)"
echo "  - .fasta   (FASTA)"
echo "  - .fna     (FASTA nucleotide)"
echo "  - .fa.gz   (Gzipped FASTA)"
echo "  - .fasta.gz (Gzipped FASTA)"
echo "  - .fna.gz  (Gzipped FASTA nucleotide)"
echo ""
echo "-------------------------------------------"
echo "AUTO-DETECTION LOGIC:"
echo "-------------------------------------------"
echo "1. If argument 1 is 'denovo' → perform de novo assembly"
echo "2. If argument 1 is 'none' → no genome, output cleaned reads"
echo "3. If argument 1 is a file that exists → use --genome"
echo "4. If argument 1 matches GC[AF]_*.* pattern → use --accession"
echo "5. Otherwise → assume it's an accession (let Nextflow validate)"
echo ""
echo "==========================================="

# Validate arguments
if [ -z "$genome_or_accession" ]; then
    echo "Error: First argument (genome file, NCBI accession, 'denovo', or 'none') is required"
    echo "Usage: sbatch run_qc.sbatch <genome_file|accession|denovo|none> <sequencing_type> <contamination_config>"
    echo "Example (local genome): sbatch run_qc.sbatch /path/to/genome.fasta ddrad configs/contam_db.conf"
    echo "Example (NCBI): sbatch run_qc.sbatch GCA_042920385.1 whole_genome configs/contam_db.conf"
    echo "Example (de novo): sbatch run_qc.sbatch denovo ddrad configs/contam_db.conf"
    echo "Example (no genome): sbatch run_qc.sbatch none whole_genome configs/contam_db.conf"
    exit 1
fi

if [ -z "$sequencing_type" ]; then
    echo "Error: Second argument (sequencing_type) is required"
    echo "Options: ddrad or whole_genome"
    exit 1
fi

if [ -z "$contamination_file" ]; then
    echo "Error: Third argument (contamination config file) is required"
    exit 1
fi

# Validate sequencing type
if [[ "$sequencing_type" != "ddrad" && "$sequencing_type" != "whole_genome" ]]; then
    echo "Error: sequencing_type must be either 'ddrad' or 'whole_genome'"
    echo "You provided: $sequencing_type"
    exit 1
fi

# Validate contamination file exists
if [ ! -f "$contamination_file" ]; then
    echo "Error: Contamination config file not found: $contamination_file"
    exit 1
fi

# Determine mode and set parameters
if [ "$genome_or_accession" == "denovo" ]; then
    # De novo assembly mode
    echo "==========================================="
    echo "De novo assembly mode selected"
    echo "==========================================="
    
    genome_param=""
    assembly_param="--assembly_mode denovo"
    genome_info="Mode: De novo assembly from cleaned reads"
    
    # Additional message for ddRAD
    if [ "$sequencing_type" == "ddrad" ]; then
        echo "Using Rainbow assembly pipeline optimized for ddRAD data"
    fi

elif [ "$genome_or_accession" == "none" ]; then
    # No genome mode
    echo "==========================================="
    echo "No-genome mode selected"
    echo "==========================================="
    
    genome_param=""
    assembly_param="--assembly_mode none"
    genome_info="Mode: No reference genome (output cleaned reads only)"

elif [ -f "$genome_or_accession" ]; then
    # It's a file that exists - use as local genome
    echo "==========================================="
    echo "Detected local genome file: $genome_or_accession"
    echo "==========================================="

    # Get absolute path to ensure it works from any directory
    genome_path=$(realpath "$genome_or_accession")

    # Validate it looks like a genome file
    if [[ ! "$genome_path" =~ \.(fa|fasta|fna|fa\.gz|fasta\.gz|fna\.gz)$ ]]; then
        echo "Warning: File doesn't have typical genome extension (.fa, .fasta, .fna, or .gz)"
        echo "Proceeding anyway..."
    fi

    genome_param="--genome ${genome_path}"
    assembly_param="--assembly_mode none"
    genome_info="Local genome: ${genome_path}"

elif [[ "$genome_or_accession" =~ ^GC[AF]_[0-9]+\.[0-9]+$ ]]; then
    # It matches NCBI accession pattern (GCA_ or GCF_ followed by numbers)
    echo "==========================================="
    echo "Detected NCBI accession: $genome_or_accession"
    echo "==========================================="

    genome_param="--accession ${genome_or_accession}"
    assembly_param="--assembly_mode none"
    genome_info="NCBI accession: ${genome_or_accession}"

elif [ -e "$genome_or_accession" ]; then
    # Path exists but is not a regular file (might be directory, symlink issue, etc.)
    echo "Error: Path exists but is not a regular file: $genome_or_accession"
    exit 1

else
    # Doesn't exist as file and doesn't match accession pattern
    # Assume it's an accession (let Nextflow handle validation)
    echo "==========================================="
    echo "File not found, assuming NCBI accession: $genome_or_accession"
    echo "==========================================="

    genome_param="--accession ${genome_or_accession}"
    assembly_param="--assembly_mode none"
    genome_info="NCBI accession (assumed): ${genome_or_accession}"
fi

# Print run configuration
echo "==========================================="
echo "QC Pipeline Configuration:"
echo "-------------------------------------------"
echo "$genome_info"
echo "Sequencing type: $sequencing_type"
echo "Contamination config: $contamination_file"
echo "Output directory: results"

# Print de novo assembly parameters if applicable
if [ "$genome_or_accession" == "denovo" ]; then
    echo "-------------------------------------------"
    echo "De novo assembly parameters:"
    echo "  Minimum reads per individual: 4"
    echo "  Minimum individuals: 4"
    echo "  CD-HIT similarity: 0.8"
    echo "  Rainbow div -f: 0.5"
    echo "  Rainbow div -K: 10"
    echo "  (Modify in nextflow.config to adjust)"
fi

echo "==========================================="
echo ""

# Run the pipeline with appropriate genome parameter
nextflow run gcl_illumina_qc/main.nf \
    -profile slurm \
    -resume \
    --reads "data/fq_raw/*.{1,2}.fq.gz" \
    ${genome_param} \
    ${assembly_param} \
    --sequencing_type ${sequencing_type} \
    --decontam_conffile ${contamination_file} \
    --outdir "results"

# Capture exit code
exit_code=$?

# Print completion status
if [ $exit_code -eq 0 ]; then
    echo ""
    echo "==========================================="
    echo "Pipeline completed successfully!"
    
    # Additional messages based on mode
    if [ "$genome_or_accession" == "denovo" ]; then
        echo "De novo assembly saved to: results/denovo_assembly/"
        echo "BAM files mapped to assembly: data/bam/"
    elif [ "$genome_or_accession" == "none" ]; then
        echo "Cleaned reads saved to: results/cleaned_reads/"
    else
        echo "BAM files saved to: data/bam/"
    fi
    
    echo "QC reports saved to: results/"
    echo "==========================================="
else
    echo ""
    echo "==========================================="
    echo "Pipeline failed with exit code: $exit_code"
    echo "Check logs for details"
    echo "==========================================="
fi

exit $exit_code