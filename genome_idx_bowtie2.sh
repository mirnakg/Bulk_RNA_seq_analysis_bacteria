#!/bin/bash
#SBATCH --job-name=Bowtie2_idx_bacterial
#SBATCH --output=logs/bowtie2_idx_bacterial_%A.out
#SBATCH --error=logs/bowtie2_idx_bacterial_%A.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mirnakg@mit.edu

# Activate conda environment
module load miniconda3/v4
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra-tools

# --------------------------------
# Input files
# -------------------------------
GENOME_FA="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/ref_genome.fa"
INDEX_DIR="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/genome_index_bowtie2"
INDEX_PREFIX="${INDEX_DIR}/reference_genome"

mkdir -p ${INDEX_DIR}
mkdir -p logs

echo "---------------------------------------"
echo "Bowtie2 Genome Index Generation - BACTERIAL"
echo "Organism: Acinetobacter baylyi ADP1"
echo "Start time: $(date)"
echo "--------------------------------------"
echo ""

# -------------------------------
# Validate files
# ---------------------------------
echo "Checking reference genome:"
echo "  ${GENOME_FA}"
if [ ! -f "${GENOME_FA}" ]; then
    echo "ERROR: Reference genome FASTA not found."
    exit 1
fi

echo ""

# -------------------------------
# Run Bowtie2 Index Build
# -------------------------------
echo "Building Bowtie2 index..."
echo "Index prefix: ${INDEX_PREFIX}"
echo ""

bowtie2-build \
    --threads 16 \
    ${GENOME_FA} \
    ${INDEX_PREFIX}

EXIT_CODE=$?

echo ""
echo "-----------------------------"
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Bowtie2 genome index generation COMPLETED"
    echo "Index directory contents:"
    ls -lh ${INDEX_DIR}
    echo ""
    echo "Index files created:"
    ls -lh ${INDEX_PREFIX}*.bt2
else
    echo "Bowtie2 genome index generation FAILED"
    echo "Check the error log."
fi
echo "End time: $(date)"
echo "-----------------------------"

conda deactivate
