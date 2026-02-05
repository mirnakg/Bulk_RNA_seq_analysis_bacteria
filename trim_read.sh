#!/bin/bash
#SBATCH --job-name=Trimmomatic_bacterial
#SBATCH --output=logs/trimmomatic_%A_%a.out
#SBATCH --error=logs/trimmomatic_%A_%a.err
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mirnakg@mit.edu

# Activate conda environment
module load miniconda3/v4
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra-tools

# Specify directories
FASTQ_DIR="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples"
OUTPUT_DIR="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/trimmed"
THREADS=8

# Trimmomatic adapter file
ADAPTERS="$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa"

# Create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Get sample for this array task
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/sample_names.txt)

echo "-----------------------------------"
echo "Processing Trimmomatic for ${sample} (line ${SLURM_ARRAY_TASK_ID})"
echo "Job ID: ${SLURM_JOB_ID}, Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Start time: $(date)"
echo "--------------------------------"

# Check if FASTQ files exist
if [[ ! -f "${FASTQ_DIR}/${sample}_R1_001.fastq.gz" ]] || [[ ! -f "${FASTQ_DIR}/${sample}_R2_001.fastq.gz" ]]; then
    echo "ERROR: FASTQ files not found for ${sample}"
    echo "  Expected: ${FASTQ_DIR}/${sample}_R1_001.fastq.gz"
    echo "           ${FASTQ_DIR}/${sample}_R2_001.fastq.gz"
    exit 0
fi

# Check if already processed
if [[ -f "${OUTPUT_DIR}/${sample}_R1_paired.fastq.gz" ]] && [[ -f "${OUTPUT_DIR}/${sample}_R2_paired.fastq.gz" ]]; then
    echo "Trimmed files already exist for ${sample}, skipping..."
    echo "  ${OUTPUT_DIR}/${sample}_R1_paired.fastq.gz"
    echo "  ${OUTPUT_DIR}/${sample}_R2_paired.fastq.gz"
    exit 0
fi

echo "Found FASTQ files for ${sample}, proceeding with trimming..."

# Check if adapter file exists
if [[ ! -f "${ADAPTERS}" ]]; then
    echo "ERROR: Adapter file not found: ${ADAPTERS}"
    echo "Searching for adapter files..."
    find $CONDA_PREFIX -name "*.fa" -path "*/adapters/*" 2>/dev/null
    exit 1
fi

echo "Using adapter file: ${ADAPTERS}"
echo ""

# Set Java memory options for Trimmomatic
export _JAVA_OPTIONS="-Xmx24g"

# Run Trimmomatic
echo "Running Trimmomatic PE..."
trimmomatic PE \
    -threads ${THREADS} \
    -phred33 \
    ${FASTQ_DIR}/${sample}_R1_001.fastq.gz \
    ${FASTQ_DIR}/${sample}_R2_001.fastq.gz \
    ${OUTPUT_DIR}/${sample}_R1_paired.fastq.gz \
    ${OUTPUT_DIR}/${sample}_R1_unpaired.fastq.gz \
    ${OUTPUT_DIR}/${sample}_R2_paired.fastq.gz \
    ${OUTPUT_DIR}/${sample}_R2_unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS}:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36 \
    2>&1 | tee ${OUTPUT_DIR}/${sample}_trimmomatic.log

# Unset Java options
unset _JAVA_OPTIONS

EXIT_CODE=${PIPESTATUS[0]}

if [ ${EXIT_CODE} -eq 0 ]; then
    echo ""
    echo "Trimmomatic completed successfully!"
    echo ""
    echo "Output files created:"
    ls -lh ${OUTPUT_DIR}/${sample}_R*_paired.fastq.gz
    echo ""
    
    # Extract summary statistics
    echo "Trimming summary:"
    grep "Input Read Pairs" ${OUTPUT_DIR}/${sample}_trimmomatic.log || echo "Check log file for details"
    
else
    echo "ERROR: Trimmomatic failed for ${sample}"
    exit 1
fi

echo "--------------------------"
echo "Trimmomatic completed for ${sample}"
echo "End time: $(date)"
echo "---------------------------"

conda deactivate
