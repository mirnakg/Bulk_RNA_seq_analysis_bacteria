#!/bin/bash
#SBATCH --job-name=FASTQ_to_BAM_bowtie2
#SBATCH --output=logs/fastq_to_bam_bowtie2_%A_%a.out
#SBATCH --error=logs/fastq_to_bam_bowtie2_%A_%a.err
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mirnakg@mit.edu

# Activate conda environment
module load miniconda3/v4
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra-tools

# Specify directories
GENOME_INDEX="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/genome_index_bowtie2/reference_genome"
FASTQ_DIR="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/trimmed"
OUTPUT_DIR="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/bam"
THREADS=8

# Create output directory
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Get sample for this array task
sample=$(sed -n "${SLURM_ARRAY_TASK_ID}p" /net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/sample_names.txt)


echo "Processing BAM generation for ${sample} (line ${SLURM_ARRAY_TASK_ID})"
echo "BACTERIAL RNA-seq with Bowtie2"
echo "Job ID: ${SLURM_JOB_ID}, Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Start time: $(date)"
echo "______________________________"

# Check if FASTQ files exist
if [[ ! -f "${FASTQ_DIR}/${sample}_R1_paired.fastq.gz" ]] || [[ ! -f "${FASTQ_DIR}/${sample}_R2_paired.fastq.gz" ]]; then
    echo "ERROR: Trimmed FASTQ files not found for ${sample}"
    echo "  Expected: ${FASTQ_DIR}/${sample}_R1_paired.fastq.gz"
    echo "           ${FASTQ_DIR}/${sample}_R2_paired.fastq.gz"
    echo "Run Trimmomatic first!"
    exit 0
fi

echo "Found trimmed FASTQ files for ${sample}, proceeding with alignment..."

# Run Bowtie2 alignment
echo "Running Bowtie2 alignment..."
bowtie2 \
    -x ${GENOME_INDEX} \
    -1 ${FASTQ_DIR}/${sample}_R1_paired.fastq.gz \
    -2 ${FASTQ_DIR}/${sample}_R2_paired.fastq.gz \
    -S ${OUTPUT_DIR}/${sample}.sam \
    --threads ${THREADS} \
    --very-sensitive \
    2> ${OUTPUT_DIR}/${sample}_bowtie2.log

if [ $? -ne 0 ]; then
    echo "ERROR: Bowtie2 alignment failed for ${sample}"
    exit 1
fi

echo "Bowtie2 alignment completed successfully!"

# Convert SAM to BAM
echo "Converting SAM to BAM..."
samtools view \
    -@ ${THREADS} \
    -b \
    -o ${OUTPUT_DIR}/${sample}.bam \
    ${OUTPUT_DIR}/${sample}.sam

if [ $? -ne 0 ]; then
    echo "ERROR: SAM to BAM conversion failed for ${sample}"
    exit 1
fi

# Remove SAM file to save space
rm ${OUTPUT_DIR}/${sample}.sam

# Sort BAM file
echo "Sorting BAM file..."
samtools sort \
    -@ ${THREADS} \
    -o ${OUTPUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam \
    ${OUTPUT_DIR}/${sample}.bam

if [ $? -ne 0 ]; then
    echo "ERROR: BAM sorting failed for ${sample}"
    exit 1
fi

# Remove unsorted BAM
rm ${OUTPUT_DIR}/${sample}.bam

# Index sorted BAM
echo "Indexing sorted BAM..."
samtools index \
    -@ ${THREADS} \
    ${OUTPUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam

if [ $? -ne 0 ]; then
    echo "ERROR: BAM indexing failed for ${sample}"
    exit 1
fi

# Generate alignment statistics
echo ""
echo "Alignment statistics:"
samtools flagstat ${OUTPUT_DIR}/${sample}_Aligned.sortedByCoord.out.bam | tee ${OUTPUT_DIR}/${sample}_flagstat.txt

echo ""
echo "Bowtie2 alignment summary:"
cat ${OUTPUT_DIR}/${sample}_bowtie2.log

echo "-------------------------------------"
echo "BAM generation completed for ${sample}"
echo "End time: $(date)"
echo "-----------------------------------"

conda deactivate
