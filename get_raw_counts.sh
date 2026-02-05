#!/bin/bash
#SBATCH --job-name=featureCounts
#SBATCH --output=logs/featurecounts_%A_%a.out
#SBATCH --error=logs/featurecounts_%A_%a.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mirnakg@mit.edu

module load miniconda3/v4
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra-tools

BAM_DIR="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/bam"
GTF_FILE="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/genome_annotations.gtf"
OUTPUT_DIR="/net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/counts"
THREADS=8

mkdir -p ${OUTPUT_DIR}
mkdir -p logs

echo ""
echo "Running featureCounts for Acinetobacter baylyi"
echo "Job ID: ${SLURM_JOB_ID}"
echo "Start time: $(date)"
echo "_________________________________________"

BAM_FILES=$(ls ${BAM_DIR}/*_Aligned.sortedByCoord.out.bam)

if [ -z "$BAM_FILES" ]; then
    echo "ERROR: No BAM files found in ${BAM_DIR}"
    exit 1
fi

echo "Found BAM files:"
echo "$BAM_FILES"
echo ""

# Primary count: Use CDS features (protein-coding genes)
echo "Running featureCounts on CDS features..."
featureCounts \
    -T ${THREADS} \
    -p \
    -B \
    -C \
    -t CDS \
    -g gene_id \
    -s 0 \
    --minOverlap 10 \
    --fracOverlap 0.2 \
    -Q 10 \
    --primary \
    -a ${GTF_FILE} \
    -o ${OUTPUT_DIR}/raw_counts_CDS.txt \
    ${BAM_FILES}

# Alternative count: Use gene features (entire gene locus)
echo ""
echo "Also running on gene features for comparison..."
featureCounts \
    -T ${THREADS} \
    -p \
    -B \
    -C \
    -t gene \
    -g gene_id \
    -s 0 \
    --minOverlap 10 \
    --fracOverlap 0.2 \
    -Q 10 \
    --primary \
    -a ${GTF_FILE} \
    -o ${OUTPUT_DIR}/raw_counts_gene.txt \
    ${BAM_FILES}

# Diagnostic: Count with multi-mapping allowed
echo ""
echo "Running diagnostic count with multi-mappers allowed..."
featureCounts \
    -T ${THREADS} \
    -p \
    -B \
    -C \
    -M \
    --fraction \
    -t CDS \
    -g gene_id \
    -s 0 \
    --minOverlap 10 \
    -Q 10 \
    -a ${GTF_FILE} \
    -o ${OUTPUT_DIR}/raw_counts_CDS_with_multimappers.txt \
    ${BAM_FILES}

if [ $? -eq 0 ]; then
    echo ""
    echo "----------------------------------"
    echo "featureCounts completed successfully!"
    echo "Main output: ${OUTPUT_DIR}/raw_counts_CDS.txt"
    echo "Alternative: ${OUTPUT_DIR}/raw_counts_gene.txt"
    echo "With multimappers: ${OUTPUT_DIR}/raw_counts_CDS_with_multimappers.txt"
    echo "----------------------------------"
    
    # Display summaries
    echo ""
    echo "CDS counts summary:"
    cat ${OUTPUT_DIR}/raw_counts_CDS.txt.summary
    
    echo ""
    echo "With multimappers summary:"
    cat ${OUTPUT_DIR}/raw_counts_CDS_with_multimappers.txt.summary
else
    echo "ERROR: featureCounts failed"
    exit 1
fi

echo ""
echo "End time: $(date)"
echo ""

conda deactivate
