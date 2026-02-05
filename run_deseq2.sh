#!/bin/bash
#SBATCH --job-name=DESeq2_analysis
#SBATCH --output=logs/deseq2_%A.out
#SBATCH --error=logs/deseq2_%A.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mirnakg@mit.edu

module load miniconda3/v4
source $(conda info --base)/etc/profile.d/conda.sh
conda activate sra-tools

echo ""
echo "DESeq2 Differential Expression Analysis"
echo "Start time: $(date)"
echo ""

# Make sure output directories exist
mkdir -p /net/bmc-pub14/data/voigt/users/mirnakg/RNA_seq_pipe/test_samples/deseq2_results
mkdir -p logs

# Run the analysis
python3 deseq2_analysis.py

EXIT_CODE=$?

echo ""
if [ ${EXIT_CODE} -eq 0 ]; then
    echo "Analysis completed successfully!"
    echo ""
    echo "Results available in:"
    echo "  test_samples/deseq2_results/"
else
    echo "Analysis failed with exit code ${EXIT_CODE}"
    echo "Check the error log for details"
fi
echo "End time: $(date)"

conda deactivate
