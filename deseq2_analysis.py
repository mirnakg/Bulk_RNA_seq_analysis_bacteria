#!/usr/bin/env python3
"""
DESeq2 pairwise analysis for all groups
Automatically detects all conditions and performs all pairwise comparisons
"""

import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations

# Set plot style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (10, 8)


print("DESeq2 Differential Expression Analysis")
print()

########
# LOAD AND PREPARE DATA
############

print("Loading count data")
counts_file = "/net/bmc-pub14/data/voigt/users/mirnakg/codon_project/bsubtilis/fastq/counts/raw_counts_gene.txt"

# Read the file, skipping the first comment line
counts_df = pd.read_csv(counts_file, sep='\t', comment='#', skiprows=1)

print(f"Raw data shape: {counts_df.shape}")
print(f"Columns: {counts_df.columns.tolist()}")
print()

# Extract count columns
gene_ids = counts_df['Geneid']
count_columns = counts_df.columns[6:]

counts_matrix = counts_df[count_columns].copy()
counts_matrix.index = gene_ids

print(f"Count matrix shape: {counts_matrix.shape}")
print(f"Genes: {counts_matrix.shape[0]}, Samples: {counts_matrix.shape[1]}")
print()

# Clean up sample names
counts_matrix.columns = [col.split('/')[-1].replace('_Aligned.sortedByCoord.out.bam', '') 
                         for col in counts_matrix.columns]

print("Sample names:")
print(counts_matrix.columns.tolist())
print()

########
# AUTO-DETECT CONDITIONS
##########

print("Creating metadata...")

# OPTION 1: Auto-detect from sample names (assumes pattern like: condition1_rep1, condition2_rep1)
# You can modify this logic based on your naming convention
sample_names = counts_matrix.columns.tolist()

# Example: If your samples follow a pattern, extract conditions
# For now, using the original pattern (alternating A/B)
# MODIFY THIS SECTION based on your actual sample naming scheme
conditions = ['15h', '12h', '12h', '12h', '10h', '10h', '10h', '8h', '24h', '24h', '24h', '18h', '18h', '18h', '15h', '15h', '8h', '8h']

# ALTERNATIVE: If you want to specify conditions manually, use this instead:
# conditions = ['control', 'treatment1', 'treatment2', 'control', 'treatment1', 'treatment2']

metadata = pd.DataFrame({
    'sample': sample_names,
    'condition': conditions
})
metadata.index = metadata['sample']

print("Metadata:")
print(metadata)
print()

# Get unique conditions
unique_conditions = metadata['condition'].unique()
print(f"Detected conditions: {list(unique_conditions)}")
print(f"Number of conditions: {len(unique_conditions)}")
print()

# Calculate number of pairwise comparisons
n_comparisons = len(list(combinations(unique_conditions, 2)))
print(f"Will perform {n_comparisons} pairwise comparisons")
print()

# Transpose for PyDESeq2 (samples as rows, genes as columns)
counts_matrix_transposed = counts_matrix.T
print(f"Transposed matrix shape (samples x genes): {counts_matrix_transposed.shape}")
print()

############
# FILTER LOW-COUNT GENES
############

print("Filtering low-count genes")
print(f"Genes before filtering: {counts_matrix_transposed.shape[1]}")

min_counts = 10
min_samples = 3
keep = (counts_matrix_transposed >= min_counts).sum(axis=0) >= min_samples
counts_filtered = counts_matrix_transposed.loc[:, keep]

print(f"Genes after filtering: {counts_filtered.shape[1]}")
print(f"Removed {counts_matrix_transposed.shape[1] - counts_filtered.shape[1]} low-count genes")
print()

###########
# RUN DESEQ2 
##########

print("Running DESeq2")
print()

# Create DESeq2 dataset with ALL samples
dds = DeseqDataSet(
    counts=counts_filtered,
    metadata=metadata,
    design_factors="condition",
    refit_cooks=True,
    n_cpus=8
)

# Run DESeq2 analysis
dds.deseq2()

print("DESeq2 analysis complete!")
print()

# #########
# PERFORM ALL PAIRWISE COMPARISONS
# #########

print("PERFORMING PAIRWISE COMPARISONS")
print()

# Create output directory
output_dir = "/net/bmc-pub14/data/voigt/users/mirnakg/codon_project/bsubtilis/fastq/deseq2_results"
import os
os.makedirs(output_dir, exist_ok=True)

# Create subdirectory for pairwise comparisons
pairwise_dir = f"{output_dir}/pairwise_comparisons"
os.makedirs(pairwise_dir, exist_ok=True)

# Store all results for summary
all_comparisons = []
all_pvalues = []  # For global FDR correction

# Generate all pairwise combinations
pairwise_combos = list(combinations(unique_conditions, 2))

for i, (cond1, cond2) in enumerate(pairwise_combos, 1):
    print(f"[{i}/{n_comparisons}] Comparing: {cond2} vs {cond1}")
    print()
    
    # Run statistical test for this comparison
    stat_res = DeseqStats(dds, contrast=["condition", cond2, cond1])
    stat_res.summary()
    
    # Get results
    results = stat_res.results_df.copy()
    results['gene_id'] = results.index
    results['comparison'] = f"{cond2}_vs_{cond1}"
    
    # Store for global correction
    all_pvalues.extend(results['pvalue'].tolist())
    
    # Sort by adjusted p-value
    results_sorted = results.sort_values('padj')
    
    # Summary statistics
    sig_05 = results[(results['padj'] < 0.05) & (abs(results['log2FoldChange']) > 1)]
    sig_01 = results[(results['padj'] < 0.01) & (abs(results['log2FoldChange']) > 1)]
    up_01 = results[(results['padj'] < 0.01) & (results['log2FoldChange'] > 1)]
    down_01 = results[(results['padj'] < 0.01) & (results['log2FoldChange'] < -1)]
    
    print(f"Significant genes (padj < 0.05, |log2FC| > 1): {len(sig_05)}")
    print(f"Significant genes (padj < 0.01, |log2FC| > 1): {len(sig_01)}")
    print(f"Upregulated in {cond2} (padj < 0.01, log2FC > 1): {len(up_01)}")
    print(f"Downregulated in {cond2} (padj < 0.01, log2FC < -1): {len(down_01)}")
    print()
    
    # Save results for this comparison
    comparison_name = f"{cond2}_vs_{cond1}"
    results_sorted.to_csv(f"{pairwise_dir}/{comparison_name}_all_genes.csv", index=True)
    
    if len(sig_05) > 0:
        sig_05.to_csv(f"{pairwise_dir}/{comparison_name}_significant_padj05.csv", index=True)
    if len(sig_01) > 0:
        sig_01.to_csv(f"{pairwise_dir}/{comparison_name}_significant_padj01.csv", index=True)
    
    # Store summary
    all_comparisons.append({
        'comparison': comparison_name,
        'condition1': cond1,
        'condition2': cond2,
        'total_genes': len(results),
        'sig_padj05': len(sig_05),
        'sig_padj01': len(sig_01),
        'upregulated': len(up_01),
        'downregulated': len(down_01)
    })
    
    # Generate MA plot
    plt.figure(figsize=(10, 8))
    plt.scatter(results['baseMean'], results['log2FoldChange'], 
               c=['red' if (p < 0.05 and abs(fc) > 1) else 'gray' 
                  for p, fc in zip(results['padj'].fillna(1), results['log2FoldChange'])],
               alpha=0.5, s=10)
    plt.xlabel('Mean Expression (log scale)')
    plt.ylabel('Log2 Fold Change')
    plt.title(f'MA Plot: {cond2} vs {cond1}')
    plt.xscale('log')
    plt.axhline(y=0, color='black', linestyle='--', linewidth=0.5)
    plt.axhline(y=1, color='blue', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.axhline(y=-1, color='blue', linestyle='--', linewidth=0.5, alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"{pairwise_dir}/{comparison_name}_MA_plot.png", dpi=300)
    plt.close()
    
    # Generate Volcano plot
    plt.figure(figsize=(10, 8))
    results_plot = results.copy()
    results_plot['-log10(padj)'] = -np.log10(results_plot['padj'].fillna(1))
    
    plt.scatter(results_plot['log2FoldChange'], results_plot['-log10(padj)'],
               c=['red' if (p < 0.05 and abs(fc) > 1) else 'gray'
                  for p, fc in zip(results['padj'].fillna(1), results['log2FoldChange'])],
               alpha=0.5, s=10)
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(adjusted p-value)')
    plt.title(f'Volcano Plot: {cond2} vs {cond1}')
    plt.axhline(y=-np.log10(0.05), color='blue', linestyle='--', linewidth=0.5)
    plt.axvline(x=1, color='blue', linestyle='--', linewidth=0.5)
    plt.axvline(x=-1, color='blue', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.savefig(f"{pairwise_dir}/{comparison_name}_volcano_plot.png", dpi=300)
    plt.close()

print("All pairwise comparisons complete!")
print()

# #######
# GLOBAL FDR CORRECTION (Optional but recommended)
########

print("APPLYING GLOBAL FDR CORRECTION")

print()

from scipy.stats import false_discovery_control

print("Applying global FDR correction across all comparisons...")
print("This adjusts p-values considering ALL tests from ALL comparisons")
print()

# Apply global correction
all_pvalues_clean = [p for p in all_pvalues if not np.isnan(p)]
global_padj = false_discovery_control(all_pvalues_clean, method='bh')

print(f"Total tests across all comparisons: {len(all_pvalues_clean)}")
print(f"Tests with global padj < 0.05: {sum(global_padj < 0.05)}")
print(f"Tests with global padj < 0.01: {sum(global_padj < 0.01)}")
print()

# Note: To use global correction, you'd need to reassign these values back to each comparison
# This is left as an exercise as it requires careful bookkeeping

#########
# CREATE SUMMARY REPORT
############

# Convert to DataFrame
summary_df = pd.DataFrame(all_comparisons)
summary_df = summary_df.sort_values('sig_padj05', ascending=False)

print()
print("SUMMARY OF ALL PAIRWISE COMPARISONS")
print()
print(summary_df.to_string(index=False))
print()

# Save summary
summary_df.to_csv(f"{output_dir}/pairwise_comparisons_summary.csv", index=False)
print(f"Saved summary: {output_dir}/pairwise_comparisons_summary.csv")

# Create a heatmap of significant genes across comparisons
print()
print("Creating comparison heatmap...")

heatmap_data = summary_df.pivot_table(
    values='sig_padj05', 
    index='condition1', 
    columns='condition2',
    fill_value=0
)

plt.figure(figsize=(12, 10))
sns.heatmap(heatmap_data, annot=True, fmt='g', cmap='YlOrRd', cbar_kws={'label': 'Significant Genes (padj<0.05)'})
plt.title('Number of Significant Genes per Pairwise Comparison')
plt.tight_layout()
plt.savefig(f"{output_dir}/comparison_heatmap.png", dpi=300)
print(f"Saved: {output_dir}/comparison_heatmap.png")
plt.close()

###########
#GENERATE PCA PLOT (All samples)
###########


print()
print("Generating PCA plot for all samples...")
from sklearn.decomposition import PCA

# Get normalized counts
normalized_counts = dds.layers['normed_counts']

# Perform PCA
log_norm = np.log2(normalized_counts + 1)
pca = PCA(n_components=2)
pca_result = pca.fit_transform(log_norm)

plt.figure(figsize=(12, 10))

# Create color map for all conditions
color_map = plt.cm.get_cmap('tab10')
condition_colors = {cond: color_map(i) for i, cond in enumerate(unique_conditions)}
colors = [condition_colors[c] for c in metadata['condition']]

plt.scatter(pca_result[:, 0], pca_result[:, 1], c=colors, s=100, alpha=0.7)
for i, sample in enumerate(metadata['sample']):
    plt.annotate(sample, (pca_result[i, 0], pca_result[i, 1]), 
                fontsize=9, alpha=0.8)

plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
plt.title('PCA Plot - All Samples')

# Create legend
from matplotlib.patches import Patch
legend_elements = [Patch(facecolor=condition_colors[cond], label=cond) 
                   for cond in unique_conditions]
plt.legend(handles=legend_elements, loc='best')

plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{output_dir}/PCA_plot_all_samples.png", dpi=300)
print(f"Saved: {output_dir}/PCA_plot_all_samples.png")
plt.close()

print()

print("PAIRWISE ANALYSIS COMPLETE!")
print()
print(f"Results saved to: {pairwise_dir}/")
print(f"Summary saved to: {output_dir}/")
print()
print("Files generated:")
print(f"  - {n_comparisons} pairwise comparison results (CSV)")
print(f"  - {n_comparisons} MA plots (PNG)")
print(f"  - {n_comparisons} Volcano plots (PNG)")
print(f"  - Summary report (CSV)")
print(f"  - Comparison heatmap (PNG)")
print(f"  - PCA plot (PNG)")
print()
