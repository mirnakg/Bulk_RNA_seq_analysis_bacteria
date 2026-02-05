# Bulk_RNA_seq_analysis_bacteria
This pipeline performs a complete RNA-seq differential expression analysis for bacterial samples. The workflow includes read trimming, alignment, quantification, and differential expression analysis using DESeq2.

## Pipeline Components

The analysis consists of five main scripts that should be run sequentially:

1. **genome_idx_bowtie2.sh** - Build genome index
2. **trim_read.sh** - Quality control and adapter trimming
3. **align_timmed_reads.sh** - Read alignment to reference genome
4. **get_raw_counts.sh** - Quantify gene expression
5. **deseq2_analysis.py** - Differential expression analysis

---

## Prerequisites

### Software Requirements

- **Conda/Miniconda** with your virtual environment (Here called SRA-tools)
- **Bowtie2** (v2.3.0 or later) - Fast and memory-efficient aligner for bacterial genomes
- **Trimmomatic** - Read quality trimming and adapter removal
- **SAMtools** (v1.10 or later) - SAM/BAM file manipulation
- **featureCounts** (from Subread package) - Gene expression quantification
- **Python 3** with the following packages:
  - pandas
  - numpy
  - pydeseq2
  - matplotlib
  - seaborn
  - scipy
  - scikit-learn

### Input Files Required (Below is a public dataset that can be used with this pipeline)

1. **Reference genome FASTA**: `GCF_000009045.1_ASM904v1_genomic.fna`
   - *Bacillus subtilis* strain 168 reference genome
   - Download from NCBI: https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=224308

2. **Genome annotation GTF**: `GCF_000009045.1_ASM904v1_genomic.gtf`
   - Gene annotations for feature counting
   - Must match the reference genome version

3. **Raw FASTQ files**: Paired-end sequencing reads
   - Format: `{sample}_1.fastq.gz` and `{sample}_2.fastq.gz`
   - Should be listed in `SRR_Acc_List.txt` (one sample per line)
   - Use this public dataset for pipeline testing: SRP422250: Time-Course Transcriptome Analysis of Bacillus subtilis DB104

---

## Step-by-Step Pipeline Execution

### Step 1: Build Bowtie2 Genome Index

**Script**: `genome_idx_bowtie2.sh`

**Purpose**: Creates a Bowtie2 index from the reference genome, required for efficient alignment.

**SLURM Parameters**:
- **--cpus-per-task=16**: Uses 16 CPU cores for parallel index building
- **--mem=32G**: Allocates 32GB RAM (bacterial genomes are small, this is conservative)

**Bowtie2-build Parameters**:
- **--threads 16**: Parallel processing using 16 threads
- **GENOME_FA**: Input reference genome in FASTA format
- **INDEX_PREFIX**: Output directory and prefix for index files

**Expected Output**:
- Six index files with `.bt2` extension
- Files: `reference_genome.1.bt2`, `reference_genome.2.bt2`, etc.

**Run Command**:
```bash
sbatch genome_idx_bowtie2.sh
```

**Runtime**: ~5-15 minutes for bacterial genome

---

### Step 2: Quality Trimming with Trimmomatic

**Script**: `trim_read.sh`

**Purpose**: Removes adapter sequences and low-quality bases from raw reads.

**SLURM Parameters**:
- **--array=1-18**: Processes 18 samples in parallel (adjust based on your sample count)
- **--cpus-per-task=8**: 8 cores per sample
- **--mem=32G**: 32GB RAM per sample

**Trimmomatic PE Parameters**:

- **-threads 8**: Uses 8 parallel threads
- **-phred33**: Input quality scores are Phred+33 encoded (Illumina 1.8+)

**Trimming Steps** (applied in order):

1. **ILLUMINACLIP:TruSeq3-PE.fa:2:30:10**
   - Removes Illumina TruSeq3 adapters
   - `2`: Seed mismatches allowed
   - `30`: Palindrome clip threshold (for paired-end adapter detection)
   - `10`: Simple clip threshold

2. **LEADING:3**
   - Removes leading bases with quality < 3

3. **TRAILING:3**
   - Removes trailing bases with quality < 3

4. **SLIDINGWINDOW:4:15**
   - Scans reads with a 4-base sliding window
   - Cuts when average quality drops below 15

5. **MINLEN:36**
   - Discards reads shorter than 36 bases after trimming

**Expected Output** (per sample):
- `{sample}_R1_paired.fastq.gz` - Forward reads that survived with their pair
- `{sample}_R2_paired.fastq.gz` - Reverse reads that survived with their pair
- `{sample}_R1_unpaired.fastq.gz` - Forward reads without surviving pair
- `{sample}_R2_unpaired.fastq.gz` - Reverse reads without surviving pair
- `{sample}_trimmomatic.log` - Trimming statistics

**Run Command**:
```bash
sbatch trim_read.sh
```

**Runtime**: ~10-30 minutes per sample (depends on file size)

---

### Step 3: Read Alignment with Bowtie2

**Script**: `align_timmed_reads.sh`

**Purpose**: Aligns trimmed reads to the reference genome and produces sorted, indexed BAM files.

**SLURM Parameters**:
- **--array=1-18**: Processes 18 samples in parallel
- **--cpus-per-task=8**: 8 cores per alignment
- **--mem=64G**: 64GB RAM (conservative for bacterial genomes)

**Bowtie2 Alignment Parameters**:

- **-x ${GENOME_INDEX}**: Path to Bowtie2 index prefix
- **-1**: Forward paired reads (R1)
- **-2**: Reverse paired reads (R2)
- **-S**: Output SAM file
- **--threads 8**: Use 8 parallel threads
- **--very-sensitive**: Preset for highly sensitive alignment
  - Equivalent to: `-D 20 -R 3 -N 0 -L 20 -i S,1,0.50`
  - `-D 20`: Max consecutive seed extension failures
  - `-R 3`: Max re-seeding attempts
  - `-N 0`: No mismatches in seed
  - `-L 20`: Seed length
  - `-i S,1,0.50`: Seeding interval function

**SAMtools Processing**:

1. **samtools view** - SAM to BAM conversion
   - **-@ 8**: Use 8 threads
   - **-b**: Output in BAM format
   - **-o**: Output file

2. **samtools sort** - Sort BAM by coordinates
   - **-@ 8**: Use 8 threads
   - **-o**: Output sorted BAM
   - Required for downstream tools like featureCounts

3. **samtools index** - Create BAM index
   - **-@ 8**: Use 8 threads
   - Creates `.bai` index file
   - Required for visualization and some tools

4. **samtools flagstat** - Alignment statistics
   - Reports mapping rates, properly paired reads, etc.

**Expected Output** (per sample):
- `{sample}_Aligned.sortedByCoord.out.bam` - Sorted BAM file
- `{sample}_Aligned.sortedByCoord.out.bam.bai` - BAM index
- `{sample}_bowtie2.log` - Alignment statistics
- `{sample}_flagstat.txt` - Detailed alignment metrics

**Run Command**:
```bash
sbatch align_timmed_reads.sh
```

**Runtime**: ~20-60 minutes per sample

**Expected Alignment Rates**:
- Good quality bacterial RNA-seq: >85% alignment rate
- If <70%: Check for contamination or wrong reference genome

---

### Step 4: Gene Expression Quantification

**Script**: `get_raw_counts.sh`

**Purpose**: Counts reads mapping to each gene using the genome annotation.

**SLURM Parameters**:
- **--cpus-per-task=8**: 8 cores for parallel counting
- **--mem=32G**: 32GB RAM

**featureCounts Parameters**:

**Three counting strategies are performed**:

#### 1. Primary Count - CDS Features (Recommended)
Counts only reads mapping to coding sequences (protein-coding genes).

- **-T 8**: Use 8 threads
- **-p**: Paired-end mode (fragments counted, not individual reads)
- **-B**: Count only properly paired fragments
- **-C**: Don't count chimeric fragments (pairs mapping to different chromosomes)
- **-t CDS**: Count reads overlapping CDS features
- **-g gene_id**: Group counts by gene_id attribute
- **-s 0**: Unstranded library (use `-s 1` for forward-stranded, `-s 2` for reverse-stranded)
- **--minOverlap 10**: Minimum overlap between read and feature (10 bp)
- **--fracOverlap 0.2**: Minimum fraction of read overlapping feature (20%)
- **-Q 10**: Minimum mapping quality score (removes low-quality alignments)
- **--primary**: Count only primary alignments (exclude secondary/supplementary)

**Output**: `raw_counts_CDS.txt`

#### 2. Alternative Count - Gene Features
Counts reads across entire gene locus (may include UTRs if annotated).

- Same parameters as above except:
- **-t gene**: Count reads overlapping gene features

**Output**: `raw_counts_gene.txt`

#### 3. Diagnostic Count - With Multi-mappers
Allows multi-mapping reads for diagnostic purposes (not recommended for DE analysis).

- Additional parameters:
- **-M**: Count multi-mapping reads
- **--fraction**: Fractionally count multi-mappers (each read weighted by 1/n, where n = number of mapping locations)

**Output**: `raw_counts_CDS_with_multimappers.txt`

**Expected Output**:
- Three count matrices with genes as rows and samples as columns
- `.summary` files with assignment statistics for each run

**Important Statistics to Check**:
- **Assigned**: Should be >60% for good RNA-seq
- **Unassigned_NoFeatures**: Reads in intergenic regions
- **Unassigned_Ambiguity**: Reads overlapping multiple genes
- If "Assigned" is very low (<40%), check:
  - GTF annotation matches genome
  - Strandedness is correct
  - Library type matches parameters

**Run Command**:
```bash
sbatch get_raw_counts.sh
```

**Runtime**: ~10-30 minutes (processes all samples together)

---

### Step 5: Differential Expression Analysis with DESeq2

**Script**: `deseq2_analysis.py`

**Purpose**: Identifies differentially expressed genes between all pairwise condition comparisons.

**Important**: Before running, verify the `conditions` variable on line 62 matches your experimental design!

**Current Configuration**:
```python
conditions = ['15h', '12h', '12h', '12h', '10h', '10h', '10h', 
              '8h', '24h', '24h', '24h', '18h', '18h', '18h', 
              '15h', '15h', '8h', '8h']
```

**How to Modify for Your Data**:
- The list must be in the same order as samples appear in the count matrix
- Check sample order by running: `head -1 raw_counts_gene.txt`
- Each entry should match the biological condition/timepoint
- Replicates should have identical condition labels

**Analysis Steps**:

#### 1. Data Loading and Preparation
- Reads count matrix from featureCounts output
- Removes metadata columns (first 6 columns)
- Cleans sample names

#### 2. Gene Filtering
- **min_counts = 10**: Minimum read count threshold
- **min_samples = 3**: Gene must have ≥10 counts in at least 3 samples
- Removes lowly expressed genes to reduce noise and multiple testing burden

#### 3. DESeq2 Normalization and Analysis
**DESeqDataSet Parameters**:
- **counts**: Filtered count matrix (samples × genes)
- **metadata**: Sample information with condition labels
- **design_factors="condition"**: Model formula
- **refit_cooks=True**: Refit outliers (recommended)
- **n_cpus=8**: Use 8 CPU cores

**Statistical Method**:
- Uses negative binomial generalized linear model
- Estimates size factors for normalization
- Estimates gene-wise dispersions
- Fits model and performs Wald test

#### 4. Pairwise Comparisons
For each pair of conditions:
- **DeseqStats Parameters**:
  - **contrast=["condition", cond2, cond1]**: Compares cond2 vs cond1
  - Calculates log2 fold changes
  - Computes p-values and adjusts for multiple testing (Benjamini-Hochberg FDR)

**Significance Thresholds**:
- **padj < 0.05**: Standard significance level
- **|log2FoldChange| > 1**: 2-fold change cutoff (biological relevance)

#### 5. Global FDR Correction (Optional, adjust code )
- Adjusts p-values across ALL comparisons using Benjamini-Hochberg
- More conservative than per-comparison FDR
- Recommended for multiple pairwise tests

#### 6. Figures Generated

**Per Comparison**:
1. **MA Plot**: Mean expression vs. log2 fold change
   - Red points: Significant genes (padj < 0.05, |log2FC| > 1)
   - Gray points: Non-significant
   - Blue dashed lines: log2FC = ±1 thresholds

2. **Volcano Plot**: Log2 fold change vs. -log10(adjusted p-value)
   - Red points: Significant genes
   - Blue lines: Significance thresholds

**Overall Analysis**:
3. **Comparison Heatmap**: Number of significant genes per pairwise comparison
4. **PCA Plot**: Sample clustering in PC1-PC2 space
   - Colors indicate conditions
   - Shows sample variability and potential outliers

**Expected Output Files**:

**Per Comparison** (in `pairwise_comparisons/` directory):
- `{cond2}_vs_{cond1}_all_genes.csv` - Complete results for all genes
- `{cond2}_vs_{cond1}_significant_padj05.csv` - Genes with padj < 0.05
- `{cond2}_vs_{cond1}_significant_padj01.csv` - Genes with padj < 0.01
- `{cond2}_vs_{cond1}_MA_plot.png`
- `{cond2}_vs_{cond1}_volcano_plot.png`

**Summary Files**:
- `pairwise_comparisons_summary.csv` - Statistics for all comparisons
- `comparison_heatmap.png` - Visual summary
- `PCA_plot_all_samples.png` - Sample clustering

**Result File Columns**:
- **gene_id**: Gene identifier
- **baseMean**: Mean normalized counts across all samples
- **log2FoldChange**: Log2(cond2/cond1)
- **lfcSE**: Standard error of log2FC
- **stat**: Wald statistic
- **pvalue**: Raw p-value
- **padj**: Benjamini-Hochberg adjusted p-value (FDR)
- **comparison**: Which comparison this result is from

**Run Command**:
```bash
python deseq2_analysis.py
```

**Runtime**: ~5-30 minutes depending on gene count and number of comparisons

---

## Interpreting Results

### Log2 Fold Change Interpretation
- **log2FC = 1**: Gene is 2× higher in condition 2
- **log2FC = 2**: Gene is 4× higher in condition 2
- **log2FC = -1**: Gene is 2× lower in condition 2 (2× higher in condition 1)

### Adjusted P-values (padj)
- **padj < 0.05**: Traditionally significant (5% FDR)
- **padj < 0.01**: Highly significant (1% FDR)
- FDR = Expected proportion of false positives among significant genes

### Quality Control Checkpoints

**After Trimming**:
- Check trimmomatic logs for >70% reads surviving

**After Alignment**:
- Check bowtie2 logs for >80% alignment rate
- Check flagstat for >85% properly paired reads

**After Counting**:
- Check featureCounts summary for >60% assigned reads
- If low assignment: verify GTF matches genome, check strandedness

**After DESeq2**:
- Check PCA plot: replicates should cluster together
- If replicates are scattered: investigate sample quality or swap

---

## Troubleshooting

### Common Issues

**Problem**: Low alignment rate (<70%)
- **Solution**: Check reference genome matches organism, verify FASTQ quality

**Problem**: Low feature assignment (<40%)
- **Solution**: 
  - Verify GTF annotation matches genome version
  - Try different `-t` parameter (gene vs CDS)
  - Check strandedness with RSeQC or similar tool

**Problem**: No replicates clustering in PCA
- **Solution**: Check sample labels, verify no batch effects, check for swaps

**Problem**: Too many/few significant genes
- **Solution**: Verify experimental design, check for outliers, adjust fold change cutoff

**Problem**: Memory errors
- **Solution**: Increase SLURM `--mem` parameter or reduce `--cpus-per-task`

---

## File Organization

Recommended directory structure:
```
project/
├── fastq/                          # Raw FASTQ files
├── trimmed/                        # Trimmed FASTQ files
├── genome_index_bowtie2/           # Bowtie2 index files
├── bam/                            # Aligned BAM files
├── counts/                         # Raw count matrices
├── deseq2_results/                 # DESeq2 output
│   ├── pairwise_comparisons/       # Individual comparison results
│   ├── pairwise_comparisons_summary.csv
│   ├── comparison_heatmap.png
│   └── PCA_plot_all_samples.png
├── logs/                           # SLURM log files
└── scripts/                        # Analysis scripts
```

---

## Additional Notes

### Bacterial RNA-seq Specifics
- Bacterial genomes are small (~4 Mb for *B. subtilis*)
- No splicing considerations (unlike eukaryotes)
- Most reads should map to CDS regions


### Choosing Between CDS and Gene Counting
- **CDS (recommended)**: More specific, protein-coding genes only
- **Gene**: May include regulatory regions if annotated
- Compare both outputs if uncertain


### Memory and Runtime Optimization
- Bacterial genomes are small; most parameters are conservative
- Can reduce memory allocation if needed
- Array jobs provide automatic parallelization

---

## References

**Tools**:
- Bowtie2: Langmead & Salzberg (2012) Nature Methods
- Trimmomatic: Bolger et al. (2014) Bioinformatics
- featureCounts: Liao et al. (2014) Bioinformatics
- DESeq2: Love et al. (2014) Genome Biology

**Study**:
- SRP422250: Time-Course Transcriptome Analysis of *Bacillus subtilis* DB104

---

## Contact

For issues with this pipeline, check:
1. SLURM error logs in `logs/` directory
2. Tool-specific log files in output directories
3. Console output for error messages

Email: mirna.kheir.gouda at gmail (dot) com
