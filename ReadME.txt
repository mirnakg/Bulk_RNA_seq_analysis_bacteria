This is an RNA-sequencing analysis pipeline to obtain differential gene expression analysis from paired-end fastq input. 


Files you need to get started: 

Reference genome (fasta file)
Genome annotations file (GTF file)
Compressed fastq files
A sample_names.txt file  containing a list of the sample names as they appear on the fastq file name. The pipeline is written for paired-end sequencing and will take into account that there are two fastq files for each sample. 

Prepping your environment: 

The following tools need to be installed, it is recommended that you create a venv to install and run everything 

python
trimmomatic 
subread
bowtie2
PyDESeq2
seaborn
pandas
matplotlib
numpy
scipy
scikit-learn

5 shell scripts are included to run the following: 
Run sbatch genome_idx_bowtie2.sh to generate a genome index 
Run sbatch trim_read.sh to remove adapters and low quality reads, tailored for illumina library prep 
Run sbatch align_trimmed_reads.sh to align trimmed reads and generate bam files
Run sbatch get_raw_counts.sh  to get raw gene counts using FeatureCounts
Run python deseq2_analysis.py (recommended to minimize python package issues) or sbatch run_deseq2.sh to obtain a differential gene expression analysis using pyDEseq2 


Editing shell scripts

In the shell scripts, update directories as well as the email address where job status would be sent. 
The trimming and alignment scripts run as parallel jobs, update the --array flag to the number of samples you have (For example, if you have 6 samples, then --array=1-6)
In the python script, deseq performs differential gene expression analysis on two groups at a time. Edit this script to indicate the group of each of your samples in the same order they appear in the sample_names.txt file . 

