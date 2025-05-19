# =========================================================
# Variant Calling Pipeline and SNP Analysis
# =========================================================

import os
import subprocess
from collections import Counter

# Ensure required tools (FastQC, Bowtie2, samtools, bcftools) are installed in the system
# (the script will call these via subprocess; modify paths if necessary)
# Define file paths for input data and reference
fastq1 = "subset_SRR099957_1.fastq.gz"   # Paired-end read 1 (subset of SRR099957 data)
fastq2 = "subset_SRR099957_2.fastq.gz"   # Paired-end read 2
reference_index = "/path/to/reference_index"  # Bowtie2 index basename (e.g., for hg38 or target genome)
reference_fasta = "/path/to/reference.fasta"  # Reference genome FASTA (needed for variant calling)


os.makedirs("fastqc", exist_ok=True)
os.makedirs("aligned", exist_ok=True)
os.makedirs("sorted", exist_ok=True)
os.makedirs("variants", exist_ok=True)

# Step 1: Quality control on raw reads using FastQC
# FastQC will produce HTML reports for each input FASTQ file in the "fastqc" directory
subprocess.run(f"fastqc {fastq1} {fastq2} -o fastqc", shell=True, check=True)

# Step 2: Align reads to reference genome using Bowtie2
# This produces a SAM file with alignments
aligned_sam = os.path.join("aligned", "subset_SRR099957.sam")
bowtie2_cmd = f"bowtie2 -x {reference_index} -1 {fastq1} -2 {fastq2} -S {aligned_sam}"
subprocess.run(bowtie2_cmd, shell=True, check=True)

# Step 3: Convert SAM to sorted BAM using samtools
sorted_bam = os.path.join("sorted", "subset_SRR099957.bam")
# The following command sorts the alignments and outputs a BAM file
subprocess.run(f"samtools sort -o {sorted_bam} {aligned_sam}", shell=True, check=True)

# Step 4: Index the BAM file (for fast random access) using samtools index
subprocess.run(f"samtools index {sorted_bam}", shell=True, check=True)

# Step 5: Call genetic variants using samtools mpileup and bcftools
# First, generate genotype likelihoods in BCF format from the aligned reads
raw_bcf = os.path.join("variants", "subset_SRR099957.bcf")
subprocess.run(f"samtools mpileup -g -f {reference_fasta} {sorted_bam} -o {raw_bcf}", shell=True, check=True)
# Next, call variants (SNPs and indels) from the BCF using bcftools
vcf_file = os.path.join("variants", "subset_SRR099957.vcf")
subprocess.run(f"bcftools call -mv {raw_bcf} -Ov -o {vcf_file}", shell=True, check=True)

# At this point, subset_SRR099957.vcf contains the variant calls in VCF format.

# Define a helper function to parse a VCF file (skipping header lines) into a pandas DataFrame
import pandas as pd
def parse_vcf(file_path):
    """Read a VCF file (skip header lines beginning with '#') into a DataFrame."""
    # VCF typically has 8 or more columns; define column names for standard fields 
    col_names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
    return pd.read_csv(file_path, sep='\t', comment='#', names=col_names)

# Function to count SNP types (transitions/transversions) in a VCF file
def count_snps(file_path):
    """Count single-nucleotide variant types (REF -> ALT substitutions) in the VCF file."""
    vcf_df = parse_vcf(file_path)
    # Filter to only single-nucleotide REF and ALT (exclude indels or multi-nucleotide variants)
    snv_mask = (vcf_df['REF'].str.len() == 1) & (vcf_df['ALT'].str.len() == 1)
    snv_df = vcf_df[snv_mask]
    # Use Counter to count occurrences of each (REF, ALT) pair
    snp_counts = Counter(zip(snv_df['REF'], snv_df['ALT']))
    return snp_counts

# Count SNP types in the called VCF
snp_counts = count_snps(vcf_file)
# Display the counts for each substitution type
print("SNP substitution counts (REF->ALT):")
for (ref_base, alt_base), count in snp_counts.items():
    print(f"{ref_base}->{alt_base}: {count}")

# Additionally, count total transitions vs transversions as a summary
transitions = 0
transversions = 0
for (ref_base, alt_base), count in snp_counts.items():
    pair = {ref_base, alt_base}
    if pair == {"A", "G"} or pair == {"C", "T"}:
        transitions += count
    else:
        transversions += count
print(f"Total transitions: {transitions}, Total transversions: {transversions}")
# Typically, transitions are more common than transversions in genomic variant data.

# Step 6: Analyze variant distribution across the genome
vcf_df = parse_vcf(vcf_file)
# Count variants per chromosome
variant_counts = vcf_df['CHROM'].value_counts()
print("\nVariant counts per chromosome:")
for chrom, count in variant_counts.items():
    print(f"{chrom}: {count}")
# (Further analysis could include examining whether variants are uniformly distributed or enriched in certain regions.)

# =========================================================
# Genomic Variant Distribution
# =========================================================

# Use the pandas DataFrame of variants (vcf_df) to analyze genomic distribution.
# For example, convert variant positions into a "GenomicRanges"-like representation.

variant_intervals = vcf_df[['CHROM', 'POS']].copy()
variant_intervals['start'] = variant_intervals['POS']
variant_intervals['end'] = variant_intervals['POS']
# The above DataFrame now holds intervals [start, end] for each variant (each SNP position).
print(f"\nTotal number of variant sites: {len(variant_intervals)}")
print("Example variant intervals:")
print(variant_intervals.head())
# We could further analyze where these variants occur (e.g., which chromosomes or genomic regions).
# (For instance, to find overlaps with gene regions, one would need gene annotation data; not shown here.)

# =========================================================
# ChIP-seq Data Processing (Peak Calling Pipeline)
# =========================================================

# we align ChIP-seq reads for a histone mark (H3K27me3) and prepare for peak calling.
fastq1_a3 = "H3K27me3_iPSC_SRA60_subset_1.fastq.gz"
fastq2_a3 = "H3K27me3_iPSC_SRA60_subset_2.fastq.gz"
# Create separate output directories 
os.makedirs("aligned_A3", exist_ok=True)
os.makedirs("sorted_A3", exist_ok=True)
# Align each FASTQ file (single-end reads) to the reference genome using Bowtie2
sam1 = os.path.join("aligned_A3", "H3K27me3_subset_1.sam")
sam2 = os.path.join("aligned_A3", "H3K27me3_subset_2.sam")
subprocess.run(f"bowtie2 -x {reference_index} -U {fastq1_a3} -S {sam1}", shell=True, check=True)
subprocess.run(f"bowtie2 -x {reference_index} -U {fastq2_a3} -S {sam2}", shell=True, check=True)
# Convert SAM to BAM and sort
bam1 = os.path.join("sorted_A3", "H3K27me3_subset_1.bam")
bam2 = os.path.join("sorted_A3", "H3K27me3_subset_2.bam")
subprocess.run(f"samtools sort -o {bam1} {sam1}", shell=True, check=True)
subprocess.run(f"samtools sort -o {bam2} {sam2}", shell=True, check=True)
# Index the BAM files
subprocess.run(f"samtools index {bam1}", shell=True, check=True)
subprocess.run(f"samtools index {bam2}", shell=True, check=True)
# At this point, the reads are aligned and ready for peak calling.
# Typically, one would use a peak-calling tool (e.g., MACS2) on these BAM files to identify enriched regions (peaks).
# (Peak calling commands would be run here, e.g., macs2 callpeak ..., but are omitted in this pure Python conversion.)

# =========================================================
# ATAC-seq Differential Accessibility Analysis
# =========================================================

# Load ATAC-seq count matrix (regions x samples)
atac_counts_file = "count_matrix_raw_atac_BRM014_ACBI1.csv.gz"
atac_df = pd.read_csv(atac_counts_file, sep=",")
# The first column is 'region' (peak identifier), remaining columns are read counts for each sample.
# Reshape the data to long format for easier manipulation
atac_long = atac_df.melt(id_vars='region', var_name='sample', value_name='read_count')
# Parse sample name into replicate, timepoint, and treatment (e.g., "R1_24h_control" -> replicate=1, timepoint=24, treatment="control")
sample_pattern = r'R(\d+)_(\d+)h_(.+)'
atac_long[['replicate', 'timepoint', 'treatment']] = atac_long['sample'].str.extract(sample_pattern)
atac_long['replicate'] = atac_long['replicate'].astype(int)
atac_long['timepoint'] = atac_long['timepoint'].astype(int)
# Calculate counts per million (CPM) for each sample to normalize for library size
atac_long['countsPerMillion'] = atac_long.groupby('sample')['read_count'].transform(lambda x: x / x.sum() * 1e6)

import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

def differential_expression(counts_matrix, group_labels, norm_method="TMM"):
    """
    Perform a basic differential analysis between two groups on the provided counts matrix.
    - counts_matrix: 2D numpy array of shape (n_regions, n_samples).
    - group_labels: list/array of length n_samples indicating the group for each sample.
    - norm_method: Normalization method to apply ("TMM" or "loess").
    Returns a DataFrame with log2 fold changes and p-values for each region.
    """
    groups = np.array(group_labels)
    unique_groups = np.unique(groups)
    assert len(unique_groups) == 2, "This function expects exactly two groups."
    group1, group2 = unique_groups[0], unique_groups[1]
    idx1 = np.where(groups == group1)[0]
    idx2 = np.where(groups == group2)[0]
    data = counts_matrix.astype(float).copy()
    # Apply normalization
    if norm_method == "TMM":
        # Approximate TMM by scaling each sample by its library size (column sum) relative to the mean library size
        lib_sizes = data.sum(axis=0)
        scale_factors = lib_sizes / lib_sizes.mean()
        data = data / scale_factors
    elif norm_method == "loess":
        # For simplicity, skip explicit loess normalization in this script
        pass
    # Compute log2 fold-change between group2 and group1 (add a small pseudocount to avoid log2 of zero)
    pseudo = 1e-6
    mean1 = data[:, idx1].mean(axis=1) + pseudo
    mean2 = data[:, idx2].mean(axis=1) + pseudo
    log2_fc = np.log2(mean2 / mean1)
    # Perform a two-sample t-test for each region
    p_values = []
    for i in range(data.shape[0]):
        vals1 = data[i, idx1]
        vals2 = data[i, idx2]
        # Welch's t-test (unequal variance) for two groups
        _, p_val = ttest_ind(vals1, vals2, equal_var=False)
        p_values.append(p_val)
    p_values = np.array(p_values)
    # Adjust p-values for multiple testing (Benjamini-Hochberg FDR)
    _, adj_p_values, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')
    # Compile results
    results = pd.DataFrame({
        'region': np.array(atac_df['region']),  # region identifiers
        'log2FC': log2_fc,
        'p_value': p_values,
        'adj_p_val': adj_p_values
    })
    results.sort_values('adj_p_val', inplace=True)
    return results

# Define comparisons (treatment pairs) and normalization methods to evaluate
comparisons = [("control", "BI-protac"), ("control", "BRM014"), ("BRM014", "BI-protac")]
normalization_methods = ["loess", "TMM"]

for comp in comparisons:
    # Filter data for the two treatments in this comparison
    mask = atac_long['treatment'].isin(comp)
    sub_df = atac_long[mask]
    # Pivot to get a matrix of CPM values (rows = regions, columns = samples)
    cpm_matrix = sub_df.pivot(index='region', columns='sample', values='countsPerMillion').fillna(0)
    # Determine group labels for each sample column (based on treatment in the sample name)
    sample_groups = [name.split('_')[-1] for name in cpm_matrix.columns]
    for norm in normalization_methods:
        de_results = differential_expression(cpm_matrix.values, sample_groups, norm_method=norm)
        # Count number of significant regions at FDR < 0.05
        sig_peaks = (de_results['adj_p_val'] < 0.05).sum()
        print(f"Comparison {comp[0]} vs {comp[1]} (Normalization: {norm}) -> Significant peaks: {sig_peaks}")

# =========================================================
# Machine Learning for Genomic Predictions
# =========================================================

# This uses DNA sequence data to predict gene expression (regression) using ML models.
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import tensorflow as tf
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Flatten, Conv1D, MaxPooling1D

# Load the dataset of sequences and expression values
sequence_data_file = "sequence_expression_data.txt"
seq_df = pd.read_csv(sequence_data_file, sep='\t', header=None, names=['sequence', 'expression'])
# One-hot encode DNA sequences (A, C, G, T) for model input
sequences = seq_df['sequence'].values
seq_length = len(sequences[0])
# Create one-hot encoding for each sequence: shape (seq_length, 4)
nuc_map = {'A': [1,0,0,0], 'C': [0,1,0,0], 'G': [0,0,1,0], 'T': [0,0,0,1]}
X_seq = np.array([[nuc_map[nuc] for nuc in seq] for seq in sequences])
y_expr = seq_df['expression'].astype(float).values
# Split into training and validation sets
X_train, X_val, y_train, y_val = train_test_split(X_seq, y_expr, test_size=0.2, random_state=42)

# Define a simple fully-connected neural network model
def create_simple_nn(input_shape):
    model = Sequential([
        Flatten(input_shape=input_shape),
        Dense(512, activation='relu'),
        Dense(1)
    ])
    return model

# Define a convolutional neural network for sequence data
def create_cnn(input_shape):
    model = Sequential([
        Conv1D(32, kernel_size=17, activation='relu', input_shape=input_shape),
        MaxPooling1D(pool_size=2),
        Conv1D(64, kernel_size=17, activation='relu'),
        MaxPooling1D(pool_size=2),
        Flatten(),
        Dense(64, activation='relu'),
        Dense(1)
    ])
    return model

# Train and evaluate the simple neural network
simple_nn_model = create_simple_nn((seq_length, 4))
simple_nn_model.compile(optimizer=tf.keras.optimizers.SGD(lr=0.1), loss='mse')
# (Training may be slow; we'll suppress verbose output)
simple_nn_model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=50, batch_size=1024, verbose=0)
y_pred_nn = simple_nn_model.predict(X_val)
mse_nn = mean_squared_error(y_val, y_pred_nn)
r2_nn = r2_score(y_val, y_pred_nn)
print(f"\nSimple NN model - MSE: {mse_nn:.4f}, R^2: {r2_nn:.4f}")

# Train and evaluate the convolutional neural network
cnn_model = create_cnn((seq_length, 4))
cnn_model.compile(optimizer=tf.keras.optimizers.Adam(lr=1e-4), loss='mse')
# Use a callback to save the best model based on validation loss
checkpoint_cb = tf.keras.callbacks.ModelCheckpoint('best_cnn_model.h5', save_best_only=True, monitor='val_loss')
cnn_model.fit(X_train, y_train, validation_data=(X_val, y_val), epochs=50, batch_size=1024, callbacks=[checkpoint_cb], verbose=0)
# Load the best weights and evaluate on validation set
cnn_model.load_weights('best_cnn_model.h5')
y_pred_cnn = cnn_model.predict(X_val)
mse_cnn = mean_squared_error(y_val, y_pred_cnn)
r2_cnn = r2_score(y_val, y_pred_cnn)
print(f"CNN model - MSE: {mse_cnn:.4f}, R^2: {r2_cnn:.4f}")
# (The CNN model, by capturing sequence motifs, is expected to perform similarly or better than the simple NN.)

# =========================================================
# GWAS and Polygenic Risk Scores (PRS)
# =========================================================

# In this, we compute a polygenic risk score manually and compare it to an automated pipeline output.
# Load GWAS summary statistics (SNP effect sizes from a GWAS)
effect_sizes_file = "Tapas_enjoyability_GWAS_sumStats.txt"
effect_df = pd.read_table(effect_sizes_file, header=0)
# Load genotype data (SNP dosages for each individual in a cohort)
genotypes_file = "MiniCohort_Tapas_SNPdosages.txt"
genotypes_df = pd.read_table(genotypes_file, header=0)
# If a family ID (FID) column is present, drop it (only need IID and SNP columns)
if 'FID' in genotypes_df.columns:
    genotypes_df.drop(columns=['FID'], inplace=True)
# Compute the Polygenic Risk Score for each individual: sum of (dosage * effect_size) for all SNPs
genotypes_df['PRS'] = 0.0
effect_map = dict(zip(effect_df['SNP'], effect_df['Effect_Size']))
for snp, beta in effect_map.items():
    if snp in genotypes_df.columns:
        genotypes_df['PRS'] += genotypes_df[snp] * beta
# Display PRS for the first few individuals
print("\nComputed Polygenic Risk Scores (PRS) for first 5 individuals:")
print(genotypes_df[['IID', 'PRS']].head())
# In practice, we would validate these PRS by correlating them with phenotypes or comparing with pipeline-generated PRS values.
