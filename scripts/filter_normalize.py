# Script to filter and normalize RNA-seq data

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import json
import os

def plot_genes_retained_by_cpm(data, min_samples=2, output_file=None):
    """
    Plot the number of genes retained as a function of different CPM thresholds.
    
    Parameters:
    -----------
    data : pandas.DataFrame
        Count matrix (without gene name column)
    min_samples : int, default=2
        Minimum number of samples that must exceed the CPM threshold
    output_file : str or None
        Path to save the plot
        
    Returns:
    --------
    matplotlib.figure.Figure
        Figure object
    """
    # Convert raw counts to CPM to normalize the data
    cpm = data.apply(lambda x: (x / x.sum()) * 1e6)
    
    # Define a range of CPM thresholds to test, from 0 to 5 with increments of 0.1
    thresholds = np.arange(0, 5, 0.1)
    
    # Initialize list to store the # of genes retained for each threshold
    genes_retained = []

    # Loop through each threshold value to determine the # of genes retained
    for min_cpm in thresholds:
        # Create mask where CPM > min_cpm in at least min_samples samples
        mask = (cpm > min_cpm).sum(axis=1) >= min_samples
        # Count # of genes that meet the criteria and append to the list
        genes_retained.append(mask.sum())

    # Plot # of genes retained as a function of CPM threshold
    fig = plt.figure(figsize=(10, 6))
    plt.plot(thresholds, genes_retained, marker='o', color='green')
    plt.axvline(x=1.0, color='red', linestyle='--', label='CPM = 1')
    plt.xlabel('Threshold (CPM)')
    plt.ylabel('Number of Genes Retained')
    plt.title('Gene Retention vs. CPM Threshold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    
    return fig

def filter_normalize(data, min_cpm=1.0, min_samples=2):
    """
    Filter low-expressed genes and normalize the count matrix.
    
    Parameters:
    -----------
    data : pandas.DataFrame
        Count matrix with gene symbols as first column
    min_cpm : float, default=1.0
        Minimum counts per million threshold
    min_samples : int, default=2
        Minimum number of samples that must exceed the CPM threshold
        
    Returns:
    --------
    tuple
        Normalized data and diagnostic statistics
    """
    # Extract structural components
    gene_names = data.iloc[:, 0]
    raw_counts = data.iloc[:, 1:]

    # Implement DESeq2-style filtering
    lib_sizes = raw_counts.sum(axis=0)
    cpm = raw_counts.div(lib_sizes, axis=1) * 1e6
    mask = (cpm > min_cpm).sum(axis=1) >= min_samples

    # Apply filtration criteria
    filtered_counts = raw_counts[mask]
    filtered_gene_names = gene_names[mask]

    # Calculate geometric means with DESeq2-inspired approach
    log_counts = np.log(filtered_counts.replace(0, np.nan))
    geometric_means = np.exp(log_counts.mean(axis=1))

    # Estimate size factors using DESeq2 methodology
    size_factor_ratios = filtered_counts.div(geometric_means, axis=0)
    size_factors = size_factor_ratios.median(axis=0)

    # Apply normalization transformation
    normalized_counts = filtered_counts.div(size_factors, axis=1)

    # Reconstruct data architecture
    normalized_data = pd.DataFrame({
        'Gene_Name': filtered_gene_names
    })
    normalized_data = pd.concat([normalized_data, normalized_counts], axis=1)

    # Generate diagnostic metrics
    diagnostics = {
        'total_genes_initial': int(len(data)),
        'genes_post_filtering': int(len(normalized_data)),
        'size_factors': {str(k): float(v) for k, v in size_factors.items()},
        'mean_size_factor': float(size_factors.mean()),
        'size_factor_variance': float(size_factors.var())
    }

    return normalized_data, diagnostics

# Get input, output and parameters from Snakemake
input_file = snakemake.input.gene_counts
output_file = snakemake.output.filtered_counts
output_stats = snakemake.output.stats
output_cpm_plot = snakemake.output.cpm_plot
min_cpm = snakemake.params.min_cpm
min_samples = snakemake.params.min_samples

# Create output directories if they don't exist
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Load the count matrix
count_matrix = pd.read_csv(input_file, index_col=0)

# Generate the CPM plot
countlist_no_name = count_matrix.iloc[:, 1:]  # Exclude gene_name column
plot_fig = plot_genes_retained_by_cpm(countlist_no_name, min_samples=min_samples, output_file=output_cpm_plot)
plt.close(plot_fig)

# Filter and normalize the data
filtered_normalized_count_matrix, stats = filter_normalize(
    count_matrix, 
    min_cpm=min_cpm, 
    min_samples=min_samples
)

# Save the filtered and normalized counts
filtered_normalized_count_matrix.to_csv(output_file)

# Save stats as JSON
with open(output_stats, 'w') as f:
    json.dump(stats, f, indent=4)

print(f"Filtering complete: {stats['genes_post_filtering']} genes retained out of {stats['total_genes_initial']}")
