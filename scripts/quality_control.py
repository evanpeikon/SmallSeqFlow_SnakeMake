# Script to perform QC on RNA-seq data

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import json
import os

def visualize_rnaseq_qc(count_matrix, figure_size=(15, 12)):
    """
    Generate quality control visualizations for RNA-seq data.
    
    Parameters:
    -----------
    count_matrix : pandas.DataFrame
        Count matrix with gene symbols as first column
    figure_size : tuple, default=(15, 12)
        Size of the main figure
        
    Returns:
    --------
    tuple
        Main figure, dendrogram figure, and QC statistics
    """
    # Drop the Gene Name column for counting
    countlist_no_name = count_matrix.iloc[:, 1:]

    # Calculate total counts and log transform
    total_counts = countlist_no_name.sum(axis=0)
    log_counts = countlist_no_name.apply(lambda x: np.log2(x + 1))

    # Create main visualization figure
    fig1, axes = plt.subplots(2, 2, figsize=figure_size)

    # Panel 1: Total counts per sample
    sns.barplot(x=countlist_no_name.columns, y=total_counts, color='skyblue', ax=axes[0,0])
    axes[0,0].set_ylabel('Total Counts')
    axes[0,0].set_title('Total Counts per Sample')
    axes[0,0].tick_params(axis='x', rotation=85)

    # Panel 2: Log transformed counts distribution
    log_counts.boxplot(ax=axes[0,1])
    axes[0,1].set_ylabel('Log2(Counts + 1)')
    axes[0,1].set_title('Log Transformed Counts per Sample')
    axes[0,1].tick_params(axis='x', rotation=85)

    # Panel 3: Sample correlation heatmap
    correlation_matrix = log_counts.corr()
    sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', center=0.5, vmin=0, vmax=1, ax=axes[1,0])
    axes[1,0].set_title('Sample Correlation Matrix')

    # Panel 4: PCA plot
    pca = PCA(n_components=2)
    scaler = StandardScaler()
    pca_result = pca.fit_transform(scaler.fit_transform(log_counts.T))
    pca_df = pd.DataFrame(data=pca_result, columns=['PC1', 'PC2'], index=log_counts.columns)
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', s=100, ax=axes[1,1])
    for idx, row in pca_df.iterrows():
        axes[1,1].annotate(idx, (row['PC1'], row['PC2']))
    axes[1,1].set_title(f'PCA Plot\nPC1 ({pca.explained_variance_ratio_[0]:.1%}) vs PC2 ({pca.explained_variance_ratio_[1]:.1%})')
    plt.tight_layout()

    # Create dendrogram figure
    fig2 = plt.figure(figsize=(8, 6))
    h_clustering = linkage(log_counts.T, 'ward')
    dendrogram(h_clustering, labels=countlist_no_name.columns)
    plt.xticks(rotation=90)
    plt.ylabel('Distance')
    plt.title('Sample Clustering Dendrogram')
    plt.tight_layout()

    # Generate QC metrics
    qc_stats = {
        'total_reads': int(total_counts.sum()),
        'mean_reads_per_sample': float(total_counts.mean()),
        'cv_reads': float(total_counts.std() / total_counts.mean()),
        'min_sample_correlation': float(correlation_matrix.min().min()),
        'max_sample_correlation': float(correlation_matrix.max().min()),
        'pc1_variance': float(pca.explained_variance_ratio_[0]),
        'pc2_variance': float(pca.explained_variance_ratio_[1])
    }
    
    print("\nRNA-seq Quality Control Metrics:")
    print(f"Total sequencing depth: {qc_stats['total_reads']:,.0f}")
    print(f"Mean reads per sample: {qc_stats['mean_reads_per_sample']:,.0f}")
    return fig1, fig2, qc_stats

# Get input and output files from Snakemake
input_file = snakemake.input.gene_counts
output_main_fig = snakemake.output.main_fig
output_dendrogram_fig = snakemake.output.dendrogram_fig
output_stats = snakemake.output.stats

# Create output directories if they don't exist
os.makedirs(os.path.dirname(output_main_fig), exist_ok=True)

# Load the count matrix
count_matrix = pd.read_csv(input_file, index_col=0)

# Perform QC
main_fig, dendrogram_fig, qc_stats = visualize_rnaseq_qc(count_matrix)

# Save figures and stats
main_fig.savefig(output_main_fig, dpi=300, bbox_inches='tight')
dendrogram_fig.savefig(output_dendrogram_fig, dpi=300, bbox_inches='tight')

# Save stats as JSON
with open(output_stats, 'w') as f:
    json.dump(qc_stats, f, indent=4)
