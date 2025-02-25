# Script to visualize differential expression results

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os

def visualize_differential_expression_matrix(results_df, filtered_degs, expression_matrix, treatment_columns, control_columns, p_adj_threshold=0.05, abs_log2fc_threshold=1.0, figure_size=(10, 8)):
    """
    Visualize differential expression analysis results.
    
    Parameters:
    -----------
    results_df : pandas.DataFrame
        Complete results from differential expression analysis
    filtered_degs : pandas.DataFrame
        Subset of significant DEGs
    expression_matrix : pandas.DataFrame
        Normalized count matrix
    treatment_columns : list
        List of treatment sample column names
    control_columns : list
        List of control sample column names
    p_adj_threshold : float, default=0.05
        Significance threshold for adjusted p-values
    abs_log2fc_threshold : float, default=1.0
        Log2 fold change threshold for significance
    figure_size : tuple, default=(10, 8)
        Size of the figure
        
    Returns:
    --------
    tuple
        Figure object and summary statistics
    """
    fig, axes = plt.subplots(2, 2, figsize=figure_size)
    scatter_params = {'alpha': 0.8, 'edgecolor': None, 'palette': 'viridis'}

    # Panel 1: Global Expression Landscape (Volcano Plot)
    sns.scatterplot(data=results_df, x='log2fc', y=-np.log10(results_df['p_adj']), hue='abs_log2fc', ax=axes[0,0], **scatter_params)
    axes[0,0].axhline(y=-np.log10(p_adj_threshold), color='red', linestyle='--', linewidth=1)
    axes[0,0].axvline(x=abs_log2fc_threshold, color='blue', linestyle='--', linewidth=1)
    axes[0,0].axvline(x=-abs_log2fc_threshold, color='blue', linestyle='--', linewidth=1)
    axes[0,0].set_xlabel('log2 Fold Change')
    axes[0,0].set_ylabel('-log10(Adjusted P-value)')
    axes[0,0].set_title('Volcano Plot')

    # Panel 2: Fold Change Distribution (All Genes)
    sns.histplot(data=results_df, x='abs_log2fc', bins=50, kde=True, ax=axes[0,1])

    # Add vertical line at fold change threshold
    axes[0,1].axvline(x=abs_log2fc_threshold, color='red', linestyle='--', linewidth=1)

    axes[0,1].set_title('Distribution of Absolute log2FC (All Genes)')
    axes[0,1].set_xlabel('Absolute log2 Fold Change')
    axes[0,1].set_ylabel('Gene Frequency')

    # Panel 3: MA Plot
    results_df['mean_expression'] = np.log2((results_df['mean_treated'] + results_df['mean_control'])/2 + 1)

    sns.scatterplot(data=results_df, x='mean_expression', y='log2fc', hue='significant' if 'significant' in results_df.columns else None, ax=axes[1,0], **scatter_params)
    axes[1,0].axhline(y=0, color='red', linestyle='--', linewidth=1)
    axes[1,0].set_title('MA Plot (Mean vs Fold Change)')
    axes[1,0].set_xlabel('Mean Expression (log2)')
    axes[1,0].set_ylabel('log2 Fold Change')

    # Panel 4: Distribution of Adjusted P-values
    sns.histplot(data=results_df, x='p_adj', bins=50, kde=True, ax=axes[1,1])

    # Add vertical line at significance threshold
    axes[1,1].axvline(x=p_adj_threshold, color='red', linestyle='--', linewidth=1)
    axes[1,1].set_title('Distribution of Adjusted P-values')
    axes[1,1].set_xlabel('Adjusted P-value')
    axes[1,1].set_ylabel('Gene Frequency')

    plt.tight_layout()

    # Generate comprehensive analytical metrics
    summary_stats = {
        'total_genes': int(len(results_df)),
        'significant_genes': int(len(filtered_degs)),
        'mean_fold_change_all': float(results_df['abs_log2fc'].mean()),
        'median_fold_change_all': float(results_df['abs_log2fc'].median()),
        'max_fold_change': float(results_df['abs_log2fc'].max()),
        'mean_fold_change_sig': float(filtered_degs['abs_log2fc'].mean()) if len(filtered_degs) > 0 else 0,
        'median_padj': float(results_df['p_adj'].median()),
        'genes_below_alpha': int(sum(results_df['p_adj'] < p_adj_threshold))
    }

    print("\nComprehensive Expression Analysis Metrics:")
    print(f"Total genes analyzed: {summary_stats['total_genes']}")
    print(f"Significant DEGs identified: {summary_stats['significant_genes']}")
    print(f"Mean absolute log2FC (all genes): {summary_stats['mean_fold_change_all']:.2f}")
    if len(filtered_degs) > 0:
        print(f"Mean absolute log2FC (significant): {summary_stats['mean_fold_change_sig']:.2f}")
    print(f"Median adjusted p-value: {summary_stats['median_padj']:.3f}")
    print(f"Genes below significance threshold: {summary_stats['genes_below_alpha']}")
    
    return fig, summary_stats

# Get input, output and parameters from Snakemake
input_deg_results = snakemake.input.deg_results
input_filtered_counts = snakemake.input.filtered_counts
output_plot = snakemake.output.de_plot
treatment_columns = snakemake.params.treatment_columns
control_columns = snakemake.params.control_columns
p_adj_threshold = snakemake.params.p_adj_threshold
abs_log2fc_threshold = snakemake.params.abs_log2fc_threshold

# Create output directories if they don't exist
os.makedirs(os.path.dirname(output_plot), exist_ok=True)

# Load the differential expression results and filtered counts
results_df = pd.read_csv(input_deg_results, index_col=0)
filtered_normalized_count_matrix = pd.read_csv(input_filtered_counts, index_col=0)

# Extract significant DEGs
filtered_degs = results_df[results_df['significant'] == True]

# Visualize the differential expression results
fig, summary_stats = visualize_differential_expression_matrix(
    results_df=results_df,
    filtered_degs=filtered_degs,
    expression_matrix=filtered_normalized_count_matrix,
    treatment_columns=treatment_columns,
    control_columns=control_columns,
    p_adj_threshold=p_adj_threshold,
    abs_log2fc_threshold=abs_log2fc_threshold
)

# Save the plot
fig.savefig(output_plot, dpi=300, bbox_inches='tight')

print(f"Visualization complete: Plot saved to {output_plot}")
