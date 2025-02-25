# Script to perform differential expression analysis

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import json
import os

def analyze_differential_expression(expression_matrix, treatment_columns, control_columns, alpha=0.05, lfc_threshold=1.0):
    """
    Perform differential expression analysis.
    
    Parameters:
    -----------
    expression_matrix : pandas.DataFrame
        Normalized count matrix with Gene_Name as first column
    treatment_columns : list
        List of treatment sample column names
    control_columns : list
        List of control sample column names
    alpha : float, default=0.05
        Significance threshold for adjusted p-values
    lfc_threshold : float, default=1.0
        Log2 fold change threshold for significance
        
    Returns:
    --------
    tuple
        Results DataFrame and summary statistics
    """
    # Input validation
    if not all(col in expression_matrix.columns for col in treatment_columns + control_columns):
        raise ValueError("Specified columns not found in expression matrix")

    # Initialize results collection
    results = []

    # Perform gene-wise differential expression analysis
    for gene in expression_matrix.index:
        try:
            # Extract and validate group-wise expression values
            treated = pd.to_numeric(expression_matrix.loc[gene, treatment_columns], errors='coerce')
            control = pd.to_numeric(expression_matrix.loc[gene, control_columns], errors='coerce')

            # Remove missing values
            treated = treated.dropna()
            control = control.dropna()

            # Validate sufficient data points
            if treated.empty or control.empty:
                continue

            # Calculate expression statistics
            mean_control = np.mean(control)
            mean_treated = np.mean(treated)

            # Compute fold change with pseudo-count
            log2fc = np.log2((mean_treated + 1) / (mean_control + 1))

            # Perform Welch's t-test (equal_var=False)
            t_stat, p_val = ttest_ind(treated, control, equal_var=False)

            # Compile gene-wise results
            results.append({
                "gene": gene,
                "Gene_Name": expression_matrix.loc[gene, "Gene_Name"] if "Gene_Name" in expression_matrix.columns else gene,
                "log2fc": log2fc,
                "mean_treated": mean_treated,
                "mean_control": mean_control,
                "t_stat": t_stat,
                "p_val": p_val,
                "var_treated": np.var(treated),
                "var_control": np.var(control)
            })

        except Exception as e:
            print(f"Warning: Error processing gene {gene}: {str(e)}")
            continue

    # Convert to DataFrame and perform quality control
    results_df = pd.DataFrame(results)
    results_df['p_val'] = pd.to_numeric(results_df['p_val'], errors='coerce')
    results_df = results_df.dropna(subset=['p_val'])

    # Apply multiple testing correction
    results_df['p_adj'] = multipletests(results_df['p_val'], method='fdr_bh')[1]

    # Calculate absolute fold change
    results_df['abs_log2fc'] = results_df['log2fc'].abs()

    # Define significance criteria
    results_df['significant'] = (results_df['p_adj'] < alpha) & \
                               (results_df['abs_log2fc'] > lfc_threshold)

    # Generate summary statistics
    summary_stats = {
        'total_genes': int(len(results_df)),
        'significant_genes': int(results_df['significant'].sum()),
        'up_regulated': int(sum((results_df['significant']) & (results_df['log2fc'] > 0))),
        'down_regulated': int(sum((results_df['significant']) & (results_df['log2fc'] < 0))),
        'mean_variance_ratio': float(np.mean(results_df['var_treated'] / results_df['var_control']))
    }

    # Sort by statistical significance
    results_df = results_df.sort_values('p_adj')

    print("\nDifferential Expression Analysis Summary:")
    print(f"Total genes analyzed: {summary_stats['total_genes']}")
    print(f"Significant genes: {summary_stats['significant_genes']}")
    print(f"Up-regulated: {summary_stats['up_regulated']}")
    print(f"Down-regulated: {summary_stats['down_regulated']}")
    print(f"Mean variance ratio (treated/control): {summary_stats['mean_variance_ratio']:.2f}")

    return results_df, summary_stats

# Get input, output and parameters from Snakemake
input_file = snakemake.input.filtered_counts
output_file = snakemake.output.deg_results
output_stats = snakemake.output.stats
treatment_columns = snakemake.params.treatment_columns
control_columns = snakemake.params.control_columns
alpha = snakemake.params.alpha
lfc_threshold = snakemake.params.lfc_threshold

# Create output directories if they don't exist
os.makedirs(os.path.dirname(output_file), exist_ok=True)

# Load the filtered and normalized count matrix
filtered_normalized_count_matrix = pd.read_csv(input_file, index_col=0)

# Perform differential expression analysis
results_df, summary_stats = analyze_differential_expression(
    expression_matrix=filtered_normalized_count_matrix,
    treatment_columns=treatment_columns,
    control_columns=control_columns,
    alpha=alpha,
    lfc_threshold=lfc_threshold
)

# Save the differential expression results
results_df.to_csv(output_file)

# Save stats as JSON
with open(output_stats, 'w') as f:
    json.dump(summary_stats, f, indent=4)

print(f"Differential expression analysis complete: {summary_stats['significant_genes']} significant genes identified")
