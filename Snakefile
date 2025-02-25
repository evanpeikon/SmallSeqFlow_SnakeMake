# Snakefile for SmallSeqFlow RNA-seq pipeline

# Configuration
configfile: "config.yaml"

# Import required libraries
import pandas as pd
import os

# Get sample information from config file
samples = config["samples"]
conditions = config["conditions"]
treatment_samples = config["treatment_samples"]
control_samples = config["control_samples"]

# Define final target rule (output files we want at the end)
rule all:
    input:
        "results/quality_control/main_qc_figure.png",
        "results/quality_control/dendrogram_figure.png",
        "results/filtered_normalized_counts.csv",
        "results/differential_expression/DEG_results.csv",
        "results/plots/differential_expression_plot.png"

# Rule to download and load RNA-seq data
rule download_data:
    output:
        count_matrix = "data/raw/count_matrix.csv"
    params:
        url = config["data_url"],
        column_filter = config.get("column_filter", None)
    script:
        "scripts/download_data.py"

# Rule to convert Ensembl IDs to gene symbols
rule convert_gene_ids:
    input:
        count_matrix = "data/raw/count_matrix.csv"
    output:
        gene_counts = "data/processed/count_matrix_gene_names.csv"
    params:
        species = config.get("species", "human")
    script:
        "scripts/convert_gene_ids.py"

# Rule to perform QC on the RNA-seq data
rule quality_control:
    input:
        gene_counts = "data/processed/count_matrix_gene_names.csv"
    output:
        main_fig = "results/quality_control/main_qc_figure.png",
        dendrogram_fig = "results/quality_control/dendrogram_figure.png",
        stats = "results/quality_control/qc_stats.json"
    script:
        "scripts/quality_control.py"

# Rule to filter low-expressed genes and normalize the counts
rule filter_normalize:
    input:
        gene_counts = "data/processed/count_matrix_gene_names.csv"
    output:
        filtered_counts = "results/filtered_normalized_counts.csv",
        stats = "results/filtering_stats.json",
        cpm_plot = "results/plots/genes_by_cpm_threshold.png"
    params:
        min_cpm = config.get("min_cpm", 1.0),
        min_samples = config.get("min_samples", 2)
    script:
        "scripts/filter_normalize.py"

# Rule to perform differential expression analysis
rule differential_expression:
    input:
        filtered_counts = "results/filtered_normalized_counts.csv"
    output:
        deg_results = "results/differential_expression/DEG_results.csv",
        stats = "results/differential_expression/de_stats.json"
    params:
        treatment_columns = treatment_samples,
        control_columns = control_samples,
        alpha = config.get("alpha", 0.05),
        lfc_threshold = config.get("lfc_threshold", 1.0)
    script:
        "scripts/differential_expression.py"

# Rule to visualize differential expression results
rule visualize_de:
    input:
        deg_results = "results/differential_expression/DEG_results.csv",
        filtered_counts = "results/filtered_normalized_counts.csv"
    output:
        de_plot = "results/plots/differential_expression_plot.png"
    params:
        treatment_columns = treatment_samples,
        control_columns = control_samples,
        p_adj_threshold = config.get("alpha", 0.05),
        abs_log2fc_threshold = config.get("lfc_threshold", 1.0)
    script:
        "scripts/visualize_de.py"
