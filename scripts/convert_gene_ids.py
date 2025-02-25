# Script to convert Ensembl IDs to gene symbols

import pandas as pd
import mygene

def convert_ensembl_to_gene_symbols(count_matrix, species='human'):
    """
    Convert Ensembl gene IDs to gene symbols.
    
    Parameters:
    -----------
    count_matrix : pandas.DataFrame
        Count matrix with Ensembl IDs as index
    species : str, default='human'
        Species for gene ID conversion
        
    Returns:
    --------
    pandas.DataFrame
        Count matrix with gene symbols added as a column
    """
    try:
        # Create a copy to avoid modifying the original
        count_matrix = count_matrix.copy()

        # Remove version numbers from Ensembl IDs
        cleaned_index = count_matrix.index.str.split('.').str[0]
        count_matrix.index = cleaned_index

        # Initialize MyGeneInfo object and query gene symbols
        mg = mygene.MyGeneInfo()
        ensembl_ids = count_matrix.index.unique().tolist()

        # Query gene information with error handling
        gene_info = mg.querymany(ensembl_ids, scopes='ensembl.gene', fields='symbol', species=species, verbose=False)

        # Convert to DataFrame and clean results
        gene_df = pd.DataFrame(gene_info)
        gene_df = gene_df.dropna(subset=['symbol'])
        gene_df = gene_df.drop_duplicates(subset='query')

        # Map gene symbols to count matrix
        symbol_map = gene_df.set_index('query')['symbol']
        count_matrix['Gene_Name'] = count_matrix.index.map(symbol_map)

        # Reorganize columns with Gene_Name first
        cols = ['Gene_Name'] + [col for col in count_matrix.columns if col != 'Gene_Name']
        count_matrix = count_matrix[cols]

        # Log conversion statistics
        total_genes = len(ensembl_ids)
        mapped_genes = len(gene_df)
        print(f"Successfully mapped {mapped_genes} out of {total_genes} genes ({mapped_genes/total_genes*100:.1f}%)")

        return count_matrix

    except Exception as e:
        raise Exception(f"Error during gene ID conversion: {str(e)}")

# Get input and output files from Snakemake
input_file = snakemake.input.count_matrix
output_file = snakemake.output.gene_counts
species = snakemake.params.species

# Load the count matrix
count_matrix = pd.read_csv(input_file, index_col=0)

# Convert gene IDs
count_matrix_gene_names = convert_ensembl_to_gene_symbols(count_matrix, species=species)

# Save the processed count matrix
count_matrix_gene_names.to_csv(output_file)
