# Script to download and load RNA-seq data

import pandas as pd
import subprocess
import gzip
import os

def download_and_load_data(url, output_filename, sep="\t", column_filter=None):
    """
    Download and load RNA-seq data from URL, with optional column filtering.
    
    Parameters:
    -----------
    url : str
        URL of the data file
    output_filename : str
        Path to save the downloaded file
    sep : str, default="\t"
        Separator in the data file
    column_filter : str or None
        If provided, only columns containing this string will be kept
        
    Returns:
    --------
    pandas.DataFrame
        Loaded data
    """
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    
    # Download the file using wget
    print(f"Downloading data from {url}...")
    subprocess.run(["wget", "-O", output_filename + ".gz", url], check=True)

    # Unzip file using gunzip
    print(f"Unzipping {output_filename}.gz...")
    with gzip.open(output_filename + ".gz", "rb") as gz_file:
        with open(output_filename, "wb") as out_file:
            out_file.write(gz_file.read())

    # Load the data into a Pandas dataframe
    print(f"Loading {output_filename} into a pandas DataFrame...")
    df = pd.read_csv(output_filename, sep=sep, index_col=0)

    # Optionally filter columns based on keyword
    if column_filter:
        print(f"Filtering columns with keyword '{column_filter}'...")
        filtered_columns = [col for col in df.columns if column_filter in col]
        df = df[filtered_columns]

    # Return pandas data frame
    return df

# Get parameters from Snakemake
url = snakemake.params.url
column_filter = snakemake.params.column_filter
output_file = snakemake.output.count_matrix

# Download and process the data
count_matrix = download_and_load_data(
    url=url,
    output_filename=output_file,
    column_filter=column_filter
)

# Save the processed count matrix
count_matrix.to_csv(output_file)
