import os
import sys, getopt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix


def ga_max(input_file, output_dir, UMI_threshold = 0):
    '''
    Guide assignment in which the most abundant gRNA per cell is assigned 
    Args:
        input_file (str): path to the stored anndata object with the gRNA counts
        output_dir (str): directory in which to store the resulting assignment
        UMI_threshold (int, optional): Additional UMI threshold for assigned cells which is applied after creating the initial assignment to remove cells with fewer UMI counts than this threshold (default: no additional UMI threshold)
    Returns:
        None
    '''
    print('Guide assignment with maximum assignment')
    
    # If output_dir doesn't exist, the output folders are created
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("The output directory " + output_dir +  " was created") 

    # Load gRNA counts data
    print('Load gRNA counts')
    adata_crispr = sc.read_h5ad(input_file)
    gRNA_list = adata_crispr.var_names.tolist()

    # Convert adata object to long df
    count_df = pd.DataFrame(adata_crispr.X.todense(), 
                              index = adata_crispr.obs_names, 
                              columns = adata_crispr.var_names)
    count_df['cell'] = count_df.index
    count_df = pd.melt(count_df, id_vars=['cell'], var_name='gRNA', value_name='UMI_counts')
    
    # Get maximum counts per cell
    count_df['max'] = count_df.groupby('cell')['UMI_counts'].transform('max')
    max_df = count_df[count_df['UMI_counts'] == count_df['max']]
    
    # Remove cells with maximum of 0
    max_df = max_df[max_df['UMI_counts'] != 0]                  
    # Remove cells where multiple gRNAs share the maximum
    occurences = max_df.groupby(['cell']).size().reset_index(name='n_max_gRNAs')
    max_df = max_df.merge(occurences, on = ['cell'])
    max_df = max_df[max_df['n_max_gRNAs'] == 1]
    
    # Optional filtering to assigned cells that have at least 'UMI_threshold' counts
    max_df = max_df[max_df['UMI_counts'] >= UMI_threshold]
     
    # Save data frames with the results
    max_df[['cell', 'gRNA', 'UMI_counts']].to_csv(output_dir + 'assignments.csv', index = False)
    print('Done: outputs are saved in ' + output_dir)
    
    
if __name__ == "__main__":
    pass
