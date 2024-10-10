import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
from itertools import chain
from scipy import sparse, stats, special
import matplotlib.pyplot as plt
from tqdm import tqdm
import time

    
def get_ratios(adata_crispr):
    '''
    Calculates the proportion of guide counts per cell
    Args:
        adata_crispr (AnnData): CRISPR gRNA counts per cell
    Returns:
        A long dataframe with the proportion of guide counts per cell
    '''
    # Normalization to get the percentages of the gRNA counts in each cell
    sc.pp.normalize_total(adata_crispr, target_sum = 1)
    
    # Convert adata object to long df
    percent_df = pd.DataFrame(adata_crispr.X.todense(), 
                              index = adata_crispr.obs_names, 
                              columns = adata_crispr.var_names)
    percent_df['cell'] = percent_df.index
    percent_df = pd.melt(percent_df, id_vars=['cell'], var_name='gRNA', value_name='percent_counts')

    # Add batch column
    percent_df['batch'] = [adata_crispr.obs.loc[cell, 'batch'] for cell in percent_df['cell']]

    return percent_df


def ga_quantiles(input_file, thresholds, output_dir, UMI_threshold = 0):
    '''
    Guide assignment in which the X% non-zero cells with highest ratios are assigned per gRNA
    
    Args:
        input_file (str): path to the stored anndata object with the gRNA counts
        thresholds (list): list of quantile thresholds for which to return the assignment
        output_dir (str): directory in which to store the resulting assignment
        UMI_threshold (int, optional): Additional UMI threshold for assigned cells which is applied after creating the initial assignment to remove cells with fewer UMI counts than this threshold (default: no additional UMI threshold)
    
    Returns:
        None
    '''   
    # If output_dir doesn't exist, the output folders are created
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("The output directory " + output_dir +  " was created")
    
    # Load gRNA counts data
    print('Load gRNA counts data')
    adata_crispr = sc.read_h5ad(input_file)
    
    # Remove cells with 0 gRNA counts
    adata_crispr.obs['total_counts'] = np.array(adata_crispr.X.sum(axis = 1)).flatten()
    adata_crispr_norm = adata_crispr[adata_crispr.obs['total_counts'] != 0].copy()
 
    # Calculate ratios per cell
    data = get_ratios(adata_crispr_norm)
        
    # Get the cells with highest ratios per gRNA
    print('Get cells with highest ratios per gRNA')
    data = data[data['percent_counts'] != 0].sort_values(by=['gRNA', 'percent_counts'], ascending=[True, False])

    # Add the UMI counts 
    umi_counts = []
    for _, row in data.iterrows():
        # Get UMI counts for selected gRNA in selected cell
        umi_count = adata_crispr[row['cell'], row['gRNA']].X.toarray().item()
        umi_counts.append(umi_count)
    data['UMI_counts'] = umi_counts

    for quantile in thresholds:
        perturbations = data.groupby('gRNA').apply(lambda x: x.head(int(quantile * len(x))))
        
        # Optional filtering to assigned cells that have at least 'UMI_threshold' counts
        perturbations = perturbations[perturbations['UMI_counts'] >= UMI_threshold]
     
        # Save the resulting dataframe
        perturbations.to_csv(output_dir + 'assignments_t' + str(quantile) + '.csv', index = False)
      
    print('Done: outputs are saved in ' + output_dir)
    
    
if __name__ == "__main__":
    pass
