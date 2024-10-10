import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from itertools import chain


def ga_ratio(input_file, thresholds, output_dir, add_UMI_counts = True, UMI_threshold = 0):
    '''
    Guide assignment in which the most abundant gRNA per cell is assigned if it comprises 
    more than X% of the total counts in a cell
    
    Args:
        input_file: (str) path to the stored anndata object with the gRNA counts
        thresholds: (list) list of ratio thresholds to use (generates one output file per ratio)
        output_dir: (str) directory in which to store the resulting assignment
        add_UMI_counts: (bool) if true, UMI counts are added to the output. To improve run time, set it to False
        UMI_threshold (int, optional): Additional UMI threshold for assigned cells which is applied after creating the initial assignment to remove cells with fewer UMI counts than this threshold (default: no additional UMI threshold)
    
    Returns:
        None
    '''
    print('Guide assignment with ratio assignment')
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("The output directory " + output_dir +  " was created")
    
    # Load gRNA counts data
    print('Load gRNA counts')
    adata_crispr = sc.read_h5ad(input_file)
    
    # Remove cells with 0 gRNA counts
    adata_crispr.obs['total_counts'] = np.array(adata_crispr.X.sum(axis = 1)).flatten()
    adata_crispr_norm = adata_crispr[adata_crispr.obs['total_counts'] != 0].copy()
    
    # Normalization to get the percentages of the gRNA counts in each cell
    sc.pp.normalize_total(adata_crispr_norm, target_sum = 1)
    
    # Convert adata object to long df
    percent_df = pd.DataFrame(adata_crispr_norm.X.todense(), 
                              index = adata_crispr_norm.obs_names, 
                              columns = adata_crispr_norm.var_names)
    percent_df['cell'] = percent_df.index
    percent_df = pd.melt(percent_df, id_vars=['cell'], var_name='gRNA', value_name='percent_counts')

    # Get maximum per cell
    max_df = percent_df.groupby('cell').agg({'percent_counts': 'max'})
    #max_df.to_csv(output_dir + 'max_df.csv')
    max_df = max_df.reset_index() 
    max_df = max_df.merge(percent_df, on = ['cell', 'percent_counts'])

    # Add the UMI counts 
    if add_UMI_counts:
        umi_counts = []
        for _, row in max_df.iterrows():
            # Get UMI counts for selected gRNA in selected cell
            umi_count = adata_crispr[row['cell'], row['gRNA']].X.toarray().item()
            umi_counts.append(umi_count)

        max_df['UMI_counts'] = umi_counts
        # Optional filtering to assigned cells that have at least 'UMI_threshold' counts
        max_df = max_df[max_df['UMI_counts'] >= UMI_threshold]
     


    # Save the results
    for t in thresholds: 
        perturbations = max_df[max_df['percent_counts'] > t]
        perturbations.to_csv(output_dir + 'assignments_t' + str(t) + '.csv', index = False)
        
    print('Done: outputs are saved in ' + output_dir)
    
    
if __name__ == "__main__":
    pass