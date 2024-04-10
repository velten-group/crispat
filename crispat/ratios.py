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


def ga_ratio(input_file, thresholds, output_dir):
    '''
    Guide assignment in which the most abundant gRNA per cell is assigned if it comprises 
    more than X% of the total counts in a cell
    
    Args:
        input_file: (str) path to the stored anndata object with the gRNA counts
        thresholds: (list) list of ratio thresholds to use (generates one output file per ratio)
        output_dir: (str) directory in which to store the resulting assignment
    
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
    adata_crispr = adata_crispr[adata_crispr.obs['total_counts'] != 0].copy()
    
    # Normalization to get the percentages of the gRNA counts in each cell
    sc.pp.normalize_total(adata_crispr, target_sum = 1)
    
    # Convert adata object to long df
    percent_df = pd.DataFrame(adata_crispr.X.todense(), 
                              index = adata_crispr.obs_names, 
                              columns = adata_crispr.var_names)
    percent_df['cell'] = percent_df.index
    percent_df = pd.melt(percent_df, id_vars=['cell'], var_name='gRNA', value_name='percent_counts')
    
    # Get maximum per cell
    max_df = percent_df.groupby('cell').agg({'percent_counts': max})
    max_df.to_csv(output_dir + 'max_df.csv')
    max_df = max_df.reset_index() 
    max_df = max_df.merge(percent_df, on = ['cell', 'percent_counts'])

    # Save the results
    for t in thresholds: 
        perturbations = max_df[max_df['percent_counts'] > t]
        perturbations.to_csv(output_dir + 'perturbations_t' + str(t) + '.csv', index = False)
        
    print('Done: outputs are saved in ' + output_dir)
    
    
if __name__ == "__main__":
    main()