import os
import sys, getopt
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix


def ga_umi(input_file, thresholds, output_dir):   
    '''
    Guide assignment with fixed UMI thresholds 
    
    Args:
        input_file (str): path to the stored anndata object with the gRNA counts
        thresholds (list): list of integers to use as thresholds (create assignment output file for each t in the list)
        output_dir (str): directory in which to store the resulting assignment
        
    Returns:
        None
    '''
    print('Guide assignment with UMI threshold')
    # If output_dir doesn't exist, the output folders are created
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print("The output directory " + output_dir +  " was created") 

    # Load gRNA counts data
    print('Load gRNA counts')
    adata_crispr = sc.read_h5ad(input_file)
    gRNA_list = adata_crispr.var_names.tolist()

    # Get perturbed cells for each gRNA based on a fixed UMI threshold
    for threshold in thresholds:
        perturbations = pd.DataFrame({'cell': [], 'gRNA': []})
        print('Get perturbed cells for each gRNA with UMI threshold = ' + str(threshold))
        for gRNA in gRNA_list:
            # Get cells with UMI counts higher than the threshold for specified gRNA
            selected_guide = adata_crispr[:,[gRNA]].X
            perturbed_cells = adata_crispr.obs_names[selected_guide.toarray().reshape(-1) >= threshold].tolist()
            UMI_counts = adata_crispr[selected_guide.toarray().reshape(-1) >= threshold, [gRNA]].X.toarray().reshape(-1)
           
            if len(perturbed_cells) != 0:
                df = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA, 'UMI_counts': UMI_counts})
                perturbations = pd.concat([perturbations, df], ignore_index = True)

        # Save data frames with the results
        perturbations.to_csv(output_dir + 'assignments_t' + str(threshold) + '.csv', index = False)
        
    print('Done: outputs are saved in ' + output_dir)
    
    
if __name__ == "__main__":
    pass
