import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
from scipy import sparse, stats
from sklearn import mixture
import time
from tqdm import tqdm


def call_presence_with_gmm_ab(umi_counts: np.ndarray, n_components: int = 2) -> np.ndarray:
    '''
    call_presence_with_gmm_ab function from cellranger github repo: https://github.com/10XGenomics/cellranger/blob/cfa9ac1a9a0fcfbd123dc574934a6a72889d2f70/lib/python/cellranger/feature/feature_assigner.py#L213C1-L235C1
    Given the UMI counts for a specific antibody, separate signal from background.
    '''
    if np.max(umi_counts) == 0 or max(umi_counts.shape) < 2:
        # there are no UMIs, or only one UMI, each barcode has 0 count
        return np.repeat(False, len(umi_counts))

    # Turn the numpy array into 2D log scale array
    umi_counts = np.reshape(umi_counts, (len(umi_counts), 1))
    log_ab = np.log10(umi_counts + 1.0)

    # Initialize and fit the gaussian Mixture Model
    gmm = mixture.GaussianMixture(n_components, n_init=10, covariance_type="tied", random_state=0)
    gmm.fit(log_ab)

    # Calculate the threshold
    umi_posterior = gmm.predict_proba(log_ab)
    high_umi_component = np.argmax(gmm.means_)

    in_high_umi_component = umi_posterior[:, high_umi_component] > 0.5

    return in_high_umi_component


def ga_cellranger(input_file, batch_list, output_dir, start_gRNA = 0, step = None):
    '''
    Guide assignment in which a Gaussian mixture model is fitted to the log-transformed UMI counts
    
    Args:
        input_file: (str) path to the stored anndata object with the gRNA counts
        batch_list: (list) list of batches for which to fit the mixture model
        output_dir: (str) directory in which to store the resulting assignment
        start_gRNA: (int, optional) index of the start gRNA when parallelizing assignment for gRNA sets
        step: (int, optional) number of gRNAs for which the assignment is done (if set to None, assignment for all gRNAs in the data)
    
    Returns:
        None
    ''' 
    print('Guide assignment according to cellranger tool')
    # Load gRNA counts data
    print('Load gRNA counts')
    adata_crispr = sc.read_h5ad(input_file)
    
    for batch in batch_list:
        print('Fit Gaussian Mixture model for batch ' + str(batch))
        # If output_dir doesn't exist, the output folders are created
        if not os.path.exists(output_dir + 'batch' + str(batch) + '/'):
            os.makedirs(output_dir + 'batch' + str(batch) + '/')
            print("The output directory " + output_dir + 'batch' + str(batch) + '/' +  " was created")

        # Subset to selected batch
        adata_crispr_batch = adata_crispr[adata_crispr.obs['batch'] == batch]
        gRNA_list = adata_crispr_batch.var_names.tolist()

        if step != None:
            end_gRNA = start_gRNA + step - 1
            if end_gRNA >= len(gRNA_list):
                end_gRNA = len(gRNA_list) - 1
            gRNA_list = gRNA_list[start_gRNA:(end_gRNA + 1)]

        # Fit Gaussian Mixture Model (GMM) for each gRNA
        perturbations = pd.DataFrame({'cell': [], 'gRNA': []})
        
        for gRNA in tqdm(gRNA_list):
            time.sleep(0.01)
            # Select data for one gRNA
            selected_guide = adata_crispr_batch[:,[gRNA]].X
            data = selected_guide.toarray() 
            
            # Fit gaussian mixture model
            perturbed = call_presence_with_gmm_ab(data)
            perturbed_cells = adata_crispr_batch.obs_names[perturbed].tolist()
            df = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA})
            perturbations = pd.concat([perturbations, df], ignore_index = True)

        # Save data frame with the perturbations assigned to each cell
        if step == None:
            perturbations.to_csv(output_dir + 'batch' + str(batch) + '/perturbations.csv', index = False)
        else:
            perturbations.to_csv(output_dir + 'batch' + str(batch) + '/perturbations_'+str(start_gRNA)+'-'+str(end_gRNA)+'.csv',
                                 index = False)
    
    # Combine the outputs of all batches
    if step == None:
        pert = pd.DataFrame()
        for batch in batch_list:
            b = pd.read_csv(output_dir + 'batch' + str(batch) + '/perturbations.csv')
            pert = pd.concat([pert, b])

        pert.to_csv(output_dir + 'perturbations.csv')
    print('Done: outputs are saved in ' + output_dir)
       
    
if __name__ == "__main__":
    main()