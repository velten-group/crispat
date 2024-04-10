import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

def create_anndata_from_csv(csv_file, save_dir):
    '''
    Creates an AnnData object from a csv file (rows: gRNAs, columns: cells; name of first column: gRNA)
    
    Args:
        csv_file (str): path of the csv file
        save_dir (str): directory in which to save the AnnData object that is used as input for guide assignment
    
    Returns:
        None
    '''
    # Load csv files with the gRNA counts
    print('Load data')
    data = pd.read_csv(csv_file).set_index('gRNA')

    # Create anndata object
    print('Create anndata object')
    data = ad.AnnData(data.transpose())
    data.obs['batch'] = 1
    data.X = csr_matrix(data.X)
    
    # Save as h5ad objects
    data.write(save_dir + 'gRNA_counts.h5ad')
    print('Done: AnnData object is saved in ' + save_dir)
    
    
def create_anndata_from_cellranger(batch_list, input_dir, save_dir):
    '''
    Creates an AnnData object from cellranger output (one folder per batch named batch1, batch2,...)
    
    Args:
        batch_list (list): list of batch numbers 
        input_dir (list): directory containing the cellranger output of each batch
        save_dir (str): directory in which to save the AnnData object that is used as input for guide assignment 
    
    Returns:
        None
    '''
    batches = []
    print('Load data')
    for batch_number in batch_list:
        # Load data of one batch 
        batch = sc.read_10x_mtx(input_dir + 'batch' + str(batch_number)+'/', 
                                var_names='gene_symbols',
                                gex_only = False, 
                                prefix = '') 
        batch.obs_names = [name.split('-')[0] for name in batch.obs_names]
        batch.obs['batch'] = batch_number
        batches.append(batch)
        
    # Concatenate batches into one AnnData object
    print('Create concatenated anndata object')
    adata = ad.concat(batches, merge = "same", label = "batch", keys = batch_list, index_unique="-")
    
    # Subset to CRISPR gRNA features
    adata_crispr = adata[:, adata.var["feature_types"] == "CRISPR Guide Capture"]
    
    # Save as h5ad object
    adata_crispr.write(save_dir + 'gRNA_counts.h5ad')
    print('Done: AnnData object is saved in ' + save_dir)
    