import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
from itertools import chain
from scipy import sparse, stats, special
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
import time

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.infer.autoguide import AutoDelta
from pyro.optim import Adam
from pyro.infer import SVI, TraceEnum_ELBO, config_enumerate
from torch.distributions.mixture_same_family import MixtureSameFamily
    
    
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


@config_enumerate
def model(data):
    '''
    Mixture Model of 2 beta distributions on the ratios
    '''
    # Global variables
    weights = pyro.sample("weights", dist.Dirichlet(torch.tensor([0.6, 0.4])))
    with pyro.plate("components", 2):
        alphas = pyro.sample("alphas", dist.LogNormal(1.0, 1.0))
        betas = pyro.sample("betas", dist.LogNormal(1.0, 1.0))

    with pyro.plate("data", len(data)):
        # Local variables
        assignment = pyro.sample("assignment", dist.Categorical(weights))
        component_distribution = dist.Beta(alphas[assignment], betas[assignment])
        pyro.sample("obs", component_distribution, obs=data)
        
        
def init_loc_fn(site):
    '''
    Define initial parameter values
    '''
    if site["name"] == "weights":
        return torch.tensor([0.6, 0.4])
    if site["name"] == "alphas":
        return torch.tensor([1.0, 10.0])
    if site["name"] == "betas":
        return torch.tensor([10.0, 1.0])
    
    raise ValueError(site["name"])

    
def initialize(seed, optim, elbo, data):
    '''
    Initialization for SVI 
    Args:
        seed (str): seed that is used in pyro
        optim: pyro optimizer
        elbo: pyro loss function
        data (tensor): observed transformed gRNA counts
    Returns:
        Initial loss
    '''
    global global_guide, svi
    pyro.set_rng_seed(seed)
    pyro.clear_param_store()
    global_guide = AutoDelta(
        poutine.block(model, expose=["weights", "alphas", "betas"]),
        init_loc_fn=init_loc_fn,
    )
    svi = SVI(model, global_guide, optim, loss=elbo)
    return svi.loss(model, global_guide, data)

    
def plot_fitted_model(data, weights, alphas, betas, threshold, batch, output_dir):
    '''
    Saves a plot of the data histogram and the fitted mixture model
    Args:
        data (tensor): observed transformed gRNA counts
        weights (np array): estimated proportion for the components
        alphas (np array): MAP estimate for the 2 alphas of the beta distributions
        betas (np array): MAP estimate for the 2 betas of the beta distributions
        threshold (int): threshold for assigning it to the Normal component
        batch (str): name of the batch used for the plot title
        output_dir (str): name of the output directory
    Returns:
        None
    '''
    X = np.arange(0, 1, 0.01)
    Y1 = weights[0] * stats.beta.pdf(X, alphas[0], betas[0]) * 0.01
    Y2 = weights[1] * stats.beta.pdf(X, alphas[1], betas[1]) * 0.01

    fig, ax = plt.subplots(figsize=(8, 3), dpi=300)
    sns.histplot(data, binwidth=0.01, color='grey', stat = "proportion")
    ax.plot(X, Y1, "r-", label = "Beta 1")
    ax.plot(X, Y2, "b-", label = "Beta 2")
    ax.plot(X, Y1 + Y2, "k--", label = "Mixture model")
    ax.axvline(threshold, c = "orange", label = "Threshold")
    plt.legend()
    plt.title("2-Beta mixture model for batch " + batch)
    plt.ylabel("Probability Density")
    plt.xlabel("Proportion highest gRNA")
    plt.savefig(output_dir + "fitted_model_plots/fitted_model_batch-" + batch + ".png", bbox_inches="tight")
    plt.close()
    

def fit_betaMM(data, batch, output_dir, seed, n_iter):
    '''
    Fits mixture model of 2 Beta distributions on the ratios of the most abundant gRNA 
    Args:
        data (tensor): tensor containing the proportions of the highest gRNA per cell 
        batch (str): number of the batch
        output_dir (str): directory in which the resulting plots will be saved
        seed (int): seed used for pyro
        n_iter (int): number of steps for training the model
    Returns:
        Fitted model parameters and threshold
    '''
    # Set optimizer and elbo parameters
    optim = pyro.optim.Adam({"lr": 0.01, "betas": [0.8, 0.99]})
    elbo = TraceEnum_ELBO(max_plate_nesting=1)
    
    # Initialization of SVI
    initialize(seed, optim, elbo, data)

    # Train the model n_iter steps with early stopping when loss doesn't change at least 0.001 for 50 steps
    losses = []
    min_loss = 1.e8
    last_step = 0

    for step in range(n_iter):
        loss = svi.step(data)
        if loss < min_loss - 0.001:
            min_loss = loss
            last_step = step
        losses.append(loss)
        if (step - last_step) > 50:
            break
        
    # MAP estimates of the model
    map_estimates = global_guide(data)
    weights = map_estimates["weights"].data.numpy()
    alphas = map_estimates["alphas"].data.numpy()
    betas = map_estimates["betas"].data.numpy()
    estimates = pd.DataFrame({'batch': [batch], 
                              'weight_Beta1': [weights[0]], 
                              'weight_Beta2': [weights[1]], 
                              'alpha_Beta1': [alphas[0]], 
                              'alpha_Beta2': [alphas[1]], 
                              'beta_Beta1': [betas[0]], 
                              'beta_Beta2': [betas[1]]})

    # find the threshold
    threshold = get_threshold(weights, alphas, betas)
    
    # plot fitted model
    plot_fitted_model(data, weights, alphas, betas, threshold, batch, output_dir)
    
    return(estimates, threshold)


def get_threshold(weights, alphas, betas):
    '''
    Finds the threshold from which on the probability of the 2nd mixture component is highest 
    Args:
        estimates (pd DataFrame): parameters estimated by the model
    Returns:
        Threshold
    '''
    X = np.arange(0.01, 1, 0.01)
    Y1 = weights[0] * stats.beta.pdf(X, alphas[0], betas[0]) 
    Y2 = weights[1] * stats.beta.pdf(X, alphas[1], betas[1]) 
    probs = pd.DataFrame({'x': X, 'Y1': Y1, 'Y2': Y2})
    probs['max'] = probs[['Y1', 'Y2']].max(axis=1)
    threshold = probs.loc[(probs['max'] == probs['Y2']), 'x'].min()
    return threshold


def ga_2beta(input_file, output_dir, n_iter = 500, batch_list = None, add_UMI_counts = True, UMI_threshold = 0):
    '''
    Guide assignment in which a mixture model of 2-Beta distributions is fitted to the ratios of 
    the most abundant gRNAs per cell to determine a batch-specific threshold on the ratios 

    Args:
        input_file (str): Path to the stored anndata object with the gRNA counts
        output_dir (str): Directory in which to store the resulting assignment
        n_iter (int, optional): Number of steps for training the model (default is 500)
        batch_list (list, opitional): List of batches for which to fit the mixture model. If none (default), all available batches are used
        add_UMI_counts (bool, optional): if true, UMI counts are added to the output. To improve run time, set it to False
        UMI_threshold (int, optional): Additional UMI threshold for assigned cells which is applied after creating the initial assignment to remove cells with fewer UMI counts than this threshold (default: no additional UMI threshold)

    Returns:
        None
    '''
    print('Guide assignment fitting a mixture model of two beta distributions')
    
    # If output_dir doesn't exist, the output folders are created
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        os.makedirs(output_dir + "fitted_model_plots/")
        print("The output directory " + output_dir +  " was created")
        
    # Load gRNA counts
    print('Load gRNA counts')
    adata_crispr = sc.read_h5ad(input_file)
    # subset to specified batches
    if batch_list != None:
        adata_crispr = adata_crispr[adata_crispr.obs['batch'].isin(batch_list)].copy()
    else:
        batch_list = adata_crispr.obs['batch'].unique()
        
    # Remove cells with 0 gRNA counts
    adata_crispr.obs['total_counts'] = np.array(adata_crispr.X.sum(axis = 1)).flatten()
    adata_crispr_norm = adata_crispr[adata_crispr.obs['total_counts'] != 0].copy()
    
    # Calculate ratios per cell
    data = get_ratios(adata_crispr_norm)
    
    perturbations = pd.DataFrame()
    thresholds = pd.DataFrame()
    estimates = pd.DataFrame()
        
    # Fit model for each batch
    print('Fit 2-beta mixture model for each batch:')
    for batch in batch_list:
        print('Batch ', str(batch))
            
        # Subset data to selected batch
        batch_data = data[data['batch'] == batch]
        
        # Restrict values to the (0, 1) interval (model fails for values that are exactly 0 or 1)
        batch_data.loc[batch_data['percent_counts'] > 0.9999, 'percent_counts'] = 0.9999
        batch_data.loc[batch_data['percent_counts'] < 0.0001, 'percent_counts'] = 0.0001
        
        # Get maximum proportion per cell
        max_df = batch_data.groupby('cell').agg({'percent_counts': "max"})
        
        # Fit 2-Beta Mixture Model
        map_estimates, threshold = fit_betaMM(torch.tensor(list(max_df['percent_counts'])), 
                                              str(batch), output_dir, 2024, n_iter)
        
        # Get perturbed cells per batch
        max_df = max_df.reset_index() 
        max_df = max_df.merge(batch_data, on = ['cell', 'percent_counts'])
        perturbed_cells = max_df[max_df['percent_counts'] > threshold]
        
        # Add batch results to the dataframes
        perturbations = pd.concat([perturbations, perturbed_cells], ignore_index = True)
        thresholds = pd.concat([thresholds, pd.DataFrame({'batch': [batch], 'threshold': [threshold]})])
        estimates = pd.concat([estimates, map_estimates])

    # Add the UMI counts 
    if add_UMI_counts:
        umi_counts = []
        for _, row in perturbations.iterrows():
            # Get UMI counts for selected gRNA in selected cell
            umi_count = adata_crispr[row['cell'], row['gRNA']].X.toarray().item()
            umi_counts.append(umi_count)

        perturbations['UMI_counts'] = umi_counts
        perturbations = perturbations[perturbations['UMI_counts'] >= UMI_threshold]

    # Save data frames with the results
    perturbations.to_csv(output_dir + 'assignments.csv', index = False)
    thresholds.to_csv(output_dir + 'batch_thresholds.csv', index = False)
    estimates.to_csv(output_dir + 'MAP_estimates.csv', index = False)
    
    print('Done: outputs are saved in ' + output_dir)
    
    
if __name__ == "__main__":
    pass
