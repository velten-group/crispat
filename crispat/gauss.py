import os
import sys, getopt
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
from scipy import sparse, stats
from scipy.sparse import csr_matrix
from sklearn import mixture
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


@config_enumerate
def model(data):
    '''
    Gaussian-Gaussian Mixture Model 
    '''
    # Global variables
    weights = pyro.sample("weights", dist.Dirichlet(torch.tensor([0.99, 0.01])))
    with pyro.plate("components", 2):
        locs = pyro.sample("locs", dist.Normal(1.0, 2.0))

    scale = pyro.sample("scales", dist.LogNormal(-3.0, 2.0))
    
    with pyro.plate("data", len(data)):
        # Local variables
        assignment = pyro.sample("assignment", dist.Categorical(weights)) 
        pyro.sample("obs", dist.Normal(locs[assignment], scale), obs=data)
        

def init_loc_fn(site):
    '''
    Define initial parameter values
    '''
    if site["name"] == "weights":
        return torch.tensor([0.99, 0.01])
    if site["name"] == "locs":
        return torch.tensor([0.0, 1.0])
    if site["name"] == "scales":
        return torch.tensor([0.01])
    raise ValueError(site["name"])

    
def initialize(seed, optim, elbo, data):
    '''
    Initialization for SVI 
    Args:
        seed: (str) seed that is used in pyro
        optim: pyro optimizer
        elbo: pyro loss function
        data: (tensor) observed transformed gRNA counts
    Returns:
        Initial loss
    '''
    global global_guide, svi
    pyro.set_rng_seed(seed)
    pyro.clear_param_store()
    global_guide = AutoDelta(
        poutine.block(model, expose=["weights", "locs", "scales"]),
        init_loc_fn = init_loc_fn,
    )
    svi = SVI(model, global_guide, optim, loss=elbo)
    return svi.loss(model, global_guide, data)


def plot_loss(losses, gRNA, output_dir):
    '''
    Saves a plot of the loss over the SVI steps
    Args:
        losses: (list) loss over the SVI steps
        gRNA: (str) name of the gRNA used for the plot title
        output_dir: (str) name of the output directory
    Returns:
        None
    '''
    plt.figure(figsize=(8, 3), dpi=300).set_facecolor("white")
    plt.plot(losses)
    plt.xlabel("iters")
    plt.ylabel("loss")
    plt.title("Convergence of SVI for " + gRNA)
    plt.savefig(output_dir+"loss_plots/loss_"+gRNA+".png", bbox_inches="tight")
    plt.close()
    
    
def plot_fitted_model(data, weights, locs, scales, threshold, gRNA, output_dir):
    '''
    Saves a plot of the data histogram and the fitted mixture model
    Args:
        data: (tensor) observed transformed gRNA counts
        weights: (np array) estimated proportion for the two Normal component
        locs: (np array) MAP estimate for the mean of the Normal distributions
        scales: (np array) MAP estimate for the scale of the Normal distributions
        threshold: (float) threshold for assigning it to the Normal component
        gRNA: (str) name of the gRNA used for the plot title
        output_dir: (str) name of the output directory
    Returns:
        None
    '''
    X = np.arange(0, max(data)+1, 0.01)
    Y1 = weights[0] * stats.norm.pdf(X, locs[0], scales)
    Y2 = weights[1] * stats.norm.pdf(X, locs[1], scales)

    fig, ax = plt.subplots(figsize=(8, 3), dpi=300)
    sns.histplot(data, binwidth=0.1, color='grey', stat = "proportion")
    ax.plot(X, Y1, "r-", label = "Normal 1")
    ax.plot(X, Y2, "b-", label = "Normal 2")
    ax.plot(X, Y1 + Y2, "k--", label = "Mixture model")
    ax.set_ylim(0, 1)
    ax.axvline(threshold, c = "green", label = "Threshold")
    plt.legend()
    plt.title("Gaussian mixture model for " + gRNA)
    plt.ylabel("Probability Density")
    plt.xlabel("Log10 " + gRNA + " UMI counts")
    plt.savefig(output_dir + "fitted_model_plots/fitted_model_" + gRNA + ".png", bbox_inches="tight")
    plt.close()
    
    
def prob_normal_component(X, weights, locs, scales):
    '''
    Calculates the probability for belonging to the Gaussian component given observations
    Args:
        X: (list) list of values for which the probability is calculated
        weights: (np array) estimated proportion for the two Normal component
        locs: (np array) MAP estimate for the mean of the Normal distributions
        scales: (np array) MAP estimate for the scale of the Normal distributions
    Returns:
        List of probabilities
    '''
    # identify which component has the higher mean
    if locs[0] > locs[1]:
        h = 0
        l = 1
    else:
        h = 1
        l = 0
      
    #pyro dist 
    high = dist.Normal(torch.tensor(locs[h]), torch.tensor(scales)).log_prob(torch.tensor(X)) + np.log(weights[h])
    low = dist.Normal(torch.tensor(locs[l]), torch.tensor(scales)).log_prob(torch.tensor(X)) + np.log(weights[l])
    
    prob = high > low
    return prob


def fit_GMM(gRNA, adata_crispr, output_dir, seed, n_iter, nonzero):
    '''
    Fits Gaussian mixture model for log2 of non-zero UMI counts of one gRNA using variational inference
    Args:
        gRNA: (str) name of the gRNA
        adata_crispr: (AnnData) anndata object with UMI counts of CRISPR Guide Capture
        output_dir: (str) directory in which the resulting plots will be saved
        seed: (int) seed used for pyro
        n_iter (int, optional): number of steps for training the model
        nonzero (bool): if true fit the GMM on the nonzero values only
    Returns:
        List of cells perturbed with the specified gRNA, as well as the inferred threshold
    '''
    # Set optimizer and elbo parameters
    optim = pyro.optim.Adam({"lr": 0.01})
    elbo = TraceEnum_ELBO(num_particles = 1, max_plate_nesting=1)
    
    # Data used to fit the model: log10 transformation of UMI counts for a given gRNA 
    selected_guide = adata_crispr[:,[gRNA]].X
    data = selected_guide.toarray() 
    data = torch.tensor(np.log10(data + 1)).reshape(-1).float() 
    if nonzero:
        data = data[data!=0]
    
    # Only fit model for gRNAs with non-zero counts in at least 2 cells and with a maximum count of at least 2
    if len(data) < 2:   
        print(gRNA + " has only " + str(len(data)) + " cells with non-zero counts, so no model is fitted for that gRNA")
        return([], 0, 0, 0)
    if max(data) < np.log10(2 + 1):
        print("Max log10 UMI count for " + gRNA + " is " + str(max(data)) + ", so no model is fitted for that gRNA")
        return([], 0, 0, 0)
      
    # Choose the best among 10 random initializations.
    loss, seed = min((initialize(seed, optim, elbo, data), seed) for seed in range(10))

    # Initialization of SVI
    initialize(seed, optim, elbo, data)

    # Train the model n_iter steps 
    losses = []
    min_loss = 1.e8
    last_step = 0

    for step in range(n_iter):
        loss = svi.step(data)
        if loss < min_loss - 0.001:
            min_loss = loss
            last_step = step
        losses.append(loss)
 
    # MAP estimates of the model
    map_estimates = global_guide(data)

    weights = map_estimates["weights"].data.numpy()
    locs = map_estimates["locs"].data.numpy()
    scales = map_estimates["scales"].data.numpy()
    estimates = pd.DataFrame({'gRNA': [gRNA], 
                              'weight_Normal1': [weights[0]], 
                              'weight_Normal2': [weights[1]], 
                              'mu1': [locs[0]], 
                              'mu2': [locs[1]], 
                              'scale1': [scales],
                              'scale2': [scales]})

    # create plot of the loss
    plot_loss(losses, gRNA, output_dir)
    
    # threshold for which probability is higher to belong to the higher normal component
    X = np.arange(1, max(selected_guide.toarray()) + 1, 1)
    log_X = np.log10(X + 1)
    df = pd.DataFrame({'t': X, 'prob_normal_component': prob_normal_component(log_X, weights, locs, scales)})
    threshold = df.loc[(df.prob_normal_component == True), 't'].min()

    # create plot of the mixture distribution
    plot_fitted_model(data, weights, locs, scales, np.log10(threshold+1), gRNA, output_dir)
    
    # get cells with gRNA counts above the threshold
    perturbed_cells = adata_crispr.obs_names[selected_guide.toarray().reshape(-1) >= threshold].tolist()
    return(perturbed_cells, threshold, losses[-1], estimates)


def fit_em(gRNA, adata_crispr, nonzero):
    '''
    Fits Gaussian mixture model for log2 of non-zero UMI counts of one gRNA using an EM algorithm
    Args:
        gRNA: (str) name of the gRNA
        adata_crispr: (AnnData) anndata object with UMI counts of CRISPR Guide Capture
        seed: (int) seed used for pyro
        nonzero (bool): if true fit the GMM on the nonzero values only
    Returns:
        Data frame of cells perturbed with the specified gRNA
    '''
    # Select data for one gRNA
    selected_guide = adata_crispr[:,[gRNA]].X
    data = selected_guide.toarray() 
    data = np.log10(data + 1)#.reshape(-1).float() 
    if nonzero:
        data = data[data!=0]
        
    # Fit gaussian mixture model
    gmm = mixture.GaussianMixture(n_components = 2, n_init = 10, covariance_type = "tied", random_state = 0)
    gmm.fit(data)

    # Calculate the threshold
    posterior = gmm.predict_proba(data)
    high_component = np.argmax(gmm.means_)
    
    # Get perturbed cells
    perturbed = posterior[:, high_component] > 0.5
    perturbed_cells = adata_crispr.obs_names[perturbed].tolist()
    perturbations = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA})
    return perturbations

    
def ga_gauss(input_file, output_dir, start_gRNA = 0, step = None, batch_list = None, UMI_threshold = 0,
                  n_iter = 250, nonzero = False, inference = "vi"):
    '''
    Guide assignment in which a Gaussian mixture model is fitted to the log-transformed UMI counts similar to the approach used in Cell Ranger. Two different inference methods are provided that can be selected with the 'inference' parameter.  
    
    Args:
        input_file (str): path to the stored anndata object with the gRNA counts
        output_dir (str): directory in which to store the resulting assignment
        start_gRNA (int, optional): index of the start gRNA when parallelizing assignment for gRNA sets
        step (int, optional): number of gRNAs for which the assignment is done (if set to None, assignment for all gRNAs in the data)
        batch_list (list, optional): list of batches for which to fit the mixture model. If none, mixture model is fited for all batches
        UMI_threshold (int, optional): Additional UMI threshold for assigned cells which is applied after creating the initial assignment to remove cells with fewer UMI counts than this threshold (default: no additional UMI threshold)
        n_iter (int, optional): number of steps for training the model
        nonzero (bool, optional): if True fit the mixture model on the nonzero values only, otherwise all values are used
        inference (str): choice of the inference method, either "vi" (default) for variational inference via pyro or "em" for using an EM algorithm
    
    Returns:
        None
    ''' 
    print('Guide assignment using a Gaussian mixture model per batch')
    # Load gRNA counts data
    print('Load gRNA counts')
    adata_crispr = sc.read_h5ad(input_file)
    
    if batch_list == None:
        batch_list = adata_crispr.obs['batch'].unique()
    
    for batch in batch_list:
        print('Fit Gaussian Mixture models for batch ' + str(batch))
        # If output_dir doesn't exist, the output folders are created
        if not os.path.exists(output_dir + 'batch' + str(batch) + '/'):
            os.makedirs(output_dir + 'batch' + str(batch) + '/')
            os.makedirs(output_dir + 'batch' + str(batch) + "/fitted_model_plots/")
            os.makedirs(output_dir + 'batch' + str(batch) + "/loss_plots/")
            #print("The output directory " + output_dir + 'batch' + str(batch) + '/' +  " was created")

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
        thresholds = pd.DataFrame({'gRNA': [], 'threshold': []})
        losses = pd.DataFrame({'gRNA': [], 'loss': []})
        estimates = pd.DataFrame()

        for gRNA in gRNA_list:
            if inference == "vi":
                perturbed_cells, threshold, loss, map_estimates = fit_GMM(gRNA, adata_crispr_batch, 
                                                                          output_dir + 'batch' + str(batch) + '/', 
                                                                          2024, n_iter, nonzero)
                if len(perturbed_cells) != 0:
                    df = pd.DataFrame({'cell': perturbed_cells, 'gRNA': gRNA})
                    thresholds = pd.concat([thresholds, pd.DataFrame({'gRNA': [gRNA], 'threshold': [threshold]})])
                    losses = pd.concat([losses, pd.DataFrame({'gRNA': [gRNA], 'loss': [loss]})])
                    estimates = pd.concat([estimates, map_estimates])
            elif inference == "em":
                df = fit_em(gRNA, adata_crispr_batch, nonzero)
                
            # get UMI_counts of assigned cells
            UMI_counts = adata_crispr_batch[df['cell'], [gRNA]].X.toarray().reshape(-1)
            df['UMI_counts'] = UMI_counts
            perturbations = pd.concat([perturbations, df], ignore_index = True)
        
        # Optional filtering to assigned cells that have at least 'UMI_threshold' counts
        if perturbations.shape[0] != 0:
            perturbations = perturbations[perturbations['UMI_counts'] >= UMI_threshold]
      
        # Save data frame with the perturbations assigned to each cell
        perturbations.to_csv(output_dir + 'batch' + str(batch) + '/assignments.csv', index = False)
        if inference == "vi":
            thresholds.to_csv(output_dir + 'batch' + str(batch) + '/gRNA_thresholds.csv', index = False)
            losses.to_csv(output_dir + 'batch' + str(batch) + '/gRNA_losses.csv', index = False)
            estimates.to_csv(output_dir + 'batch' + str(batch) + '/estimated_parameters.csv', index = False)

    # Combine the outputs of all batches
    if step == None:
        pert = pd.DataFrame()
        for batch in batch_list:
            b = pd.read_csv(output_dir + 'batch' + str(batch) + '/assignments.csv')
            pert = pd.concat([pert, b])

        pert.to_csv(output_dir + 'assignments.csv', index = False)
    print('Done: outputs are saved in ' + output_dir)
       
    
if __name__ == "__main__":
    pass