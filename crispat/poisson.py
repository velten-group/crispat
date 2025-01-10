import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import torch
from torch.distributions import constraints
from scipy import sparse, stats, special
from scipy.stats import poisson
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from tqdm import tqdm
import time
from dask.distributed import Client, LocalCluster
import dask.bag as db
from functools import partial

import pyro
import pyro.distributions as dist
from pyro import poutine
from pyro.infer.autoguide import AutoDelta
from pyro.optim import Adam
from pyro.infer import SVI, TraceEnum_ELBO, config_enumerate
  

@config_enumerate
def model(data, n_batches):
    '''
    Poisson Mixture Model 
    '''
    # Global variables for each gRNA
    beta0 = pyro.sample("beta0", dist.Normal(0.0, 3.0))
    beta1 = pyro.sample("beta1", dist.LogNormal(1.0, 2.0))
    gamma = pyro.sample("gamma", dist.MultivariateNormal(torch.full((n_batches,), 0.0),
                                                         torch.eye(n_batches) * 1.0))
    pi = pyro.sample("pi", dist.Beta(1.0, 10.0))
    
    # Get data and confounders
    batch = torch.tensor(data.obs['batch'].iloc[:].values)
    seq_depth = torch.tensor(data.obs['total_counts'].iloc[:].values)
    obs = torch.tensor(data.X.toarray()).reshape(-1)

    with pyro.plate("cells", len(data)):
        # Local variables for each cell
        perturbation = pyro.sample("perturbation", dist.Bernoulli(pi)) 
        mu = torch.exp(beta0 + beta1 * perturbation + gamma[batch - 1] + torch.log(seq_depth + 1))
        pyro.sample("obs", dist.Poisson(mu), obs=obs)
        
        
def init_loc_fn(site, n_batches):
    '''
    Define initial parameter values
    '''
    if site["name"] == "beta0":
        return torch.tensor([0.0])
    if site["name"] == "beta1":
        return torch.tensor([3.0])
    if site["name"] == "gamma":
        return torch.zeros(n_batches)
    if site["name"] == "pi":
        return torch.tensor([0.01])
    raise ValueError(site["name"])

    
def initialize(seed, optim, elbo, data, n_batches):
    '''
    Initialization for SVI 
    
    Args:
        seed (str): seed that is used in pyro
        optim: pyro optimizer
        elbo: pyro loss function
        data (tensor): observed transformed gRNA counts
        conditioned_model (pyro model): conditioned pyro mixture model
        n_batches (int): number of batches
    
    Returns:
        Initial loss
    '''
    global global_guide, svi
    pyro.set_rng_seed(seed)
    pyro.clear_param_store()
    global_guide = AutoDelta(
        poutine.block(model, expose=["beta0", "beta1", "gamma", "pi", "overdispersion"]),
        init_loc_fn=lambda site: init_loc_fn(site, n_batches)
    )
    svi = SVI(model, global_guide, optim, loss=elbo)
    return svi.loss(model, global_guide, data, n_batches)


def plot_loss(losses, gRNA, output_dir):
    '''
    Saves a plot of the loss over the SVI steps
    
    Args:
        losses (list): loss over the SVI steps
        gRNA (str): name of the gRNA used for the plot title
        output_dir (str): name of the output directory
    Returns:
        None
    '''
    plt.figure(figsize=(8, 3), dpi=300).set_facecolor("white")
    plt.plot(losses)
    plt.xlabel("iters")
    plt.ylabel("loss")
    plt.yscale("log")
    plt.title("Convergence of SVI for " + gRNA)
    plt.savefig(output_dir+"loss_plots/loss_"+gRNA+".png", bbox_inches="tight")
    plt.close()
    

def fit_SCEPTRE(gRNA, adata_crispr, batch_list, output_dir, seed, n_iter, subsample_size):
    '''
    Fits SCEPTRE model for UMI counts of one gRNA 
    
    Args:
        gRNA (str): name of the gRNA
        adata_crispr (AnnData): anndata object with UMI counts of CRISPR Guide Capture
        batch_list (list): list of all batches
        output_dir (str): directory in which the resulting plots will be saved
        seed (int): seed used for pyro
        n_iter (int): number of steps for training the model
        subsample_size (int): number of cells to use for each step
    
    Returns:
        List of cells perturbed with the specified gRNA, as well as the inferred threshold
    '''
    # Set optimizer and elbo parameters
    optim = pyro.optim.Adam({"lr": 0.01, "betas": [0.8, 0.99]})
    elbo = TraceEnum_ELBO(num_particles = 1, max_plate_nesting=1)
   
    # Data used to fit the model
    selected_gRNA = adata_crispr[:,[gRNA]]

    # Only fit model for gRNAs with non-zero counts in at least 3 cells and with a maximum count of at least 2
    data = selected_gRNA.X.toarray() 
    data = torch.tensor(data[data != 0])
    if len(data) < 2:   
        print(gRNA + " has only " + str(len(data)) + " cells with non-zero counts, so no model is fitted for that gRNA")
        return(0, pd.DataFrame())
    if max(data) < 2:
        print("Max UMI count for " + gRNA + " is " + str(max(data).item()) + ", so no model is fitted for that gRNA")
        return(0, pd.DataFrame())

    # Initialization of SVI
    n_batches = len(batch_list)
    initialize(seed, optim, elbo, selected_gRNA, n_batches)

    # Train the model n_iter steps 
    losses = []
    for step in range(n_iter):
        subsample = sc.pp.subsample(selected_gRNA, n_obs=subsample_size, random_state=2023+step, copy = True)
        loss = svi.step(subsample, n_batches)
        losses.append(loss)
        
    # MAP estimates of the model
    map_estimates = global_guide(data)
    batch_effects = pd.DataFrame({'gRNA': gRNA,
                                  'param': ['gamma_' + str(batch) for batch in batch_list],
                                  'value': map_estimates['gamma'].detach().numpy()})
    beta0 = map_estimates["beta0"].item()
    beta1 = map_estimates["beta1"].item()
    pi = map_estimates["pi"].item()
    estimates = pd.concat([pd.DataFrame({'gRNA': gRNA, 
                                         'param': ['beta0', 'beta1', 'pi'], 
                                         'value': [map_estimates["beta0"].item(), map_estimates["beta1"].item(),
                                                   map_estimates["pi"].item()]}), 
                           batch_effects])

    # create plot of the loss
    plot_loss(losses, gRNA, output_dir)
    
    return(losses[-1], estimates)


def get_perturbed_cells(adata_crispr, estimates, gRNA):
    '''
    Gets perturbed cells for one gRNA 
    
    Args:
        adata_crispr (AnnData object): UMI counts of CRISPR gRNAs in all cells
        estimates (pd DataFrame): parameter estimates from the Negative Binomial model
        gRNA (str): gRNA name
        
    Returns:
        Perturbed cells
    '''
    # Get data and confounders
    selected_gRNA = adata_crispr[:,[gRNA]]
    data = pd.DataFrame({'UMI_counts': selected_gRNA.X.toarray().reshape(-1), 
                         'batch': selected_gRNA.obs['batch'], 
                         'seq_depth': selected_gRNA.obs['total_counts']})
    data = data[data['UMI_counts'] != 0].copy()
    
    # get inferred parameters
    beta0 = estimates.loc[estimates['param'] == 'beta0', 'value'].item()
    beta1 = estimates.loc[estimates['param'] == 'beta1', 'value'].item()
    gamma = list(estimates.loc[estimates['param'].str.startswith('gamma_'), 'value'])
    gamma_cells = np.array([gamma[batch - 1] for batch in data['batch']])
    pi = estimates.loc[estimates['param'] == 'pi', 'value'].item()
    
    # calculate mean per cell
    data['mu0'] = np.exp(beta0 + gamma_cells + np.log(data['seq_depth']))
    data['mu1'] = np.exp(beta0 + beta1 + gamma_cells + np.log(data['seq_depth']))
    
    # get probabilities for the two mixture components
    data['p0'] = poisson.pmf(data['UMI_counts'], data['mu0']) * (1 - pi)
    data['p1'] = poisson.pmf(data['UMI_counts'], data['mu1']) * pi
    data['pert_probability'] = data['p1'] / (data['p0'] + data['p1'])
    data['perturbation'] =  np.where(data['pert_probability'] >= 0.8, 1, 0)
    
    # filter for perturbed cells
    perturbed_cells = data[data['perturbation'] == 1].copy()
    perturbed_cells['gRNA'] = gRNA
    perturbed_cells.index.name = 'cell'
    perturbed_cells = perturbed_cells.reset_index()
    perturbed_cells = perturbed_cells[['cell', 'gRNA', 'UMI_counts']]

    return perturbed_cells


def parallel_assignment(gRNA, adata_crispr, batch_list, output_dir, seed, n_iter, subsample_size):
    '''
    Call assignment function for one gRNA 
    
    Args:
        gRNA (str): name of the gRNA
        adata_crispr (AnnData): anndata object with UMI counts of CRISPR Guide Capture
        batch_list (list): list of all batches
        output_dir (str): directory in which the resulting plots will be saved
        seed (int): seed used for pyro
        n_iter (int): number of steps for training the model
        subsample_size (int): number of cells to use for each step
    
    Returns:
        gRNA, loss, MAP_estimates, perturbed_cells
    '''
    loss, map_estimates = fit_SCEPTRE(gRNA, adata_crispr, batch_list, output_dir, seed, n_iter, subsample_size)
    if map_estimates.shape[0] != 0:
        perturbed_cells = get_perturbed_cells(adata_crispr, map_estimates, gRNA)
    else:
        perturbed_cells = pd.DataFrame()
    return gRNA, loss, map_estimates, perturbed_cells
            
        
def ga_poisson(input_file, output_dir, start_gRNA = 0, gRNA_step = None, batch_list = None, UMI_threshold = 0, n_iter = 2500, 
               subsample_size = 15000, parallelize = True, n_processes = None, mem_limit = '10GB'):
    '''
    Guide assignment with a Poisson mixture model based on SCEPTRE approach
    
    Args:
        input_file (str): path to the stored anndata object with the gRNA counts
        output_dir (str): directory in which to store the resulting assignment
        start_gRNA (int, optional): index of the start gRNA when parallelizing assignment for gRNA sets
        gRNA_step (int, optional): number of gRNAs for which the assignment is done (if set to None, assignment for all gRNAs in the data)
        batch_list (list, optional): list of batches for which to fit the mixture model. If none (default), all available batches are used. 
        UMI_threshold (int, optional): Additional UMI threshold for assigned cells which is applied after creating the initial assignment to remove cells with fewer UMI counts than this threshold (default: no additional UMI threshold)
        n_iter (int, optional): number of steps for training the model
        subsample_size (int, optional): number of cells to use for each step
        parallelize (bool, optional): whether to parallelize the computation over the gRNA (default = True)
        n_processes (int, optional): specifies number of processes to use for parallelization if parallelize = True. If set to None (default), all available CPUs will be used (if this number is not higher than the number of gRNAs). 
        mem_limit (str, optional): set memory limit for the dask cluster (default: 10GB)

    Returns:
        None
    '''   
    print('Guide assignment based on SCEPTRE Poisson mixture model')
    # If output_dir doesn't exist, the output folders are created
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs(output_dir + "loss_plots/", exist_ok=True)
        print("The output directory " + output_dir +  " was created")

    # Load gRNA counts data
    print('Load gRNA counts')
    adata_crispr = sc.read_h5ad(input_file)
    # subset to specified batches
    if batch_list != None:
        adata_crispr = adata_crispr[adata_crispr.obs['batch'].isin(batch_list)].copy()
    else:
        batch_list = adata_crispr.obs['batch'].unique()
    
    # Add total number of gRNA counts per cell 
    adata_crispr.obs['total_counts'] = np.array(adata_crispr.X.sum(axis = 1)).flatten()

    gRNA_list = adata_crispr.var_names.tolist()
    
    if gRNA_step != None:
        end_gRNA = start_gRNA + gRNA_step - 1
        if end_gRNA >= len(gRNA_list):
            end_gRNA = len(gRNA_list) - 1
        gRNA_list = gRNA_list[start_gRNA:(end_gRNA + 1)]

    # Fit SCEPTRE Model for each gRNA
    perturbations = pd.DataFrame({})
    losses = pd.DataFrame({'gRNA': [], 'loss': []})
    estimates = pd.DataFrame()
    print('Fit SCEPTRE Poisson Model for each gRNA')
    
    # Fit SCEPTRE model in parallel to all gRNAs
    if parallelize:  
        if n_processes == None:
            n_processes = os.cpu_count() # if not specified all available CPUs will be used
        if n_processes > len(gRNA_list):
            n_processes = len(gRNA_list)
        print(str(n_processes) + ' parallel processes used')
          
        partial_function = partial(parallel_assignment, adata_crispr = adata_crispr, batch_list = batch_list, 
                                   output_dir = output_dir, seed = 2024, n_iter = n_iter, subsample_size = subsample_size)
        
        # start a local Dask cluster
        cluster = LocalCluster(n_workers = n_processes, threads_per_worker = 1, memory_limit = mem_limit)
        client = Client(cluster, heartbeat_interval='5s', timeout='30s')

        # run in parallel over gRNAs
        bag = db.from_sequence(gRNA_list)
        results = bag.map(partial_function).compute()
        client.close()
        cluster.close()
        
        # combine the results per gRNA
        for gRNA, loss, map_estimates, perturbed_cells in results:
            losses = pd.concat([losses, pd.DataFrame({'gRNA': [gRNA], 'loss': [loss]})])
            estimates = pd.concat([estimates, map_estimates])
            perturbations = pd.concat([perturbations, perturbed_cells])
    
    # Fit SCEPTRE model without parallelization
    else:
        for gRNA in tqdm(gRNA_list):
            time.sleep(0.01)
            # fit the model to infer the parameters
            loss, map_estimates = fit_SCEPTRE(gRNA, adata_crispr, batch_list, output_dir, 2024, n_iter, subsample_size)
            losses = pd.concat([losses, pd.DataFrame({'gRNA': [gRNA], 'loss': [loss]})])
            estimates = pd.concat([estimates, map_estimates])

            # get the perturbed cells
            if map_estimates.shape[0] != 0:
                perturbed_cells = get_perturbed_cells(adata_crispr, map_estimates, gRNA)
            else:
                perturbed_cells = pd.DataFrame()
            perturbations = pd.concat([perturbations, perturbed_cells])
    
    # Optional filtering to assigned cells that have at least 'UMI_threshold' counts
    if perturbations.shape[0] != 0:
        perturbations = perturbations[perturbations['UMI_counts'] >= UMI_threshold]
    
    # Save data frames with the results
    if gRNA_step == None:
        perturbations.to_csv(output_dir + 'assignments.csv', index = False)
        losses.to_csv(output_dir + 'gRNA_losses.csv', index = False)
        estimates.to_csv(output_dir + 'gRNA_estimates.csv', index = False)
    else:
        perturbations.to_csv(output_dir + 'assignments_'+str(start_gRNA)+'-'+str(end_gRNA)+'.csv', index = False)
        losses.to_csv(output_dir + 'gRNA_losses_'+str(start_gRNA)+'-'+str(end_gRNA)+'.csv', index = False)
        estimates.to_csv(output_dir + 'gRNA_estimates_'+str(start_gRNA)+'-'+str(end_gRNA)+'.csv', index = False)
    
    print('Done: outputs are saved in ' + output_dir)


if __name__ == "__main__":
    pass
