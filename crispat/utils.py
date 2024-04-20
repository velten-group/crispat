import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def combine_assignments(data_dir):
    '''
    Combines assignment files for methods that have been performed on gRNA subsets
    
    Args:
        data_dir (str): directory containing the assignment output
    Returns:
        None
    '''
    # check which perturbation files are in the directory
    all_files = os.listdir(data_dir)
    csv_files = [file for file in all_files if file.startswith('perturbations_') and file.endswith('.csv')]

    perturbations = pd.DataFrame()
    
    # read in data for all gRNA subsets
    for file in csv_files:
        perturbations_g = pd.read_csv(data_dir + file)
        perturbations = pd.concat([perturbations, perturbations_g], ignore_index = True)
    
    # remove duplicates (in case there were overlapping indices)
    perturbations = perturbations.drop_duplicates()
    
    # save file containing all gRNAs
    perturbations.to_csv(data_dir + "perturbations.csv")
    
    
def load_assignments(gc_dict, data_dir):
    '''
    Loads and combines assignment files of specified methods
    
    Args:
        data_dir (str): directory containing the assignment output
        start_gRNAs (list): list of gRNA start indices
        end_gRNAs (list): list of gRNA end indices
        
    Returns:
        A pd DataFrame containing the gRNA-cell assignments of each method
    '''
    perturbations = pd.DataFrame()
    for name, method_params in gc_dict.items():
        method_dir, t = method_params

        if t == None: 
            assignment = pd.read_csv(data_dir + method_dir + '/perturbations.csv')[['cell', 'gRNA']]
        else:
            assignment = pd.read_csv(data_dir + method_dir + '/perturbations_t' + str(t) + '.csv')[['cell', 'gRNA']]
        assignment['method'] = name

        # Combine with all others
        perturbations = pd.concat([perturbations, assignment])
    return perturbations


def get_intersections(pert_dict):
    '''
    Calculates pairwise intersection proportions
    
    Args:
        pert_dict (dictionary): dictionary which contains assigned cell_gRNAs for each method
        
    Returns:
        A pd DataFrame with the pairwise similarity scores
    '''
    matrix = np.zeros((len(pert_dict), len(pert_dict)))
    results = pd.DataFrame(matrix, index=pert_dict.keys(), columns=pert_dict.keys())
    
    # loop over all pair of methods
    for i in list(pert_dict.keys()):
        for j in list(pert_dict.keys()):
            # calculate number of intersecting assignments divided by the number of assignments in set 1
            set1 = set(pert_dict[i])
            set2 = set(pert_dict[j])    
            intersection_size = len(set1.intersection(set2))
            set1_size = len(set1)
            similarity = intersection_size / set1_size
            results.loc[i, j] = similarity
    return results


def plot_intersection_heatmap(perturbations, colors = None):
    '''
    Plot a heatmap with intersection proportions and a barplot with the number of assignments per method
    
    Args:
        perturbations (pd DataFrame): df with the assigned perturbations (needed columns: method, cell, gRNA)
        colors (dictionary, optional): specifies the colors to use for each method in the barplot
        
    Returns:
        A matplotlib plot consisting of two subfigures (intersection heatmap and number of assignments barplot)
    '''
    
    # create dictionary with the cell-perturbation pairs per method
    pert_dict = {}
    methods = perturbations['method'].unique()
    for method in methods:
        subset = perturbations[perturbations['method'] == method]
        cell_gRNA_list = subset['cell'] + '_' + subset['gRNA'] 
        pert_dict[method] = cell_gRNA_list.unique()
    
    # calculate matrix with intersection / assignment in row
    matrix = get_intersections(pert_dict)
    intersections = pd.DataFrame(matrix, columns = methods, index = methods)
    
    # calculate total number of assignments per method
    n_assignments_per_method = perturbations[['cell', 'method', 
                                         'gRNA']].drop_duplicates().groupby(['method']).size().reset_index(name='count')
    n_assignments_per_method = n_assignments_per_method.sort_values(by='count')
    # sort the intersections according to the total number of assignments
    order = list(n_assignments_per_method['method'])
    intersections = intersections.loc[order, order]
    
    # plot heatmap
    fig, axs = plt.subplots(figsize = (14, 5), nrows = 1, ncols = 2, sharey = True, width_ratios = (0.6, 0.4))
    cax = axs[0].matshow(intersections, cmap='YlOrRd')
    axs[0].set_xticks(range(len(methods)))
    axs[0].set_yticks(range(len(methods)))
    axs[0].set_xticklabels(order, rotation=45, ha='left')
    axs[0].set_yticklabels(order)
    plt.colorbar(cax, label='Intersection / N Row')
    
    # plot barplot with set size on the right
    if colors == None:
        colors = sns.color_palette("husl", n_colors = n_assignments_per_method.shape[0])
    
    sns.barplot(data = n_assignments_per_method, y = "method", x = "count", ax = axs[1],
               palette = colors)
    axs[1].set_xlabel('')
    axs[1].set_ylabel('')
    axs[1].set_title('Number of cells\nwith single assignment')
    