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
        gc_dict (dict): dictionary with method directory and threshold for each assignment
        data_dir (str): directory containing all assignment output folders
        
    Returns:
        A pd DataFrame containing the gRNA-cell assignments of each method
    '''
    perturbations = pd.DataFrame()
    for name, method_params in gc_dict.items():
        method_dir, t = method_params

        if t == None: 
            assignment = pd.read_csv(data_dir + method_dir + '/assignments.csv')[['cell', 'gRNA']]
        else:
            assignment = pd.read_csv(data_dir + method_dir + '/assignments_t' + str(t) + '.csv')[['cell', 'gRNA']]
        assignment['method'] = name

        # Combine with all others
        perturbations = pd.concat([perturbations, assignment])
    return perturbations


def calculate_jaccard(pert_dict):
    '''
    Calculates pairwise Jaccard index
    
    Args:
        pert_dict (dictionary): dictionary which contains assigned cell_gRNAs for each method
        
    Returns:
        A pd DataFrame with the pairwise similarity scores
    '''
    matrix = np.zeros((len(pert_dict), len(pert_dict)))
    results = pd.DataFrame(matrix, index=pert_dict.keys(), columns=pert_dict.keys())

    for i in list(pert_dict.keys()):
        for j in list(pert_dict.keys()):
            set1 = set(pert_dict[i])
            set2 = set(pert_dict[j])    
            intersection_size = len(set1.intersection(set2))
            union_size = len(set1.union(set2))
            similarity = intersection_size / union_size
            results.loc[i, j] = similarity
    return results


def plot_intersection_heatmap(perturbations, method_order = None):
    '''
    Plot a heatmap with Jaccard index showing the intersecting assignments
    
    Args:
        perturbations (pd DataFrame): df with the assigned perturbations (needed columns: method, cell, gRNA)
        method_order (list): list defining the order of the rows and columns (default: alphabetic order)
        
    Returns:
        A matplotlib plot 
    '''
    # create dictionary with the cell-perturbation pairs per method
    pert_dict = {}
    
    if method_order != None:
        methods = method_order
    else:
        methods = perturbations['method'].unique()
        
    for method in methods:
        subset = perturbations[perturbations['method'] == method]
        cell_gRNA_list = subset['cell'] + '_' + subset['gRNA'] 
        pert_dict[method] = cell_gRNA_list.unique()
    
    # calculate matrix with intersection / total in row
    matrix = calculate_jaccard(pert_dict)
    intersections = pd.DataFrame(matrix, columns = methods, index = methods)
    
    # plot heatmap
    fig, axs = plt.subplots(figsize = (6, 3), nrows = 1, ncols = 1)
    cax = axs.matshow(intersections, cmap='YlOrRd')
    axs.set_xticks(range(len(methods)))
    axs.set_yticks(range(len(methods)))
    axs.set_xticklabels(methods, rotation=45, ha='left')
    axs.set_yticklabels(methods)
    plt.colorbar(cax, label='Jaccard index')
      
    
def plot_n_assigned_cells(perturbations, colors = None):
    '''
    Plots a barplot with the number of assigned cells and uniquely assigned cells per method
    
    Args:
        perturbations (pd DataFrame): df with the assigned perturbations (needed columns: method, cell, gRNA)
        colors (dictionary, optional): specifies the colors to use for each method in the barplot 
        
    Returns:
        A matplotlib plot 
    '''
    perturbations = perturbations[['cell', 'method', 'gRNA']].drop_duplicates()
    
    # calculate number of assignments per method
    n_total_assignments = perturbations[['method','cell']].drop_duplicates().groupby(['method']).size().reset_index(name='count')
    n_total_assignments['group'] = 'All assigned cells'

    # calculate number of cells with single gRNA assigned
    grouped_df = perturbations.groupby(['method', 'cell']).size().reset_index(name='grna_count')
    filtered_df = perturbations.merge(grouped_df, on=['method', 'cell'])
    filtered_df = filtered_df[filtered_df['grna_count'] == 1]
    n_single_assignments = filtered_df.groupby(['method']).size().reset_index(name='count')
    n_single_assignments['group'] = 'Uniquely assigned cells'

    combined_df = pd.concat([n_total_assignments, n_single_assignments])
    combined_df['group'] = pd.Categorical(combined_df['group'], categories=['Uniquely assigned cells', 'All assigned cells'], ordered=True)
    combined_df.sort_values(by='group', inplace=True)
    
    if colors == None:
        colors = sns.color_palette("husl", n_colors = n_total_assignments.shape[0])
    else:
        # order by guide assignment groups
        method_order = list(colors.keys())
        combined_df['method'] = pd.Categorical(combined_df['method'], categories=method_order, ordered=True)
        combined_df.sort_values(by='method', inplace=True)
    
    plt.figure(figsize = (3,3))
    scatter = sns.scatterplot(
        data = combined_df, x = 'count', y='method', 
        palette=colors, hue='method', s=50, style='group', 
        markers={'All assigned cells': 'X', 'Uniquely assigned cells': 'o'}
    )
    plt.grid(True, which='both', axis='y', color='grey', linestyle='-', linewidth=0.3)
    plt.grid(False, which='both', axis='x')

    # manually create legend
    handles, labels = scatter.get_legend_handles_labels()
    new_handles = [handles[idx] for idx, label in enumerate(labels) if label in ['All assigned cells', 'Uniquely assigned cells']]
    new_labels = [label for label in labels if label in ['All assigned cells', 'Uniquely assigned cells']]
    scatter.legend(new_handles, new_labels, title='', loc = "upper center", bbox_to_anchor=(0.5, 1.25))

    plt.xlabel('Number of assigned cells')
    plt.ylabel('')
    plt.xlim(0, None)