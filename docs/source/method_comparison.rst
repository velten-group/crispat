Method comparison
=================

Load assignment results
-----------------------
In case some of the methods were manually run on distinct gRNA subsets, we provide the function `combine_assignments` to combine the outputs into one csv file. After running all guide assignment methods of interest, we can create a data frame which contains all gRNA cell assignments per method using the `load_assignments` function and filter for cells that got exactly one gRNA assigned.  

.. autofunction:: crispat.combine_assignments
    :noindex:
    
.. autofunction:: crispat.load_assignments
    :noindex:
    
Number of assigned cells and intersection
-----------------------------------------
To plot the number of total assigned cells, as well as the number of uniquely assigned cells per method, crispat contains the function `plot_n_assigned_cells`. To compare how similar various assignments are to each other, crispat also includes a function (`plot_intersection_heatmap`) which creates  a heatmap of the pairwise Jaccard index of uniquely assigned cells per method. 
    
.. autofunction:: crispat.plot_n_assigned_cells
    :noindex:

.. autofunction:: crispat.plot_intersection_heatmap
    :noindex:

Effects on downstream analysis
------------------------------
To investigate the consequences of the assignment differences for discovery analysis, the assignments from crispat can serve as input to different differential expression tests. For the analysis shown in our paper, we use the crispat output as an input for the R package `SCEPTRE`_, which is tailored to single-cell CRISPR screen analysis and combines multiple analyses and control checks in a statistical rigorous fashion. First, SCEPTRE calculates the log2 fold changes and p-values for the target genes (power check). Next, it calculates the number of false discoveries which are genes that are significantly differentially expressed comparing the cells of one non-targeting gRNA against all other non-targeting control cells (calibration). And finally, it calculates the differentially expressed genes for each perturbation (discovery analysis). However, it is also possible to input the assignments obtained by crispat into other tools for differential expression testing such as scanpy or Seurat. Therefore, we provide tutorials on how to do downstream analysis with SCEPTRE, scanpy and Seurat in our github `repository`_. 

.. _SCEPTRE: https://timothy-barry.github.io/sceptre-book/sceptre.html
.. _repository: https://github.com/velten-group/crispat/blob/main/tutorials/downstream_analyses/