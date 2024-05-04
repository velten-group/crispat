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
To compare the assignments with each other, crispat contains the function `plot_intersection_heatmap` which creates two figures next to each other. On the left, a heatmap of the pairwise overlap of assignments divided by the number of assignments from the method shown in the row is plotted. On the right, the number of cells with single assignments is shown and the y-axis of both plots is sorted according to an increasing number of assigned cells.
    
.. autofunction:: crispat.plot_intersection_heatmap
    :noindex:

Effects on downstream analysis
------------------------------
To investigate the consequences of the assignment differences for discovery analysis, the assignments from crispat can serve as input to different differential expression tests. For the analysis shown in our paper, we use the crispat output as an input for the R package `SCEPTRE`_, which is tailored to single-cell CRISPR screen analysis and combines multiple analyses and control checks in a statistical rigorous fashion. First, SCEPTRE calculates the log2 fold changes and p-values for the target genes (power check). Next, it calculates the number of false discoveries which are genes that are significantly differentially expressed comparing the cells of one non-targeting gRNA against all other non-targeting control cells (calibration). And finally, it calculates the differentially expressed genes for each perturbation (discovery analysis). In our github `repository`_, we provide a tutorial on how the assignments from crispat can be used as an inut for SCEPTRE. 

.. _SCEPTRE: https://timothy-barry.github.io/sceptre-book/sceptre.html
.. _repository: https://github.com/velten-group/crispat/blob/main/tutorials/downstream_analyses/SCEPTRE_discovery.R