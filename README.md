# crispat: CRISPR guide Assignment Toolkit

Pooled single-cell CRISPR screens are a powerful tool for systematically gaining new insights into the functional consequences of genetic perturbations in high-throughput analyses. To allow for meaningful downstream analyses and biological insights about gene regulatory mechanisms from single-cell CRISPR screen experiments a first crucial step is the assignment of cells to specific perturbations that are defined during guide assignment. In guide assignment, thresholds on the measured gRNA counts per cell are used to distinguish between background noise and the actual perturbation state of a cell. However, lots of different guide assignment strategies and thresholds are used by different labs without any guidance on what model or threshold to choose when. Focusing on **low MOI** CRISPRi screens we systematically compared various methods for guide assignment in terms of the number of assigned cells, the effectiveness of target gene downregulation in assigned cells, as well as the number of false and total discoveries. The methods in this benchmark include simple approaches such as a threshold on the UMI counts or assigning the gRNA with highest counts per cell, as well as more advanced models taking into account the variability per cell, the variability per gRNA, or both. Even though there is a high overlap for the assigned cells across methods, the total number of assigned cells varies strongly. Furthermore, we observed that it is crucial to determine reasonable dataset dependent thresholds per method since there is a trade-off between a gain in statistical power and a less effective average target gene downregulation when assigning more cells. Therefore, we provide our code as a tool to easily compare various guide assignment approaches and find a suitable assignment for new data sets. 

crispat contains 11 guide assignment methods, which are grouped into 4 main categories based on the information that is used during assignment:
- independent (`ga_UMI`)
- across gRNAs (`ga_max`, `ga_ratio`, `ga_2beta`, `ga_3beta`)
- across cells (`ga_cellranger`, `ga_replogle`)
- across gRNAs and across cells (`ga_SCEPTRE`, `ga_negative_binomial`, `ga_binomial`, `ga_quantiles`)

In addition to the guide assignment functions, the package includes some additional helper functions. First, the package provides two functions for importing the data starting from either a csv file containing the gRNA count matrix or from the cellranger count output. Next, the `combine_assignments` functions allow to create a combined csv file when running some of the methods in parallel on subsets of the data and `load_assignments` creates a pandas DataFrame of the output from all specified methods. Finally, the function `plot_intersection_heatmap` creates a heatmap with the pair-wise intersection score (intersection of the assignments / N assignments in method of the row), as well as a barplot showing the total number of cells with a single gRNA assigned. An example for the usage of all functions in our python package is shown for `example_data` in the `guide_assignment.ipynb` script. Tutorials on how we performed our downstream analyses in R using the SCEPTRE package and comparing the outputs across methods, can be found in the `tutorials` directory. 
