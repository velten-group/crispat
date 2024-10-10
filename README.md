<p align="center">
    <img src="https://github.com/velten-group/crispat/blob/main/crispat_logo.png" alt="logo" width="200"/>
</p>

# crispat: CRISPR guide assignment tool

Pooled single-cell CRISPR screens are a powerful tool for systematically gaining new insights into the functional consequences of genetic perturbations in high-throughput analyses. To allow for meaningful downstream analyses and biological insights about gene regulatory mechanisms from single-cell CRISPR screen experiments a first crucial step is guide assignment, where cells are assigned to specific guides and corresponding genetic targets. For this, thresholds on the measured gRNA counts per cell are used to distinguish between background contamination and the actual guide presence in a cell. However, lots of different guide assignment strategies and thresholds are used by different labs without any guidance on what model or threshold to choose when. 

As demonstrated on *low MOI CRISPRi* screens in our paper [Braunger et al, 2024](https://academic.oup.com/bioinformatics/article/40/9/btae535/7750392) the choice of guide assignment strategy strongly influences the results, highlighting the need to choose a suitable strategy for the data at hand for a reliable and powerful analysis of the data. To help with this choice the **crispat** package implements 11 different assignment methods and facilitates their comparison. 

## Guide assignment methods
Methods include simple approaches such as a threshold on the UMI counts or assigning the gRNA with highest counts per cell, as well as more advanced models taking into account the variability per cell, the variability per gRNA, or both. Methods are grouped into 4 main categories based on the information that is used during assignment:
- independent (`ga_umi`)
- across gRNAs (`ga_max`, `ga_ratio`)
- across cells (`ga_gauss`, `ga_poisson_gauss`)
- across gRNAs and across cells (`ga_2beta`, `ga_3beta`, `ga_poisson`, `ga_negative_binomial`, `ga_binomial`, `ga_quantiles`)

For details on the individual methods please refer to our our paper [Braunger et al, 2024](https://academic.oup.com/bioinformatics/article/40/9/btae535/7750392).

## Tools for method comparison 
In addition to the guide assignment functions, the package includes some additional helper functions incl. for
* data import (starting from either a csv file containing the gRNA count matrix or from the cellranger count output)
* running and combining results from different guide assignment methods (`combine_assignments`, `load_assignments`)
* visualization and comparison of different methods (`plot_intersection_heatmap`, `plot_n_assigned_cells`, tutorials)

For details on the downstream analysis and comparison on different guide assignments obtained by crispat take a look at our tutorial page on [downstream analyses](tutorials/downstream_analyses).

## Installation
To install crispat, run `pip install crispat` (see [PyPI page](https://pypi.org/project/crispat/)). To avoid dependency conflicts, we recommend to first create a new conda environment with python 3.10. and then using `pip install crispat`. Alternatively, the current version of crispat can also be installed by cloning this repository and then running `python setup.py install`.

## Getting started
An example use case of crispat is shown in the [`guide_assignment.ipynb`](tutorials/guide_assignment.ipynb) script. 

Tutorials on how to evaluate differences between methods for downstream analyses can be found in the [`tutorials/downstream_analyses`](tutorials/downstream_analyses) directory. 

## Documentation

For details on individual functions in the package refer to the [crispat documentation](https://crispat.readthedocs.io/en/latest/index.html).
