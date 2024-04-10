.. crispat documentation master file, created by
   sphinx-quickstart on Fri Apr  5 16:03:56 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for crispat 
=========================

Pooled single-cell CRISPR screens are a powerful tool for systematically gaining new insights into the functional consequences of genetic perturbations in high-throughput analyses. To allow for meaningful downstream analyses and biological insights about gene regulatory mechanisms from single-cell CRISPR screen experiments a first crucial step is the assignment of cells to specific perturbations that are defined during guide assignment. In guide assignment, thresholds on the measured gRNA counts per cell are used to distinguish between background noise and the actual perturbation state of a cell. However, lots of different guide assignment strategies and thresholds are used by different labs without any guidance on what model or threshold to choose when. Focusing on low MOI CRISPRi screens we systematically compared various methods for guide assignment in terms of the number of assigned cells, the effectiveness of target gene downregulation in assigned cells, as well as the number of false and total discoveries. The methods in this benchmark include simple approaches such as a threshold on the UMI counts or assigning the gRNA with highest counts per cell, as well as more advanced models taking into account the variability per cell, the variability per gRNA, or both. Even though there is a high overlap for the assigned cells across methods, the total number of assigned cells varies strongly. Furthermore, we observed that it is crucial to determine reasonable dataset dependent thresholds per method since there is a trade-off between a gain in statistical power and a less effective average target gene downregulation when assigning more cells. Therefore, we provide our CRISPR guide Assignment Tool (**crispat**) as a python package which allows easy comparison of various guide assignment approaches to find a suitable assignment for low MOI CRISPR screen data sets. 


.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   crispat


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
