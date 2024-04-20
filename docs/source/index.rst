.. crispat documentation master file, created by
   sphinx-quickstart on Fri Apr  5 16:03:56 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Documentation for crispat 
=========================

Pooled single-cell CRISPR screens are a powerful tool for systematically gaining new insights into the functional consequences of genetic perturbations in high-throughput analyses. To allow for meaningful downstream analyses and biological insights about gene regulatory mechanisms from single-cell CRISPR screen experiments a first crucial step is guide assignment, where cells are assigned to specific guides and corresponding genetic targets. For this, thresholds on the measured gRNA counts per cell are used to distinguish between background contamination and the actual guide presence in a cell. However, lots of different guide assignment strategies and thresholds are used by different labs without any guidance on what model or threshold to choose when.

As demonstrated on low MOI CRISPRi screens in our preprint [Braunger et al, 2024](link to the preprint will be added soon) the choice of guide assignment strategy strongly influences the results, highlighting the need to choose a suitable strategy for the data at hand for a reliable and powerful analysis of the data. To help with this choice the crispat package implements 11 different assignment methods and facilitates their comparison.


.. toctree::
   :maxdepth: 2
   :caption: List of all functions:
   
   crispat


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
