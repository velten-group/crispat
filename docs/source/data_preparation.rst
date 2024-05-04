Data preparation
================

To apply one of the 11 guide assignment functions of the crispat package, we first have to create and save an AnnData object containing the gRNA counts for all cells. Depending on the format of your input data, this can either be done by using one of our import functions or for more customizability directly using the anndata package. In our package, we provide two functions that either read in cellranger output (`create_anndata_from_cellranger`) or a csv file with the gRNA-cell matrix (`create_anndata_from_csv`) to create the .h5ad file. If creating an AnnData object yourself, make sure to add a batch column to adata.obs even if you only have one batch (use same value for all cells in this case).

.. autofunction:: crispat.create_anndata_from_cellranger
   :noindex:

.. autofunction:: crispat.create_anndata_from_csv
   :noindex: