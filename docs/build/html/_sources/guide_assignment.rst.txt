Guide assignment
================
This package contains implementations of 11 different guide assignment methods which we grouped into four main categories based on whether information is shared across gRNAs, cells or both. Details on all functions can be found in our manuscript (link will be added very soon).

Independent
-----------
**UMI threshold (UMI_t):** The simplest approach is to not share any information across the gRNA-cell matrix and check for each value separately whether it is at least as high than a fixed threshold. To find a suitable thresholds, a list of thresholds can be passed as one argument of the function

.. autofunction:: crispat.ga_umi


Across gRNAs
------------

**Maximum:** This method assigns each cell the gRNA with highest UMI count in this cell.

.. autofunction:: crispat.ga_max

**Relative frequency threshold (Ratio_X%):** This method assigns for each cell the gRNA with highest counts in this cell if its counts comprise at least X% of the total gRNA counts in this cell. 

.. autofunction:: crispat.ga_ratio


Across cells
------------

**Poisson-Gaussian mixture model (Poisson-Gauss):** For every gRNA, this method fits a Poisson-Gaussian mixture model on the log2-transformed non-zero UMI counts of this gRNA over all cells across all batches. Next, all cells for which the probability of observing the guide counts from the Gaussian component is higher than for the Poisson (background) component are assigned to this gRNA.

.. autofunction:: crispat.ga_poisson_gauss

**Gaussian-Gaussian mixture model (Cell Ranger):** For every gRNA, this method fits a Gaussian-Gaussian mixture model on the log10-transformed UMI counts of this gRNA with a pseudocount of 1 over all cells in a batch. 

.. autofunction:: crispat.ga_cellranger


Across gRNAs and cells
----------------------
In the last group of methods, information is shared across cells and across gRNAs. Since `ga_poisson`, `ga_negative_binomial` and `ga_binomial` have the longest run time, these functions automatically are parallelized to run over all available CPUs. If you want to change this default behaviour, you can set `parallelize` to False or specify the number of processes (`n_processes`) that should be used (instead of all available CPUs).

**2-Beta mixture model (2-Beta):** Like the Ratio_X% approach, this method calculates for every cell the relative frequency of a gRNA as the ratio of its counts over the total number of gRNA counts. Using the highest ratio for every cell the method then fits a mixture model of two Beta distributions across all cells from a given batch to determine a threshold on the ratio based on where the two Beta distributions intersect. This results in one threshold per batch without distinguishing between gRNAs. This threshold is then used as X in the Ratio_X% approach. 

.. autofunction:: crispat.ga_2beta

**3-Beta mixture model (3-Beta):** This method uses the same approach as the 2-Beta model but using a 3-component Beta mixture and defining the threshold X as the intersection of the two highest components. The third component might thus capture cells infected with two gRNAs.

.. autofunction:: crispat.ga_3beta

**Latent variable Poisson generalized linear model (Poisson):** For every gRNA, this method fits a Poisson mixture on the UMI counts of this gRNA across all cells with mean :math:`\lambda = e^{\beta_0+\beta_1p_c+\beta_2b_c+log(s_c)}` with :math:`\beta_0 \in R, \beta_1 \in R^+, \beta_2 \in R^n`, perturbation state :math:`p_c \in {0,1}` and cell-specific covariates (sequencing depth :math:`s_c` and one-hot encoded batch :math:`b_c` for n batches). This approach is based on the R package SCEPTRE. 

.. autofunction:: crispat.ga_poisson

**Latent variable Negative Binomial generalized linear model (Negative Binomial):** This method is a modified version of the Poisson method using a Negative Binomial distribution instead of Poisson distribution. The overdispersion is learnt as an additional parameter in the model.

.. autofunction:: crispat.ga_negative_binomial

**Latent variable Binomial generalized linear model (Binomial):** For every gRNA, this method fits a binomial distribution :math:`B(N_c, \theta_c)` with :math:`N_c` being the total number of gRNA counts per cell and :math:`\theta_c=sigmoid(e^\beta_0+\beta_1p_c+\beta_2b_c)` with :math:`\beta_0 \in R, \beta_1 \in R^+, \beta_2 \in R^n`, perturbation state :math:`p_c \in {0,1}` and one-hot encoded batch :math:`b_c`.

.. autofunction:: crispat.ga_binomial

**Quantile approach (Top_X% cells):** For every gRNA, this method chooses the top X% of cells with the highest gRNA, excluding cells with zero counts for the gRNA.

.. autofunction:: crispat.ga_quantiles