Downstream Analyses 
--------------------------

To investigate the consequences of the assignment differences for discovery analysis, the assignments from crispat can serve as input to differential expression tests performed by various different R and python packages. In this directory we provide scripts that show how the assignments obtained by crispat can be used for differential expression analysis with

 - **scanpy** (python)
 - **Seurat** (R)
 - **SCEPTRE** (R)

We recommend to use the R package **SCEPTRE** for downstream analyses. It is tailored to single-cell CRISPR screen analysis and combines multiple analyses and control checks in a statistical rigorous fashion. It consists of the following steps:

1) power check: calculates the log2 fold changes and p-values for the target genes
2) calibration check: calculates the number of false discoveries based on genes that are significantly differentially expressed in cells assigned to one non-targeting gRNA against all other non-targeting control cells
3) discovery analysis: calculates differentially expressed genes for each perturbation

Details on SCEPTRE and an explanation of all options can be found [here](https://timothy-barry.github.io/sceptre-book/sceptre.html). 

All code used for the analyses shown in our paper such as the comparison of the differential expression analysis output obtained by SCEPTRE across methods can be found in a separate [github repository](https://github.com/velten-group/crispat_analysis). 

## Practical considerations for method comparison
For small screens the whole downstream analyses can be run for all guide assignment methods in a reasonable time. However, we note  that differential expression analysis with SCEPTRE can be time consuming for large screens such as the screen by Replogle et al. (around 2200 gRNAs and 600k cells) used in our paper. For such screens, we recommend to first perform the differential expression analysis on a subset of the data only (e.g. subsetting to 40 random targeting gRNAs and all non-targeting gRNAs as done in our analysis). 

Based on the results of steps 1) and 2) users of crispat can select a suitable crispat assignment method based on true and false positive rate, calculated as the number of target genes with significantly differential expression compared to the control cells over the total number of all targeted genes (TPR) and the number of significantly differentially expressed genes in the calibration check out of all tested non-targeting gRNA - gene pairs (FPR).
