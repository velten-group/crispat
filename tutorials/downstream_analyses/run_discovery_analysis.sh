#!/bin/bash
#SBATCH --partition=single
#SBATCH --ntasks=48
#SBATCH --time=80:00:00
#SBATCH --mem=50gb
#SBATCH --export=NONE
export OMP_NUM_THREADS=${SLURM_NTASKS}

module load math/R/4.2.1

data_dir="../data/Replogle/Replogle_K562_essential/"
data="Replogle_K562"
batches="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48"
annotation_file="../data/Replogle/gRNA_pair_selection_50_sceptre.csv"
gc_file="../data/guide_calling/Replogle_K562/2-BetaMM/perturbations.csv"
save_dir="../data/sceptre_pipeline/discovery_analysis/replogle_K562/2-BetaMM"

echo '2-Beta Mixture Model'
Rscript SCEPTRE_discovery.R "$data_dir" "$data" "$batches" "$annotation_file" "$gc_file" "$save_dir" 
