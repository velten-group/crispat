library(tidyverse)
library(sceptre)
library(Seurat)
library(Matrix)
source('utils_SCEPTRE.R')

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if there are any arguments
if (length(args) > 0) {
  data_dir <- args[1]
  data_name <- args[2]
  batches <- args[3]
  batches <- as.integer(strsplit(batches, ",")[[1]])
  cat("Batches:", batches, "\n")
  annotation_file <- args[4] # path to the annotation file specifying the target of each gRNA
  gc_file <- args[5] # path to file containing the guide assignment results
  cat("Guide calling file:", gc_file, "\n")
  save_dir <- args[6]
} else {
  cat("Please provide the following 6 arguments: data directory, data name, batch numbers, annoatation file, guide calling result file and output directory \n")
}

# Load annotation data specifying the target gene name of each gRNA
# csv with column names 'gRNA' and 'target_gene'
print('Load data')
t0 <- Sys.time()
gRNA_target_df <- read_csv(annotation_file, show_col_types = FALSE) %>%
  rename(grna_id = gRNA, grna_target = target_gene) 

# Read in guide assignment results
PGMM <- read_csv(gc_file, show_col_types = FALSE) 

## subset to cells with single perturbations
single_perturbations <- PGMM %>% subset(select = c(cell, gRNA)) %>%
  distinct() %>%
  group_by(cell) %>%
  mutate(n = n()) %>%
  filter(n == 1)
PGMM <- PGMM %>%
  filter(cell %in% single_perturbations$cell, gRNA %in% gRNA_target_df$grna_id)

assigned_cells <- unique(PGMM$cell)

# Create sceptre object for assigned cells
out <- create_sceptre(gRNA_target_df, assigned_cells, batches, "low", data_dir, data_name)
sceptre_object <- out$sceptre
cell_indices <- out$cell_indices
t1 <- Sys.time()
print(round(t1 - t0,2))

# Set analysis parameters
print('Set analysis parameters')
positive_control_pairs <- construct_positive_control_pairs(
  sceptre_object = sceptre_object
)
discovery_pairs <- construct_trans_pairs(sceptre_object = sceptre_object)

# set data set specific parameters
if (data_name == 'Replogle_K562' | data_name == 'Replogle_RPE1'){
  selected_genes <- read_csv(paste0(data_dir, 'selected_response_genes.csv'), show_col_types = FALSE)
  discovery_pairs <- filter(discovery_pairs, response_id %in% selected_genes$genes)
  formula_object = formula(~ log(response_n_nonzero) + log(response_n_umis) +
                             log(grna_n_nonzero + 1) + log(grna_n_umis + 1) +
                             response_p_mito + Batch + S_score + G2M_score)
  strategy = "union"
} else if (data_name == 'Schraivogel_whole'){
  selected_genes <- read_csv(paste0(data_dir, 'selected_response_genes.csv'), show_col_types = FALSE)
  discovery_pairs <- filter(discovery_pairs, response_id %in% selected_genes$genes)
  strategy = "singleton"
  formula_object = formula(~ log(response_n_nonzero) + log(response_n_umis) +
                             log(grna_n_nonzero + 1) + log(grna_n_umis + 1) +
                             response_p_mito + S_score + G2M_score)
} else{
  strategy = "singleton"
  formula_object = formula(~ log(response_n_nonzero) + log(response_n_umis) +
                             log(grna_n_nonzero + 1) + log(grna_n_umis + 1))
}


sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = positive_control_pairs,
  side = "both", 
  grna_integration_strategy = strategy, 
  control_group = 'nt_cells',
  formula_object = formula_object,
  multiple_testing_alpha = 0.05
)
t2 <- Sys.time()
print(round(t2 - t1,2))

# Add guide assignment
print('Add own guide assignment to SCEPTRE')
## add cell indices
PGMM <- PGMM %>%
  inner_join(cell_indices) %>%
  mutate(cell_index = as.integer(cell_index))

# convert it to a list
PGMM_list <- split(PGMM$cell_index, PGMM$gRNA)
# add gRNAs without any assigned cells
null_grnas <- filter(gRNA_target_df, ! grna_id %in% PGMM$gRNA)$grna_id
PGMM_list <- c(PGMM_list, 
               setNames(lapply(rep(list(integer(0)), length(null_grnas)), function(x) x), null_grnas))

# get processed assignment outputs
processed_assignment_out <- process_initial_assignment_list(initial_assignment_list = PGMM_list,
                                                            grna_target_data_frame = sceptre_object@grna_target_data_frame,
                                                            n_cells = ncol(sceptre_object@grna_matrix[[1]]), 
                                                            low_moi = sceptre_object@low_moi,
                                                            maximum_assignment = FALSE, 
                                                            perturbation_df = PGMM)

# write guide assignment results to the sceptre object
sceptre_object@grna_assignments_raw <- processed_assignment_out$grna_assignments_raw
sceptre_object@grnas_per_cell <- as.integer(processed_assignment_out$grnas_per_cell)
sceptre_object@cells_w_multiple_grnas <- processed_assignment_out$cells_w_multiple_grnas
sceptre_object@initial_grna_assignment_list <- PGMM_list 
sceptre_object@last_function_called <- "assign_grnas"
sceptre_object@functs_called["assign_grnas"] <- TRUE
sceptre_object@grna_assignment_method = "PGMM"
sceptre_object@cells_in_use <- sort(unique(PGMM$cell_index))
t3 <- Sys.time()
print(round(t3 - t2,2))

# Run quality control
print('Run quality control')

## set quality thresholds depending on the data set 
if (data_name == 'Replogle_K562' | data_name == 'Replogle_RPE1'){
  n_umis_range = c(0,1)
  n_nonzero_range = c(0,1)
} else{
  n_umis_range = c(0.01,0.99)
  n_nonzero_range = c(0.01,0.99)
}

sceptre_object <- run_qc(sceptre_object,
                         response_n_umis_range = n_umis_range, 
                         response_n_nonzero_range = n_nonzero_range, 
                         p_mito_threshold = 0.2,
                         n_nonzero_trt_thresh = 0L,
                         n_nonzero_cntrl_thresh = 7L)
t4 <- Sys.time()
print(round(t4 - t3,2))

# Run calibration check
print('Run calibration check')
sceptre_object <- run_calibration_check(
  sceptre_object = sceptre_object, print_progress = FALSE
)
t5 <- Sys.time()
print(round(t5 - t4,2))

# Run power check
print('Run power check')
sceptre_object <- run_power_check(
  sceptre_object = sceptre_object, print_progress = FALSE
)
t6 <- Sys.time()
print(round(t6 - t5,2))

# Run discovery analysis
print('Run discovery analysis')
sceptre_object <- run_discovery_analysis(
  sceptre_object = sceptre_object, print_progress = FALSE
)
t7 <- Sys.time()
print(round(t7 - t6,2))

# Save the results
write_outputs_to_directory(
  sceptre_object = sceptre_object, 
  directory = save_dir
)

times <- data.frame('step' = c('create SCEPTRE', 'set analysis parameters', 'read in guide assignment',
                               'quality control', 'calibration check', 'power check', 'discovery analysis'), 
                    'time' = c(t1-t0, t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6))
write_csv(times, paste0(save_dir, '/run_times.csv'))

print(paste0('Done! Outputs are saved in ', save_dir))
