# This script contains functions used in the SCEPTRE analyses

get_cov_per_batch <- function(data, batch){
  # create seurat object
  seurat_object = CreateSeuratObject(counts = data$`Gene Expression`)
  seurat_object[["CRISPR"]] = CreateAssayObject(counts = data$`CRISPR Guide Capture`)
  
  # calculate percent mitochondrial gene expression
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  
  # preprocessing
  seurat_object_proc <- NormalizeData(seurat_object)
  seurat_object_proc <- FindVariableFeatures(seurat_object_proc, selection.method = "vst")
  seurat_object_proc <- ScaleData(seurat_object_proc, features = rownames(seurat_object_proc))
  
  # get intersection of cell cycle genes and measured gene expression
  cell_cycle <- TRUE
  measured_genes <- rownames(seurat_object_proc@assays$RNA)
  if (length(intersect(measured_genes, cc.genes$s.genes)) < 5 | length(intersect(measured_genes, cc.genes$g2m.genes)) < 5){
    cell_cycle <- FALSE
  }
  
  # calculate cell cycle scores if the expression of enough genes was measured
  if (cell_cycle == TRUE) {
    seurat_object_proc <- CellCycleScoring(seurat_object_proc, s.features = cc.genes$s.genes, 
                                           g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
    covariates <- data.frame('cell' = colnames(seurat_object_proc), 
                             'S_score' = seurat_object_proc$S.Score, 
                             'G2M_score' = seurat_object_proc$G2M.Score, 
                             'Batch' = batch, 
                             'total_UMI_counts' = seurat_object_proc@meta.data$nCount_RNA,
                             'percent_mito' = seurat_object_proc@meta.data$percent.mt)
  } else {
    covariates <- data.frame('cell' = colnames(seurat_object_proc), 
                             'Batch' = batch, 
                             'total_UMI_counts' = seurat_object_proc@meta.data$nCount_RNA,
                             'percent_mito' = seurat_object_proc@meta.data$percent.mt)
  }
  
  return(covariates)
}

get_data_per_batch_replogle <- function(gRNA_target_df, assigned_cells, batch_number, data_dir, data_name){
  # load data
  batch_data <- Read10X(data.dir = paste0(data_dir, 'batch', batch_number),
                        strip.suffix = TRUE)
  
  # sum gRNA counts over the two gRNAs per pair
  guide_matrix <- as.matrix(batch_data$`CRISPR Guide Capture`)
  summed_df <- guide_matrix %>% 
    data.frame() %>% 
    rownames_to_column(var = 'individual_gRNA') %>% 
    inner_join(subset(gRNA_target_df, select = c(individual_gRNA, grna_id))) %>% 
    group_by(grna_id) %>%
    summarize(across(all_of(colnames(guide_matrix)), sum)) %>%
    column_to_rownames(var = 'grna_id') 
  
  # remove cells that are not assigned (to save memory)
  assignment_df = data.frame('cell' = sub("-.*", "", assigned_cells),
                             'batch' = sub(".*-", "", assigned_cells))
  assigned_cells_batch <- filter(assignment_df, batch == batch_number)$cell
  summed_df <- summed_df[, assigned_cells_batch]
  batch_data$`CRISPR Guide Capture` <- Matrix(as.matrix(summed_df), sparse = TRUE)
  batch_data$`Gene Expression` <- batch_data$`Gene Expression`[, assigned_cells_batch]
  
  # add batch to the cell names
  colnames(batch_data$`Gene Expression`) <- paste0(colnames(batch_data$`Gene Expression`), '-', batch_number)
  colnames(batch_data$`CRISPR Guide Capture`) <- paste0(colnames(batch_data$`CRISPR Guide Capture`), '-', batch_number)
  
  # calculate cell cycle score covariates, proportion mitochondrial genes and total UMIs
  batch_covariates <- get_cov_per_batch(batch_data, as.character(batch_number))
  
  # remove cells with less than 3000 counts and more than 20%/11% mitochondrial genes
  if (data_name == 'Replogle_RPE1'){
    batch_covariates <- filter(batch_covariates, total_UMI_counts > 3000, percent_mito < 11)
  } else{
    batch_covariates <- filter(batch_covariates, total_UMI_counts > 3000, percent_mito < 20)
  }
  batch_data$`Gene Expression` <- batch_data$`Gene Expression`[, batch_covariates$cell]
  batch_data$`CRISPR Guide Capture` <- batch_data$`CRISPR Guide Capture`[, batch_covariates$cell]
  
  return(list(grna_matrix = batch_data$`CRISPR Guide Capture`,
              response_matrix = batch_data$`Gene Expression`,
              covariates = batch_covariates))
}


get_data_csv <- function(gRNA_target_df, assigned_cells, batch_number, data_dir){
  # load data
  ge_df = read_csv(paste0(data_dir, 'gene_expression_counts.csv')) %>%
    column_to_rownames(var = 'gene')
  gRNA_df = read_csv(paste0(data_dir, 'gRNA_counts.csv')) %>%
    column_to_rownames(var = 'gRNA')
  batch_data = list('Gene Expression' = ge_df, 
                    'CRISPR Guide Capture' = gRNA_df)
  
  # calculate cell cycle score covariates, proportion mitochondrial genes and total UMIs
  batch_covariates <- get_cov_per_batch(batch_data, "1")
  
  return(list(grna_matrix = data.matrix(batch_data$`CRISPR Guide Capture`),
              response_matrix = data.matrix(batch_data$`Gene Expression`),
              covariates = batch_covariates))
}


create_sceptre <- function(gRNA_target_df, assigned_cells, batches, moi, data_dir, name){
  # load data for the first batch
  if (name == 'Replogle_K562' | name == 'Replogle_RPE1'){
    batch_data <- get_data_per_batch_replogle(gRNA_target_df, assigned_cells, batches[1], data_dir, name) 
  } else {
    batch_data <- get_data_csv(gRNA_target_df, assigned_cells, batches[1], data_dir) 
  }
  
  grna_matrix <- batch_data$grna_matrix
  response_matrix <- batch_data$response_matrix
  all_covariates <- batch_data$covariates
  
  if (name == 'Replogle_K562' | name == 'Replogle_RPE1'){
    for (batch in batches[-1]){
      # get data per batch
      batch_data <- get_data_per_batch_replogle(gRNA_target_df, assigned_cells, batch, data_dir, name) 
      
      # merge with the other batches
      grna_matrix = cbind(grna_matrix, batch_data$grna_matrix)
      response_matrix = cbind(response_matrix, batch_data$response_matrix)
      all_covariates <- rbind(all_covariates, batch_data$covariates)
    }
  }
  
  # save the cell indices
  cell_indices = data.frame(cell_index = seq(1, ncol(grna_matrix)),
                            cell = colnames(grna_matrix))
  
  gRNA_target_df <- filter(gRNA_target_df, grna_id %in% rownames(grna_matrix)) %>% 
    subset(select = c(grna_id, grna_target)) %>%
    distinct()

  # create SCEPTRE object
  sceptre_object = import_data(response_matrix = response_matrix,
                               grna_matrix = grna_matrix,
                               grna_target_data_frame = gRNA_target_df,
                               moi = moi,
                               extra_covariates = all_covariates)
  return(list(sceptre = sceptre_object, cell_indices = cell_indices))
}

# the following function is a slightly modified version of the SCETPRE function 
# (modified such that it doesn't have to use a C++ function but computes the grnas_per_cell based on the perturbation_df)
process_initial_assignment_list <- function(initial_assignment_list, grna_target_data_frame, n_cells, low_moi, maximum_assignment, perturbation_df) {
  # 1. compute the vector of grnas per cell
  cell_indices <- data.frame('cell_index' = seq(n_cells))
  perturbations <- perturbation_df %>% 
    subset(select = c('cell_index', 'gRNA')) %>% 
    distinct() %>% 
    group_by(cell_index) %>% 
    summarise('n_perturbations' = n()) %>%
    right_join(cell_indices) %>% 
    mutate(n_perturbations = ifelse(is.na(n_perturbations), 0, n_perturbations)) %>%
    arrange(cell_index)
  grnas_per_cell <- perturbations$n_perturbations
  
  # 2. determine the cells that contain multiple grnas (if in low MOI and not using maximum_assignment)
  if (low_moi && !maximum_assignment) {
    cells_w_multiple_grnas <- which(grnas_per_cell >= 2L)
  } else {
    cells_w_multiple_grnas <- integer()
  }
  # 3. pool together targeting gRNAs via the or operation
  targeting_grna_group_data_table <- grna_target_data_frame |>
    dplyr::filter(grna_group != "non-targeting") |> data.table::as.data.table()
  targeting_grna_groups <- targeting_grna_group_data_table$grna_group |> unique()
  grna_group_idxs <- lapply(targeting_grna_groups, function(targeting_grna_group) {
    curr_grna_ids <- targeting_grna_group_data_table[
      targeting_grna_group_data_table$grna_group == targeting_grna_group,]$grna_id
    initial_assignment_list[curr_grna_ids] |> unlist() |> unique()
  }) |> stats::setNames(targeting_grna_groups)
  # 4. obtain the individual non-targeting grna idxs
  nontargeting_grna_ids <- grna_target_data_frame |>
    dplyr::filter(grna_group == "non-targeting") |> dplyr::pull(grna_id)
  indiv_nt_grna_idxs <- initial_assignment_list[nontargeting_grna_ids]
  # 5. construct the grna_group_idxs list
  grna_assignments_raw <- list(grna_group_idxs = grna_group_idxs,
                               indiv_nt_grna_idxs = indiv_nt_grna_idxs)
  # 6. compute the number of targeting grna groups per cell
  out <- list(grna_assignments_raw = grna_assignments_raw,
              grnas_per_cell = grnas_per_cell,
              cells_w_multiple_grnas = cells_w_multiple_grnas)
  return(out)
}
