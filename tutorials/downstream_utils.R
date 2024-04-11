library(tidyverse)

get_power_check_results <- function(data_dir, name, method_names){
  # load results from the power check step of SCEPTRE for multiple methods
  method_name = method_names[method_names$dir_name == name, 'method']
  df = readRDS(paste0(data_dir, 'sceptre_output/',  name, 
                      '/results_run_power_check.rds')) %>%
    filter(pass_qc == TRUE) %>%
    mutate(method = method_name)
  df <- df %>% mutate(p_adj = p.adjust(df$p_value, method="BH"))
  return(df)
}

get_calibration_results <- function(data_dir, name, method_names){
  # load results from the calibration step of SCEPTRE for multiple methods
  method_name = method_names[method_names$dir_name == name, 'method']
  df = readRDS(paste0(data_dir, 'sceptre_output/', name, 
                      '/results_run_calibration_check.rds')) %>%
    filter(pass_qc == TRUE) %>%
    mutate(method = method_name)
  df <- df %>% mutate(p_adj = p.adjust(df$p_value, method="BH"))
  return(df)
}

get_discovery_analysis_results <- function(data_dir, name, method_names){
  # load results from the discovery analysis step of SCEPTRE for multiple methods
  method_name = method_names[method_names$dir_name == name, 'method']
  df <- readRDS(paste0(data_dir, 'sceptre_output/', name, 
                       '/results_run_discovery_analysis.rds')) %>%
    filter(pass_qc == TRUE, significant == TRUE) %>%
    mutate(method = method_name)
}

