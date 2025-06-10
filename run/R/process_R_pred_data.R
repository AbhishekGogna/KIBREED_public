library(targets)

# Set basic directory structure
if (dir.exists("/proj")) {
  project_path <- "/proj"
} else {
  project_path <- getwd() # assumes that this script is run from KIBREED_public
}

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", "ggpubr", "gghalves", "ggh4x",
                            "tidyr", "reshape2", "stringr",
                            "lubridate", "patchwork", "hablar", "readxl",
                            "foreach", "doParallel",
                            "Metrics", "BBmisc",
                            "qs", "feather"),
               format = "qs")

run_name <- "process_R_pred_data"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))

# Load exiting data

data <- list(
  "cv_acr_5f" = tar_read(run_scripts_acr_5f , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  , "cv_acr_sce" = tar_read(run_scripts_acr_sce , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  , "cv_wtn_tra" = tar_read(run_scripts_wtn_tra , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
)

# define pipeline
list(
  tar_target(
    name = pred_file_paths,
    command = get_pred_file_paths(existing_data = data, 
                                  log_at = sprintf("%s/logs/R/%s", project_path, run_name),
                                  tmp_at = sprintf("%s/tmp_data/R/%s", project_path, run_name))
  )
  #, tar_target(
  #  name = pred_data,
  #  command = load_pred_data(data = pred_file_paths,
  #                           write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name),
  #                           geno_mapping_from = sprintf("%s/results/R/process_cgm_data/BLUEs_within_env_cgm.qs", project_path))
  #)
  #, tar_target(
  #  name = overviews,
  #  command = make_plots(data = pred_data,
  #                       write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name))
  #                           
  #)
)
