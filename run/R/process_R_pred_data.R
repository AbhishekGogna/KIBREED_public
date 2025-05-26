library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

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
  #, "cv_acr_str" = tar_read(run_scripts_acr_str , 
  #                        store = sprintf("%s/store/generate_prediction_data", project_path))
  , "cv_wtn_tra" = tar_read(run_scripts_wtn_tra , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  , "cv_wtn_LoO" = tar_read(run_scripts_wtn_LoO , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  , "cv_wtn_cvL" = tar_read(run_scripts_wtn_cvL , 
                          store = sprintf("%s/store/generate_prediction_data", project_path))
  , "pred_paths_py" = jsonlite::read_json(sprintf("%s/create_slurm_scripts/pred_dirs_paths.json", core_paths[["results_Py"]]))
  , "data_bahareh_at" = sprintf("%s/ext_dir/KIBREED/source_data/ZALF_Bonn/21_11_23_data", project_path)
)

# define pipeline
list(
  tar_target(
    name = pred_file_paths,
    command = get_pred_file_paths(existing_data = data, 
                                  log_at = sprintf("%s/%s", core_paths[["logs_R"]], run_name),
                                  tmp_at = sprintf("%s/%s", core_paths[["tmp_data_R"]], run_name))
  )
  , tar_target(
    name = pred_data,
    command = load_pred_data(data = pred_file_paths,
                             write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name),
                             geno_mapping_from = sprintf("%s/results/R/process_cgm_data/BLUEs_within_env_cgm.qs", project_path))
  )
  , tar_target(
    name = overviews,
    command = make_plots(data = pred_data,
                         write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name))
                             
  )
  #, tar_target(
  #  name = reset_job_meta_data,
  #  command = change_job_meta_data(data = pred_file_paths,
  #                                 new_time = "2-0:0:0",
  #                                 cv_type_vec = c("cv_wtn_tra"),
  #                                 model_name_vec = paste0("M_", c()))
  #)
)
