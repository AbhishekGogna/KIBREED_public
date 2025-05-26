library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data
options(java.parameters = "-Xmx100G")

tar_option_set(packages = c("dplyr", "readr", "ggplot2", "ggmap", 
                            "tidyr", "reshape2", "stringr",
                            "lubridate", "hablar", "readxl",
                            "qs", "feather", "grid", "patchwork",
                            "corehunter", "rJava", "BGLR", "cluster",
                            "foreach", "doParallel",
                            "asreml", "RColorBrewer", "dendextend",
                            "ggpubr"),
               format = "qs")

run_name <- "feature_importance"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))

# Load exiting data
ext_parse <- "results/R"
data_paths <- list("g_data" = sprintf("/proj/%s/KIBREED_data_generation/GNdata_comb_add.qs", ext_parse), # based on genetic data
                   "p_wtn" = sprintf("/proj/%s/KIBREED_data_generation/BLUES_within_env.qs", ext_parse), # raw p_data
                   #"p_wtn" = sprintf("/proj/%s/process_cgm_data/BLUEs_within_env_cgm.qs", ext_parse), # raw p_data
                   "G_a_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_a.qs", ext_parse),  # based on genetic data
                   "G_d_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_d.qs", ext_parse),  # based on genetic data
                   "G_aa_RM" = sprintf("/proj/%s/preprocessing_geno_data/kin_aa.qs", ext_parse),  # based on genetic data
                   "ERM_l" = sprintf("/proj/%s/KIBREED_data_generation/ERM_data_linear.qs", ext_parse), # based on environment data
                   "ERM_nl" = sprintf("/proj/%s/KIBREED_data_generation/ERM_data_non_linear.qs", ext_parse), # based on environment data
                   "SRM" = sprintf("/proj/%s/KIBREED_data_generation/SRM.qs", ext_parse), # based on environment data
                   "YRM" = sprintf("/proj/%s/KIBREED_data_generation/YRM.qs", ext_parse), # based on environment data
                   "G_a_S" = sprintf("/proj/%s/process_cgm_data/g_s_mat.qs", ext_parse) # based on CGM output
)


list(
  tar_target(
    name = RD_data,
    command = get_RD_mat(existing_data = data_paths,
                         write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name),
                         log_at = sprintf("%s/%s", core_paths[["logs_R"]], run_name),
                         tmp_at = sprintf("%s/%s", core_paths[["tmp_data_R"]], run_name),
                         key = run_name)
  ) 
  , tar_target(
    name = core_set,
    command = get_core_set(existing_data = data_paths,
                           data = RD_data)
  )
  , tar_target(
    name = training_data,
    command = get_training_data(existing_data = data_paths,
                                data = core_set)
  )
  , tar_target(
    name = training_data_2,
    command = get_training_data_2(existing_data = data_paths,
                                data = core_set)
  )
  , tar_target(
    name = predicted_data,
    command = predict_missing(existing_data = data_paths,
                              data_1 = training_data,
                              data_2 = training_data_2,
                              model = "ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_nl@RKHS",
                              debug = FALSE)
    
  )
  , tar_target(
    name = residual_vals,
    command = get_residuals(existing_data = data_paths,
                            data = predicted_data)
    
  )
  , tar_target(
    name = env_clusters,
    command = get_env_clusters(existing_data = data_paths,
                               data = residual_vals,
                               ec_at = "/proj/store/KIBREED_data_generation/")
  )
  , tar_target(
    name = feature_importance,
    command = get_feature_imp_plots(existing_data = data_paths,
                                    data = residual_vals,
                                    scores_at = "/proj/results/Py/feature_importance/feature_imp_scores.csv",
                                    ref_time = tar_timestamp(env_clusters))
  
  )
  , tar_target(
    name = yield_gain,
    command = get_yield_gain(existing_data = env_clusters,
                             data = predicted_data)
    
  )
)
