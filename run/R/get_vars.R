library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", "gghalves",
                            "tidyr", "reshape2", "stringr", "tibble",
                            "lubridate", "hablar", "readxl",
                            "ggpubr",
                            "foreach", "doParallel", "logger",
                            "qs", "feather", "RColorBrewer"),
               format = "qs")

run_name <- "get_vars"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))

# Load exiting data

# acr_models
models_acr <- c("G_a@RKHS",  
                "G_a@RKHS&G_d@RKHS", 
                "G_a@RKHS&G_d@RKHS&G_aa@RKHS")
acr_models <- data.frame(models = paste0("M_", 1:length(models_acr)),
                         model_specs = models_acr)
# wtn models
models_wtn <- c("E_i@BRR&G_i@BRR", # M_1
                "E_i@BRR&G_a@BRR&G_d@BRR&G_aa@BRR", #M_2
                "ERM_l@BRR&G_a@BRR&G_d@BRR&G_aa@BRR", #M_3
                "ERM_l@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_l@RKHS", #M_4
                "ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR", #M_5
                "ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_nl@RKHS", #M_6
                "S@BRR&Y@BRR&G_a@BRR&G_d@BRR&G_a_S_i@BRR&G_a_Y@RKHS&G_a_ERM_l@RKHS", #M_7
                "S@BRR&Y@BRR&G_a@BRR&G_d@BRR&G_a_S@BRR&G_a_Y@RKHS&G_a_ERM_l@RKHS" #M_8
                #,"S@BRR&Y@BRR&G_a@BRR&G_d@BRR&G_a_S_p@BRR&G_a_Y@RKHS&G_a_ERM_l@RKHS" #M_9
                #,"S_al@BRR&Y@BRR&G_a@BRR&G_d@BRR&G_a_S_al@RKHS&G_a_Y@RKHS&G_a_ERM_l@RKHS" #M_10
)
wtn_models <- data.frame(models =  paste0("M_", 1:length(models_wtn)),
                         model_specs = models_wtn)

# input data
ext_parse <- "results/R"
data_paths <- list("pheno_wtn" = sprintf("/proj/%s/process_cgm_data/BLUEs_within_env_cgm.qs", ext_parse),
                   "pheno_acr" = sprintf("/proj/%s/KIBREED_data_generation/BLUES_acr_env.qs", ext_parse),
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
    name = var_acr,
    command = get_vars(existing_data = data_paths,
                       write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name), 
                       log_at = sprintf("%s/%s/acr", core_paths[["logs_R"]], run_name),
                       tmp_at = sprintf("%s/%s/acr", core_paths[["tmp_data_R"]], run_name),
                       model_info = acr_models,
                       key = "pheno_data_acr")
  ),
  tar_target(
    name = var_wtn,
    command = get_vars(existing_data = data_paths,
                       write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name), 
                       log_at = sprintf("%s/%s/wtn", core_paths[["logs_R"]], run_name),
                       tmp_at = sprintf("%s/%s/wtn", core_paths[["tmp_data_R"]], run_name),
                       model_info = wtn_models,
                       key = "pheno_data_wtn")
  ),
  tar_target(
    name = acr_overviews,
    command = make_acr_plots(linked_proj = "process_R_pred_data",
                             data = var_acr)
  )
)
