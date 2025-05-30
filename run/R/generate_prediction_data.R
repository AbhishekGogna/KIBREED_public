library(targets)

# Set basic directory structure
if (dir.exists("/proj")) {
  project_path <- "/proj"
} else {
  project_path <- getwd() # assumes that this script is run from KIBREED_public
}

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "ggplot2", "tidyr", "stringr",
                            "hablar", "qs", "jsonlite", "cvTools",
                            "tibble", "BGLR"),
               format = "qs")

run_name <- "generate_prediction_data"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))

# acr_models
models_acr <- c("G_a@BRR",  
                "G_a@BRR&G_d@BRR", 
                "G_a@BRR&G_d@BRR&G_aa@BRR")
acr_models <- data.frame(models = paste0("M_", 1:length(models_acr)),
                         model_specs = models_acr)

# wtn models
models_wtn <- c("E_i@BRR&G_i@BRR", # M_1
                "E_i@BRR&G_a@BRR&G_d@BRR&G_aa@BRR", #M_2
                "ERM_l@BRR&G_a@BRR&G_d@BRR&G_aa@BRR", #M_3
                "ERM_l@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_l@RKHS", #M_4
                "ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR", #M_5
                "ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_nl@RKHS" #M_6
                )

wtn_models <- data.frame(models =  paste0("M_", 1:length(models_wtn)),
                         model_specs = models_wtn)

# define pipeline #to run fresh try deleting sub directories in results, tmp_data and logs
list(
  tar_target(
    name = pred_acr_objects,
    command = create_pred_acr_objects(existing_data_path = sprintf("%s/%s", project_path, "data"),
                                      write_at = sprintf("%s/results/R/%s", project_path, run_name),
                                      log_at = sprintf("%s/logs/R/%s", project_path, run_name),
                                      tmp_at = sprintf("%s/tmp_data/R/%s", project_path, run_name))
  ),
  tar_target(
    name = cv_acr_5f_data,
    command = cv_acr_5f(data = pred_acr_objects, runs = 10, folds = 5) # should also provide eigen info
  ),
  tar_target(
    name = run_scripts_acr_5f,
    command = generate_run_scripts(data = cv_acr_5f_data,
                                   run_script_at = sprintf("%s/src/R/scr_genomic_prediction_acr.R", project_path),
                                   input_data_at = sprintf("%s/results/R", project_path),
                                   model_info = acr_models,
                                   project_path = project_path)
  ),
  tar_target(
    name = cv_acr_sce_data,
    command = cv_acr_sce(data = pred_acr_objects, take_parts = FALSE)
  ),
  tar_target(
    name = run_scripts_acr_sce,
    command = generate_run_scripts(data = cv_acr_sce_data,
                                   run_script_at = sprintf("%s/src/R/scr_genomic_prediction_acr.R", project_path),
                                   input_data_at = sprintf("%s/results/R", project_path),
                                   model_info = acr_models,
                                   project_path = project_path)
  ),
  tar_target(
    name = pred_wtn_objects,
    command = create_pred_wtn_objects(existing_data_path = sprintf("%s/%s", project_path, "data"),
                                      write_at = sprintf("%s/results/R/%s", project_path, run_name),
                                      log_at = sprintf("%s/logs/R/%s", project_path, run_name),
                                      tmp_at = sprintf("%s/tmp_data/R/%s", project_path, run_name))
  ),
  tar_target(
    name = cv_wtn_tra_data,
    command = cv_wtn_tra(data = pred_wtn_objects, runs = 50, test_prop = 0.33)
  ),
  tar_target(
    name = run_scripts_wtn_tra,
    command = generate_run_scripts(data = cv_wtn_tra_data,
                                   run_script_at = sprintf("%s/src/R/scr_genomic_prediction_wtn.R", project_path),
                                   input_data_at = sprintf("%s/results/R", project_path),
                                   model_info = wtn_models,
                                   wtn = TRUE,
                                   project_path = project_path)
  )
)
