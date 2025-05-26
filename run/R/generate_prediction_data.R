library(targets)

project_path <- "/proj"
core_paths <- lapply(jsonlite::read_json(sprintf("%s/results/core_paths.json", project_path)), function(x) x[[1]])

# Define packages needed and default storage format for intermittent data

tar_option_set(packages = c("dplyr", "readr", "ggplot2", 
                            "tidyr", "reshape2", "stringr",
                            "lubridate", "hablar", "readxl",
                            "AGHmatrix",
                            "qs", "feather", "jsonlite", "cvTools"),
               format = "qs")

run_name <- "generate_prediction_data"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))

# Load exiting data

data <- list(#"pheno_data_acr" = tar_read(KIBREED_data_full,
             #                            store = sprintf("%s/store/KIBREED_data_generation", project_path))[["BLUES_acr_env"]],  
             "geno_data_a" = tar_read(KIBREED_data_full,
                                      store = sprintf("%s/store/KIBREED_data_generation", project_path))[["GNdata_comb_add"]],
             "geno_data_d" = tar_read(KIBREED_data_full,
                                      store = sprintf("%s/store/KIBREED_data_generation", project_path))[["GNdata_comb_dom"]],
             "ec_data" = tar_read(KIBREED_data_full,
                                  store = sprintf("%s/store/KIBREED_data_generation", project_path))[["ERM_data"]],
             "climate_data" = tar_read(KIBREED_data_full,
                                       store = sprintf("%s/store/KIBREED_data_generation", project_path))[["env_data_kibreed_raw"]],
             #"pheno_data_wtn" = tar_read(BLUES_within_env_cgm ,
             #                            store = sprintf("%s/store/process_cgm_data", project_path))[["BLUEs_within_env_cgm"]],
             "param_data" = tar_read(BLUES_within_env_cgm ,
                                     store = sprintf("%s/store/process_cgm_data", project_path))[["params_data"]]
)

# acr_models
models_acr <- c("G_a@BRR",  
                "G_a@BRR&G_d@BRR", 
                "G_a@BRR&G_d@BRR&G_aa@BRR")
acr_models <- data.frame(models = paste0("M_", 1:length(models_acr)),
                         model_specs = models_acr,
                         cpu = sapply(models_acr, function(x) length(strsplit(x, "&")[[1]]), USE.NAMES = F),
                         memory = c(50, 50, 70),
                         time = "1-0:0:0")

#"ERM_l@BRR&ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_l@RKHS&G_a_ERM_nl@RKHS&G_d_ERM_l@RKHS&G_d_ERM_nl@RKHS",#omit
#"ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_nl@RKHS&G_d_ERM_nl@RKHS", #omit
#"ERM_l@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_l@RKHS&G_d_ERM_l@RKHS", #omit
#"E_i@BRR&G_a@BRR&G_d@BRR", #omit

# wtn models
models_wtn <- c("E_i@BRR&G_i@BRR", # M_1
                "E_i@BRR&G_a@BRR&G_d@BRR&G_aa@BRR", #M_2
                "ERM_l@BRR&G_a@BRR&G_d@BRR&G_aa@BRR", #M_3
                "ERM_l@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_l@RKHS", #M_4
                "ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR", #M_5
                "ERM_nl@BRR&G_a@BRR&G_d@BRR&G_aa@BRR&G_a_ERM_nl@RKHS", #M_6
                "S@BRR&Y@BRR&G_a@BRR&G_d@BRR&G_a_S_i@BRR&G_a_Y@RKHS&G_a_ERM_l@RKHS", #M_7
                "S@BRR&Y@BRR&G_a@BRR&G_d@BRR&G_a_S@BRR&G_a_Y@RKHS&G_a_ERM_l@RKHS", #M_8
                "S@BRR&Y@BRR&G_a@BRR&G_d@BRR&G_a_S_p@BRR&G_a_Y@RKHS&G_a_ERM_l@RKHS", #M_9
                "S_al@BRR&Y@BRR&G_a@BRR&G_d@BRR&G_a_S_al@RKHS&G_a_Y@RKHS&G_a_ERM_l@RKHS" #M_10
                )

#to_run <- data.frame(model = models_wtn[1:6], cv = rep(paste0("cv", 1:4, "_tra"), each = 6))
#to_run_2 <- data.frame(model = models_wtn[c(1:7, 8)], cv = rep(paste0("cv1_tra"), times = 2))
#to_run_3 <- data.frame(model = models_wtn[c(1:7, 9)], cv = rep(paste0("cv3_tra"), times = 2))
#to_run_5 <- data.frame(model = models_wtn[c(1:7)], cv = rep(paste0("cv2_Lo0"), times = 7)) #cv2
#to_run_4 <- data.frame(model = models_wtn[c(1:7, 10)], cv = rep(paste0("cv4_cvL"), times = 4)) #cv4
#
#to_run<- bind_rows(to_run, to_run_2, to_run_3, to_run_4, to_run_5) %>% distinct()

wtn_models <- data.frame(models =  paste0("M_", 1:length(models_wtn)),
                         model_specs = models_wtn,
                         cpu = sapply(models_wtn, function(x) length(strsplit(x, "&")[[1]]), USE.NAMES = F),
                         mem = c(30, 50, 100, 200, 100, 200, 200, 200, 200,200),
                         time = "2-0:0:0")

# todo: rewrite the log file if the process is rerun from the top
# define pipeline
list(
  tar_target(
    name = pheno_data_acr_path,
    sprintf("%s/results/R/KIBREED_data_generation/BLUES_acr_env.qs", project_path), 
    format = "file"
  ),
  tar_target(
    name = data_acr,
    command = get_data(pheno_data_path = pheno_data_acr_path,
                       existing_data = data,
                       key = "pheno_data_acr")
  ),
  tar_target(
    name = pred_acr_objects,
    command = create_pred_acr_objects(existing_data = data_acr,
                                      write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name), 
                                      log_at = sprintf("%s/%s/acr", core_paths[["logs_R"]], run_name),
                                      tmp_at = sprintf("%s/%s/acr", core_paths[["tmp_data_R"]], run_name),
                                      run_name = run_name)
  ),
  tar_target(
    name = cv_acr_5f_data,
    command = cv_acr_5f(data = pred_acr_objects, runs = 10, folds = 5) # should also provide eigen info
  ),
  tar_target(
    name = run_scripts_acr_5f,
    command = generate_run_scripts(data = cv_acr_5f_data,
                                   run_script_at = c(sprintf("%s/scr_genomic_prediction_acr.R", core_paths[["src_R"]])),
                                   input_data_at = c(sprintf("%s", core_paths[["results_R"]])),
                                   model_info = acr_models)
  ),
  #tar_target(
  #  name = cv_acr_str_data,
  #  command = cv_acr_str(data = pred_acr_objects, runs = 50, take_parts = TRUE)
  #),
  #tar_target(
  #  name = run_scripts_acr_str,
  #  command = generate_run_scripts(data = cv_acr_str_data,
  #                                 run_script_at = c(sprintf("%s/scr_genomic_prediction_acr.R", core_paths[["src_R"]])),
  #                                 input_data_at = c(sprintf("%s", core_paths[["results_R"]])),
  #                                 model_info = acr_models)
  #),
  tar_target(
    name = cv_acr_sce_data,
    command = cv_acr_sce(data = pred_acr_objects, take_parts = FALSE)
  ),
  tar_target(
    name = run_scripts_acr_sce,
    command = generate_run_scripts(data = cv_acr_sce_data,
                                   run_script_at = c(sprintf("%s/scr_genomic_prediction_acr.R", core_paths[["src_R"]])),
                                   input_data_at = c(sprintf("%s", core_paths[["results_R"]])),
                                   model_info = acr_models)
  ),
  tar_target(
    name = pheno_data_wtn_path,
    sprintf("%s/results/R/process_cgm_data/BLUEs_within_env_cgm.qs", project_path), 
    format = "file"
  ),
  tar_target(
    name = data_wtn,
    command = get_data(pheno_data_path = pheno_data_wtn_path,
                       #existing_data = data,
                       key = "pheno_data_wtn")
  ),
  tar_target(
    name = pred_wtn_objects,
    command = create_pred_wtn_objects(existing_data = data_wtn,
                                      write_at = sprintf("%s/%s", core_paths[["results_R"]], run_name),
                                      log_at = sprintf("%s/%s/wtn", core_paths[["logs_R"]], run_name),
                                      tmp_at = sprintf("%s/%s/wtn", core_paths[["tmp_data_R"]], run_name))
  ),
  tar_target(
    name = cv_wtn_tra_data,
    command = cv_wtn_tra(data = pred_wtn_objects, runs = 50, test_prop = 0.33)
  ),
  tar_target(
    name = run_scripts_wtn_tra,
    command = generate_run_scripts(data = cv_wtn_tra_data,
                                   run_script_at = c(sprintf("%s/scr_genomic_prediction_wtn.R", core_paths[["src_R"]])),
                                   input_data_at = c(sprintf("%s", core_paths[["results_R"]])),
                                   model_info = wtn_models[1:9, ], # here cv1 and cv3 will run with common M1 to M7 but different M8 and M9
                                   wtn = TRUE,
                                   reservation = NULL)
  ),
  tar_target(
    name = cv_wtn_LoO_data,
    command = cv_wtn_LoO(data = pred_wtn_objects)
  ),
  tar_target(
    name = run_scripts_wtn_LoO,
    command = generate_run_scripts(data = cv_wtn_LoO_data,
                                   run_script_at = c(sprintf("%s/scr_genomic_prediction_wtn.R", core_paths[["src_R"]])),
                                   input_data_at = c(sprintf("%s", core_paths[["results_R"]])),
                                   model_info = wtn_models[c(1:7), ], # cv2 alt will get inputs from bahareh
                                   wtn = TRUE,
                                   reservation = NULL)
  ),
  tar_target(
    name = cv_wtn_cvL_data,
    command = cv_wtn_cvL(data = pred_wtn_objects)
  ),
  tar_target(
    name = run_scripts_wtn_cvL,
    command = generate_run_scripts(data = cv_wtn_cvL_data,
                                   run_script_at = c(sprintf("%s/scr_genomic_prediction_wtn.R", core_paths[["src_R"]])),
                                   input_data_at = c(sprintf("%s", core_paths[["results_R"]])),
                                   model_info = wtn_models[c(1:7, 10), ], # cv4 alt with parameter derived location kinship 
                                   wtn = TRUE,
                                   reservation = NULL)
  )
  #,
  #tar_target(
  #  name = cv_wtn_trsz_data, #tr-training, sz-size
  #  command = cv_wtn_trsz(data = pred_wtn_objects, runs = 20, test_prop = 0.33)
  #),
  #tar_target(
  #  name = run_scripts_wtn_trsz,
  #  command = generate_run_scripts(data = cv_wtn_trsz_data,
  #                                 run_script_at = c(sprintf("%s/scr_genomic_prediction_wtn.R", core_paths[["src_R"]])),
  #                                 input_data_at = c(sprintf("%s", core_paths[["results_R"]])),
  #                                 model_info = wtn_models[c(1, 2), ],
  #                                 wtn = TRUE,
  #                                 reservation = NULL)
  #)
)
#paths <- unique(sapply(run_scripts_wtn_5f$run_data$run_data, function(x) x[["run_scripts"]]))
#path_files <- unique((sapply(paths, function(x) list.files(x, full.names = T))))
#for (i in path_files) {
#  file <- readLines(i)
#  file[grepl('ext_lib_blas="/qg-10/', file)] <- "ext_lib_blas=\"/qg-10/data/AGR-QG/Gogna/computing_containers/openblas_3.23/inst/qg-10.ipk-gatersleben.de/lib/libopenblas.so\""
#  cat(file, file = i, sep = "\n")
#  system(sprintf("chmod +x %s", i))
#}
