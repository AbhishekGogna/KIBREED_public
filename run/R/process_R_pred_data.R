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

# For t-tests ------------------------------------------------------------
#to_check <- "plot_wtn_with_cgm"
##to_check <- "plot_wtn_no_cgm"
#
#data <- tar_read(pred_data)[[to_check]]$data
#models <- as.vector(unique(data$model_idx))
#cv <- as.vector(unique(data$cv_idx))
#types <- as.vector(unique(data$type))
#check_against <- "M_1"
#
#data_wide <- data %>% 
#  mutate(value = ifelse(value == "NaN", NA, value)) %>%  
#  pivot_wider(id_cols = c(run_type, cv_idx, run_idx, type), 
#              names_from = model_idx, values_from = value) %>%
#  as.data.frame()
#
#cases <- data %>% 
#  distinct(model_idx, cv_idx) %>% 
#  filter(model_idx != check_against) %>%
#  as.data.frame()
#
#out <- NULL
#for (row in 1:nrow(cases)){
#  idx <- as.character(cases[row, "model_idx"])
#  idx_cv <- as.character(cases[row, "cv_idx"])
#  for(typ in types){
#    data_1 <- data_wide %>% filter(cv_idx == idx_cv, type == typ) %>% pull(all_of(check_against))
#    data_2 <- data_wide %>% filter(cv_idx == idx_cv, type == typ) %>% pull(all_of(idx))
#    val <- t.test(data_1, data_2, paired = TRUE)$p.value
#    
#    res <- data.frame("cv" = idx_cv,
#                      "idx" = idx,
#                      "type" = typ,
#                      "p_val" = ifelse(val < 0.05, "sig", "non_sig"))
#    if(is.null(out)){
#      out <- res
#    } else {
#      out <- rbind(out, res)
#    }
#  }
#}
#
#out_wide <- out %>% 
#  pivot_wider(id_cols = c("cv", "idx"), 
#              names_from = "type", 
#              values_from = "p_val") %>%
#  arrange(cv)
#
#write.table(out_wide, 
#            sprintf("/proj/results/R/process_R_pred_data/t_tests_%s_%s.txt", 
#                    to_check, check_against), row.names = FALSE)
#

# For acr_data ------------------------------------------------------------
# run wise comparison between M_3 and CNN
#model_comp <- pred_data$plot_acr_sce_2$data %>%
#  filter(model_idx %in% c("acr_CNN", "M_3")) %>% 
#  distinct() %>% 
#  pivot_wider(id_cols = c("run_idx", "type", "scenario"), names_from = "model_idx", values_from = "value") %>% 
#  group_by(type, scenario) %>% 
#  mutate(delta = acr_CNN - M_3) %>%
#  summarize(total_runs = n(),
#            superior_cnn = sum(delta > 0),
#            competetive_cnn = sum(delta > -0.06 & delta < 0), 
#            performance_CNN = round((superior_cnn + competetive_cnn)/total_runs, 2),
#            .groups = "drop")
#
#write.table(model_comp, 
#            "/proj/results/R/process_R_pred_data/compare_models.txt", row.names = FALSE)
#
## mean prediction abilities per scenario
#pred_data$plot_acr_sce_2$data %>%
#  filter(model_idx %in% c("acr_CNN", "M_3")) %>% 
#  distinct() %>%
#  group_by(scenario, model_idx, type) %>% 
#  summarize(mean_val = mean(value), .groups = "drop") %>% 
#  pivot_wider(id_cols = c("scenario", "type"), names_from = "model_idx", values_from = "mean_val") %>%
#  arrange(type)
#
## correlation of prediction ability with training set
#pred_data$plot_acr_sce_2$data %>%
#  filter(model_idx %in% c("acr_CNN", "M_3")) %>% 
#  distinct() %>% group_by(model_idx, scenario, type) %>% 
#  summarize(cor_val = cor(train, value) , obs = n(), .groups = "drop") %>% 
#  pivot_wider(id_cols = c(scenario, type, obs), names_from = "model_idx", values_from = "cor_val") %>% 
#  arrange(type) %>% 
#  mutate(cnn_better = ifelse(acr_CNN > M_3, TRUE, FALSE))
#
## compare means across 5f and sce
#mean_5f <- pred_data$plot_acr_cv$data %>% 
#  group_by(name, type) %>%
#  summarize(mean_val = mean(value), .groups = "drop") %>%
#  filter(name %in% c("E-GBLUP_D", "acr_CNN"))
#
#mean_sce <- pred_data$plot_acr_sce_2$data %>% 
#  filter(model_idx %in% c("acr_CNN", "M_3")) %>% 
#  distinct() %>% group_by(model_idx, scenario, type) %>%
#  summarize(mean_val = mean(value), .groups = "drop") %>%
#  filter(scenario == 7)
#
# M_3 for hybrid
#100*((0.308 - 0.840)/0.840) ##-64
# M_3 for lines
#100*((0.382 - 0.755)/0.755) ## -50

# CNN for hybrid
#100*((0.190 - 0.737)/0.737) ## -74
# CNN for lines
#100*((0.415 - 0.680)/0.680) ## 39

# For wtn_data ------------------------------------------------------------
#pred_data <- tar_read(pred_data)
#"%!in%" <- Negate("%in%") 
#calculate_percent_differences <- function(data, cols_to_calc, base_col) {
#  base_col_sym <- sym(base_col)
#  data %>%
#    mutate(across(all_of(cols_to_calc), ~ (. - !!base_col_sym) / !!base_col_sym * 100, .names = "{col}_pct_diff"))
#}
#
#get_comparision_matrix <- function(data_in,
#                                   #cv_idx_in = "cv1",
#                                   #label_idx = "no_cgm_input",
#                                   what = "mean",
#                                   idx_1 = c("run_type", "run_idx", "cv_idx", "type"), 
#                                   idx_2 = c("run_type", "cv_idx", "type")){
#  data <- data_in %>% 
#    #filter(cv_idx == cv_idx_in, label == label_idx) %>%
#    pivot_wider(id_cols = idx_1, names_from = "model_idx", values_from = "value")
#  
#  models <- colnames(data)[which(colnames(data) %!in% idx_1)]
#  mat_line <- mat_hybrid <- matrix(NA, length(models), length(models),
#                                   dimnames = list(models, models))
#  
#  for (model_1 in 1:length(models)) {
#    for(model_2 in 1:length(models)){
#      if(model_2 > model_1){
#        out_0 <- calculate_percent_differences(data, cols_to_calc = models[model_2], 
#                                             base_col = models[model_1]) %>%
#          select(all_of(idx_1), contains("diff")) %>%
#          pivot_longer(cols = contains("diff"), names_to = "name", values_to = "value") 
#        if(what == "mean"){
#          out <- out_0 %>%
#            aggregate(as.formula(paste("value", paste(idx_2, collapse="+"), sep="~")), ., 
#                      function(i) round(mean(i), 3)) %>%
#            mutate(value_mod = ifelse(value > 1, round(value), NA))
#        } else if (what == "median"){
#          out <- out_0 %>%
#            aggregate(as.formula(paste("value", paste(idx_2, collapse="+"), sep="~")), ., 
#                      function(i) round(median(i), 3)) %>%
#            mutate(value_mod = ifelse(value > 1, round(value), NA))
#        }
#        
#        # check significance
#        p_lines <- t.test(data[data$type == "Lines", models[model_1]] %>% pull(models[model_1]),
#                          data[data$type == "Lines", models[model_2]] %>% pull(models[model_2]),
#                          paired = TRUE)$p.value
#        p_hybrid <- t.test(data[data$type == "Hybrid", models[model_1]] %>% pull(models[model_1]),
#                           data[data$type == "Hybrid", models[model_2]] %>% pull(models[model_2]),
#                           paired = TRUE)$p.value
#        
#        mat_line[model_1, model_2] <- sprintf("%s%s", 
#                                              round(out %>% filter(type == "Lines") %>% pull(value)),
#                                              ifelse(p_lines < 0.05, "*", "*ns"))
#        mat_hybrid[model_1, model_2] <- sprintf("%s%s",
#                                                round(out %>% filter(type == "Hybrid") %>% pull(value)),
#                                                ifelse(p_hybrid < 0.05, "*", "*ns"))
#      }
#    }
#  }
#
# return(list("line" = mat_line,
#             "hybrid" = mat_hybrid)) 
#}
#
#cv1 <- get_comparision_matrix(pred_data$wtn_cor_data %>%
#                                filter(cv_idx == "cv1") %>% 
#                                select(-label) %>% 
#                                distinct()
#                              #, "median"
#                              )
#cv2 <- get_comparision_matrix(pred_data$wtn_cor_data %>%
#                                filter(cv_idx == "cv2", label == "no_cgm_input") %>% 
#                                select(-label) %>% 
#                                distinct()
#                              #, "median"
#                              )
#cv3 <- get_comparision_matrix(pred_data$wtn_cor_data %>%
#                                filter(cv_idx == "cv3") %>% 
#                                select(-label) %>% 
#                                distinct()
#                              #, "median"
#                              )
#cv4 <- get_comparision_matrix(pred_data$wtn_cor_data %>%
#                                filter(cv_idx == "cv4", label == "no_cgm_input") %>% 
#                                select(-label) %>% 
#                                distinct()
#                              #, "median"
#                              )
#LoO <- get_comparision_matrix(pred_data$wtn_cor_data %>%
#                                filter(cv_idx == "cv2", label == "with_cgm_input") %>% 
#                                select(-label) %>% 
#                                distinct()
#                              #, "median"
#                              )
#cvL <- get_comparision_matrix(pred_data$wtn_cor_data %>%
#                                filter(cv_idx == "cv4", label == "with_cgm_input") %>% 
#                                select(-label) %>% 
#                                distinct()
#                              #, "median"
#                              )
#
# For data size  ----------------------------------------------------------
#cv_schema <- grid::rasterGrob(png::readPNG("/proj/results/R/generate_prediction_data/cross_validations.png"))
#cv_table <- grid::rasterGrob(png::readPNG("/proj/results/R/generate_prediction_data/cv_table.png"))
#data_tra <- do.call(rbind, lapply(data[grep("wtn", names(data), value = TRUE)],
#                                  function(x) x[["run_meta"]])) %>%
#  mutate(cv = gsub("\\S+_(cv\\d)|(LoO)_\\S+|(cvL)_\\S+", "\\1\\2\\3", run_id),
#         run = gsub("(\\S+)_cv\\d|LoO_(\\S+)|cvL_(\\S+)", "\\1\\2\\3", run_id),
#         run = gsub("run_", "", run)) 
#
#plot_cvs <- data_tra %>%
#  filter(cv %!in% c("cvL")) %>%
#  distinct(cv, run, train, test) %>%
#  pivot_longer(cols = c("train", "test"),
#               names_to = "type",
#               values_to = "value") %>%
#  convert(num(run, value),
#          fct(cv),
#          fct(type, .args = list(levels = c("train", "test")))) %>%
#  ggplot(aes(x = run, y = value, color = type)) +
#  geom_line() +
#  facet_wrap(~ cv) +
#  labs(y = "Data points", x = "Runs") +
#  scale_color_manual(name = "Point type", 
#                     values = c("train" = "#fc6020",
#                                "test" = "#e9c01c"),
#                     labels = c("train" = "Train",
#                                "test" = "Test")) + 
#  theme_classic(base_size = 10) +
#  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))
#
#layout <- "
#AAAAAA
#BBBBCC
#"
#joint_plot <- wrap_elements(cv_schema) +
#  wrap_elements(plot_cvs) + 
#  wrap_elements(cv_table) +
#  plot_layout(design = layout
#              #, heights = c(1, 1.25)
#              )
#write_at <- "/proj/tmp_data"
#ggsave(plot = joint_plot, filename = sprintf("%s/cv_plot.png", write_at), 
#       width = 17, height = 17, units = "cm", dpi = 600)
# Reruns ------------------------------------------------------------------
#neg_cnn <- pred_data$wtn_cor_data %>% filter(cv_idx == "cv1", model_idx == "CNN_EC") %>%
#  filter(value < 0.1) %>%
#  count(run_idx, cv_idx, model_idx, label)
#miss <- neg_cnn
#
#miss <- pred_data$wtn_cor_data %>% filter(is.na(value)) %>% count(run_idx, cv_idx, model_idx, label)
#no_file <- pred_file_paths$pred_data_meta %>% filter(!pred_check)
#
#source_dir <- "/proj/results/Py/create_slurm_scripts/wtn_CNN_EC_wtn_tra"
#file_name <- "/pred/output.csv"
#master_script <- sprintf("%s/master_script.sh", source_dir)
#
#for(row in 1:nrow(miss)){
#  run <- miss %>% filter(row_number() == row) %>% 
#    pull(run_idx) %>% as.vector()
#  cv <- miss %>% filter(row_number() == row) %>% 
#    pull(cv_idx) %>% as.vector()
#  dir <- paste0(run, "_", cv)
#  path <- sprintf("%s/%s/%s", source_dir, dir, file_name)
#  path_model <- sprintf("%s/%s/model", source_dir, dir)
#  path_cb <- sprintf("%s/%s/callback_data", source_dir, dir)
#  # for the pred file
#  if(file.exists(path)){
#    file.remove(path)
#  } else {
#    sprintf("%s does not exist", dir)
#  }
#  # for the model directory
#  if(dir.exists(path_model)){
#    system(sprintf("rm -rf %s", path_model))
#  } else {
#    sprintf("model dir is not present in %s", dir)
#  }
#  if(!dir.exists(path_model)){
#    system(sprintf("mkdir -p %s", path_model))
#  } else {
#    sprintf("model dir is already here %s", path_model)
#  }
#  # for the callback directory
#  if(dir.exists(path_cb)){
#    system(sprintf("rm -rf %s", path_cb))
#  } else {
#    sprintf("callback dir is not present in %s", dir)
#  }
#}
