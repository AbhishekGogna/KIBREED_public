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
  #, tar_target(
  #  name = overviews,
  #  command = make_overview_plots(existing_data = data_paths,
  #                                data1 = RD_data,
  #                                data2 = core_set,
  #                                data3 = env_clusters,
  #                                data4 = feature_importance,
  #                                data5 = training_data_2)
  #)
)

# identify best and worst performers per environment cluster
#core_dist_series <- core_set$core_dist_series
##write.csv2(core_dist_series, "/proj/tmp_data/core_dist_series.csv", row.names = F, quote = T)
#per_series <- core_dist_series %>% select(-freq_env_class, -total) %>% as.matrix()
#per_series_sum <- colSums(per_series, na.rm = T)
#per_series_sum[order(per_series_sum, decreasing = T)]

## how many of the core set are in the extremes
#dist_plot <- core_set$core_ext_plot
#ggsave(plot = dist_plot,
#       filename = sprintf("%s/dist_plot.png", "/proj/results/R/feature_importance"), 
#       width = 8.4, height = 8.4, dpi = 600, units = "cm")
#
# a possible line plot
#data <- predicted_data[["predicted_data"]][["pred"]]%>% filter(is.na(set)) 
#to_color <- data %>% filter(geno %in% c("zf003zm027", "hystarhystar"))
#data %>% 
#  ggplot(aes(x = env, y = pred, group = geno)) +
#  geom_line() +
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#  geom_line(aes(x = env, y = pred, group = geno, color = "red"), data = to_color)  