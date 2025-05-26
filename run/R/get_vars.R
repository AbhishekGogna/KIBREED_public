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

# pcent increases
#acr_cv_plot$data %>% 
#  pivot_wider(id_cols = c("type", "run_idx"), names_from = "model_idx", values_from = "value") %>% 
#  group_by(type) %>% summarize(mean_diff = mean((M_3- M_2)/M_2))


#mod_plot <- var_acr$plot_vars_pheno_data_acr + 
#  scale_x_discrete(labels=c("M_1" = "GBLUP", "M_2" = "GBLUP_D", "M_3" = "E-GBLUP_D")) + 
#  ylab("Proportion explained") + 
#  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#ggsave(plot = mod_plot, filename = "/proj/results/R/get_vars/pheno_data_acr_var_plot_mod.png", 
#       width = 8.4, height = 8.4, units = "cm", dpi = 600)
#to_table <- var_wtn$plot_vars_pheno_data_wtn$data %>% 
#  mutate(value = round(value, 2)) %>%
#  pivot_wider(id_cols = "type", names_from = "model", values_from = "value")
#baseline <- to_table %>% filter(type == "G_i@BRR_var") %>% pull(M_1)
#var_wtn$plot_vars_pheno_data_wtn$data %>%
#  mutate(value = round(value, 2)) %>%
#  pivot_wider(id_cols = "model", names_from = "type", values_from = "value") %>%
#  filter(model %in% c("M_2", "M_3", "M_5")) %>%
#  group_by(model) %>%
#  summarise(var_g = (`G_a@BRR_var` + `G_d@BRR_var` + `G_aa@BRR_var`)/baseline)
#var_wtn$plot_vars_pheno_data_wtn$data %>%
#  mutate(value = round(value, 2)) %>%
#  pivot_wider(id_cols = "model", names_from = "type", values_from = "value") %>%
#  filter(model %in% c("M_4", "M_6")) %>%
#  group_by(model) %>%
#  summarise(var_g = sum(`G_a@BRR_var`, `G_d@BRR_var`, `G_aa@BRR_var`,
#                        `G_a_ERM_l@RKHS_var`, `G_a_ERM_nl@RKHS_var`, na.rm = TRUE)/baseline)
#write.csv(to_table, "/proj/tmp_data/table_var.csv", row.names = F)

# for institute days poster
#data_ins <- p2_data %>% filter(type == "Lines") %>%
#  convert(chr(name)) %>%
#  mutate(name = case_when(
#    name == "E-GBLUP_D" ~ "GBLUP",
#    name == "acr_CNN" ~ "CNN",
#    .default = name
#  )) %>%
#  convert(fct(name, .args = list(levels = c("GBLUP", "CNN"))))
#
#p2 <- ggplot(aes(x = train_val, y = value), data = data_ins) +
#  geom_point() +
#  coord_cartesian(ylim = c(-0.25, 1)) +
#  facet_wrap(~name, ncol = 1) +
#  #ylab("Prediction ability") +
#  #xlab("Training set size") +
#  geom_smooth(method = "loess", color = "blue") +
#  stat_cor(aes(label =paste(after_stat(r.label), cut(after_stat(p), 
#                                                     breaks = c(-Inf, 0.05, Inf),
#                                                     labels = c("'*'", "''")), 
#                            sep = "~")),
#           p.accuracy = 0.01, r.accuracy = 0.01, 
#           label.x = 0, label.y = 0.9, size = 5) +
#  theme_classic(base_size = 12) +
#  scale_x_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE)) +
#  theme(legend.position = "none") +
#  scale_fill_manual(values = c("grey", "white"))
#
#ggsave(p2, filename = sprintf("%s/%s", "/proj/tmp_data", "pred_acr_res.png"),
#       width = 8.4, height = 12.6, units = "cm", dpi = 600, bg = "white")
#
#

