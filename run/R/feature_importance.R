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

list(
  tar_target(
    name = env_clusters,
    command = get_env_clusters(existing_data_path = sprintf("%s/%s", project_path, "data/data_pub.qs"))
  )
  , tar_target(
    name = feature_importance,
    command = get_feature_imp_plots(existing_data_path = sprintf("%s/%s", project_path, "data/data_pub.qs"),
                                    scores_at = "/proj/results/Py/feature_importance/feature_imp_scores.csv")
  
  )
)