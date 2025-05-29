library(targets)

# Set basic directory structure
if (dir.exists("/proj")) {
  project_path <- "/proj"
} else {
  project_path <- getwd() # assumes that this script is run from KIBREED_public
}

# Define packages needed and default storage format for intermittent data
tar_option_set(packages = c("dplyr", "ggplot2", "qs", 
                            "cluster", "RColorBrewer", 
                            "ggpubr"),
               format = "qs")

run_name <- "env_clusters"

tar_source(sprintf("/%s/src/R/fun_%s.R", project_path, run_name))

list(
  tar_target(
    name = env_clusters,
    command = get_env_clusters(existing_data_path = sprintf("%s/%s", project_path, "data"),
                               write_at = sprintf("%s/results/R/%s", project_path, run_name),
                               log_at = sprintf("%s/logs/R/%s", project_path, run_name))
  )
)