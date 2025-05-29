# Get the project directory
project_dir <- getwd()
lockfile_path <- file.path(project_dir, "renv.lock")

cat("Checking for renv.lock in:", lockfile_path, "\n")

# Check if renv.lock exists
if (!file.exists(lockfile_path)) {
  stop("Error: renv.lock file not found. Please provide a valid renv.lock to proceed.\n")
}

# Install renv if not already installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv in bare mode if not already done
if (!file.exists(file.path(project_dir, "renv"))) {
  renv::init(project = project_dir, bare = TRUE, force = TRUE)
}

# Restore the environment from the lock file
cat("Restoring R packages from renv.lock...\n")
Sys.setenv(LC_ALL = "C")
options(install.packages.compile.from.source = "never")
renv::restore(lockfile = lockfile_path, prompt = FALSE)

cat("All packages restored successfully from renv.lock.\n")

# Some house keeping
if (dir.exists("/proj")) {
  project_path <- "/proj"
} else {
  project_path <- getwd()
}

# Create YAML content with correct paths
yaml_content <- paste0(
  "generate_prediction_data:\n",
  "  script: ", project_path, "/run/R/generate_prediction_data.R\n",
  "  store: ", project_path, "/store/generate_prediction_data\n",
  "process_R_pred_data:\n",
  "  script: ", project_path, "/run/R/process_R_pred_data.R\n",
  "  store: ", project_path, "/store/process_R_pred_data\n",
  "env_clusters:\n",
  "  script: ", project_path, "/run/R/env_clusters.R\n",
  "  store: ", project_path, "/store/env_clusters\n",
  "\n"  # Empty line at the end
)

# Write to YAML file at project path
writeLines(yaml_content, file.path(project_path, "_targets.yaml"))