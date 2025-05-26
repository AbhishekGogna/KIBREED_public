# Portable R Project Setup Script
# This script sets up a complete R project environment without container dependencies

# Get current working directory as project directory
project_dir <- getwd()
cat("Setting up project in:", project_dir, "\n")

# Create project structure
cat("Creating project directory structure...\n")

# Define core directories
core_dirs <- c("src", "run", "results", "tmp_data", "logs", "raw_data", "store")

# Create main directories and R/Python subdirectories
my_dirs <- list()
for (dir_name in core_dirs) {
  # Create main directory
  dir_path <- file.path(project_dir, dir_name)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Create R and Python subdirectories for relevant folders
  if (dir_name %in% c("src", "run")) {
    r_path <- file.path(dir_path, "R")
    py_path <- file.path(dir_path, "Py")
    
    if (!dir.exists(r_path)) dir.create(r_path, recursive = TRUE)
    if (!dir.exists(py_path)) dir.create(py_path, recursive = TRUE)
    
    my_dirs[[paste0(dir_name, "_R")]] <- r_path
    my_dirs[[paste0(dir_name, "_Py")]] <- py_path
  }
}

# Store directory paths for later use
jsonlite::write_json(my_dirs, file.path(project_dir, "results", "core_paths.json"))

# Initialize renv for dependency management
cat("Initializing renv for dependency management...\n")

# Check if renv is installed, install if not
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Initialize renv in the project
renv::init(project = project_dir, 
           bare = TRUE,
           settings = list(use.cache = TRUE),  # Enable cache for faster installs
           force = TRUE)

# Modify .Rprofile to activate renv when in project directory
rprofile_content <- sprintf('if(getwd() == "%s") source("renv/activate.R")\n', project_dir)
cat(rprofile_content, file = file.path(project_dir, ".Rprofile"))

# Define required packages
cat("Installing required R packages...\n")

# Core packages with specific versions where needed
required_packages <- c(
  "remotes",           # for version-specific installs
  "tidyverse",         # data wrangling suite
  "hablar",            # additional data manipulation
  "readxl",            # Excel file reading
  "reshape2",          # data reshaping
  "targets",           # reproducible analysis pipeline
  "tarchetypes",       # additional targets functions
  "visNetwork",        # network visualization
  "usethis",           # project setup utilities
  "rmarkdown",         # dynamic documents
  "knitr",             # document generation
  "pbkrtest",          # statistical testing
  "ggpubr",            # publication-ready plots
  "foreach",           # parallel processing
  "doParallel",        # parallel backend
  "RhpcBLASctl",       # BLAS control
  "qs",                # fast serialization
  "jsonlite",          # JSON handling
  "feather",           # fast data frames
  "BGLR",              # Bayesian analysis
  "corehunter"         # core collection optimization
)

# Install packages that aren't already installed
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if (length(new_packages) > 0) {
  cat("Installing packages:", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, Ncpus = parallel::detectCores() - 1)
}

# Special handling for packages with version requirements
version_specific <- list(
  "tidyverse" = "1.3.2",
  "pbkrtest" = "0.5.1",
  "ggpubr" = "0.5.0",
  "corehunter" = "3.2.2"
)

# Install version-specific packages using remotes
library(remotes)
for (pkg in names(version_specific)) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing", pkg, "version", version_specific[[pkg]], "\n")
    install_version(pkg, version = version_specific[[pkg]], upgrade = "never")
  }
}

# Setup for ASReml (optional commercial package)
cat("Setting up ASReml installation directory...\n")
cellar_dir <- file.path(project_dir, "renv", "cellar")
if (!dir.exists(cellar_dir)) {
  dir.create(cellar_dir, recursive = TRUE)
}

cat("Note: To install ASReml, place the ASReml tar.gz file in:", cellar_dir, "\n")
cat("The file should be named 'asreml_*.tar.gz'\n")

# Check for and install ASReml if tar file is present
asreml_files <- list.files(cellar_dir, pattern = "asreml.*\\.tar\\.gz", full.names = TRUE)
if (length(asreml_files) > 0) {
  cat("Found ASReml file:", basename(asreml_files[1]), "\n")
  cat("Installing ASReml...\n")
  install.packages(asreml_files[1], repos = NULL, type = "source")
} else {
  cat("No ASReml installation file found. Skipping ASReml installation.\n")
}

# Setup Python environment
cat("Setting up Python environment...\n")

# Create requirements.txt for Python dependencies
python_requirements <- "# Deep Learning packages
tensorflow==2.8.4
tensorboard==2.8.4
pyarrow==5.0.0
matplotlib==3.5.1
pandas==1.4
scikit-learn==1.0.2
patsy==0.5.2
protobuf==3.19.6
keras-tuner==1.1.3
ipykernel==6.22.0

# Machine Learning packages
xgboost==2.0.0

# Workflow management
graphviz==0.20
doit==0.36.0
pygraphviz==1.9
import_deps==0.2.0
"

writeLines(python_requirements, file.path(project_dir, "requirements.txt"))

# Create Python setup instructions
python_setup <- "# Python Environment Setup Instructions

# 1. Create a virtual environment:
python3 -m venv py_env

# 2. Activate the environment:
# On Linux/Mac:
source py_env/bin/activate
# On Windows:
# py_env\\Scripts\\activate

# 3. Install required packages:
pip install -r requirements.txt

# 4. (Optional) Add to Jupyter if you use it:
python -m ipykernel install --user --name=project_env

# 5. To activate automatically, add this to your shell profile:
# echo 'source $(pwd)/py_env/bin/activate' >> ~/.bashrc
"

writeLines(python_setup, file.path(project_dir, "python_setup_instructions.txt"))

# Initialize targets workflow
cat("Setting up targets workflow...\n")

if (!requireNamespace("targets", quietly = TRUE)) {
  install.packages("targets")
}

library(targets)
use_targets(open = FALSE)

# Create subdirectories for specific analysis types
analysis_dirs <- c(
  "generate_prediction_data",
  "process_R_pred_data", 
  "get_vars",
  "feature_importance"
)

# Create subdirectories within run and src
for (analysis in analysis_dirs) {
  # Create subdirectory in run/R
  run_subdir <- file.path(my_dirs[["run_R"]], analysis)
  if (!dir.exists(run_subdir)) {
    dir.create(run_subdir, recursive = TRUE)
  }
  
  # Create subdirectory in src/R  
  src_subdir <- file.path(my_dirs[["src_R"]], analysis)
  if (!dir.exists(src_subdir)) {
    dir.create(src_subdir, recursive = TRUE)
  }
  
  # Create store subdirectory
  store_subdir <- file.path(project_dir, "store", analysis)
  if (!dir.exists(store_subdir)) {
    dir.create(store_subdir, recursive = TRUE)
  }
  
  cat("Created directories for:", analysis, "\n")
}

# Create .gitignore
cat("Creating .gitignore file...\n")
gitignore_content <- "# R specific
.Rproj.user/
.Rhistory
.RData
.Ruserdata
*.Rproj

# renv
renv/library/
renv/python/
renv/staging/

# Python
py_env/
__pycache__/
*.pyc
*.pyo
*.egg-info/

# Temporary data
tmp_data/
*.tmp

# Logs
logs/*.log

# OS specific
.DS_Store
Thumbs.db

# IDE
.vscode/
.idea/

# Large data files (uncomment if needed)
# raw_data/
# results/*.rds
# results/*.qs
"

writeLines(gitignore_content, file.path(project_dir, ".gitignore"))

# Create README
cat("Creating README file...\n")
readme_content <- paste0("# ", basename(project_dir), " Project

This project uses renv for R dependency management and includes setup for both R and Python environments.

## Setup Instructions

### R Environment
1. Open R/RStudio in this directory
2. Run `source('setup_project.R')` to initialize everything
3. Use `renv::restore()` to install exact package versions

### Python Environment  
1. Follow instructions in `python_setup_instructions.txt`
2. Or simply run: `python3 -m venv py_env && source py_env/bin/activate && pip install -r requirements.txt`

### Targets Workflow
- Analysis scripts are in `run/R/`
- Helper functions are in `src/R/`
- Use `targets::tar_make()` to run analyses
- Use `targets::tar_visnetwork()` to visualize workflow

## Project Structure
- `src/`: Source code and functions
- `run/`: Analysis scripts and workflows  
- `results/`: Analysis outputs
- `raw_data/`: Input data files
- `tmp_data/`: Temporary files
- `logs/`: Log files
- `store/`: Targets data store

## Notes
- ASReml installation requires placing the tar.gz file in `renv/cellar/`
- All paths are relative to project directory for portability
")

writeLines(readme_content, file.path(project_dir, "README.md"))

# Update renv snapshot
cat("Creating renv snapshot...\n")
renv::snapshot(force = TRUE)

# Final status check
cat("\n=== Project Setup Complete ===\n")
cat("Project directory:", project_dir, "\n")
cat("Required packages installed:", length(required_packages), "\n")
cat("Python requirements file created\n")
cat("Targets workflow initialized\n")
cat("Analysis directories created:", length(analysis_dirs), "\n")

cat("\nNext steps:\n")
cat("1. Add your data files to 'raw_data/'\n") 
cat("2. Create your analysis scripts in the appropriate subdirectories:\n")
cat("   - 'run/R/generate_prediction_data/'\n")
cat("   - 'run/R/process_R_pred_data/'\n") 
cat("   - 'run/R/get_vars/'\n")
cat("   - 'run/R/feature_importance/'\n")
cat("3. Add your functions to corresponding 'src/R/' subdirectories\n")
cat("4. Set up Python environment using 'python_setup_instructions.txt'\n")
cat("5. Use targets::tar_make() to run your analyses\n")

cat("\nProject setup completed successfully!\n")