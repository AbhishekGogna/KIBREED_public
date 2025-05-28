# KIBREED

A genomic prediction project for plant breeding applications using R and Python environments with specific version requirements.

## Prerequisites

- **R**: Version 4.0.5
- **Python**: Version 3.8.11

## Environment Setup

### 1. Verify Language Versions

Ensure R 4.0.5 and Python 3.8.11 are installed and properly configured in your PATH (likely by adding their paths to `.bashrc`).

```bash
./check_Langs.sh
```

### 2. Python Environment Setup

Create a virtual environment named `py_env` and install all required Python packages from `pyenv.lock`:

```bash
./setup_py.sh
```

**Activate the environment:**
```bash
source /path/to/KIBREED_public/py_env/bin/activate
```

### 3. R Environment Setup

Restore the R package environment using `renv` and the `renv.lock` file:

```bash
Rscript setup.R
```

Once completed, your R environment will be automatically activated when you start R.

### Setup Notes

- Run all commands from the project root directory
- If you encounter version mismatch errors, ensure the correct versions are installed and available in your `$PATH`
- The input data is available at Dryad. Download those files and put them at `/path/to/KIBREED_public/data`

## Running R Subprojects with {targets}

After environment setup, use the {targets} pipeline system to run subprojects:

```r
# Start R
R

# Set basic directory structure
project_path <- "/proj" 

# Note: Each subproject runs from /proj root directory, mapped to KIBREED_public
# You may need to modify run scripts accordingly

# Create required directories
dirs <- c("results/R", "results/Py", "logs/R", "logs/Py", "tmp_data/R", "tmp_data/Py")
sapply(file.path(project_path, dirs), function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
})

# View available subprojects
proj_list <- names(yaml::read_yaml("_targets.yaml"))
# Available projects: "generate_prediction_data", "process_R_pred_data", "get_vars", "feature_importance"

# Set the subproject to run
Sys.setenv("TAR_PROJECT" = proj_list[1])

# Load targets package
library(targets)

# View planned targets
tar_manifest() # Each name is a target to be run

# Run complete pipeline
tar_make() 

# Run specific target
tar_make(names = "target_name") # Use names from tar_manifest()
```

## Analysis Examples

### Example 1: Genomic Prediction Across Environments

**Project:** `generate_prediction_data`

Generates files for genomic predictions with genotypic values averaged across environments.

**Output Location:** `/proj/results/R/generate_prediction_data/cv_acr*`

#### Output Files

- **`cv_acr_5f.json`** - Train/test splits data
- **`eigen_data/`** - Eigen decompositions used for predictions
- **`master_files/`** - R scripts required for running predictions (must be executed to generate output)
- **`run_data/`** - Directories for storing prediction results

**Related Projects:**
- **`process_R_pred_data`** - Processes prediction output and generates visualization figures
- **`get_vars`** - Calculates variances and generates corresponding figures

### Example 2: Genomic Prediction Within Environments

**Project:** `generate_prediction_data`

Generates files for genomic predictions with genotypic values within specific environments.

**Output Location:** `/proj/results/R/generate_prediction_data/cv_wtn_tra`

#### Output Files

- **`cv_wtn_tra.json`** - Train/test splits data
- **`cv_wtn_tra_meta.txt`** - Metadata file containing:
  - `connect` - CV identifier (combination of "cv" + "run")
  - `train` - Data points in training set
  - `test` - Data points in test set
  - `train_env` - Environments in training set
  - `train_geno` - Genotypes in training set
  - `test_env` - Environments in test set
  - `test_geno` - Genotypes in test set
- **`cv_wtn_tra_sizes.png`** - Visualization of training and test set data points
- **`master_files/`** - R scripts required for running predictions
- **`run_data/`** - Directories for storing prediction results

#### Important Notes

Console messages like *"run_50_cv1 has 6 environments with low number of genotypes"* indicate that after the 80:20 quadrant 1 split, some test environments contain fewer than 50 genotypes. Check logs for detailed counts per environment. In testing, these typically contained more than 40 genotypes, which is acceptable.

**Related Project:**
- **`process_R_pred_data`** - Processes prediction output and generates visualization figures

### Example 3: Environment Clustering Analysis

**Project:** `feature_importance`

Uses predicted values of the core set to derive clusters of environments based on predicted GxE (Genotype × Environment) patterns.

**Output Location:** `/proj/results`

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for full details.
