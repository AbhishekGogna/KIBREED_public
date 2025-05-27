# KIBREED

This project requires R 4.0.5 and Python 3.8.11 with specific package environments managed by lock files.

## Steps to set up the working environments. 

1. Check R and Python versions
Make sure you have R 4.0.5 and Python 3.8.11 installed and configured, likely by adding their paths to .bashrc

```bash
./check_Langs.sh
```

2. Set up Python environment
This script creates a virtual environment named py_env and installs all required Python packages from pyenv.lock.

```bash
./setup_py.sh
```

Once this is done you can activate your environment with source /path/to/KIBREED_public/py_env/bin/activate

3. Set up R environment
This script restores the R package environment using renv and the renv.lock file.

```bash
Rscript setup.R
```

Once this is done your environment should be activated when you enter R

Notes
1. Run all commands from the project root directory.
2. If you get version mismatch errors, ensure the correct versions are installed and available in your $PATH.

4. Running R subprojects using {targets}

Once the environment is set up, you can start R and run subprojects using the {targets} pipeline system.

```r
# Start R
R

# Set basic directory structure
project_path <- "/proj" 

# Please note that each sub projects tries to run from /proj root directory. This 
#is mapped to KIBREED_public. You may need to modify run scripts accordingly. 

dirs <- c("results/R", "results/Py", "logs/R", "logs/Py", "tmp_data/R", "tmp_data/Py")
sapply(file.path(project_path, dirs), function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
})

# Run sub projects
proj_list <- names(yaml::read_yaml("_targets.yaml"))
# [1] "generate_prediction_data" "process_R_pred_data"
# [3] "get_vars"  [4] "feature_importance"

# Set the subproject you want to run
Sys.setenv("TAR_PROJECT" = proj_list[1])

# Load the targets package
library(targets)

# View planned targets
tar_manifest() # each name is a target to be run

# Run the sub project pipeline completely
tar_make() 

# Run the sub project pipeline one target at a time
tar_make(names = "target_name") # from names in tar_manifest

```

## Example code 1: Genomic prediction for genotypic values averaged across environments

**generate_prediction_data** - Generates files for genomic predictions. Output stored in `/proj/results`.
    you need to run the files to generate the output, which is stored in /proj/results

**process_R_pred_data** - Processes prediction output and generates figures.

**get_vars** - Calculates variances and generates corresponding figures.

## Example code 2: Genomic prediction for genotypic values within environments

## generate_prediction_data

Generates files for genomic predictions with output stored in `/proj/results/R/generate_prediction_data/cv_wtn_tra`

### Output Files

**cv_wtn_tra.json** - Stores train/test splits data

**cv_wtn_tra_meta.txt** - Contains metadata for cv_wtn_tra.json with the following fields:
- `connect` - CV identifier (combination of "cv" + "run")  
- `train` - Data points in training set
- `test` - Data points in test set
- `train_env` - Environments in training set
- `train_geno` - Genotypes in training set
- `test_env` - Environments in test set
- `test_geno` - Genotypes in test set

**cv_wtn_tra_sizes.png** - Visualization of training and test set data points from cv_wtn_tra_meta

**master_files** - R scripts required for running predictions (must be executed to generate output)

**run_data** - Directories for storing prediction results

Note: Console messages like "run_50_cv1 has 6 environments with low number of genotypes" indicate that after the 80:20 quadrant 1 split, some test environments contain fewer than 50 genotypes. Check logs for detailed counts per environment. In testing, these typically contained more than 40 genotypes, which is acceptable.

## process_R_pred_data

Processes prediction output and generates visualization figures.

## Example code 3: Clustering environments based on predicted GxE patterns

**feature_importance** - Uses the predicted values of the core set to derive 
clusters of environments. Results are stored at `/proj/results`

## License

This project is licensed under the MIT License.  
See the [LICENSE](LICENSE) file for full details.
