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
proj_list
# [1] "generate_prediction_data" "process_R_pred_data"
# [3] "get_vars"

# Set the subproject you want to run
Sys.setenv("TAR_PROJECT" = proj_list[1])

# Load the targets package
library(targets)

# Confirm the project name
Sys.getenv("TAR_PROJECT")

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

**generate_prediction_data** - Generates files for genomic predictions. Output stored in `/proj/results`
    you need to run the files to generate the output, which is stored in `/proj/results`

**process_R_pred_data** - Processes prediction output and generates figures.

## Example code 3: Clustering environments based on predicted GxE patterns

**feature_importance** - Uses the predicted values of the core set to derive 
clusters of environments. Results are stored at `/proj/results`

## License

This project is licensed under the MIT License.  
See the [LICENSE](LICENSE) file for full details.
