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

# Inside R
proj_list <- names(yaml::read_yaml("_targets.yaml"))
proj_list
# [1] "generate_prediction_data" "process_R_pred_data"
# [3] "get_vars"                 "feature_importance"

# Set the subproject you want to run
Sys.setenv("TAR_PROJECT" = proj_list[1])

# Load the targets package
library(targets)

# Confirm the project name
Sys.getenv("TAR_PROJECT")

# Run the pipeline
tar_make()

```

## Example code 1: Genomic prediction for genotypic values averaged across environments

## Example code 2: Genomic prediction for genotypic values within environments

## Example code 3: Generating feature importance scores for core set

## License

This project is licensed under the MIT License.  
See the [LICENSE](LICENSE) file for full details.
