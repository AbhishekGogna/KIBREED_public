# KIBREED_public

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

Once completed, your R environment will be automatically activated when you start R. If that does not happen check the .Rprofile file.

### Setup Notes

- Run all commands from the project root directory
- If you encounter version mismatch errors, ensure the correct versions are installed and available in your `$PATH`
- The input data is available at [FigShare](https://figshare.com/). Download those files and put them at `/path/to/KIBREED_public/data`

## Running Genomic Prediction Scripts

After environment setup, you can run the genomic prediction analyses directly using the provided R scripts.

### Required Directory Structure

```r
# Start R from the project root
R

# Create required directories
dirs <- c("results", "tmp", "logs")
sapply(dirs, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE)
})
```

### Available Scripts

1. **`scr_genomic_prediction_acr.R`** - Across environments genomic prediction
2. **`scr_genomic_prediction_wtn.R`** - Within environments G×E genomic prediction

## Analysis Examples

### Example 1: Genomic Prediction Across Environments (ACR)

**Script:** `scr_genomic_prediction_acr.R`

Performs genomic predictions with genotypic values averaged across environments using a simple 80/20 random cross-validation split.

```bash
# Run the script
Rscript scr_genomic_prediction_acr.R
```

#### Models Included
- **Model 1:** Additive + Dominance effects
- **Model 2:** Additive + Dominance + Epistatic effects

#### Output Files
- **`results/prediction_results_acr.qs`** - Detailed prediction results
- **`results/accuracy_summary_acr.qs`** - Accuracy summary by genotype type
- **`logs/genomic_prediction_acr.log`** - Analysis log file

### Example 2: Genomic Prediction Within Environments (WTN)

**Script:** `scr_genomic_prediction_wtn.R`

Performs complex multi-environment G×E genomic predictions using CV3 cross-validation scenario with modular ETA components.

```bash
# Run the script
Rscript scr_genomic_prediction_wtn.R
```

#### Models Available
The script defines 6 different model specifications:
- **M1:** Basic environment + genotype effects
- **M2:** Environment + genomic relationships (A+D+AA)
- **M3:** Linear ERM + genomic relationships
- **M4:** Linear ERM + genomics + G×E (additive×linear)
- **M5:** Nonlinear ERM + genomic relationships
- **M6:** Nonlinear ERM + genomics + G×E (additive×nonlinear)

*Note: The script runs M1 and M6 as examples. Model specifications are saved for running additional models.*

#### Output Files
- **`results/prediction_results_wtn.qs`** - Detailed prediction results
- **`results/accuracy_summary_wtn.qs`** - Accuracy summary by genotype type
- **`results/model_specifications_wtn.qs`** - All model definitions for future use
- **`logs/genomic_prediction_wtn.log`** - Analysis log file

### Data Requirements

Both scripts expect the following data files in the `data/` directory:
- `pheno_acr.qs` / `pheno_wtn.qs` - Phenotype data
- `g_data_add.qs` - Additive genomic data
- `g_data_dom.qs` - Dominance genomic data
- `ev_data.qs` - Environmental data (WTN only)

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for full details.
