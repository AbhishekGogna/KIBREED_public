# R Statistical Analysis

Traditional genomic prediction using GBLUP and Bayesian models with R/BGLR.

## Available Scripts

1. **`scr_genomic_prediction_acr.R`** - Across environments genomic prediction
2. **`scr_genomic_prediction_wtn.R`** - Within environments G×E genomic prediction

## Example 1: Genomic Prediction Across Environments (ACR)

**Script:** `scr_genomic_prediction_acr.R`

Performs genomic predictions with genotypic values averaged across environments using a simple 80/20 random cross-validation split.

```bash
# Run the script
Rscript scr_genomic_prediction_acr.R
```

### Models Included
- **Model 1:** Additive + Dominance effects
- **Model 2:** Additive + Dominance + Epistatic effects

### Output Files
- **`results/prediction_results_acr.qs`** - Detailed prediction results
- **`results/accuracy_summary_acr.qs`** - Accuracy summary by genotype type
- **`logs/genomic_prediction_acr.log`** - Analysis log file

## Example 2: Genomic Prediction Within Environments (WTN)

**Script:** `scr_genomic_prediction_wtn.R`

Performs complex multi-environment G×E genomic predictions using CV3 cross-validation scenario with modular ETA components.

```bash
# Run the script
Rscript scr_genomic_prediction_wtn.R
```

### Models Available
The script defines 6 different model specifications:
- **M1:** Basic environment + genotype effects
- **M2:** Environment + genomic relationships (A+D+AA)
- **M3:** Linear ERM + genomic relationships
- **M4:** Linear ERM + genomics + G×E (additive×linear)
- **M5:** Nonlinear ERM + genomic relationships
- **M6:** Nonlinear ERM + genomics + G×E (additive×nonlinear)

*Note: The script runs M1 and M6 as examples. Model specifications are saved for running additional models.*

### Output Files
- **`results/prediction_results_wtn.qs`** - Detailed prediction results
- **`results/accuracy_summary_wtn.qs`** - Accuracy summary by genotype type
- **`logs/genomic_prediction_wtn.log`** - Analysis log file

## Data Requirements for R Scripts

R scripts expect the following data files in the `data/` directory:
- `pheno_acr.qs` / `pheno_wtn.qs` - Phenotype data
- `g_data_add.qs` - Additive genomic data
- `g_data_dom.qs` - Dominance genomic data
- `ev_data.qs` - Environmental data (WTN only)
