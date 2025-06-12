# Python CNN Analysis

Deep learning genomic prediction using Convolutional Neural Networks with TensorFlow/Keras.

## Main Script

**`scr_genomic_prediction.Py`** - CNN-based genomic prediction with automatic hyperparameter tuning

This script provides deep learning-based genomic prediction using Convolutional Neural Networks with the following features:

## CNN Models Available

The `Py/` directory contains modular CNN architectures:

- **`model_acr_CNN.py`** - CNN for across-environment predictions
- **`model_CNN_EV.py`** - CNN with environmental covariates for within-environment predictions

## Python Module Structure

```
Py/
├── libs.py           # Core libraries and imports
├── func.py           # Utility functions (data loading, scaling, etc.)
├── model_acr_CNN.py  # ACR CNN architecture
└── model_CNN_EV.py   # WTN CNN with environmental data
```

## Running CNN Analysis

1. **Activate Python environment:**
```bash
source py_env/bin/activate
```

2. **Configure prediction type in script:**
```python
# Edit scr_genomic_prediction.Py
PREDICTION_TYPE = "acr"  # or "wtn"
```

3. **Run CNN analysis:**
```bash
python scr_genomic_prediction.Py
```

## CNN Output Structure

Results are organized under `results/CNN/[prediction_type]/`:

```
results/CNN/
├── acr/
│   ├── model/
│   │   ├── model_tuned/     # Best hyperparameter model
│   │   ├── model_fitted/    # Final trained model
│   │   └── best_params.json # Optimal hyperparameters
│   ├── pred/
│   │   └── output.csv       # Predictions with accuracy metrics
│   ├── callback_data/       # TensorBoard logs and checkpoints
│   └── predictions.log      # Detailed analysis log
└── wtn/                     # Similar structure for within-environment
```

## Data Requirements for Python Scripts

Python scripts expect the following data files in the `data/` directory:
- `pheno_final_acr.feather` / `pheno_final_wtn.feather` - Preprocessed phenotype data
- `g_data_add_acr.feather` / `g_data_add_wtn.feather` - Preprocessed genomic data
- `ev_data_wtn.feather` - Environmental data (WTN only)

*Note: Data preprocessing and format conversion between R and Python is handled by the respective setup scripts.*
