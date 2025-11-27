# KIBREED_public

[![DOI](https://zenodo.org/badge/990447884.svg)](https://doi.org/10.5281/zenodo.17737517)

A genomic prediction project for plant breeding applications using both **R statistical models** and **Python deep learning (CNN)** approaches with specific version requirements.

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

## Documentation

| Guide | Description |
|-------|-------------|
| **[GBLUP-based models](docs/r_analysis.md)** | Traditional genomic prediction using R/BGLR |
| **[CNN-based models](docs/py_analysis.md)** | Deep learning genomic prediction with TensorFlow |

## Data Requirements

Both R and Python scripts expect specific data files in the `data/` directory. See the individual analysis guides for detailed data format requirements.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for full details.

