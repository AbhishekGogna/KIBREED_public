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

3. Set up R environment
This script restores the R package environment using renv and the renv.lock file.

```bash
Rscript setup.R
```

Notes
1. Run all commands from the project root directory.
2. If you get version mismatch errors, ensure the correct versions are installed and available in your $PATH.
