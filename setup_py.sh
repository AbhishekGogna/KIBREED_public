#!/bin/bash
set -e  # Exit on first error

PROJECT_DIR="$(pwd)"
VENV_DIR="$PROJECT_DIR/py_env"
FREEZE_FILE="$PROJECT_DIR/pyenv.lock"
PIP="$VENV_DIR/bin/pip"

echo "📁 Working directory: $PROJECT_DIR"

# Step 0: Check Python version
PYTHON_VERSION_FULL=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
PYTHON_MAJOR=$(python3 -c 'import sys; print(sys.version_info[0])')
PYTHON_MINOR=$(python3 -c 'import sys; print(sys.version_info[1])')

if (( PYTHON_MAJOR < 3 )) || (( PYTHON_MAJOR == 3 && PYTHON_MINOR < 8 )); then
  echo "❌ Python 3.8 or higher is required. Found Python $PYTHON_VERSION_FULL"
  exit 1
elif (( PYTHON_MAJOR > 3 )) || (( PYTHON_MAJOR == 3 && PYTHON_MINOR > 8 )); then
  echo "⚠️  Warning: Detected Python $PYTHON_VERSION_FULL. This project was tested with Python 3.8.11"
fi

# Step 1: Check if freeze file exists
if [[ ! -f "$FREEZE_FILE" ]]; then
  echo "❌ Error: pip_freeze.txt not found in $PROJECT_DIR"
  exit 1
fi

# Step 2: Create virtual environment if not present
if [[ ! -d "$VENV_DIR" ]]; then
  echo "Creating virtual environment in py_env/"
  python3 -m venv "$VENV_DIR"
else
  echo "✅ Found existing virtual environment."
fi

# Step 3: Activate and install packages
echo "Installing packages from pyenv.lock..."
"$PIP" install --upgrade pip
"$PIP" install -r "$FREEZE_FILE"

echo "✅ Python environment setup complete!"
