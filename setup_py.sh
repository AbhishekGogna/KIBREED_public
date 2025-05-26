set -e  # Exit on first error

PROJECT_DIR="$(pwd)"
VENV_DIR="$PROJECT_DIR/py_env"
FREEZE_FILE="$PROJECT_DIR/pip_freeze.txt"
PIP="$VENV_DIR/bin/pip"

echo "📁 Working directory: $PROJECT_DIR"

# Step 1: Check if freeze file exists
if [[ ! -f "$FREEZE_FILE" ]]; then
  echo "❌ Error: pip_freeze.txt not found in $PROJECT_DIR"
  echo "Please generate it with: pip freeze > pip_freeze.txt"
  exit 1
fi

# Step 2: Create virtual environment if not present
if [[ ! -d "$VENV_DIR" ]]; then
  echo "🐍 Creating virtual environment in py_env/"
  python3 -m venv "$VENV_DIR"
else
  echo "✅ Found existing virtual environment."
fi

# Step 3: Activate and install packages
echo "📦 Installing packages from pip_freeze.txt..."
"$PIP" install --upgrade pip
"$PIP" install -r "$FREEZE_FILE"

echo "✅ Python environment setup complete!"
