#!/bin/bash

# Check R version
R_VERSION=$(R --version 2>/dev/null | head -n1 | grep -o "R version [0-9.]*" | cut -d' ' -f3)
if [[ "$R_VERSION" != "4.0.5" ]]; then
    echo "R version 4.0.5 not found (found: ${R_VERSION:-'none'})"
    echo "Please check if R 4.0.5 is added to your .bashrc"
    exit 1
fi

# Check Python version
PYTHON_VERSION=$(python --version 2>/dev/null | cut -d' ' -f2)
if [[ "$PYTHON_VERSION" != "3.8.11" ]]; then
    echo "Python version 3.8.11 not found (found: ${PYTHON_VERSION:-'none'})"
    echo "Please check if Python 3.8.11 is added to your .bashrc"
    exit 1
fi

echo "Both R 4.0.5 and Python 3.8.11 found. Stopping script."
exit 0
