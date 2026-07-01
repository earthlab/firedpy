#!/bin/bash
# Usage: source cyverse_setup.sh
# Run this at the start of every CyVerse session.

FIREDPY_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Source conda for use in this script
source /opt/conda/etc/profile.d/conda.sh

# --- Environment Setup (skipped if env already exists) ---
if conda env list | grep -q "^fired "; then
    echo "'fired' environment found, skipping setup."
else
    echo "=== First-time setup ==="

    echo "Creating conda environment with Python 3.12..."
    conda create -n fired python=3.12 -c conda-forge -y

    conda activate fired

    echo "Installing GDAL..."
    conda install gdal libgdal-hdf4 -c conda-forge -y

    echo "Installing firedpy..."
    pip install -e "$FIREDPY_DIR"

    echo "Setup complete."
fi

# --- Activate environment ---
conda activate fired

# --- Earthdata Authentication ---
echo ""
echo "Authenticating with NASA Earthdata (credentials stored in ~/.netrc for this session)..."
python -c "import earthaccess; earthaccess.login(persist=True)"

echo ""
echo "Ready to run firedpy."
echo "Example: firedpy -p ./runs -c finland -y1 2023 -y2 2023 -sp 5 -tp 11 -d"