#!/bin/bash
# Wrapper script to run transpoly pipeline with conda environment

set -euo pipefail

# Configuration
CONFIG_FILE="${1:-examples/config_pah_example.yaml}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Activate conda environment
echo "Activating conda environment..."
source /data/srinivab/miniconda/etc/profile.d/conda.sh
conda activate ambertools

# Set GROMACS library path
export GMXLIB=/data/srinivab/miniconda/envs/ambertools/share/gromacs/top
export OMP_NUM_THREADS=8

# Verify tools
echo "Checking tools..."
which antechamber >/dev/null || { echo "ERROR: antechamber not found"; exit 1; }
which gmx_mpi >/dev/null || { echo "ERROR: gmx_mpi not found"; exit 1; }
which packmol >/dev/null || { echo "ERROR: packmol not found"; exit 1; }

echo "✓ All tools available"
echo ""

# Run pipeline
cd "$REPO_ROOT"
echo "Running transpoly pipeline with config: $CONFIG_FILE"
python3 workflows/run_pipeline.py --config "$CONFIG_FILE"
