#!/bin/bash
# Quick environment validation script

set -euo pipefail

echo "=== Transpoly Environment Validation ==="
echo ""

# Check required commands
commands=("antechamber" "parmchk2" "tleap" "gmx_mpi" "packmol" "python3")
missing=()

for cmd in "${commands[@]}"; do
    if command -v "$cmd" &>/dev/null; then
        echo "✓ $cmd found: $(which $cmd)"
    else
        echo "✗ $cmd NOT found"
        missing+=("$cmd")
    fi
done

echo ""

# Check Python packages
python3 -c "import numpy; print('✓ numpy')" 2>/dev/null || echo "✗ numpy"
python3 -c "import pandas; print('✓ pandas')" 2>/dev/null || echo "✗ pandas"
python3 -c "import matplotlib; print('✓ matplotlib')" 2>/dev/null || echo "✗ matplotlib"
python3 -c "import scipy; print('✓ scipy')" 2>/dev/null || echo "✗ scipy"
python3 -c "import MDAnalysis; print('✓ MDAnalysis')" 2>/dev/null || echo "✗ MDAnalysis"

echo ""

if [[ ${#missing[@]} -gt 0 ]]; then
    echo "Missing commands: ${missing[@]}"
    echo "Please install AmberTools, GROMACS, and Packmol (see docs/SETUP.md)"
    exit 1
else
    echo "✓ All required tools found!"
fi

# Check GMXLIB
if [[ -z "${GMXLIB-}" ]]; then
    echo "⚠ GMXLIB not set. Set it manually or it will be set automatically."
else
    echo "✓ GMXLIB = $GMXLIB"
fi

echo ""
echo "=== Environment validation complete ==="
