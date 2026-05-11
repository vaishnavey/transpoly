#!/bin/bash
#SBATCH --job-name=transpoly_equil
#SBATCH --partition=defq
#SBATCH --output=transpoly_equil_%j.out
#SBATCH --error=transpoly_equil_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Equilibration stage for transpoly

set -euo pipefail

# Load environment
source /data/srinivab/miniconda/etc/profile.d/conda.sh
export OLD_GMXLDLIB="${OLD_GMXLDLIB-}"
set +u
conda activate ambertools
set -u

# Set GROMACS library path
export GMXLIB=/data/srinivab/miniconda/envs/ambertools/share/gromacs/top
export OMP_NUM_THREADS=8

cd "$(dirname "$0")"

echo "Host: $(hostname)"
echo "Start: $(date)"
echo "Working directory: $(pwd)"

which gmx_mpi
which python3

# Run equilibration
python3 -m transpoly.equilibration

echo "End: $(date)"
