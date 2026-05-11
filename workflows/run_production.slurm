#!/bin/bash
#SBATCH --job-name=transpoly_prod
#SBATCH --partition=gpuq
#SBATCH --output=transpoly_prod_%j.out
#SBATCH --error=transpoly_prod_%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Production stage for transpoly (10 ns NVT)

set -euo pipefail

# Load environment
source /data/srinivab/miniconda/etc/profile.d/conda.sh
export OLD_GMXLDLIB="${OLD_GMXLDLIB-}"
set +u
conda activate ambertools
set -u

export GMXLIB=/data/srinivab/miniconda/envs/ambertools/share/gromacs/top
export OMP_NUM_THREADS=8

cd "$(dirname "$0")"

echo "Host: $(hostname)"
echo "Start: $(date)"

which gmx_mpi

# Run production
python3 -m transpoly.production

echo "End: $(date)"
