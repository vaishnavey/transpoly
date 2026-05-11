# Quick Start Guide

## Setup Environment

1. Install dependencies (see docs/SETUP.md)
2. Clone this repository
3. Create a configuration file (see `examples/config_example.yaml`)

## Basic Usage

### Option 1: Configuration File (Recommended)

```bash
cd transpoly/workflows
python run_pipeline.py --config ../examples/config_example.yaml
```

### Option 2: Python Script

```python
from transpoly.config import SimulationConfig
from transpoly.pipeline import TranspolyPipeline

config = SimulationConfig(
    single_chain_pdb="my_polymer.pdb",
    n_chains=50,
    box_x=30.0, box_y=30.0, box_z=90.0,
    output_dir="my_simulation"
)

pipeline = TranspolyPipeline(config)
pipeline.run_all()
```

### Option 3: SLURM Submission

```bash
cd transpoly/workflows
bash submit_pipeline.sh
```

This submits equilibration and production as separate jobs with dependencies.

## Input Files

- `my_polymer.pdb`: Single polymer chain structure (e.g., 10-mer)

## Output Structure

```
my_simulation/
├── 01_parameterization/     # AmberTools outputs
├── 02_packing/              # Packed structure
├── 03_gromacs_prep/         # MDP files and topology
├── 04_solvation/            # Solvated system
├── 05_equilibration/        # EM and equilibration trajectories
├── 06_production/           # 10 ns NVT production
└── 07_analysis/             # MSD, RDF, diffusivity plots
```

## Density-Based Chain Count

Instead of specifying `n_chains`, you can specify `target_density` and the pipeline will estimate the number of chains:

```yaml
target_density: 1.2  # g/cm³
box_x: 30.0
box_y: 30.0
box_z: 90.0
```

## Customization

Edit `src/transpoly/config.py` or pass parameters to override defaults:

- Temperature: `temperature: 310` (for 310 K)
- Ion counts: `kcl_count: 10`, `nh4cl_count: 5`
- Equilibration times (ps): adjust `equil_time_*` parameters
- Force field: currently supports OPLS-AA

## Troubleshooting

- **Missing antechamber/gmx**: Ensure conda environment is activated (see SETUP.md)
- **packmol errors**: Check PDB format; remove non-standard atoms
- **LINCS failures**: Increase NVT Berendsen time or reduce initial velocities

## Analysis Output

After production completes, check:
- `07_analysis/msd_ions.png`: Ion mean square displacement
- `07_analysis/rdf_ions.png`: Radial distribution functions
- `07_analysis/diffusivity_summary.txt`: Diffusivity estimates
- `pipeline.log`: Full execution log
