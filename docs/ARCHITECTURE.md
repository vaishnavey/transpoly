# Pipeline Architecture

## Workflow Stages

### 1. Parameterization (01_parameterization/)

Uses AmberTools to convert PDB to GROMACS-compatible format:
- `antechamber`: PDB → GAFF2 force field assignment
- `parmchk2`: Validate/complete missing parameters
- `tleap`: Build Amber topology and coordinates
- `acpype`: Convert to GROMACS ITP/TOP/GRO

**Inputs**: `single_chain.pdb`
**Outputs**: `polymer_gaff2.mol2`, `polymer.prmtop`, `polymer_GMX.itp`

### 2. Packing (02_packing/)

Pack multiple chains into simulation box:
- `packmol`: Fill box with requested number of chains
- Support density-based estimation: `n_chains = estimate_n_chains(target_density)`

**Inputs**: Single chain PDB, box dimensions, chain count
**Outputs**: `packed_box.pdb`

### 3. GROMACS Preparation (03_gromacs_prep/)

Generate all MDP files and set up topology:
- Energy minimization: SD (restrained) + CG (unrestrained)
- NVT equilibration: Berendsen → v-rescale
- NPT equilibration: Berendsen → C-rescale
- Production: 10 ns NVT

**Outputs**: `em_sd.mdp`, `em_cg.mdp`, `nvt_berendsen.mdp`, `nvt_vrescale.mdp`, `npt_berendsen.mdp`, `npt_crescale.mdp`, `nvt_prod.mdp`, `POLY.itp`

### 4. Solvation & Ionization (04_solvation/)

Prepare solvated, ionized system:
- Convert packed PDB to GRO with box dimensions
- Insert K⁺ and Cl⁻ ions (5 KCl + 5 NH₄Cl by default)
- Solvate with TIP4P water
- Update topology: `topol_ions.top`

**Outputs**: `solvated.gro`, `topol_ions.top`

### 5. Equilibration (05_equilibration/)

Multi-stage equilibration to avoid LINCS failures:
1. Energy minimization (SD + CG)
2. NVT Berendsen (restrained, 500 ps)
3. NVT v-rescale (restrained, 500 ps)
4. NPT Berendsen (pre-equilibration, 200 ps)
5. NPT C-rescale (main equilibration, 500 ps)

Extract final density from NPT trajectory.

**Outputs**: `npt_crescale.gro`, `npt_crescale.edr`, `npt_density.xvg`

### 6. Production (06_production/)

10 ns NVT ensemble production run:
- Uses final coordinates and velocities from equilibration
- Generates trajectory (XTC/TRR format)

**Outputs**: `nvt_prod.xtc`, `nvt_prod.gro`, `nvt_prod.edr`

### 7. Analysis (07_analysis/)

Post-production analysis and visualization:
- **MSD Analysis**: Mean square displacement for each ion
  - Linear fit to estimate diffusivity
  - Plot: `msd_ions.png`
  
- **RDF Analysis**: Radial distribution functions
  - Requires `gmx rdf` pre-processing
  - Plot: `rdf_ions.png`
  
- **Coordination**: Coordination numbers over time (framework)
  
- **Energy Tracking**: LJ-SR and Coul-SR terms
  
- **Diffusivity**: 
  - From MSD linear fit: D = slope / 6
  - Optional Stokes-Einstein refinement
  - Summary: `diffusivity_summary.txt`

## Module Design

```
src/transpoly/
├── config.py              # Configuration and defaults
├── utils.py               # Logging, file I/O, checkpointing
├── parameterize.py        # AmberTools orchestration
├── packing.py             # Packmol + density estimation
├── gromacs_prep.py        # MDP generation
├── solvation.py           # Water + ion addition
├── equilibration.py       # Multi-stage equilibration
├── production.py          # Production MD
├── pipeline.py            # Main orchestrator
└── main.py                # Entry point

analysis/
├── analysis.py            # MSD, RDF, diffusivity, energy tracking
└── __init__.py

workflows/
├── run_pipeline.py        # Main script
├── run_equilibration.slurm
├── run_production.slurm
└── submit_pipeline.sh     # SLURM submission with dependencies
```

## Key Features

### Checkpointing

Each stage saves checkpoint files. If a stage is interrupted, re-running automatically skips completed steps.

### SLURM Integration

- Separate SLURM scripts for equilibration and production
- Automatic job dependency submission: production waits for equilibration
- Configurable partitions and resources

### Reproducibility

- All configurations saved to YAML
- Unique random seeds for production (deterministic but independent)
- Detailed logging to `pipeline.log`

### Density Estimation

If `target_density` is provided instead of `n_chains`:
```
n_chains = (target_density * box_volume * avogadro) / chain_mass
```

Assumes average chain mass ~ 1200 g/mol (customizable).

## Configuration Example

```yaml
single_chain_pdb: polymer_10mer.pdb
n_chains: 50              # or target_density: 1.2
box_x: 30.0
box_y: 30.0
box_z: 90.0
kcl_count: 5
nh4cl_count: 5
temperature: 300
equil_time_nvt: 500       # ps
equil_time_npt_berendsen: 200
equil_time_npt_crescale: 500
prod_time: 10000          # 10 ns
output_dir: my_simulation
```

## Extensibility

To add custom analysis:
1. Extend `analysis.py` with new analysis class
2. Call from `pipeline.run_analysis()`
3. Follow existing pattern: extract data → fit → plot → save

To modify MDP parameters:
1. Edit `gromacs_prep.py::generate_mdp_files()`
2. Or override via configuration + template system (future enhancement)
