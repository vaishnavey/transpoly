# Transpoly
<img width="1536" height="1024" alt="image" src="https://github.com/user-attachments/assets/fd208124-82b1-40eb-959c-fb1edfb61b46" />

A reproducible molecular simulation pipeline for polymer/salt systems.

Transpoly automates the workflow from a single-chain protein structure to a packed, solvated GROMACS system, through multistage equilibration, production MD, and comprehensive analysis including ion transport metrics and coordination statistics.

# What you will need
- A .pdb file for a single chain

# User inputs
- Box dimensions
- Ions to consider
- Target density/ Number of chains to fill

## Quick Start

See [docs/SETUP.md](docs/SETUP.md) for environment setup and [docs/USAGE.md](docs/USAGE.md) for pipeline execution.

## Requirements

- AmberTools (antechamber, parmchk2, tleap)
- GROMACS
- Packmol
- Acpype
- Python 3.8+

## Pipeline Overview

1. **Parameterization**: Convert single-chain PDB to Amber topology using GAFF2
2. **Packing**: Fill simulation box with multiple chains using Packmol
3. **GROMACS Prep**: Convert to GROMACS with OPLS-AA force field
4. **Solvation & Ions**: Add water (TIP4P), KCl, NH4Cl, and neutralize
5. **Equilibration**: Multistage equilibration (EM → NVT → NPT)
6. **Production**: 10 ns NVT ensemble dynamics
7. **Analysis**: Ion transport, coordination numbers, RDF, diffusivity estimates

## Repository Structure

```
transpoly/
├── src/transpoly/          # Core Python modules
├── workflows/              # Stage-specific runners
├── templates/              # .mdp, tleap, packmol templates
├── analysis/               # Post-processing and plotting
├── examples/               # Minimal runnable examples
├── docs/                   # Documentation
├── requirements.txt        # Python dependencies
└── README.md
```

## License

See LICENSE file.
