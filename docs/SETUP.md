# Environment Setup

## Required Tools

### AmberTools
Install AmberTools (includes antechamber, parmchk2, tleap):
```bash
conda install -c conda-forge ambertools
```

### GROMACS
Install GROMACS (with MPI support recommended):
```bash
conda install -c conda-forge gromacs
```

### Packmol
Install Packmol:
```bash
conda install -c conda-forge packmol
```

### Acpype
Install Acpype for Amber→GROMACS conversion:
```bash
pip install acpype
```

## Python Environment

Create a virtual environment and install dependencies:
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Verify Installation

```bash
which antechamber
which gmx
which packmol
python -c "import MDAnalysis; print(MDAnalysis.__version__)"
```

All commands should return valid paths/versions.

## Force Field Files

OPLS-AA parameters are typically included with GROMACS. Set the environment variable:
```bash
export GMXLIB=${GROMACS_INSTALL_DIR}/share/gromacs/top
```

TIP4P water model files will be referenced during solvation.
