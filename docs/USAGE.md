# Usage Guide

## Input Requirements

1. **single_chain.pdb**: A single polymer chain structure (e.g., a 10-mer protein or polyelectrolyte)
2. **n_chains** or **target_density**: Number of chains to pack, or target density (g/cm³) for automatic estimation
3. **box_dimensions**: x, y, z dimensions in angstroms

## Example Workflow

```python
from transpoly import Pipeline

pipeline = Pipeline(
    single_chain_pdb='polymer_10mer.pdb',
    n_chains=50,  # or density=1.2 for density-based estimation
    box_dims=(30.0, 30.0, 90.0),
    temperature=300,  # K, optional
    output_dir='my_simulation'
)

# Run full pipeline
pipeline.run_all()

# Or run stages individually
pipeline.parameterize()
pipeline.pack()
pipeline.prepare_gromacs()
pipeline.solvate_and_ionize()
pipeline.equilibrate()
pipeline.production()
pipeline.analyze()
```

## Output Structure

```
my_simulation/
├── 01_parameterization/
│   ├── polymer.mol2
│   ├── polymer.frcmod
│   └── polymer_tleap.pdb
├── 02_packing/
│   └── packed_box.pdb
├── 03_gromacs_prep/
│   ├── topol.top
│   └── system.gro
├── 04_solvation/
│   ├── solvated.gro
│   └── ions_added.top
├── 05_equilibration/
│   ├── em.trr
│   ├── nvt_eq.trr
│   └── npt_eq.trr
├── 06_production/
│   ├── prod_10ns.trr
│   └── prod_10ns.xtc
└── 07_analysis/
    ├── msd_ions.png
    ├── rdf_ions.png
    ├── coordination.csv
    └── summary.txt
```

## Customization

Edit temperature, ion counts, box size, and force field parameters in the config file or pass them as arguments. See examples/ for minimal working setups.
