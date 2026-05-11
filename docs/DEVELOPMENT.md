# Development and Testing

## Running Tests

```bash
# Validate environment
bash check_env.sh

# Run minimal example
cd examples
python3 minimal_run.py
```

## Adding Custom Analysis

1. Create a new analysis class in `analysis/analysis.py`:

```python
class MyAnalysis:
    @staticmethod
    def run(output_dir, logger):
        # Your analysis code
        pass
```

2. Call from `pipeline.py::run_analysis()`:

```python
def run_analysis(self):
    # ... existing code ...
    MyAnalysis.run(self.output_dir, self.logger)
```

## Extending the Pipeline

### Add a new stage

1. Create `src/transpoly/my_stage.py`
2. Implement class inheriting from pattern in existing stages
3. Add to `pipeline.py::TranspolyPipeline.run_all()`

### Override MDP parameters

Currently, MDP generation is in `gromacs_prep.py`. To make it configurable:

1. Add fields to `SimulationConfig` in `config.py`
2. Reference in `gromacs_prep.py::generate_mdp_files()`
3. Document in example config

## Debugging

Set logging level in `utils.py::setup_logger()`:

```python
logger.setLevel(logging.DEBUG)  # More verbose
```

Check log file in output directory:
```bash
tail -f transpoly_output/pipeline.log
```

## Common Issues

### LINCS blow-up

- Increase NVT Berendsen equilibration time
- Lower temperature or reduce initial velocities
- Use more aggressive initial minimization

### Packmol fails

- Check PDB format (no QM atoms, hydrogens may need removal)
- Verify box dimensions are reasonable
- Try smaller chain count first

### AmberTools hangs

- Kill and restart: `killall antechamber`
- Check for permission issues on temp files
- Verify PDB has standard amino acids / atoms

### GPU not detected

- Check SLURM partition (use `gpuq` for GPU partitions)
- Set `export CUDA_VISIBLE_DEVICES=0` if specific GPU needed
- Verify GROMACS was compiled with GPU support

## Future Enhancements

- [ ] Template-based MDP generation (allow user templates)
- [ ] Automatic charge neutralization logic
- [ ] Support for multiple force fields (AMBER, CHARMM, etc.)
- [ ] GPU acceleration detection and setup
- [ ] Online analysis (compute properties during production)
- [ ] Restart from checkpoint (mid-run resume)
- [ ] Trajectory compression / storage optimization
- [ ] Integration with analysis frameworks (e.g., MDAnalysis pipelines)

## Code Style

Follow PEP 8. Use type hints where practical:

```python
def my_function(param: str, count: int) -> Path:
    """Do something.
    
    Args:
        param: A string parameter
        count: Integer count
    
    Returns:
        Path to output file
    """
    pass
```

Document all public functions and classes with docstrings.
