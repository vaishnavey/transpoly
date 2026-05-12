"""
Configuration and defaults for transpoly pipeline.
"""
from dataclasses import dataclass, asdict
import yaml
from pathlib import Path
from typing import Optional


@dataclass
class SimulationConfig:
    """Core simulation parameters."""
    workflow_mode: str = "polymer"
    single_chain_pdb: Optional[str] = None
    n_chains: Optional[int] = None
    target_density: Optional[float] = None
    box_x: float = 30.0  # angstrom
    box_y: float = 30.0
    box_z: float = 90.0
    
    # Ion counts
    kcl_count: int = 5
    nh4cl_count: int = 5
    ammonium_itp: str = "ammonium.itp"
    ammonium_gro: str = "NH4.gro"

    # Replicate controls
    independent_runs: bool = False
    independent_run_count: int = 5
    
    # Temperature and pressure
    temperature: int = 300  # K
    pressure: float = 1.0  # bar
    
    # Force field
    force_field: str = "opls_aa"
    water_model: str = "tip4p"
    
    # Polymer parameterization
    charge_method: str = "bcc"  # bcc or gas
    
    # Equilibration
    equil_time_nvt: int = 500  # ps
    equil_time_npt_berendsen: int = 200  # ps
    equil_time_npt_crescale: int = 500  # ps
    
    # Production
    prod_time: int = 10000  # ps (10 ns)
    
    # SLURM
    equil_partition: str = "defq"
    prod_partition: str = "gpuq"
    conda_env: str = "ambertools"
    
    # Output
    output_dir: str = "transpoly_output"
    
    def to_dict(self) -> dict:
        return asdict(self)
    
    @classmethod
    def from_yaml(cls, path: Path) -> "SimulationConfig":
        """Load configuration from YAML file."""
        with open(path) as f:
            data = yaml.safe_load(f)
        return cls(**data)
    
    def to_yaml(self, path: Path) -> None:
        """Save configuration to YAML file."""
        with open(path, 'w') as f:
            yaml.dump(self.to_dict(), f, default_flow_style=False)
    
    def estimate_n_chains(self) -> int:
        """
        Estimate number of chains from target density.
        Assumes average chain mass and box volume.
        """
        if self.n_chains is not None:
            return self.n_chains
        if self.target_density is None:
            raise ValueError("Either n_chains or target_density must be specified")
        
        # Box volume in nm^3
        box_vol_nm3 = (self.box_x / 10.0) * (self.box_y / 10.0) * (self.box_z / 10.0)
        box_vol_cm3 = box_vol_nm3 * 1e-24  # Convert nm^3 to cm^3
        
        # Assume average polymer chain mass ~ 1200 g/mol (e.g., 10-mer at ~120 g/mol per residue)
        # Avogadro's number
        avogadro = 6.022e23
        chain_mass = 1200.0  # g/mol
        
        # mass = density * volume * n_chains / avogadro
        # n_chains = mass * avogadro / (density * volume * chain_mass)
        n_chains = int(
            (self.target_density * box_vol_cm3 * avogadro) / chain_mass
        )
        return max(1, n_chains)


# Defaults
DEFAULT_CONFIG = SimulationConfig(
    workflow_mode="polymer",
    single_chain_pdb="chain.pdb",
    n_chains=50,
    box_x=30.0,
    box_y=30.0,
    box_z=90.0,
)
