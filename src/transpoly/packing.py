"""
Packing stage: packmol integration and density-based chain estimation.
"""
import logging
from pathlib import Path
from .utils import run_command, checkpoint_file, write_file
from .config import SimulationConfig


class PackingStage:
    """Handle polymer chain packing with packmol."""
    
    def __init__(self, config: SimulationConfig, output_dir: Path, logger: logging.Logger):
        self.config = config
        self.output_dir = output_dir
        self.logger = logger
        self.stage_dir = output_dir / "02_packing"
        self.stage_dir.mkdir(parents=True, exist_ok=True)
    
    def get_n_chains(self) -> int:
        """Determine number of chains: use n_chains if given, else estimate from density."""
        if self.config.n_chains is not None:
            self.logger.info(f"Using specified chain count: {self.config.n_chains}")
            return self.config.n_chains
        
        n_chains = self.config.estimate_n_chains()
        self.logger.info(
            f"Estimated chain count from density {self.config.target_density} g/cm³: {n_chains}"
        )
        return n_chains
    
    def generate_packmol_input(self, pdb_file: Path, n_chains: int) -> Path:
        """Generate packmol .inp file."""
        packmol_inp = self.stage_dir / "packmol_box.inp"
        
        bx, by, bz = self.config.box_x, self.config.box_y, self.config.box_z
        
        content = f"""tolerance 2.0
filetype pdb
output packed_box.pdb

structure {pdb_file.name}
  number {n_chains}
  inside box 0.0 0.0 0.0 {bx:.1f} {by:.1f} {bz:.1f}
end structure
"""
        
        write_file(packmol_inp, content)
        self.logger.info(f"Generated packmol input: {packmol_inp}")
        return packmol_inp
    
    def run_packmol(self, pdb_file: Path, n_chains: int) -> None:
        """Run packmol to pack chains into box."""
        packed_pdb = self.stage_dir / "packed_box.pdb"
        
        if checkpoint_file(packed_pdb, f"Packmol packing ({n_chains} chains)", self.logger):
            return
        
        # Copy input PDB to stage directory if not already there
        pdb_src = pdb_file
        if not (self.stage_dir / pdb_src.name).exists():
            import shutil
            shutil.copy(pdb_src, self.stage_dir / pdb_src.name)
        
        # Generate packmol input
        self.generate_packmol_input(self.stage_dir / pdb_src.name, n_chains)
        
        # Run packmol
        cmd = f"packmol < packmol_box.inp > packmol_box.log 2>&1"
        run_command(
            cmd,
            self.stage_dir,
            self.logger,
            description=f"Packmol (packing {n_chains} chains)"
        )
        
        if not packed_pdb.exists():
            raise RuntimeError("Packmol failed to generate packed_box.pdb")
    
    def run_all(self, pdb_file: Path) -> Path:
        """Execute packing pipeline."""
        self.logger.info("="*60)
        self.logger.info("PACKING STAGE")
        self.logger.info("="*60)
        
        n_chains = self.get_n_chains()
        self.run_packmol(pdb_file, n_chains)
        
        packed_pdb = self.stage_dir / "packed_box.pdb"
        self.logger.info(f"Packing complete: {packed_pdb}")
        return packed_pdb
