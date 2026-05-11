"""
Production MD stage: 10 ns NVT ensemble.
"""
import logging
from pathlib import Path
from .utils import run_command, checkpoint_file
from .config import SimulationConfig


class ProductionStage:
    """Handle 10 ns NVT production MD."""
    
    def __init__(self, config: SimulationConfig, output_dir: Path, logger: logging.Logger):
        self.config = config
        self.output_dir = output_dir
        self.logger = logger
        self.stage_dir = output_dir / "06_production"
        self.stage_dir.mkdir(parents=True, exist_ok=True)
    
    def copy_equil_files(self) -> None:
        """Copy necessary files from equilibration stage."""
        import shutil
        
        equil_dir = self.output_dir / "05_equilibration"
        if equil_dir.exists():
            # Copy MDP
            for f in equil_dir.glob("nvt_prod.mdp"):
                shutil.copy(f, self.stage_dir / f.name)
            
            # Copy final coordinates and checkpoint
            for f in equil_dir.glob("npt_crescale.gro"):
                shutil.copy(f, self.stage_dir / f.name)
            for f in equil_dir.glob("npt_crescale.cpt"):
                shutil.copy(f, self.stage_dir / f.name)
            
            # Copy topology and ITP files
            for f in equil_dir.glob("*.top"):
                shutil.copy(f, self.stage_dir / f.name)
            for f in equil_dir.glob("*.itp"):
                shutil.copy(f, self.stage_dir / f.name)
    
    def run_production(self) -> None:
        """Run 10 ns NVT production."""
        
        start_gro = self.stage_dir / "npt_crescale.gro"
        topol = self.stage_dir / "topol_ions.top"
        prod_gro = self.stage_dir / "nvt_prod.gro"
        
        if not start_gro.exists() or not topol.exists():
            raise FileNotFoundError("Missing start coordinates or topology")
        
        if checkpoint_file(prod_gro, "10 ns NVT production", self.logger):
            return
        
        start_cpt = self.stage_dir / "npt_crescale.cpt"
        
        cmd = (
            f"gmx_mpi grompp -f nvt_prod.mdp -c {start_gro.name} "
            f"-p {topol.name} -o nvt_prod.tpr -maxwarn 1"
        )
        if start_cpt.exists():
            cmd += f" -t {start_cpt.name}"
        
        run_command(cmd, self.stage_dir, self.logger, description="grompp production")
        
        cmd = "gmx_mpi mdrun -deffnm nvt_prod"
        run_command(cmd, self.stage_dir, self.logger, description="mdrun production (10 ns)")
    
    def run_all(self) -> Path:
        """Execute production MD pipeline."""
        self.logger.info("="*60)
        self.logger.info("PRODUCTION MD STAGE (10 ns NVT)")
        self.logger.info("="*60)
        
        self.copy_equil_files()
        self.run_production()
        
        prod_xtc = self.stage_dir / "nvt_prod.xtc"
        prod_gro = self.stage_dir / "nvt_prod.gro"
        
        self.logger.info(f"Production complete: {prod_gro}")
        return prod_xtc
