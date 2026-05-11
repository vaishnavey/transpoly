"""
Parameterization using AmberTools: antechamber, parmchk2, tleap, acpype.
"""
import logging
from pathlib import Path
from .utils import run_command, checkpoint_file, write_file, get_topology_resname
from .config import SimulationConfig


class ParameterizeStage:
    """Handle AmberTools parameterization workflow."""
    
    def __init__(self, config: SimulationConfig, output_dir: Path, logger: logging.Logger):
        self.config = config
        self.output_dir = output_dir
        self.logger = logger
        self.stage_dir = output_dir / "01_parameterization"
        self.stage_dir.mkdir(parents=True, exist_ok=True)
    
    def get_base_name(self) -> str:
        """Get base name for parameterized files (e.g., 'polymer' from 'polymer.pdb')."""
        return Path(self.config.single_chain_pdb).stem
    
    def run_antechamber(self) -> None:
        """Run antechamber to convert PDB to GAFF2 .mol2."""
        base = self.get_base_name()
        pdb_src = Path(self.config.single_chain_pdb)
        
        if not pdb_src.exists():
            raise FileNotFoundError(f"Input PDB not found: {pdb_src}")
        
        # Copy PDB to stage directory
        pdb_file = self.stage_dir / pdb_src.name
        if not pdb_file.exists():
            import shutil
            shutil.copy(pdb_src, pdb_file)
        
        mol2_file = self.stage_dir / f"{base}_gaff2.mol2"
        frcmod_file = self.stage_dir / f"{base}_gaff2.frcmod"
        
        if checkpoint_file(mol2_file, "GAFF2 parameterization", self.logger):
            return
        
        # Try with BCC first, then fall back to gas-phase
        charge_method = self.config.charge_method
        
        for attempt, method in enumerate([charge_method, "gas"], start=1):
            cmd = (
                f"antechamber -i {pdb_file} -fi pdb "
                f"-o {mol2_file} -fo mol2 "
                f"-at gaff2 -c {method} -s 2"
            )
            rc, stdout, stderr = run_command(
                cmd,
                self.stage_dir,
                self.logger,
                check=False,
                description=f"Antechamber (attempt {attempt}, {method})"
            )
            
            if rc == 0 and mol2_file.exists():
                self.logger.info(f"Antechamber succeeded with {method} charge method")
                break
            elif attempt == 1:
                self.logger.warning(f"Antechamber failed with {method}; retrying with gas")
            else:
                raise RuntimeError(f"Antechamber failed: {stderr}")
        
        # Run parmchk2
        if not frcmod_file.exists():
            cmd = f"parmchk2 -i {mol2_file} -f mol2 -o {frcmod_file}"
            run_command(
                cmd,
                self.stage_dir,
                self.logger,
                description="parmchk2 (generate missing parameters)"
            )
    
    def run_tleap(self) -> None:
        """Use tleap to build Amber topology and coordinates."""
        base = self.get_base_name()
        mol2_file = self.stage_dir / f"{base}_gaff2.mol2"
        frcmod_file = self.stage_dir / f"{base}_gaff2.frcmod"
        prmtop_file = self.stage_dir / f"{base}.prmtop"
        inpcrd_file = self.stage_dir / f"{base}.inpcrd"
        
        if checkpoint_file(prmtop_file, "tleap topology building", self.logger):
            return
        
        # Create tleap input script
        tleap_script = self.stage_dir / "tleap.in"
        write_file(
            tleap_script,
            f"""source leaprc.gaff2
mol = loadmol2 {mol2_file.name}
loadamberparams {frcmod_file.name}
saveamberparm mol {prmtop_file.name} {inpcrd_file.name}
savepdb mol {base}_tleap.pdb
quit
"""
        )
        
        cmd = f"tleap -f {tleap_script.name}"
        run_command(
            cmd,
            self.stage_dir,
            self.logger,
            description="tleap (Amber topology)"
        )
    
    def run_acpype(self) -> None:
        """Convert Amber topology to GROMACS format using acpype."""
        base = self.get_base_name()
        prmtop_file = self.stage_dir / f"{base}.prmtop"
        inpcrd_file = self.stage_dir / f"{base}.inpcrd"
        
        # Check if acpype output already exists
        for ac_dir in [f"{base}.acpype", f"{base}.amb2gmx"]:
            if (self.stage_dir / ac_dir).exists():
                self.logger.info(f"Skipping acpype: found {ac_dir}")
                return
        
        self.logger.info("Running acpype (Amber→GROMACS)")
        
        cmd = f"acpype -p {prmtop_file.name} -x {inpcrd_file.name} -b {base}"
        rc, stdout, stderr = run_command(
            cmd,
            self.stage_dir,
            self.logger,
            check=False,
            description="acpype (Amber→GROMACS conversion)"
        )
        
        if rc != 0:
            self.logger.warning("acpype with topology failed; trying with mol2")
            mol2_file = self.stage_dir / f"{base}_gaff2.mol2"
            cmd = f"acpype -i {mol2_file.name} -b {base}"
            run_command(
                cmd,
                self.stage_dir,
                self.logger,
                description="acpype (mol2 input)"
            )
    
    def extract_gmx_files(self) -> tuple[Path, Path, Path]:
        """
        Extract .gro, .itp, and .top from acpype output.
        Returns (gro_file, itp_file, top_file)
        """
        base = self.get_base_name()
        
        # Find acpype output directory
        ac_dir = None
        for d in [f"{base}.acpype", f"{base}.amb2gmx"]:
            potential_dir = self.stage_dir / d
            if potential_dir.exists():
                ac_dir = potential_dir
                break
        
        if not ac_dir:
            raise FileNotFoundError(f"No acpype output directory found for {base}")
        
        gro_file = ac_dir / f"{base}_GMX.gro"
        itp_file = ac_dir / f"{base}_GMX.itp"
        top_file = ac_dir / f"{base}_GMX.top"
        
        if not gro_file.exists():
            raise FileNotFoundError(f"Missing coordinate file: {gro_file}")
        
        return gro_file, itp_file, top_file
    
    def run_all(self) -> None:
        """Execute full parameterization pipeline."""
        self.logger.info("="*60)
        self.logger.info("PARAMETERIZATION STAGE")
        self.logger.info("="*60)
        
        self.run_antechamber()
        self.run_tleap()
        self.run_acpype()
        
        gro_file, itp_file, top_file = self.extract_gmx_files()
        self.logger.info(f"Parameterization complete:")
        self.logger.info(f"  Coordinates: {gro_file}")
        self.logger.info(f"  Topology ITP: {itp_file}")
        self.logger.info(f"  Topology TOP: {top_file}")
