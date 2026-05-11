"""
Solvation and ionization stage: water + counter-ions.
"""
import logging
from pathlib import Path
from .utils import run_command, checkpoint_file, write_file
from .config import SimulationConfig


class SolvationStage:
    """Handle solvation and ion addition."""
    
    def __init__(self, config: SimulationConfig, output_dir: Path, logger: logging.Logger):
        self.config = config
        self.output_dir = output_dir
        self.logger = logger
        self.stage_dir = output_dir / "04_solvation"
        self.stage_dir.mkdir(parents=True, exist_ok=True)
    
    def prepare_system_box(self, packed_pdb: Path, resname: str) -> Path:
        """Convert PDB to GRO with box dimensions."""
        system_gro = self.stage_dir / "system.gro"
        
        if checkpoint_file(system_gro, "Convert packed PDB to GRO", self.logger):
            return system_gro
        
        # Copy packed PDB to stage
        if not (self.stage_dir / packed_pdb.name).exists():
            import shutil
            shutil.copy(packed_pdb, self.stage_dir / packed_pdb.name)
        
        # editconf to set box dimensions
        bx_nm = self.config.box_x / 10.0
        by_nm = self.config.box_y / 10.0
        bz_nm = self.config.box_z / 10.0
        
        cmd = (
            f"gmx_mpi editconf -f {packed_pdb.name} -o {system_gro.name} "
            f"-box {bx_nm:.3f} {by_nm:.3f} {bz_nm:.3f} -bt triclinic"
        )
        
        run_command(
            cmd,
            self.stage_dir,
            self.logger,
            description="editconf (set box dimensions)"
        )
        
        return system_gro
    
    def add_ions(
        self,
        system_gro: Path,
        resname: str,
        topol_top: Path
    ) -> tuple[Path, int, int]:
        """
        Add K+ and Cl- ions using gmx insert-molecules.
        Returns (system_with_ions, n_k, n_cl)
        """
        
        # Calculate ion counts
        # First, estimate net charge by trying a grompp pass
        # (This is simplified; in production, parse actual charge)
        
        n_kcl = self.config.kcl_count
        n_nh4cl = self.config.nh4cl_count
        
        # For now, use direct counts. Can expand to neutralization logic.
        n_k = n_kcl + n_nh4cl  # K+ from both sources
        n_cl = n_kcl + 2 * n_nh4cl  # Cl- to neutralize
        
        self.logger.info(f"Adding ions: K+ ({n_k}), Cl- ({n_cl})")
        
        # Create single ion .gro files
        self._create_ion_gro(self.stage_dir / "K_single.gro", "K", "K")
        self._create_ion_gro(self.stage_dir / "CL_single.gro", "CL", "CL")
        
        # Add K+ ions
        system_k_gro = self.stage_dir / "system_k.gro"
        if not checkpoint_file(system_k_gro, "Add K+ ions", self.logger):
            cmd = (
                f"gmx_mpi insert-molecules -f {system_gro.name} "
                f"-ci K_single.gro -nmol {n_k} -radius 0.08 -o {system_k_gro.name}"
            )
            run_command(cmd, self.stage_dir, self.logger, description="Insert K+ ions")
        
        # Add Cl- ions
        system_kcl_gro = self.stage_dir / "system_kcl.gro"
        if not checkpoint_file(system_kcl_gro, "Add Cl- ions", self.logger):
            cmd = (
                f"gmx_mpi insert-molecules -f {system_k_gro.name} "
                f"-ci CL_single.gro -nmol {n_cl} -radius 0.08 -o {system_kcl_gro.name}"
            )
            run_command(cmd, self.stage_dir, self.logger, description="Insert Cl- ions")
        
        return system_kcl_gro, n_k, n_cl
    
    def _create_ion_gro(self, path: Path, resname: str, atom_name: str) -> None:
        """Create a minimal .gro file for a single ion."""
        content = f"""{atom_name} ion
     1
     1{resname:>5s}{atom_name:>5s}    1    0.000    0.000    0.000
   1.0   1.0   1.0
"""
        write_file(path, content)
    
    def create_ionized_topology(
        self,
        topol_base: Path,
        resname: str,
        n_polymer: int,
        n_k: int,
        n_cl: int
    ) -> Path:
        """Create topology file for ionized system."""
        topol_ions = self.stage_dir / "topol_ions.top"
        
        content = f"""#include "oplsaa.ff/forcefield.itp"
#include "POLY.itp"
#include "oplsaa.ff/ions.itp"

[ system ]
{resname} packed box with ions

[ molecules ]
{resname}   {n_polymer}
K     {n_k}
CL    {n_cl}
"""
        
        write_file(topol_ions, content)
        self.logger.info(f"Created ionized topology: {topol_ions}")
        return topol_ions
    
    def solvate(self, system_kcl_gro: Path, topol_ions: Path) -> Path:
        """Solvate with TIP4P water."""
        
        solvated_gro = self.stage_dir / "solvated.gro"
        
        if checkpoint_file(solvated_gro, "Solvate with TIP4P", self.logger):
            return solvated_gro
        
        # Use genbox (or gmx solvate in newer versions)
        cmd = (
            f"gmx_mpi solvate -cp {system_kcl_gro.name} -cs tip4p.gro "
            f"-o {solvated_gro.name} -p {topol_ions.name}"
        )
        
        rc, _, _ = run_command(
            cmd,
            self.stage_dir,
            self.logger,
            check=False,
            description="Solvate with TIP4P water"
        )
        
        if rc != 0:
            self.logger.warning("solvate failed; trying genbox")
            cmd = (
                f"genbox -cp {system_kcl_gro.name} -cs tip4p.gro "
                f"-o {solvated_gro.name} -p {topol_ions.name}"
            )
            run_command(cmd, self.stage_dir, self.logger, description="genbox solvation")
        
        return solvated_gro
    
    def run_all(
        self,
        packed_pdb: Path,
        acpype_dir: Path,
        resname: str,
        n_polymer: int
    ) -> tuple[Path, Path]:
        """Execute solvation pipeline."""
        self.logger.info("="*60)
        self.logger.info("SOLVATION & IONIZATION STAGE")
        self.logger.info("="*60)
        
        # Copy POLY.itp from gromacs_prep stage
        gromacs_prep_dir = self.output_dir / "03_gromacs_prep"
        for itp_file in [gromacs_prep_dir / "POLY.itp", gromacs_prep_dir / "POLY_posre.itp"]:
            if itp_file.exists():
                import shutil
                shutil.copy(itp_file, self.stage_dir / itp_file.name)
        
        # Copy acpype-generated ITP file
        import shutil
        acpype_itp = acpype_dir / f"{resname}_GMX.itp"
        if acpype_itp.exists():
            shutil.copy(acpype_itp, self.stage_dir / acpype_itp.name)
        
        # Prepare system
        system_gro = self.prepare_system_box(packed_pdb, resname)
        
        # Add ions
        system_kcl_gro, n_k, n_cl = self.add_ions(system_gro, resname, None)
        
        # Create topology
        topol_ions = self.create_ionized_topology(
            None, resname, n_polymer, n_k, n_cl
        )
        
        # Solvate
        solvated_gro = self.solvate(system_kcl_gro, topol_ions)
        
        self.logger.info("Solvation complete")
        return solvated_gro, topol_ions
