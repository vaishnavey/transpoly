"""
Equilibration stage: energy minimization + multistage equilibration.
"""
import logging
from pathlib import Path
from .utils import run_command, checkpoint_file
from .config import SimulationConfig


class EquilibrationStage:
    """Handle energy minimization and multistage equilibration."""
    
    def __init__(self, config: SimulationConfig, output_dir: Path, logger: logging.Logger):
        self.config = config
        self.output_dir = output_dir
        self.logger = logger
        self.stage_dir = output_dir / "05_equilibration"
        self.stage_dir.mkdir(parents=True, exist_ok=True)
    
    def copy_prep_files(self) -> None:
        """Copy MDP and ITP files from previous stages."""
        import shutil
        
        prep_dir = self.output_dir / "03_gromacs_prep"
        if prep_dir.exists():
            for f in prep_dir.glob("*.mdp"):
                shutil.copy(f, self.stage_dir / f.name)
            for f in prep_dir.glob("*.itp"):
                shutil.copy(f, self.stage_dir / f.name)
        
        solv_dir = self.output_dir / "04_solvation"
        if solv_dir.exists():
            for f in solv_dir.glob("solvated.gro"):
                shutil.copy(f, self.stage_dir / f.name)
            for f in solv_dir.glob("topol_ions.top"):
                shutil.copy(f, self.stage_dir / f.name)
    
    def run_energy_minimization(self) -> None:
        """Run two-stage energy minimization (SD + CG)."""
        
        solvated_gro = self.stage_dir / "solvated.gro"
        topol = self.stage_dir / "topol_ions.top"
        
        if not solvated_gro.exists() or not topol.exists():
            raise FileNotFoundError("Missing solvated.gro or topol_ions.top")
        
        # SD minimization
        em_sd_gro = self.stage_dir / "em_sd.gro"
        if not checkpoint_file(em_sd_gro, "Steepest descent minimization", self.logger):
            cmd = (
                f"gmx_mpi grompp -f em_sd.mdp -c {solvated_gro.name} "
                f"-r {solvated_gro.name} -p {topol.name} -o em_sd.tpr -maxwarn 1"
            )
            run_command(cmd, self.stage_dir, self.logger, description="grompp EM-SD")
            
            cmd = "gmx_mpi mdrun -deffnm em_sd"
            run_command(cmd, self.stage_dir, self.logger, description="mdrun EM-SD")
        
        # CG minimization
        em_cg_gro = self.stage_dir / "em_cg.gro"
        if not checkpoint_file(em_cg_gro, "Conjugate gradient minimization", self.logger):
            cmd = (
                f"gmx_mpi grompp -f em_cg.mdp -c {em_sd_gro.name} "
                f"-p {topol.name} -o em_cg.tpr -maxwarn 1"
            )
            run_command(cmd, self.stage_dir, self.logger, description="grompp EM-CG")
            
            cmd = "gmx_mpi mdrun -deffnm em_cg"
            run_command(cmd, self.stage_dir, self.logger, description="mdrun EM-CG")
    
    def run_nvt_berendsen(self) -> None:
        """Run NVT equilibration with Berendsen thermostat (restrained)."""
        
        em_cg_gro = self.stage_dir / "em_cg.gro"
        topol = self.stage_dir / "topol_ions.top"
        nvt_ber_gro = self.stage_dir / "nvt_berendsen.gro"
        
        if checkpoint_file(nvt_ber_gro, "NVT Berendsen equilibration", self.logger):
            return
        
        cmd = (
            f"gmx_mpi grompp -f nvt_berendsen.mdp -c {em_cg_gro.name} "
            f"-r {em_cg_gro.name} -p {topol.name} -o nvt_berendsen.tpr -maxwarn 1"
        )
        run_command(cmd, self.stage_dir, self.logger, description="grompp NVT-Berendsen")
        
        cmd = "gmx_mpi mdrun -deffnm nvt_berendsen"
        run_command(cmd, self.stage_dir, self.logger, description="mdrun NVT-Berendsen")
    
    def run_nvt_vrescale(self) -> None:
        """Run NVT equilibration with v-rescale thermostat (restrained)."""
        
        nvt_ber_gro = self.stage_dir / "nvt_berendsen.gro"
        nvt_ber_cpt = self.stage_dir / "nvt_berendsen.cpt"
        topol = self.stage_dir / "topol_ions.top"
        nvt_vr_gro = self.stage_dir / "nvt_vrescale.gro"
        
        if checkpoint_file(nvt_vr_gro, "NVT v-rescale equilibration", self.logger):
            return
        
        cmd = (
            f"gmx_mpi grompp -f nvt_vrescale.mdp -c {nvt_ber_gro.name} "
            f"-r {nvt_ber_gro.name} -p {topol.name} -o nvt_vrescale.tpr -maxwarn 1"
        )
        if nvt_ber_cpt.exists():
            cmd += f" -t {nvt_ber_cpt.name}"
        
        run_command(cmd, self.stage_dir, self.logger, description="grompp NVT-vrescale")
        
        cmd = "gmx_mpi mdrun -deffnm nvt_vrescale"
        run_command(cmd, self.stage_dir, self.logger, description="mdrun NVT-vrescale")
    
    def run_npt_berendsen(self) -> None:
        """Run NPT pre-equilibration with Berendsen barostat."""
        
        nvt_vr_gro = self.stage_dir / "nvt_vrescale.gro"
        nvt_vr_cpt = self.stage_dir / "nvt_vrescale.cpt"
        topol = self.stage_dir / "topol_ions.top"
        npt_ber_gro = self.stage_dir / "npt_berendsen.gro"
        
        if checkpoint_file(npt_ber_gro, "NPT Berendsen equilibration", self.logger):
            return
        
        cmd = (
            f"gmx_mpi grompp -f npt_berendsen.mdp -c {nvt_vr_gro.name} "
            f"-p {topol.name} -o npt_berendsen.tpr -maxwarn 1"
        )
        if nvt_vr_cpt.exists():
            cmd += f" -t {nvt_vr_cpt.name}"
        
        run_command(cmd, self.stage_dir, self.logger, description="grompp NPT-Berendsen")
        
        cmd = "gmx_mpi mdrun -deffnm npt_berendsen"
        run_command(cmd, self.stage_dir, self.logger, description="mdrun NPT-Berendsen")
    
    def run_npt_crescale(self) -> None:
        """Run NPT main equilibration with C-rescale barostat."""
        
        npt_ber_gro = self.stage_dir / "npt_berendsen.gro"
        npt_ber_cpt = self.stage_dir / "npt_berendsen.cpt"
        topol = self.stage_dir / "topol_ions.top"
        npt_cres_gro = self.stage_dir / "npt_crescale.gro"
        
        if checkpoint_file(npt_cres_gro, "NPT C-rescale equilibration", self.logger):
            return
        
        cmd = (
            f"gmx_mpi grompp -f npt_crescale.mdp -c {npt_ber_gro.name} "
            f"-p {topol.name} -o npt_crescale.tpr -maxwarn 1"
        )
        if npt_ber_cpt.exists():
            cmd += f" -t {npt_ber_cpt.name}"
        
        run_command(cmd, self.stage_dir, self.logger, description="grompp NPT-Crescale")
        
        cmd = "gmx_mpi mdrun -deffnm npt_crescale"
        run_command(cmd, self.stage_dir, self.logger, description="mdrun NPT-Crescale")
    
    def extract_final_density(self) -> float:
        """Extract and report final density from NPT."""
        
        npt_edr = self.stage_dir / "npt_crescale.edr"
        if not npt_edr.exists():
            self.logger.warning("No NPT edr file found")
            return 0.0
        
        density_xvg = self.stage_dir / "npt_density.xvg"
        
        cmd = (
            f"printf 'Density\\n0\\n' | gmx_mpi energy -f {npt_edr.name} "
            f"-o {density_xvg.name}"
        )
        
        run_command(cmd, self.stage_dir, self.logger, check=False, description="Extract density")
        
        # Parse final density
        try:
            from .utils import read_xvg
            times, densities = read_xvg(density_xvg)
            if densities:
                final_density = densities[-1]
                self.logger.info(f"Final density: {final_density:.2f} kg/m³")
                return final_density
        except Exception as e:
            self.logger.warning(f"Could not extract density: {e}")
        
        return 0.0
    
    def run_all(self) -> Path:
        """Execute full equilibration pipeline."""
        self.logger.info("="*60)
        self.logger.info("EQUILIBRATION STAGE")
        self.logger.info("="*60)
        
        self.copy_prep_files()
        self.run_energy_minimization()
        self.run_nvt_berendsen()
        self.run_nvt_vrescale()
        self.run_npt_berendsen()
        self.run_npt_crescale()
        
        final_density = self.extract_final_density()
        
        final_gro = self.stage_dir / "npt_crescale.gro"
        self.logger.info(f"Equilibration complete: {final_gro}")
        
        return final_gro
