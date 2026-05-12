"""
GROMACS preparation: MDP generation, topology setup, conversions.
"""
import logging
import random
from pathlib import Path
from .utils import run_command, write_file, timestamp, checkpoint_file
from .config import SimulationConfig


class GromacsPrepStage:
    """Handle GROMACS file preparation and MDP generation."""
    
    def __init__(self, config: SimulationConfig, output_dir: Path, logger: logging.Logger):
        self.config = config
        self.output_dir = output_dir
        self.logger = logger
        self.stage_dir = output_dir / "03_gromacs_prep"
        self.stage_dir.mkdir(parents=True, exist_ok=True)
    
    def generate_mdp_files(self, resname: str) -> None:
        """Generate all required .mdp files."""
        if checkpoint_file(self.stage_dir / "em_sd.mdp", "Generate MDP files", self.logger):
            return
        
        T = self.config.temperature
        
        # Energy minimization - steepest descent (restrained)
        write_file(
            self.stage_dir / "em_sd.mdp",
            f"""title                   = SD minimization (restrained)
define                  = -DPOSRES
integrator              = steep
nsteps                  = 50000
emtol                   = 1000.0
emstep                  = 0.01

nstenergy               = 100
nstlog                  = 100

cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz
constraints             = none
"""
        )
        
        # Energy minimization - conjugate gradient (unrestrained)
        write_file(
            self.stage_dir / "em_cg.mdp",
            f"""title                   = CG minimization (unrestrained)
integrator              = cg
nsteps                  = 50000
emtol                   = 500.0
emstep                  = 0.01

nstenergy               = 100
nstlog                  = 100

cutoff-scheme           = Verlet
ns_type                 = grid
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz
constraints             = none
"""
        )
        
        # NVT equilibration (Berendsen, then v-rescale in actual equilibration)
        nvt_steps = int((self.config.equil_time_nvt * 1000) / 1)  # fs to steps at dt=1fs
        write_file(
            self.stage_dir / "nvt_berendsen.mdp",
            f"""title                   = NVT equilibration (Berendsen, restrained)
define                  = -DPOSRES
integrator              = md
dt                      = 0.001
nsteps                  = {nvt_steps}

nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 1000
nstenergy               = 1000
nstlog                  = 1000

continuation            = no
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 2
lincs_order             = 4

cutoff-scheme           = Verlet
nstlist                 = 20
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz

tcoupl                  = Berendsen
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = {T}

pcoupl                  = no

gen_vel                 = yes
gen_temp                = {T}
gen_seed                = -1
"""
        )
        
        # NVT equilibration (v-rescale, restrained)
        write_file(
            self.stage_dir / "nvt_vrescale.mdp",
            f"""title                   = NVT equilibration (v-rescale, restrained)
define                  = -DPOSRES
integrator              = md
dt                      = 0.001
nsteps                  = {nvt_steps}

nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 1000
nstenergy               = 1000
nstlog                  = 1000

continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 2
lincs_order             = 4

cutoff-scheme           = Verlet
nstlist                 = 20
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz

tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = {T}

pcoupl                  = no

gen_vel                 = no
"""
        )
        
        # NPT equilibration - Berendsen (pre-stage)
        npt_ber_steps = int((self.config.equil_time_npt_berendsen * 1000) / 1)
        write_file(
            self.stage_dir / "npt_berendsen.mdp",
            f"""title                   = NPT equilibration (Berendsen, pre-stage)
integrator              = md
dt                      = 0.001
nsteps                  = {npt_ber_steps}

nstxout-compressed      = 1000
nstenergy               = 500
nstlog                  = 500

continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 2
lincs_order             = 4

cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz

tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = {T}

pcoupl                  = Berendsen
pcoupltype              = isotropic
tau_p                   = 5.0
ref_p                   = {self.config.pressure}
compressibility         = 4.5e-5

DispCorr                = EnerPres
gen_vel                 = no
"""
        )
        
        # NPT equilibration - C-rescale (main stage)
        npt_crescale_steps = int((self.config.equil_time_npt_crescale * 1000) / 1)
        write_file(
            self.stage_dir / "npt_crescale.mdp",
            f"""title                   = NPT equilibration (C-rescale, main)
integrator              = md
dt                      = 0.001
nsteps                  = {npt_crescale_steps}

nstxout-compressed      = 1000
nstenergy               = 500
nstlog                  = 500

continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 2
lincs_order             = 4

cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz

tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 0.1
ref_t                   = {T}

pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 5.0
ref_p                   = {self.config.pressure}
compressibility         = 4.5e-5

DispCorr                = EnerPres
gen_vel                 = no
"""
        )
        
        # NVT production
        prod_steps = int((self.config.prod_time * 1000) / 2)  # dt=2fs
        random_seed = random.randint(1, 2**31 - 1)
        energy_groups = f"K CL SOL {resname}"
        if self.config.nh4cl_count > 0:
            energy_groups = f"{energy_groups} NH4P"
        write_file(
            self.stage_dir / "nvt_prod.mdp",
            f"""title                   = NVT production ({self.config.prod_time} ps)
integrator              = md
dt                      = 0.002
nsteps                  = {prod_steps}

nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 5000
nstenergy               = 5000
nstlog                  = 5000

continuation            = yes
constraint_algorithm    = lincs
constraints             = h-bonds
lincs_iter              = 2
lincs_order             = 4

cutoff-scheme           = Verlet
nstlist                 = 20
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz

; Energy groups for ion interaction analysis
energygrps              = {energy_groups}

tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 0.5
ref_t                   = {T}

pcoupl                  = no

gen_vel                 = no
"""
        )
        
        self.logger.info(f"Generated MDP files in {self.stage_dir}")
    
    def create_polymer_itp(
        self,
        acpype_dir: Path,
        resname: str,
        posre_fc: float = 1000.0
    ) -> Path:
        """Create polymer ITP file from acpype output."""
        gmx_itp = acpype_dir / f"{resname}_GMX.itp"
        poly_itp = self.stage_dir / "POLY.itp"

        # If the canonical ITP from acpype isn't present, try to locate any ITP in the provided directory
        if not gmx_itp.exists():
            itp_candidates = list(acpype_dir.glob("*.itp")) if acpype_dir.exists() else []
            # prefer file containing the residue name
            chosen = None
            for p in itp_candidates:
                if resname.lower() in p.name.lower():
                    chosen = p
                    break
            if not chosen and itp_candidates:
                chosen = itp_candidates[0]
            if chosen:
                gmx_itp = chosen
            else:
                raise FileNotFoundError(f"No ITP file found in {acpype_dir} or {acpype_dir} contains no .itp files")

        if checkpoint_file(poly_itp, "Create polymer ITP with posre", self.logger):
            return poly_itp

        # Read original ITP
        original_content = gmx_itp.read_text()
        
        # Create wrapper with position restraints
        wrapper_content = f"""#define {resname}_POSRES
#include "{gmx_itp.name}"

#ifdef POSRES
#include "POLY_posre.itp"
#endif
"""
        
        write_file(poly_itp, wrapper_content)
        
        # Generate position restraint ITP (simplified)
        posre_content = f"""[ position_restraints ]
; atom  type      fx      fy      fz
1       1    {posre_fc}    {posre_fc}    {posre_fc}
"""
        
        write_file(self.stage_dir / "POLY_posre.itp", posre_content)
        
        self.logger.info(f"Created polymer ITP: {poly_itp}")
        return poly_itp
    
    def run_all(self, acpype_dir: Path, resname: str) -> None:
        """Execute GROMACS prep pipeline."""
        self.logger.info("="*60)
        self.logger.info("GROMACS PREPARATION STAGE")
        self.logger.info("="*60)
        
        self.generate_mdp_files(resname)
        self.create_polymer_itp(acpype_dir, resname)
        
        self.logger.info("GROMACS prep complete")
