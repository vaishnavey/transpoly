"""
Solvent-only workflow for water + ions systems.
"""
from __future__ import annotations

import logging
import re
from pathlib import Path

import matplotlib.pyplot as plt

from .config import SimulationConfig
from .utils import run_command, write_file


class SolventOnlyStage:
    """Build and simulate a water + ions box without a polymer input."""

    def __init__(self, config: SimulationConfig, output_dir: Path, logger: logging.Logger):
        self.config = config
        self.output_dir = output_dir
        self.logger = logger
        self.stage_dir = output_dir / "04_solvation"
        self.stage_dir.mkdir(parents=True, exist_ok=True)
        self.analysis_dir = output_dir / "07_analysis"
        self.analysis_dir.mkdir(parents=True, exist_ok=True)

    def _write_ion_gro(self, path: Path, resname: str, atom_name: str) -> None:
        content = (
            f"{atom_name} ion\n"
            f"    1\n"
            f"{1:5d}{resname:<5}{atom_name:>5}{1:5d}{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}\n"
            f"   0.00000   0.00000   0.00000\n"
        )
        write_file(path, content)

    def _plot_xvg(self, xvg: Path, output_png: Path, title: str, xlabel: str, ylabel: str) -> None:
        xs: list[float] = []
        ys: list[float] = []
        with xvg.open(encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("@"):  # comment/header
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    xs.append(float(parts[0]))
                    ys.append(float(parts[1]))

        plt.figure(figsize=(7.2, 4.6))
        plt.plot(xs, ys, linewidth=1.8)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.title(title)
        plt.tight_layout()
        plt.savefig(output_png, dpi=220)
        plt.close()

    def _write_mdp_files(self) -> None:
        temperature = self.config.temperature
        pressure = self.config.pressure
        prod_steps = int((self.config.prod_time * 1000) / 2)

        write_file(
            self.stage_dir / "em.mdp",
            """integrator              = steep
nsteps                  = 2000
emtol                   = 1000.0
emstep                  = 0.01

cutoff-scheme           = Verlet
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz
constraints             = none
""",
        )
        write_file(
            self.stage_dir / "nvt.mdp",
            f"""integrator              = md
dt                      = 0.001
nsteps                  = 50000

nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 5000
nstenergy               = 1000
nstlog                  = 1000

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
tau_t                   = 0.5
ref_t                   = {temperature}

pcoupl                  = no
gen_vel                 = yes
gen_temp                = {temperature}
gen_seed                = -1
""",
        )
        write_file(
            self.stage_dir / "npt.mdp",
            f"""integrator              = md
dt                      = 0.001
nsteps                  = 50000

nstxout                 = 0
nstvout                 = 0
nstfout                 = 0
nstxout-compressed      = 5000
nstenergy               = 1000
nstlog                  = 1000

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
tau_t                   = 0.5
ref_t                   = {temperature}

pcoupl                  = C-rescale
pcoupltype              = isotropic
tau_p                   = 5.0
ref_p                   = {pressure}
compressibility         = 4.5e-5

DispCorr                = EnerPres
gen_vel                 = no
""",
        )
        write_file(
            self.stage_dir / "prod.mdp",
            f"""integrator              = md
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
coulombtype             = PME
rcoulomb                = 1.0
vdwtype                 = Cut-off
rvdw                    = 1.0
pbc                     = xyz

tcoupl                  = v-rescale
tc-grps                 = System
tau_t                   = 0.5
ref_t                   = {temperature}

pcoupl                  = no
gen_vel                 = no
""",
        )

    @staticmethod
    def _count_solvent_molecules(gro_path: Path) -> int:
        """Count distinct SOL residue IDs in a GRO file."""
        residues = set()
        with gro_path.open(encoding="utf-8") as handle:
            lines = handle.readlines()
        for line in lines[2:-1]:
            if len(line) >= 10 and line[5:10].strip() == "SOL":
                residues.add(line[:5].strip())
        return len(residues)

    @staticmethod
    def _rewrite_topology(topol_path: Path, n_solvent: int, n_k: int, n_cl: int) -> None:
        topol_path.write_text(
            f'''#include "oplsaa.ff/forcefield.itp"
#include "oplsaa.ff/tip4p.itp"
#include "oplsaa.ff/ions.itp"

[ system ]
Pure water and ions box

[ molecules ]
SOL     {n_solvent}
K       {n_k}
CL      {n_cl}
''',
            encoding="utf-8",
        )

    def _make_index(self, tpr_name: str) -> None:
        run_command(
            f"printf 'keep 0\nr K\nname 1 K_ions\nr CL\nname 2 CL_ions\nq\n' | gmx_mpi make_ndx -f {tpr_name} -o ions.ndx",
            self.stage_dir,
            self.logger,
            description="build ion index",
        )

    def _collect_msd(self, prod_tpr: str) -> None:
        self._make_index(prod_tpr)
        run_command(
            f"printf '1\n' | gmx_mpi msd -f prod.xtc -s {prod_tpr} -n ions.ndx -o msd_k.xvg -tu ps",
            self.stage_dir,
            self.logger,
            description="MSD K ions",
        )
        run_command(
            f"printf '2\n' | gmx_mpi msd -f prod.xtc -s {prod_tpr} -n ions.ndx -o msd_cl.xvg -tu ps",
            self.stage_dir,
            self.logger,
            description="MSD Cl ions",
        )

    def build_and_run(self) -> None:
        self.logger.info("=" * 60)
        self.logger.info("SOLVENT-ONLY WORKFLOW")
        self.logger.info("=" * 60)

        seed_gro = self.stage_dir / "seed.gro"
        box_gro = self.stage_dir / "box.gro"
        solv_gro = self.stage_dir / "solv.gro"
        solv_k_gro = self.stage_dir / "solv_k.gro"
        solv_kcl_gro = self.stage_dir / "solv_kcl.gro"
        topol_top = self.stage_dir / "topol.top"

        write_file(
            seed_gro,
            """Pure water and ions seed
    1
    1DUM    DUM    1   0.000   0.000   0.000
   0.00000   0.00000   0.00000
""",
        )
        self._write_ion_gro(self.stage_dir / "K_single.gro", "K", "K")
        self._write_ion_gro(self.stage_dir / "CL_single.gro", "CL", "CL")
        write_file(
            topol_top,
            """#include "oplsaa.ff/forcefield.itp"
#include "oplsaa.ff/tip4p.itp"
#include "oplsaa.ff/ions.itp"

[ system ]
Pure water and ions box

[ molecules ]
SOL     0
K       0
CL      0
""",
        )

        bx_nm = self.config.box_x / 10.0
        by_nm = self.config.box_y / 10.0
        bz_nm = self.config.box_z / 10.0

        run_command(
            f"gmx_mpi editconf -f {seed_gro.name} -o {box_gro.name} -box {bx_nm:.3f} {by_nm:.3f} {bz_nm:.3f} -bt triclinic",
            self.stage_dir,
            self.logger,
            description="build solvent box",
        )
        run_command(
            f"gmx_mpi solvate -cp {box_gro.name} -cs tip4p.gro -o {solv_gro.name} -p {topol_top.name}",
            self.stage_dir,
            self.logger,
            description="solvate box",
        )
        n_solvent = self._count_solvent_molecules(solv_gro)
        if not n_solvent:
            self.logger.warning("Could not derive solvent count from solvated GRO; using 0")

        n_k = int(self.config.kcl_count)
        n_cl = int(self.config.nh4cl_count)

        run_command(
            f"gmx_mpi insert-molecules -f {solv_gro.name} -ci K_single.gro -nmol {n_k} -radius 0.08 -o {solv_k_gro.name}",
            self.stage_dir,
            self.logger,
            description="insert K ions",
        )
        run_command(
            f"gmx_mpi insert-molecules -f {solv_k_gro.name} -ci CL_single.gro -nmol {n_cl} -radius 0.08 -o {solv_kcl_gro.name}",
            self.stage_dir,
            self.logger,
            description="insert Cl ions",
        )

        self._rewrite_topology(topol_top, n_solvent, n_k, n_cl)
        self._write_mdp_files()

        run_command(
            f"gmx_mpi grompp -f em.mdp -c {solv_kcl_gro.name} -p {topol_top.name} -o em.tpr -maxwarn 1",
            self.stage_dir,
            self.logger,
            description="grompp EM",
        )
        run_command(
            "gmx_mpi mdrun -deffnm em",
            self.stage_dir,
            self.logger,
            description="mdrun EM",
        )
        run_command(
            "printf 'Potential\n0\n' | gmx_mpi energy -f em.edr -o em_potential.xvg",
            self.stage_dir,
            self.logger,
            description="extract EM potential",
        )

        run_command(
            f"gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p {topol_top.name} -o nvt.tpr -maxwarn 1",
            self.stage_dir,
            self.logger,
            description="grompp NVT",
        )
        run_command(
            "gmx_mpi mdrun -deffnm nvt",
            self.stage_dir,
            self.logger,
            description="mdrun NVT",
        )

        run_command(
            f"gmx_mpi grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p {topol_top.name} -o npt.tpr -maxwarn 1",
            self.stage_dir,
            self.logger,
            description="grompp NPT",
        )
        run_command(
            "gmx_mpi mdrun -deffnm npt",
            self.stage_dir,
            self.logger,
            description="mdrun NPT",
        )
        run_command(
            "printf 'Density\n0\n' | gmx_mpi energy -f npt.edr -o npt_density.xvg",
            self.stage_dir,
            self.logger,
            description="extract NPT density",
        )

        run_command(
            f"gmx_mpi grompp -f prod.mdp -c npt.gro -t npt.cpt -p {topol_top.name} -o prod.tpr -maxwarn 1",
            self.stage_dir,
            self.logger,
            description="grompp production",
        )
        run_command(
            "gmx_mpi mdrun -deffnm prod",
            self.stage_dir,
            self.logger,
            description="mdrun production",
        )

        self._collect_msd("prod.tpr")

        self._plot_xvg(
            self.stage_dir / "em_potential.xvg",
            self.analysis_dir / "em_potential.png",
            "Minimization Potential Energy",
            "Time (ps)",
            "Potential Energy (kJ/mol)",
        )
        self._plot_xvg(
            self.stage_dir / "npt_density.xvg",
            self.analysis_dir / "npt_density.png",
            "NPT Density vs Time",
            "Time (ps)",
            "Density (kg/m^3)",
        )
        self._plot_xvg(
            self.stage_dir / "msd_k.xvg",
            self.analysis_dir / "msd_k.png",
            "K+ MSD vs Time",
            "Time (ps)",
            "MSD (nm^2)",
        )
        self._plot_xvg(
            self.stage_dir / "msd_cl.xvg",
            self.analysis_dir / "msd_cl.png",
            "Cl- MSD vs Time",
            "Time (ps)",
            "MSD (nm^2)",
        )

        summary = self.stage_dir / "summary.txt"
        write_file(
            summary,
            f"""Solvent-only build and MD complete.
Box (nm): {bx_nm:.3f} {by_nm:.3f} {bz_nm:.3f}
Water molecules: {n_solvent}
Inserted K ions: {n_k}
Inserted Cl ions: {n_cl}
Final build GRO: {solv_kcl_gro.name}
""",
        )

        self.logger.info("Solvent-only workflow complete")
