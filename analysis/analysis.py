"""
Analysis modules: MSD, RDF, coordination, diffusivity, energy tracking.
"""
import logging
import numpy as np
from pathlib import Path
from typing import Tuple, List
import matplotlib.pyplot as plt
from scipy import stats
import subprocess
import os


def read_xvg(path: Path) -> Tuple[List[float], List[float]]:
    """Parse GROMACS XVG file."""
    x, y = [], []
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    x.append(float(parts[0]))
                    y.append(float(parts[1]))
    except Exception as e:
        raise ValueError(f"Error reading {path}: {e}")
    return x, y


def run_gmx_command(cmd: str, cwd: Path, logger: logging.Logger, description: str = "") -> bool:
    """Execute a GROMACS command and log output."""
    try:
        logger.info(f"Running: {description}")
        logger.info(f"  Command: {cmd}")
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=str(cwd),
            capture_output=True,
            text=True
        )
        if result.returncode != 0:
            logger.warning(f"  Warning: Command returned exit code {result.returncode}")
            if result.stderr:
                logger.warning(f"  stderr: {result.stderr[:200]}")
            return False
        logger.info(f"  ✓ {description} complete")
        return True
    except Exception as e:
        logger.error(f"  Error running command: {e}")
        return False


class MSDAnalysis:
    """Mean Square Displacement analysis."""
    
    @staticmethod
    def clean_trajectory(output_dir: Path, logger: logging.Logger) -> bool:
        """
        Clean PBC and rotation from production trajectory.
        
        Steps:
        1. Center system (to remove whole-molecule rotation)
        2. Remove PBC artifacts (noPBC option for molecules to not wrap)
        """
        prod_dir = output_dir / "06_production"
        
        tpr = prod_dir / "nvt_prod.tpr"
        xtc_raw = prod_dir / "nvt_prod.xtc"
        xtc_center = prod_dir / "nvt_prod_center.xtc"
        xtc_clean = prod_dir / "nvt_prod_clean.xtc"
        
        if not xtc_raw.exists():
            logger.error(f"Production trajectory not found: {xtc_raw}")
            return False
        
        # Step 1: Center on polymer to remove rotation
        if not xtc_center.exists():
            cmd = f"printf 'Protein\nSystem\n' | gmx_mpi trjconv -f {xtc_raw.name} -s {tpr.name} -center -o {xtc_center.name} -pbc none"
            if not run_gmx_command(cmd, prod_dir, logger, "Center trajectory on polymer"):
                logger.warning("Could not center trajectory; continuing anyway")
        
        # Step 2: Remove PBC wrapping artifacts
        if not xtc_clean.exists():
            cmd = f"printf 'System\n' | gmx_mpi trjconv -f {xtc_center.name if xtc_center.exists() else xtc_raw.name} -s {tpr.name} -o {xtc_clean.name} -pbc nojump"
            if not run_gmx_command(cmd, prod_dir, logger, "Remove PBC wrapping"):
                logger.warning("Could not clean trajectory; trying without centering")
                xtc_clean_alt = prod_dir / "nvt_prod_nojump.xtc"
                cmd = f"printf 'System\n' | gmx_mpi trjconv -f {xtc_raw.name} -s {tpr.name} -o {xtc_clean_alt.name} -pbc nojump"
                run_gmx_command(cmd, prod_dir, logger, "Remove PBC wrapping (direct from raw)")
        
        logger.info("Trajectory cleaning complete")
        return True
    
    @staticmethod
    def extract_msd(output_dir: Path, logger: logging.Logger) -> dict:
        """
        Extract MSD for K+ and Cl-.
        Requires cleaned trajectory (PBC and rotation removed).
        """
        prod_dir = output_dir / "06_production"
        
        # First clean the trajectory
        MSDAnalysis.clean_trajectory(output_dir, logger)
        
        # Determine which cleaned trajectory exists
        xtc_clean = prod_dir / "nvt_prod_clean.xtc"
        xtc_nojump = prod_dir / "nvt_prod_nojump.xtc"
        xtc_center = prod_dir / "nvt_prod_center.xtc"
        
        if xtc_clean.exists():
            xtc_to_use = xtc_clean
        elif xtc_nojump.exists():
            xtc_to_use = xtc_nojump
        elif xtc_center.exists():
            xtc_to_use = xtc_center
        else:
            xtc_to_use = prod_dir / "nvt_prod.xtc"
        
        logger.info(f"Using trajectory for MSD: {xtc_to_use.name}")
        
        results = {}
        
        # Make ion index if needed
        tpr = prod_dir / "nvt_prod.tpr"
        ndx = prod_dir / "index_ions.ndx"
        
        if not ndx.exists():
            cmd = "printf 'keep 0\nr K\nname 1 K_ions\nr CL\nname 2 CL_ions\nq\n' | gmx_mpi make_ndx -f nvt_prod.tpr -o index_ions.ndx"
            run_gmx_command(cmd, prod_dir, logger, "Build ion index")
        
        # Extract MSD for each ion
        for ion_name, ion_idx in [("K", "1"), ("Cl", "2")]:
            msd_xvg = prod_dir / f"msd_{ion_name.lower()}.xvg"
            
            if not msd_xvg.exists():
                cmd = f"printf '{ion_idx}\n' | gmx_mpi msd -f {xtc_to_use.name} -s nvt_prod.tpr -n index_ions.ndx -o msd_{ion_name.lower()}.xvg -tu ps"
                if run_gmx_command(cmd, prod_dir, logger, f"MSD for {ion_name}+ ions"):
                    logger.info(f"  MSD file created: {msd_xvg.name}")
            
            try:
                times, msd_vals = read_xvg(msd_xvg)
                results[ion_name] = (np.array(times), np.array(msd_vals))
                logger.info(f"  Loaded MSD: {len(times)} time points for {ion_name}+")
            except Exception as e:
                logger.warning(f"Could not read MSD for {ion_name}: {e}")
        
        return results
    
    @staticmethod
    def estimate_diffusivity(times: np.ndarray, msd: np.ndarray, logger: logging.Logger) -> float:
        """
        Estimate diffusivity from linear region of MSD-time plot.
        D = MSD / (6 * t)
        """
        if len(times) < 10:
            logger.warning("Not enough data points for diffusivity fit")
            return 0.0
        
        # Use last 50% of trajectory for linear fit
        n_start = len(times) // 2
        t_fit = times[n_start:]
        msd_fit = msd[n_start:]
        
        # Linear regression
        slope, intercept, r_value, p_value, std_err = stats.linregress(t_fit, msd_fit)
        
        # D = slope / 6 (in ps and nm^2)
        # slope is nm^2/ps
        D_msd = slope / 6.0  # nm^2/ps
        
        logger.info(f"MSD linear fit: slope={slope:.4f}, R²={r_value**2:.4f}")
        logger.info(f"Diffusivity from MSD: {D_msd:.4e} nm²/ps")
        
        return D_msd
    
    @staticmethod
    def plot_msd(output_dir: Path, results: dict, logger: logging.Logger) -> None:
        """Plot MSD for all ions."""
        analysis_dir = output_dir / "07_analysis"
        analysis_dir.mkdir(exist_ok=True)
        
        plt.figure(figsize=(8, 5))
        
        for ion_name, (times, msd) in results.items():
            plt.plot(times, msd, label=f"{ion_name}⁺", linewidth=1.5)
        
        plt.xlabel("Time (ps)")
        plt.ylabel("MSD (nm²)")
        plt.title("Ion Mean Square Displacement vs Time")
        plt.legend()
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig(analysis_dir / "msd_ions.png", dpi=150)
        plt.close()
        
        logger.info(f"Saved: msd_ions.png")


class RDFAnalysis:
    """Radial Distribution Function analysis."""
    
    @staticmethod
    def extract_rdf(output_dir: Path, ion_type: str, logger: logging.Logger) -> Tuple[np.ndarray, np.ndarray]:
        """
        Extract RDF for a given ion (e.g., 'K', 'Cl').
        Runs gmx rdf if the file doesn't exist.
        """
        prod_dir = output_dir / "06_production"
        
        # Use cleaned trajectory
        xtc_clean = prod_dir / "nvt_prod_clean.xtc"
        xtc_nojump = prod_dir / "nvt_prod_nojump.xtc"
        xtc_center = prod_dir / "nvt_prod_center.xtc"
        
        if xtc_clean.exists():
            xtc_to_use = xtc_clean
        elif xtc_nojump.exists():
            xtc_to_use = xtc_nojump
        elif xtc_center.exists():
            xtc_to_use = xtc_center
        else:
            xtc_to_use = prod_dir / "nvt_prod.xtc"
        
        rdf_xvg = prod_dir / f"rdf_{ion_type.lower()}_water.xvg"
        tpr = prod_dir / "nvt_prod.tpr"
        
        # Map ion names to GROMACS residue names
        ion_resname = {"K": "K", "Cl": "CL"}.get(ion_type, ion_type)
        water_resname = "SOL"
        
        if not rdf_xvg.exists():
            logger.info(f"Computing RDF for {ion_type} - Water oxygen...")
            # gmx rdf: selection pairs are "residue_name and name"
            cmd = f"gmx_mpi rdf -f {xtc_to_use.name} -s {tpr.name} -ref 'resname {ion_resname}' -sel 'resname {water_resname} and name OW' -o {rdf_xvg.name} -bin 0.01"
            if not run_gmx_command(cmd, prod_dir, logger, f"RDF {ion_type}-Water"):
                logger.warning(f"Could not compute RDF for {ion_type}; checking if file was created")
        
        if not rdf_xvg.exists():
            logger.warning(f"RDF file not found: {rdf_xvg}")
            return np.array([]), np.array([])
        
        try:
            r, rdf = read_xvg(rdf_xvg)
            logger.info(f"  Loaded RDF: {len(r)} distance points for {ion_type}-Water")
            return np.array(r), np.array(rdf)
        except Exception as e:
            logger.error(f"Error reading RDF: {e}")
            return np.array([]), np.array([])
    
    @staticmethod
    def plot_rdf(output_dir: Path, logger: logging.Logger) -> None:
        """Plot RDF profiles."""
        analysis_dir = output_dir / "07_analysis"
        analysis_dir.mkdir(exist_ok=True)
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        
        has_data = False
        for ion_type, ax in [("K", ax1), ("Cl", ax2)]:
            r, rdf = RDFAnalysis.extract_rdf(output_dir, ion_type, logger)
            
            if len(r) > 0:
                ax.plot(r, rdf, linewidth=1.5, color="steelblue")
                ax.set_xlabel("r (nm)")
                ax.set_ylabel("g(r)")
                ax.set_title(f"RDF: {ion_type} - Water Oxygen")
                ax.grid(alpha=0.3)
                has_data = True
            else:
                ax.text(0.5, 0.5, "No data", ha='center', va='center', transform=ax.transAxes)
        
        if has_data:
            plt.tight_layout()
            plt.savefig(analysis_dir / "rdf_ions.png", dpi=150)
            plt.close()
            logger.info("Saved: rdf_ions.png")
        else:
            plt.close()
            logger.warning("No RDF data to plot")


class CoordinationAnalysis:
    """Ion coordination number tracking."""
    
    @staticmethod
    def compute_coordination_numbers(output_dir: Path, logger: logging.Logger) -> dict:
        """
        Compute average coordination numbers.
        In full implementation, would use gmx select or MDAnalysis.
        """
        results = {}
        
        # Placeholder: coordination number computation would require
        # trajectory analysis. Here we document the approach.
        
        logger.info("Coordination number analysis would track neighbors over time")
        
        return results


class EnergyAnalysis:
    """Track LJ-SR and Coulomb interactions over time."""

    ION_NAMES = {"K", "CL", "NH4", "NH4P"}

    @staticmethod
    def _list_energy_terms(edr_file: Path, cwd: Path, logger: logging.Logger) -> list[str]:
        """Return available energy term labels from gmx energy."""
        proc = subprocess.Popen(
            ["gmx", "energy", "-f", str(edr_file), "-xvg", "none"],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=str(cwd),
            text=True,
        )
        stdout, stderr = proc.communicate(input="0\n")
        combined = f"{stdout}\n{stderr}"
        terms = []
        for line in combined.splitlines():
            line = line.strip()
            if not line or not line[0].isdigit():
                continue
            parts = line.split(maxsplit=1)
            if len(parts) == 2:
                terms.append(parts[1].strip())
        if not terms:
            logger.warning("Could not parse energy term list from gmx energy")
        return terms

    @staticmethod
    def _classify_pair(term: str) -> tuple[str, str] | None:
        """Map a term like Coulomb-(SR):K-SOL to (series_kind, legend_label)."""
        if ":" not in term:
            return None
        head, pair = term.split(":", 1)
        if head == "Coulomb-(SR)":
            series_kind = "coul_sr"
        elif head == "LJ-(SR)":
            series_kind = "lj_sr"
        else:
            return None

        if "-" not in pair:
            return None
        left, right = pair.split("-", 1)
        left = left.strip().upper()
        right = right.strip().upper()

        ion = None
        partner = None
        if left in EnergyAnalysis.ION_NAMES:
            ion = left
            partner = right
        elif right in EnergyAnalysis.ION_NAMES:
            ion = right
            partner = left

        if ion is None:
            return None

        ion_label = "NH4" if ion in {"NH4", "NH4P"} else ion
        if partner == "SOL":
            return series_kind, f"{ion_label}-Water"
        return series_kind, f"{ion_label}-Polymer"
    
    @staticmethod
    def extract_energy_terms(output_dir: Path, logger: logging.Logger) -> dict:
        """Extract ion-water and ion-polymer Coul-SR/LJ-SR energy terms."""
        prod_dir = output_dir / "06_production"
        analysis_dir = output_dir / "07_analysis"
        analysis_dir.mkdir(exist_ok=True)

        results: dict[str, list[tuple[str, Path]]] = {"coul_sr": [], "lj_sr": []}

        # Try both prod.edr and nvt_prod.edr
        edr_file = None
        if (prod_dir / "prod.edr").exists():
            edr_file = prod_dir / "prod.edr"
        elif (prod_dir / "nvt_prod.edr").exists():
            edr_file = prod_dir / "nvt_prod.edr"
        else:
            logger.warning("No .edr file found in production directory, skipping energy analysis")
            return results

        try:
            terms = EnergyAnalysis._list_energy_terms(edr_file, prod_dir, logger)
            extracted_labels: set[str] = set()
            for term in terms:
                classified = EnergyAnalysis._classify_pair(term)
                if classified is None:
                    continue
                series_kind, legend_label = classified
                key = f"{series_kind}:{legend_label}"
                if key in extracted_labels:
                    continue

                safe_label = legend_label.lower().replace("-", "_")
                out_file = analysis_dir / f"{series_kind}_{safe_label}.xvg"
                if not out_file.exists():
                    proc = subprocess.Popen(
                        ["gmx", "energy", "-f", str(edr_file), "-o", str(out_file)],
                        stdin=subprocess.PIPE,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE,
                        cwd=str(prod_dir),
                        text=True,
                    )
                    _, stderr = proc.communicate(input=f"{term}\n0\n")
                    if proc.returncode != 0:
                        logger.warning(f"Failed to extract {term}: {stderr}")
                        continue

                results[series_kind].append((legend_label, out_file))
                extracted_labels.add(key)

        except Exception as e:
            logger.warning(f"Energy extraction error: {e}")

        return results
    
    @staticmethod
    def plot_energies(energy_files: dict, output_dir: Path, logger: logging.Logger) -> None:
        """Plot all ion-water and ion-polymer curves on one Coul-SR/LJ-SR panel each."""
        analysis_dir = output_dir / "07_analysis"
        analysis_dir.mkdir(exist_ok=True)

        for base_type in ("coul_sr", "lj_sr"):
            entries = energy_files.get(base_type, [])
            if not entries:
                logger.warning(f"No {base_type} energy terms were extracted")
                continue

            try:
                plt.figure(figsize=(12, 6))

                plotted = 0
                for legend_label, xvg_file in entries:
                    if not xvg_file.exists():
                        continue
                    data = []
                    with open(xvg_file, "r") as f:
                        for line in f:
                            if not line.startswith(("@", "#")):
                                try:
                                    parts = line.split()
                                    if len(parts) >= 2:
                                        data.append([float(parts[0]), float(parts[1])])
                                except ValueError:
                                    continue

                    if not data:
                        logger.warning(f"No valid data in {xvg_file}")
                        continue

                    data = np.array(data)
                    plt.plot(data[:, 0], data[:, 1], linewidth=1.8, label=legend_label, alpha=0.9)
                    plotted += 1

                if plotted == 0:
                    plt.close()
                    logger.warning(f"No plottable traces for {base_type}")
                    continue

                if base_type == "coul_sr":
                    plt.title("Coulomb Short-Range Energy: Ion-Water and Ion-Polymer")
                else:
                    plt.title("Lennard-Jones Short-Range Energy: Ion-Water and Ion-Polymer")

                plt.xlabel("Time (ps)")
                plt.ylabel("Energy (kJ/mol)")
                plt.legend()
                plt.grid(alpha=0.25)
                plt.tight_layout()

                output_file = analysis_dir / f"{base_type}.png"
                plt.savefig(output_file, dpi=150)
                plt.close()
                logger.info(f"Energy plot saved: {output_file}")

            except Exception as e:
                logger.warning(f"Error plotting {base_type}: {e}")


class DiffusivityEstimator:
    """Estimate diffusivity with optional Stokes-Einstein refinement."""
    
    @staticmethod
    def stokes_einstein_diffusivity(
        temperature: float,
        viscosity: float,
        radius_nm: float
    ) -> float:
        """
        D = kT / (6πηr)
        
        Args:
            temperature: K
            viscosity: Pa·s (= g/(cm·s))
            radius_nm: nm
        
        Returns:
            D in nm²/ps
        """
        k_b = 1.380649e-23  # J/K
        T = temperature
        eta = viscosity  # Pa·s
        r = radius_nm * 1e-9  # m
        
        D_si = k_b * T / (6 * np.pi * eta * r)  # m²/s
        D_nm_ps = D_si * 1e12 * 1e-18  # nm²/ps
        
        return D_nm_ps
    
    @staticmethod
    def summary_report(output_dir: Path, results: dict, logger: logging.Logger) -> None:
        """Generate summary report of all diffusivity estimates."""
        analysis_dir = output_dir / "07_analysis"
        analysis_dir.mkdir(exist_ok=True)
        
        report_file = analysis_dir / "diffusivity_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write("="*60 + "\n")
            f.write("ION DIFFUSIVITY SUMMARY\n")
            f.write("="*60 + "\n\n")
            
            for ion_name, data in results.items():
                f.write(f"Ion: {ion_name}\n")
                if 'D_msd' in data:
                    f.write(f"  Diffusivity (from MSD): {data['D_msd']:.4e} nm²/ps\n")
                if 'D_se' in data:
                    f.write(f"  Diffusivity (Stokes-Einstein): {data['D_se']:.4e} nm²/ps\n")
                f.write("\n")
        
        logger.info(f"Wrote diffusivity summary: {report_file}")
