"""
Analysis modules: MSD, RDF, coordination, diffusivity, energy tracking.
"""
import logging
import numpy as np
from pathlib import Path
from typing import Tuple, List
import matplotlib.pyplot as plt
from scipy import stats


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


class MSDAnalysis:
    """Mean Square Displacement analysis."""
    
    @staticmethod
    def extract_msd(output_dir: Path, logger: logging.Logger) -> dict:
        """Extract MSD for K+ and Cl-."""
        prod_dir = output_dir / "06_production"
        
        results = {}
        
        # Make ion index if needed
        tpr = prod_dir / "nvt_prod.tpr"
        make_ndx_script = prod_dir / "make_ndx_ions.sh"
        
        if not (prod_dir / "index_ions.ndx").exists():
            script = """#!/bin/bash
printf "keep 0\\nr K\\nname 1 K_ions\\nr CL\\nname 2 CL_ions\\nq\\n" | gmx_mpi make_ndx -f nvt_prod.tpr -o index_ions.ndx >/dev/null
"""
            make_ndx_script.write_text(script)
            logger.info("Generated make_ndx script")
        
        # Extract MSD for each ion
        for ion_name, ion_idx in [("K", "1"), ("Cl", "2")]:
            msd_xvg = prod_dir / f"msd_{ion_name.lower()}.xvg"
            
            if not msd_xvg.exists():
                logger.info(f"Computing MSD for {ion_name}...")
                # This would require running gmx msd in actual workflow
            
            try:
                times, msd_vals = read_xvg(msd_xvg)
                results[ion_name] = (np.array(times), np.array(msd_vals))
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
        Requires gmx rdf to be run separately.
        """
        prod_dir = output_dir / "06_production"
        rdf_xvg = prod_dir / f"rdf_{ion_type.lower()}_water.xvg"
        
        if not rdf_xvg.exists():
            logger.warning(f"RDF file not found: {rdf_xvg}")
            return np.array([]), np.array([])
        
        try:
            r, rdf = read_xvg(rdf_xvg)
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
        
        for ion_type, ax in [("K", ax1), ("Cl", ax2)]:
            r, rdf = RDFAnalysis.extract_rdf(output_dir, ion_type, logger)
            
            if len(r) > 0:
                ax.plot(r, rdf, linewidth=1.5, color="steelblue")
                ax.set_xlabel("r (nm)")
                ax.set_ylabel("g(r)")
                ax.set_title(f"RDF: {ion_type} - Water Oxygen")
                ax.grid(alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(analysis_dir / "rdf_ions.png", dpi=150)
        plt.close()
        
        logger.info("Saved: rdf_ions.png")


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
    
    @staticmethod
    def extract_energy_terms(output_dir: Path, ion_type: str, logger: logging.Logger) -> dict:
        """Extract LJ-SR and Coul-SR for ion-polymer and ion-water."""
        prod_dir = output_dir / "06_production"
        
        results = {}
        
        # Would use gmx energy to extract specific terms
        # e.g., "LJ-SR:K-Water", "Coul-SR:K-Polymer"
        
        logger.info(f"Energy analysis for {ion_type} interactions (placeholder)")
        
        return results
    
    @staticmethod
    def plot_energies(output_dir: Path, logger: logging.Logger) -> None:
        """Plot energy terms over time."""
        analysis_dir = output_dir / "07_analysis"
        analysis_dir.mkdir(exist_ok=True)
        
        logger.info("Energy tracking plot generation (placeholder)")


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
