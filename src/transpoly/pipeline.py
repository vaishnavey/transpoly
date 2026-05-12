"""
Main pipeline orchestrator: coordinates all stages.
"""
import logging
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from .config import SimulationConfig
from .utils import setup_logger, get_topology_resname
from .parameterize import ParameterizeStage
from .packing import PackingStage
from .gromacs_prep import GromacsPrepStage
from .solvation import SolvationStage
from .solvent_only import SolventOnlyStage
from .equilibration import EquilibrationStage
from .production import ProductionStage


class TranspolyPipeline:
    """Main orchestrator for transpoly workflow."""
    
    def __init__(self, config: SimulationConfig):
        self.config = config
        self.output_dir = Path(config.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger = setup_logger("transpoly", self.output_dir)
        self.logger.info("TranspolyPipeline initialized")
        self.logger.info(f"Output directory: {self.output_dir}")

    def _is_solvent_only(self) -> bool:
        return self.config.workflow_mode.strip().lower() in {"solvent_only", "water_only", "solvent", "build_only"}
    
    def run_parameterization(self) -> tuple[Path, str, Path]:
        """Run parameterization stage."""
        stage = ParameterizeStage(self.config, self.output_dir, self.logger)
        stage.run_all()

        gro_file, itp_file, top_file = stage.extract_gmx_files()
        # Directory to use for GROMACS prep: prefer directory containing the ITP
        if itp_file is not None:
            acpype_dir = itp_file.parent
        else:
            acpype_dir = top_file.parent

        # Infer residue name
        resname = get_topology_resname(stage.stage_dir / f"{stage.get_base_name()}_gaff2.mol2", self.logger)

        return gro_file, resname, acpype_dir
    
    def run_packing(self, pdb_file: Path) -> Path:
        """Run packing stage."""
        stage = PackingStage(self.config, self.output_dir, self.logger)
        packed_pdb = stage.run_all(pdb_file)
        return packed_pdb
    
    def run_gromacs_prep(self, acpype_dir: Path, resname: str) -> None:
        """Run GROMACS preparation stage."""
        stage = GromacsPrepStage(self.config, self.output_dir, self.logger)
        stage.run_all(acpype_dir, resname)
    
    def run_solvation(
        self,
        packed_pdb: Path,
        acpype_dir: Path,
        resname: str,
        n_polymer: int
    ) -> tuple[Path, Path]:
        """Run solvation and ionization stage."""
        stage = SolvationStage(self.config, self.output_dir, self.logger)
        solvated_gro, topol = stage.run_all(packed_pdb, acpype_dir, resname, n_polymer)
        return solvated_gro, topol
    
    def run_equilibration(self) -> Path:
        """Run equilibration stage."""
        stage = EquilibrationStage(self.config, self.output_dir, self.logger)
        final_gro = stage.run_all()
        return final_gro
    
    def run_production(self) -> Path:
        """Run production MD stage."""
        stage = ProductionStage(self.config, self.output_dir, self.logger)
        prod_xtc = stage.run_all()
        return prod_xtc

    def _analyze_output_dir(self, output_dir: Path, logger: logging.Logger) -> dict:
        """Run analysis against an arbitrary output directory and return diffusivity values."""
        from analysis.analysis import MSDAnalysis, RDFAnalysis, DiffusivityEstimator, EnergyAnalysis

        logger.info("="*60)
        logger.info("ANALYSIS STAGE")
        logger.info("="*60)

        energy_files = EnergyAnalysis.extract_energy_terms(output_dir, logger)
        if energy_files:
            EnergyAnalysis.plot_energies(energy_files, output_dir, logger)

        msd_results = MSDAnalysis.extract_msd(output_dir, logger)
        diffusivity_estimates = {}
        for ion_name, (times, msd) in msd_results.items():
            D = MSDAnalysis.estimate_diffusivity(times, msd, logger)
            diffusivity_estimates[ion_name] = {"D_msd": D}

        if msd_results:
            MSDAnalysis.plot_msd(output_dir, msd_results, logger)

        RDFAnalysis.plot_rdf(output_dir, logger)
        DiffusivityEstimator.summary_report(output_dir, diffusivity_estimates, logger)
        logger.info("Analysis complete")
        return diffusivity_estimates

    def _prepare_replicate_inputs(self, replicate_dir: Path) -> None:
        """Copy prepared inputs so each replicate can run equilibration+production independently."""
        for stage_name in ("03_gromacs_prep", "04_solvation"):
            src = self.output_dir / stage_name
            dst = replicate_dir / stage_name
            if not src.exists():
                raise FileNotFoundError(f"Required stage directory missing: {src}")
            if dst.exists():
                shutil.rmtree(dst)
            shutil.copytree(src, dst)

    def _run_single_replicate(self, replicate_idx: int) -> dict:
        """Run one equilibration+production+analysis replicate."""
        replicate_dir = self.output_dir / f"replicate_{replicate_idx:02d}"
        replicate_dir.mkdir(parents=True, exist_ok=True)
        rep_logger = setup_logger(f"transpoly_rep{replicate_idx:02d}", replicate_dir)

        self._prepare_replicate_inputs(replicate_dir)

        equil_stage = EquilibrationStage(self.config, replicate_dir, rep_logger)
        equil_stage.run_all()

        prod_stage = ProductionStage(self.config, replicate_dir, rep_logger)
        prod_stage.run_all()

        return self._analyze_output_dir(replicate_dir, rep_logger)

    def _run_independent_replicates(self) -> None:
        """Launch N independent equilibration+production runs in parallel and average diffusivity."""
        n_runs = max(1, int(self.config.independent_run_count))
        self.logger.info(f"Launching {n_runs} independent replicate runs in parallel")

        replicate_results: list[tuple[int, dict]] = []
        with ThreadPoolExecutor(max_workers=n_runs) as executor:
            futures = {
                executor.submit(self._run_single_replicate, idx): idx
                for idx in range(1, n_runs + 1)
            }
            for future in as_completed(futures):
                idx = futures[future]
                result = future.result()
                replicate_results.append((idx, result))
                self.logger.info(f"Replicate {idx:02d} finished")

        analysis_dir = self.output_dir / "07_analysis"
        analysis_dir.mkdir(parents=True, exist_ok=True)

        per_ion_values: dict[str, list[float]] = {}
        for _, result in replicate_results:
            for ion_name, data in result.items():
                D = data.get("D_msd")
                if D is None:
                    continue
                per_ion_values.setdefault(ion_name, []).append(float(D))

        import numpy as np

        summary_file = analysis_dir / "diffusivity_summary.txt"
        with open(summary_file, "w", encoding="utf-8") as handle:
            handle.write("="*60 + "\n")
            handle.write("ION DIFFUSIVITY SUMMARY (INDEPENDENT RUNS)\n")
            handle.write("="*60 + "\n\n")
            handle.write(f"Replicate count: {n_runs}\n\n")

            for ion_name in sorted(per_ion_values):
                values = np.array(per_ion_values[ion_name], dtype=float)
                mean_val = float(np.mean(values))
                std_val = float(np.std(values, ddof=0))
                handle.write(f"Ion: {ion_name}\n")
                handle.write(f"  Mean diffusivity (MSD): {mean_val:.4e} nm²/ps\n")
                handle.write(f"  Std diffusivity (MSD):  {std_val:.4e} nm²/ps\n")
                handle.write(f"  Samples: {len(values)}\n\n")

        self.logger.info(f"Wrote independent-runs diffusivity summary: {summary_file}")
    
    def run_analysis(self) -> None:
        """Run post-production analysis."""
        self._analyze_output_dir(self.output_dir, self.logger)
    
    def run_all(self) -> None:
        """Execute full pipeline from start to finish."""
        self.logger.info("="*70)
        self.logger.info(" TRANSPOLY: POLYMER/SALT MD SIMULATION PIPELINE")
        self.logger.info("="*70)
        
        try:
            if self._is_solvent_only():
                self.logger.info("Workflow mode: solvent_only")
                stage = SolventOnlyStage(self.config, self.output_dir, self.logger)
                stage.build_and_run()
                self.logger.info("="*70)
                self.logger.info("PIPELINE COMPLETED SUCCESSFULLY")
                self.logger.info("="*70)
                self.logger.info(f"Results in: {self.output_dir}")
                return

            if not self.config.single_chain_pdb:
                raise ValueError("single_chain_pdb is required for polymer workflow mode")

            # 1. Parameterize
            self.logger.info("\n[1/6] Parameterization")
            pdb_input = Path(self.config.single_chain_pdb)
            gro_file, resname, acpype_dir = self.run_parameterization()
            
            # 2. Pack
            self.logger.info("\n[2/6] Packing")
            packed_pdb = self.run_packing(pdb_input)
            n_chains = self.config.n_chains or self.config.estimate_n_chains()
            
            # 3. GROMACS prep
            self.logger.info("\n[3/6] GROMACS Preparation")
            self.run_gromacs_prep(acpype_dir, resname)
            
            # 4. Solvation & Ions
            self.logger.info("\n[4/6] Solvation & Ionization")
            solvated_gro, topol = self.run_solvation(packed_pdb, acpype_dir, resname, n_chains)
            
            # 5. Equilibration
            if self.config.independent_runs:
                self.logger.info("\n[5/6] Parallel Independent Equilibration + Production")
                self._run_independent_replicates()
            else:
                self.logger.info("\n[5/6] Equilibration")
                self.run_equilibration()

                # 6. Production
                self.logger.info("\n[6/6] Production MD")
                self.run_production()

                # 7. Analysis
                self.logger.info("\n[7/7] Analysis")
                self.run_analysis()
            
            self.logger.info("\n" + "="*70)
            self.logger.info("PIPELINE COMPLETED SUCCESSFULLY")
            self.logger.info("="*70)
            self.logger.info(f"Results in: {self.output_dir}")
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {e}", exc_info=True)
            raise
