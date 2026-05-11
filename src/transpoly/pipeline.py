"""
Main pipeline orchestrator: coordinates all stages.
"""
import logging
from pathlib import Path
from .config import SimulationConfig
from .utils import setup_logger, get_topology_resname
from .parameterize import ParameterizeStage
from .packing import PackingStage
from .gromacs_prep import GromacsPrepStage
from .solvation import SolvationStage
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
    
    def run_parameterization(self) -> tuple[Path, str]:
        """Run parameterization stage."""
        stage = ParameterizeStage(self.config, self.output_dir, self.logger)
        stage.run_all()
        
        gro_file, itp_file, top_file = stage.extract_gmx_files()
        acpype_dir = stage.stage_dir / f"{stage.get_base_name()}.acpype"
        
        # Infer residue name
        resname = get_topology_resname(stage.stage_dir / f"{stage.get_base_name()}_gaff2.mol2", self.logger)
        
        return gro_file, resname
    
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
    
    def run_analysis(self) -> None:
        """Run post-production analysis."""
        from analysis.analysis import MSDAnalysis, RDFAnalysis, DiffusivityEstimator
        
        self.logger.info("="*60)
        self.logger.info("ANALYSIS STAGE")
        self.logger.info("="*60)
        
        # MSD analysis
        msd_results = MSDAnalysis.extract_msd(self.output_dir, self.logger)
        for ion_name, (times, msd) in msd_results.items():
            D = MSDAnalysis.estimate_diffusivity(times, msd, self.logger)
        
        MSDAnalysis.plot_msd(self.output_dir, msd_results, self.logger)
        
        # RDF analysis
        RDFAnalysis.plot_rdf(self.output_dir, self.logger)
        
        # Diffusivity summary
        results = {
            "K": {"D_msd": 1.0e-5},  # Placeholder
            "Cl": {"D_msd": 2.0e-5}
        }
        DiffusivityEstimator.summary_report(self.output_dir, results, self.logger)
        
        self.logger.info("Analysis complete")
    
    def run_all(self) -> None:
        """Execute full pipeline from start to finish."""
        self.logger.info("="*70)
        self.logger.info(" TRANSPOLY: POLYMER/SALT MD SIMULATION PIPELINE")
        self.logger.info("="*70)
        
        try:
            # 1. Parameterize
            self.logger.info("\n[1/6] Parameterization")
            pdb_input = Path(self.config.single_chain_pdb)
            gro_file, resname = self.run_parameterization()
            
            # 2. Pack
            self.logger.info("\n[2/6] Packing")
            packed_pdb = self.run_packing(pdb_input)
            n_chains = self.config.n_chains or self.config.estimate_n_chains()
            
            # 3. GROMACS prep
            self.logger.info("\n[3/6] GROMACS Preparation")
            acpype_dir = Path(self.output_dir) / "01_parameterization" / f"{pdb_input.stem}.acpype"
            self.run_gromacs_prep(acpype_dir, resname)
            
            # 4. Solvation & Ions
            self.logger.info("\n[4/6] Solvation & Ionization")
            solvated_gro, topol = self.run_solvation(packed_pdb, acpype_dir, resname, n_chains)
            
            # 5. Equilibration
            self.logger.info("\n[5/6] Equilibration")
            equil_gro = self.run_equilibration()
            
            # 6. Production
            self.logger.info("\n[6/6] Production MD")
            prod_xtc = self.run_production()
            
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
