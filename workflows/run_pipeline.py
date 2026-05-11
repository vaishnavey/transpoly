#!/usr/bin/env python3
"""
Main entry point for transpoly pipeline.
Usage: python run_pipeline.py --config config.yaml
"""
import argparse
import sys
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

from transpoly.config import SimulationConfig
from transpoly.pipeline import TranspolyPipeline


def main():
    parser = argparse.ArgumentParser(description="Transpoly: Polymer/Salt MD Simulation Pipeline")
    parser.add_argument("--config", required=True, help="Configuration YAML file")
    
    args = parser.parse_args()
    config_file = Path(args.config)
    
    if not config_file.exists():
        print(f"Error: Config file not found: {config_file}")
        sys.exit(1)
    
    # Load configuration
    config = SimulationConfig.from_yaml(config_file)
    
    # Run pipeline
    pipeline = TranspolyPipeline(config)
    pipeline.run_all()


if __name__ == "__main__":
    main()
