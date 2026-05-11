#!/usr/bin/env python3
"""
Verify transpoly installation and structure.
"""
import sys
from pathlib import Path


def check_structure():
    """Verify repo structure."""
    repo_root = Path(__file__).parent.parent
    
    required_dirs = [
        "src/transpoly",
        "workflows",
        "templates",
        "analysis",
        "examples",
        "docs",
    ]
    
    required_files = [
        "README.md",
        "requirements.txt",
        ".copilot-instructions.md",
        "LICENSE",
        "src/transpoly/__init__.py",
        "src/transpoly/config.py",
        "src/transpoly/utils.py",
        "src/transpoly/parameterize.py",
        "src/transpoly/packing.py",
        "src/transpoly/gromacs_prep.py",
        "src/transpoly/solvation.py",
        "src/transpoly/equilibration.py",
        "src/transpoly/production.py",
        "src/transpoly/pipeline.py",
        "analysis/analysis.py",
        "workflows/run_pipeline.py",
        "docs/SETUP.md",
        "docs/USAGE.md",
        "docs/ARCHITECTURE.md",
        "examples/config_example.yaml",
    ]
    
    print("Checking repository structure...")
    print()
    
    missing = []
    
    for d in required_dirs:
        p = repo_root / d
        if p.is_dir():
            print(f"✓ {d}/")
        else:
            print(f"✗ {d}/ (missing)")
            missing.append(d)
    
    print()
    
    for f in required_files:
        p = repo_root / f
        if p.is_file():
            print(f"✓ {f}")
        else:
            print(f"✗ {f} (missing)")
            missing.append(f)
    
    print()
    
    if missing:
        print(f"ERROR: {len(missing)} items missing")
        return False
    else:
        print("✓ All required files and directories present")
        return True


def check_imports():
    """Verify Python modules can be imported."""
    print("Checking Python imports...")
    print()
    
    sys.path.insert(0, str(Path(__file__).parent / "src"))
    
    modules_to_check = [
        "transpoly",
        "transpoly.config",
        "transpoly.utils",
        "transpoly.parameterize",
        "transpoly.packing",
        "transpoly.gromacs_prep",
        "transpoly.solvation",
        "transpoly.equilibration",
        "transpoly.production",
        "transpoly.pipeline",
    ]
    
    failed = []
    
    for mod in modules_to_check:
        try:
            __import__(mod)
            print(f"✓ {mod}")
        except ImportError as e:
            print(f"✗ {mod}: {e}")
            failed.append(mod)
    
    print()
    
    if failed:
        print(f"ERROR: {len(failed)} modules failed to import")
        return False
    else:
        print("✓ All modules import successfully")
        return True


def main():
    print("="*60)
    print("TRANSPOLY INSTALLATION VERIFICATION")
    print("="*60)
    print()
    
    struct_ok = check_structure()
    print()
    import_ok = check_imports()
    
    print()
    print("="*60)
    
    if struct_ok and import_ok:
        print("✓ Transpoly installation verified successfully!")
        print()
        print("Next steps:")
        print("  1. Check environment: bash check_env.sh")
        print("  2. Edit config: examples/config_example.yaml")
        print("  3. Run pipeline: workflows/run_pipeline.py --config examples/config_example.yaml")
        print()
        return 0
    else:
        print("✗ Verification failed. See errors above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
