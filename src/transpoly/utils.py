"""
Utility functions: logging, file I/O, subprocess management, checkpointing.
"""
import logging
import os
import subprocess
from pathlib import Path
from typing import Optional, Tuple
import sys
import time


def setup_logger(name: str, output_dir: Optional[Path] = None) -> logging.Logger:
    """Set up a logger that writes to both stdout and file."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    
    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    
    # Formatter
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    
    # File handler
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(output_dir / "pipeline.log")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    
    return logger


def write_file(path: Path, content: str) -> None:
    """Write content to file, creating parent directories if needed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")


def run_command(
    cmd: str,
    cwd: Path,
    logger: logging.Logger,
    check: bool = True,
    description: str = "",
    env: Optional[dict[str, str]] = None,
) -> Tuple[int, str, str]:
    """
    Run a shell command and return (returncode, stdout, stderr).
    
    Args:
        cmd: Shell command to run
        cwd: Working directory
        logger: Logger instance
        check: Raise exception if returncode != 0
        description: Optional description for logging
    
    Returns:
        (returncode, stdout, stderr)
    """
    if description:
        logger.info(f"Running: {description}")
    logger.debug(f"Command: {cmd}")
    logger.debug(f"CWD: {cwd}")
    
    cwd.mkdir(parents=True, exist_ok=True)
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    tmp_dir = cwd / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    run_env.setdefault("TMPDIR", str(tmp_dir))
    run_env.setdefault("TMP", str(tmp_dir))
    run_env.setdefault("TEMP", str(tmp_dir))
    run_env.setdefault("OMPI_MCA_orte_tmpdir_base", str(tmp_dir))
    run_env.setdefault("OMPI_MCA_shmem_mmap_enable_nfs_warning", "0")
    
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=str(cwd),
            env=run_env,
            capture_output=True,
            text=True,
            check=False
        )
        
        if result.stdout:
            logger.debug(f"STDOUT:\n{result.stdout}")
        if result.stderr:
            logger.debug(f"STDERR:\n{result.stderr}")
        
        if result.returncode != 0 and check:
            logger.error(f"Command failed with return code {result.returncode}")
            logger.error(f"STDOUT:\n{result.stdout}")
            logger.error(f"STDERR:\n{result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, cmd)
        
        return result.returncode, result.stdout, result.stderr
    
    except subprocess.CalledProcessError:
        raise
    except Exception as e:
        logger.error(f"Exception running command: {e}")
        raise


def checkpoint_exists(path: Path, logger: logging.Logger) -> bool:
    """Check if a checkpoint file exists and log it."""
    exists = path.exists()
    if exists:
        logger.info(f"Checkpoint found: {path}")
    return exists


def checkpoint_file(
    path: Path,
    description: str,
    logger: logging.Logger
) -> bool:
    """
    Check if checkpoint exists. If yes, log and return True.
    If no, log that we're starting this stage and return False.
    """
    if checkpoint_exists(path, logger):
        logger.info(f"Skipping: {description} (checkpoint exists)")
        return True
    logger.info(f"Starting: {description}")
    return False


def get_topology_resname(pdb_file: Path, logger: logging.Logger) -> str:
    """
    Infer residue name from PDB file.
    Returns the first ATOM residue name found.
    """
    try:
        with open(pdb_file) as f:
            for line in f:
                if line.startswith("ATOM"):
                    resname = line[17:20].strip()
                    logger.info(f"Inferred residue name: {resname}")
                    return resname
    except Exception as e:
        logger.error(f"Error reading PDB file {pdb_file}: {e}")
    
    # Fallback
    logger.warning("Could not infer residue name; using default 'POL'")
    return "POL"


def read_xvg(path: Path) -> Tuple[list, list]:
    """
    Parse GROMACS .xvg file and return (x_values, y_values).
    Skips comment lines starting with # or @.
    """
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
        raise ValueError(f"Error reading XVG file {path}: {e}")
    return x, y


def estimate_density_from_gro(gro_file: Path) -> float:
    """
    Rough estimate of system density from .gro file.
    Assumes box dimensions and atom count.
    """
    try:
        with open(gro_file) as f:
            lines = f.readlines()
        
        # Line 0: title
        # Line 1: number of atoms
        n_atoms = int(lines[1].strip())
        
        # Last line: box dimensions (x, y, z)
        box_line = lines[-1].strip().split()
        if len(box_line) >= 3:
            x, y, z = float(box_line[0]), float(box_line[1]), float(box_line[2])
            vol_nm3 = x * y * z
            vol_cm3 = vol_nm3 * 1e-24
            
            # Rough average mass per atom ~ 10 g/mol / avogadro
            avg_mass_per_atom = 10.0 / 6.022e23
            total_mass = n_atoms * avg_mass_per_atom
            density = total_mass / vol_cm3
            return density
    except Exception:
        pass
    
    return 1.0  # Default if calculation fails


def timestamp() -> str:
    """Return current timestamp for unique identifiers."""
    return time.strftime("%Y%m%d_%H%M%S")
