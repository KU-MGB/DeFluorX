#!/usr/bin/env python3

"""
===============================================================================
DeFluorX Pipeline  |  Step 00  |  Environment Management
===============================================================================
Exports the active Conda environment to 'PFAS.yml' and 'requirements.txt'
for reproducibility. Provides an automated installation routine to
re-synchronise environments across compute nodes.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 20 August 2026 <─────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 00_03_Environment_DeFluorX.py
  Role          : Infrastructure - environment export and sync.
  Imports from  : None (standalone sys/os/subprocess).
  Reads         : Active conda environment.
  Writes        : PFAS.yml, requirements.txt.
  Upstream      : None.
  Downstream    : 00_00_run_pipeline_DeFluorX.sh (Step 00).
───────────────────────────────────────────────────────────────────────────────

── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Conda Dependency: Assumes `conda` is available on the system PATH; will
     fail if using alternative package managers (mamba/micromamba) without alias.
  2. Internet Connectivity: Resyncing environments requires an active internet
     connection to reach Conda/PyPI repositories.
  3. Environment Naming: Hard-coded to export the 'PFAS' environment; requires
     manual modification if a different environment name is used.
  4. Platform Specificity: YAML exports include platform-specific builds;
     syncing between Linux and macOS/Windows may require `--from-history`.
───────────────────────────────────────────────────────────────────────────────
"""
# =============================================================================
# SECTION 1: IMPORTS & CONSOLE INFRASTRUCTURE
# =============================================================================

import subprocess
import sys
import argparse
from datetime import datetime
from pathlib import Path
# --- consolidated top-level imports (optional/heavy + Schrodinger stay function-local) ---
import importlib.util as _ilu
import shutil
import time as _time


'''
ConsoleColours is defined locally because this script runs BEFORE the PFAS
conda environment is guaranteed to exist - importing 00_02_Project_Utils is not
safe here. The canonical definition lives in 00_02; if it changes, sync this copy
manually (the drift assertion in main() guards the shared codes).
'''
class ConsoleColours:
    OKGREEN = "\033[92m"  # Green text designating success
    WARNING = "\033[93m"  # Yellow text designating caution
    FAIL    = "\033[91m"  # Red text designating failure
    OKBLUE  = "\033[94m"  # Blue text designating information
    MAGENTA = "\033[95m"  # Magenta text designating script banners
    BOLD    = "\033[1m"   # Bold text designed for headers
    ENDC    = "\033[0m"   # Reset colour formatting

SEPARATOR_HEAVY = "═" * 80


# =============================================================================
# SECTION 2: ENVIRONMENT MANAGEMENT
# =============================================================================

# -----------------------------------------------------------------------------
# Step 2.1: Environment Export
# -----------------------------------------------------------------------------

def export_environment():
    """Exports the current active Conda environment to PFAS.yml and requirements.txt.

    Overwrites the canonical files with the current host's exact versions; the export
    timestamp is written into each file's header.
    """
    timestamp = datetime.now().strftime("%d %B %Y, %H:%M")
    _yml_path, _req_path, _mode = Path("PFAS.yml"), Path("requirements.txt"), "--export"

    print(f"Exporting Conda environment to {_yml_path}...")
    try:
        '''
        Pin the export to the PFAS environment explicitly (-n PFAS) - never the
        ACTIVE env: running this from `base` would otherwise export base and
        silently overwrite PFAS.yml. --no-builds omits OS-specific build hashes.
        '''
        result = subprocess.run(
            ["conda", "env", "export", "-n", "PFAS", "--no-builds"],
            capture_output=True, text=True, check=True
        )
        '''
        Remove local environment prefix line if present
        Drop the local prefix line and the obsolete `dataclasses` backport
        (built into the Python 3.7+ stdlib; the pip 0.6 backport shadows it
        and breaks imports). _is_dropped filters both yaml `- dataclasses==`
        and bare `dataclasses==` forms.
        '''
        def _is_dropped(line: str) -> bool:
            s = line.strip()
            return s.startswith("prefix:") or s.lstrip("- ").startswith("dataclasses==")
        cleaned_lines = [line for line in result.stdout.splitlines() if not _is_dropped(line)]
        cleaned_yaml = "\n".join(cleaned_lines) + "\n"
        header = (
            f"# PFAS Conda Environment\n"
            f"# Exported : {timestamp}\n"
            f"# Script   : 00_03_Environment_DeFluorX.py {_mode}\n"
            f"#\n"
        )
        with open(_yml_path, "w") as f:
            f.write(header + cleaned_yaml)
        print(f"{ConsoleColours.OKGREEN}✔ Successfully created {_yml_path}{ConsoleColours.ENDC}")
    except subprocess.CalledProcessError as e:
        print(f"{ConsoleColours.FAIL}✘ Error exporting Conda environment: {e}{ConsoleColours.ENDC}")
        sys.exit(1)

    print(f"\nExporting pip requirements to {_req_path}...")
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "list", "--format=freeze"],
            capture_output=True, text=True, check=True
        )
        header = (
            f"# PFAS pip Requirements\n"
            f"# Exported : {timestamp}\n"
            f"# Script   : 00_03_Environment_DeFluorX.py {_mode}\n"
            f"#\n"
        )
        # Filter the obsolete `dataclasses` backport (stdlib since Python 3.7).
        _pip_lines = [ln for ln in result.stdout.splitlines()
                      if not ln.strip().startswith("dataclasses==")]
        with open(_req_path, "w") as f:
            f.write(header + "\n".join(_pip_lines) + "\n")
        print(f"{ConsoleColours.OKGREEN}✔ Successfully created {_req_path}{ConsoleColours.ENDC}")
    except subprocess.CalledProcessError as e:
        print(f"{ConsoleColours.FAIL}✘ Error exporting pip requirements: {e}{ConsoleColours.ENDC}")
        sys.exit(1)


# -----------------------------------------------------------------------------
# Step 2.2: Environment Installation
# -----------------------------------------------------------------------------

def install_environment(env_name: str):
    """Creates a fresh Conda environment from the exported PFAS.yml."""
    print(f"Installing fresh Conda environment '{env_name}' from PFAS.yml...")

    if not Path("PFAS.yml").exists():
        print(f"{ConsoleColours.FAIL}✘ Error: PFAS.yml not found. "
              f"Ensure the file is present in the working directory.{ConsoleColours.ENDC}")
        sys.exit(1)

    try:
        subprocess.run(
            ["conda", "env", "create", "-f", "PFAS.yml", "-n", env_name],
            check=True
        )
        print(f"\n{ConsoleColours.OKGREEN}✔ Successfully created environment '{env_name}'.{ConsoleColours.ENDC}")
        print("-" * 80)
        print(f"To activate your new environment, run:\n    conda activate {env_name}")
        print("-" * 80)
    except subprocess.CalledProcessError as e:
        print(f"{ConsoleColours.FAIL}✘ Error creating environment. "
              f"Verify your Conda installation: {e}{ConsoleColours.ENDC}")
        sys.exit(1)   # propagate failure so the pipeline runner halts instead of running on a broken env


# -----------------------------------------------------------------------------
# Step 2.3: Environment Verification
# -----------------------------------------------------------------------------

def verify_environment() -> None:
    """Verify the installed pipeline packages and HALT if a mandatory one is missing.

    This is a non-optional gate (00_00 runs it before Step 01), so a missing boltz/colabfold/MDAnalysis/
    rdkit/gemmi/torch must fail HERE with a non-zero exit - not print a red ✘ and let the run report PASS,
    only to die hours later inside Step 01/02. PyMOL/PLIP are the exception: Step 05 auto-installs them, so
    they are a ⚠, not a failure."""
    print("Verifying installed pipeline packages:")
    _missing: list[str] = []

    # Python
    py_ver = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    print(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} Python           : {py_ver}")

    def check_pkg(name: str, import_name: str | None = None) -> None:
        import_name = import_name or name
        try:
            mod = __import__(import_name)
            ver = getattr(mod, "__version__", "Available")
            print(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} {name:<16} : {ver}")
        except ImportError:
            print(f"  {ConsoleColours.FAIL}✘{ConsoleColours.ENDC} {name:<16} : Missing")
            _missing.append(name)

    check_pkg("boltz")
    check_pkg("colabfold")
    check_pkg("MDAnalysis")
    check_pkg("rdkit")
    check_pkg("gemmi")

    # Torch (CUDA)
    try:
        import torch
        cuda_avail = "CUDA available" if torch.cuda.is_available() else "CUDA NOT available"
        print(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} torch (CUDA)     : {torch.__version__} - {cuda_avail}")
    except ImportError:
        print(f"  {ConsoleColours.FAIL}✘{ConsoleColours.ENDC} torch (CUDA)     : Missing")
        _missing.append("torch")

    # PyMOL
    try:
        import pymol
        ver = getattr(pymol, "__version__", "Available")
        print(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} PyMOL            : {ver}")
    except ImportError:
        p = shutil.which("pymol")
        if p:
            print(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} PyMOL            : Available (binary located)")
        else:
            print(f"  {ConsoleColours.WARNING}⚠{ConsoleColours.ENDC} PyMOL            : Missing (will be auto-installed in Step 05)")

    # PLIP
    try:
        import plip
        ver = getattr(plip, "__version__", "Available")
        print(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} PLIP             : {ver}")
    except ImportError:
        p = shutil.which("plip")
        if p:
            print(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} PLIP             : Available (binary located)")
        else:
            print(f"  {ConsoleColours.WARNING}⚠{ConsoleColours.ENDC} PLIP             : Missing (will be auto-installed in Step 05)")

    if _missing:
        print(f"\n  {ConsoleColours.FAIL}✘ Environment check FAILED - mandatory package(s) missing: "
              f"{', '.join(_missing)}.{ConsoleColours.ENDC} Fix the conda env before running the pipeline.")
        sys.exit(1)


# =============================================================================
# SECTION 3: MAIN EXECUTION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Environment Manager for the DeFluorX Pipeline"
    )
    parser.add_argument(
        "--export",  action="store_true",
        help="Export the current environment to PFAS.yml and requirements.txt (overwrite)"
    )
    parser.add_argument(
        "--install", action="store_true",
        help="Install a fresh environment from PFAS.yml"
    )
    parser.add_argument(
        "--name", type=str, default="PFAS",
        help="Name of the environment to create (default: PFAS)"
    )

    args = parser.parse_args()

    _now = _time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"\n{SEPARATOR_HEAVY}", flush=True)
    print(
        f"  {ConsoleColours.MAGENTA}{ConsoleColours.BOLD}▶  00_03_Environment_DeFluorX.py"
        f"{ConsoleColours.ENDC}  │  DeFluorX Pipeline",
        flush=True,
    )
    print("  Conda Environment Export & Installation Manager", flush=True)
    print(f"  Started : {_now}", flush=True)
    print(f"{SEPARATOR_HEAVY}\n", flush=True)

    try:
        _spec = _ilu.spec_from_file_location("utils", Path(__file__).parent / "00_02_Project_Utils_DeFluorX.py")
        _u = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_u)
        for _code in ("OKGREEN", "WARNING", "FAIL", "OKBLUE", "MAGENTA", "BOLD", "ENDC"):
            assert getattr(ConsoleColours, _code) == getattr(_u.ConsoleColours, _code), \
                f"ConsoleColours.{_code} drift: sync this local copy in 00_03 to match canonical 00_02"
    except AssertionError as _drift:
        # A drift is a real (but non-fatal) maintenance issue - warn, don't crash the installer.
        print(f"{ConsoleColours.WARNING}⚠ {_drift}{ConsoleColours.ENDC}")
    except (FileNotFoundError, ModuleNotFoundError, AttributeError):
        pass

    if args.export:
        export_environment()
    elif args.install:
        install_environment(args.name)
    else:
        verify_environment()


if __name__ == "__main__":
    _t0 = _time.perf_counter()
    main()
    _el = _time.perf_counter() - _t0
    _h, _rem = divmod(int(_el), 3600)
    _m, _s   = divmod(_rem, 60)
    _fmt = (f"{_h}h {_m:02d}m {_s:02d}s" if _h else
            f"{_m}m {_s:02d}s"            if _m else
            f"{_s}s")
    print(f"\n{SEPARATOR_HEAVY}", flush=True)
    print(
        f"  {ConsoleColours.OKGREEN}✔  00_03_Environment_DeFluorX.py  -  Pipeline Phase Complete"
        f"  │  Total Elapsed: {_fmt}{ConsoleColours.ENDC}",
        flush=True,
    )
    print(f"{SEPARATOR_HEAVY}\n", flush=True)
