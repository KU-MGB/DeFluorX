#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 00  |  Environment Management
===============================================================================
Exports the active Conda environment to 'PFAS.yml' and 'requirements.txt'
for reproducibility. Provides an automated installation routine to
re-synchronise environments across compute nodes.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 June 2026 <─────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 00_01_Environment_Installation_FAcDs.py
  Role          : Infrastructure — environment export and sync.
  Imports from  : None (standalone sys/os/subprocess).
  Reads         : Active conda environment.
  Writes        : PFAS.yml, requirements.txt.
  Upstream      : None.
  Downstream    : 00_00_run_pipeline_FAcDs.sh (Step 00).
───────────────────────────────────────────────────────────────────────────────

# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
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
# ===============================================================================
# SECTION 1: IMPORTS & CONSOLE INFRASTRUCTURE
# ===============================================================================

import subprocess
import sys
import argparse
from datetime import datetime
from pathlib import Path


# ConsoleColours is defined locally here because this script runs BEFORE the
# PFAS conda environment is guaranteed to exist — importing 00_02_Project_Utils
# is not safe at this point. If the canonical definition in 00_02 is updated,
# this local copy must be synchronised manually.
class ConsoleColours:
    OKGREEN = '\033[92m'  # Green text designating success
    WARNING = '\033[93m'  # Yellow text designating caution
    FAIL    = '\033[91m'  # Red text designating failure
    OKBLUE  = '\033[94m'  # Blue text designating information
    BOLD    = '\033[1m'   # Bold text designed for headers
    ENDC    = '\033[0m'   # Reset colour formatting

SEPARATOR_HEAVY = '═' * 80


# ===============================================================================
# SECTION 2: ENVIRONMENT MANAGEMENT
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 2.1: Environment Export
# -------------------------------------------------------------------------------

def export_environment():
    """Exports the current active Conda environment to .yml and .txt files."""
    timestamp = datetime.now().strftime("%d %B %Y, %H:%M")

    print("Exporting Conda environment to PFAS.yml...")
    try:
        # --no-builds ensures cross-platform compatibility by omitting
        # OS-specific build hashes from the exported specification.
        result = subprocess.run(
            ["conda", "env", "export", "--no-builds"],
            capture_output=True, text=True, check=True
        )
        # Remove local environment prefix line if present
        cleaned_lines = [line for line in result.stdout.splitlines() if not line.strip().startswith("prefix:")]
        cleaned_yaml = "\n".join(cleaned_lines) + "\n"
        header = (
            f"# PFAS Conda Environment\n"
            f"# Exported : {timestamp}\n"
            f"# Script   : 00_01_Environment_Installation_FAcDs.py --export\n"
            f"#\n"
        )
        with open("PFAS.yml", "w") as f:
            f.write(header + cleaned_yaml)
        print(f"{ConsoleColours.OKGREEN}✔ Successfully created PFAS.yml{ConsoleColours.ENDC}")
    except subprocess.CalledProcessError as e:
        print(f"{ConsoleColours.FAIL}✘ Error exporting Conda environment: {e}{ConsoleColours.ENDC}")
        sys.exit(1)

    print("\nExporting pip requirements to requirements.txt...")
    try:
        result = subprocess.run(
            [sys.executable, "-m", "pip", "list", "--format=freeze"],
            capture_output=True, text=True, check=True
        )
        header = (
            f"# PFAS pip Requirements\n"
            f"# Exported : {timestamp}\n"
            f"# Script   : 00_01_Environment_Installation_FAcDs.py --export\n"
            f"#\n"
        )
        with open("requirements.txt", "w") as f:
            f.write(header + result.stdout)
        print(f"{ConsoleColours.OKGREEN}✔ Successfully created requirements.txt{ConsoleColours.ENDC}")
    except subprocess.CalledProcessError as e:
        print(f"{ConsoleColours.FAIL}✘ Error exporting pip requirements: {e}{ConsoleColours.ENDC}")
        sys.exit(1)


# -------------------------------------------------------------------------------
# Step 2.2: Environment Installation
# -------------------------------------------------------------------------------

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


# ===============================================================================
# SECTION 3: MAIN EXECUTION
# ===============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Environment Manager for the PFAS Pipeline"
    )
    parser.add_argument(
        "--export",  action="store_true",
        help="Export the current environment to PFAS.yml and requirements.txt"
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

    import time as _t
    _now = _t.strftime("%Y-%m-%d %H:%M:%S")
    print(f"\n{SEPARATOR_HEAVY}", flush=True)
    print(
        f"  \033[95m\033[1m▶  00_01_Environment_Installation_FAcDs.py"
        f"\033[0m  │  FAcDs Pipeline",
        flush=True,
    )
    print(f"  Conda Environment Export & Installation Manager", flush=True)
    print(f"  Started : {_now}", flush=True)
    print(f"{SEPARATOR_HEAVY}\n", flush=True)

    try:
        import importlib.util as _ilu
        _spec = _ilu.spec_from_file_location("utils", Path(__file__).parent / "00_02_Project_Utils_FAcDs.py")
        _u = _ilu.module_from_spec(_spec); _spec.loader.exec_module(_u)
        assert ConsoleColours.OKGREEN == _u.ConsoleColours.OKGREEN, \
            "ConsoleColours.OKGREEN drift: update 00_01 to match 00_02"
        assert ConsoleColours.FAIL == _u.ConsoleColours.FAIL, \
            "ConsoleColours.FAIL drift: update 00_01 to match 00_02"
    except (FileNotFoundError, ModuleNotFoundError, AttributeError):
        pass

    if args.export:
        export_environment()
    elif args.install:
        install_environment(args.name)
    else:
        print(f"{ConsoleColours.FAIL}Please specify an action: --export or --install.{ConsoleColours.ENDC}\n")
        parser.print_help()


if __name__ == "__main__":
    import time as _time
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
        f"  {ConsoleColours.OKGREEN}✔  00_01_Environment_Installation_FAcDs.py  —  Pipeline Phase Complete"
        f"  │  Total Elapsed: {_fmt}{ConsoleColours.ENDC}",
        flush=True,
    )
    print(f"{SEPARATOR_HEAVY}\n", flush=True)
