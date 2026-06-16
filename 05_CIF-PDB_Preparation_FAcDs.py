#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 05  |  PDB Generation & Protein Preparation
===============================================================================
CIF → PDB conversion (Gemmi) and Schrödinger PrepWizard preparation for all
best-complex structures selected by the Production script (Step 02).

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 June 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 05_CIF-PDB_Preparation_FAcDs.py
  Role          : "Builder" — converts Boltz-2 CIF outputs to analysis-ready PDB.
  Imports from  : 00_02_Project_Config_FAcDs.py  (CFG — pH values)
                  00_03_Project_Utils_FAcDs.py   (ConsoleColours, setup_logging,
                                            console_title, console_info,
                                            console_separator)
  Reads         : <Run>/2_Best_Complexes_CIFs/*.cif
                  <Run>/1_Boltz2_Production/7_Boltz2_FAcDs_Ranked_*.csv
  Writes        : <Run>/5_PDB_Generation_Preparation/1_Converted_Raw_PDB/*.pdb
                  <Run>/5_PDB_Generation_Preparation/2_Prepared_PDBs/*.pdb
                  <Run>/5_PDB_Generation_Preparation/3_PDB_prep_master.log
  Upstream      : 02_Production_FAcDs.py  → writes Best_Complexes_CIFs and ranked CSV
  Downstream    : 06_Top-N_Extraction_FAcDs.py  → reads prepared PDBs
                  08_MD_Thermodynamics_QMMM_Engine_FAcDs.py → reads prepared PDBs for MD/QM-MM
───────────────────────────────────────────────────────────────────────────────

# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Schrödinger Dependency: Hard requirement for Schrödinger's `prepwizard`
     binary; will fall back to raw (unprepared) PDBs if missing, which may compromise downstream MD/QM-MM quality.
  2. Ligand Naming: CIF-to-PDB conversion using Gemmi may struggle with highly
     non-standard ligands; relies on Chain L isolation for downstream detection.
  3. Parallel Overhead: Spawns multiple `prepwizard` instances; requires
     sufficient license tokens and CPU cores to avoid job starvation.
  4. Header Sensitivity: Relies on the SOURCE_CIF tag in PDB headers for caching;
     manual header modification will trigger redundant re-processing.
───────────────────────────────────────────────────────────────────────────────

Usage:
    conda activate PFAS
    python 05_CIF-PDB_Preparation_FAcDs.py Boltz-2_Run_20260309T085406Z

── Key features ───────────────────────────────────────────────────────────────
  • Smart Validation: SOURCE_CIF tag in PDB header detects stale files from
    previous runs and triggers automatic regeneration.
  • Stale-file Detection: re-prepares if Raw PDB is newer than Prepared PDB.
  • Gemmi Conversion: robust CIF → PDB; moves non-protein residues to Chain L.
  • Smart Ligand Management: metals (Zn, Mg) and modified residues (MSE) stay
    with the protein chain to preserve topology for Maestro.
  • PrepWizard: fills side chains, PropKa pH 8.0, Epik pH 8.0, 0.3 Å RMSD
    restrained minimisation, disulfide detection.
───────────────────────────────────────────────────────────────────────────────
"""

# ===============================================================================
# SECTION 1: SYSTEM CONFIGURATION & IMPORTS
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 1.1: Standard Library Imports
# -------------------------------------------------------------------------------
import os
import sys

# CPU usage cap (total cores − 2; mirrors CFG.PREP_CPU_RESERVE). Reserve 2 cores
# for OS/desktop stability by limiting the thread-pool maths libraries (BLAS /
# MKL / OpenMP / NumExpr). Must precede numpy/pandas import to take effect;
# setdefault() preserves any value exported by the caller or pipeline runner.
_CPU_CAP = str(max(1, (os.cpu_count() or 4) - 2))
for _tv in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_tv, _CPU_CAP)

import shutil
import subprocess
import tempfile
import time
import argparse
import re
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# -------------------------------------------------------------------------------
# Step 1.2: Scientific Stack
# -------------------------------------------------------------------------------
import gemmi

# -------------------------------------------------------------------------------
# Step 1.3: Pipeline modules (00_02 (config & utils)) via importlib
# Filenames begin with digits and cannot be imported with standard `import`.
# -------------------------------------------------------------------------------
import importlib.util as _ilu

def _load_module(name: str, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod  = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_REPO_DIR   = Path(__file__).resolve().parent
_cfg_mod    = _load_module("ProjectConfig", _REPO_DIR / "00_02_Project_Config_FAcDs.py")
_utils_mod  = _load_module("ProjectUtils",  _REPO_DIR / "00_03_Project_Utils_FAcDs.py")

CFG             = _cfg_mod.CFG()
ConsoleColours     = _utils_mod.ConsoleColours
SEPARATOR_HEAVY    = _utils_mod.SEPARATOR_HEAVY
SEPARATOR_LIGHT    = _utils_mod.SEPARATOR_LIGHT

# -------------------------------------------------------------------------------
# Step 1.4: Global Constants & Paths
# -------------------------------------------------------------------------------

# Ensure SCHRODINGER is set in the process environment so child processes
# (PrepWizard subprocess calls) inherit it without requiring a prior `export`.
if "SCHRODINGER" not in os.environ:
    os.environ["SCHRODINGER"] = "/opt/schrodinger"

SCHRODINGER_PATH = Path(os.environ["SCHRODINGER"])
PREPWIZARD_BIN   = SCHRODINGER_PATH / "utilities" / "prepwizard"

# Concurrency: reserve 2 cores for OS/desktop stability

MAX_PREP_JOBS = CFG.GLOBAL_MAX_WORKERS

# Formatting
SEPARATOR          = "-" * 80
BATCH_SIZE         = CFG.PROC_BATCH_SIZE
HEADER_CHECK_LINES = CFG.PROC_HEADER_CHECK_LINES
logger             = None

STANDARD_AA        = CFG.PREP_STANDARD_AA
PROTEIN_ASSOCIATED = CFG.PREP_PROTEIN_ASSOCIATED


# ===============================================================================
# SECTION 2: LOGGING INFRASTRUCTURE
# ===============================================================================
# Logging and console functions are provided by 00_03_Project_Utils.
# Script-level wrappers capture the module-global `logger` so existing call
# sites require no modification.

def setup_logging(prep_base_dir: Path) -> Path:
    """Initialises the preparation log via the shared utility."""
    global logger
    prep_base_dir.mkdir(parents=True, exist_ok=True)
    master_log = prep_base_dir / "3_PDB_prep_master.log"
    logger = _utils_mod.setup_logging(master_log, logger_name="pdb_prep", mode="a", timestamp=True)
    return master_log

def console_title(msg: str) -> None:
    _utils_mod.console_title(msg, logger)

def console_info(msg: str) -> None:
    _utils_mod.console_info(msg, logger)

def console_separator() -> None:
    _utils_mod.console_separator(logger, heavy=True)


# ===============================================================================
# SECTION 3: CORE LOGIC & HELPERS
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 3.1: File Indexing & Metadata
# -------------------------------------------------------------------------------
def index_existing_files(directory: Path, suffix: str) -> dict:
    """Creates a fast lookup dictionary {JobIndex: Path} for existing files.
    Special handling for Controls (ID 0000000) to avoid collisions.
    """
    index = {}
    if not directory.exists(): return index
    
    for f in directory.glob(f"*{suffix}"):
        try:
            parts = f.name.split("_")
            if parts and parts[0].isdigit():
                idx = parts[0]
                # If it's a control ID, it is not indexed for renaming
                # because multiple controls share this ID.
                if idx == "0000000":
                    continue
                index[idx] = f
        except Exception as e:
            if logger: logger.debug(f"Error parsing file name {f.name}: {e}")
    return index

def collect_best_cifs(best_cifs_dir: Path) -> list:
    """
    Scans 2_Best_Complexes_CIFs for all CIF files and derives job names.
    Returns list of (job_name, cif_path) tuples.
    Job name is CIF stem with '_model_N' suffix stripped.
    Ensures uniqueness: if multiple models exist for the same job name,
    keeps the most recently modified one.
    """
    latest_cifs = {}
    for cif_path in best_cifs_dir.glob("*.cif"):
        job_name = re.sub(r'_model_\d+$', '', cif_path.stem)
        mtime = cif_path.stat().st_mtime
        if job_name not in latest_cifs or mtime > latest_cifs[job_name][1]:
            latest_cifs[job_name] = (cif_path, mtime)
    
    # Return sorted list for deterministic order
    return [(jn, path) for jn, (path, _) in sorted(latest_cifs.items())]

def load_rank_map(prod_dir: Path) -> dict:
    """Loads the master CSV to annotate PDB headers with a score/rank.
    Prefers the FAcDs Ranked CSV (7_Boltz2_FAcDs_Ranked_*), then any ranked CSV,
    then falls back to any master CSV.
    Uses Scientific_Rank if present, otherwise confidence_score (rounded to 4dp).
    """
    rank_csvs = sorted(prod_dir.glob("7_Boltz2_FAcDs_Ranked_*.csv"))
    if not rank_csvs:
        rank_csvs = sorted(prod_dir.glob("*_Ranked_*.csv"))
    if not rank_csvs:
        console_info("Warning: No Ranked CSV found in production directory.")
        return {}

    latest_csv = sorted(rank_csvs)[-1]
    console_info(f"Loaded Rank File : {latest_csv.name}")

    try:
        df = pd.read_csv(latest_csv, low_memory=False)
        if 'job_name' not in df.columns:
            return {}
        if 'Scientific_Rank' in df.columns:
            return dict(zip(df['job_name'], df['Scientific_Rank']))
        if 'confidence_score' in df.columns:
            # Use rounded confidence score as a rank proxy
            return {jn: f"{sc:.4f}" for jn, sc in zip(df['job_name'], df['confidence_score'])}
    except Exception as e:
        if logger: logger.debug(f"Error loading Rank Map: {e}")
    return {}

# -------------------------------------------------------------------------------
# Step 3.2: PDB Manipulation (Header & Source Tracking)
# -------------------------------------------------------------------------------
def get_source_tag(pdb_path: Path) -> str:
    """Reads the REMARK 999 SOURCE_CIF tag from a PDB file."""
    try:
        with open(pdb_path, 'r') as f:
            for _ in range(HEADER_CHECK_LINES): # Check first HEADER_CHECK_LINES lines
                line = f.readline()
                if not line: break
                if "REMARK 999 SOURCE_CIF:" in line:
                    return line.split(":", 1)[1].strip()
    except Exception as e:
        logger.debug(f"Header tag not found: {e}") # Not critical if header is missing, return empty string
    return ""

def update_pdb_header(pdb_path: Path, job_name: str, rank: str, source_cif: str):
    """
    Injects metadata into the PDB Header.
    Crucially, adds 'SOURCE_CIF' to track which model generated this PDB.
    """
    try:
        lines = pdb_path.read_text().splitlines()
        # Remove old custom headers to avoid duplication
        lines = [l for l in lines if not l.startswith(("HEADER", "TITLE", "REMARK 999"))]
        
        header_block = [
            f"HEADER    {job_name[:40].ljust(40)}",
            f"TITLE     Scientific Rank: {rank} | FAcD Pipeline",
            f"REMARK 999 SOURCE_CIF: {source_cif}"
        ]
        
        new_content = "\n".join(header_block + lines) + "\n"
        pdb_path.write_text(new_content)
        return True
    except Exception as e:
        if logger: logger.warning(f"Could not update header for {job_name}: {e}")
        return False

# -------------------------------------------------------------------------------
# Step 3.3: Conversion Engine (Gemmi)
# -------------------------------------------------------------------------------
def cif_to_pdb_gemmi(cif_path: Path, pdb_path: Path, job_name: str, rank: str):
    """
    Converts mmCIF to PDB using Gemmi.
    - Moves Ligands (Non-Standard Residues) to Chain 'L' for Maestro compatibility.
    - Preserves Metals/Modified AAs with protein using PROTEIN_ASSOCIATED whitelist.
    - Preserves atom names and residue numbering.
    """
    try:
        if cif_path.stat().st_size == 0: return False, "Empty CIF File"

        doc = gemmi.cif.read_file(str(cif_path))
        block = doc.sole_block()
        st = gemmi.make_structure_from_block(block)
        
        if len(st) > 0:
            model = st[0] 
            lig_chain = gemmi.Chain(CFG.PREP_LIGAND_CHAIN)
            has_ligands = False

            # Iterate chains to find and move ligands
            for chain in model:
                if chain.name == CFG.PREP_LIGAND_CHAIN: continue
                indices_to_move = []
                
                for i, res in enumerate(chain):
                    # Move to Chain L only if it's NOT standard, NOT water, and NOT protein-associated
                    if res.name not in STANDARD_AA and res.name != "HOH" and res.name not in PROTEIN_ASSOCIATED:
                        indices_to_move.append(i)
                        
                # Move from back to front to avoid index shift issues
                for i in reversed(indices_to_move):
                    res = chain[i]
                    lig_chain.add_residue(res) 
                    del chain[i]
                    has_ligands = True
                    
            if has_ligands: model.add_chain(lig_chain)
        
        st.setup_entities() 
        st.write_pdb(str(pdb_path))
        
        # Inject validation metadata
        update_pdb_header(pdb_path, job_name, rank, cif_path.name)
        
        return True, "Success"
    except Exception as e:
        return False, f"Gemmi Error: {str(e)}"

# -------------------------------------------------------------------------------
# Step 3.4: Preparation Engine (Schrödinger PrepWizard)
# -------------------------------------------------------------------------------
def run_prepwizard(raw_pdb: Path, final_dest: Path):
    """
    Wraps Schrödinger's prepwizard for headless protein-ligand preparation.
    - Adds hydrogens (default behaviour; -nohtreat NOT used)
    - Fills missing side chains with Prime
    - PropKa protonation at pH 8.0 (FAcD physiological context)
    - Epik protonation of PFAS ligands at pH 8.0
    - Restrained minimisation (0.3 Å RMSD) to resolve clashes
    - Runs inline via -NOJOBID (no Schrödinger job-server dependency)
    - Uses /tmp for work dir — nothing written to the output folder
    """
    job_name = raw_pdb.stem.replace("_RAW", "")
    if not PREPWIZARD_BIN.exists(): return False

    job_work_dir = Path(tempfile.mkdtemp(prefix=f"prepwiz_{job_name}_"))
    try:
        input_pdb_local = job_work_dir / f"{job_name}_RAW.pdb"
        output_pdb_name = f"{job_name}_Prepared.pdb"
        shutil.copy2(raw_pdb, input_pdb_local)

        # Headless execution environment — strip Python-env overrides that confuse
        # Schrödinger's bundled Python interpreter.
        env = os.environ.copy()
        env["SCHRODINGER"]     = str(SCHRODINGER_PATH)
        env["QT_QPA_PLATFORM"] = "offscreen"
        for k in ("PYTHONPATH", "PYTHONHOME", "LD_LIBRARY_PATH"):
            env.pop(k, None)

        cmd = [
            str(PREPWIZARD_BIN),
            "-fillsidechains",       # Fill truncated side chains (common in AI models)
            "-propka_pH", str(CFG.PREPWIZARD_PROPKA_PH),  # protein protonation pH
            "-epik_pH",   str(CFG.PREPWIZARD_EPIK_PH),    # PFAS ligand protonation via Epik
            "-r",         str(CFG.PREPWIZARD_RMSD_RESTRAIN),  # RMSD-restrained minimisation
            "-disulfides",           # Detect and bond proximal Cys pairs
            "-NOJOBID",              # Run inline — no Schrödinger Job Control layer
            input_pdb_local.name,
            output_pdb_name,
        ]

        log_file = job_work_dir / f"{job_name}_prepwizard.log"

        with open(log_file, "w") as log_f:
            result = subprocess.run(
                cmd,
                cwd=job_work_dir,
                stdout=log_f,
                stderr=subprocess.STDOUT,
                env=env,
                check=False,
            )

        if result.returncode != 0:
            if logger:
                logger.error(f"PrepWizard non-zero exit ({result.returncode}) for {job_name}. "
                             f"Log: {log_file}")
            return False

        prepared_path = job_work_dir / output_pdb_name
        if prepared_path.exists() and prepared_path.stat().st_size > 0:
            shutil.copy2(prepared_path, final_dest)
            return True

        if logger:
            logger.error(f"PrepWizard returned 0 but output missing for {job_name}")
    except Exception as e:
        if logger: logger.error(f"PrepWizard execution failed for {job_name}: {e}")
    finally:
        shutil.rmtree(job_work_dir, ignore_errors=True)
    return False


# ===============================================================================
# SECTION 4: TASK EXECUTION & VALIDATION LOGIC
# ===============================================================================

def generate_raw_step(job_name: str, best_cif: Path, dir_raw: Path, rank: str, raw_index: dict):
    """
    Task 1: Generate Raw PDB from CIF.
    CIF is passed directly from 2_Best_Complexes_CIFs (pre-collected by Production script).
    Logic: Checks if existing PDB matches the source CIF. If not, regenerates.
    """
    raw_pdb_path = dir_raw / f"{job_name}_RAW.pdb"

    # --- Validation Logic ---
    if raw_pdb_path.exists() and raw_pdb_path.stat().st_size > 0:
        # Check if the file on disk matches the current best CIF name
        recorded_source = get_source_tag(raw_pdb_path)
        is_valid = (recorded_source == best_cif.name)
        
        if is_valid:
            # Skip this job, but signal 'Skipped' so main loop knows it exists
            return {"job": job_name, "status": "Skipped", "rank": rank}
        else:
            # Source Mismatch (New model selected in Production) -> Force Regen
            try: raw_pdb_path.unlink()
            except Exception as e:
                if logger: logger.debug(f"Failed to unlink stale raw PDB {raw_pdb_path.name}: {e}")
            # Fall through to generation below

    # Check if renaming an old file from the index is required
    job_index = job_name.split("_")[0]
    if job_index in raw_index:
        old_file = raw_index[job_index]
        if old_file.exists() and old_file.name != raw_pdb_path.name:
            try:
                old_file.rename(raw_pdb_path)
                # Check validity after rename
                recorded_source = get_source_tag(raw_pdb_path)
                if recorded_source == best_cif.name:
                    return {"job": job_name, "status": "Renamed", "rank": rank}
                else:
                    raw_pdb_path.unlink() # Renamed file was stale
            except Exception as e:
                if logger: logger.debug(f"Failed to rename raw PDB {old_file.name}: {e}") 

    # --- Generation Logic ---
    ok, error_msg = cif_to_pdb_gemmi(best_cif, raw_pdb_path, job_name, rank)
    if not ok: return {"job": job_name, "status": "Conversion_Failed", "reason": error_msg, "rank": rank}
    
    # Return Regened status if it existed before but was invalid, else Success
    return {"job": job_name, "status": "Success", "rank": rank}

def check_prep_needed(job_name: str, dir_raw: Path, dir_prep_clean: Path):
    """
    Pre-check for Phase 2: Returns True if preparation is needed, False if skipped.
    """
    final_prep_path = dir_prep_clean / f"{job_name}_Prepared.pdb"
    raw_pdb_path = dir_raw / f"{job_name}_RAW.pdb"
    
    if not raw_pdb_path.exists():
        return True # Can't validate, so let it fail/handle in step

    if final_prep_path.exists() and final_prep_path.stat().st_size > 0:
        raw_source = get_source_tag(raw_pdb_path)
        prep_source = get_source_tag(final_prep_path)
        
        # Check Source Tag Consistency
        if raw_source and prep_source and raw_source != prep_source:
            try: final_prep_path.unlink() # Delete bad file
            except Exception: pass
            return True # Needed
            
        # Check Timestamp (Backup)
        elif raw_pdb_path.stat().st_mtime > final_prep_path.stat().st_mtime:
            try: final_prep_path.unlink()
            except Exception: pass
            return True # Needed
            
        return False # Skipped (Valid)
        
    return True # Needed (Missing)

def preparation_step(job_name: str, dir_raw: Path, dir_prep_clean: Path, rank: str, prep_index: dict):
    """
    Task 2: Protein Preparation.
    Logic: Checks consistency between Raw and Prep PDBs (Timestamp + Source Tag).
    PrepWizard runs in a /tmp work dir that is deleted automatically on completion.
    """
    final_prep_path = dir_prep_clean / f"{job_name}_Prepared.pdb"
    raw_pdb_path = dir_raw / f"{job_name}_RAW.pdb"

    if not raw_pdb_path.exists():
        return {"job": job_name, "status": "Missing_Raw", "rank": rank}

    # Get source info from Raw
    raw_source = get_source_tag(raw_pdb_path)

    # Rename Check (Handle renames before execution)
    job_index = job_name.split("_")[0]
    if job_index in prep_index:
        old_file = prep_index[job_index]
        if old_file.exists() and old_file.name != final_prep_path.name:
            try:
                old_file.rename(final_prep_path)
            except Exception as e:
                if logger: logger.debug(f"Failed to rename prepared PDB {old_file.name}: {e}")

    # --- Execution Logic ---
    if not PREPWIZARD_BIN.exists():
        if logger: logger.warning(f"PrepWizard binary not found at {PREPWIZARD_BIN}")
        return {"job": job_name, "status": "No_PrepWizard", "rank": rank}

    prep_success = run_prepwizard(raw_pdb_path, final_prep_path)
    if prep_success:
        update_pdb_header(final_prep_path, job_name, rank, raw_source or "Unknown")

    status = "Success" if prep_success else "Prep_Failed"
    return {"job": job_name, "status": status, "rank": rank}


# ===============================================================================
# SECTION 5: MAIN EXECUTION
# ===============================================================================

def _check_residue_identity_guard(prepared_pdb_path: Path, job_name: str, cfg) -> dict:
    """Parse prepared PDB; verify catalytic residues at expected positions. Write index_offset.json."""
    import json

    nuc_ref  = cfg.DREAM_TEAM_REFS["Nuc"]   # 110
    acid_ref = cfg.DREAM_TEAM_REFS["Acid"]  # 134
    base_ref = cfg.DREAM_TEAM_REFS["Base"]  # 277

    ASP_TYPES = {"ASP", "ASH"}
    HIS_TYPES = {"HIS", "HIE", "HID", "HIP"}

    resnum_to_resname = {}
    try:
        with open(prepared_pdb_path) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    resname = line[17:20].strip()
                    resnum  = int(line[22:26].strip())
                    resnum_to_resname[resnum] = resname
    except Exception:
        return {"offset": 0}

    def find_type_near(ref, expected_types, window=CFG.SMART_LOCK_RESNUM_WINDOW):
        if resnum_to_resname.get(ref) in expected_types:
            return ref, 0
        for delta in range(1, window + 1):
            for sign in (-1, 1):
                r = ref + sign * delta
                if resnum_to_resname.get(r) in expected_types:
                    return r, sign * delta
        return None, None

    nuc_found,  nuc_offset  = find_type_near(nuc_ref,  ASP_TYPES)
    acid_found, acid_offset = find_type_near(acid_ref, ASP_TYPES)
    base_found, base_offset = find_type_near(base_ref, HIS_TYPES)

    # Use the most common non-zero offset (consensus across Nuc/Acid/Base);
    # fall back to 0 when all three agree on the reference numbering.
    from collections import Counter as _Counter
    offsets = [o for o in (nuc_offset, acid_offset, base_offset) if o is not None and o != 0]
    offset = _Counter(offsets).most_common(1)[0][0] if offsets else 0

    record = {
        "offset":     offset,
        "nuc_found":  nuc_found,
        "acid_found": acid_found,
        "base_found": base_found,
    }

    out_json = prepared_pdb_path.parent / f"{job_name}_index_offset.json"
    with open(out_json, "w") as fh:
        json.dump(record, fh, indent=2)

    return record


def _print_run_delta_table(n_raw_cached: int, n_raw_new: int, n_prep_cached: int, n_prep_new: int):
    """Print a 4-row delta table summarising what this run generated vs. cached."""
    _w = 45
    rows = [
        ("RAW PDBs already on disk (cached)",      n_raw_cached),
        ("RAW PDBs generated this run",             n_raw_new),
        ("Prepared PDBs already on disk (cached)",  n_prep_cached),
        ("Prepared PDBs generated this run",        n_prep_new),
    ]
    console_info(f"┌{'─' * _w}┬{'─' * 10}┐")
    for label, value in rows:
        console_info(f"│  {label:<{_w - 2}}│{value:>9,} │")
    console_info(f"└{'─' * _w}┴{'─' * 10}┘")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("run_folder_name", help="Run Folder Name (e.g. Boltz-2_Run_...)")
    parser.add_argument("--quick", action="store_true", help="Quick resume: skip deep validation of existing files")
    args = parser.parse_args()

    root = Path.cwd()
    run_path = root / args.run_folder_name
    if not run_path.exists(): sys.exit(f"Run path missing: {run_path}")

    # Paths
    analysis_dir   = run_path / "5_PDB_Generation_Preparation"
    dir_raw        = analysis_dir / "1_Converted_Raw_PDB"
    dir_prep_clean = analysis_dir / "2_Prepared_PDBs"
    prod_dir       = run_path / "1_Boltz2_Production"

    # Init output directories and log (log goes directly in analysis_dir)
    for d in [dir_raw, dir_prep_clean]: d.mkdir(parents=True, exist_ok=True)
    setup_logging(analysis_dir)

    best_cifs_dir = run_path / "2_Best_Complexes_CIFs"
    if not best_cifs_dir.exists(): sys.exit(f"Critical: 2_Best_Complexes_CIFs directory not found in {run_path}")
    jobs = collect_best_cifs(best_cifs_dir)

    _utils_mod.print_script_banner(
        "05_CIF-PDB_Preparation_FAcDs.py",
        "CIF → PDB Conversion  ·  PrepWizard Refinement  ·  Chain Assignment",
    )
    console_info(f"  Run Name : {args.run_folder_name}")
    console_info(f"Found {len(jobs)} Best Complex CIFs in 2_Best_Complexes_CIFs")
    if args.quick:
        console_info(f"{ConsoleColours.OKBLUE}Quick Resume Mode: Deep validation disabled.{ConsoleColours.ENDC}")
    
    rank_map = load_rank_map(prod_dir)

    print("")
    console_info(f"Outputs      : {analysis_dir}")
    console_info(f"Schrödinger  : {SCHRODINGER_PATH}")
    console_info(f"PrepWizard   : {'Found' if PREPWIZARD_BIN.exists() else 'NOT FOUND — Step 5.2 will be skipped'}")
    console_info(f"Parallel Jobs: {MAX_PREP_JOBS} (Globally managed by CFG)")
    console_separator()
    
    # -------------------------------------------------------------------------------
    # Step 5.1: Raw PDB Generation (With Source Validation)
    # -------------------------------------------------------------------------------
    print(SEPARATOR_LIGHT, flush=True)
    console_info("Generating Raw PDBs")
    
    console_info("Indexing existing Raw PDBs...")
    raw_index = index_existing_files(dir_raw, "_RAW.pdb")
    # Quick index of full filenames for exact matches
    existing_raw_files = {f.name for f in dir_raw.glob("*_RAW.pdb")}
    
    t0 = time.time()
    
    # Re-scan for missing or invalid files
    jobs_to_generate = []
    final_valid_raw_names = []
    
    for job_name, best_cif in jobs:
        raw_name = f"{job_name}_RAW.pdb"
        raw_pdb_path = dir_raw / raw_name

        needs_run = True
        try:
            _raw_sz = raw_pdb_path.stat().st_size if raw_pdb_path.exists() else 0
        except FileNotFoundError:
            _raw_sz = 0
        if raw_name in existing_raw_files and _raw_sz > 0:
            if args.quick:
                needs_run = False
            elif get_source_tag(raw_pdb_path) == best_cif.name:
                needs_run = False

        if needs_run:
            jobs_to_generate.append((job_name, best_cif))
        else:
            final_valid_raw_names.append(job_name)
            
    n_raw_cached = len(final_valid_raw_names)
    n_raw_new    = len(jobs_to_generate)
    console_info(f"Existing Valid Raw PDBs : {n_raw_cached}")
    console_info(f"Jobs to Generate/Update : {n_raw_new}")

    if jobs_to_generate:
        total_gen = len(jobs_to_generate)
        console_info(f"Queuing {total_gen} generation tasks (batch mode)...")
        count = 0
        with ThreadPoolExecutor(max_workers=MAX_PREP_JOBS) as executor:
            for batch_start in range(0, total_gen, BATCH_SIZE):
                batch = jobs_to_generate[batch_start: batch_start + BATCH_SIZE]
                futures = {
                    executor.submit(
                        generate_raw_step,
                        job_name, best_cif, dir_raw, str(rank_map.get(job_name, "N/A")), raw_index
                    ): job_name for job_name, best_cif in batch
                }
                for future in as_completed(futures):
                    res = future.result()
                    count += 1
                    if res['status'] in ("Success", "Renamed", "Skipped"):
                        final_valid_raw_names.append(res['job'])
                    if count % 100 == 0 or count == total_gen:
                        if sys.stdout.isatty():
                            print(f"\rGenerated: {count}/{total_gen}", end="", flush=True)
                        elif count == total_gen:
                            print(f"Generated: {count}/{total_gen}", flush=True)
        if sys.stdout.isatty():
            print("")

    console_info(f"Raw PDB generation stage completed in {time.time() - t0:.1f}s")
    console_separator()

    # -------------------------------------------------------------------------------
    # Step 5.2: Protein Preparation via PrepWizard (With Stale Check)
    # -------------------------------------------------------------------------------
    print(SEPARATOR_LIGHT, flush=True)
    console_info("Protein Preparation")

    console_info("Indexing existing Prepared PDBs...")
    prep_index = index_existing_files(dir_prep_clean, "_Prepared.pdb")
    existing_prep_files = {f.name for f in dir_prep_clean.glob("*_Prepared.pdb")}

    # Filter Prep Jobs
    prep_jobs_to_run = []
    prep_already_done_count = 0
    
    for job_name in final_valid_raw_names:
        prep_name = f"{job_name}_Prepared.pdb"
        if args.quick and prep_name in existing_prep_files:
            prep_already_done_count += 1
        elif check_prep_needed(job_name, dir_raw, dir_prep_clean):
            prep_jobs_to_run.append(job_name)
        else:
            prep_already_done_count += 1
            
    console_info(f"Existing Valid Prepared : {prep_already_done_count}")
    console_info(f"Jobs to Prepare         : {len(prep_jobs_to_run)}")

    count = 0
    total_prep = len(prep_jobs_to_run)
    
    if total_prep > 0:
        console_info(f"Queuing {total_prep} preparation tasks (batch mode)...")
        with ThreadPoolExecutor(max_workers=MAX_PREP_JOBS) as executor:
            for batch_start in range(0, total_prep, BATCH_SIZE):
                batch = prep_jobs_to_run[batch_start: batch_start + BATCH_SIZE]
                futures = {
                    executor.submit(preparation_step, job_name, dir_raw, dir_prep_clean,
                                    str(rank_map.get(job_name, "N/A")), prep_index): job_name
                    for job_name in batch
                }
                for future in as_completed(futures):
                    count += 1
                    res = future.result()
                    symbol = "✓"
                    stat_text = res['status']
                    if stat_text != "Success": symbol = "✗"
                    try:
                        job_id_num = res['job'].split('_')[0]
                    except Exception:
                        job_id_num = "?"
                    print(f"({count}/{total_prep} | ID:{job_id_num}) {symbol} {res['job']} | [{stat_text}]", flush=True)
                    if logger: logger.info(f"{symbol} {res['job']} | [{stat_text}]")
                    # Residue identity guard: run after every successful PrepWizard job
                    if stat_text == "Success":
                        _prepared_pdb = dir_prep_clean / f"{res['job']}_Prepared.pdb"
                        _check_residue_identity_guard(_prepared_pdb, res['job'], CFG)

    print("")
    console_separator()

    n_prep_cached = prep_already_done_count
    n_prep_new    = len(prep_jobs_to_run)
    _print_run_delta_table(n_raw_cached, n_raw_new, n_prep_cached, n_prep_new)


if __name__ == "__main__":
    import time as _time
    _t0 = _time.perf_counter()
    main()
    _utils_mod.print_elapsed(_t0, "05_CIF-PDB_Preparation_FAcDs.py")