#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 05  |  Top-N Selection + PDB Generation & Preparation
===============================================================================
Two phases on the MD-ready cohort only (MD_Selected column from Step 02, §18):
  Phase 1 — CIF → PDB conversion (Gemmi) + Schrödinger PrepWizard preparation.
  Phase 2 — Top-N extraction, Ramachandran validation, molecular handover files,
            and PyMOL/PLIP interaction figures.
Restricting both phases to the ~10 MD-ready complexes keeps this step cheap
instead of converting/preparing the entire predicted library.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 30 June 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 05_TopN_and_PDB_Preparation_FAcDs.py
  Role          : "Builder + Selector" — gates on the MD-ready cohort (MD_Selected),
                  converts CIF outputs to analysis-ready PDB, prepares them with
                  PrepWizard, then extracts and renders that cohort for handover.
  Imports from  : 00_02_Project_Config_FAcDs.py  (CFG — pH values, MD-selection §18)
                  00_03_Project_Utils_FAcDs.py   (ConsoleColours, setup_logging,
                                            console helpers, Ramachandran helpers)
  Reads         : <Run>/2_Best_Complexes_CIFs/*.cif
                  <Run>/1_Boltz2_Production/7_Boltz2_FAcDs_Ranked_*.csv  (MD_Selected)
                  <Run>/1_Boltz2_Production/1_Input_FASTA_and_SMILES/*
  Writes        : <Run>/5_TopN_and_Preparation/  (one consolidated folder)
                    1_Converted_Raw_PDB/   (raw PDBs + Figures/)
                    2_Prepared_PDBs/       (prepared PDBs + Figures/)
                    3_Comparative_Analysis/ (Ramachandran, Controls, handover, combined CSV)
                    00_TopN_and_Preparation_Log.txt  (single log for both phases)
  Upstream      : 02_Production_FAcDs.py  → writes Best_Complexes_CIFs and ranked CSV
  Downstream    : 06_SID_Prime-MMGBSA_FAcDs.py         → reads prepared PDBs / handover
                  07_MD_Thermodynamics_QMMM_Engine_FAcDs.py → reads prepared PDBs for MD/QM-MM
───────────────────────────────────────────────────────────────────────────────

── The Critic's Corner: Known Limitations & Failure Points ──────────────────
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
    python 05_TopN_and_PDB_Preparation_FAcDs.py Boltz-2_Run_20260309T085406Z

── Key features ───────────────────────────────────────────────────────────────
  • Smart Validation: SOURCE_CIF tag in PDB header detects stale files from
    previous runs and triggers automatic regeneration.
  • Stale-file Detection: re-prepares if Raw PDB is newer than Prepared PDB.
  • Gemmi Conversion: robust CIF → PDB; moves non-protein residues to Chain L.
  • Smart Ligand Management: metals (Zn, Mg) and modified residues (MSE) stay
    with the protein chain to preserve topology for Maestro.
  • PrepWizard: fills side chains, PropKa pH 8.0, Epik pH 8.0, 0.3 Å RMSD
    restrained minimisation, disulfide detection.
  • MD-ready gate: both phases process only the MD_Selected cohort (CFG §18) plus
    the control jobs (ID 0000000_*), so the whole predicted library is never
    converted or prepared.
  • Top-N extraction: prepared-complex handover (FASTA/SMILES/SDF), Ramachandran
    validation, and PyMOL/PLIP/matplotlib interaction figures for each complex.
───────────────────────────────────────────────────────────────────────────────

-------------------------------------------------------------------------------
Scientific References:
    1. mmCIF/PDB structure handling (Gemmi):
       - Wojdyr, M. (2022) GEMMI: a library for structural biology. J Open Source Softw 7:4200.
       - DOI: https://doi.org/10.21105/joss.04200
    2. Protein preparation (Schrödinger Protein Preparation Wizard / prepwizard):
       - Sastry, G.M., Adzhigirey, M., Day, T., Annabhimoju, R. & Sherman, W. (2013)
         J Comput-Aided Mol Des 27:221–234. DOI: https://doi.org/10.1007/s10822-013-9644-8
       - Schrödinger Release: Protein Preparation Wizard, Schrödinger, LLC, New York, NY. https://www.schrodinger.com
    3. Protonation-state assignment:
       - PropKa3 (protein pKa): Olsson, M.H.M., Søndergaard, C.R., Rostkowski, M. & Jensen, J.H.
         (2011) J Chem Theory Comput 7:525–537. DOI: https://doi.org/10.1021/ct100578z
       - Epik (ligand protonation/tautomers): Shelley, J.C. et al. (2007)
         J Comput-Aided Mol Des 21:681–691. DOI: https://doi.org/10.1007/s10822-007-9133-z
    4. Data handling:
       - McKinney, W. (2010) Data Structures for Statistical Computing in Python. Proc 9th Python in Science Conf 56–61. DOI: https://doi.org/10.25080/Majora-92bf1922-00a
-------------------------------------------------------------------------------
"""

# ===============================================================================
# SECTION 1: SYSTEM CONFIGURATION & IMPORTS
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 1.1: Standard Library Imports
# -------------------------------------------------------------------------------
import os
import sys

"""
CPU usage cap (total cores − 2; mirrors CFG.PREP_CPU_RESERVE). Reserve 2 cores
for OS/desktop stability by limiting the thread-pool maths libraries (BLAS /
MKL / OpenMP / NumExpr). Must precede numpy/pandas import to take effect;
setdefault() preserves any value exported by the caller or pipeline runner.
"""
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

# --- Additional stack for the Top-N extraction + figure phase ---
import tempfile
import xml.etree.ElementTree as _ET
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend; must be set before pyplot import
import matplotlib.patches as _mpatches
from matplotlib.patches import FancyBboxPatch as _FancyBboxPatch
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.lines import Line2D as _Line2D
from Bio import SeqIO
from Bio.PDB import PDBParser as _PDBParser
from rdkit import Chem
DEFAULT_BASE_PATH = Path.cwd()


# -------------------------------------------------------------------------------
# Step 1.2: Scientific Stack
# -------------------------------------------------------------------------------
import gemmi

# -------------------------------------------------------------------------------
# Step 1.3: Pipeline modules (00_02 (config & utils)) via importlib
# -------------------------------------------------------------------------------
"""
Filenames begin with digits and cannot be imported with standard `import`.
"""
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
# Utility bindings used bare by the Top-N extraction + figure phase.
compute_ramachandran_angles  = _utils_mod.compute_ramachandran_angles
_rama_stats                  = _utils_mod._rama_stats
save_ramachandran_comparison = _utils_mod.save_ramachandran_comparison
safe_name                    = _utils_mod.safe_name
SEPARATOR                    = "-" * 80

# -------------------------------------------------------------------------------
# Step 1.4: Global Constants & Paths
# -------------------------------------------------------------------------------

"""
Ensure SCHRODINGER is set in the process environment so child processes
(PrepWizard subprocess calls) inherit it without requiring a prior `export`.
"""
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

"""
Logging and console functions are provided by 00_03_Project_Utils.
Script-level wrappers capture the module-global `logger` so existing call
sites require no modification.
"""

def setup_logging(prep_base_dir: Path) -> Path:
    """Initialises the preparation log via the shared utility."""
    global logger
    prep_base_dir.mkdir(parents=True, exist_ok=True)
    master_log = prep_base_dir / "00_TopN_and_Preparation_Log.txt"
    logger = _utils_mod.setup_logging(master_log, logger_name="pdb_prep", mode="w", timestamp=True)
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
                """
                If it's a control ID, it is not indexed for renaming
                because multiple controls share this ID.
                """
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
        job_name = re.sub(r"_model_\d+$", "", cif_path.stem)
        mtime = cif_path.stat().st_mtime
        if job_name not in latest_cifs or mtime > latest_cifs[job_name][1]:
            latest_cifs[job_name] = (cif_path, mtime)

    # Return sorted list for deterministic order
    return [(jn, path) for jn, (path, _) in sorted(latest_cifs.items())]

def load_rank_map(prod_dir: Path) -> dict:
    """Loads the master CSV to annotate PDB headers with a score/rank.
    Prefers the FAcDs Ranked CSV (7_Boltz2_FAcDs_Ranked_*), then any ranked CSV,
    then falls back to any master CSV.
    Uses Scientific_Rank if present, otherwise Boltz_Model_Confidence (rounded to 4dp).
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
        if "job_name" not in df.columns:
            return {}
        if "Scientific_Rank" in df.columns:
            return dict(zip(df["job_name"], df["Scientific_Rank"]))
        if "Boltz_Model_Confidence" in df.columns:
            # Use rounded Boltz model confidence as a rank proxy
            return {jn: f"{sc:.4f}" for jn, sc in zip(df["job_name"], df["Boltz_Model_Confidence"])}
    except Exception as e:
        if logger: logger.debug(f"Error loading Rank Map: {e}")
    return {}


def load_md_selected_jobs(prod_dir: Path) -> set:
    """Return the set of job_names flagged MD_Selected in the ranked CSV.

    This is the SECTION 18 gate: only the MD-ready cohort is converted and
    prepared here, so the stage processes ~10 complexes instead of the full
    predicted library. Returns an empty set if the column is absent (older CSV),
    in which case the caller falls back to processing every best-complex CIF.
    """
    md_col = getattr(CFG, "MD_SELECTED_COL", "MD_Selected")
    rank_csvs = (sorted(prod_dir.glob("7_Boltz2_FAcDs_Ranked_*.csv")) or
                 sorted(prod_dir.glob("*_Ranked_*.csv")))
    if not rank_csvs:
        return set()
    try:
        df = pd.read_csv(sorted(rank_csvs)[-1], low_memory=False)
        if md_col not in df.columns or "job_name" not in df.columns:
            return set()
        _sel = df[md_col].astype(str).str.lower().isin(["true", "1", "1.0"])
        return set(df.loc[_sel, "job_name"].astype(str))
    except Exception as e:
        if logger: logger.debug(f"Error loading MD selection: {e}")
        return set()

# -------------------------------------------------------------------------------
# Step 3.2: PDB Manipulation (Header & Source Tracking)
# -------------------------------------------------------------------------------
def get_source_tag(pdb_path: Path) -> str:
    """Reads the REMARK 999 SOURCE_CIF tag from a PDB file."""
    try:
        with open(pdb_path, "r") as f:
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

            # Re-number Chain L residues sequentially so multiple ligand copies
            # (homodimers / multiple identical PFAS) never collide on res.seqid,
            # which would break PyMOL/PLIP selections downstream.
            for _new_seq, _lres in enumerate(lig_chain, start=1):
                _lres.seqid = gemmi.SeqId(_new_seq, " ")

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
    _success = False
    log_file = None
    try:
        input_pdb_local = job_work_dir / f"{job_name}_RAW.pdb"
        output_pdb_name = f"{job_name}_Prepared.pdb"
        shutil.copy2(raw_pdb, input_pdb_local)

        """
        Headless execution environment — strip Python-env overrides that confuse
        Schrödinger's bundled Python interpreter.
        """
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
            _success = True
            return True

        if logger:
            logger.error(f"PrepWizard returned 0 but output missing for {job_name}")
    except Exception as e:
        if logger: logger.error(f"PrepWizard execution failed for {job_name}: {e}")
    finally:
        # Preserve the PrepWizard log on failure (the temp work dir is wiped below).
        if not _success and log_file is not None and log_file.exists():
            try:
                shutil.copy2(log_file, final_dest.parent / f"{job_name}_prepwizard_FAILED.log")
            except Exception:
                pass
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

def load_catalytic_anchor_map(prod_dir: Path) -> dict:
    """Per-job catalytic residue numbers taken from 02's dynamic global alignment
    (Mapped_Nucleophile / Mapped_Acid / Mapped_Base columns of the ranked CSV).
    Lets the identity guard anchor on each homolog's actual aligned positions —
    robust to insertions/deletions — instead of a static ±window around the
    canonical reference numbers. Returns {job_name: {"Nuc"/"Acid"/"Base": int}}."""
    import re as _re
    rank_csvs = sorted(prod_dir.glob("7_Boltz2_FAcDs_Ranked_*.csv")) or sorted(prod_dir.glob("*_Ranked_*.csv"))
    if not rank_csvs:
        return {}
    try:
        df = pd.read_csv(sorted(rank_csvs)[-1], low_memory=False)
    except Exception:
        return {}
    _cols = {"Nuc": "Mapped_Nucleophile", "Acid": "Mapped_Acid", "Base": "Mapped_Base"}
    if "job_name" not in df.columns or not all(c in df.columns for c in _cols.values()):
        return {}
    def _num(v):
        m = _re.search(r"(\d+)\s*$", str(v))
        return int(m.group(1)) if m else None
    return {row["job_name"]: {k: _num(row[c]) for k, c in _cols.items()}
            for _, row in df.iterrows()}


def _check_residue_identity_guard(prepared_pdb_path: Path, job_name: str, cfg, anchors: dict = None) -> dict:
    """Parse prepared PDB; verify catalytic residues at expected positions, writing
    a {job_name}_index_offset.json QC record. NOTE: this is a Step-05 QC diagnostic
    only — Step 07 maps catalytic residues from its OWN dynamic global alignment
    (Full_Sequence_Alignment_Map), so it does not consume this file.

    When `anchors` (02's dynamic Mapped_* positions for this job) is supplied the
    search is anchored there — robust to homolog insertions/deletions — rather than
    on a static window around the canonical reference numbers."""
    import json

    anchors = anchors or {}
    nuc_ref  = anchors.get("Nuc")  or cfg.DREAM_TEAM_REFS["Nuc"]   # 02 dynamic align, else canonical 110
    acid_ref = anchors.get("Acid") or cfg.DREAM_TEAM_REFS["Acid"]  # else 134
    base_ref = anchors.get("Base") or cfg.DREAM_TEAM_REFS["Base"]  # else 277

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

    def find_type_near(ref, expected_types, window=cfg.SMART_LOCK_RESNUM_WINDOW):
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

    """
    Use the most common non-zero offset (consensus across Nuc/Acid/Base);
    fall back to 0 when all three agree on the reference numbering.
    """
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


def load_reference_data(input_data_dir: Path):
    """
    Scans the 1_Input_Data folder to map names to Sequences and SMILES.
    Returns: (seq_map, smi_map, fasta_count, smi_count)
    """
    seq_map = {}
    smi_map = {}
    fasta_count = 0
    smi_count = 0

    if not input_data_dir.exists():
        return seq_map, smi_map, 0, 0

    # 1. Load FASTA
    if SeqIO:
        for f in list(input_data_dir.glob("*.fasta")) + list(input_data_dir.glob("*.fa")):
            try:
                for r in SeqIO.parse(str(f), "fasta"):
                    fasta_count += 1
                    clean_id = safe_name(re.split(r"[ \t|()-]", r.id)[0])
                    _seq = str(r.seq)
                    seq_map[clean_id] = _seq
                    '''
                    Also key by the merge-index-prefix-stripped name, which is the form
                    carried in the CSV Protein_Name column (e.g. "1_DeHa4" → "DeHa4");
                    without this the per-candidate sequence lookup fails and the handover
                    FASTA is written as "Sequence_Not_Found".
                    '''
                    _stripped = re.sub(r"^\d+_", "", clean_id)
                    if _stripped in seq_map and seq_map[_stripped] != _seq:
                        # Two prefixed entries collapse to the same bare name but carry
                        # different sequences; overwriting would silently hand the wrong
                        # sequence to the earlier candidate. Keep the first, warn loudly.
                        console_info(
                            f"Warning: stripped sequence key '{_stripped}' collides with a "
                            f"differing sequence (from '{clean_id}'); keeping the first occurrence.")
                    else:
                        seq_map[_stripped] = _seq
                    # Store original ID as fallback
                    seq_map[r.id] = _seq
            except Exception as e:
                console_info(f"Warning: Could not parse FASTA {f.name}: {e}")

    # 2. Load SMILES
    for s in input_data_dir.glob("*.smi"):
        try:
            with open(s) as f:
                for line in f:
                    if line.strip(): # Check empty lines
                        smi_count += 1
                        parts = line.strip().split()
                        if len(parts) >= 1:
                            smi = parts[0]
                            # If name exists, use it, else generic
                            name = parts[1] if len(parts) > 1 else "LIG"
                            clean_name = safe_name(name)
                            smi_map[clean_name] = smi
        except Exception as e:
            console_info(f"Warning: Could not parse SMILES {s.name}: {e}")

    return seq_map, smi_map, fasta_count, smi_count

# -------------------------------------------------------------------------------
# Step 4.2: Molecule Handling
# -------------------------------------------------------------------------------
def extract_chain_l_mol(pdb_path: Path):
    """Extracts Chain L (Ligand) lines from PDB and returns RDKit Mol."""
    ligand_lines = []
    try:
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    if len(line) > 21 and line[21] == "L":
                        ligand_lines.append(line)
                elif line.startswith("END"):
                    break
        block = "".join(ligand_lines)
        if not block: return None
        return Chem.MolFromPDBBlock(block, removeHs=False, sanitize=False)
    except Exception as e:
        if logger: logger.debug(f"Failed to extract Chain L mol from {pdb_path.name}: {e}")
        return None


# ===============================================================================
# SECTION 5: MAIN EXECUTION LOGIC
# ===============================================================================


def _fig_constants():
    """Return a dict of figure-generation constants from CFG."""
    return {
        "IMG_WIDTH":      CFG.VIS_IMG_WIDTH,
        "IMG_HEIGHT":     CFG.VIS_IMG_HEIGHT,
        "RAY_TRACE":      CFG.VIS_RAY_TRACE,
        "DIST_HBOND":     CFG.THRESHOLD_HB_DIST_MAX,   # heavy-atom D···A proxy (H stripped from PDB)
        "DIST_SALT":      CFG.THRESHOLD_SALT_BRIDGE,   # Maestro salt-bridge cutoff (single source)
        "DIST_F_CONTACT": CFG.VIS_F_CONTACT_RADIUS,
        "DIST_POCKET":    CFG.VIS_POCKET_RADIUS,
        "TIMEOUT_PYMOL":  CFG.VIS_TIMEOUT_PYMOL,
        "TIMEOUT_PLIP":   CFG.VIS_TIMEOUT_PLIP,
    }

"""
Per-interaction-type distance ceilings for filtering PLIP output down to the
project's Schrödinger-Maestro criteria (CFG §3 — single source of truth).
PLIP's own internal cutoffs are looser (e.g. H-bond 4.1 Å, salt 5.5 Å); any
PLIP contact whose reported distance exceeds the matching CFG ceiling is
dropped so all three figure engines agree on what counts as a bond.
"""
def _plip_cfg_cutoffs():
    return {
        "hbond":       CFG.THRESHOLD_HB_DIST_MAX,    # D···A proxy (PLIP reports D-A distance)
        "halogen":     CFG.HALOGEN_BOND_DIST_MAX,
        "salt":        CFG.THRESHOLD_SALT_BRIDGE,
        "pistack":     CFG.PI_STACK_EDGE_DIST_MAX,   # widest π-π ceiling (edge–face)
        "pication":    CFG.THRESHOLD_PI_CATION_MAX,
        "hydrophobic": CFG.THRESHOLD_HYDROPHOBIC_MAX,
        "water":       CFG.WATER_BRIDGE_DIST_MAX,   # Maestro water-bridge ceiling
    }

# Metadata cache: job_name → {p: protein_name, l: ligand_name}
METADATA_CACHE: dict = {}

def load_metadata(ext_dir: Path):
    """Populates METADATA_CACHE from any *_Scientific_Data.csv found under ext_dir."""
    global METADATA_CACHE
    try:
        import pandas as pd
        for csv in ext_dir.rglob("*_Scientific_Data.csv"):
            try:
                df = pd.read_csv(csv)
                for _, row in df.iterrows():
                    jn = str(row.get("job_name", ""))
                    if jn:
                        METADATA_CACHE[jn] = {
                            "p": str(row.get("Protein_Name", "Protein")),
                            "l": str(row.get("Ligand_Name", "PFAS")),
                        }
            except Exception:
                pass
    except ImportError:
        pass  # pandas optional; filenames fall back to base_name

# -------------------------------------------------------------------------------
# --- Figure engine: logger (separate from the extraction logger) ---
# -------------------------------------------------------------------------------
_fig_logger = None

def _fig_log(text=""):
    """Unified logging for figure generation — routes through console_info."""
    _utils_mod.console_info(text, _fig_logger)

def _append_auxiliary_log(tool_name, complex_name, log_path):
    """Reads a specific tool's output log and appends it to the main report."""
    if not _fig_logger or not log_path.exists():
        return
    try:
        content = log_path.read_text()
        if content.strip():
            border = "=" * 60
            _fig_logger.info(f"\n{border}\nFIGURE LOG: {tool_name} | {complex_name}\n{border}\n")
            _fig_logger.info(content + "\n")
    except Exception as e:
        _fig_logger.debug(f"Failed to append figure log: {e}")

# -------------------------------------------------------------------------------
# --- Figure engine: subprocess runner ---
# -------------------------------------------------------------------------------
def _run_cmd(cmd_list, cwd=None, env=None, timeout=60, log_file=None):
    """Run external command; explicitly kill the subprocess on timeout so no orphans remain."""
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    run_env["QT_QPA_PLATFORM"]      = "offscreen"
    run_env["LIBGL_ALWAYS_SOFTWARE"] = "1"
    if "DISPLAY" not in run_env:
        run_env["DISPLAY"] = ":99"
    """
    Cap internal threading in numpy/OpenBabel/MKL so each subprocess uses
    exactly 1 thread; concurrency is controlled at the pool level instead.
    """
    run_env.setdefault("OMP_NUM_THREADS",      "1")
    run_env.setdefault("MKL_NUM_THREADS",      "1")
    run_env.setdefault("OPENBLAS_NUM_THREADS", "1")
    run_env.setdefault("NUMEXPR_NUM_THREADS",  "1")

    f_handle = None
    proc     = None
    try:
        if log_file:
            f_handle = open(log_file, "w")
        out_dest = f_handle              if f_handle else subprocess.DEVNULL
        err_dest = subprocess.STDOUT    if f_handle else subprocess.DEVNULL
        proc = subprocess.Popen(cmd_list, cwd=cwd, env=run_env,
                                stdout=out_dest, stderr=err_dest)
        proc.communicate(timeout=timeout)
        rc = proc.returncode
        if f_handle:
            f_handle.close()
        return rc == 0
    except subprocess.TimeoutExpired:
        if proc:
            proc.kill()
            proc.communicate()
        if f_handle:
            f_handle.write(f"\n[TIMEOUT after {timeout}s]\n")
            f_handle.close()
        return False
    except Exception as e:
        if proc:
            try:
                proc.kill()
                proc.wait()          # reap the killed child so it does not linger as a zombie
            except Exception: pass
        if f_handle:
            f_handle.write(f"\n[ERROR]: {e}\n")
            f_handle.close()
        return False
    except BaseException:
        # KeyboardInterrupt / SIGTERM: reap the PyMOL/PLIP child before the
        # exception unwinds, so it is not orphaned, then re-raise to abort.
        if proc and proc.poll() is None:
            try:
                proc.kill(); proc.wait()
            except Exception: pass
        if f_handle:
            try: f_handle.close()
            except Exception: pass
        raise

# -------------------------------------------------------------------------------
# --- Figure engine: visualisation engines (PyMOL / PLIP) ---
# -------------------------------------------------------------------------------
def _run_pymol(pdb_path, fig_root, log_dir, base_name, lig_name, lig_num, has_f, sw, C):
    if not sw["PyMOL"]: return "Missing"
    out_dir = fig_root / "PyMOL"
    out_dir.mkdir(exist_ok=True, parents=True)

    out_opaque      = out_dir / f"{base_name}_Interaction_Opaque.png"
    out_transparent = out_dir / f"{base_name}_Interaction_Transparent.png"
    out_pse         = out_dir / f"{base_name}_Session.pse"
    pml             = out_dir / f"{base_name}_render.pml"
    log_path        = log_dir / f"{base_name}_PyMOL.log"

    """
    Precise ligand selection: chain L residue number is authoritative.
    Using AND (not OR) prevents UNK/non-standard protein residues also placed
    in chain L by Schrödinger prep from being selected as ligand atoms.
    """
    lig_sel = f"(chain L and resi {lig_num})"

    """
    Cap PyMOL ray-tracing threads: 2 PyMOL workers run concurrently, so each
    should use at most (cpu-2)//2 ray threads to stay within the cpu-2 budget.
    """
    _max_rt = max(1, CFG.GLOBAL_MAX_WORKERS // 2)

    script = [
        "reinitialize",
        f'load "{pdb_path.resolve()}", complex',
        "remove solvent",
        # Keep polar (N–H / O–H) hydrogens so PyMOL's H-bond detection uses real
        # donor–H···acceptor geometry; drop non-polar C–H to reduce visual clutter.
        "remove (hydro and not (neighbor (elem N+O)))",
        # Lighting — two_sided illuminates inside faces of the binding pocket,
        # preventing the cavity from rendering as a solid black void.
        f"bg_color {CFG.VIS_PYMOL_BG_COLOR}",
        "set ray_trace_mode, 1", "set ray_trace_gain, 0.12",
        "set ray_shadows, 1",
        "set ambient, 0.45", "set direct, 0.55", "set specular, 0.35",
        "set two_sided_lighting, 1",   # illuminate interior surface faces
        "set antialias, 2", "set orthoscopic, on",
        "set ambient_occlusion_mode, 0",  # AO darkens cavities — disable to keep pocket visible
        "set surface_quality, 1",
        f"set max_threads, {_max_rt}",
        # Selections
        f"select ligand, {lig_sel}",
        f"select interacting, byres polymer.protein within {C['DIST_F_CONTACT']} of ligand",
        "hide all",
        # Full protein surface — grey, 50% transparent so ligand sticks show through
        "show surface, polymer.protein",
        "color gray65, polymer.protein",
        f"set transparency, {CFG.VIS_PYMOL_SURFACE_TRANSPARENCY}",
        # Protein cartoon underneath surface
        "show cartoon, polymer.protein",
        "color gray70, polymer.protein",
        f"set cartoon_transparency, {CFG.VIS_PYMOL_CARTOON_TRANSPARENCY}",
        "set cartoon_fancy_helices, 1", "set cartoon_fancy_sheets, 1",
        # Interacting residues — thin grey sticks, element-coloured heteroatoms
        "show sticks, interacting",
        "color gray85, interacting and elem C",
        "color red,    interacting and elem O",
        "color blue,   interacting and elem N",
        "color yellow, interacting and elem S",
        "set stick_radius, 0.15, interacting",
        # Ligand — thick bright orange sticks; highly visible against grey surface
        # orange chosen for maximum contrast on grey without clashing with interaction colours
        "show sticks, ligand",
        "color tv_orange, ligand",
        "color tv_green, ligand and elem F",
        "color red,      ligand and elem O",
        "color blue,     ligand and elem N",
        "set stick_radius, 0.32, ligand",
        # H-bonds — yellow dashes
        f"dist hbonds, ligand, interacting, mode=2, cutoff={C['DIST_HBOND']}",
        "color yellow, hbonds",
        "set dash_gap, 0.20", "set dash_width, 3.5", "set dash_color, yellow",
        # Salt bridges — cyan dashes (cutoff from CFG.THRESHOLD_SALT_BRIDGE — Maestro)
        f"dist salt_bridge, (ligand and (name O*,N*)), (interacting and (name N*,O*)), mode=2, cutoff={C['DIST_SALT']}",
        "color cyan, salt_bridge",
        # Residue labels — white, bold, readable on grey background
        "set label_size, 9", "set label_font_id, 7",
        f"set label_color, {CFG.VIS_PYMOL_LABEL_COLOR}",
        "label (interacting and n. CA), '%s%s' % (resn.capitalize(), resi)",
        # Orient to face the binding pocket, clip front protein to expose pocket interior
        "orient (interacting or ligand)",
        "zoom ligand, 13.0",
        "clip near, 4, ligand",  # clip 4 Å in front of ligand centre, exposing the pocket
        # Opaque render
        "set ray_opaque_background, on",
        f"png {out_opaque.resolve()}, {C['IMG_WIDTH']}, {C['IMG_HEIGHT']}, ray=1",
        # Transparent render
        "set ray_opaque_background, off",
        f"png {out_transparent.resolve()}, {C['IMG_WIDTH']}, {C['IMG_HEIGHT']}, ray=1",
        f"save {out_pse.resolve()}",
        "quit",
    ]

    pml.write_text("\n".join(script))
    timeout = CFG.VIS_TIMEOUT_PYMOL_HEAVY if has_f else CFG.VIS_TIMEOUT_PYMOL
    try:
        success = _run_cmd([sw["PyMOL"], "-c", "-q", str(pml)], timeout=timeout, log_file=log_path)
    finally:
        if pml.exists():
            pml.unlink(missing_ok=True)
    if not success:
        _append_auxiliary_log("PyMOL", base_name, log_path)
        if log_path.exists():
            tail = log_path.read_text(errors="replace")[-600:]
            _fig_log(f"  [!] PyMOL log tail for {base_name}:\n{tail}")
    return "Success" if success else "Failed"

def _plip_coo(el):
    """Parse a PLIP <ligcoo> or <protcoo> element → numpy array [x, y, z]."""
    return np.array([float(el.find("x").text),
                     float(el.find("y").text),
                     float(el.find("z").text)])


def _parse_plip_xml(xml_path, lig, pro):
    """
    Parse a PLIP XML report into a contacts list compatible with _im_project /
    _im_render_diagram.

    Each contact dict carries:
      key, resname, resnum, chain, dist, itype, is_hbond, is_salt,
      lig_atom (nearest ligand atom), prot_atom (''), center (3-D protcoo)
    """
    _ITYPE_ORDER = ["hbond", "halogen", "salt", "water", "pistack", "pication",
                    "hydrophobic", "contact"]

    tree = _ET.parse(str(xml_path))
    root = tree.getroot()

    lig_pos = np.array([a["pos"] for a in lig])

    contacts_by_key = {}   # key → best contact dict (deduplicate by residue)

    def _closest_lig_atom(ligcoo):
        dists = np.linalg.norm(lig_pos - ligcoo, axis=1)
        idx   = int(dists.argmin())
        return lig[idx], float(dists[idx])

    _cfg_cut = _plip_cfg_cutoffs()   # Maestro/CFG distance ceilings per interaction type

    def _register(resnr, restype, reschain, dist, itype, protcoo_3d, ligcoo_3d):
        key = (reschain, int(resnr), restype)
        """
        Clamp PLIP's looser internal cutoffs to the project's CFG (Maestro) criteria
        so PLIP, PyMOL, and InteractionMap all agree on what counts as a bond.
        """
        ceil = _cfg_cut.get(itype)
        if ceil is not None and dist > ceil:
            return
        la, _  = _closest_lig_atom(ligcoo_3d)
        # C–F bonds are leaving groups in FAcDs SN2; exclude F from H-bond / halogen contacts
        if la["elem"] == "F" and itype in ("hbond", "halogen"):
            return
        is_hbond  = itype == "hbond"
        is_salt   = itype == "salt"
        entry = {"key": key, "resname": restype, "resnum": int(resnr),
                 "chain": reschain, "dist": dist, "itype": itype,
                 "is_hbond": is_hbond, "is_salt": is_salt,
                 "lig_atom": la, "prot_atom": "", "center": protcoo_3d}
        # Keep highest-priority interaction type per residue
        if key not in contacts_by_key:
            contacts_by_key[key] = entry
        else:
            existing = contacts_by_key[key]
            if _ITYPE_ORDER.index(itype) < _ITYPE_ORDER.index(existing["itype"]):
                contacts_by_key[key] = entry

    for bs in root.findall("bindingsite"):
        iact = bs.find("interactions")
        if iact is None:
            continue

        for hb in iact.findall("./hydrogen_bonds/hydrogen_bond"):
            protcoo = _plip_coo(hb.find("protcoo"))
            ligcoo  = _plip_coo(hb.find("ligcoo"))
            dist    = float(hb.find("dist_d-a").text)
            _register(hb.find("resnr").text, hb.find("restype").text,
                      hb.find("reschain").text, dist, "hbond", protcoo, ligcoo)

        for hx in iact.findall("./halogen_bonds/halogen_bond"):
            protcoo = _plip_coo(hx.find("protcoo"))
            ligcoo  = _plip_coo(hx.find("ligcoo"))
            dist    = float(hx.find("dist").text)
            _register(hx.find("resnr").text, hx.find("restype").text,
                      hx.find("reschain").text, dist, "halogen", protcoo, ligcoo)

        for sb in iact.findall("./salt_bridges/salt_bridge"):
            protcoo = _plip_coo(sb.find("protcoo"))
            ligcoo  = _plip_coo(sb.find("ligcoo"))
            dist    = float(sb.find("dist").text)
            _register(sb.find("resnr").text, sb.find("restype").text,
                      sb.find("reschain").text, dist, "salt", protcoo, ligcoo)

        for wb in iact.findall("./water_bridges/water_bridge"):
            protcoo = _plip_coo(wb.find("protcoo"))
            ligcoo  = _plip_coo(wb.find("ligcoo"))
            dist_d  = float(wb.find("dist_d-w").text)
            dist_a  = float(wb.find("dist_a-w").text)
            _register(wb.find("resnr").text, wb.find("restype").text,
                      wb.find("reschain").text, max(dist_d, dist_a), "water",
                      protcoo, ligcoo)

        for ps in iact.findall("./pi_stacks/pi_stack"):
            protcoo = _plip_coo(ps.find("protcoo"))
            ligcoo  = _plip_coo(ps.find("ligcoo"))
            dist    = float(ps.find("dist").text)
            _register(ps.find("resnr").text, ps.find("restype").text,
                      ps.find("reschain").text, dist, "pistack", protcoo, ligcoo)

        for pc in iact.findall("./pi_cation_interactions/pi_cation_interaction"):
            protcoo = _plip_coo(pc.find("protcoo"))
            ligcoo  = _plip_coo(pc.find("ligcoo"))
            dist    = float(pc.find("dist").text)
            _register(pc.find("resnr").text, pc.find("restype").text,
                      pc.find("reschain").text, dist, "pication", protcoo, ligcoo)

        for hp in iact.findall("./hydrophobic_interactions/hydrophobic_interaction"):
            protcoo = _plip_coo(hp.find("protcoo"))
            ligcoo  = _plip_coo(hp.find("ligcoo"))
            dist    = float(hp.find("dist").text)
            _register(hp.find("resnr").text, hp.find("restype").text,
                      hp.find("reschain").text, dist, "hydrophobic", protcoo, ligcoo)

        # bs_residues with contact="True" and no typed interaction → "contact"
        for bsr in bs.findall("./bs_residues/bs_residue"):
            if bsr.get("contact", "False") != "True":
                continue
            text   = bsr.text.strip()           # e.g. "150A"
            chain  = text[-1]
            resnr  = text[:-1]
            restype = bsr.get("aa", "UNK")
            key    = (chain, int(resnr), restype)
            if key not in contacts_by_key:
                # No typed interaction — approximate position from PDB
                pkey = next((k for k in pro if k[0] == chain and k[1] == int(resnr)), None)
                if pkey:
                    protcoo = pro[pkey]["center"]
                    la, _ = _closest_lig_atom(protcoo)
                    contacts_by_key[key] = {
                        "key": key, "resname": restype, "resnum": int(resnr),
                        "chain": chain, "dist": float(bsr.get("min_dist", CFG.PLIP_CONTACT_FALLBACK_DIST)),
                        "itype": "contact", "is_hbond": False, "is_salt": False,
                        "lig_atom": la, "prot_atom": "", "center": protcoo,
                    }

    return sorted(contacts_by_key.values(), key=lambda x: x["dist"])


def _im_render_diagram(lig_2d, lig, res_2d, contacts, out_png, mode="distance"):
    """
    Shared matplotlib renderer for both InteractionMap (distance-based) and
    PLIP (XML-based) diagrams.  mode='distance' or 'plip'.
    """
    # Interaction type → (linewidth, linestyle, colour, show_dist_label)
    _ITYPE_STYLE = {
        "hbond":       (2.2, (0, (6, 3)),  "#E67E22", True),
        "arom_hbond":  (1.9, (0, (5, 2, 1, 2)), "#16A085", True),
        "halogen":     (2.0, (0, (4, 2)),  "#1D8348", True),
        "salt":        (2.0, (0, (3, 2)),  "#C0392B", True),
        "water":       (1.5, (0, (2, 2)),  "#5DADE2", True),
        "pistack":     (1.8, (0, (5, 2)),  "#2471A3", False),
        "pication":    (1.8, (0, (4, 2)),  "#7D3C98", False),
        "hydrophobic": (0.9, (0, (2, 4)),  "#BDC3C7", False),
        "contact":     (0.9, (0, (2, 4)),  "#BDC3C7", False),
    }
    _DIST_COL = {"hbond": ("#884400","#FDEBD0","#E67E22"),
                 "arom_hbond": ("#0B5345","#D1F2EB","#16A085"),
                 "halogen": ("#0A3D0A","#D5F5E3","#1D8348"),
                 "salt": ("#7B241C","#FADBD8","#C0392B"),
                 "water": ("#1A5276","#D6EAF8","#5DADE2")}

    # Dynamic axis bounds — zoom in when few residues to eliminate blank space
    _all_x = [p[0] for p in lig_2d] + [v[0] for v in res_2d.values()]
    _all_y = [p[1] for p in lig_2d] + [v[1] for v in res_2d.values()]
    if _all_x and _all_y:
        _pad = max(1.5, (max(_all_x) - min(_all_x)) * 0.20, (max(_all_y) - min(_all_y)) * 0.20)
        _xlo = min(_all_x) - _pad - 1.0   # 1.0 for residue box half-width
        _xhi = max(_all_x) + _pad + 1.0
        _ylo = min(_all_y) - _pad - 0.6
        _yhi = max(_all_y) + _pad + 0.6
        # Square view required for aspect='equal'
        _span = max(_xhi - _xlo, _yhi - _ylo)
        _cx = (_xlo + _xhi) / 2;  _cy = (_ylo + _yhi) / 2
        _xlo, _xhi = _cx - _span / 2, _cx + _span / 2
        _ylo, _yhi = _cy - _span / 2, _cy + _span / 2
    else:
        _xlo, _xhi, _ylo, _yhi = -7.5, 7.5, -7.5, 7.5

    fig = Figure(figsize=(12, 12), facecolor="white")
    canvas = FigureCanvasAgg(fig)
    ax  = fig.add_subplot(111, aspect="equal")
    ax.axis("off")
    ax.set_xlim(_xlo, _xhi)
    ax.set_ylim(_ylo, _yhi)

    # Binding pocket background
    ax.add_patch(_mpatches.Circle((0, 0), 4.2, color="#EAF2FF", zorder=0, alpha=0.6))
    ax.add_patch(_mpatches.Circle((0, 0), 4.2, color="#AED6F1", fill=False,
                            linewidth=1.2, linestyle="--", zorder=0, alpha=0.4))
    ax.text(0, -3.7, "Binding Pocket", ha="center", fontsize=7.5,
            color="#85929E", style="italic", zorder=1)

    # Interaction lines
    name2idx = {a["name"]: i for i, a in enumerate(lig)}
    for c in contacts:
        rpos = res_2d[c["key"]]
        li   = name2idx.get(c["lig_atom"]["name"], 0) if c.get("lig_atom") else 0
        lap  = lig_2d[li]
        itype = c.get("itype", "hbond" if c.get("is_hbond") else
                      "salt" if c.get("is_salt") else "contact")
        lw, ls, col, show_dist = _ITYPE_STYLE.get(itype, _ITYPE_STYLE["contact"])
        ax.plot([lap[0], rpos[0]], [lap[1], rpos[1]],
                lw=lw, ls=ls, color=col, alpha=0.9, zorder=2,
                solid_capstyle="round")
        if show_dist:
            mx, my = (lap[0]+rpos[0])/2, (lap[1]+rpos[1])/2
            tc, fc, ec = _DIST_COL.get(itype, ("#555","#EEE","#999"))
            ax.text(mx, my, f"{c['dist']:.1f} Å",
                    fontsize=8.5, ha="center", va="center", fontweight="bold",
                    color=tc,
                    bbox=dict(fc=fc, ec=ec, alpha=0.88,
                              boxstyle="round,pad=0.22", linewidth=0.8),
                    zorder=8)

    # Ligand bonds
    for i, a in enumerate(lig):
        for j, b in enumerate(lig):
            if j <= i: continue
            if np.linalg.norm(a["pos"] - b["pos"]) < CFG.LIG_COVALENT_BOND_DIST:
                p1, p2 = lig_2d[i], lig_2d[j]
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                        color="#2C3E50", lw=2.8, solid_capstyle="round",
                        zorder=4, alpha=0.85)

    # Ligand atoms
    for i, a in enumerate(lig):
        xy  = lig_2d[i]
        col = _IM_ELEM_COLORS.get(a["elem"], _IM_ELEM_COLORS["other"])
        r   = {"C":0.17,"N":0.19,"O":0.19,"F":0.17,"S":0.22}.get(a["elem"], 0.15)
        ax.add_patch(_mpatches.Circle(xy, r, color=col, zorder=5, ec="white", lw=1.4))
        if a["elem"] not in ("C",):
            ax.text(xy[0], xy[1], a["elem"],
                    ha="center", va="center", fontsize=6.5,
                    color="white", fontweight="bold", zorder=6)

    # Residue boxes
    bw, bh = 1.5, 0.68
    for c in contacts:
        rpos  = res_2d[c["key"]]
        col   = _IM_RES_COLORS.get(c["resname"], "#717D7E")
        label = f"{c['resname']} {c['resnum']}"
        itype = c.get("itype", "contact")
        _BADGE = {
            "hbond": "H-bond", "arom_hbond": "arom. H-bond",
            "halogen": "halogen bond", "salt": "salt bridge",
            "water": "water bridge", "pistack": "π-stack", "pication": "π-cation",
            "hydrophobic": "hydrophobic", "contact": "contact",
        }
        badge = _BADGE.get(itype, "contact")
        # Annotate steric-contact quality (Bondi vdW ratio: good/bad/ugly)
        if itype == "contact" and c.get("quality") in ("bad", "ugly"):
            badge = f"contact ({c['quality']})"
        ax.add_patch(_FancyBboxPatch(
            (rpos[0]-bw/2+0.04, rpos[1]-bh/2-0.04), bw, bh,
            boxstyle="round,pad=0.1", facecolor="#C0C0C0",
            alpha=0.22, zorder=5, linewidth=0))
        ax.add_patch(_FancyBboxPatch(
            (rpos[0]-bw/2, rpos[1]-bh/2), bw, bh,
            boxstyle="round,pad=0.1", facecolor=col,
            edgecolor="white", linewidth=1.6, alpha=0.95, zorder=6))
        ax.text(rpos[0], rpos[1]+0.10, label,
                ha="center", va="center", fontsize=9.5,
                fontweight="bold", color="white", zorder=7)
        ax.text(rpos[0], rpos[1]-0.18, badge,
                ha="center", va="center", fontsize=7,
                color="white", alpha=0.9, zorder=7)

    # Unified legend (bottom, all entries in one block)
    _leg = [
        _mpatches.Patch(color="#E67E22", label="H-bond"),
        _mpatches.Patch(color="#16A085", label="Aromatic H-bond"),
        _mpatches.Patch(color="#1D8348", label="Halogen bond"),
        _mpatches.Patch(color="#C0392B", label="Salt bridge"),
        _mpatches.Patch(color="#5DADE2", label="Water bridge"),
        _mpatches.Patch(color="#2471A3", label="π-stack"),
        _mpatches.Patch(color="#7D3C98", label="π-cation"),
        _mpatches.Patch(color="#BDC3C7", label="Contact"),
        _Line2D([],[],color="none", label=""),
        _mpatches.Patch(color="#C0392B", label="ASP/GLU"),
        _mpatches.Patch(color="#2471A3", label="ARG/LYS"),
        _mpatches.Patch(color="#1E8449", label="HIS"),
        _mpatches.Patch(color="#7D3C98", label="TRP/PHE"),
        _mpatches.Patch(color="#BA4A00", label="TYR"),
        _mpatches.Patch(color="#117A65", label="SER/THR/ASN/GLN"),
        _mpatches.Patch(color="#626567", label="Hydrophobic"),
        _Line2D([],[],color="none", label=""),
    ]
    # Append atom entries inline so everything sits in one box
    for elem, ec in [("C","#2C3E50"),("N","#1A5276"),("O","#A93226"),("F","#1D8348")]:
        _leg.append(_mpatches.Patch(color=ec, label=f"Lig {elem}"))

    """
    Dynamic legend: anchor just below the lowest content point (circle bottom
    or lowest residue box, whichever is further down), with a small gap.
    Uses data-space coordinates so placement adapts to content density.
    """
    _y_res_bottom = (min(pos[1] - bh / 2 for pos in res_2d.values())
                     if res_2d else -4.2)
    _y_legend_top = min(_y_res_bottom, -4.2) - 0.35
    _cx_data      = (_xlo + _xhi) / 2
    ax.legend(handles=_leg, loc="upper center", ncol=4,
              bbox_to_anchor=(_cx_data, _y_legend_top),
              bbox_transform=ax.transData,
              fontsize=8.5, frameon=True, framealpha=0.95,
              edgecolor="#CCCCCC",
              title=f'{"PLIP interactions" if mode=="plip" else "Interactions"}  |  Residue type  |  Ligand atoms',
              title_fontsize=8.5, columnspacing=0.8, handlelength=1.1,
              borderpad=0.6)

    try:
        canvas.print_figure(str(out_png), dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight",
                            pad_inches=0.04, facecolor="white", edgecolor="none")
    except Exception as _render_err:
        raise _render_err
    finally:
        del fig, canvas


def _run_plip(pdb_path, fig_root, log_dir, base_name, sw, C):
    """Run PLIP (XML only), parse typed interactions, render matplotlib diagram."""
    if not sw["PLIP"]: return "Missing"
    out_dir = fig_root / "PLIP"
    out_dir.mkdir(exist_ok=True, parents=True)
    log_path = log_dir / f"{base_name}_PLIP.log"

    tmp_dir = Path(tempfile.mkdtemp(prefix="plip_"))
    try:
        # Handle both binary and module-based PLIP invocation
        if sw["PLIP"].endswith(":plip_module"):
            cmd_list = [sys.executable, "-m", "plip", "-f", str(pdb_path), "-o", str(tmp_dir),
                        "--name", base_name, "-x"]
        else:
            cmd_list = [sw["PLIP"], "-f", str(pdb_path), "-o", str(tmp_dir),
                        "--name", base_name, "-x"]

        if not _run_cmd(cmd_list, timeout=C["TIMEOUT_PLIP"], log_file=log_path):
            _append_auxiliary_log("PLIP", base_name, log_path)
            return "Failed"

        xml_files = list(tmp_dir.glob("*.xml"))
        if not xml_files:
            return "Failed"

        lig, pro, _lig_hb = _im_parse_pdb(pdb_path)
        if not lig:
            return "Failed"

        rendered = False
        for xml_path in xml_files:
            contacts = _parse_plip_xml(xml_path, lig, pro)
            if not contacts:
                continue
            lig_2d, res_2d = _im_project(lig, contacts)
            lig_2d = _im_separate_atoms(lig_2d)
            out_png = out_dir / f"{base_name}_PLIP.png"
            _im_render_diagram(lig_2d, lig, res_2d, contacts, out_png, mode="plip")
            rendered = True

        return "Success" if rendered else "Failed"
    except Exception as _e:
        _fig_log(f"  [!] PLIP render failed for {base_name}: {type(_e).__name__}: {_e}")
        return "Failed"
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# -------------------------------------------------------------------------------
# --- Figure engine: matplotlib interaction diagram (pure Python, no PyMOL) ---
# -------------------------------------------------------------------------------

_IM_CONTACT_DIST  = CFG.CATALYTIC_DIST_CUTOFF
_IM_HBOND_DIST    = CFG.THRESHOLD_HB_DIST_MAX
_IM_SALT_DIST     = CFG.THRESHOLD_SALT_BRIDGE   # Maestro salt-bridge cutoff (single source)
_IM_HBOND_DONORS  = {"N","NH1","NH2","NE","NE2","ND1","ND2","NZ","OG","OG1","OH","NE1"}
_IM_HBOND_ACC     = {"O","OD1","OD2","OE1","OE2","OG","OG1","OH","NE2","ND1",
                     "O1","O2","O3","O4","O5","O6"}
# Maestro H-bond / aromatic-H-bond / contact criteria — all sourced from CFG (single source).
_IM_HB_HA_DIST    = CFG.THRESHOLD_HB_DIST_HA      # Å  H···A maximum (explicit-H structures)
_IM_HB_DON_ANGLE  = CFG.THRESHOLD_HB_ANGLE_MIN    # °  donor minimum angle (D–H···A)
_IM_DH_BOND_MAX   = CFG.HB_DH_BOND_MAX            # Å  X–H covalent ceiling (assign H to donor)
_IM_AROM_O        = CFG.AROM_HB_DIST_O_ACC        # Å  aromatic H-bond, O donor
_IM_AROM_N        = CFG.AROM_HB_DIST_N_ACC        # Å  aromatic H-bond, N donor
_IM_AROM_DON_ANGLE = CFG.AROM_HB_DON_ANGLE_O      # °  aromatic H-bond donor minimum angle
_IM_VDW           = CFG.VDW_RADII                 # Bondi (1964) vdW radii (contact ratio)
_IM_VDW_DEF       = CFG.VDW_RADIUS_DEFAULT
_IM_CONTACT_GOOD  = CFG.CONTACT_RATIO_GOOD        # ratio ≥ this → no/weak overlap
_IM_CONTACT_BAD   = CFG.CONTACT_RATIO_BAD
_IM_CONTACT_UGLY  = CFG.CONTACT_RATIO_UGLY
# Protein aromatic ring atom names (π acceptor for aromatic H-bonds). HIS variants normalised.
_IM_RING_ATOMS = {
    "PHE": {"CG","CD1","CD2","CE1","CE2","CZ"},
    "TYR": {"CG","CD1","CD2","CE1","CE2","CZ"},
    "TRP": {"CG","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"},
    "HIS": {"CG","ND1","CD2","CE1","NE2"},
}
_IM_RES_COLORS = {
    "ASP":"#C0392B","GLU":"#C0392B",
    "ARG":"#2471A3","LYS":"#2471A3",
    "HIS":"#1E8449","TRP":"#7D3C98","TYR":"#BA4A00",
    "SER":"#117A65","THR":"#117A65","ASN":"#117A65","GLN":"#117A65",
    "PHE":"#6C3483",
    "LEU":"#626567","ILE":"#626567","VAL":"#626567",
    "ALA":"#626567","GLY":"#626567","PRO":"#626567",
    "MET":"#7E5109","CYS":"#7E5109",
}
_IM_ELEM_COLORS = {
    "C":"#2C3E50","N":"#1A5276","O":"#A93226",
    "F":"#1D8348","S":"#D4AC0D","other":"#717D7E",
}
_IM_HYDROPHOBIC = {"LEU","ILE","VAL","PHE","TRP","PRO","MET","ALA","GLY","CYS"}
_IM_AA3 = {
    "ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
    "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
    "UNK","MSE","SEC","PYL","HIE","HID","HIP","ASH","GLH","LYN",
    "CYX","CYM",
}


def _im_angle(a, b, c):
    """Angle (degrees) at vertex b for points a-b-c."""
    v1, v2 = a - b, c - b
    cosv = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9)
    return float(np.degrees(np.arccos(np.clip(cosv, -1.0, 1.0))))


def _im_don_acc(heavy, hpos):
    """Split heavy N/O atoms into donors (with bonded H) and acceptors.

    Returns (donors, acceptors) where donors = [(Dpos, elem, [Hpos,...])] and
    acceptors = [Dpos]. A heavy N/O is a donor when ≥1 H sits within
    CFG.HB_DH_BOND_MAX (covalent X–H); every heavy N/O is an acceptor candidate.
    """
    donors, acceptors = [], []
    harr = np.array(hpos) if hpos else None
    for nm, el, ps in heavy:
        if el not in ("N", "O"):
            continue
        acceptors.append(ps)
        if harr is not None and len(harr):
            d  = np.linalg.norm(harr - ps, axis=1)
            hs = [harr[i] for i in range(len(d)) if d[i] <= _IM_DH_BOND_MAX]
            if hs:
                donors.append((ps, el, hs))
    return donors, acceptors


def _im_ring(heavy_dict, resname):
    """Aromatic ring (centroid, normal) for an aromatic residue, else None."""
    rn  = CFG.PROTONATION_MAP.get(resname, resname)
    nms = _IM_RING_ATOMS.get(rn)
    if not nms:
        return None
    pts = np.array([heavy_dict[n] for n in nms if n in heavy_dict])
    if len(pts) < 5:
        return None
    cent = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - cent, full_matrices=False)
    return (cent, vt[2])


def _im_parse_pdb(path):
    """Return (lig_atoms, pro_residues, lig_hb) from a prepared (protonated) PDB.

    Hydrogens are RETAINED (PrepWizard adds them) so true Maestro H-bond geometry
    (H···A + donor angle) and aromatic H-bonds can be detected. `lig` holds heavy
    atoms only (for the 2D projection); ligand donor/acceptor-H geometry is
    returned separately as `lig_hb`. Protein donors/acceptors/rings are attached
    to each `pro` entry.

    Schrödinger prep places non-standard backbone residues as HETATM in chain L.
    Guard: exclude AA/UNK residues from ligand atoms even when chain == 'L'.
    """
    parser = _PDBParser(QUIET=True)
    struct = parser.get_structure("X", str(path))
    model  = struct[0]
    lig, pro = [], {}
    lig_heavy, lig_h = [], []
    for chain in model:
        for res in chain:
            hetflag = res.id[0].strip()
            resname = res.resname.strip()
            allat = [(a.name.strip(),
                      (a.element.strip().upper() or a.name.strip()[0].upper()),
                      a.coord.copy()) for a in res]
            heavy = [t for t in allat if t[1] != "H"]
            hyd   = [t[2] for t in allat if t[1] == "H"]
            is_lig_chain = chain.id == CFG.PREP_LIGAND_CHAIN and resname not in _IM_AA3
            is_hetatm    = bool(hetflag) and hetflag != "W" and resname not in _IM_AA3
            if is_lig_chain or is_hetatm:
                for nm, el, ps in heavy:
                    lig.append({"name": nm, "elem": el, "pos": ps})
                    lig_heavy.append((nm, el, ps))
                lig_h.extend(hyd)
            elif not hetflag or resname in _IM_AA3:
                if not heavy:
                    continue
                key = (chain.id, res.id[1], resname)
                heavy_dict = {nm: ps for nm, el, ps in heavy}
                elem_dict  = {nm: el for nm, el, ps in heavy}   # true parsed element per atom (ZN/SE/MG safe)
                donors, acceptors = _im_don_acc(heavy, hyd)
                pro[key] = {"resname": resname, "resnum": res.id[1],
                            "chain": chain.id, "atoms": heavy_dict, "elems": elem_dict,
                            "center": np.mean([ps for _, _, ps in heavy], axis=0),
                            "donors": donors, "acceptors": acceptors,
                            "ring": _im_ring(heavy_dict, resname)}
    l_don, l_acc = _im_don_acc(lig_heavy, lig_h)
    lig_hb = {"donors": l_don, "acceptors": l_acc}
    return lig, pro, lig_hb


def _im_hbond_check(res, lig_hb):
    """True Maestro H-bond (H···A ≤ cutoff, D–H···A ≥ angle), both directions.

    Returns the donor–acceptor heavy-atom distance of the best qualifying pair,
    or None. Uses CFG.THRESHOLD_HB_DIST_HA / THRESHOLD_HB_ANGLE_MIN.
    """
    best = None
    pairs = [(res["donors"], lig_hb["acceptors"]),     # protein donor → ligand acceptor
             (lig_hb["donors"], res["acceptors"])]     # ligand donor  → protein acceptor
    for donors, acceptors in pairs:
        for (Dpos, _el, Hs) in donors:
            for A in acceptors:
                for H in Hs:
                    if (np.linalg.norm(H - A) <= _IM_HB_HA_DIST and
                            _im_angle(Dpos, H, A) >= _IM_HB_DON_ANGLE):
                        dDA = float(np.linalg.norm(Dpos - A))
                        best = dDA if best is None else min(best, dDA)
    return best


def _im_arom_hbond_check(ring, lig_donors):
    """Aromatic H-bond: ligand donor-H → protein π-ring centroid (acceptor).

    Distance ceiling by donor element (CFG.AROM_HB_DIST_O_ACC / _N_ACC); donor
    angle ≥ CFG.AROM_HB_DON_ANGLE_O. Returns H···centroid distance or None.
    """
    cent, _normal = ring
    best = None
    for (Dpos, el, Hs) in lig_donors:
        ceil = _IM_AROM_O if el == "O" else _IM_AROM_N
        for H in Hs:
            d = float(np.linalg.norm(H - cent))
            if d <= ceil and _im_angle(Dpos, H, cent) >= _IM_AROM_DON_ANGLE:
                best = d if best is None else min(best, d)
    return best


def _im_find_contacts(lig, pro, lig_hb=None):
    """Detect interactions using the Maestro criteria defined in CFG §3.

    Priority per residue: H-bond > aromatic H-bond > salt bridge > steric
    contact. H-bonds and aromatic H-bonds use true H···A geometry (explicit H);
    salt bridges use heavy-atom charged-group distance; remaining proximal
    residues are classified as steric contacts with a Bondi vdW ratio quality
    tag (H-bonds and salt bridges are excluded from contacts per CFG flags).
    """
    lig_hb = lig_hb or {"donors": [], "acceptors": []}
    lig_pos = np.array([a["pos"] for a in lig])
    contacts = []
    for key, res in pro.items():
        prot_pos = np.array(list(res["atoms"].values()))
        prot_nms = list(res["atoms"].keys())
        diffs = lig_pos[:, None, :] - prot_pos[None, :, :]
        dists = np.linalg.norm(diffs, axis=-1)
        mind  = float(dists.min())
        if mind > _IM_CONTACT_DIST:
            continue
        li, pi = np.unravel_index(dists.argmin(), dists.shape)
        la, pname = lig[li], prot_nms[pi]

        itype, disp = None, mind
        quality = None

        # 1) True H-bond (Maestro H···A + donor angle)
        hb = _im_hbond_check(res, lig_hb)
        if hb is not None:
            itype, disp = "hbond", hb

        # 2) Aromatic H-bond (ligand donor-H → protein ring)
        if itype is None and res["ring"] is not None:
            ah = _im_arom_hbond_check(res["ring"], lig_hb["donors"])
            if ah is not None:
                itype, disp = "arom_hbond", ah

        # 3) Salt bridge (charged heavy groups within CFG cutoff)
        if itype is None:
            charged = pname[:2] in {"NH", "NZ", "NE", "ND", "OD", "OE"}
            if mind <= _IM_SALT_DIST and la["elem"] in {"O", "N"} and charged:
                itype, disp = "salt", mind

        # 4) Steric contact — Bondi vdW ratio (excludes H-bonds / salt above)
        if itype is None:
            ra = _IM_VDW.get(la["elem"], _IM_VDW_DEF)
            # True parsed element (ZN/SE/MG/FE safe); fall back to first letter only if absent.
            rb = _IM_VDW.get(res["elems"].get(pname, pname[0].upper()).upper(), _IM_VDW_DEF)
            ratio = mind / (ra + rb)
            quality = ("ugly" if ratio < _IM_CONTACT_UGLY else
                       "bad"  if ratio < _IM_CONTACT_BAD  else
                       "good" if ratio < _IM_CONTACT_GOOD else "far")
            itype = "contact"

        contacts.append({"key": key, "resname": res["resname"], "resnum": res["resnum"],
                         "chain": res["chain"], "dist": disp, "itype": itype,
                         "quality": quality,
                         "is_hbond": itype in ("hbond", "arom_hbond"),
                         "is_salt": itype == "salt",
                         "lig_atom": la, "prot_atom": pname,
                         "center": res["center"]})
    return sorted(contacts, key=lambda x: x["dist"])


def _im_project(lig, contacts):
    """SVD project onto ligand principal plane; place residues on fixed ring."""
    lig_c3 = np.mean([a["pos"] for a in lig], axis=0)
    L = np.array([a["pos"] - lig_c3 for a in lig])
    if len(L) >= 2:
        _, _, Vt = np.linalg.svd(L, full_matrices=False)
        u1, u2 = Vt[0], Vt[1]
    else:
        u1, u2 = np.array([1.,0.,0.]), np.array([0.,1.,0.])
    lig_2d = np.array([[np.dot(a["pos"]-lig_c3, u1),
                        np.dot(a["pos"]-lig_c3, u2)] for a in lig])
    span = np.max(np.linalg.norm(lig_2d, axis=1)) if len(lig_2d) else 1.0
    scale = 3.0 / max(span, 0.5)
    lig_2d *= scale
    ring_r = 5.5
    res_2d = {}
    for c in contacts:
        v = c["center"] - lig_c3
        x2 = np.dot(v, u1) * scale
        y2 = np.dot(v, u2) * scale
        d  = np.sqrt(x2**2 + y2**2)
        if d < 0.1: d = ring_r
        res_2d[c["key"]] = np.array([ring_r * x2/d, ring_r * y2/d])
    # Angular collision separation
    keys = list(res_2d.keys())
    for _ in range(60):
        moved = False
        for i in range(len(keys)):
            for j in range(i+1, len(keys)):
                a, b = res_2d[keys[i]], res_2d[keys[j]]
                d = np.linalg.norm(a - b)
                if d < 2.0:
                    push = (a - b) / max(d, 0.01) * (2.0 - d) * 0.5
                    res_2d[keys[i]] = a + push
                    res_2d[keys[j]] = b - push
                    for k in (keys[i], keys[j]):
                        r = np.linalg.norm(res_2d[k])
                        if r > 0.1:
                            res_2d[k] = res_2d[k] * ring_r / r
                    moved = True
        if not moved:
            break
    return lig_2d, res_2d


def _im_separate_atoms(lig_2d, min_sep=0.28):
    """Push overlapping ligand atoms apart (e.g. CF3 group after SVD projection)."""
    coords = lig_2d.copy().astype(float)
    n = len(coords)
    for _ in range(40):
        moved = False
        for i in range(n):
            for j in range(i + 1, n):
                d = np.linalg.norm(coords[i] - coords[j])
                if d < min_sep:
                    if d < 1e-6:
                        ang = 2 * np.pi * i / max(n, 1)
                        coords[i] += np.array([np.cos(ang), np.sin(ang)]) * min_sep * 0.7
                    else:
                        push = (coords[i] - coords[j]) / d * (min_sep - d) * 0.55
                        coords[i] += push
                        coords[j] -= push
                    moved = True
        if not moved:
            break
    return coords


def _draw_interaction_diagram(pdb_path, fig_root, base_name, C):
    """Render a publication-quality 2D protein–ligand interaction map (no PyMOL)."""
    out_dir = fig_root / "InteractionMap"
    out_dir.mkdir(exist_ok=True, parents=True)
    out_png = out_dir / f"{base_name}_interaction.png"
    try:
        lig, pro, lig_hb = _im_parse_pdb(pdb_path)
        if not lig:
            return "Failed"
        contacts = _im_find_contacts(lig, pro, lig_hb)
        if not contacts:
            return "Failed"
        """
        Detector already assigns c['itype'] (hbond / arom_hbond / salt / contact).
        Only refine a plain steric contact on a hydrophobic residue for display.
        """
        for c in contacts:
            if c["itype"] == "contact" and c["resname"] in _IM_HYDROPHOBIC:
                c["itype"] = "hydrophobic"
        lig_2d, res_2d = _im_project(lig, contacts)
        lig_2d = _im_separate_atoms(lig_2d)
        _im_render_diagram(lig_2d, lig, res_2d, contacts, out_png, mode="distance")
        return "Success"
    except Exception as _e:
        _fig_log(f"  [!] InteractionMap failed for {base_name}: {_e}")
        return "Failed"


# -------------------------------------------------------------------------------
# --- Figure engine: software manager ---
# -------------------------------------------------------------------------------
class SoftwareManager:
    """Handles verification and automated installation of visual tools."""

    def __init__(self):
        self.status = {}

    def check_all(self):
        _fig_log(f"\n{ConsoleColours.BOLD}>>> Checking Visualisation Engine Availability{ConsoleColours.ENDC}")

        """
        1. PyMOL — prefer the binary in the running conda env over the system one
        (system /usr/bin/pymol on Ubuntu 24.04 uses python3-pymol 2.5 which
        calls `from imp import find_module`; imp was removed in Python 3.12)
        """
        _conda_pymol = Path(sys.executable).parent / "pymol"
        if _conda_pymol.exists():
            self.status["PyMOL"] = str(_conda_pymol)
        else:
            self.status["PyMOL"] = shutil.which("pymol")
        if not self.status["PyMOL"]:
            self.status["PyMOL"] = self._attempt_install_pymol()

        # 2. PLIP — try binary first, then check if module is importable
        self.status["PLIP"] = shutil.which("plip")
        if not self.status["PLIP"]:
            # Try to import PLIP as a module; if successful, use python -m plip
            if self._check_plip_module():
                self.status["PLIP"] = f"{sys.executable}:plip_module"
            else:
                self.status["PLIP"] = self._attempt_install_plip()

        self._print_summary()
        return self.status

    def _check_plip_module(self):
        """Check if PLIP is installed as a Python module."""
        try:
            import plip  # noqa: F401 — availability probe; module presence is the signal
            return True
        except ImportError:
            return False

    def _attempt_install_pymol(self):
        _fig_log(f"  {ConsoleColours.WARNING}⚠{ConsoleColours.ENDC} PyMOL      : Missing. Attempting Conda install...")
        try:
            subprocess.run(["conda", "install", "-y", "-c", "conda-forge", "pymol-open-source"],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            p = shutil.which("pymol")
            if p:
                _fig_log(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} PyMOL      : Successfully installed.")
            return p
        except Exception:
            _fig_log(f"  {ConsoleColours.FAIL}✗{ConsoleColours.ENDC} PyMOL      : Conda install failed.")
            return None

    def _attempt_install_plip(self):
        _fig_log(f"  {ConsoleColours.WARNING}⚠{ConsoleColours.ENDC} PLIP       : Missing. Attempting Pip install...")
        try:
            subprocess.run([sys.executable, "-m", "pip", "install", "plip"],
                           check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            p = shutil.which("plip")
            if p:
                _fig_log(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} PLIP       : Successfully installed.")
            return p
        except Exception:
            _fig_log(f"  {ConsoleColours.FAIL}✗{ConsoleColours.ENDC} PLIP       : Pip install failed.")
            return None

    def _print_summary(self):
        _fig_log(f"\n{ConsoleColours.BOLD}--- Final Tool Status ---{ConsoleColours.ENDC}")
        missing = False
        for t, v in self.status.items():
            if v:
                _fig_log(f"  {ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC} {t.ljust(10)} : Available")
            else:
                _fig_log(f"  {ConsoleColours.WARNING}⚠{ConsoleColours.ENDC} {t.ljust(10)} : Unavailable (Skipping)")
                missing = True

        if missing:
            _fig_log(f"\n{ConsoleColours.BOLD}>>> Software Installation Guide for Missing Tools:{ConsoleColours.ENDC}")
            if not self.status.get("PyMOL"):
                _fig_log("  • PyMOL    : https://pymol.org/")
            if not self.status.get("PLIP"):
                _fig_log("  • PLIP     : pip install plip")
            _fig_log(f"\n{ConsoleColours.OKBLUE}Note: Missing tools will be gracefully skipped.{ConsoleColours.ENDC}")

    def get_ligand_info(self, pdb_path, smi_map=None):
        """Return (res_name, res_num, has_f) for the ligand in chain L.

        Skips AA/UNK residues that Schrödinger prep places in chain L.
        Falls back to filename-pattern lookup in smi_map when available.
        """
        res_name, res_num, has_f = "LIG", "1", False
        _found_lig = False
        # Try filename-based name lookup: e.g. '25_TFA' in stem → name='TFA'
        if smi_map:
            _parts = pdb_path.stem.split("_")
            for _i, _p in enumerate(_parts[:-1]):
                _candidate = f"{_p}_{_parts[_i + 1]}"
                if _candidate in smi_map:
                    res_name = _parts[_i + 1]   # e.g. 'TFA'
                    break
        try:
            with open(pdb_path, "r") as f:
                for line in f:
                    if (line.startswith("HETATM") or line.startswith("ATOM")) and len(line) > 21 and line[21] == "L":
                        _rn = line[17:20].strip()
                        if _rn not in _IM_AA3 and not _found_lig:
                            if res_name == "LIG":
                                res_name = _rn
                            res_num  = line[22:26].strip()
                            _found_lig = True
                        elem = line[76:78].strip() if len(line) > 76 else ""
                        if elem == "F" or " F " in line[60:] or line.rstrip().endswith(" F"):
                            has_f = True
        except Exception:
            pass
        return res_name, res_num, has_f

# -------------------------------------------------------------------------------
# --- Figure engine: phase-2 entry point ---
# -------------------------------------------------------------------------------
def run_figure_generation(run_dir: Path, ext_dir: Path):
    """
    Phase 2: Run PyMOL and PLIP rendering on all PDB files under ext_dir.
    Called at the end of main() after Phase 1 (extraction) has completed.
    """
    global _fig_logger, logger

    if not ext_dir.exists():
        console_info("Phase 2: Extraction directory not found — skipping figure generation.")
        return

    # Use the existing global logger for consolidation
    _fig_logger = logger

    # Render figures in-place on the canonical PDB folders (siblings of ext_dir);
    # controls already live there too, so no separate control pass is needed.
    prep_parent = ext_dir.parent
    dirs = [d for d in (prep_parent / "1_Converted_Raw_PDB",
                        prep_parent / "2_Prepared_PDBs") if d.exists()]

    load_metadata(ext_dir)

    # Load SMILES map for ligand name lookup (used in get_ligand_info)
    _input_data_dir = next(
        (run_dir / n for n in ["1_Boltz2_Production/1_Input_FASTA_and_SMILES",
                                "1_Boltz2_Production/1_Input_Data"]
         if (run_dir / n).exists()), None)
    _, _smi_map_fig, _, _ = load_reference_data(_input_data_dir) if _input_data_dir else ({}, {}, 0, 0)

    # Resolve figure constants
    C = _fig_constants()

    sw_mgr    = SoftwareManager()
    sw_status = sw_mgr.check_all()
    workers   = CFG.GLOBAL_MAX_WORKERS

    _fig_log(f"\n{ConsoleColours.BOLD}Multi-Engine Figure Generation Pipeline (Concurrency: {workers}){ConsoleColours.ENDC}")

    """
    Build separate task queues:
    - subprocess_tasks: PyMOL + PLIP (spawn external processes — thread-safe)
    - imap_tasks: InteractionMap (matplotlib Agg — NOT thread-safe, run serially)
    """
    subprocess_tasks = []
    imap_tasks       = []
    for d in dirs:
        pdbs = sorted(list(d.glob("*.pdb")))
        if not pdbs:
            continue

        fig_root = d / "Figures"
        log_dir  = fig_root / "Figure_Logs"
        log_dir.mkdir(parents=True, exist_ok=True)

        # Wipe stale tool subdirs so regeneration starts clean
        for _tool_subdir in ("PyMOL", "PLIP", "InteractionMap"):
            _td = fig_root / _tool_subdir
            if _td.exists():
                shutil.rmtree(_td)

        for pdb in pdbs:
            lig_name, lig_num, has_f = sw_mgr.get_ligand_info(pdb, smi_map=_smi_map_fig)
            task_base = (pdb, fig_root, log_dir, pdb.stem)
            subprocess_tasks.append(("PyMOL", _run_pymol, task_base + (lig_name, lig_num, has_f, sw_status, C)))
            subprocess_tasks.append(("PLIP",  _run_plip,  task_base + (sw_status, C)))
            imap_tasks.append((pdb.stem, pdb, fig_root))

    _plip_tasks  = [t for t in subprocess_tasks if t[0] == "PLIP"]
    _pymol_tasks = [t for t in subprocess_tasks if t[0] == "PyMOL"]
    """
    PyMOL workers capped at 2: ray-tracing is CPU-bound; high concurrency causes
    CPU saturation, massively inflated render times, and timeout-related orphan processes.
    """
    _pymol_workers = min(2, max(1, workers))
    _fig_log(
        f"Subprocess tasks: PLIP={len(_plip_tasks)} ({workers} workers) | "
        f"PyMOL={len(_pymol_tasks)} ({_pymol_workers} workers)  |  "
        f"InteractionMap={len(imap_tasks)} (serial)"
    )

    results    = {}
    completed  = 0

    def _collect(fut_map, pool_label):
        nonlocal completed
        for f in as_completed(fut_map):
            engine, cname = fut_map[f]
            status = f.result()
            if cname not in results:
                results[cname] = {}
            results[cname][engine] = status
            completed += 1
            total_sub = len(subprocess_tasks)
            if completed % 5 == 0 or completed == total_sub:
                if sys.stdout.isatty():
                    print(f"\r  [{pool_label}] {completed}/{total_sub} tasks ...", end="", flush=True)
                elif completed == total_sub:
                    print(f"  [{pool_label}] {completed}/{total_sub} tasks finished.", flush=True)

    # Phase A1: PLIP — I/O-bound, run at full worker count
    with ThreadPoolExecutor(max_workers=workers) as executor:
        _collect({executor.submit(t[1], *t[2]): (t[0], t[2][3]) for t in _plip_tasks}, "PLIP")

    # Phase A2: PyMOL — CPU-bound ray-tracing, limited to 2 concurrent processes
    with ThreadPoolExecutor(max_workers=_pymol_workers) as executor:
        _collect({executor.submit(t[1], *t[2]): (t[0], t[2][3]) for t in _pymol_tasks}, "PyMOL")

    if sys.stdout.isatty():
        print("\n")

    # Phase B: InteractionMap — serial (matplotlib Agg is not thread-safe)
    _fig_log(f"\n  Running {len(imap_tasks)} InteractionMap diagrams (serial)...")
    for idx, (cname, pdb, fig_root) in enumerate(imap_tasks, 1):
        status = _draw_interaction_diagram(pdb, fig_root, cname, C)
        if cname not in results:
            results[cname] = {}
        results[cname]["InteractionMap"] = status
        if sys.stdout.isatty():
            print(f"\r  [InteractionMap] {idx}/{len(imap_tasks)}", end="", flush=True)
    if sys.stdout.isatty():
        print("\n")

    # Formal reporting
    for d in dirs:
        pdbs = [p.stem for p in sorted(list(d.glob("*.pdb")))]
        if not pdbs:
            continue

        _fig_log(f"\n{ConsoleColours.BOLD}>>> Folder: {d.relative_to(run_dir)}{ConsoleColours.ENDC}")
        _fig_log(f"  {'Complex Name':<40} | PyMOL | PLIP | IntMap")
        _fig_log(f"  {'-'*70}")

        for cn in pdbs:
            r = results.get(cn, {})
            def icon(s):
                if s == "Success":
                    return f"{ConsoleColours.OKGREEN}✔{ConsoleColours.ENDC}"
                return f"{ConsoleColours.OKBLUE}-{ConsoleColours.ENDC}"

            _fig_log(f"  {cn[:37]+'..' if len(cn)>39 else cn:<40} |   {icon(r.get('PyMOL'))}   |  {icon(r.get('PLIP'))}  |    {icon(r.get('InteractionMap'))}")

    # Final cleanup — remove tool-side non-image files and Figure_Logs
    _keep_ext = {".png", ".pse", ".jpg", ".jpeg", ".svg", ".html"}
    for d in dirs:
        fig_dir = d / "Figures"
        if not fig_dir.exists():
            continue
        for tool_dir in fig_dir.iterdir():
            if not tool_dir.is_dir() or tool_dir.name == "Figure_Logs":
                continue
            for item in list(tool_dir.iterdir()):
                if item.is_file() and item.suffix.lower() not in _keep_ext:
                    item.unlink(missing_ok=True)
        logs = fig_dir / "Figure_Logs"
        if logs.exists():
            shutil.rmtree(logs, ignore_errors=True)

    _fig_log(f"\n\n{ConsoleColours.BOLD}{ConsoleColours.OKGREEN}                ✔ Figure Generation Complete                {ConsoleColours.ENDC}\n")




def setup_logging_extraction(output_dir: Path) -> Path:
    """Initialise the Top-N extraction log (separate handler from the prep log)."""
    global logger
    log_file = output_dir / "00_TopN_and_Preparation_Log.txt"
    logger = _utils_mod.setup_logging(log_file, logger_name="extraction", mode="a", timestamp=True)
    return log_file


def prep_and_convert_phase(args):
    """Phase 1 - CIF->PDB conversion + PrepWizard preparation of the MD-ready cohort."""

    root = Path.cwd()
    run_path = root / args.run_folder_name
    if not run_path.exists(): sys.exit(f"Run path missing: {run_path}")

    # Paths — one consolidated Step-05 folder for both phases (prep + extraction)
    analysis_dir   = run_path / "5_TopN_and_Preparation"
    dir_raw        = analysis_dir / "1_Converted_Raw_PDB"
    dir_prep_clean = analysis_dir / "2_Prepared_PDBs"
    prod_dir       = run_path / "1_Boltz2_Production"

    # Init output directories and log (log goes directly in analysis_dir)
    for d in [dir_raw, dir_prep_clean]: d.mkdir(parents=True, exist_ok=True)
    setup_logging(analysis_dir)

    best_cifs_dir = run_path / "2_Best_Complexes_CIFs"
    if not best_cifs_dir.exists(): sys.exit(f"Critical: 2_Best_Complexes_CIFs directory not found in {run_path}")
    jobs = collect_best_cifs(best_cifs_dir)
    _n_all_cifs = len(jobs)

    # SECTION 18 gate: restrict to the MD-ready cohort so only ~10 complexes are
    # converted and prepared, not the whole predicted library. Control jobs
    # (ID 0000000_*) are always kept alongside the cohort so their raw/prepared
    # PDBs exist for the downstream control comparison. Falls back to all CIFs
    # when the ranked CSV carries no MD_Selected column.
    _md_jobs = load_md_selected_jobs(prod_dir)
    if _md_jobs:
        jobs = [(jn, cif) for (jn, cif) in jobs
                if jn in _md_jobs or jn.startswith("0000000")]

    _utils_mod.print_script_banner(
        "05_TopN_and_PDB_Preparation_FAcDs.py",
        "CIF → PDB Conversion  ·  PrepWizard Refinement  ·  Chain Assignment",
    )
    console_info(f"  Run Name : {args.run_folder_name}")
    if _md_jobs:
        console_info(f"MD-ready gate: preparing {len(jobs)} of {_n_all_cifs} Best Complex CIFs (MD_Selected)")
    else:
        console_info(f"Found {len(jobs)} Best Complex CIFs in 2_Best_Complexes_CIFs")
    if args.quick:
        console_info(f"{ConsoleColours.OKBLUE}Quick Resume Mode: Deep validation disabled.{ConsoleColours.ENDC}")

    rank_map = load_rank_map(prod_dir)
    cat_anchor_map = load_catalytic_anchor_map(prod_dir)   # 02 dynamic-alignment anchors for the identity guard

    print("")
    console_info(f"Outputs      : {analysis_dir}")
    console_info(f"Schrödinger  : {SCHRODINGER_PATH}")
    console_info(f"PrepWizard   : {'Found' if PREPWIZARD_BIN.exists() else 'NOT FOUND — preparation will be skipped'}")
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
        _raw_failures = []
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
                    if res["status"] in ("Success", "Renamed", "Skipped"):
                        final_valid_raw_names.append(res["job"])
                    else:
                        # Conversion failures are otherwise silent — record them so
                        # an all-fail run is surfaced rather than queuing zero
                        # PrepWizard jobs and exiting 0.
                        _raw_failures.append((res["job"], res.get("reason", res["status"])))
                    if count % 100 == 0 or count == total_gen:
                        if sys.stdout.isatty():
                            print(f"\rGenerated: {count}/{total_gen}", end="", flush=True)
                        elif count == total_gen:
                            print(f"Generated: {count}/{total_gen}", flush=True)
        if sys.stdout.isatty():
            print("")
        if _raw_failures:
            console_info(f"Raw PDB conversion failures: {len(_raw_failures)}")
            for _fjob, _freason in _raw_failures[:20]:
                console_info(f"  ! {_fjob}: {_freason}")
            if len(_raw_failures) > 20:
                console_info(f"  ! ... and {len(_raw_failures) - 20} more")

    console_info(f"Raw PDB generation stage completed in {time.time() - t0:.1f}s")
    console_separator()

    # Abort if every input CIF failed to convert — downstream stages require at
    # least one raw PDB, and a silent exit-0 here masks a total conversion failure.
    if jobs and not final_valid_raw_names:
        print("CRITICAL: all CIF->PDB conversions failed; no raw PDBs produced. "
              "Aborting before PrepWizard.", file=sys.stderr, flush=True)
        sys.exit(1)

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
    _prep_ok = 0
    _prep_fail = 0

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
                    stat_text = res["status"]
                    if stat_text != "Success": symbol = "✗"
                    try:
                        job_id_num = res["job"].split("_")[0]
                    except Exception:
                        job_id_num = "?"
                    print(f"({count}/{total_prep} | ID:{job_id_num}) {symbol} {res['job']} | [{stat_text}]", flush=True)
                    if logger: logger.info(f"{symbol} {res['job']} | [{stat_text}]")
                    # Residue identity guard: run after every successful PrepWizard job
                    if stat_text == "Success":
                        _prep_ok += 1
                        _prepared_pdb = dir_prep_clean / f"{res['job']}_Prepared.pdb"
                        _check_residue_identity_guard(_prepared_pdb, res["job"], CFG,
                                                      anchors=cat_anchor_map.get(res["job"]))
                    else:
                        _prep_fail += 1

    print("")
    console_separator()

    n_prep_cached = prep_already_done_count
    n_prep_new    = len(prep_jobs_to_run)
    _print_run_delta_table(n_raw_cached, n_raw_new, n_prep_cached, n_prep_new)

    '''
    Surface PrepWizard failures so they are not silently masked. A total failure
    (no structure prepared this run) aborts with a non-zero status so the pipeline
    runner halts before Steps 07/08 try to run on missing structures; partial
    failures are reported but allowed through (per-structure resilience).
    '''
    if _prep_fail:
        console_info(f"  [WARNING] PrepWizard failed for {_prep_fail}/{total_prep} structure(s) — see *_prepwizard_FAILED.log.")
        if _prep_ok == 0:
            console_info("  [FATAL] No structures were prepared this run — aborting (downstream MD/QM would run on missing files).")
            sys.exit(1)




def topn_extraction_phase(args):
    """Phase 2 - extract, validate and render the MD-ready complexes for handover."""

    # -------------------------------------------------------------------------------
    # Step 5.3: Path validation
    # -------------------------------------------------------------------------------
    run_dir = DEFAULT_BASE_PATH / args.run_folder_name
    if not run_dir.exists():
        # Fallback to check relative path if not in default base
        run_dir = Path.cwd() / args.run_folder_name
        if not run_dir.exists():
            console_info(f"Error: Run folder not found at: {run_dir}")
            sys.exit(1)

    prod_dir = run_dir / "1_Boltz2_Production"           # Step 02 Output
    prep_dir = run_dir / "5_TopN_and_Preparation"        # Step 05 consolidated folder
    input_data_dir = next((prod_dir / n for n in ["1_Input_FASTA_and_SMILES", "1_Input_Data"] if (prod_dir / n).exists()), prod_dir / "1_Input_FASTA_and_SMILES")

    # Comparison + handover outputs (Ramachandran, Controls, handover files, combined
    # CSV) live under 3_Comparative_Analysis/. The raw and prepared PDBs themselves
    # are NOT duplicated here — they stay in the sibling 1_Converted_Raw_PDB/ and
    # 2_Prepared_PDBs/ folders, where their figures are also rendered.
    final_dir = prep_dir / "3_Comparative_Analysis"
    final_dir.mkdir(parents=True, exist_ok=True)
    setup_logging_extraction(prep_dir)   # append to the single Step-05 log in the parent

    raw_pdb_dir  = prep_dir / "1_Converted_Raw_PDB"
    prep_pdb_dir = prep_dir / "2_Prepared_PDBs"

    # -------------------------------------------------------------------------------
    # Step 5.4: Header output
    # -------------------------------------------------------------------------------
    _utils_mod.print_script_banner(
        "05_TopN_and_PDB_Preparation_FAcDs.py",
        "Top-N Delivery  ·  Tier Extraction  ·  Ramachandran Validation  ·  PyMOL/PLIP Figure Generation",
    )
    console_info(f"  Run Name : {args.run_folder_name}")
    console_info(f"  Output   : {final_dir}")
    console_separator()

    # -------------------------------------------------------------------------------
    # Step 5.5: Load ranking logic
    # -------------------------------------------------------------------------------
    # Primary: FAcDs Ranked CSV written by 02_Production_FAcDs.py
    rank_csvs = (sorted(prod_dir.glob("7_Boltz2_FAcDs_Ranked_*.csv")) or
                 sorted(prod_dir.glob("*_Ranked_*.csv")))
    if not rank_csvs:
        console_info("Error: No Ranked CSV found in production directory.")
        sys.exit(1)

    rank_csv = rank_csvs[-1]
    console_info(f"Loading Logic : {rank_csv.name}")

    try:
        df = pd.read_csv(rank_csv, low_memory=False)
    except Exception as e:
        console_info(f"Error reading CSV file: {e}")
        sys.exit(1)

    # Filter for successful jobs
    if "status" in df.columns:
        df = df[df["status"] == "Success"]

    # Sort by Rank
    if "Scientific_Rank" in df.columns:
        df = df.sort_values("Scientific_Rank")
    elif "Production_Rank" in df.columns:
        df = df.sort_values("Production_Rank")
    else:
        console_info("Warning: Rank column missing. Using default sort.")

    # -------------------------------------------------------------------------------
    # Step 5.6: Load reference data (FASTA/SMILES)
    # -------------------------------------------------------------------------------
    console_info("Loading reference sequences and SMILES...")
    seq_map, smi_map, fasta_count, smi_count = load_reference_data(input_data_dir)
    console_info(f"Loaded {fasta_count} Sequences, {smi_count} SMILES.")

    # -------------------------------------------------------------------------------
    # Step 5.7: Tier / MD-ready selection
    # -------------------------------------------------------------------------------
    import select as _select

    TIER_ORDER = CFG.TIER_ORDER

    # Separate controls from candidates (controls always extracted independently)
    def _is_control(row):
        jn = str(row.get("job_name", ""))
        return jn.startswith("0000000")

    df_controls  = df[df.apply(_is_control, axis=1)].copy()
    df_candidates = df[~df.apply(_is_control, axis=1)].copy()

    # SECTION 18 gate: when the ranked CSV carries MD_Selected, extract exactly that
    # cohort and skip the interactive tier prompt (non-interactive and reproducible,
    # matching the CIF->PDB preparation stage).
    _md_col = getattr(CFG, "MD_SELECTED_COL", "MD_Selected")
    _md_subset = None
    if _md_col in df_candidates.columns:
        _mask = df_candidates[_md_col].astype(str).str.lower().isin(["true", "1", "1.0"])
        if _mask.any():
            _md_subset = df_candidates[_mask].copy()

    # Build tier availability table from candidates only
    tier_col = next((c for c in CFG.VIS_PHYLO_COLUMN_MAP["Tier"] if c in df_candidates.columns), None)
    available_tiers = []
    if tier_col:
        # Standardise tiers using utility logic
        df_candidates = _utils_mod.standardise_dataframe_tiers(df_candidates, CFG)
        _tier_series = df_candidates[CFG.COL_TIER]
        counts = _tier_series.value_counts()
        for t in TIER_ORDER:
            if t in counts.index and counts[t] > 0:
                available_tiers.append((t, int(counts[t])))

    # Dynamically count actual files on disk
    actual_raw_count  = len(list(raw_pdb_dir.glob("*.pdb")))  if raw_pdb_dir.exists()  else 0
    actual_prep_count = len(list(prep_pdb_dir.glob("*.pdb"))) if prep_pdb_dir.exists() else 0
    console_info(f"\n  Total Raw Structures Available      : {actual_raw_count}")
    console_info(f"  Total Prepared Structures Available : {actual_prep_count}")
    console_info(f"  Total Candidate Jobs                : {len(df_candidates)}")
    console_info(f"  Control Jobs (auto-extracted)       : {len(df_controls)}")
    console_separator()

    if _md_subset is not None:
        subset = _md_subset
        selected_tiers = ["MD_Selected"]
        console_info(f"  MD-ready gate: extracting {len(subset)} MD_Selected complexes (tier prompt skipped)")
    elif available_tiers:
        print("\n  Select Degrader Tiers to Extract")
        print(  "  ┌───────┬──────────────────┬─────────────┐")
        print(  "  │  No.  │  Tier            │     Count   │")
        print(  "  ├───────┼──────────────────┼─────────────┤")
        for idx, (tier, count) in enumerate(available_tiers, 1):
            print(f"  │  [{idx:2d}] │  {tier:<16}│  {count:>9,}  │")
        print(  "  └───────┴──────────────────┴─────────────┘")
        print(  "\n  Enter tier numbers separated by commas (e.g. 1,2,3).")
        print(  f"  Auto-selecting top tier [{available_tiers[0][0]}] in {CFG.VIS_TOP_N_AUTO_SELECT_TIMEOUT} seconds...\n")
        sys.stdout.flush()

        selected_indices = []
        try:
            rlist, _, _ = _select.select([sys.stdin], [], [], CFG.VIS_TOP_N_AUTO_SELECT_TIMEOUT)
            if rlist:
                raw_input = sys.stdin.readline().strip()
                for part in raw_input.split(","):
                    try:
                        idx = int(part.strip()) - 1
                        if 0 <= idx < len(available_tiers):
                            selected_indices.append(idx)
                    except ValueError:
                        pass
        except Exception:
            pass

        if not selected_indices:
            selected_indices = [0]
            console_info(f"  -> No selection — auto-selecting: [{available_tiers[0][0]}]")
        else:
            console_info(f"  -> User selected tier(s): {', '.join(available_tiers[i][0] for i in selected_indices)}")

        selected_tiers = [available_tiers[i][0] for i in selected_indices]
        subset = df_candidates[df_candidates[tier_col].isin(selected_tiers)].copy()
    else:
        # No tier column or no data — fall back to top-N by rank
        console_info("  Warning: No tier data found. Falling back to top-50 by rank.")
        top_n_fallback = min(args.top or 50, len(df_candidates))
        subset = df_candidates.head(top_n_fallback).copy()
        selected_tiers = ["All"]

    top_n = len(subset)
    tier_label = "_".join(selected_tiers)
    console_info(f"  Extracting {top_n} structures from tier(s): {', '.join(selected_tiers)}")
    console_separator()

    # -------------------------------------------------------------------------------
    # Step 5.8: Prepare output folders
    # -------------------------------------------------------------------------------
    folder_tag = f"{tier_label}_{top_n}hits" if top_n else "Selected"

    # Comparative Ramachandran plots (raw vs prepared) for the selected cohort.
    out_rama = final_dir / f"1_{folder_tag}_Comparative_Ramachandran_Plots"
    out_rama.mkdir(parents=True, exist_ok=True)

    # Control output folders (always created)
    out_ctrl = final_dir / "2_Controls"
    out_ctrl_raw  = out_ctrl / "1_Raw"
    out_ctrl_prep = out_ctrl / "2_Prepared"
    out_ctrl_raw.mkdir(parents=True, exist_ok=True)
    out_ctrl_prep.mkdir(parents=True, exist_ok=True)

    # Molecular Handover Folder
    out_handover = final_dir / f"3_{folder_tag}_Molecular_Handover_Files"
    out_handover.mkdir(parents=True, exist_ok=True)

    # Initialise SDF Writer
    sdf_path = out_handover / f"{folder_tag}_Ligands.sdf"
    sdf_writer = Chem.SDWriter(str(sdf_path))
    count_sdf = 0
    try:

        # -------------------------------------------------------------------------------
        # Step 5.8.1: Control-case extraction (always automatic)
        # -------------------------------------------------------------------------------
        ctrl_extracted_raw  = 0
        ctrl_extracted_prep = 0
        ctrl_csv_rows = []

        if not df_controls.empty:
            console_info(f"Extracting {len(df_controls)} control structure(s) → Controls/")
            for _, crow in df_controls.iterrows():
                ctrl_name  = str(crow.get("job_name", ""))
                # Control jobs: 0000000_02–3 (DeHa4) / 0000000_4–6 (3R3U)
                fname_raw  = f"{ctrl_name}_RAW.pdb"
                fname_prep = f"{ctrl_name}_Prepared.pdb"

                src_raw  = raw_pdb_dir  / fname_raw
                src_prep = prep_pdb_dir / fname_prep

                if src_raw.exists():
                    shutil.copy2(src_raw,  out_ctrl_raw  / fname_raw)
                    ctrl_extracted_raw += 1
                else:
                    console_info(f"    {ConsoleColours.WARNING}⚠ Missing Raw Control: {ctrl_name}{ConsoleColours.ENDC}")

                if src_prep.exists():
                    shutil.copy2(src_prep, out_ctrl_prep / fname_prep)
                    ctrl_extracted_prep += 1
                else:
                    console_info(f"    {ConsoleColours.WARNING}⚠ Missing Prepared Control: {ctrl_name}{ConsoleColours.ENDC}")

                ctrl_csv_rows.append(crow)

            # Write controls CSV
            if ctrl_csv_rows:
                ctrl_csv_path = out_ctrl / "Controls_Scientific_Data.csv"
                pd.DataFrame(ctrl_csv_rows).to_csv(ctrl_csv_path, index=False)

            console_info(f"  -> Controls: {ctrl_extracted_raw} raw, {ctrl_extracted_prep} prepared structures saved.")
        else:
            console_info("  Note: No control jobs found in the CSV (FAcDs_Control / 3R3U_Control).")
        console_separator()

        # -------------------------------------------------------------------------------
        # Step 5.9: Extraction loop with unique aggregation
        # -------------------------------------------------------------------------------
        extracted_raw_count = 0
        extracted_prep_count = 0
        extracted_fa = []

        # Dictionaries for ligand aggregation
        ligand_ranks = defaultdict(list)  # Key: SMILES,    Value: list of ranks
        ligand_data  = {}                 # Key: SMILES,    Value: (LigandName, Mol)
        extracted_smi_lines = []

        # Dictionaries for protein aggregation (mirrors ligand logic)
        protein_ranks = defaultdict(list) # Key: sequence,  Value: list of ranks
        protein_data  = {}                # Key: sequence,  Value: protein_name

        for i, row in subset.iterrows():
            # Metadata
            rank = row.get("Scientific_Rank", i+1)
            name = row["job_name"]

            # Determine Protein and Ligand names for Lookup (column aliases sourced from CFG)
            p_name = next((str(row[c]) for c in CFG.VIS_PHYLO_COLUMN_MAP["Protein_Name"]
                           if c in row.index and pd.notna(row[c])), "Unknown")
            l_name = next((str(row[c]) for c in CFG.VIS_PHYLO_COLUMN_MAP["Ligand_Name"]
                           if c in row.index and pd.notna(row[c])), "Unknown")

            # File names
            fname_raw = f"{name}_RAW.pdb"
            fname_prep = f"{name}_Prepared.pdb"

            # Parse critical residues from Active_Site_Triad_Map
            triad_map = str(row.get("Active_Site_Triad_Map", ""))
            critical_res = {}
            if triad_map and triad_map != "nan":
                for item in triad_map.split("|"):
                    item = item.strip()
                    if ":" in item:
                        role, res = item.split(":", 1)
                        res = res.strip()
                        resname = res[:3]
                        _m = re.search(r"\d+", res)
                        resnum_str = _m.group() if _m else ""
                        if resnum_str:
                            critical_res[int(resnum_str)] = (resname, role.strip())

            # --- A. Ramachandran comparison (read canonical PDBs in-place; no copies) ---
            # The raw and prepared PDBs already live in 1_Converted_Raw_PDB/ and
            # 2_Prepared_PDBs/ — they are read directly rather than duplicated here.
            src_raw = raw_pdb_dir / fname_raw
            if src_raw.exists():
                extracted_raw_count += 1

            src_prep = prep_pdb_dir / fname_prep
            if src_prep.exists():
                extracted_prep_count += 1
                try:
                    angles_raw = []
                    angles_prep = []
                    if src_raw.exists():
                        st_raw = gemmi.read_structure(str(src_raw))
                        angles_raw = compute_ramachandran_angles(st_raw)
                    if src_prep.exists():
                        st_prep = gemmi.read_structure(str(src_prep))
                        angles_prep = compute_ramachandran_angles(st_prep)

                    if angles_raw or angles_prep:
                        rama_path = out_rama / f"Rank_{rank}_{fname_prep.replace('.pdb', '_Rama_Comparison.png')}"
                        save_ramachandran_comparison(
                            angles_raw, angles_prep,
                            f"Rank {rank}: {p_name} (Raw)",
                            f"Rank {rank}: {p_name} (Prepared)",
                            rama_path,
                            critical_res=critical_res,
                            dpi=CFG.VIS_FIGURE_DPI
                        )
                except Exception as e:
                    console_info(f"    ! Ramachandran generation failed for {name}: {e}")

            # --- B. Ligand Data Collection (RAW source preferred) ---
            # Prioritise Raw for the structure data if available
            mol_source = src_raw if src_raw.exists() else (src_prep if src_prep.exists() else None)

            smi = row.get("smiles", row.get("Ligand_Smile", smi_map.get(l_name, "SMILES_Not_Found")))

            if smi and smi != "SMILES_Not_Found":
                # Store Rank
                ligand_ranks[smi].append(rank)

                """
                Store Mol if this SMILES has not been seen yet
                Only attempt extraction if a valid source file exists
                """
                if smi not in ligand_data and mol_source:
                    mol = extract_chain_l_mol(mol_source)
                    ligand_data[smi] = (l_name, mol)

            # --- C. Extract Sequence (individual entry + unique-protein aggregation) ---
            seq = (row.get("target_sequence") or
                   row.get("sequence") or
                   row.get("Protein_Sequence") or
                   seq_map.get(p_name) or
                   "Sequence_Not_Found")
            # Full per-complex FASTA (one entry per rank — always written)
            extracted_fa.append(f">Rank_{rank}_{name} | {p_name}\n{seq}")
            # Unique-protein aggregation: key on sequence string; fall back to name when absent
            seq_key = seq if (seq and seq != "Sequence_Not_Found") else f"NAME:{p_name}"
            protein_ranks[seq_key].append(rank)
            if seq_key not in protein_data:
                protein_data[seq_key] = p_name

        # -------------------------------------------------------------------------------
        # Step 5.10: Unique ligand writing (SDF & SMILES)
        # -------------------------------------------------------------------------------

        # Iterate through unique SMILES found
        for smi, ranks in ligand_ranks.items():
            # Sort ranks numerically
            sorted_ranks = sorted(list(set(ranks)))

            """
            Construct Composite Name: Rank-2-4-5___LigandName
            1. Join ranks with hyphens
            """
            rank_str = "-".join(map(str, sorted_ranks))
            l_name, mol = ligand_data.get(smi, ("Unknown", None))

            # 2. Create formatted string
            composite_name = f"Rank-{rank_str}___{l_name}"

            # 1. Write SMILES
            extracted_smi_lines.append(f"{smi}\t{composite_name}")

            # 2. Write SDF (if 3D mol exists)
            if mol:
                mol.SetProp("_Name", composite_name)
                mol.SetProp("Ranks", rank_str)
                mol.SetProp("Ligand_Name", l_name)
                sdf_writer.write(mol)
                count_sdf += 1

        # -------------------------------------------------------------------------------
        # Step 5.11: Save data & summary
        # -------------------------------------------------------------------------------

        # Save CSV Data (hits + controls merged; controls appended with is_control flag)
        subset_csv_path = final_dir / f"4_{folder_tag}_Combined_Scientific_Data.csv"
        _ctrl_df = pd.DataFrame(ctrl_csv_rows) if ctrl_csv_rows else pd.DataFrame()
        if not _ctrl_df.empty:
            _ctrl_df = _ctrl_df.reindex(columns=subset.columns)
        _combined = pd.concat([subset, _ctrl_df], ignore_index=True)
        _combined.to_csv(subset_csv_path, index=False)

        # Save full per-complex FASTA (one entry per rank)
        fa_path = out_handover / f"{folder_tag}_All_Sequences.fasta"
        with open(fa_path, "w") as f:
            f.write("\n".join(extracted_fa))

        # Save SMILES
        smi_path = out_handover / f"{folder_tag}_Ligands.smi"
        with open(smi_path, "w") as f:
            f.write("\n".join(extracted_smi_lines))
    finally:
        sdf_writer.close()

    # Build and save unique-proteins FASTA (mirrors ligand composite naming)
    unique_fa_lines = []
    for seq_key, ranks in protein_ranks.items():
        sorted_ranks = sorted(list(set(ranks)))
        rank_str     = "-".join(map(str, sorted_ranks))
        p_name_u     = protein_data[seq_key]
        composite    = f"Rank-{rank_str}___{p_name_u}"
        actual_seq   = seq_key if not seq_key.startswith("NAME:") else "Sequence_Not_Found"
        unique_fa_lines.append(f">{composite}\n{actual_seq}")

    unique_fa_path = out_handover / f"{folder_tag}_Unique_Proteins.fasta"
    with open(unique_fa_path, "w") as f:
        f.write("\n".join(unique_fa_lines))

    # --- Delivery Package Summary ---
    _dp_rows = [
        ("1", "Raw PDBs + Figures",            str(raw_pdb_dir.resolve())),
        ("2", "Prepared PDBs + Figures",       str(prep_pdb_dir.resolve())),
        ("3", "Comparative Ramachandran Plots", str(out_rama.resolve())),
        ("4", "Controls (Raw/Prep/CSV)",       str(out_ctrl.resolve())),
        ("5", "Molecular Handover (FASTA/SDF)", str(out_handover.resolve())),
        ("6", "Combined Scientific Data (CSV)", str(subset_csv_path.resolve())),
        ("7", "Step-05 Log (prep + extraction)", str((prep_dir / "00_TopN_and_Preparation_Log.txt").resolve())),
    ]
    _lw = max(max(len(r[1]) for r in _dp_rows), max(len(r[2]) for r in _dp_rows))
    console_info(f"\n{SEPARATOR_LIGHT}")
    console_info(f"  Delivery Package  │  Tiers: {', '.join(selected_tiers)}")
    console_info(SEPARATOR_LIGHT)
    console_info(f"  ┌───┬{'─'*(_lw+2)}┐")
    for _n, _label, _path in _dp_rows:
        console_info(f"  │ {_n} │  {_label:<{_lw}}│")
        console_info(f"  │   │  {_path:<{_lw}}│")
    console_info(f"  └───┴{'─'*(_lw+2)}┘")

    # Detailed extraction stats
    _stat_rows = [
        ("Selected-tier raw PDBs",    extracted_raw_count),
        ("Selected-tier prepared PDBs", extracted_prep_count),
        ("Sequences in FASTA",        len(extracted_fa)),
        ("Unique proteins",           len(unique_fa_lines)),
        ("Unique ligands",            len(extracted_smi_lines)),
        ("CSV shape",                 f"{subset.shape[0]:,} rows × {subset.shape[1]} cols"),
        ("Control raw PDBs",          ctrl_extracted_raw),
        ("Control prepared PDBs",     ctrl_extracted_prep),
    ]
    _sk = max(len(r[0]) for r in _stat_rows)
    _sv = max(len(str(r[1])) for r in _stat_rows)
    console_info("")
    console_info(f"  ┌{'─'*(_sk+2)}┬{'─'*(_sv+4)}┐")
    for _k, _v in _stat_rows:
        console_info(f"  │  {_k:<{_sk}}│  {str(_v):>{_sv}}  │")
    console_info(f"  └{'─'*(_sk+2)}┴{'─'*(_sv+4)}┘")

    # ===============================================================================
    # Phase 2: Figure Generation
    # ===============================================================================
    run_figure_generation(run_dir, final_dir)


# ===============================================================================
# SECTION 6: PHASE 2 — FIGURE GENERATION (PyMOL / PLIP)
# ===============================================================================

# -------------------------------------------------------------------------------
# --- Figure engine: constants (sourced from CFG) ---
# -------------------------------------------------------------------------------
"""
These are resolved lazily (inside run_figure_generation) so that CFG is
available; they are module-level names for clarity.
"""


def main():
    parser = argparse.ArgumentParser(
        description="Merged Top-N selection + CIF->PDB generation/preparation (MD-ready cohort only)")
    parser.add_argument("run_folder_name", help="Run Folder Name (e.g. Boltz-2_Run_...)")
    parser.add_argument("--quick", action="store_true", help="Quick resume: skip deep validation")
    parser.add_argument("--top", type=int, default=None, help="Fallback top-N when no MD_Selected column")
    args = parser.parse_args()
    prep_and_convert_phase(args)
    topn_extraction_phase(args)


if __name__ == "__main__":
    main()
