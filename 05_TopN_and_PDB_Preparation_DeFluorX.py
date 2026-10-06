#!/usr/bin/env python3

"""
===============================================================================
DeFluorX Pipeline  |  Step 05  |  Top-N Selection + PDB Generation & Preparation
===============================================================================
Two phases on the MD-ready cohort only (MD_Selected column from Step 02, §18):
  Phase 1 - CIF → PDB conversion (Gemmi) + Schrödinger PrepWizard preparation.
  Phase 2 - Top-N extraction, Ramachandran validation, molecular handover files,
            and PyMOL/PLIP interaction figures.
Restricting both phases to the ~10 MD-ready complexes keeps this step cheap
instead of converting/preparing the entire predicted library.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 09 October 2026 <───────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 05_TopN_and_PDB_Preparation_DeFluorX.py
  Role          : "Builder + Selector" - gates on the MD-ready cohort (MD_Selected),
                  converts CIF outputs to analysis-ready PDB, prepares them with
                  PrepWizard, then extracts and renders that cohort for handover.
  Imports from  : 00_01_Project_Config_DeFluorX.py  (CFG - pH values, MD-selection §18)
                  00_02_Project_Utils_DeFluorX.py   (ConsoleColours, setup_logging,
                                            console helpers, Ramachandran helpers)
  Reads         : <Run>/2_Best_Complexes_CIFs/*.cif
                  <Run>/1_Boltz2_Production/6_Boltz2_DeFluorX_Ranked_*.csv  (MD_Selected)
                  <Run>/1_Boltz2_Production/1_Input_Data/*
  Writes        : <Run>/5_TopN_and_Preparation/  (one consolidated folder)
                    1_Converted_Raw_PDB/   (raw PDBs + Figures/)
                    2_Prepared_PDBs/       (prepared PDBs + Figures/)
                    3_Comparative_Analysis/ (Ramachandran, Controls, handover, combined CSV)
                    4_Ligand_ESP_Charges/  (<stem>_ESP.mae + 00_ESP_Charges_Summary.csv; gated on --esp)
                    00_TopN_and_Preparation.log  (single log for both phases)
  Upstream      : 02_Production_DeFluorX.py  → writes Best_Complexes_CIFs and ranked CSV
  Downstream    : 06_Physics_Validation_DeFluorX.py       → reads the MD-selected handover (R{N}_*.pdb
                                                          + *_ESP.mae) → WaterMap · System Builder · MD · SID · MM-GBSA
                  07_QMMM_Defluorination_DeFluorX.py   → reads the MD/WaterMap outputs for QM/MM defluorination
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
    python 05_TopN_and_PDB_Preparation_DeFluorX.py Boltz-2_Run_20260309T085406Z

── Key features ───────────────────────────────────────────────────────────────
  • Smart Validation: SOURCE_CIF tag in PDB header detects stale files from
    previous runs and triggers automatic regeneration.
  • Stale-file Detection: re-prepares if Raw PDB is newer than Prepared PDB.
  • Gemmi Conversion: robust CIF → PDB; moves non-protein residues to Chain L.
  • Smart Ligand Management: metals (Zn, Mg) and modified residues (MSE) stay
    with the protein chain to preserve topology for Maestro.
  • PrepWizard: fills side chains, PropKa pH 8.0, Epik pH 8.0, 0.15 Å RMSD
    restrained minimisation, disulfide detection.
  • MD-ready gate: both phases process only the MD_Selected cohort (CFG §18) plus
    the control jobs (ID 0000000_*), so the whole predicted library is never
    converted or prepared.
  • Top-N extraction: prepared-complex handover (FASTA/SMILES/SDF), Ramachandran
    validation, and PyMOL/PLIP/matplotlib interaction figures for each complex.
  • Interaction diagrams: PLIP (validated typing) plus a distance-based InteractionMap,
    both role-coloured from CFG (nucleophile / acid-base / clamp / fluoride pocket) with
    adaptive residue placement - upper-arc when few contacts, full ring when many.
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

# =============================================================================
# SECTION 1: SYSTEM CONFIGURATION & IMPORTS
# =============================================================================

# -----------------------------------------------------------------------------
# Step 1.1: Standard Library Imports
# -----------------------------------------------------------------------------
import os
import sys

'''
CPU usage cap (total cores − 2; mirrors CFG.PREP_CPU_RESERVE). Reserve 2 cores
for OS/desktop stability by limiting the thread-pool maths libraries (BLAS /
MKL / OpenMP / NumExpr). Must precede numpy/pandas import to take effect;
setdefault() preserves any value exported by the caller or pipeline runner.
'''
_CPU_CAP = str(max(1, (os.cpu_count() or 4) - 2))
for _tv in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_tv, _CPU_CAP)

import shutil
import csv
import warnings
import subprocess
import tempfile
import time
import argparse
import re
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- Additional stack for the Top-N extraction + figure phase ---
import xml.etree.ElementTree as _ET
from collections import defaultdict
import math
import numpy as np
import matplotlib
matplotlib.use("Agg")  # non-interactive backend; must be set before pyplot import
# --- consolidated matplotlib imports (after backend selection) ---
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import matplotlib.patches as _mpatches
from matplotlib.patches import FancyBboxPatch as _FancyBboxPatch, Rectangle as _Rectangle
from matplotlib.offsetbox import (HPacker as _HPacker, TextArea as _TextArea,
                                  DrawingArea as _DrawingArea, AnchoredOffsetbox as _AnchoredOffsetbox)
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.lines import Line2D as _Line2D
from Bio import SeqIO
from Bio.PDB import PDBParser as _PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')
DEFAULT_BASE_PATH = Path.cwd()


# -----------------------------------------------------------------------------
# Step 1.2: Scientific Stack
# -----------------------------------------------------------------------------
import gemmi

# -----------------------------------------------------------------------------
# Step 1.3: Pipeline modules (00_01 config, 00_02 utils) via importlib
# -----------------------------------------------------------------------------
'''
Filenames begin with digits and cannot be imported with standard `import`.
'''
import importlib.util as _ilu
# --- consolidated top-level imports (optional/heavy + Schrodinger stay function-local) ---
from collections import Counter as _Counter
import re as _re
import select as _select

def _load_module(name: str, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod  = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_REPO_DIR   = Path(__file__).resolve().parent
_cfg_mod    = _load_module("ProjectConfig", _REPO_DIR / "00_01_Project_Config_DeFluorX.py")
_utils_mod  = _load_module("ProjectUtils",  _REPO_DIR / "00_02_Project_Utils_DeFluorX.py")

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

# -----------------------------------------------------------------------------
# Step 1.4: Global Constants & Paths
# -----------------------------------------------------------------------------

'''
Ensure SCHRODINGER is set in the process environment so child processes
(PrepWizard subprocess calls) inherit it without requiring a prior `export`.
'''
def _latest_schrodinger() -> str:
    """Newest /opt/schrodinger* install that carries the `run` binary, so the pipeline follows a
    Schrodinger version upgrade (e.g. 2026-3 -> 2026-4) without an env edit. Falls back to the bare
    path only when nothing matches."""
    import glob as _glob, re as _re
    _cands = [d for d in _glob.glob("/opt/schrodinger*") if os.path.isdir(d)]
    _vkey = lambda p: [int(x) if x.isdigit() else x for x in _re.split(r"(\d+)", p)]
    for d in sorted(_cands, key=_vkey, reverse=True):
        if os.path.exists(os.path.join(d, "run")):
            return d
    return "/opt/schrodinger"

if "SCHRODINGER" not in os.environ:
    os.environ["SCHRODINGER"] = _latest_schrodinger()

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


# =============================================================================
# SECTION 2: LOGGING INFRASTRUCTURE
# =============================================================================

'''
Logging and console functions are provided by 00_02_Project_Utils.
Script-level wrappers capture the module-global `logger` so existing call
sites require no modification.
'''

def setup_logging(prep_base_dir: Path) -> Path:
    """Initialises the preparation log via the shared utility."""
    global logger
    prep_base_dir.mkdir(parents=True, exist_ok=True)
    primary_log = prep_base_dir / "00_TopN_and_Preparation.log"
    logger = _utils_mod.setup_logging(primary_log, logger_name="pdb_prep", mode="w", timestamp=True)
    return primary_log

def console_title(msg: str) -> None:
    _utils_mod.console_title(msg, logger)

def console_info(msg: str) -> None:
    _utils_mod.console_info(msg, logger)

def console_separator() -> None:
    _utils_mod.console_separator(logger, heavy=True)


# =============================================================================
# SECTION 3: CORE LOGIC & HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# Step 3.1: File Indexing & Metadata
# -----------------------------------------------------------------------------
def index_existing_files(directory: Path, suffix: str) -> dict:
    """Creates a fast lookup dictionary {JobIndex: Path} for existing files.
    Special handling for Controls (ID 0000000) to avoid collisions.
    """
    index = {}
    if not directory.exists(): return index

    for f in sorted(directory.glob(f"*{suffix}")):
        try:
            parts = f.name.split("_")
            if parts and parts[0].isdigit():
                idx = parts[0]
                '''
                If it's a control ID, it is not indexed for renaming
                because multiple controls share this ID.
                '''
                if idx == CFG.CONTROL_JOB_PREFIX:
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
    for cif_path in sorted(best_cifs_dir.glob("*.cif")):
        job_name = re.sub(r"_model_\d+$", "", cif_path.stem)
        mtime = cif_path.stat().st_mtime
        if job_name not in latest_cifs or mtime > latest_cifs[job_name][1]:
            latest_cifs[job_name] = (cif_path, mtime)

    # Return sorted list for deterministic order
    return [(jn, path) for jn, (path, _) in sorted(latest_cifs.items())]

def load_rank_map(prod_dir: Path) -> dict:
    """Loads the primary CSV to annotate PDB headers with a score/rank.
    Prefers the DeFluorX Ranked CSV (6_Boltz2_DeFluorX_Ranked_*), then any ranked CSV,
    then falls back to any primary CSV.
    Uses Scientific_Rank if present, otherwise Boltz_Model_Confidence (rounded to 4dp).
    """
    rank_csvs = sorted(prod_dir.glob(CFG.GLOB_RANKED_CSV))
    if not rank_csvs:
        rank_csvs = sorted(prod_dir.glob("*_Ranked_*.csv"))
    if not rank_csvs:
        console_info("Warning: No Ranked CSV found in production directory.")
        return {}

    latest_csv = _utils_mod.latest_by_mtime(rank_csvs)   # newest by mtime, not name (a name sort can rank an older file last when the leading number differs)
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
    rank_csvs = (sorted(prod_dir.glob(CFG.GLOB_RANKED_CSV)) or
                 sorted(prod_dir.glob("*_Ranked_*.csv")))
    if not rank_csvs:
        return set()
    try:
        df = pd.read_csv(_utils_mod.latest_by_mtime(rank_csvs), low_memory=False)
        if md_col not in df.columns or "job_name" not in df.columns:
            return set()
        _sel = df[md_col].astype(str).str.lower().isin(["true", "1", "1.0"])
        return set(df.loc[_sel, "job_name"].astype(str))
    except Exception as e:
        if logger: logger.debug(f"Error loading MD selection: {e}")
        return set()

# -----------------------------------------------------------------------------
# Step 3.2: PDB Manipulation (Header & Source Tracking)
# -----------------------------------------------------------------------------
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
            f"TITLE     Scientific Rank: {rank} | DeFluorX Pipeline",
            f"REMARK 999 SOURCE_CIF: {source_cif}"
        ]

        new_content = "\n".join(header_block + lines) + "\n"
        pdb_path.write_text(new_content)
        return True
    except Exception as e:
        if logger: logger.warning(f"Could not update header for {job_name}: {e}")
        return False

# -----------------------------------------------------------------------------
# Step 3.3: Conversion Engine (Gemmi)
# -----------------------------------------------------------------------------
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

            '''
            Re-number Chain L residues sequentially so multiple ligand copies
            (homodimers / multiple identical PFAS) never collide on res.seqid,
            which would break PyMOL/PLIP selections downstream.
            '''
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

# -----------------------------------------------------------------------------
# Step 3.4: Preparation Engine (Schrödinger PrepWizard)
# -----------------------------------------------------------------------------
def run_prepwizard(raw_pdb: Path, final_dest: Path):
    """
    Wraps Schrödinger's prepwizard for headless protein-ligand preparation.
    - Adds hydrogens (default behaviour; -nohtreat NOT used)
    - Fills missing side chains with Prime
    - PropKa protonation at pH 8.0 (FAcD physiological context)
    - Epik protonation of PFAS ligands at pH 8.0
    - Restrained minimisation (0.15 Å RMSD) to resolve clashes
    - Runs inline via -NOJOBID (no Schrödinger job-server dependency)
    - Uses /tmp for work dir - nothing written to the output folder
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

        '''
        Headless execution environment - strip Python-env overrides that confuse
        Schrödinger's bundled Python interpreter.
        '''
        env = os.environ.copy()
        env["SCHRODINGER"]     = str(SCHRODINGER_PATH)
        env["QT_QPA_PLATFORM"] = "offscreen"
        for k in ("PYTHONPATH", "PYTHONHOME", "LD_LIBRARY_PATH"):
            env.pop(k, None)

        cmd = [
            str(PREPWIZARD_BIN),
            "-fillsidechains",       # Fill truncated side chains (common in AI models)
            "-fillloops",            # Rebuild missing loops - experimental control PDBs
                                     # can have unresolved loops; a physical chain break
                                     # otherwise crashes the downstream Amber/Desmond build
            "-propka_pH", str(CFG.PREPWIZARD_PROPKA_PH),  # protein protonation pH
            "-epik_pH",   str(CFG.PREPWIZARD_EPIK_PH),    # PFAS ligand protonation via Epik
            "-r",         str(CFG.PREPWIZARD_RMSD_RESTRAIN),  # RMSD-restrained minimisation
            # Minimisation force field, stated rather than defaulted: PrepWizard's own default is
            # OPLS_2005, while every downstream stage (Desmond build, MD, WaterMap, Prime MM-GBSA)
            # runs the OPLS4 family. Minimising under one force field and simulating under another
            # means the prepared geometry is not a minimum of the potential the MD uses.
            "-f",         str(CFG.PREPWIZARD_FORCEFIELD),
            "-disulfides",           # Detect and bond proximal Cys pairs
            "-NOJOBID",              # Run inline - no Schrödinger Job Control layer
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


# =============================================================================
# SECTION 4: TASK EXECUTION & VALIDATION LOGIC
# =============================================================================

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
                    raw_pdb_path.unlink() # source tag mismatch after rename - drop the stale file
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

def _is_reference_control(job) -> bool:
    """The single control retained in the comparative figures - the CFG-designated MD control:
    the reference structure (CFG.REFERENCE_PDB_ID) paired with a CFG.CONTROL_MD_LIGANDS ligand
    (3R3U × fluoroacetate). The other reference systems (3R3U-DFA/TFA, DeHa4-*) are scored/tiered
    but not simulated, so they are dropped here. Nothing about the control is hardcoded - both the
    protein and the ligand come from CFG."""
    j = str(job)
    return CFG.REFERENCE_PDB_ID in j and any(str(lig) in j for lig in CFG.CONTROL_MD_LIGANDS)


def _lig_short_token(s: str) -> str:
    """Collapse the long PFAS acid name inside a compound job token to its CFG short form.
    The map is CFG.VIS_LIGAND_SHORT (single source of truth); the longest acid name is matched
    first so 'trifluoroacetate' is not shadowed by the 'fluoroacetate' substring."""
    for _full in ("trifluoroacetate", "difluoroacetate", "fluoroacetate"):
        s = s.replace(_full.title(), CFG.VIS_LIGAND_SHORT[_full])
    return s


def plot_pose_drift(geom_rows: list, out_dir: Path) -> Path | None:
    """What preparation does to the two numbers the tier is decided on.

    One row per complex, drawn as the pose\'s journey through three stages: the Boltz CIF the SCREEN
    scored (open circle), the RAW PDB after gemmi conversion (small dot), and the MINIMISED pose MD
    actually STARTS from (arrowhead). The gap between the circle and the arrowhead is the whole point.

    CIF and RAW are drawn as SEPARATE markers even though the conversion is lossless (it moves the angle
    by at most 0.04°): they land on top of each other, and that coincidence IS the result - it shows,
    per structure rather than as an aggregate footnote, that none of the drift is in the conversion and
    all of it is in the minimisation.

    Arrows are coloured by DIRECTION, because the direction is the finding: the poses the screen ranked
    highest move AWAY from the gate, the poses it ranked lowest move TOWARDS it. That is regression to
    the mean - the best-of-5 pick is partly luck, and minimisation takes the luck back. A figure that
    coloured by tier would hide it; a figure that coloured by direction cannot.
    """
    if not geom_rows:
        return None
    _utils_mod.apply_figure_style(CFG)

    df = pd.DataFrame(geom_rows)
    for _c in ("cif_sn2_angle", "cif_dist_nuc", "raw_sn2_angle", "prep_sn2_angle",
               "raw_dist_nuc", "prep_dist_nuc"):
        if _c in df.columns:
            df[_c] = pd.to_numeric(df[_c], errors="coerce")
    # CIF falls back to RAW for older rows written before the CIF stage was measured (lossless, so equal).
    for _cc, _rc in (("cif_sn2_angle", "raw_sn2_angle"), ("cif_dist_nuc", "raw_dist_nuc")):
        if _cc not in df.columns:
            df[_cc] = df[_rc]
        else:
            df[_cc] = df[_cc].fillna(df[_rc])
    df = df.dropna(subset=["raw_sn2_angle", "prep_sn2_angle"])
    if df.empty:
        return None
    '''
    Group the two cohorts on the y-axis: controls together at the foot, MD-selected together above
    them, each block sorted by prepared angle. Reading the MD candidates as one contiguous set is
    clearer than interleaving them with the controls by angle. The 'job' secondary key keeps equal
    angles deterministic (rows otherwise arrive in thread-completion order).
    '''
    df["_is_md"] = ~df["job"].astype(str).str.startswith(CFG.CONTROL_JOB_PREFIX)

    def _short(j):
        return (_lig_short_token("_".join(str(j).split("_")[2:]).replace("_Control", ""))
                .replace("_26", "").replace("_27", "").replace("_25", ""))

    # Keep the MD-selected candidates + only the CFG-designated control (3R3U × fluoroacetate); the
    # other reference systems (3R3U-DFA/TFA, DeHa4-*) are dropped so one canonical control is shown.
    _keep = df["_is_md"] | df["job"].map(_is_reference_control)
    df = df[_keep].reset_index(drop=True)
    if df.empty:
        return None
    df = df.sort_values(["_is_md", "prep_sn2_angle", "job"],
                        ascending=[True, True, True], kind="mergesort").reset_index(drop=True)

    _label = [_short(j) for j in df["job"]]
    _is_md = [not str(j).startswith(CFG.CONTROL_JOB_PREFIX) for j in df["job"]]
    _y = np.arange(len(df))

    fig, (ax_a, ax_d) = plt.subplots(1, 2, figsize=(17.0, 0.62 * len(df) + 3.0), sharey=True,
                                     gridspec_kw={"wspace": 0.06})

    def _panel(ax, ccif, c0, c1, gate, gate_lbl, relaxed, relaxed_lbl, xlab, worse_is, strict, dfmt="{:+.1f}"):
        """
        The two zones are SHADED, not merely ruled. A dashed line tells the reader where the gate is; a
        filled band tells them which side of it means competent, without translating a number first. For
        the angle the gate is a FLOOR (pass = to the right of it); for the distance a CEILING (pass = to
        the left). The fills carry that reversal so the eye does not have to. A dash-dot strict-NAC line
        (CFG.NAC_*_STRICT) is drawn as a third reference (keyed in the legend, not the margin).
        """
        _lo = min(df[ccif].min(), df[c0].min(), df[c1].min(), gate, relaxed, strict)
        _hi = max(df[ccif].max(), df[c0].max(), df[c1].max(), gate, relaxed, strict)
        _m = (_hi - _lo) * 0.13
        _x0, _x1 = _lo - _m, _hi + _m
        if worse_is == "down":                      # angle: elite is HIGH
            ax.axvspan(gate, _x1, color=CFG.VIS_TINT["green"], alpha=0.55, lw=0, zorder=0)
            ax.axvspan(_x0, relaxed, color=CFG.VIS_TINT["red"], alpha=0.55, lw=0, zorder=0)
            _zones = [((gate + _x1) / 2, "clears Tier_1A gate", CFG.VIS_BAND["high"]),
                      ((relaxed + gate) / 2, "relaxed NAC,\nbelow Tier_1A gate", CFG.VIS_INK["muted"]),
                      ((_x0 + relaxed) / 2, "outside NAC", CFG.VIS_BAND["low"])]
        else:                                       # distance: elite is LOW
            ax.axvspan(_x0, gate, color=CFG.VIS_TINT["green"], alpha=0.55, lw=0, zorder=0)
            ax.axvspan(relaxed, _x1, color=CFG.VIS_TINT["red"], alpha=0.55, lw=0, zorder=0)
            _zones = [((_x0 + gate) / 2, "clears Tier_1A gate", CFG.VIS_BAND["high"]),
                      ((gate + relaxed) / 2, "relaxed NAC,\nbelow Tier_1A gate", CFG.VIS_INK["muted"]),
                      ((relaxed + _x1) / 2, "outside NAC", CFG.VIS_BAND["low"])]
        # The three zones are labelled WHERE THEY ARE, at the foot of the panel, so the reader never has
        # to translate a legend swatch into a region. Each label is centred in its own band.
        for _zx, _zt, _zc in _zones:
            ax.annotate(_zt, xy=(_zx, 0.0), xycoords=("data", "axes fraction"),
                        xytext=(0, 4), textcoords="offset points", ha="center", va="bottom",
                        fontsize=CFG.VIS_FONT_ANNOT - 1.0, color=_zc, style="italic",
                        alpha=0.9, zorder=1, clip_on=True)
        for _i, (_cf, _a, _b) in enumerate(zip(df[ccif], df[c0], df[c1])):
            '''
            Three stages on one row: CIF (open circle) → RAW (small dot, ≈CIF) → MINIMISED (arrowhead).
            The arrow runs RAW → MINIMISED, because that leg carries all the movement; the CIF→RAW leg
            is the lossless conversion and is shown only as the two coincident markers.
            'worse' = away from the gate. For the angle the gate is a floor; for the distance a ceiling.
            '''
            _worse = (_b < _a) if worse_is == "down" else (_b > _a)
            _col = CFG.VIS_BAND["low"] if _worse else CFG.VIS_BAND["high"]
            ax.annotate("", xy=(_b, _i), xytext=(_a, _i),
                        arrowprops=dict(arrowstyle="-|>,head_width=0.28,head_length=0.6",
                                        color=_col, lw=2.2, shrinkA=0, shrinkB=0), zorder=4)
            # RAW: small filled dot at the conversion pose.
            ax.plot([_a], [_i], "o", ms=4.5, color=CFG.VIS_INK["muted"], zorder=6)
            # CIF: open circle - the pose the screen scored. Coincident with RAW when lossless.
            ax.plot([_cf], [_i], "o", ms=8, color=CFG.VIS_INK["white"],
                    markeredgecolor=CFG.VIS_INK["dark"], markeredgewidth=1.3, zorder=5)
            # The value rides BEYOND the arrowhead, in the direction of travel, so it can never sit on
            # top of the gate line the arrow is crossing.
            _dir = 1 if _b >= _a else -1
            # Endpoint value + the signed minimisation drift (RAW→minimised), e.g. "177.0 (+3)".
            ax.annotate(f"{_b:.1f} ({dfmt.format(_b - _a)})", xy=(_b, _i), xytext=(9 * _dir, 0),
                        textcoords="offset points", ha="left" if _dir > 0 else "right", va="center",
                        fontsize=CFG.VIS_FONT_ANNOT - 0.5, color=_col, zorder=6)
        ax.axvline(gate, ls="--", lw=1.4, color=CFG.VIS_INK["dark"], alpha=0.9, zorder=2)
        ax.axvline(relaxed, ls=":", lw=1.2, color=CFG.VIS_INK["ghost"], zorder=2)
        ax.axvline(strict, ls="-.", lw=1.3, color=CFG.VIS_ACCENT["blue"], alpha=0.9, zorder=2)
        '''
        The gate labels sit in the MARGIN ABOVE the panel, not inside it. Rotated 90° across the data
        they crossed arrows, values and rows - the reader had to decode the label before reading the
        plot. Anchored to the axis in data-x and figure-y, they stay attached to their own line and
        collide with nothing.
        '''
        for _gv, _gl, _gc in ((gate, gate_lbl, CFG.VIS_INK["dark"]),
                              (relaxed, relaxed_lbl, CFG.VIS_INK["ghost"])):
            ax.annotate(_gl, xy=(_gv, 1.0), xycoords=("data", "axes fraction"),
                        xytext=(0, 6), textcoords="offset points",
                        ha="center", va="bottom", rotation=0,
                        fontsize=CFG.VIS_FONT_ANNOT, color=_gc, clip_on=False)
        ax.set_xlabel(xlab)
        ax.grid(True, axis="x", alpha=CFG.VIS_GRID_ALPHA, color=CFG.VIS_GRID_COLOUR)
        ax.set_axisbelow(True)
        ax.set_xlim(_x0, _x1)

    _panel(ax_a, "cif_sn2_angle", "raw_sn2_angle", "prep_sn2_angle",
           CFG.TIER_ANGLE_MIN[CFG.TIER_TOP], f"Tier_1A gate  {CFG.TIER_ANGLE_MIN[CFG.TIER_TOP]:.0f}°",
           CFG.NAC_ANGLE_RELAXED, f"relaxed NAC  {CFG.NAC_ANGLE_RELAXED:.0f}°",
           "SN2 attack angle at the mapped nucleophile  (°)", worse_is="down",
           strict=CFG.NAC_ANGLE_STRICT, dfmt="{:+.0f}")
    _panel(ax_d, "cif_dist_nuc", "raw_dist_nuc", "prep_dist_nuc",
           CFG.TIER_NUC_DIST[CFG.TIER_TOP], f"Tier_1A gate  {CFG.TIER_NUC_DIST[CFG.TIER_TOP]:.1f} Å",
           CFG.NAC_DIST_RELAXED, f"relaxed NAC  {CFG.NAC_DIST_RELAXED:.1f} Å",
           "Nucleophile distance  (Å)", worse_is="up",
           strict=CFG.NAC_DIST_STRICT, dfmt="{:+.1f}")

    ax_a.set_yticks(_y)
    ax_a.set_yticklabels(_label)
    for _t, _md in zip(ax_a.get_yticklabels(), _is_md):
        _t.set_color(CFG.VIS_ACCENT["blue"] if _md else CFG.VIS_INK["muted"])
    ax_a.set_ylim(-0.8, len(df) - 0.2)

    '''
    The movement handles carry an arrowhead, because the figure draws arrows - the legend must show
    the same mark. The three zone bands are NOT in the legend: they are labelled in place, inside their
    own fills at the foot of each panel, so the legend only explains the marks the reader cannot
    otherwise decode (the two pose markers and the two drift directions).
    '''
    _h = [_Line2D([0], [0], marker="o", ls="", markerfacecolor=CFG.VIS_INK["white"],
                  markeredgecolor=CFG.VIS_INK["dark"], ms=8, label="CIF  (Boltz - what the screen scored)"),
          _Line2D([0], [0], marker="o", ls="", markerfacecolor=CFG.VIS_INK["muted"],
                  markeredgecolor="none", ms=5, label="RAW  (converted - lossless, sits on CIF)"),
          _Line2D([0], [0], color=CFG.VIS_BAND["low"], lw=2.2, marker=">", markersize=9,
                  markerfacecolor=CFG.VIS_BAND["low"], markeredgecolor=CFG.VIS_BAND["low"],
                  label="minimised pose moved AWAY from the gate"),
          _Line2D([0], [0], color=CFG.VIS_BAND["high"], lw=2.2, marker=">", markersize=9,
                  markerfacecolor=CFG.VIS_BAND["high"], markeredgecolor=CFG.VIS_BAND["high"],
                  label="minimised pose moved TOWARDS the gate"),
          _Line2D([0], [0], color=CFG.VIS_ACCENT["blue"], ls="-.", lw=1.4,
                  label=f"strict NAC  ({CFG.NAC_ANGLE_STRICT:.0f}° / {CFG.NAC_DIST_STRICT:.1f} Å)")]
    ax_a.legend(handles=_h, loc="lower center", bbox_to_anchor=(1.03, 1.045), ncol=len(_h),
                frameon=False, fontsize=CFG.VIS_FONT_LEGEND)

    _da = pd.to_numeric(df["prep_d_angle"], errors="coerce").dropna()
    _dd = pd.to_numeric(df["prep_d_dist"], errors="coerce").dropna()
    fig.text(0.5, 0.012,
             f"CIF (Boltz) → RAW (gemmi, lossless) → Minimised (PrepWizard, MD start).   "
             f"Minimisation drift: mean |Δangle| {_da.abs().mean():.1f}°, mean Δdist {_dd.mean():+.2f} Å.   "
             f"Blue = MD-selected.",
             ha="center", fontsize=CFG.VIS_FONT_ANNOT, color=CFG.VIS_INK["muted"])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(rect=(0, 0.045, 1, 0.90))

    out_dir.mkdir(parents=True, exist_ok=True)
    _path = out_dir / "02_Pose_Drift_CIF_to_Prepared.svg"
    fig.savefig(_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return _path


'''
The eight catalytic residues, in the order the figure stacks them: the reactive relay that does the
chemistry first, then the machinery that positions the substrate. Each entry is
(role_key, ranked-CSV column, short label, functional group, partner). The PARTNER is what the residue
actually engages mechanistically - and it is NOT the substrate for every residue. FAcD's catalysis is a
RELAY: the nucleophile attacks the substrate carbon, but the base (His) engages the NUCLEOPHILE, and the
acid (Asp) engages the BASE (the Asp–His dyad that polarises it, "not a ligand contact" per Step 02).
Measuring the acid or the base against the substrate would report them as "disengaged" when the dyad is
in fact intact, which is exactly backwards. So each residue is measured against its own partner:
"substrate" for the direct contacts, or another role_key for the relay links.
'''
_MACHINERY_ROLES = [
    ("Nuc",     "Mapped_Nucleophile",  "Nucleophile → warhead C", "reactive",   "substrate"),
    ("Base",    "Mapped_Base",         "Base → nucleophile",      "reactive",   "Nuc"),
    ("Acid",    "Mapped_Acid",         "Acid → base (dyad)",      "reactive",   "Base"),
    ("StabH",   "Mapped_Stabiliser_H", "F⁻ stabiliser → substrate", "stabiliser", "substrate"),
    ("Clamp1",  "Mapped_Clamp1",       "Carboxylate clamp (Arg) → substrate", "clamp", "substrate"),
    ("Clamp2",  "Mapped_Clamp2",       "Carboxylate clamp (Arg) → substrate", "clamp", "substrate"),
    ("CradleW", "Mapped_Stabiliser_W", "Aromatic cradle (Trp) → substrate", "cradle", "substrate"),
    ("CradleY", "Mapped_Stabiliser_Y", "Aromatic cradle (Tyr) → substrate", "cradle", "substrate"),
]


def load_machinery_map(prod_dir: Path) -> dict:
    """Per-job residue numbers for ALL EIGHT catalytic roles, from 02's dynamic alignment
    (the Mapped_* columns of the ranked CSV). Returns {job_name: {role_key: int|None}}.

    This is the 8-residue superset of load_catalytic_anchor_map, which returns only the three
    protonation-relevant roles. The extra five (the fluoride stabiliser, the two carboxylate
    clamps, the two aromatic-cradle residues) are read only for the machinery-engagement figure,
    never for protonation, so they live in their own loader rather than widening the anchor map."""
    _rank = sorted(prod_dir.glob(CFG.GLOB_RANKED_CSV)) or sorted(prod_dir.glob("*_Ranked_*.csv"))
    if not _rank:
        return {}
    try:
        _df = pd.read_csv(_utils_mod.latest_by_mtime(_rank), low_memory=False)
    except Exception:
        return {}
    if "job_name" not in _df.columns:
        return {}
    def _num(v):
        _m = _re.search(r"(\d+)\s*$", str(v))
        return int(_m.group(1)) if _m else None
    _out = {}
    for _, _row in _df.iterrows():
        _out[_row["job_name"]] = {rk: (_num(_row[col]) if col in _df.columns else None)
                                  for rk, col, _lbl, _grp, _partner in _MACHINERY_ROLES}
    return _out


def measure_machinery_engagement(struct_path: Path, role_resnums: dict) -> dict:
    """How close each catalytic residue sits to the substrate - the minimum heavy-atom distance from
    the residue's SIDECHAIN to the nearest ligand heavy atom, for every mapped role.

    'Engagement' is deliberately the sidechain-to-ligand contact distance, not the Cα distance: the
    chemistry is done by the sidechains (the Asp carboxylate attacks, the His imidazole relays a proton,
    the Arg guanidinium clamps the substrate carboxylate, the Trp/Tyr ring cradles it). A residue whose
    backbone is in place but whose sidechain has swung away is DISENGAGED, and only the sidechain metric
    catches that. Returns {role_key: min_distance_A}; a role with no mapped residue or no sidechain atoms
    present is NaN so the figure can mark it 'not resolved' rather than draw a false contact."""
    _AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO",
           "SER","THR","TRP","TYR","VAL","HID","HIE","HIP","ASH","GLH","LYN","HOH","WAT","NA","CL","SPC","T3P"}
    _bb = {"N", "CA", "C", "O", "OXT", "H", "HA"}
    _out = {rk: float("nan") for rk, *_ in _MACHINERY_ROLES}
    try:
        st = gemmi.read_structure(str(struct_path)); st.setup_entities(); st.remove_hydrogens()
    except Exception:
        return _out
    _lig, _prot = [], []
    for _ch in st[0]:
        for _r in _ch:
            (_prot if _r.name.strip().upper() in _AA else _lig).append(_r)
    _lig = [r for r in _lig if len(r) > 2]
    if not _lig:
        return _out
    _L = max(_lig, key=len)
    _lig_pos = [a.pos for a in _L if a.element.name != "H"]
    if not _lig_pos:
        return _out
    _by_num = {}
    for _r in _prot:
        _by_num.setdefault(_r.seqid.num, _r)

    def _sidechain(_num):
        if _num is None:
            return None
        _res = _by_num.get(int(_num))
        if _res is None:
            return None
        _side = [a.pos for a in _res if a.element.name != "H" and a.name.strip().upper() not in _bb]
        if not _side:                                    # glycine-like: no sidechain, use all heavy atoms
            _side = [a.pos for a in _res if a.element.name != "H"]
        return _side or None

    _partner_of = {rk: partner for rk, _c, _l, _g, partner in _MACHINERY_ROLES}
    for _rk, _num in (role_resnums or {}).items():
        _side = _sidechain(_num)
        if not _side:
            continue
        _partner = _partner_of.get(_rk, "substrate")
        if _partner == "substrate":
            _target = _lig_pos                            # residue → nearest substrate heavy atom
        else:
            _target = _sidechain((role_resnums or {}).get(_partner))   # residue → partner residue sidechain
        if not _target:
            continue
        _out[_rk] = float(min(_sp.dist(_tp) for _sp in _side for _tp in _target))
    return _out


def plot_machinery_distribution(rows: list, out_dir: Path) -> Path | None:
    """All eight catalytic residues in one mixed figure: for each residue, the cohort's engagement to
    its mechanistic partner is drawn as a raincloud - a violin (the distribution across complexes at the
    minimised, MD-start pose), a box (median and quartiles) inside it, and a strip of the individual
    complexes on top (each MD-selected complex in its own colour, controls grey). A thin stem runs from
    the CIF median down/up to the minimised median, so the drift of the whole cohort through preparation
    is visible per residue. The three shaded bands read reactive contact → in contact → out of contact.

    This is the cohort-level companion to the per-complex figures: it answers whether, across every
    candidate, each piece of the machinery still sits where the chemistry needs it after minimisation -
    the reactive relay (nucleophile → warhead C, base → nucleophile, acid → base dyad) drawn first and
    heavier, then the substrate-positioning residues (fluoride stabiliser, carboxylate clamps, aromatic
    cradle)."""
    if not rows:
        return None
    # Keep the MD-selected candidates + only the CFG-designated control (3R3U × fluoroacetate); drop the
    # other reference systems (3R3U-DFA/TFA, DeHa4-*) from every layer (violin, strip, drift stem, legend).
    rows = [r for r in rows if not str(r["job"]).startswith(CFG.CONTROL_JOB_PREFIX)
            or _is_reference_control(r["job"])]
    if not rows:
        return None
    _utils_mod.apply_figure_style(CFG)

    _react = float(CFG.NAC_DIST_STRICT)
    _outer = float(getattr(CFG, "THRESHOLD_SALT_BRIDGE", 5.0))
    _grp_col = {"reactive": CFG.VIS_BAND["high"], "stabiliser": CFG.VIS_BAND["moderate"],
                "clamp": CFG.VIS_ACCENT["blue"], "cradle": CFG.VIS_INK["mid"]}
    _pal = CFG.VIS_TREND_SERIES
    _md_jobs = sorted(r["job"] for r in rows if not str(r["job"]).startswith(CFG.CONTROL_JOB_PREFIX))
    _md_colour = {j: _pal[_i % len(_pal)] for _i, j in enumerate(_md_jobs)}

    _n = len(_MACHINERY_ROLES)
    fig, ax = plt.subplots(figsize=(1.55 * _n + 2.0, 8.2))
    _allvals = [v for r in rows for _d in (r.get("cif_eng", {}), r.get("prep_eng", {}))
                for v in _d.values() if v is not None and v == v]
    _ymax = min(max(_allvals) if _allvals else _outer, _outer + 3.0) + 0.4
    _ymin = max(0.0, (min(_allvals) if _allvals else 0.0) - 0.3)
    ax.axhspan(_ymin, _react, color=CFG.VIS_BAND["high_fill"], lw=0, zorder=0)
    ax.axhspan(_react, _outer, color=CFG.VIS_BAND["mod_fill"], lw=0, zorder=0)
    ax.axhspan(_outer, _ymax, color=CFG.VIS_BAND["low_fill"], lw=0, zorder=0)
    ax.axhline(_react, ls="--", lw=1.1, color=CFG.VIS_INK["dark"], alpha=0.75, zorder=1)
    ax.axhline(_outer, ls=":", lw=1.0, color=CFG.VIS_INK["ghost"], zorder=1)

    _rng = np.random.default_rng(0)
    for _i, (_rk, _col_name, _lbl, _grp, _partner) in enumerate(_MACHINERY_ROLES):
        _gc = _grp_col.get(_grp, CFG.VIS_INK["dark"])
        _prep = [r.get("prep_eng", {}).get(_rk) for r in rows]
        _prep = [v for v in _prep if v is not None and v == v]
        _cif = [r.get("cif_eng", {}).get(_rk) for r in rows]
        _cif = [v for v in _cif if v is not None and v == v]
        if len(_prep) >= 3:                      # violin needs a few points to form a shape
            _vp = ax.violinplot([_prep], positions=[_i], widths=0.8, showextrema=False)
            for _b in _vp["bodies"]:
                _b.set_facecolor(_gc); _b.set_alpha(0.20); _b.set_edgecolor(_gc); _b.set_linewidth(0.9); _b.set_zorder(2)
        # strip: each complex as a point, MD coloured, controls grey; small horizontal jitter
        for _r in rows:
            _v = _r.get("prep_eng", {}).get(_rk)
            if _v != _v or _v is None:
                continue
            _md = not str(_r["job"]).startswith(CFG.CONTROL_JOB_PREFIX)
            _pc = _md_colour[_r["job"]] if _md else CFG.VIS_INK["silver"]
            ax.plot(_i + float(_rng.uniform(-0.14, 0.14)), _v, "o", ms=6.5 if _md else 4.5,
                    color=_pc, markeredgecolor=CFG.VIS_INK["dark"], markeredgewidth=0.6,
                    alpha=0.95 if _md else 0.7, zorder=6 if _md else 5)
        '''
        Drift stem: the MD-SELECTED cohort's median engagement, CIF → minimised. Only the MD-selected
        complexes are pooled here - those are the poses that go to MD, and pooling the controls (a
        different enzyme in the DeHa4 case) into one median would not be a meaningful number.
        '''
        _cif_md = [r.get("cif_eng", {}).get(_rk) for r in rows
                   if not str(r["job"]).startswith(CFG.CONTROL_JOB_PREFIX)]
        _cif_md = [v for v in _cif_md if v is not None and v == v]
        _prep_md = [r.get("prep_eng", {}).get(_rk) for r in rows
                    if not str(r["job"]).startswith(CFG.CONTROL_JOB_PREFIX)]
        _prep_md = [v for v in _prep_md if v is not None and v == v]
        if _cif_md and _prep_md:
            _mc, _mp = float(np.median(_cif_md)), float(np.median(_prep_md))
            ax.plot([_i - 0.36, _i - 0.36], [_mc, _mp], "-", color=_gc, lw=1.6, alpha=0.8, zorder=4)
            ax.plot([_i - 0.36], [_mc], "o", ms=6, markerfacecolor=CFG.VIS_INK["white"],
                    markeredgecolor=_gc, markeredgewidth=1.4, zorder=5)
            ax.plot([_i - 0.36], [_mp], "o", ms=6, color=_gc, markeredgecolor=CFG.VIS_INK["dark"],
                    markeredgewidth=0.5, zorder=5)

    # separate the reactive relay (first 3) from the positioning machinery
    _n_react = sum(1 for _rk, _c, _l, _g, _p in _MACHINERY_ROLES if _g == "reactive")
    # dotted separators between every residue column (same style as the per-tier separators in 03)
    for _sx in range(_n - 1):
        if abs((_sx + 0.5) - (_n_react - 0.5)) > 1e-6:      # the relay|positioning divider is drawn solid next
            ax.axvline(_sx + 0.5, color=CFG.VIS_INK["mid"], ls=":", alpha=0.5, lw=1.0, zorder=1)
    ax.axvline(_n_react - 0.5, ls="-", lw=0.8, color=CFG.VIS_INK["pale"], zorder=1)
    ax.annotate("reactive relay", xy=((_n_react - 1) / 2, 1.0), xycoords=("data", "axes fraction"),
                xytext=(0, 6), textcoords="offset points", ha="center", va="bottom",
                fontsize=CFG.VIS_FONT_ANNOT, color=CFG.VIS_BAND["high"], fontweight="bold")
    ax.annotate("substrate-positioning machinery", xy=((_n_react + _n - 1) / 2, 1.0),
                xycoords=("data", "axes fraction"), xytext=(0, 6), textcoords="offset points",
                ha="center", va="bottom", fontsize=CFG.VIS_FONT_ANNOT, color=CFG.VIS_INK["muted"])

    ax.set_xticks(range(_n))
    _xt = ax.set_xticklabels([_lbl for _rk, _c, _lbl, _g, _p in _MACHINERY_ROLES],
                             rotation=30, ha="right", fontsize=CFG.VIS_FONT_ANNOT)
    for _t, (_rk, _c, _lbl, _g, _p) in zip(_xt, _MACHINERY_ROLES):
        _t.set_color(_grp_col.get(_g, CFG.VIS_INK["dark"]))
        if _g == "reactive":
            _t.set_fontweight("bold")
    ax.set_xlim(-0.7, _n - 0.3)
    ax.set_ylim(_ymin, _ymax)
    ax.set_ylabel("Distance to mechanistic partner  (Å)", )
    ax.yaxis.set_major_locator(MultipleLocator(0.2))     # finer gridlines for reading the tight spread
    ax.grid(True, axis="y", alpha=CFG.VIS_GRID_ALPHA, color=CFG.VIS_GRID_COLOUR)
    ax.set_axisbelow(True)
    '''
    Zone labels at the right edge, each centred in its VISIBLE band. The reactive label rides just
    ABOVE its dashed line on the top layer (high zorder) so the line never sits over it. The "in
    contact" band is only drawn up to the visible top (_outer may exceed the axis), and the "out of
    contact" label is drawn only when that band is actually on screen (_ymax > _outer) - otherwise
    the two upper labels would collide at the top.
    '''
    _mod_top = min(_outer, _ymax)
    _zlabels = [(_react, f"reactive ≤{_react:g} Å", CFG.VIS_BAND["high"], "bottom", 3),
                ((_react + _mod_top) / 2, f"in contact ≤{_outer:g} Å", CFG.VIS_INK["muted"], "center", 0)]
    if _ymax > _outer + 0.05:
        _zlabels.append(((_outer + _ymax) / 2, "out of contact", CFG.VIS_BAND["low"], "center", 0))
    for _yv, _yt, _yc, _va, _dy in _zlabels:
        ax.annotate(_yt, xy=(1.0, _yv), xycoords=("axes fraction", "data"), xytext=(-4, _dy),
                    textcoords="offset points", ha="right", va=_va, style="italic",
                    fontsize=CFG.VIS_FONT_ANNOT - 1.0, color=_yc, alpha=0.95, zorder=7)

    _h = ([_Line2D([0], [0], marker="o", ls="", markerfacecolor=CFG.VIS_INK["white"],
                   markeredgecolor=CFG.VIS_INK["dark"], ms=6, label="MD-selected median @ CIF"),
           _Line2D([0], [0], marker="o", ls="", markerfacecolor=CFG.VIS_INK["dark"],
                   markeredgecolor="none", ms=6, label="MD-selected median @ minimised")]
          + [_Line2D([0], [0], marker="o", ls="", markerfacecolor=_md_colour[_j],
                     markeredgecolor=CFG.VIS_INK["dark"], ms=6, label=_short_job(_j)) for _j in _md_jobs]
          + [_Line2D([0], [0], marker="o", ls="", markerfacecolor=CFG.VIS_INK["silver"],
                     markeredgecolor=CFG.VIS_INK["dark"], ms=5, label=f"{CFG.REFERENCE_PDB_ID} control")])
    ax.legend(handles=_h, loc="upper left", frameon=False, fontsize=CFG.VIS_FONT_LEGEND - 0.5, ncol=2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout(rect=(0, 0.0, 1, 0.95))
    out_dir.mkdir(parents=True, exist_ok=True)
    _path = out_dir / "03_Machinery_Engagement_Distribution.svg"
    fig.savefig(_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return _path


def _short_job(j):
    return (_lig_short_token("_".join(str(j).split("_")[2:]).replace("_Control", ""))
            .replace("_26", "").replace("_27", "").replace("_25", ""))


def measure_sn2_geometry(struct_path: Path, nuc_resnum: int | None = None) -> dict:
    """Measure the SN2 geometry of a structure - CIF or PDB - with the pipeline's own joint-NAC rule.

    WHY THIS EXISTS
    ---------------
    The tier is decided on the Boltz CIF. The MD starts from the PREPARED PDB. Those are not the same
    structure, and the difference is neither small nor random: measured across the MD picks and the
    controls, PrepWizard's restrained minimisation moves the attack angle by 6.3 deg on average (up to
    18.5) and pushes the nucleophile 0.30 A further away, in 8 of 9 structures. Worse, the shift is
    DIRECTIONAL - the poses the screen selected degrade (-9.7, -7.3 deg) while the poses it rejected
    improve (+18.5, +8.7). That is regression to the mean: the best-of-5 pick is partly luck, and
    minimisation takes the luck back.

    So the tier ladder - 5 deg rungs, a 0.2 A distance step - discriminates more finely than the
    structure it screens is reproducible. Nothing here changes the ranking: the screen over 58,000
    complexes can only ever run on the CIF, because preparation costs minutes per structure. What this
    does is MEASURE the pose that is actually handed to MD, so the drift is visible instead of silent.

    THE NUCLEOPHILE IS THE NUCLEOPHILE
    ----------------------------------
    The attack is made by ONE residue - the catalytic aspartate that 02 mapped for this protein
    (Mapped_Nucleophile, e.g. ASP110), found by global alignment and robust to insertions. It is passed
    in as nuc_resnum. Scanning every Asp/Glu in the protein and taking whichever scores best would let a
    surface carboxylate 6 A away win a geometry contest it has no business entering, and report an
    attack that no enzyme makes. Without the mapped residue the measurement REFUSES to guess: it returns
    status "no_mapped_nucleophile" rather than a plausible number from the wrong atom.

    Between the two oxygens OF THAT residue the choice is the pipeline's own JOINT-NAC rule: the
    attacking oxygen must satisfy the distance AND the angle on the SAME atom, chosen by maximising
    sigmoid(d) x sigmoid(angle). Taking the nearest oxygen, or the best-angle oxygen, pairs a distance
    from one atom with an angle from the other and describes a nucleophile that does not exist.
    """
    _AA = {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO",
           "SER","THR","TRP","TYR","VAL","HID","HIE","HIP","ASH","GLH","LYN","HOH","WAT","NA","CL","SPC","T3P"}
    _out = {"sn2_angle": float("nan"), "dist_nuc": float("nan"), "attack_o": "", "status": "ok"}
    try:
        st = gemmi.read_structure(str(struct_path))
        st.setup_entities()
        st.remove_hydrogens()
    except Exception as _e:                                       # noqa: BLE001
        _out["status"] = f"unreadable:{type(_e).__name__}"
        return _out

    _lig, _prot = [], []
    for _ch in st[0]:
        for _r in _ch:
            (_prot if _r.name.strip().upper() in _AA else _lig).append(_r)
    _lig = [r for r in _lig if len(r) > 2]
    if not _lig:
        _out["status"] = "no_ligand"
        return _out
    _L = max(_lig, key=len)
    _F = [a for a in _L if a.element.name == "F"]
    _C = [a for a in _L if a.element.name == "C"]
    if not _F or not _C:
        _out["status"] = "no_C_F"
        return _out
    # the scissile C-F: a bonded C-F pair. The alpha carbon carries the leaving fluorine.
    _pairs = [(c, f) for c in _C for f in _F if c.pos.dist(f.pos) < CFG.CF_BOND_MAX_A]
    if not _pairs:
        _out["status"] = "no_CF_bond"
        return _out

    if nuc_resnum is None:
        _out["status"] = "no_mapped_nucleophile"
        return _out
    '''
    ASP/ASH only, and only the aspartate's own carboxylate oxygens. FAcD attacks with an aspartate;
    admitting GLU here would let a glutamate at the mapped position pass as the nucleophile, which is
    the same mistake - by residue type instead of by residue number - that put a serine in Step 03.
    '''
    _nucs = [(r, a) for r in _prot
             if r.seqid.num == int(nuc_resnum) and r.name.strip().upper() in ("ASP", "ASH")
             for a in r if a.name.strip() in ("OD1", "OD2")]
    if not _nucs:
        _out["status"] = f"nucleophile_{nuc_resnum}_not_found"
        return _out

    def _sig(x, k, x0):
        try:
            return 1.0 / (1.0 + math.exp(-k * (x - x0)))
        except OverflowError:
            return 0.0 if k * (x - x0) < 0 else 1.0

    _best = None
    for _c, _f in _pairs:
        for _r, _o in _nucs:
            _d = _o.pos.dist(_c.pos)
            if _d > CFG.NUC_SEARCH_RADIUS_A:
                continue   # the ligand has left the catalytic site entirely
            _v1 = np.array([_o.pos.x - _c.pos.x, _o.pos.y - _c.pos.y, _o.pos.z - _c.pos.z])
            _v2 = np.array([_f.pos.x - _c.pos.x, _f.pos.y - _c.pos.y, _f.pos.z - _c.pos.z])
            _cs = float(np.dot(_v1, _v2) / (np.linalg.norm(_v1) * np.linalg.norm(_v2) + 1e-12))
            _ang = math.degrees(math.acos(max(-1.0, min(1.0, _cs))))
            '''
            JOINT: the same oxygen must satisfy the distance AND the angle. The distance sigmoid is
            written EXACTLY as Step 02 writes it - plain d, with the negative SOFT_K_NUC doing the
            inversion. Negating d as well would invert it twice and score the FARTHEST oxygen best.
            '''
            _score = _sig(_d, CFG.SOFT_K_NUC, CFG.NAC_DIST_STRICT) * _sig(_ang, CFG.SOFT_K_ANG, CFG.NAC_ANGLE_STRICT)
            if _best is None or _score > _best[0]:
                _best = (_score, _ang, _d, f"{_r.name}{_r.seqid.num}:{_o.name.strip()}")
    if _best is None:
        _out["status"] = "no_nucleophile_in_range"
        return _out
    _out.update({"sn2_angle": round(_best[1], 2), "dist_nuc": round(_best[2], 2), "attack_o": _best[3]})
    return _out


def _build_hd1_line(res_lines: list) -> str | None:
    """Build the PDB line for a histidine Nδ1-H (HD1), needed to convert an HIE base to HID.

    HD1 lies in the imidazole plane along the external bisector of the CG–ND1–CE1 angle, one N–H bond
    length (1.01 Å) from ND1. The line is cloned from an existing hydrogen of the same residue (HE2,
    the proton being removed) so the column layout, chain, resnum and residue name match exactly; only
    the atom name and coordinates change. Returns None if ND1/CG/CE1 are not all present.
    """
    def _xyz(_l):
        return np.array([float(_l[30:38]), float(_l[38:46]), float(_l[46:54])], float)
    _by = {_l[12:16].strip().upper(): _l for _l in res_lines}
    if not {"ND1", "CG", "CE1"}.issubset(_by) or "HE2" not in _by:
        return None
    _nd1, _cg, _ce1 = _xyz(_by["ND1"]), _xyz(_by["CG"]), _xyz(_by["CE1"])
    _u = _cg - _nd1;  _u /= (np.linalg.norm(_u) or 1.0)
    _v = _ce1 - _nd1; _v /= (np.linalg.norm(_v) or 1.0)
    _d = -(_u + _v)
    _n = np.linalg.norm(_d)
    if _n < 1e-6:                        # CG/ND1/CE1 collinear - cannot place HD1
        return None
    _hd1 = _nd1 + (_d / _n) * 1.01
    _tmpl = _by["HE2"]
    _name = _tmpl[12:16].replace("HE2", "HD1")            # keep the template's justification
    return (_tmpl[:12] + _name + _tmpl[16:30]
            + f"{_hd1[0]:8.3f}{_hd1[1]:8.3f}{_hd1[2]:8.3f}" + _tmpl[54:])


def enforce_catalytic_protonation(pdb_path: Path, anchors: dict, job_name: str) -> dict:
    """Give each catalytic residue the protonation its ROLE requires (CFG §15).

    Every rule is read from CFG.CATALYTIC_PROTONATION_POLICY - one entry per catalytic role, each tied to
    ONE canonical residue via `role_key` (CFG.ROLE_EXPECTED_RESIDUES). There is no residue-type guess
    anywhere: no Asp-or-Glu, no Asp/Glu/Ser branch. Only the three roles PropKa cannot infer - the
    aspartate nucleophile, the dyad aspartate and the histidine base - are carried in `anchors`
    (load_catalytic_anchor_map builds it from the Nuc/Acid/Base Mapped_* columns) and are ENFORCED here.
    The other five CATALYTIC_PROTONATION_POLICY roles are enforce=False (never modified) and are NOT loaded
    into `anchors`, so this function neither re-verifies their identity nor logs them.

    The cost of getting the trio wrong is measured, not theoretical: a controlled experiment (two systems
    identical but for one hydrogen) collapses the SN2 attack angle 150° → 98° when the base is HIP (+1),
    while the neutral HID tautomer holds at 170° / 3.77 Å. So:
      · the aspartate nucleophile and dyad aspartate are deprotonated by removing the HD2 carboxyl proton
        (Desmond reads the charge state from the hydrogens present; the residue keeps the name ASP);
      · the histidine base is driven to HID and GUARANTEED there regardless of the input tautomer -
        HIP → remove HE2; HID → unchanged; HIE → remove HE2 and build the Nδ1-H (never left proton-less,
        never left as an unconverted HIE).

    Returns {role: (resname, resnum, action)} for the run log. The structure is rewritten in place only
    when something actually changes.
    """
    _policy = CFG.CATALYTIC_PROTONATION_POLICY
    if not anchors:
        return {}
    try:
        lines = pdb_path.read_text().splitlines()
    except Exception as e:
        if logger:
            logger.warning(f"Protonation enforcement skipped for {job_name}: {e}")
        return {}

    _drop, _insert, _report = set(), {}, {}     # _insert: line index → [new lines to add after it]
    for _role, _num in anchors.items():
        _pol = _policy.get(_role)
        if _pol is None or not _num:
            continue
        _expected = {r.upper() for r in CFG.ROLE_EXPECTED_RESIDUES.get(_pol.get("role_key", ""), set())}
        _cand = [i for i, l in enumerate(lines)
                 if l.startswith(("ATOM", "HETATM")) and l[22:26].strip() == str(_num)]
        if not _cand:
            _report[_role] = ("?", _num, "residue not found")
            continue
        # A residue number alone can collide across chains, or with a water/ligand sharing it. Key on the
        # chain that actually carries this role's residue so only the intended residue is ever modified.
        _chain = next((lines[i][21] for i in _cand
                       if lines[i][17:20].strip().upper() in _expected), lines[_cand[0]][21])
        _idx = [i for i in _cand if lines[i][21] == _chain]
        _rname = lines[_idx[0]][17:20].strip()
        _names = {lines[i][12:16].strip().upper() for i in _idx}

        # CFG-driven identity guard: fail closed if the alignment mapped this role onto a residue of a
        # type CFG.ROLE_EXPECTED_RESIDUES does not permit - stripping/adding atoms there would corrupt it.
        if _expected and _rname.upper() not in _expected:
            _report[_role] = (_rname, _num,
                              f"REFUSED - {_rname} is not a valid {_role} residue "
                              f"(expected {sorted(_expected)}); protonation not enforced")
            if logger:
                logger.warning(f"[QC] {job_name}: {_role} mapped to {_rname}{_num}, which is not in "
                               f"CFG.ROLE_EXPECTED_RESIDUES[{_pol.get('role_key')}] - not modified.")
            continue

        # Roles that carry their standard state (clamp Arg, His155 stabiliser, Trp/Tyr cradle):
        # declared for function + identity check only, never stripped.
        if not _pol.get("enforce"):
            _report[_role] = (_rname, _num, f"{_pol['state']} - not enforced ({_pol['why']})")
            continue

        _strip = {h.upper() for h in _pol.get("strip_H", ())}

        '''
        Histidine base → GUARANTEE HID (Nδ1-H present, Nε2 free) for any input tautomer. Keyed on the
        ROLE's canonical residue (CFG), so a CHARMM-named histidine (HSD/HSE/HSP) - admitted by the
        identity set - takes the His path, not the aspartate carboxyl strip.
        '''
        if _pol.get("residue") == "HIS":
            _has_hd1, _has_he2 = "HD1" in _names, "HE2" in _names
            if _has_hd1 and _has_he2:                       # HIP → HID
                _drop.update(i for i in _idx if lines[i][12:16].strip().upper() == "HE2")
                _act = "HIP → HID (removed HE2)"
            elif _has_hd1:                                  # already HID
                _act = "HID (already correct)"
            elif _has_he2:                                  # HIE → HID: strip HE2, build HD1
                _new = _build_hd1_line([lines[i] for i in _idx])
                if _new is None:
                    _report[_role] = (_rname, _num, "HIE - could not build HD1 (ring atoms missing); left as-is")
                    if logger:
                        logger.warning(f"[QC] {job_name}: base His{_num} is HIE and HD1 could not be "
                                       f"built - left unconverted, verify by hand.")
                    continue
                _drop.update(i for i in _idx if lines[i][12:16].strip().upper() == "HE2")
                _insert.setdefault(min(_idx), []).append(_new)   # anchor to the residue's first atom
                _act = "HIE → HID (removed HE2, added HD1)"
            else:
                _report[_role] = (_rname, _num, "no imidazole protons found; left as-is")
                continue
            _report[_role] = (_rname, _num, _act)
            continue

        # Aspartate nucleophile / dyad acid → deprotonated: remove the HD2 carboxyl proton.
        _removed = [lines[i][12:16].strip() for i in _idx
                    if lines[i][12:16].strip().upper() in _strip]
        _drop.update(i for i in _idx if lines[i][12:16].strip().upper() in _strip)
        _report[_role] = (_rname, _num,
                          f"{_pol['state']} (removed {', '.join(_removed)})" if _removed
                          else f"{_pol['state']} (already correct)")

    if _drop or _insert:
        _out = []
        for i, l in enumerate(lines):
            if i not in _drop:
                _out.append(l)
            _out.extend(_insert.get(i, ()))    # emit inserts even if the anchor line was dropped
        pdb_path.write_text("\n".join(_out) + "\n")
    return _report


def preparation_step(job_name: str, dir_raw: Path, dir_prep_clean: Path, rank: str, prep_index: dict,
                     anchors: dict | None = None):
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
    _prot, _geo = {}, {}
    if prep_success:
        if CFG.PREPWIZARD_ENFORCE_PROTONATION:
            _prot = enforce_catalytic_protonation(final_prep_path, anchors or {}, job_name)
        update_pdb_header(final_prep_path, job_name, rank, raw_source or "Unknown")

        '''
        Measure the pose that is actually handed to MD, and the pose it came from.

        The screen tiers on the Boltz CIF; MD starts here. Preparation moves the geometry - measurably,
        and in the direction that undoes the selection - so the two are recorded side by side and the
        drift is reported. This is instrumentation, not a gate: nothing is rejected on these numbers.
        A candidate whose prepared pose has fallen out of the relaxed NAC envelope is FLAGGED, so that
        it is known before a week of GPU time is spent on it, and never silently dropped.
        '''
        _nuc_res = (anchors or {}).get("Nuc")
        '''
        Three stages: the Boltz CIF (what the screen scored), the RAW PDB (gemmi conversion), and the
        minimised PDB (what MD starts from). The CIF is measured directly rather than assumed equal to
        RAW: the conversion is lossless, but MEASURING it proves that per structure instead of citing an
        aggregate, and lets the figure show CIF and RAW as coincident points - the visual proof that all
        the drift is in minimisation. The CIF mirror is flat: 2_Best_Complexes_CIFs/<job>_model_*.cif.
        '''
        _cif_g = {"sn2_angle": float("nan"), "dist_nuc": float("nan")}
        try:
            _cif_dir = dir_raw.parent.parent / "2_Best_Complexes_CIFs"
            _cif_hit = next(iter(sorted(_cif_dir.glob(f"{job_name}_model_*.cif"))), None) if _cif_dir.is_dir() else None
            if _cif_hit is not None:
                _cif_g = measure_sn2_geometry(_cif_hit, _nuc_res)
        except Exception:
            pass
        _raw_g  = measure_sn2_geometry(raw_pdb_path, _nuc_res)
        _prep_g = measure_sn2_geometry(final_prep_path, _nuc_res)
        _geo = {
            "cif_sn2_angle":   _cif_g["sn2_angle"],
            "cif_dist_nuc":    _cif_g["dist_nuc"],
            "raw_sn2_angle":   _raw_g["sn2_angle"],
            "raw_dist_nuc":    _raw_g["dist_nuc"],
            "prep_sn2_angle":  _prep_g["sn2_angle"],
            "prep_dist_nuc":   _prep_g["dist_nuc"],
            "prep_attack_o":   _prep_g["attack_o"],
            "prep_geom_status": _prep_g["status"],
        }
        _da = _prep_g["sn2_angle"] - _raw_g["sn2_angle"]
        _dd = _prep_g["dist_nuc"] - _raw_g["dist_nuc"]
        _geo["prep_d_angle"] = round(_da, 2) if _da == _da else float("nan")
        _geo["prep_d_dist"]  = round(_dd, 2) if _dd == _dd else float("nan")
        # Did minimisation move the pose out of the envelope the tier was granted on?
        _left_nac = bool(_prep_g["sn2_angle"] == _prep_g["sn2_angle"]
                         and (_prep_g["sn2_angle"] < CFG.NAC_ANGLE_RELAXED
                              or _prep_g["dist_nuc"] > CFG.NAC_DIST_RELAXED))
        _geo["prep_left_nac"] = int(_left_nac)
        if _prep_g["status"] == "ok":
            _msg = (f"  [geom] {job_name}: angle {_raw_g['sn2_angle']:.1f}° → {_prep_g['sn2_angle']:.1f}° "
                    f"({_geo['prep_d_angle']:+.1f}°)   nuc {_raw_g['dist_nuc']:.2f} → "
                    f"{_prep_g['dist_nuc']:.2f} Å ({_geo['prep_d_dist']:+.2f})   attack O {_prep_g['attack_o']}")
            console_info(_msg + ("   [!] LEFT the relaxed NAC envelope" if _left_nac else ""))
        else:
            console_info(f"  [geom] {job_name}: prepared geometry not measurable ({_prep_g['status']})")

    status = "Success" if prep_success else "Prep_Failed"
    return {"job": job_name, "status": status, "rank": rank, "protonation": _prot, **_geo}


# #############################################################################
# SUBSECTION 4B: QM (Jaguar ESP) LIGAND CHARGES  -  opt-in, CFG.ESP_CHARGES_ENABLE / --esp
# #############################################################################
'''
WHY THIS IS HERE, AND NOT IN 07.

The ligand's charges reach the physics through exactly one door: the Desmond SYSTEM BUILD, which Step 06
runs. So the charge set has to exist BEFORE that build. Put this in Step 07 and it would fire after the
trajectory it was meant to influence had already been run.

It is also irrelevant to Step 07 on its own terms: QSite puts the ligand INSIDE the QM region, where DFT
computes its electron density directly and never consults a point charge. ESP charges matter only for
the CLASSICAL regions - the Desmond trajectory and Prime MM-GBSA.

WHAT IT DOES NOT DO. It does not build a system, and it does not run MD or WaterMap: that is Step 06's
work. It writes <ligand>_ESP.mae and stops. Step 06 name-matches each complex's own _ESP.mae by stem,
writes those charges into the built .cms force field, and re-reads them through msys to prove they
reached the MD engine - a charge set that silently failed to apply would be worse than none.
'''

_ESP_BUILD = r"""
import sys
from pathlib import Path
from schrodinger import structure
from schrodinger.application.jaguar.input import JaguarInput

prepared, out_dir, stem, basis, dft = sys.argv[1:6]
out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
_AA = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE',
       'PRO','SER','THR','TRP','TYR','VAL','HID','HIE','HIP','ASH','HOH','NA','CL','SPC','T3P'}

with structure.StructureReader(prepared) as _r:
    st = next(iter(_r), None)
if st is None:
    print('EMPTY_STRUCTURE %s' % prepared); sys.exit(5)
lig = None
for mol in st.molecule:
    if not ({a.pdbres.strip() for a in mol.atom} & _AA) and len(mol.atom) > 2:
        lig = mol.extractStructure(); break
if lig is None:
    print('NO_LIGAND'); sys.exit(2)

chg = sum(a.formal_charge for a in lig.atom)
lig.write(str(out_dir / (stem + '_lig.mae')))

# The charges are fitted on the PREPARED geometry - the structure the MD actually starts from - not on
# an idealised gas-phase optimum. A charge set derived from a different conformer is a charge set for a
# different molecule. icfit=1 fits to the electrostatic potential; the fitted charges sum to the formal
# charge, which is the check that the fit converged.
ji = JaguarInput(name=stem + '_ESP')
ji.setStructure(lig)
ji.setValues({'basis': basis, 'dftname': dft, 'molchg': int(chg), 'multip': 1, 'icfit': 1})
ji.saveAs(str(out_dir / (stem + '_ESP.in')))
print('WROTE_INPUT charge', chg, 'atoms', len(lig.atom))
"""

_ESP_PARSE = r"""
import sys, csv
from pathlib import Path
from schrodinger import structure

out_file, lig_mae, stem, out_dir = sys.argv[1:5]
out_dir = Path(out_dir)
lines = Path(out_file).read_text(errors='ignore').splitlines()

# Jaguar prints the fitted charges as label/charge row PAIRS under one header. The block is read to its
# own end rather than to a fixed line budget: a larger ligand wraps the table over several pairs and a
# fixed cut would silently drop the tail.
labels, charges, capture = [], [], False
for i, ln in enumerate(lines):
    if 'Atomic charges from electrostatic potential' in ln:
        capture = True; continue
    if capture:
        t = ln.strip()
        if t.startswith('Atom') and i + 1 < len(lines) and lines[i + 1].strip().startswith('Charge'):
            labels += ln.split()[1:]; charges += lines[i + 1].split()[1:]
        elif t and not t.startswith('Charge') and labels:
            break
if not labels:
    print('NO_ESP_BLOCK'); sys.exit(3)

vals = [float(c) for c in charges]
with (out_dir / (stem + '_ESP_charges.csv')).open('w', newline='') as fh:
    w = csv.writer(fh); w.writerow(['atom_label', 'esp_charge']); w.writerows(zip(labels, vals))

with structure.StructureReader(lig_mae) as _r:
    st = next(iter(_r), None)
if st is None:
    print('EMPTY_STRUCTURE %s' % lig_mae); sys.exit(5)
if len(st.atom) != len(vals):
    print('ATOM_COUNT_MISMATCH %d vs %d' % (len(st.atom), len(vals))); sys.exit(4)
# The count check is permutation-invariant, so it cannot catch a Jaguar atom reorder relative to the
# .mae; the parsed labels (C1/F5/O6 …) are the only reorder-sensitive evidence. Assert the element of
# each label matches st.atom[i] before binding the charge, so a reorder fails loudly instead of silently
# assigning every atom the wrong charge (PFAS ligands are C/F/O/H/N - single-character element symbols).
for i, (a, q) in enumerate(zip(st.atom, vals)):
    if labels[i][0].upper() != a.element.upper():
        print('ESP_LABEL_MISMATCH idx=%d label=%s elem=%s' % (i, labels[i], a.element)); sys.exit(6)
    a.partial_charge = float(q)
st.write(str(out_dir / (stem + '_ESP.mae')))
print('OK atoms=%d sum=%+.4f' % (len(vals), sum(vals)))
"""

# Jaguar drops its scratch beside the results. Kept only when CFG.ESP_KEEP_SCRATCH is set: 26 files per
# ligand of babel/symtry/restart noise buries the four that are actually the answer.
_ESP_SCRATCH_GLOBS = ("babel*.com", "*.dat", "*.prm", "*.ark", "symtry.*", "restart*.in",
                      "default.mass", "*_tmp.mae", "*_trunc.mae", "*_ESP.01.*", "nbyn.dat", "pbf*",
                      "PID", "runflags", ".write_dir.json")   # Jaguar names these two with no suffix


def _esp_run(script: str, *args: str, timeout: int = 1800):
    """Run a snippet under $SCHRODINGER/run - its interpreter, not ours."""
    return subprocess.run([str(SCHRODINGER_PATH / "run"), "python3", "-c", script, *args],
                          capture_output=True, text=True, timeout=timeout)


def generate_esp_charges(prep_dir: Path, out_dir: Path) -> Path | None:
    """QM (Jaguar ESP) partial charges for every prepared MD ligand. Returns the summary CSV."""
    _jag = SCHRODINGER_PATH / "jaguar"
    if not _jag.exists():
        console_info(f"  ! ESP charges skipped: jaguar not found at {_jag}")
        return None
    _pdbs = sorted(prep_dir.glob("*_Prepared.pdb"))
    if not _pdbs:
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    _basis, _dft = str(CFG.ESP_CHARGE_BASIS), str(CFG.ESP_CHARGE_DFT)
    console_info(f"QM ligand charges - {len(_pdbs)} ligand(s)  |  {_dft}/{_basis}  (Jaguar, icfit=1)")

    _summary, _ok = [], 0
    for _pdb in _pdbs:
        _stem = _pdb.name.replace("_Prepared.pdb", "")
        _r = _esp_run(_ESP_BUILD, str(_pdb), str(out_dir), _stem, _basis, _dft)
        if "WROTE_INPUT" not in _r.stdout:
            console_info(f"  ! {_stem}: input build failed - {(_r.stdout + _r.stderr).strip()[:120]}")
            continue
        _jr = subprocess.run([str(_jag), "run", "-WAIT", "-NOJOBID", f"{_stem}_ESP.in"],
                             cwd=str(out_dir), capture_output=True, text=True, timeout=7200)
        _out = out_dir / f"{_stem}_ESP.out"
        if _jr.returncode != 0 or not _out.exists():
            console_info(f"  ! {_stem}: Jaguar failed (rc={_jr.returncode})")
            continue
        _pr = _esp_run(_ESP_PARSE, str(_out), str(out_dir / f"{_stem}_lig.mae"), _stem, str(out_dir))
        if not _pr.stdout.startswith("OK"):
            console_info(f"  ! {_stem}: {(_pr.stdout + _pr.stderr).strip()[:120]}")
            continue
        _ok += 1
        console_info(f"  ✔ {_stem}  {_pr.stdout.strip()}")
        with (out_dir / f"{_stem}_ESP_charges.csv").open() as _fh:
            for _row in csv.DictReader(_fh):
                _summary.append({"structure": _stem, "atom_label": _row["atom_label"],
                                 "esp_charge": float(_row["esp_charge"])})

    if not _summary:
        return None
    _sum_path = out_dir / "00_ESP_Charges_Summary.csv"
    _utils_mod.atomic_write_csv(pd.DataFrame(_summary), _sum_path)

    if not CFG.ESP_KEEP_SCRATCH:
        _n = 0
        for _g in _ESP_SCRATCH_GLOBS:
            for _f in sorted(out_dir.glob(_g)):
                try:
                    _f.unlink(); _n += 1
                except OSError:
                    pass
        if _n:
            console_info(f"  Removed {_n} Jaguar scratch file(s); the .in/.out/.mae/.csv are kept for audit.")

    console_info(f"  \u2714 {_ok}/{len(_pdbs)} ligand(s) charged  \u2192  {_sum_path.name}")
    console_info("  NEXT (Step 06 - applied automatically, no Maestro step):")
    console_info(f"    1. In {out_dir.name}/ there is one *_ESP.mae per prepared complex ({_ok} total).")
    console_info("    2. Step 06 name-matches each complex's own _ESP.mae by stem and merges the charges.")
    console_info("    3. They are written into the built .cms force field and verified against the MD engine.")
    try:
        _f = plot_esp_alpha_carbon(_summary, out_dir)
        if _f:
            console_info(f"  ESP figure  →  {_utils_mod.deflx_fig_name(_f.name)}")
    except Exception as _e:                                       # noqa: BLE001
        console_info(f"  ! ESP figure skipped: {type(_e).__name__}: {_e}")
    return _sum_path


def plot_esp_alpha_carbon(summary_rows: list, out_dir: Path) -> Path | None:
    """The charge on the carbon the nucleophile actually attacks - the number OPLS4 cannot see.

    An SN2 rate turns on how electrophilic the α-carbon is. The QM charge on that atom runs +0.023 (FA)
    → +0.126 (DFA) → +0.297 (TFA): a thirteen-fold spread across the three substrates. OPLS4 assigns by
    atom type, so in the classical trajectory those three carbons look much alike. The figure exists to
    put that gap where it cannot be missed, because it is buried in a per-atom CSV otherwise.

    The α-carbon is taken as the most positive carbon AFTER the carboxylate carbon (which is always the
    most positive of all, ~+0.6 e, and is not attacked). This assumes exactly two carbons - the
    carboxylate and the α - which holds for FA/DFA/TFA; a ligand with a third carbon is skipped rather
    than mislabelled, because the second-most-positive carbon is then not guaranteed to be the α.
    """
    if not summary_rows:
        return None
    _utils_mod.apply_figure_style(CFG)

    df = pd.DataFrame(summary_rows)
    _rows = []
    for _st, _g in df.groupby("structure"):
        _c = _g[_g["atom_label"].str.match(r"^C\d+$")].sort_values("esp_charge", ascending=False)
        if len(_c) != 2:                          # carboxylate C + α-C only; see docstring
            if logger:
                logger.warning(f"ESP α-carbon figure: {_st} has {len(_c)} carbons (expected 2); "
                               f"skipping - the most-positive-after-carboxylate rule assumes n_C=2.")
            continue
        _alpha = _c.iloc[1]                        # 0 = carboxylate C, 1 = the α-carbon
        _lig = CFG.VIS_LIGAND_SHORT.get(_st.split("_")[-1].lower(), _st.split("_")[-1])   # ligand short via CFG SSOT
        _rows.append({"ligand": _lig, "structure": str(_st), "q_alpha": float(_alpha["esp_charge"]),
                      "n_F": int(_g["atom_label"].str.match(r"^F\d+$").sum())})
    if not _rows:
        return None
    # Dedup per ligand BEFORE sorting for display, preferring a ranked hit over the 0000000_* control
    # (input order is lexicographic groupby order, so a plain dedup would always keep the control).
    d = pd.DataFrame(_rows)
    d["_ctrl"] = d["structure"].str.startswith("0000000")
    # 'structure' (unique per row) as the final tie-break makes the kept row deterministic: without it,
    # two structures sharing a ligand + control-status would resolve by pandas' unstable sort order.
    d = (d.sort_values(["ligand", "_ctrl", "structure"]).drop_duplicates("ligand", keep="first")
           .sort_values("q_alpha"))

    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    _cols = [CFG.VIS_ACCENT["blue"], CFG.VIS_ACCENT["amber"], CFG.VIS_ACCENT["vermillion"]]
    _bar_cols = [_cols[min(i, 2)] for i in range(len(d))]
    _b = ax.barh(d["ligand"], d["q_alpha"],
                 color=_bar_cols,
                 edgecolor=CFG.VIS_INK["dark"], linewidth=0.8, height=0.55, zorder=3)
    # Tint each ligand tick label to its own bar colour so the FA/DFA/TFA axis reads with the bars.
    for _tick, _tc in zip(ax.get_yticklabels(), _bar_cols):
        _tick.set_color(_tc)
    for _r, _v, _nf in zip(_b, d["q_alpha"], d["n_F"]):
        _lx = _v + 0.006 if _v >= 0 else 0.006          # keep negative-bar labels inside the plot
        ax.text(_lx, _r.get_y() + _r.get_height() / 2, f"{_v:+.3f}   ({_nf} F)",
                va="center", ha="left", fontsize=CFG.VIS_FONT_ANNOT, color=CFG.VIS_INK["dark"], zorder=4)
    ax.set_xlabel("QM (Jaguar ESP) charge on the α-carbon - the atom the nucleophile attacks  (e)")
    # Include negatives: FA's α-carbon is slightly negative, so an axis starting at 0 renders its bar
    # entirely off-plot (only the floating label survives).
    _qmin, _qmax = float(d["q_alpha"].min()), float(d["q_alpha"].max())
    ax.set_xlim(min(0.0, _qmin * 1.32), max(_qmax * 1.32, 0.01))
    ax.grid(True, axis="x", alpha=CFG.VIS_GRID_ALPHA, color=CFG.VIS_GRID_COLOUR)
    ax.set_axisbelow(True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        fig.tight_layout()
    _p = out_dir / "01_ESP_Alpha_Carbon_Charge.svg"
    fig.savefig(_p, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    return _p


# =============================================================================
# SECTION 5: STRUCTURE & MOLECULE HANDLING
# =============================================================================

def load_catalytic_anchor_map(prod_dir: Path) -> dict:
    """Per-job catalytic residue numbers taken from 02's dynamic global alignment
    (Mapped_Nucleophile / Mapped_Acid / Mapped_Base columns of the ranked CSV).
    Lets the identity guard anchor on each homolog's actual aligned positions -
    robust to insertions/deletions - instead of a static ±window around the
    canonical reference numbers. Returns {job_name: {"Nuc"/"Acid"/"Base": int}}."""
    rank_csvs = sorted(prod_dir.glob(CFG.GLOB_RANKED_CSV)) or sorted(prod_dir.glob("*_Ranked_*.csv"))
    if not rank_csvs:
        return {}
    try:
        df = pd.read_csv(_utils_mod.latest_by_mtime(rank_csvs), low_memory=False)
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
    only - Step 07 maps catalytic residues from its OWN dynamic global alignment
    (Full_Sequence_Alignment_Map), so it does not consume this file.

    When `anchors` (02's dynamic Mapped_* positions for this job) is supplied the
    search is anchored there - robust to homolog insertions/deletions - rather than
    on a static window around the canonical reference numbers."""
    anchors = anchors or {}
    nuc_ref  = anchors.get("Nuc")  or cfg.DREAM_TEAM_REFS["Nuc"]   # 02 dynamic align, else canonical 110
    acid_ref = anchors.get("Acid") or cfg.DREAM_TEAM_REFS["Acid"]  # else 134
    base_ref = anchors.get("Base") or cfg.DREAM_TEAM_REFS["Base"]  # else 277

    ASP_TYPES = {r.upper() for r in cfg.ROLE_EXPECTED_RESIDUES["Nucleophile"]}
    HIS_TYPES = {r.upper() for r in cfg.ROLE_EXPECTED_RESIDUES["Base_Catalyst"]}

    resnum_to_resname = {}
    resnum_to_atoms   = defaultdict(set)
    try:
        with open(prepared_pdb_path) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    resname = line[17:20].strip()
                    resnum  = int(line[22:26].strip())
                    resnum_to_resname[resnum] = resname
                    resnum_to_atoms[resnum].add(line[12:16].strip().upper())
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

    '''
    QC: the catalytic nucleophile Asp must be DEPROTONATED for the SN2 attack; a protonated
    carboxylic acid cannot attack, so a protonated nucleophile is a catalytically DEAD enzyme.

    The state is read from the HYDROGENS, not from the residue NAME. Schrödinger's preparation keeps
    the name ASP whether or not the carboxyl carries its proton - ASH is an AMBER convention it never
    writes - so a name test never fires and the guard it was supposed to provide is silently absent.
    The observable is the aspartate carboxyl proton itself (HD2, taken from the CFG policy so this and
    the enforcement agree on exactly one atom); there is no Glu HE2 here, the nucleophile is an ASP.

    This is a WARNING, not an edit: enforce_catalytic_protonation strips exactly this hydrogen from the
    nucleophile, so the chemistry is already imposed. The warning's job is to report when the preparation
    has handed over a dead enzyme in the first place, which the enforcement step alone never surfaces.
    '''
    _CARBOXYL_H = {h.upper() for h in cfg.CATALYTIC_PROTONATION_POLICY["Nuc"]["strip_H"]}
    if nuc_found is not None:
        _nuc_atoms = resnum_to_atoms.get(nuc_found, set())
        _present = _nuc_atoms & _CARBOXYL_H
        if _present:
            _msg = (f"[QC] Nucleophile Asp{nuc_found} is PROTONATED "
                    f"(carboxyl H present: {sorted(_present)}) - a deprotonated ASP is required for the "
                    f"SN2 defluorination. The protonation enforcement strips it, but check PropKa/Epik "
                    f"pH: the preparation produced a catalytically dead nucleophile.")
            (logger.warning if logger else print)(_msg)

    '''
    Use the most common non-zero offset (consensus across Nuc/Acid/Base);
    fall back to 0 when all three agree on the reference numbering.
    '''
    offsets = [o for o in (nuc_offset, acid_offset, base_offset) if o is not None and o != 0]
    offset = _Counter(offsets).most_common(1)[0][0] if offsets else 0

    record = {
        "offset":     offset,
        "nuc_found":  nuc_found,
        "acid_found": acid_found,
        "base_found": base_found,
    }

    out_json = prepared_pdb_path.parent / f"{job_name}_index_offset.json"
    _utils_mod.write_json_atomic(out_json, record)   # atomic (temp + rename), like every other data write

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
        '''
        Sorted: two FASTAs carrying the same sequence ID would otherwise overwrite each other in
        seq_map in whatever order the filesystem listed them, so the surviving sequence could differ
        between runs on identical input.
        '''
        for f in sorted(input_data_dir.glob("*.fasta")) + sorted(input_data_dir.glob("*.fa")):
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
                        '''
                        Two prefixed entries collapse to the same bare name but carry
                        different sequences; overwriting would silently hand the wrong
                        sequence to the earlier candidate. Keep the first, warn loudly.
                        '''
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
    for s in sorted(input_data_dir.glob("*.smi")):
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

# -----------------------------------------------------------------------------
# Step 5.1: Molecule Handling
# -----------------------------------------------------------------------------
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


# =============================================================================
# SECTION 6: FIGURE GENERATION ENGINE (PyMOL / PLIP / InteractionMap)
# =============================================================================


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

'''
Per-interaction-type distance ceilings for filtering PLIP output down to the
project's Schrödinger-Maestro criteria (CFG §3 - single source of truth).
PLIP's own internal cutoffs are looser (e.g. H-bond 4.1 Å, salt 5.5 Å); any
PLIP contact whose reported distance exceeds the matching CFG ceiling is
dropped so all three figure engines agree on what counts as a bond.
'''
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


'''
ANGLE CRITERIA - the half of the definition that a distance-only clamp cannot enforce.

A hydrogen bond is not "an N or O within 3.5 Å"; a halogen bond is not "a halogen within 3.5 Å". Both
are defined by GEOMETRY: a donor angle that says the interaction points the right way, and an acceptor
angle that says the lone pair is oriented to receive it. CFG §3-4 carries those angles - the
Schrödinger-Maestro values - and this gate enforces them on PLIP's reported geometry. A distance-only
clamp cannot: it lets PLIP's looser internal angles decide which contacts survive and leaves the CFG
angle constants unapplied. Reading each PLIP angle and holding it to the CFG criterion closes that half.

PLIP reports each angle in its XML, so the criteria can be enforced on its output rather than
re-implemented. Each entry is (xml_tag, minimum, maximum); None means unbounded on that side.

Not listed, deliberately:
  · pi-stacking / pi-cation - PLIP applies its own offset and angle tests and does not expose them in a
    form that can be re-gated without re-deriving the ring geometry, so the distance clamp stands alone.
  · the aromatic-H-bond family (AROM_HB_*) - a Maestro interaction class PLIP does not model at all.
    Those constants describe a criterion no engine in this pipeline can evaluate; they are reference,
    not configuration, and are marked as such in CFG.
'''
_PLIP_ANGLE_RULES: dict = {
    "hbond":   [("don_angle", CFG.THRESHOLD_HB_ANGLE_MIN, None)],
    "halogen": [("don_angle", CFG.HALOGEN_DON_ANGLE_MIN, None),
                ("acc_angle", CFG.HALOGEN_DON_ACC_ANGLE_MIN, None)],
    "water":   [("don_angle", CFG.WATER_BRIDGE_DON_ANGLE_MIN, None),
                ("water_angle", CFG.WATER_BRIDGE_OMEGA_MIN, CFG.WATER_BRIDGE_OMEGA_MAX)],
}


def _plip_angles_ok(node, itype: str) -> bool:
    """True when every CFG angle criterion for this interaction type is satisfied.

    A missing angle in the XML is NOT treated as a pass: the criterion exists because the geometry
    decides whether the contact is real, and an unverifiable contact is not a verified one.
    """
    _rules = _PLIP_ANGLE_RULES.get(itype)
    if not _rules:
        return True
    for _tag, _lo, _hi in _rules:
        _el = node.find(_tag)
        if _el is None or _el.text is None:
            return False
        try:
            _v = float(_el.text)
        except ValueError:
            return False
        if _lo is not None and _v < float(_lo):
            return False
        if _hi is not None and _v > float(_hi):
            return False
    return True

# Metadata cache: job_name → {p: protein_name, l: ligand_name}
METADATA_CACHE: dict = {}

def load_metadata(ext_dir: Path):
    """Populates METADATA_CACHE from any *_Scientific_Data.csv found under ext_dir."""
    try:
        for _csv_path in sorted(ext_dir.rglob("*_Scientific_Data.csv")):
            try:
                df = pd.read_csv(_csv_path)
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

# -----------------------------------------------------------------------------
# --- Figure engine: logger (separate from the extraction logger) ---
# -----------------------------------------------------------------------------
_fig_logger = None

def _fig_log(text=""):
    """Unified logging for figure generation - routes through console_info."""
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

# -----------------------------------------------------------------------------
# --- Figure engine: subprocess runner ---
# -----------------------------------------------------------------------------
def _run_cmd(cmd_list, cwd=None, env=None, timeout=60, log_file=None):
    """Run external command; explicitly kill the subprocess on timeout so no orphans remain."""
    run_env = os.environ.copy()
    if env:
        run_env.update(env)
    run_env["QT_QPA_PLATFORM"]      = "offscreen"
    run_env["LIBGL_ALWAYS_SOFTWARE"] = "1"
    if "DISPLAY" not in run_env:
        run_env["DISPLAY"] = ":99"
    '''
    Cap internal threading in numpy/OpenBabel/MKL so each subprocess uses
    exactly 1 thread; concurrency is controlled at the pool level instead.
    '''
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
        return proc.returncode == 0
    except subprocess.TimeoutExpired:
        if proc:
            proc.kill()
            proc.communicate()
        if f_handle:
            try: f_handle.write(f"\n[TIMEOUT after {timeout}s]\n")
            except Exception: pass
        return False
    except Exception as e:
        if proc:
            try:
                proc.kill()
                proc.wait()          # reap the killed child so it does not linger as a zombie
            except Exception: pass
        if f_handle:
            try: f_handle.write(f"\n[ERROR]: {e}\n")
            except Exception: pass
        return False
    except BaseException:
        # KeyboardInterrupt / SIGTERM: reap the PyMOL/PLIP child before the
        # exception unwinds, so it is not orphaned, then re-raise to abort.
        if proc and proc.poll() is None:
            try:
                proc.kill(); proc.wait()
            except Exception: pass
        raise
    finally:
        # Guarantee the log handle is released on every path, so a failing write in a
        # handler above cannot leak the descriptor under parallel execution.
        if f_handle:
            try: f_handle.close()
            except Exception: pass

# -----------------------------------------------------------------------------
# --- Figure engine: visualisation engines (PyMOL / PLIP) ---
# -----------------------------------------------------------------------------
def _run_pymol(pdb_path, fig_root, log_dir, base_name, lig_name, lig_num, has_f, sw, C):
    if not sw["PyMOL"]: return "Missing"
    out_dir = fig_root / "PyMOL"
    out_dir.mkdir(exist_ok=True, parents=True)

    out_opaque      = out_dir / f"{base_name}_Interaction_Opaque.png"
    out_transparent = out_dir / f"{base_name}_Interaction_Transparent.png"
    out_pse         = out_dir / f"{base_name}_Session.pse"
    pml             = out_dir / f"{base_name}_render.pml"
    log_path        = log_dir / f"{base_name}_PyMOL.log"

    '''
    Precise ligand selection: chain L residue number is authoritative.
    Using AND (not OR) prevents UNK/non-standard protein residues also placed
    in chain L by Schrödinger prep from being selected as ligand atoms.
    '''
    lig_sel = f"(chain L and resi {lig_num})"

    '''
    Cap PyMOL ray-tracing threads: 2 PyMOL workers run concurrently, so each
    should use at most (cpu-2)//2 ray threads to stay within the cpu-2 budget.
    '''
    _max_rt = max(1, CFG.GLOBAL_MAX_WORKERS // 2)

    script = [
        "reinitialize",
        f'load "{pdb_path.resolve()}", complex',
        "remove solvent",
        # Keep polar (N–H / O–H) hydrogens so PyMOL's H-bond detection uses real
        # donor–H···acceptor geometry; drop non-polar C–H to reduce visual clutter.
        "remove (hydro and not (neighbor (elem N+O)))",
        # Lighting - two_sided illuminates inside faces of the binding pocket,
        # preventing the cavity from rendering as a solid black void.
        f"bg_color {CFG.VIS_PYMOL_BG_COLOR}",
        "set ray_trace_mode, 1", "set ray_trace_gain, 0.12",
        "set ray_shadows, 1",
        "set ambient, 0.45", "set direct, 0.55", "set specular, 0.35",
        "set two_sided_lighting, 1",   # illuminate interior surface faces
        "set antialias, 2", "set orthoscopic, on",
        "set ambient_occlusion_mode, 0",  # AO darkens cavities - disable to keep pocket visible
        "set surface_quality, 1",
        f"set max_threads, {_max_rt}",
        # Selections
        f"select ligand, {lig_sel}",
        f"select interacting, byres polymer.protein within {C['DIST_F_CONTACT']} of ligand",
        "hide all",
        # Pocket-local surface - grey, 50% transparent so ligand sticks show through.
        # Restricted to residues near the ligand: the far protein is off-frame after
        # the zoom/clip below, so rendering its surface is wasted ray time.
        f"show surface, byres (polymer.protein within {CFG.VIS_PYMOL_SURFACE_RADIUS} of ligand)",
        "color gray65, polymer.protein",
        f"set transparency, {CFG.VIS_PYMOL_SURFACE_TRANSPARENCY}",
        # Protein cartoon underneath surface
        "show cartoon, polymer.protein",
        "color gray70, polymer.protein",
        f"set cartoon_transparency, {CFG.VIS_PYMOL_CARTOON_TRANSPARENCY}",
        "set cartoon_fancy_helices, 1", "set cartoon_fancy_sheets, 1",
        # Interacting residues - thin grey sticks, element-coloured heteroatoms
        "show sticks, interacting",
        "color gray85, interacting and elem C",
        "color red,    interacting and elem O",
        "color blue,   interacting and elem N",
        "color yellow, interacting and elem S",
        "set stick_radius, 0.15, interacting",
        # Ligand - thick bright orange sticks; highly visible against grey surface
        # orange chosen for maximum contrast on grey without clashing with interaction colours
        "show sticks, ligand",
        "color tv_orange, ligand",
        "color tv_green, ligand and elem F",
        "color red,      ligand and elem O",
        "color blue,     ligand and elem N",
        "set stick_radius, 0.32, ligand",
        # H-bonds - yellow dashes
        f"dist hbonds, ligand, interacting, mode=2, cutoff={C['DIST_HBOND']}",
        "color yellow, hbonds",
        "set dash_gap, 0.20", "set dash_width, 3.5", "set dash_color, yellow",
        # Salt bridges - cyan dashes (cutoff from CFG.THRESHOLD_SALT_BRIDGE - Maestro)
        f"dist salt_bridge, (ligand and (name O*,N*)), (interacting and (name N*,O*)), mode=2, cutoff={C['DIST_SALT']}",
        "color cyan, salt_bridge",
        # Residue labels - white, bold, readable on grey background
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
      lig_atom (nearest ligand atom), prot_atom (''), centre (3-D protcoo)
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
    _ANGLE_REJECTS: dict = {}        # contacts PLIP found but the CFG geometry criteria reject

    def _register(resnr, restype, reschain, dist, itype, protcoo_3d, ligcoo_3d, node=None):
        key = (reschain, int(resnr), restype)
        '''
        Clamp PLIP's looser internal cutoffs to the project's CFG (Maestro) criteria - on DISTANCE and
        on ANGLE - so PLIP, PyMOL and InteractionMap all agree on what counts as a bond, and so the
        geometry half of each definition is actually applied rather than left to PLIP's own defaults.
        '''
        ceil = _cfg_cut.get(itype)
        if ceil is not None and dist > ceil:
            return
        if node is not None and not _plip_angles_ok(node, itype):
            _ANGLE_REJECTS[itype] = _ANGLE_REJECTS.get(itype, 0) + 1
            return
        la, _  = _closest_lig_atom(ligcoo_3d)
        # C–F bonds are leaving groups in FAcD SN2; exclude F from H-bond / halogen contacts
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
                      hb.find("reschain").text, dist, "hbond", protcoo, ligcoo, node=hb)

        for hx in iact.findall("./halogen_bonds/halogen_bond"):
            protcoo = _plip_coo(hx.find("protcoo"))
            ligcoo  = _plip_coo(hx.find("ligcoo"))
            dist    = float(hx.find("dist").text)
            _register(hx.find("resnr").text, hx.find("restype").text,
                      hx.find("reschain").text, dist, "halogen", protcoo, ligcoo, node=hx)

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
                      protcoo, ligcoo, node=wb)

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

        '''
        bs_residues with contact="True" and no typed interaction → "contact"
        Index protein residues by (chain, resnum) once so the lookup below is O(1)
        rather than a linear scan of `pro` per contact residue.
        '''
        _pro_by_cr = {(k[0], k[1]): k for k in pro}
        for bsr in bs.findall("./bs_residues/bs_residue"):
            if bsr.get("contact", "False") != "True":
                continue
            text   = bsr.text.strip()           # e.g. "150A"
            chain  = text[-1]
            resnr  = text[:-1]
            restype = bsr.get("aa", "UNK")
            key    = (chain, int(resnr), restype)
            if key not in contacts_by_key:
                # No typed interaction - approximate position from PDB
                pkey = _pro_by_cr.get((chain, int(resnr)))
                if pkey:
                    protcoo = pro[pkey]["center"]
                    la, _ = _closest_lig_atom(protcoo)
                    contacts_by_key[key] = {
                        "key": key, "resname": restype, "resnum": int(resnr),
                        "chain": chain, "dist": float(bsr.get("min_dist", CFG.PLIP_CONTACT_FALLBACK_DIST)),
                        "itype": "contact", "is_hbond": False, "is_salt": False,
                        "lig_atom": la, "prot_atom": "", "center": protcoo,
                    }

    '''
    A filter that removes contacts silently is indistinguishable from a filter that is not running. The
    count of contacts PLIP found and the CFG geometry rejected is reported, so tightening a criterion in
    CFG has a visible consequence rather than an invisible one.
    '''
    if _ANGLE_REJECTS:
        _msg = ", ".join(f"{_n} {_t}" for _t, _n in sorted(_ANGLE_REJECTS.items()))
        console_info(f"    PLIP contacts rejected on the CFG angle criterion: "
                     f"{sum(_ANGLE_REJECTS.values())} ({_msg})")

    return sorted(contacts_by_key.values(), key=lambda x: x["dist"])


def _im_render_diagram(lig_2d, lig, res_2d, contacts, out_png, mode="distance", role_resnums=None):
    """Shared 2D interaction renderer for InteractionMap (mode='distance') and PLIP (mode='plip').

    The pocket circle hugs the ligand; residues sit outside it (upper arc when few, full ring when
    many - see _im_project). Each residue box shows its name and interaction type and is coloured by
    catalytic role (CFG.ACTIVE_SITE_ROLE_GROUP_COLOUR) when its number is in role_resnums, else the
    neutral contact colour. The measured distance rides each line - for InteractionMap always, for
    PLIP only where CFG.INTERACTION_DIAGRAM_STYLE marks it. Ligand atoms are labelled on the atom; the
    legend is one row of role swatches. Emitted in CFG.VIS_FIGURE_FORMAT (true-vector SVG by default).
    """
    _chrome = CFG.INTERACTION_DIAGRAM_CHROME
    _bond = CFG.BOND_TYPE_COLOUR
    _itype_col = {
        "hbond": _bond["H-Bond"], "salt": _bond["Salt Bridge"], "halogen": _bond["Halogen"],
        "hydrophobic": _bond["Hydrophobic"], "contact": _bond["Hydrophobic"],
        **CFG.INTERACTION_DIAGRAM_EXTRA_COLOUR,
    }
    _ITYPE_STYLE = {_k: (_lw, _ls, _itype_col[_k], _sd)
                    for _k, (_lw, _ls, _sd) in CFG.INTERACTION_DIAGRAM_STYLE.items()}
    _DIST_COL = {_k: (_itype_col[_k], CFG.VIS_INK["white"], _itype_col[_k])
                 for _k in ("hbond", "arom_hbond", "halogen", "salt", "water")}
    _BADGE = {
        "hbond": "H-bond", "arom_hbond": "arom. H-bond", "halogen": "halogen bond",
        "salt": "salt bridge", "water": "water bridge", "pistack": "π-stack",
        "pication": "π-cation", "hydrophobic": "hydrophobic", "contact": "contact",
    }
    role_resnums = role_resnums or {}

    circ_r = (max(float(np.linalg.norm(p)) for p in lig_2d) + 0.55) if len(lig_2d) else 4.2
    _rc = list(res_2d.values())
    _bhw, _bhh = 0.46, 0.24
    if _rc:
        _xlo = min(min(v[0] - _bhw for v in _rc), -circ_r) - 0.12
        _xhi = max(max(v[0] + _bhw for v in _rc),  circ_r) + 0.12
        _yhi = max(max(v[1] + _bhh for v in _rc),  circ_r) + 0.12
        _content_bot = min(min(v[1] - _bhh for v in _rc), -circ_r)
    else:
        _xlo, _xhi, _yhi, _content_bot = -circ_r - 0.12, circ_r + 0.12, circ_r + 0.12, -circ_r
    _ylo = _content_bot - 0.10

    fig = Figure(figsize=(12, 12), facecolor="white")
    canvas = FigureCanvasAgg(fig)
    ax = fig.add_subplot(111, aspect="equal")
    ax.axis("off")
    ax.set_xlim(_xlo, _xhi)
    ax.set_ylim(_ylo, _yhi)

    ax.add_patch(_mpatches.Circle((0, 0), circ_r, color=_chrome["pocket_fill"], zorder=0, alpha=0.6))
    ax.add_patch(_mpatches.Circle((0, 0), circ_r, color=_chrome["pocket_edge"], fill=False,
                                  linewidth=1.2, linestyle="--", zorder=0, alpha=0.4))
    ax.text(0, -(circ_r - 0.35), "Binding Pocket", ha="center", fontsize=CFG.VIS_FONT_ANNOT,
            color=_chrome["annotation"], style="italic", zorder=1)

    name2idx = {a["name"]: i for i, a in enumerate(lig)}
    for c in contacts:
        rpos = res_2d[c["key"]]
        li = name2idx.get(c["lig_atom"]["name"], 0) if c.get("lig_atom") else 0
        lap = lig_2d[li]
        itype = c.get("itype", "hbond" if c.get("is_hbond") else "salt" if c.get("is_salt") else "contact")
        lw, ls, col, show_dist = _ITYPE_STYLE.get(itype, _ITYPE_STYLE["contact"])
        ax.plot([lap[0], rpos[0]], [lap[1], rpos[1]], lw=lw, ls=ls, color=col, alpha=0.9,
                zorder=2, solid_capstyle="round")
        if (show_dist or mode == "distance") and c.get("dist"):
            mx, my = (lap[0] + rpos[0]) / 2, (lap[1] + rpos[1]) / 2
            tc, fc, ec = _DIST_COL.get(itype, (col, CFG.VIS_INK["white"], col))
            ax.text(mx, my, f"{c['dist']:.1f} Å", fontsize=CFG.VIS_FONT_ANNOT, ha="center", va="center",
                    fontweight="bold", color=tc,
                    bbox=dict(fc=fc, ec=ec, alpha=0.88, boxstyle="round,pad=0.2", linewidth=0.8), zorder=8)

    for i, a in enumerate(lig):
        for j in range(i + 1, len(lig)):
            if np.linalg.norm(a["pos"] - lig[j]["pos"]) < CFG.LIG_COVALENT_BOND_DIST:
                p1, p2 = lig_2d[i], lig_2d[j]
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=CFG.VIS_INK["ink_deep"], lw=2.8,
                        solid_capstyle="round", zorder=4, alpha=0.85)
    for i, a in enumerate(lig):
        xy = lig_2d[i]
        col = _IM_ELEM_COLORS.get(a["elem"], _IM_ELEM_COLORS["other"])
        r = {"C": 0.17, "N": 0.19, "O": 0.19, "F": 0.17, "S": 0.22}.get(a["elem"], 0.15)
        ax.add_patch(_mpatches.Circle(xy, r, color=col, zorder=5, ec="white", lw=1.4))
        ax.text(xy[0], xy[1], a["elem"], ha="center", va="center", fontsize=CFG.VIS_FONT_ANNOT - 1,
                color="white", fontweight="bold", zorder=6)

    bw, bh = 0.82, 0.38
    for c in contacts:
        rpos = res_2d[c["key"]]
        _grp = CFG.ACTIVE_SITE_ROLE_GROUP.get(role_resnums.get(_im_resnum(c.get("resnum")), ""), "")
        col = CFG.ACTIVE_SITE_ROLE_GROUP_COLOUR.get(_grp, _chrome["residue_default"])
        label = f"{c['resname']}{c['resnum']}"
        itype = c.get("itype", "contact")
        badge = _BADGE.get(itype, "contact")
        if itype == "contact" and c.get("quality") in ("bad", "severe"):
            badge = f"contact ({c['quality']})"
        ax.add_patch(_FancyBboxPatch((rpos[0] - bw / 2 + 0.025, rpos[1] - bh / 2 - 0.025), bw, bh,
                     boxstyle="round,pad=0.035", facecolor=_chrome["label_box"], alpha=0.22, zorder=5, linewidth=0))
        ax.add_patch(_FancyBboxPatch((rpos[0] - bw / 2, rpos[1] - bh / 2), bw, bh, boxstyle="round,pad=0.035",
                     facecolor=col, edgecolor="white", linewidth=1.2, alpha=0.95, zorder=6))
        ax.text(rpos[0], rpos[1] + 0.075, label, ha="center", va="center", fontsize=CFG.VIS_FONT_TICK - 1,
                fontweight="bold", color="white", zorder=7)
        ax.text(rpos[0], rpos[1] - 0.085, badge, ha="center", va="center", fontsize=CFG.VIS_FONT_ANNOT,
                color="white", alpha=0.92, zorder=7)

    def _swatch(colour, lbl):
        _da = _DrawingArea(16, 11, 0, 0)
        _da.add_artist(_Rectangle((0, 1), 13, 9, fc=colour, ec="none"))
        return _HPacker(children=[_da, _TextArea(lbl, textprops=dict(fontsize=CFG.VIS_FONT_LEGEND))],
                        align="center", pad=0, sep=3)
    _roles = [("Nucleophile", "Nucleophile"), ("Acid/base catalysis", "Acid / base catalysis"),
              ("Carboxylate clamp", "Carboxylate clamp"), ("Fluoride pocket", "Fluoride pocket")]
    _sw = [_swatch(CFG.ACTIVE_SITE_ROLE_GROUP_COLOUR[g], lbl) for g, lbl in _roles]
    _sw.append(_swatch(_chrome["residue_default"], "Other contact"))
    _leg = _AnchoredOffsetbox(loc="upper center",
                              child=_HPacker(children=_sw, pad=0, sep=14, align="center"),
                              frameon=True, pad=0.4, borderpad=0,
                              bbox_to_anchor=((_xlo + _xhi) / 2, _ylo), bbox_transform=ax.transData)
    _leg.patch.set_edgecolor(CFG.VIS_INK["palest"])
    _leg.patch.set_linewidth(1.0)
    _leg.set_clip_on(False)
    ax.add_artist(_leg)

    _out = out_png.with_suffix("." + CFG.VIS_FIGURE_FORMAT)
    fig.savefig(str(_out), dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight", pad_inches=0.05,
                facecolor="white", edgecolor="none")
    del fig, canvas


def _run_plip(pdb_path, fig_root, log_dir, base_name, sw, C, role_resnums=None):
    """Run PLIP (XML only), parse typed interactions, render matplotlib diagram."""
    if not sw["PLIP"]: return "Missing"
    out_dir = fig_root / "PLIP"
    out_dir.mkdir(exist_ok=True, parents=True)
    log_path = log_dir / f"{base_name}{CFG.SUFFIX_PLIP_LOG}"

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

        xml_files = sorted(tmp_dir.glob("*.xml"))
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
            out_png = out_dir / f"{base_name}_PLIP.svg"
            _im_render_diagram(lig_2d, lig, res_2d, contacts, out_png, mode="plip", role_resnums=role_resnums)
            rendered = True

        return "Success" if rendered else "Failed"
    except Exception as _e:
        _fig_log(f"  [!] PLIP render failed for {base_name}: {type(_e).__name__}: {_e}")
        return "Failed"
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# -----------------------------------------------------------------------------
# --- Figure engine: matplotlib interaction diagram (pure Python, no PyMOL) ---
# -----------------------------------------------------------------------------

_IM_CONTACT_DIST  = CFG.CATALYTIC_DIST_CUTOFF
_IM_HBOND_DIST    = CFG.THRESHOLD_HB_DIST_MAX
_IM_SALT_DIST     = CFG.THRESHOLD_SALT_BRIDGE   # Maestro salt-bridge cutoff (single source)
_IM_HBOND_DONORS  = {"N","NH1","NH2","NE","NE2","ND1","ND2","NZ","OG","OG1","OH","NE1"}
_IM_HBOND_ACC     = {"O","OD1","OD2","OE1","OE2","OG","OG1","OH","NE2","ND1",
                     "O1","O2","O3","O4","O5","O6"}
# Maestro H-bond / aromatic-H-bond / contact criteria - all sourced from CFG (single source).
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
_IM_ELEM_COLORS = CFG.LIGAND_ELEMENT_COLOUR
_IM_HYDROPHOBIC = {"LEU","ILE","VAL","PHE","TRP","PRO","MET","ALA","GLY","CYS"}

# Catalytic-role colouring for the interaction diagrams: residue number -> role -> CFG role-group colour.
_IM_ROLE_MAPPED_COLS = {
    "Mapped_Nucleophile": "Nuc", "Mapped_Acid": "Acid", "Mapped_Base": "Base",
    "Mapped_Clamp1": "Carb1", "Mapped_Clamp2": "Carb2",
    "Mapped_Stabiliser_H": "Stab_H", "Mapped_Stabiliser_W": "Stab_W", "Mapped_Stabiliser_Y": "Stab_Y",
}


def _im_resnum(value) -> "int | None":
    _m = _re.search(r"(\d+)\s*$", str(value))
    return int(_m.group(1)) if _m else None


def _im_role_resnums(run_dir: Path) -> dict:
    """{job_name: {residue_number: role_key}} from the ranked CSV Mapped_* columns, so each
    interaction-diagram residue box can be coloured by its catalytic role."""
    prod = run_dir / "1_Boltz2_Production"
    rank_csvs = sorted(prod.glob(CFG.GLOB_RANKED_CSV)) or sorted(prod.glob("*_Ranked_*.csv"))
    if not rank_csvs:
        return {}
    try:
        df = pd.read_csv(_utils_mod.latest_by_mtime(rank_csvs), low_memory=False)
    except Exception:
        return {}
    if "job_name" not in df.columns:
        return {}
    _cols = [(c, r) for c, r in _IM_ROLE_MAPPED_COLS.items() if c in df.columns]
    out = {}
    for _, row in df.iterrows():
        _m = {}
        for col, role in _cols:
            _n = _im_resnum(row.get(col))
            if _n is not None:
                _m[_n] = role
        out[str(row["job_name"])] = _m
    return out
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

        '''
        3) Salt bridge (formally charged groups within CFG cutoff). Residue-aware:
        only ASP/GLU carboxylate O and ARG/LYS/HIS(+) cationic N count - atom-name
        prefix alone would mislabel neutral ASN/GLN amide (ND2/OD1, NE2/OE1) and
        backbone atoms as ionic. Ligand partner must be an O/N (elem-gated).
        '''
        if itype is None:
            _rn = str(res.get("resname", "")).upper()
            _anion  = _rn in {"ASP", "GLU"} and pname in {"OD1", "OD2", "OE1", "OE2"}
            _cation = _rn in {"ARG", "LYS", "HIS", "HIP", "HIE", "HID"} and \
                      pname in {"NH1", "NH2", "NE", "NZ", "ND1", "NE2"}
            charged = _anion or _cation
            if mind <= _IM_SALT_DIST and la["elem"] in {"O", "N"} and charged:
                itype, disp = "salt", mind

        # 4) Steric contact - Bondi vdW ratio (excludes H-bonds / salt above)
        if itype is None:
            ra = _IM_VDW.get(la["elem"], _IM_VDW_DEF)
            # True parsed element (ZN/SE/MG/FE safe); fall back to first letter only if absent.
            rb = _IM_VDW.get(res["elems"].get(pname, pname[0].upper()).upper(), _IM_VDW_DEF)
            ratio = mind / (ra + rb)
            quality = ("severe" if ratio < _IM_CONTACT_UGLY else
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
    """SVD-project the ligand onto its principal plane and place each interacting residue outside the
    pocket. Few residues (<= CFG.INTERACTION_DIAGRAM_UPPER_ARC_MAX) spread across the upper arc so the
    bottom stays clear for the legend; more take the full ring with collision separation."""
    lig_c3 = np.mean([a["pos"] for a in lig], axis=0)
    L = np.array([a["pos"] - lig_c3 for a in lig])
    if len(L) >= 2:
        _, _, Vt = np.linalg.svd(L, full_matrices=False)
        u1, u2 = Vt[0], Vt[1]
    else:
        u1, u2 = np.array([1., 0., 0.]), np.array([0., 1., 0.])
    lig_2d = np.array([[np.dot(a["pos"] - lig_c3, u1),
                        np.dot(a["pos"] - lig_c3, u2)] for a in lig])
    span = np.max(np.linalg.norm(lig_2d, axis=1)) if len(lig_2d) else 1.0
    lig_2d *= CFG.INTERACTION_DIAGRAM_LIGAND_SCALE / max(span, 0.5)
    lig_r = float(np.max(np.linalg.norm(lig_2d, axis=1))) if len(lig_2d) else 1.0

    res_2d = {}
    if not contacts:
        return lig_2d, res_2d

    upper = len(contacts) <= CFG.INTERACTION_DIAGRAM_UPPER_ARC_MAX
    # Residues hug the pocket edge; the full ring only grows if too many to fit at that radius.
    ring_r = lig_r + 1.1 if upper else max(
        lig_r + 1.1, len(contacts) * CFG.INTERACTION_DIAGRAM_MIN_RES_SEP / (2 * np.pi) + 0.3)
    for c in contacts:
        v = c["center"] - lig_c3
        x2, y2 = np.dot(v, u1), np.dot(v, u2)
        d = np.hypot(x2, y2) or ring_r
        res_2d[c["key"]] = np.array([ring_r * x2 / d, ring_r * y2 / d])

    keys = list(res_2d.keys())
    if upper:
        nat = {k: np.degrees(np.arctan2(res_2d[k][1], res_2d[k][0])) for k in keys}
        order = sorted(keys, key=lambda k: nat[k], reverse=True)
        _hi, _lo, _n = 200.0, -20.0, len(order)
        for i, k in enumerate(order):
            ang = np.radians(_hi - ((i + 0.5) / _n) * (_hi - _lo))
            res_2d[k] = ring_r * np.array([np.cos(ang), np.sin(ang)])
    else:
        _sep = CFG.INTERACTION_DIAGRAM_MIN_RES_SEP
        for _ in range(90):
            moved = False
            for i in range(len(keys)):
                for j in range(i + 1, len(keys)):
                    a, b = res_2d[keys[i]], res_2d[keys[j]]
                    dd = float(np.linalg.norm(a - b))
                    if dd < _sep:
                        push = (a - b) / max(dd, 0.01) * (_sep - dd) * 0.5
                        res_2d[keys[i]] = a + push
                        res_2d[keys[j]] = b - push
                        for k in (keys[i], keys[j]):
                            r = float(np.linalg.norm(res_2d[k]))
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


def _draw_interaction_diagram(pdb_path, fig_root, base_name, C, role_resnums=None):
    """Render a publication-quality 2D protein–ligand interaction map (no PyMOL)."""
    out_dir = fig_root / "InteractionMap"
    out_dir.mkdir(exist_ok=True, parents=True)
    out_png = out_dir / f"{base_name}_interaction.svg"
    try:
        lig, pro, lig_hb = _im_parse_pdb(pdb_path)
        if not lig:
            return "Failed"
        contacts = _im_find_contacts(lig, pro, lig_hb)
        if not contacts:
            return "Failed"
        '''
        Detector already assigns c['itype'] (hbond / arom_hbond / salt / contact).
        Only refine a plain steric contact on a hydrophobic residue for display.
        '''
        for c in contacts:
            if c["itype"] == "contact" and c["resname"] in _IM_HYDROPHOBIC:
                c["itype"] = "hydrophobic"
        lig_2d, res_2d = _im_project(lig, contacts)
        lig_2d = _im_separate_atoms(lig_2d)
        _im_render_diagram(lig_2d, lig, res_2d, contacts, out_png, mode="distance", role_resnums=role_resnums)
        return "Success"
    except Exception as _e:
        _fig_log(f"  [!] InteractionMap failed for {base_name}: {_e}")
        return "Failed"


# -----------------------------------------------------------------------------
# --- Figure engine: software manager ---
# -----------------------------------------------------------------------------
class SoftwareManager:
    """Handles verification and automated installation of visual tools."""

    def __init__(self):
        self.status = {}

    def check_all(self):
        _fig_log(f"\n{ConsoleColours.BOLD}>>> Checking Visualisation Engine Availability{ConsoleColours.ENDC}")

        '''
        1. PyMOL - prefer the binary in the running conda env over the system one
        (system /usr/bin/pymol on Ubuntu 24.04 uses python3-pymol 2.5 which
        calls `from imp import find_module`; imp is absent in Python 3.12+)
        '''
        _conda_pymol = Path(sys.executable).parent / "pymol"
        if _conda_pymol.exists():
            self.status["PyMOL"] = str(_conda_pymol)
        else:
            self.status["PyMOL"] = shutil.which("pymol")
        if not self.status["PyMOL"]:
            self.status["PyMOL"] = self._attempt_install_pymol()

        # 2. PLIP - try binary first, then check if module is importable
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
            import plip  # noqa: F401 - availability probe; module presence is the signal
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
            _fig_log(f"  {ConsoleColours.FAIL}✘{ConsoleColours.ENDC} PyMOL      : Conda install failed.")
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
            _fig_log(f"  {ConsoleColours.FAIL}✘{ConsoleColours.ENDC} PLIP       : Pip install failed.")
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
            '''
            The stem is scanned for EVERY ligand key it contains, not just the first. Taking the first
            match and breaking silently would pick one ligand out of an ambiguous filename and give no
            sign that another was equally valid - the kind of choice that is only ever discovered when
            the wrong ligand turns up in a result table.
            '''
            _parts = pdb_path.stem.split("_")
            _matches = [_parts[_i + 1] for _i, _p in enumerate(_parts[:-1])
                        if f"{_p}_{_parts[_i + 1]}" in smi_map]
            if _matches:
                res_name = _matches[0]          # e.g. 'TFA'
                if len(_matches) > 1:
                    console_info(f"  [!] {pdb_path.name}: filename matches {len(_matches)} ligand keys "
                                 f"({', '.join(_matches)}); using '{_matches[0]}'. Check the naming.")
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

# -----------------------------------------------------------------------------
# --- Figure engine: phase-2 entry point ---
# -----------------------------------------------------------------------------
def run_figure_generation(run_dir: Path, ext_dir: Path):
    """
    Phase 2: Run PyMOL and PLIP rendering on all PDB files under ext_dir.
    Called during Phase 2, after the Top-N extraction has completed.
    """
    global _fig_logger

    if not ext_dir.exists():
        console_info("Phase 2: Extraction directory not found - skipping figure generation.")
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
        (run_dir / n for n in ["1_Boltz2_Production/1_Input_Data",
                                "1_Boltz2_Production/1_Input_FASTA_and_SMILES"]
         if (run_dir / n).exists()), None)
    _, _smi_map_fig, _, _ = load_reference_data(_input_data_dir) if _input_data_dir else ({}, {}, 0, 0)

    # Resolve figure constants
    C = _fig_constants()

    # Residue-number -> catalytic-role map per job, for role-coloured interaction-diagram boxes.
    _role_by_job = _im_role_resnums(run_dir)

    def _roles_for(stem: str) -> dict:
        _j = max((j for j in _role_by_job if stem.startswith(j)), key=len, default=None)
        return _role_by_job.get(_j, {}) if _j else {}

    sw_mgr    = SoftwareManager()
    sw_status = sw_mgr.check_all()
    workers   = CFG.GLOBAL_MAX_WORKERS

    _fig_log(f"\n{ConsoleColours.BOLD}Multi-Engine Figure Generation Pipeline (Concurrency: {workers}){ConsoleColours.ENDC}")

    '''
    Build separate task queues:
    - subprocess_tasks: PyMOL + PLIP (spawn external processes - thread-safe)
    - imap_tasks: InteractionMap (matplotlib Agg - NOT thread-safe, run serially)
    '''
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
            subprocess_tasks.append(("PLIP",  _run_plip,  task_base + (sw_status, C, _roles_for(pdb.stem))))
            imap_tasks.append((pdb.stem, pdb, fig_root))

    _plip_tasks  = [t for t in subprocess_tasks if t[0] == "PLIP"]
    _pymol_tasks = [t for t in subprocess_tasks if t[0] == "PyMOL"]
    '''
    PyMOL workers capped at 2: ray-tracing is CPU-bound; high concurrency causes
    CPU saturation, massively inflated render times, and timeout-related orphan processes.
    '''
    _pymol_workers = min(2, max(1, workers))
    _fig_log(
        f"Subprocess tasks: PLIP={len(_plip_tasks)} ({workers} workers) | "
        f"PyMOL={len(_pymol_tasks)} ({_pymol_workers} workers)  |  "
        f"InteractionMap={len(imap_tasks)} (serial)"
    )

    results    = {}

    def _collect(fut_map, pool_label, total):
        # Per-pool progress: PLIP and PyMOL each report against their own task
        # count (not the combined subprocess total), so the denominator is not inflated.
        done = 0
        for f in as_completed(fut_map):
            engine, cname = fut_map[f]
            status = f.result()
            if cname not in results:
                results[cname] = {}
            results[cname][engine] = status
            done += 1
            if done % 5 == 0 or done == total:
                if sys.stdout.isatty():
                    print(f"\r  [{pool_label}] {done}/{total} tasks ...", end="", flush=True)
                elif done == total:
                    print(f"  [{pool_label}] {done}/{total} tasks finished.", flush=True)

    # Phase A1: PLIP - I/O-bound, run at full worker count
    with ThreadPoolExecutor(max_workers=workers) as executor:
        _collect({executor.submit(t[1], *t[2]): (t[0], t[2][3]) for t in _plip_tasks}, "PLIP", len(_plip_tasks))

    # Phase A2: PyMOL - CPU-bound ray-tracing, limited to 2 concurrent processes
    with ThreadPoolExecutor(max_workers=_pymol_workers) as executor:
        _collect({executor.submit(t[1], *t[2]): (t[0], t[2][3]) for t in _pymol_tasks}, "PyMOL", len(_pymol_tasks))

    if sys.stdout.isatty():
        print("\n")

    # Phase B: InteractionMap - serial (matplotlib Agg is not thread-safe)
    _fig_log(f"\n  Running {len(imap_tasks)} InteractionMap diagrams (serial)...")
    for idx, (cname, pdb, fig_root) in enumerate(imap_tasks, 1):
        status = _draw_interaction_diagram(pdb, fig_root, cname, C, _roles_for(cname))
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

    # Final cleanup - remove tool-side non-image files and Figure_Logs
    _keep_ext = {".png", ".pse", ".jpg", ".jpeg", ".svg", ".html"}
    for d in dirs:
        fig_dir = d / "Figures"
        if not fig_dir.exists():
            continue
        for tool_dir in sorted(fig_dir.iterdir()):
            if not tool_dir.is_dir() or tool_dir.name == "Figure_Logs":
                continue
            for item in list(sorted(tool_dir.iterdir())):
                if item.is_file() and item.suffix.lower() not in _keep_ext:
                    item.unlink(missing_ok=True)
        logs = fig_dir / "Figure_Logs"
        if logs.exists():
            shutil.rmtree(logs, ignore_errors=True)

    _fig_log(f"\n\n{ConsoleColours.BOLD}{ConsoleColours.OKGREEN}                ✔ Figure Generation Complete                {ConsoleColours.ENDC}\n")




def setup_logging_extraction(output_dir: Path) -> Path:
    """Initialise the Top-N extraction log (separate handler from the prep log)."""
    global logger
    log_file = output_dir / "00_TopN_and_Preparation.log"
    logger = _utils_mod.setup_logging(log_file, logger_name="extraction", mode="a", timestamp=True)
    return log_file


def prep_and_convert_phase(args):
    """Phase 1 - CIF->PDB conversion + PrepWizard preparation of the MD-ready cohort."""

    root = Path.cwd()
    run_path = root / args.run_folder_name
    if not run_path.exists(): sys.exit(f"Run path missing: {run_path}")

    # Paths - one consolidated Step-05 folder for both phases (prep + extraction)
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

    '''
    SECTION 18 gate: restrict to the MD-ready cohort so only the MD_Selected complexes are converted
    and prepared, not the whole predicted library. Of the control jobs (ID 0000000_*) only the single
    CFG-designated MD control is kept - the reference structure paired with a CFG.CONTROL_MD_LIGANDS
    ligand (3R3U × fluoroacetate); the other reference systems (3R3U-DFA/TFA, DeHa4-*) are scored and
    tiered upstream but never simulated, so they are not converted here. Falls back to all CIFs when
    the ranked CSV carries no MD_Selected column.
    '''
    _md_jobs = load_md_selected_jobs(prod_dir)
    if _md_jobs:
        jobs = [(jn, cif) for (jn, cif) in jobs
                if jn in _md_jobs or (jn.startswith(CFG.CONTROL_JOB_PREFIX) and _is_reference_control(jn))]

    _utils_mod.print_script_banner(
        "05_TopN_and_PDB_Preparation_DeFluorX.py",
        "PDB Conversion & Preparation  ·  Top-N Selection & Delivery  ·  Figure Generation",
    )
    console_info(f"  Run Name : {args.run_folder_name}")
    if _md_jobs:
        _n_ctrl = sum(1 for (jn, _c) in jobs if jn.startswith("0000000"))
        _n_sel  = len(jobs) - _n_ctrl
        console_info(f"MD-ready gate: preparing {len(jobs)} of {_n_all_cifs} Best Complex CIFs "
                     f"= {_n_sel} MD_Selected + {_n_ctrl} control(s)")
    else:
        console_info(f"Found {len(jobs)} Best Complex CIFs in 2_Best_Complexes_CIFs")
    if args.quick:
        console_info(f"{ConsoleColours.OKBLUE}Quick Resume Mode: Deep validation disabled.{ConsoleColours.ENDC}")

    rank_map = load_rank_map(prod_dir)
    cat_anchor_map = load_catalytic_anchor_map(prod_dir)   # 02 dynamic-alignment anchors for the identity guard

    print("")
    console_info(f"Outputs      : {analysis_dir}")
    console_info(f"Schrödinger  : {SCHRODINGER_PATH}")
    console_info(f"PrepWizard   : {'Found' if PREPWIZARD_BIN.exists() else 'NOT FOUND - preparation will be skipped'}")
    console_info(f"Parallel Jobs: {MAX_PREP_JOBS} (Globally managed by CFG)")
    console_separator()

    # -----------------------------------------------------------------------------
    # Step 6.1: Raw PDB Generation (With Source Validation)
    # -----------------------------------------------------------------------------
    print(SEPARATOR_LIGHT, flush=True)
    console_info("Generating Raw PDBs")

    console_info("Indexing existing Raw PDBs...")
    raw_index = index_existing_files(dir_raw, CFG.SUFFIX_RAW_PDB)
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
                        '''
                        Conversion failures are otherwise silent - record them so
                        an all-fail run is surfaced rather than queuing zero
                        PrepWizard jobs and exiting 0.
                        '''
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
            # Sorted: the list is filled in thread-completion order, so WHICH 20 of them get
            # printed would otherwise change between identical runs.
            _raw_failures = sorted(_raw_failures)
            for _fjob, _freason in _raw_failures[:20]:
                console_info(f"  ! {_fjob}: {_freason}")
            if len(_raw_failures) > 20:
                console_info(f"  ! ... and {len(_raw_failures) - 20} more")

    console_info(f"Raw PDB generation stage completed in {time.time() - t0:.1f}s")
    console_separator()

    # Abort if every input CIF failed to convert - downstream stages require at
    # least one raw PDB, and a silent exit-0 here masks a total conversion failure.
    if jobs and not final_valid_raw_names:
        print("CRITICAL: all CIF->PDB conversions failed; no raw PDBs produced. "
              "Aborting before PrepWizard.", file=sys.stderr, flush=True)
        sys.exit(1)

    # -----------------------------------------------------------------------------
    # Step 6.2: Protein Preparation via PrepWizard (With Stale Check)
    # -----------------------------------------------------------------------------
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

    _prep_geom_rows: list = []
    if total_prep > 0:
        console_info(f"Queuing {total_prep} preparation tasks (batch mode)...")
        with ThreadPoolExecutor(max_workers=MAX_PREP_JOBS) as executor:
            for batch_start in range(0, total_prep, BATCH_SIZE):
                batch = prep_jobs_to_run[batch_start: batch_start + BATCH_SIZE]
                futures = {
                    executor.submit(preparation_step, job_name, dir_raw, dir_prep_clean,
                                    str(rank_map.get(job_name, "N/A")), prep_index,
                                    cat_anchor_map.get(job_name)): job_name
                    for job_name in batch
                }
                for future in as_completed(futures):
                    count += 1
                    res = future.result()
                    symbol = "✔"
                    stat_text = res["status"]
                    if stat_text != "Success": symbol = "✘"
                    try:
                        job_id_num = res["job"].split("_")[0]
                    except Exception:
                        job_id_num = "?"
                    print(f"({count}/{total_prep} | ID:{job_id_num}) {symbol} {res['job']} | [{stat_text}]", flush=True)
                    if logger: logger.info(f"{symbol} {res['job']} | [{stat_text}]")
                    '''
                    The catalytic protonation is REPORTED, never silent. A residue whose charge was
                    changed - or one left as it was - decides whether the SN2 can happen at all, and
                    the run log is where that decision has to be visible.
                    '''
                    for _role, (_rn, _num, _act) in (res.get("protonation") or {}).items():
                        _msg = f"      protonation · {_role:<5} {_rn}{_num} → {_act}"
                        print(_msg, flush=True)
                        if logger: logger.info(_msg)
                    # Residue identity guard: run after every successful PrepWizard job
                    if stat_text == "Success":
                        _prep_ok += 1
                        _prepared_pdb = dir_prep_clean / f"{res['job']}_Prepared.pdb"
                        _check_residue_identity_guard(_prepared_pdb, res["job"], CFG,
                                                      anchors=cat_anchor_map.get(res["job"]))
                        if "prep_sn2_angle" in res:
                            _prep_geom_rows.append({k: res.get(k) for k in (
                                "job", "rank", "cif_sn2_angle", "cif_dist_nuc",
                                "raw_sn2_angle", "raw_dist_nuc",
                                "prep_sn2_angle", "prep_dist_nuc", "prep_d_angle", "prep_d_dist",
                                "prep_attack_o", "prep_left_nac", "prep_geom_status")})
                    else:
                        _prep_fail += 1

    '''
    The prepared-pose geometry is written out as its own table, because it answers a question no other
    file in the run can: does the structure that MD actually starts from still hold the geometry the
    tier was granted on? The screen tiers the Boltz CIF; preparation then moves the angle by ~6° and the
    nucleophile ~0.3 Å outward, in the direction that undoes the selection. Recording both poses side by
    side makes that drift auditable instead of invisible. Nothing is gated on it.

    EVERY prepared structure on disk is measured, not merely the ones prepared on this run. Preparation
    is cached - a second run queues zero PrepWizard jobs - so measuring only the newly-prepared ones
    would leave the table empty on exactly the runs where the structures already exist, which is most of
    them. The measurement is cheap; the cache is not a reason to under-report.
    '''
    _measured = {r["job"] for r in _prep_geom_rows}
    for _pp in sorted(dir_prep_clean.glob("*_Prepared.pdb")):
        _job = _pp.name.replace("_Prepared.pdb", "")
        if _job in _measured:
            continue
        _rp = dir_raw / f"{_job}_RAW.pdb"
        if not _rp.exists():
            continue
        _nr = (cat_anchor_map.get(_job) or {}).get("Nuc")
        _cg2 = {"sn2_angle": float("nan"), "dist_nuc": float("nan")}
        try:
            _cdir = dir_raw.parent.parent / "2_Best_Complexes_CIFs"
            _chit = next(iter(sorted(_cdir.glob(f"{_job}_model_*.cif"))), None) if _cdir.is_dir() else None
            if _chit is not None:
                _cg2 = measure_sn2_geometry(_chit, _nr)
        except Exception:
            pass
        _rg, _pg2 = measure_sn2_geometry(_rp, _nr), measure_sn2_geometry(_pp, _nr)
        _da2 = _pg2["sn2_angle"] - _rg["sn2_angle"]
        _dd2 = _pg2["dist_nuc"] - _rg["dist_nuc"]
        _prep_geom_rows.append({
            "job": _job, "rank": str(rank_map.get(_job, "N/A")),
            "cif_sn2_angle": _cg2["sn2_angle"], "cif_dist_nuc": _cg2["dist_nuc"],
            "raw_sn2_angle": _rg["sn2_angle"], "raw_dist_nuc": _rg["dist_nuc"],
            "prep_sn2_angle": _pg2["sn2_angle"], "prep_dist_nuc": _pg2["dist_nuc"],
            "prep_d_angle": round(_da2, 2) if _da2 == _da2 else float("nan"),
            "prep_d_dist": round(_dd2, 2) if _dd2 == _dd2 else float("nan"),
            "prep_attack_o": _pg2["attack_o"],
            "prep_left_nac": int(_pg2["sn2_angle"] == _pg2["sn2_angle"]
                                 and (_pg2["sn2_angle"] < CFG.NAC_ANGLE_RELAXED
                                      or _pg2["dist_nuc"] > CFG.NAC_DIST_RELAXED)),
            "prep_geom_status": _pg2["status"],
        })

    if _prep_geom_rows:
        _pg = pd.DataFrame(_prep_geom_rows).sort_values(
            ["prep_sn2_angle", "job"], ascending=[False, True], kind="mergesort")
        _pg_path = dir_prep_clean.parent / "3_Comparative_Analysis" / "01_Prepared_Pose_Geometry.csv"
        _pg_path.parent.mkdir(parents=True, exist_ok=True)
        _utils_mod.atomic_write_csv(_pg, _pg_path)
        _n_left = int(pd.to_numeric(_pg["prep_left_nac"], errors="coerce").fillna(0).sum())
        _da = pd.to_numeric(_pg["prep_d_angle"], errors="coerce").dropna()
        _dd = pd.to_numeric(_pg["prep_d_dist"], errors="coerce").dropna()
        console_separator()
        console_info(f"  Prepared-pose geometry  (CIF \u2192 minimised)  \u2192  {_pg_path.name}")
        console_info(SEPARATOR_LIGHT)

        def _short_name(_j):
            _p = str(_j).split("_")
            return "_".join(_p[2:]) if len(_p) > 2 and _p[0].isdigit() else str(_j)
        _c1 = max(30, min(40, max((len(_short_name(r["job"])) for _, r in _pg.iterrows()), default=30)))
        _hn, _h22, _h11, _h4 = "\u2500" * _c1, "\u2500" * 22, "\u2500" * 11, "\u2500" * 4
        console_info("  \u250c\u2500" + _hn + "\u2500\u252c\u2500" + _h22 + "\u2500\u252c\u2500" + _h22 + "\u2500\u252c\u2500" + _h11 + "\u2500\u252c\u2500" + _h4 + "\u2510")
        _h_ang, _h_nuc = "SN2 angle CIF\u2192prep", "Nuc dist CIF\u2192prep"
        console_info(f"  \u2502 {'Complex':<{_c1}} \u2502 {_h_ang:>22} \u2502 {_h_nuc:>22} \u2502 {'Attack O':<11} \u2502 NAC \u2502")
        console_info("  \u251c\u2500" + _hn + "\u2500\u253c\u2500" + _h22 + "\u2500\u253c\u2500" + _h22 + "\u2500\u253c\u2500" + _h11 + "\u2500\u253c\u2500" + _h4 + "\u2524")

        def _row_group_key(_j):
            """Group the rows so each family prints together: identified hits first, then the control
            families (3R3U, DeHa4) as blocks; within a block, ordered by short name."""
            _s = str(_j)
            _m = re.search(r"_([A-Za-z0-9]+)_Control", _s)
            return (1 if _m else 0, _m.group(1) if _m else "", _short_name(_j))
        for _r in sorted(_pg.to_dict("records"), key=lambda _rr: _row_group_key(_rr["job"])):
            _a0, _a1 = _r.get("raw_sn2_angle"), _r.get("prep_sn2_angle")
            _d0, _d1 = _r.get("raw_dist_nuc"), _r.get("prep_dist_nuc")
            _ang = (f"{_a0:.1f}\u2192{_a1:.1f} ({_r.get('prep_d_angle'):+.1f}\u00b0)"
                    if pd.notna(_a0) and pd.notna(_a1) else "-")
            _nuc = (f"{_d0:.2f}\u2192{_d1:.2f} ({_r.get('prep_d_dist'):+.2f})"
                    if pd.notna(_d0) and pd.notna(_d1) else "-")
            _nac = "OUT" if int(_r.get("prep_left_nac", 0) or 0) else "in"
            console_info(f"  \u2502 {_short_name(_r['job']):<{_c1}} \u2502 {_ang:>22} \u2502 {_nuc:>22} \u2502 {str(_r.get('prep_attack_o','')):<11} \u2502 {_nac:>3} \u2502")
        console_info("  \u2514\u2500" + _hn + "\u2500\u2534\u2500" + _h22 + "\u2500\u2534\u2500" + _h22 + "\u2500\u2534\u2500" + _h11 + "\u2500\u2534\u2500" + _h4 + "\u2518")
        if len(_da):
            _sfx = (f"  \u00b7  [!] {_n_left} pose(s) left the relaxed NAC envelope "
                    f"(angle < {CFG.NAC_ANGLE_RELAXED:.0f}\u00b0 or nuc > {CFG.NAC_DIST_RELAXED:.1f} \u00c5) - flagged, not dropped"
                    if _n_left else "")
            console_info(f"  Drift  \u00b7  angle mean|\u0394| {_da.abs().mean():.1f}\u00b0 (max {_da.abs().max():.1f}\u00b0)  \u00b7  "
                         f"nucleophile mean {_dd.mean():+.2f} \u00c5 (max {_dd.max():+.2f}){_sfx}")
        try:
            _fig_p = plot_pose_drift(_prep_geom_rows, _pg_path.parent)
            if _fig_p:
                console_info(f"  \u2714 Comparative figure  \u2192  {_utils_mod.deflx_fig_name(_fig_p.name)}")
        except Exception as _e:                                   # noqa: BLE001
            console_info(f"  ! Pose-drift figure skipped: {type(_e).__name__}: {_e}")

        # ── All eight catalytic residues, cohort distribution CIF → minimised (raincloud) ─
        try:
            _mach_map = load_machinery_map(prod_dir)
            _cif_dir_m = dir_raw.parent.parent / "2_Best_Complexes_CIFs"
            _mach_rows = []
            for _pp in sorted(dir_prep_clean.glob("*_Prepared.pdb")):
                _job = _pp.name.replace("_Prepared.pdb", "")
                _roles = _mach_map.get(_job)
                if not _roles:
                    continue
                _chit = (next(iter(sorted(_cif_dir_m.glob(f"{_job}_model_*.cif"))), None)
                         if _cif_dir_m.is_dir() else None)
                _mach_rows.append({"job": _job,
                                   "cif_eng": measure_machinery_engagement(_chit, _roles) if _chit else {},
                                   "prep_eng": measure_machinery_engagement(_pp, _roles)})
            if _mach_rows:
                _dfig = plot_machinery_distribution(_mach_rows, _pg_path.parent)
                if _dfig:
                    console_info(f"  \u2714 Comparative figure  \u2192  {_utils_mod.deflx_fig_name(_dfig.name)}  ({len(_mach_rows)} complexes \u00d7 8 residues)")
        except Exception as _e:                                   # noqa: BLE001
            console_info(f"  ! Machinery-engagement figure skipped: {type(_e).__name__}: {_e}")

    '''
    QM ligand charges - only when asked for. The step is minutes of DFT per ligand, so it must never run
    by surprise. Step 06 consumes the resulting .mae when it builds and charges the MD system.
    '''
    if bool(getattr(CFG, "ESP_CHARGES_ENABLE", False)) or bool(globals().get("_ESP_REQUESTED", False)):
        console_separator()
        try:
            generate_esp_charges(dir_prep_clean, dir_prep_clean.parent / "4_Ligand_ESP_Charges")
        except Exception as _e:                                   # noqa: BLE001
            console_info(f"  ! ESP charges skipped: {type(_e).__name__}: {_e}")

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
        console_info(f"  [WARNING] PrepWizard failed for {_prep_fail}/{total_prep} structure(s) - see *_prepwizard_FAILED.log.")
        if _prep_ok == 0:
            console_info("  [FATAL] No structures were prepared this run - aborting (downstream MD/QM would run on missing files).")
            sys.exit(1)




def topn_extraction_phase(args):
    """Phase 2 - extract, validate and render the MD-ready complexes for handover."""

    # -----------------------------------------------------------------------------
    # Step 6.3: Path validation
    # -----------------------------------------------------------------------------
    run_dir = DEFAULT_BASE_PATH / args.run_folder_name
    if not run_dir.exists():
        # Fallback to check relative path if not in default base
        run_dir = Path.cwd() / args.run_folder_name
        if not run_dir.exists():
            console_info(f"Error: Run folder not found at: {run_dir}")
            sys.exit(1)

    prod_dir = run_dir / "1_Boltz2_Production"           # Step 02 Output
    prep_dir = run_dir / "5_TopN_and_Preparation"        # Step 05 consolidated folder
    input_data_dir = next((prod_dir / n for n in ["1_Input_Data", "1_Input_FASTA_and_SMILES"] if (prod_dir / n).exists()), prod_dir / "1_Input_Data")

    '''
    Comparison + handover outputs (Ramachandran, Controls, handover files, combined
    CSV) live under 3_Comparative_Analysis/. The raw and prepared PDBs themselves
    are NOT duplicated here - they stay in the sibling 1_Converted_Raw_PDB/ and
    2_Prepared_PDBs/ folders, where their figures are also rendered.
    '''
    final_dir = prep_dir / "3_Comparative_Analysis"
    final_dir.mkdir(parents=True, exist_ok=True)
    setup_logging_extraction(prep_dir)   # append to the single Step-05 log in the parent

    raw_pdb_dir  = prep_dir / "1_Converted_Raw_PDB"
    prep_pdb_dir = prep_dir / "2_Prepared_PDBs"

    # -----------------------------------------------------------------------------
    # Step 6.4: Load ranking logic
    # -----------------------------------------------------------------------------
    # Primary: DeFluorX Ranked CSV written by 02_Production_DeFluorX.py
    rank_csvs = (sorted(prod_dir.glob(CFG.GLOB_RANKED_CSV)) or
                 sorted(prod_dir.glob("*_Ranked_*.csv")))
    if not rank_csvs:
        console_info("Error: No Ranked CSV found in production directory.")
        sys.exit(1)

    rank_csv = _utils_mod.latest_by_mtime(rank_csvs)   # newest by mtime, not name
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

    # -----------------------------------------------------------------------------
    # Step 6.5: Load reference data (FASTA/SMILES)
    # -----------------------------------------------------------------------------
    console_info("Loading reference sequences and SMILES...")
    seq_map, smi_map, fasta_count, smi_count = load_reference_data(input_data_dir)
    console_info(f"Loaded {fasta_count} Sequences, {smi_count} SMILES.")

    # -----------------------------------------------------------------------------
    # Step 6.6: Tier / MD-ready selection
    # -----------------------------------------------------------------------------

    TIER_ORDER = CFG.TIER_ORDER

    '''
    Separate controls from candidates. Only the CFG-designated control (3R3U × fluoroacetate) is
    extracted - the other reference systems (3R3U-DFA/TFA, DeHa4-*) are never converted/prepared
    (SECTION 18 gate), so listing them here only produced spurious "Missing Raw/Prepared Control"
    warnings. The prefix test isolates control rows; _is_reference_control keeps the one real control.
    '''
    _is_control = df["job_name"].astype(str).str.startswith(CFG.CONTROL_JOB_PREFIX)
    df_controls   = df[_is_control & df["job_name"].map(_is_reference_control)].copy()
    df_candidates = df[~_is_control].copy()

    '''
    SECTION 18 gate: when the ranked CSV carries MD_Selected, extract exactly that
    cohort and skip the interactive tier prompt (non-interactive and reproducible,
    matching the CIF->PDB preparation stage).
    '''
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
    _raw_ctrl  = len([p for p in raw_pdb_dir.glob("*.pdb")  if "Control" in p.name]) if raw_pdb_dir.exists()  else 0
    _prep_ctrl = len([p for p in prep_pdb_dir.glob("*.pdb") if "Control" in p.name]) if prep_pdb_dir.exists() else 0
    console_info(f"\n  Total Raw Structures Available      : {actual_raw_count}  "
                 f"({actual_raw_count - _raw_ctrl} MD_Selected + {_raw_ctrl} control)")
    console_info(f"  Total Prepared Structures Available : {actual_prep_count}  "
                 f"({actual_prep_count - _prep_ctrl} MD_Selected + {_prep_ctrl} control)")
    console_info(f"  Total Candidate Jobs (full library) : {len(df_candidates)}")
    console_info(f"  Control Jobs (auto-extracted)       : {len(df_controls)}")
    console_separator()

    if _md_subset is not None:
        subset = _md_subset
        '''
        Route the MD-selected 3R3U reference control (Step 02 forces 3R3U × fluoroacetate into
        MD_Selected) through the SAME handover + MD machinery as the candidates, so Step 06
        simulates it as a first-class job. It carries its own Scientific_Rank, so its handover
        PDB is R{rank}_..._3R3U_Control_...; downstream flags it via is_control (looked up by
        job_name from the ranked CSV), and the run gets a real WT-FAcD·FA reference trajectory.
        '''
        if _md_col in df_controls.columns:
            _ctrl_md = df_controls[df_controls[_md_col].astype(str).str.lower().isin(["true", "1", "1.0"])]
            if len(_ctrl_md):
                subset = pd.concat([subset, _ctrl_md], ignore_index=False)
                console_info(f"  + {len(_ctrl_md)} MD-selected 3R3U reference control(s) added to the MD handover")
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
            console_info(f"  -> No selection - auto-selecting: [{available_tiers[0][0]}]")
        else:
            console_info(f"  -> User selected tier(s): {', '.join(available_tiers[i][0] for i in selected_indices)}")

        selected_tiers = [available_tiers[i][0] for i in selected_indices]
        subset = df_candidates[df_candidates[tier_col].isin(selected_tiers)].copy()
    else:
        # No tier column or no data - fall back to top-N by rank
        console_info("  Warning: No tier data found. Falling back to top-50 by rank.")
        top_n_fallback = min(args.top or 50, len(df_candidates))
        subset = df_candidates.head(top_n_fallback).copy()
        selected_tiers = ["All"]

    top_n = len(subset)
    tier_label = "_".join(selected_tiers)
    console_info(f"  Extracting {top_n} structures from tier(s): {', '.join(selected_tiers)}")
    console_separator()

    # -----------------------------------------------------------------------------
    # Step 6.7: Prepare output folders
    # -----------------------------------------------------------------------------
    folder_tag = f"{tier_label}_{top_n}hits" if top_n else "Selected"

    # Comparative Ramachandran plots (raw vs prepared) for the selected cohort.
    out_rama = final_dir / f"04_{folder_tag}_Comparative_Ramachandran_Plots"
    out_rama.mkdir(parents=True, exist_ok=True)

    # Control output folders (always created)
    out_ctrl = final_dir / "05_Controls"
    out_ctrl_raw  = out_ctrl / "1_Converted_Raw_PDB"
    out_ctrl_prep = out_ctrl / "2_Prepared_PDBs"
    out_ctrl_raw.mkdir(parents=True, exist_ok=True)
    out_ctrl_prep.mkdir(parents=True, exist_ok=True)

    # Molecular Handover Folder
    out_handover = final_dir / f"06_{folder_tag}_Molecular_Handover_Files"
    out_handover.mkdir(parents=True, exist_ok=True)

    # Initialise SDF Writer
    sdf_path = out_handover / f"{folder_tag}_Ligands.sdf"
    sdf_writer = Chem.SDWriter(str(sdf_path))
    try:

        #  -- Sub-step 6.7.1: Control-case extraction (always automatic) --
        ctrl_extracted_raw  = 0
        ctrl_extracted_prep = 0
        ctrl_csv_rows = []

        if not df_controls.empty:
            console_info(f"Extracting {len(df_controls)} control structure(s) → Controls/")
            for _, crow in df_controls.iterrows():
                ctrl_name  = str(crow.get("job_name", ""))   # the single CFG control (3R3U × fluoroacetate)
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
                _utils_mod.atomic_write_csv(pd.DataFrame(ctrl_csv_rows), ctrl_csv_path)

            console_info(f"  -> Controls: {ctrl_extracted_raw} raw, {ctrl_extracted_prep} prepared structures saved.")
        else:
            console_info("  Note: No control jobs found in the CSV (DeHa4_Control / 3R3U_Control).")
        console_separator()

        # -----------------------------------------------------------------------------
        # Step 6.8: Extraction loop with unique aggregation
        # -----------------------------------------------------------------------------
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

        for idx, (i, row) in enumerate(subset.iterrows()):
            # Metadata - handover naming keys on the true Scientific_Rank
            # (not selection order) so R_N matches every downstream MD/WaterMap artefact.
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

            '''
            --- A. Ramachandran comparison (read canonical PDBs in-place; no copies) ---
            The raw and prepared PDBs already live in 1_Converted_Raw_PDB/ and
            2_Prepared_PDBs/ - they are read directly rather than duplicated here.
            '''
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
                        rama_path = out_rama / f"Rank_{rank}_{fname_prep.replace('.pdb', '_Rama_Comparison.svg')}"
                        save_ramachandran_comparison(
                            angles_raw, angles_prep,
                            f"Rank {rank}: {p_name} (Raw)",
                            f"Rank {rank}: {p_name} (Prepared)",
                            rama_path,
                            critical_res=critical_res,
                            dpi=CFG.VIS_FIGURE_DPI, cfg=CFG
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

                '''
                Store Mol if this SMILES has not been seen yet
                Only attempt extraction if a valid source file exists
                '''
                if smi not in ligand_data and mol_source:
                    mol = extract_chain_l_mol(mol_source)
                    '''
                    The mol was parsed from PDB with sanitize=False (all-single-bond
                    topology). Restore true bond orders from the reference SMILES so the
                    handover SDF carries correct PFAS connectivity. NOTE: the SDF then
                    reflects the reference-SMILES protonation, NOT the pH-adjusted
                    (Epik) MD state - this SDF is a connectivity handover record; QM/MM
                    (Step 07) parametrises from the Schrödinger .maegz, not this SDF.
                    Best-effort: on any template-match failure keep the geometry-only
                    mol rather than dropping the ligand.
                    '''
                    if mol is not None and smi:
                        try:
                            _tmpl = Chem.MolFromSmiles(smi)
                            if _tmpl is not None:
                                mol = AllChem.AssignBondOrdersFromTemplate(_tmpl, mol)
                        except Exception as _be:
                            if logger:
                                logger.debug(f"Bond-order assignment failed for {smi}: {_be}")
                    ligand_data[smi] = (l_name, mol)

            # --- C. Extract Sequence (individual entry + unique-protein aggregation) ---
            seq = (row.get("target_sequence") or
                   row.get("sequence") or
                   row.get("Protein_Sequence") or
                   seq_map.get(p_name) or
                   "Sequence_Not_Found")
            # Full per-complex FASTA (one entry per rank - always written)
            extracted_fa.append(f">Rank_{rank}_{name} | {p_name}\n{seq}")
            # Unique-protein aggregation: key on sequence string; fall back to name when absent
            seq_key = seq if (seq and seq != "Sequence_Not_Found") else f"NAME:{p_name}"
            protein_ranks[seq_key].append(rank)
            if seq_key not in protein_data:
                protein_data[seq_key] = p_name

            # --- D. Export MD-ready PDB file to the handover folder ---
            src_prep = prep_pdb_dir / fname_prep
            if src_prep.exists():
                try:
                    title_line = f"TITLE     R{rank}_{name}\n"
                    with open(src_prep, "r") as f:
                        pdb_lines = f.readlines()
                    pdb_lines = [line for line in pdb_lines if not line.startswith("TITLE")]
                    pdb_lines.insert(0, title_line)
                    
                    dest_pdb = out_handover / f"R{rank}_{name}.pdb"
                    with open(dest_pdb, "w") as f:
                        f.writelines(pdb_lines)
                except Exception as e:
                    console_info(f"    ! Failed to copy MD-ready PDB to handover: {e}")

        # -----------------------------------------------------------------------------
        # Step 6.9: Unique ligand writing (SDF & SMILES)
        # -----------------------------------------------------------------------------

        # Iterate through unique SMILES found
        for smi, ranks in ligand_ranks.items():
            # Sort ranks numerically
            sorted_ranks = sorted(list(set(ranks)))

            '''
            Construct Composite Name: Rank-2-4-5___LigandName
            1. Join ranks with hyphens
            '''
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

        # -----------------------------------------------------------------------------
        # Step 6.10: Save data & summary
        # -----------------------------------------------------------------------------

        # Save CSV Data (hits + controls merged; controls appended with is_control flag). It lives in the
        # Molecular Handover folder - it is the scientific metadata for exactly those handover complexes.
        subset_csv_path = out_handover / f"{folder_tag}_Combined_Scientific_Data.csv"
        _ctrl_df = pd.DataFrame(ctrl_csv_rows) if ctrl_csv_rows else pd.DataFrame()
        if not _ctrl_df.empty:
            _ctrl_df = _ctrl_df.reindex(columns=subset.columns)
        '''
        Provenance flag so controls stay distinguishable from candidate hits in the merged CSV (the
        reindex above drops any control marker carried in the rows). The MD control (3R3U-FA) is both
        is_control and MD_Selected, so it was folded into `subset` for extraction AND is carried in
        _ctrl_df - drop it from the candidate side here so the merged CSV lists it once (from _ctrl_df,
        is_control=True) instead of twice.
        '''
        subset = subset[~subset["job_name"].astype(str).str.startswith(CFG.CONTROL_JOB_PREFIX)].copy()
        subset["is_control"] = False
        if not _ctrl_df.empty:
            _ctrl_df["is_control"] = True
        _combined = pd.concat([subset, _ctrl_df], ignore_index=True)
        _utils_mod.atomic_write_csv(_combined, subset_csv_path)

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

    '''
    --- Delivery Package Summary ---
    Lists everything Step 05 produced, in the order it sits on disk, so the log is a complete map of
    the deliverables (including the analysis figures and, when generated, the ESP charges).
    '''
    _esp_dir = prep_dir / "4_Ligand_ESP_Charges"
    _dp_rows = [
        ("1", "Raw PDBs + PyMOL/PLIP figures",   str(raw_pdb_dir.resolve())),
        ("2", "Prepared PDBs + PyMOL/PLIP figures", str(prep_pdb_dir.resolve())),
        ("3", "Comparative Analysis · geometry CSV + pose-drift + machinery figures", str(final_dir.resolve())),
        ("4", "Comparative Ramachandran plots",  str(out_rama.resolve())),
        ("5", "Controls (raw/prep/CSV)",         str(out_ctrl.resolve())),
        ("6", "Molecular Handover (FASTA/SDF/PDB + Combined Scientific Data CSV)", str(out_handover.resolve())),
    ]
    if _esp_dir.is_dir():
        _dp_rows.append(("7", "Ligand ESP charges (Jaguar) + summary + figure", str(_esp_dir.resolve())))
    _dp_rows.append((str(len(_dp_rows) + 1), "Step-05 log (prep + extraction)",
                     str((prep_dir / "00_TopN_and_Preparation.log").resolve())))
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

    # -----------------------------------------------------------------------------
    # Phase 2: Figure Generation
    # -----------------------------------------------------------------------------
    run_figure_generation(run_dir, final_dir)


# =============================================================================
# SECTION 7: MAIN EXECUTION
# =============================================================================


def main():
    # Collapse stacked separator rules: a caller prints a rule, a helper prints its own,
    # and the log grows triple bars with nothing between them.
    _utils_mod.install_console_rule_filter()
    _t0 = time.perf_counter()
    parser = argparse.ArgumentParser(
        description="Merged Top-N selection + CIF->PDB generation/preparation (MD-ready cohort only)")
    parser.add_argument("run_folder_name", help="Run Folder Name (e.g. Boltz-2_Run_...)")
    parser.add_argument("--quick", action="store_true", help="Quick resume: skip deep validation")
    parser.add_argument("--top", type=int, default=None, help="Fallback top-N when no MD_Selected column")
    parser.add_argument("--esp", action="store_true",
                        help="Compute QM (Jaguar ESP) partial charges for the prepared MD ligands and "
                             "write <ligand>_ESP.mae. Step 06 applies them to the built system "
                             "automatically; this step never touches the MD or WaterMap setup. Also "
                             "settable as CFG.ESP_CHARGES_ENABLE.")
    args = parser.parse_args()

    # --esp turns the QM charge step on for this run; CFG.ESP_CHARGES_ENABLE turns it on permanently.
    global _ESP_REQUESTED
    _ESP_REQUESTED = bool(args.esp)
    prep_and_convert_phase(args)
    topn_extraction_phase(args)
    _utils_mod.print_elapsed(_t0, "05_TopN_and_PDB_Preparation_DeFluorX.py")


if __name__ == "__main__":
    main()
