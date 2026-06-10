#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 06  |  Top-N Extraction & Figure Generation
===============================================================================
Phase 1 — Extraction:
  Reads the ranked CSV from Step 02, packages the best structures into a curated
  delivery folder, generates Ramachandran comparison plots, and exports FASTA,
  SMILES, and SDF files for downstream experimental validation.

Phase 2 — Figure Generation:
  Orchestrates PyMOL and PLIP to produce publication-quality figures for the
  Top-N curated complexes extracted in Phase 1.
  Uses high-concurrency task-level parallelism for maximum performance.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 June 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 06_Top-N_Extraction_FAcDs.py
  Role          : "Delivery & Photographer" — packages top-ranked structures
                  for handover, then renders publication-quality figures.
  Imports from  : 00_02_Project_Config_FAcDs.py  (CFG — threshold documentation,
                                            rendering parameters)
                  00_03_Project_Utils_FAcDs.py   (ConsoleColours, setup_logging,
                                            console_info, console_separator)
  Reads         : <Run>/1_Boltz2_Production/7_Boltz2_FAcDs_Ranked_*.csv
                  <Run>/5_PDB_Generation_Preparation/1_Converted_Raw_PDB/*.pdb
                  <Run>/5_PDB_Generation_Preparation/2_Prepared_PDBs/*.pdb
                  <Run>/1_Boltz2_Production/1_Input_FASTA_and_SMILES/*.fasta *.smi
                  <Run>/6_Top_N_Extracted/1_<Tier>_Raw_Complexes/*.pdb
                  <Run>/6_Top_N_Extracted/2_<Tier>_Prepared_Complexes/*.pdb
                  <Run>/6_Top_N_Extracted/Controls/*/*.pdb
  Writes        : <Run>/6_Top_N_Extracted/1_<Tier>_Raw_Complexes/
                  <Run>/6_Top_N_Extracted/2_<Tier>_Prepared_Complexes/
                  <Run>/6_Top_N_Extracted/3_<Tier>_Comparative_Ramachandran_Plots/
                  <Run>/6_Top_N_Extracted/4_Controls/
                  <Run>/6_Top_N_Extracted/5_<Tier>_Molecular_Handover_Files/
                  <Run>/6_Top_N_Extracted/6_Top-N_Extraction_Report.log
                  <Run>/6_Top_N_Extracted/7_<Tier>_Combined_Scientific_Data.csv
  Upstream      : 02_Production_FAcDs.py → writes ranked CSV
                  05_CIF-PDB_Preparation_FAcDs.py → writes the PDB files extracted here
  Downstream    : 08_MD_Thermodynamics_QMMM_Engine_FAcDs.py → uses extracted structures for MD
───────────────────────────────────────────────────────────────────────────────

# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Tool Dependency: Hard requirement for `pymol` and `plip` binaries; figure
     generation will be skipped if these are not available on the PATH.
  2. Parallel Rendering: Spawns multiple PyMOL instances; on some systems (e.g.,
     macOS/Wayland), this may trigger GUI focus-stealing or X11 errors if not properly configured.
  3. Interactive Timeout: The 30-second tier-selection prompt requires an
     active terminal; will auto-select the elite tier (CFG.T_PA) if no input is detected.
  4. PLIP Sensitivity: PLIP interaction detection is highly sensitive to PDB
     formatting; prepared PDBs from Step 05 are required for reliable signal.
───────────────────────────────────────────────────────────────────────────────

Usage:
    conda activate PFAS
    python 06_Top-N_Extraction_FAcDs.py Boltz-2_Run_20260309T085406Z [--top N]

── Key features ───────────────────────────────────────────────────────────────
  Phase 1 — Extraction:
  • Tier-driven extraction: interactive 30-second timeout with auto-select.
  • Ramachandran comparison: Raw vs Prepared overlaid plots with triad markers.
  • Unique-ligand aggregation: composite rank labels (Rank-1-5-10___LigandName).
  • Control extraction: separate Controls/ folder, always automatic.
  • SDF writer: 3-D ligand structures with rank and SMILES metadata.

  Phase 2 — Figure Generation:
  • PyMOL rendering: opaque + transparent interaction images, .pse session files.
  • PLIP analysis: XML, text, and PyMOL-rendered interaction reports.
  • Automated tool installation via Conda/Pip if tools are absent.
  • Task-level parallelism across all complexes for maximum throughput.

  Scientific references (PLIP/interaction profiling):
  • Salentin et al. (2015) Nucleic Acids Res 43:W443–W447.
  • McGaughey et al. (1998) J Biol Chem 273:15458–15463.
  • Gallivan & Dougherty (1999) PNAS 96:9459–9464.
  • Wilcken et al. (2013) J Med Chem 56:1363–1388.
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
import shutil
import argparse
import logging
import re
import subprocess
import multiprocessing
import time
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

# -------------------------------------------------------------------------------
# Step 1.2: Scientific Stack Imports
# -------------------------------------------------------------------------------
import math
import gemmi
import tempfile
import xml.etree.ElementTree as _ET
import numpy as np
import matplotlib
matplotlib.use('Agg')  # non-interactive backend; must be set before pyplot import
import matplotlib.pyplot as plt
import matplotlib.patches as _mpatches
from matplotlib.patches import FancyBboxPatch as _FancyBboxPatch
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from matplotlib.lines import Line2D as _Line2D

import pandas as pd
from typing import List, Tuple, Dict, Any
from Bio import SeqIO
from Bio.PDB import PDBParser as _PDBParser
from rdkit import Chem

# -------------------------------------------------------------------------------
# Step 1.3: Pipeline modules (00_02 (config & utils)) via importlib
# -------------------------------------------------------------------------------
import importlib.util as _ilu

def _load_module(name: str, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod  = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_REPO_DIR  = Path(__file__).resolve().parent
_cfg_mod   = _load_module("ProjectConfig", _REPO_DIR / "00_02_Project_Config_FAcDs.py")
_utils_mod = _load_module("ProjectUtils",  _REPO_DIR / "00_03_Project_Utils_FAcDs.py")

CFG            = _cfg_mod.CFG()
ConsoleColours     = _utils_mod.ConsoleColours
SEPARATOR_HEAVY    = _utils_mod.SEPARATOR_HEAVY
SEPARATOR_LIGHT    = _utils_mod.SEPARATOR_LIGHT

compute_ramachandran_angles = _utils_mod.compute_ramachandran_angles
_rama_stats = _utils_mod._rama_stats
save_ramachandran_comparison = _utils_mod.save_ramachandran_comparison
safe_name = _utils_mod.safe_name

# -------------------------------------------------------------------------------
# Step 1.4: Global Configuration
# -------------------------------------------------------------------------------
DEFAULT_BASE_PATH = Path.cwd()
SEPARATOR = "-" * 80
logger = None


# ===============================================================================
# SECTION 2: LOGGING & UTILITIES
# ===============================================================================

def setup_logging(output_dir: Path) -> Path:
    """Initialises the extraction log via the shared utility."""
    global logger
    log_file = output_dir / "6_Top-N_Extraction_Report.log"
    logger = _utils_mod.setup_logging(log_file, logger_name="extraction", mode="w", timestamp=False)
    return log_file

def console_info(msg: str) -> None:
    _utils_mod.console_info(msg, logger)

def console_separator() -> None:
    _utils_mod.console_separator(logger, heavy=True)




# ===============================================================================
# SECTION 3: RAMACHANDRAN PLOTTING ENGINE
# ===============================================================================

# --- Ramachandran helpers (relocated to 00_03_Project_Utils_FAcDs.py) ---------


# ===============================================================================
# SECTION 4: CORE LOGIC & HELPERS
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 4.1: Data Loading
# -------------------------------------------------------------------------------
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
                    seq_map[clean_id] = str(r.seq)
                    # Store original ID as fallback
                    seq_map[r.id] = str(r.seq)
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
        with open(pdb_path, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    if len(line) > 21 and line[21] == 'L':
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

def main():
    # -------------------------------------------------------------------------------
    # Step 5.1: Argument Parsing
    # -------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description="Boltz-2 Analysis & Delivery Engine")
    
    # Argument 1: Run Folder Name (Required)
    parser.add_argument("run_folder_name", help="Run Folder Name (e.g. Boltz-2_Run_...)")
    
    # Argument 2: Top N (Optional via CLI)
    parser.add_argument("--top", type=int, default=None, help="Number of top compounds to extract (Optional)")
    
    args = parser.parse_args()

    # -------------------------------------------------------------------------------
    # Step 5.2: Path Validation
    # -------------------------------------------------------------------------------
    run_dir = DEFAULT_BASE_PATH / args.run_folder_name
    if not run_dir.exists():
        # Fallback to check relative path if not in default base
        run_dir = Path.cwd() / args.run_folder_name
        if not run_dir.exists():
            console_info(f"Error: Run folder not found at: {run_dir}")
            sys.exit(1)

    prod_dir = run_dir / "1_Boltz2_Production"           # Step 02 Output
    prep_dir = run_dir / "5_PDB_Generation_Preparation"  # Step 05 Output
    input_data_dir = next((prod_dir / n for n in ["1_Input_FASTA_and_SMILES", "1_Input_Data"] if (prod_dir / n).exists()), prod_dir / "1_Input_FASTA_and_SMILES")

    final_dir = run_dir / "6_Top_N_Extracted"
    final_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(final_dir)

    raw_pdb_dir  = prep_dir / "1_Converted_Raw_PDB"
    prep_pdb_dir = prep_dir / "2_Prepared_PDBs"

    # -------------------------------------------------------------------------------
    # Step 5.3: Header Output
    # -------------------------------------------------------------------------------
    _utils_mod.print_script_banner(
        "06_Top-N_Extraction_FAcDs.py",
        "Top-N Delivery  ·  Tier Extraction  ·  Ramachandran Validation  ·  PyMOL/PLIP Figure Generation",
    )
    console_info(f"  Run Name : {args.run_folder_name}")
    console_info(f"  Output   : {final_dir}")
    console_separator()

    # -------------------------------------------------------------------------------
    # Step 5.4: Load Ranking Logic
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
    if 'status' in df.columns:
        df = df[df['status'] == 'Success']

    # Sort by Rank
    if 'Scientific_Rank' in df.columns:
        df = df.sort_values('Scientific_Rank')
    elif 'Production_Rank' in df.columns:
        df = df.sort_values('Production_Rank')
    else:
        console_info("Warning: Rank column missing. Using default sort.")

    # -------------------------------------------------------------------------------
    # Step 5.5: Load Reference Data (FASTA/SMILES)
    # -------------------------------------------------------------------------------
    console_info("Loading reference sequences and SMILES...")
    seq_map, smi_map, fasta_count, smi_count = load_reference_data(input_data_dir)
    console_info(f"Loaded {fasta_count} Sequences, {smi_count} SMILES.")

    # -------------------------------------------------------------------------------
    # Step 5.6: Interactive Tier Selection (30-Second Timeout)
    # -------------------------------------------------------------------------------
    import select as _select

    TIER_ORDER = CFG.TIER_ORDER

    # Separate controls from candidates (controls always extracted independently)
    def _is_control(row):
        jn = str(row.get('job_name', ''))
        return jn.startswith('0000000')

    df_controls  = df[df.apply(_is_control, axis=1)].copy()
    df_candidates = df[~df.apply(_is_control, axis=1)].copy()

    # Build tier availability table from candidates only
    tier_col = next((c for c in ['degrader_tier', 'Degrader_Tier', 'tier'] if c in df_candidates.columns), None)
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

    if available_tiers:
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
                for part in raw_input.split(','):
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
    # Step 5.7: Prepare Output Folders
    # -------------------------------------------------------------------------------
    folder_tag = f"{tier_label}_{top_n}hits" if top_n else "Selected"
    out_raw  = final_dir / f"1_{folder_tag}_Raw_Complexes"
    out_prep = final_dir / f"2_{folder_tag}_Prepared_Complexes"
    out_raw.mkdir(parents=True, exist_ok=True)
    out_prep.mkdir(parents=True, exist_ok=True)
    out_rama = final_dir / f"3_{folder_tag}_Comparative_Ramachandran_Plots"
    out_rama.mkdir(parents=True, exist_ok=True)

    # Control output folders (always created)
    out_ctrl = final_dir / "4_Controls"
    out_ctrl_raw  = out_ctrl / "1_Raw"
    out_ctrl_prep = out_ctrl / "2_Prepared"
    out_ctrl_raw.mkdir(parents=True, exist_ok=True)
    out_ctrl_prep.mkdir(parents=True, exist_ok=True)

    # Molecular Handover Folder
    out_handover = final_dir / f"5_{folder_tag}_Molecular_Handover_Files"
    out_handover.mkdir(parents=True, exist_ok=True)

    # Initialise SDF Writer
    sdf_path = out_handover / f"{folder_tag}_Ligands.sdf"
    sdf_writer = Chem.SDWriter(str(sdf_path))
    count_sdf = 0
    try:

        # -------------------------------------------------------------------------------
        # Step 5.7.1: Control Case Extraction (Always Automatic)
        # -------------------------------------------------------------------------------
        ctrl_extracted_raw  = 0
        ctrl_extracted_prep = 0
        ctrl_csv_rows = []

        if not df_controls.empty:
            console_info(f"Extracting {len(df_controls)} control structure(s) → Controls/")
            for _, crow in df_controls.iterrows():
                ctrl_name  = str(crow.get('job_name', ''))
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
        # Step 5.8: Extraction Loop with Unique Aggregation
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
            rank = row.get('Scientific_Rank', i+1)
            name = row['job_name']
        
            # Determine Protein and Ligand names for Lookup
            p_name = str(row.get('Protein_Name', row.get('protein', 'Unknown')))
            l_name = str(row.get('Ligand_Name', row.get('ligand', 'Unknown')))

            # File names
            fname_raw = f"{name}_RAW.pdb"
            fname_prep = f"{name}_Prepared.pdb"
        
            # Parse critical residues from Active_Site_Triad_Map
            triad_map = str(row.get('Active_Site_Triad_Map', ''))
            critical_res = {}
            if triad_map and triad_map != 'nan':
                for item in triad_map.split('|'):
                    item = item.strip()
                    if ':' in item:
                        role, res = item.split(':', 1)
                        res = res.strip()
                        resname = res[:3]
                        _m = re.search(r'\d+', res)
                        resnum_str = _m.group() if _m else ''
                        if resnum_str:
                            critical_res[int(resnum_str)] = (resname, role.strip())

            # --- A. Extract Structures (Always extract the PDBs individually) ---
            src_raw = raw_pdb_dir / fname_raw
            if src_raw.exists():
                shutil.copy2(src_raw, out_raw / f"Rank_{rank}_{fname_raw}")
                extracted_raw_count += 1
        
            src_prep = prep_pdb_dir / fname_prep
            if src_prep.exists():
                shutil.copy2(src_prep, out_prep / f"Rank_{rank}_{fname_prep}")
                extracted_prep_count += 1
                try:
                    angles_raw = []
                    angles_prep = []
                    if src_raw.exists():
                        extracted_raw_path = out_raw / f"Rank_{rank}_{fname_raw}"
                        st_raw = gemmi.read_structure(str(extracted_raw_path))
                        angles_raw = compute_ramachandran_angles(st_raw)
                    if src_prep.exists():
                        extracted_prep_path = out_prep / f"Rank_{rank}_{fname_prep}"
                        st_prep = gemmi.read_structure(str(extracted_prep_path))
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

            # --- B. Ligand Data Collection (Prefer RAW as requested) ---
            # Prioritize Raw for the structure data if available
            mol_source = src_raw if src_raw.exists() else (src_prep if src_prep.exists() else None)
        
            smi = row.get('smiles', row.get('Ligand_Smile', smi_map.get(l_name, "SMILES_Not_Found")))
        
            if smi and smi != "SMILES_Not_Found":
                # Store Rank
                ligand_ranks[smi].append(rank)
            
                # Store Mol if this SMILES has not been seen yet
                # Only attempt extraction if a valid source file exists
                if smi not in ligand_data and mol_source:
                    mol = extract_chain_l_mol(mol_source)
                    ligand_data[smi] = (l_name, mol) 
        
            # --- C. Extract Sequence (individual entry + unique-protein aggregation) ---
            seq = (row.get('target_sequence') or
                   row.get('sequence') or
                   row.get('Protein_Sequence') or
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
        # Step 5.9: Unique Ligand Writing (SDF & SMILES)
        # -------------------------------------------------------------------------------
    
        # Iterate through unique SMILES found
        for smi, ranks in ligand_ranks.items():
            # Sort ranks numerically
            sorted_ranks = sorted(list(set(ranks)))
        
            # Construct Composite Name: Rank-2-4-5___LigandName
            # 1. Join ranks with hyphens
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
        # Step 5.10: Save Data & Summary
        # -------------------------------------------------------------------------------

        # Save CSV Data (hits + controls merged; controls appended with is_control flag)
        subset_csv_path = final_dir / f"7_{folder_tag}_Combined_Scientific_Data.csv"
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
        ("1", "Raw PDBs (selected)",           str(out_raw.resolve())),
        ("2", "Prepared PDBs (selected)",      str(out_prep.resolve())),
        ("3", "Comparative Ramachandran Plots", str(out_rama.resolve())),
        ("4", "Controls (Raw/Prep/CSV)",       str(out_ctrl.resolve())),
        ("5", "Molecular Handover (FASTA/SDF)", str(out_handover.resolve())),
        ("6", "Top-N Extraction Report (Log)", str((final_dir / "6_Top-N_Extraction_Report.log").resolve())),
        ("7", "Combined Scientific Data (CSV)", str(subset_csv_path.resolve())),
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
# Step 6.1: Figure-generation constants (sourced from CFG)
# -------------------------------------------------------------------------------
# These are resolved lazily (inside run_figure_generation) so that CFG is
# available; they are module-level names for clarity.
def _fig_constants():
    """Return a dict of figure-generation constants from CFG."""
    return {
        "IMG_WIDTH":      CFG.VIS_IMG_WIDTH,
        "IMG_HEIGHT":     CFG.VIS_IMG_HEIGHT,
        "RAY_TRACE":      CFG.VIS_RAY_TRACE,
        "DIST_HBOND":     CFG.THRESHOLD_HB_DIST_MAX,
        "DIST_F_CONTACT": CFG.VIS_F_CONTACT_RADIUS,
        "DIST_POCKET":    CFG.VIS_POCKET_RADIUS,
        "TIMEOUT_PYMOL":  CFG.VIS_TIMEOUT_PYMOL,
        "TIMEOUT_PLIP":   CFG.VIS_TIMEOUT_PLIP,
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

def get_meta(base_name: str) -> dict:
    """Returns {'p': protein, 'l': ligand} for a complex base_name."""
    for jn, val in METADATA_CACHE.items():
        if jn in base_name:
            return val
    return {"p": base_name, "l": "PFAS"}

# -------------------------------------------------------------------------------
# Step 6.2: Figure-generation logger (separate from extraction logger)
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
# Step 6.3: Subprocess runner
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
    # Cap internal threading in numpy/OpenBabel/MKL so each subprocess uses
    # exactly 1 thread; concurrency is controlled at the pool level instead.
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
            try: proc.kill()
            except Exception: pass
        if f_handle:
            f_handle.write(f"\n[ERROR]: {e}\n")
            f_handle.close()
        return False

# -------------------------------------------------------------------------------
# Step 6.4: Visualisation engines
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

    # Precise ligand selection: chain L residue number is authoritative.
    # Using AND (not OR) prevents UNK/non-standard protein residues also placed
    # in chain L by Schrödinger prep from being selected as ligand atoms.
    lig_sel = f"(chain L and resi {lig_num})"

    # Cap PyMOL ray-tracing threads: 2 PyMOL workers run concurrently, so each
    # should use at most (cpu-2)//2 ray threads to stay within the cpu-2 budget.
    _max_rt = max(1, (multiprocessing.cpu_count() - 2) // 2)

    script = [
        "reinitialize",
        f'load "{pdb_path.resolve()}", complex',
        "remove solvent", "remove hydrogens",
        # Lighting — two_sided illuminates inside faces of the binding pocket,
        # preventing the cavity from rendering as a solid black void.
        "bg_color gray60",
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
        "set transparency, 0.50",
        # Protein cartoon underneath surface
        "show cartoon, polymer.protein",
        "color gray70, polymer.protein",
        "set cartoon_transparency, 0.35",
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
        # Salt bridges — cyan dashes
        "dist salt_bridge, (ligand and (name O*,N*)), (interacting and (name N*,O*)), mode=2, cutoff=4.5",
        "color cyan, salt_bridge",
        # Residue labels — white, bold, readable on grey background
        "set label_size, 9", "set label_font_id, 7",
        "set label_color, white",
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
            tail = log_path.read_text(errors='replace')[-600:]
            _fig_log(f"  [!] PyMOL log tail for {base_name}:\n{tail}")
    return "Success" if success else "Failed"

def _plip_coo(el):
    """Parse a PLIP <ligcoo> or <protcoo> element → numpy array [x, y, z]."""
    return np.array([float(el.find('x').text),
                     float(el.find('y').text),
                     float(el.find('z').text)])


def _parse_plip_xml(xml_path, lig, pro):
    """
    Parse a PLIP XML report into a contacts list compatible with _im_project /
    _im_render_diagram.

    Each contact dict carries:
      key, resname, resnum, chain, dist, itype, is_hbond, is_salt,
      lig_atom (nearest ligand atom), prot_atom (''), center (3-D protcoo)
    """
    _ITYPE_ORDER = ['hbond', 'halogen', 'salt', 'water', 'pistack', 'pication',
                    'hydrophobic', 'contact']

    tree = _ET.parse(str(xml_path))
    root = tree.getroot()

    lig_pos = np.array([a['pos'] for a in lig])

    contacts_by_key = {}   # key → best contact dict (deduplicate by residue)

    def _closest_lig_atom(ligcoo):
        dists = np.linalg.norm(lig_pos - ligcoo, axis=1)
        idx   = int(dists.argmin())
        return lig[idx], float(dists[idx])

    def _register(resnr, restype, reschain, dist, itype, protcoo_3d, ligcoo_3d):
        key = (reschain, int(resnr), restype)
        la, _  = _closest_lig_atom(ligcoo_3d)
        # C–F bonds are leaving groups in FAcDs SN2; exclude F from H-bond / halogen contacts
        if la['elem'] == 'F' and itype in ('hbond', 'halogen'):
            return
        is_hbond  = itype == 'hbond'
        is_salt   = itype == 'salt'
        entry = {'key': key, 'resname': restype, 'resnum': int(resnr),
                 'chain': reschain, 'dist': dist, 'itype': itype,
                 'is_hbond': is_hbond, 'is_salt': is_salt,
                 'lig_atom': la, 'prot_atom': '', 'center': protcoo_3d}
        # Keep highest-priority interaction type per residue
        if key not in contacts_by_key:
            contacts_by_key[key] = entry
        else:
            existing = contacts_by_key[key]
            if _ITYPE_ORDER.index(itype) < _ITYPE_ORDER.index(existing['itype']):
                contacts_by_key[key] = entry

    for bs in root.findall('bindingsite'):
        iact = bs.find('interactions')
        if iact is None:
            continue

        for hb in iact.findall('./hydrogen_bonds/hydrogen_bond'):
            protcoo = _plip_coo(hb.find('protcoo'))
            ligcoo  = _plip_coo(hb.find('ligcoo'))
            dist    = float(hb.find('dist_d-a').text)
            _register(hb.find('resnr').text, hb.find('restype').text,
                      hb.find('reschain').text, dist, 'hbond', protcoo, ligcoo)

        for hx in iact.findall('./halogen_bonds/halogen_bond'):
            protcoo = _plip_coo(hx.find('protcoo'))
            ligcoo  = _plip_coo(hx.find('ligcoo'))
            dist    = float(hx.find('dist').text)
            _register(hx.find('resnr').text, hx.find('restype').text,
                      hx.find('reschain').text, dist, 'halogen', protcoo, ligcoo)

        for sb in iact.findall('./salt_bridges/salt_bridge'):
            protcoo = _plip_coo(sb.find('protcoo'))
            ligcoo  = _plip_coo(sb.find('ligcoo'))
            dist    = float(sb.find('dist').text)
            _register(sb.find('resnr').text, sb.find('restype').text,
                      sb.find('reschain').text, dist, 'salt', protcoo, ligcoo)

        for wb in iact.findall('./water_bridges/water_bridge'):
            protcoo = _plip_coo(wb.find('protcoo'))
            ligcoo  = _plip_coo(wb.find('ligcoo'))
            dist_d  = float(wb.find('dist_d-w').text)
            dist_a  = float(wb.find('dist_a-w').text)
            _register(wb.find('resnr').text, wb.find('restype').text,
                      wb.find('reschain').text, max(dist_d, dist_a), 'water',
                      protcoo, ligcoo)

        for ps in iact.findall('./pi_stacks/pi_stack'):
            protcoo = _plip_coo(ps.find('protcoo'))
            ligcoo  = _plip_coo(ps.find('ligcoo'))
            dist    = float(ps.find('dist').text)
            _register(ps.find('resnr').text, ps.find('restype').text,
                      ps.find('reschain').text, dist, 'pistack', protcoo, ligcoo)

        for pc in iact.findall('./pi_cation_interactions/pi_cation_interaction'):
            protcoo = _plip_coo(pc.find('protcoo'))
            ligcoo  = _plip_coo(pc.find('ligcoo'))
            dist    = float(pc.find('dist').text)
            _register(pc.find('resnr').text, pc.find('restype').text,
                      pc.find('reschain').text, dist, 'pication', protcoo, ligcoo)

        for hp in iact.findall('./hydrophobic_interactions/hydrophobic_interaction'):
            protcoo = _plip_coo(hp.find('protcoo'))
            ligcoo  = _plip_coo(hp.find('ligcoo'))
            dist    = float(hp.find('dist').text)
            _register(hp.find('resnr').text, hp.find('restype').text,
                      hp.find('reschain').text, dist, 'hydrophobic', protcoo, ligcoo)

        # bs_residues with contact="True" and no typed interaction → "contact"
        for bsr in bs.findall('./bs_residues/bs_residue'):
            if bsr.get('contact', 'False') != 'True':
                continue
            text   = bsr.text.strip()           # e.g. "150A"
            chain  = text[-1]
            resnr  = text[:-1]
            restype = bsr.get('aa', 'UNK')
            key    = (chain, int(resnr), restype)
            if key not in contacts_by_key:
                # No typed interaction — approximate position from PDB
                pkey = next((k for k in pro if k[0] == chain and k[1] == int(resnr)), None)
                if pkey:
                    protcoo = pro[pkey]['center']
                    la, dist_la = _closest_lig_atom(protcoo)
                    contacts_by_key[key] = {
                        'key': key, 'resname': restype, 'resnum': int(resnr),
                        'chain': chain, 'dist': float(bsr.get('min_dist', 5.0)),
                        'itype': 'contact', 'is_hbond': False, 'is_salt': False,
                        'lig_atom': la, 'prot_atom': '', 'center': protcoo,
                    }

    return sorted(contacts_by_key.values(), key=lambda x: x['dist'])


def _im_render_diagram(lig_2d, lig, res_2d, contacts, out_png, mode='distance'):
    """
    Shared matplotlib renderer for both InteractionMap (distance-based) and
    PLIP (XML-based) diagrams.  mode='distance' or 'plip'.
    """
    # Interaction type → (linewidth, linestyle, colour, show_dist_label)
    _ITYPE_STYLE = {
        'hbond':       (2.2, (0, (6, 3)),  '#E67E22', True),
        'halogen':     (2.0, (0, (4, 2)),  '#1D8348', True),
        'salt':        (2.0, (0, (3, 2)),  '#C0392B', True),
        'water':       (1.5, (0, (2, 2)),  '#5DADE2', True),
        'pistack':     (1.8, (0, (5, 2)),  '#2471A3', False),
        'pication':    (1.8, (0, (4, 2)),  '#7D3C98', False),
        'hydrophobic': (0.9, (0, (2, 4)),  '#BDC3C7', False),
        'contact':     (0.9, (0, (2, 4)),  '#BDC3C7', False),
    }
    _DIST_COL = {'hbond': ('#884400','#FDEBD0','#E67E22'),
                 'halogen': ('#0A3D0A','#D5F5E3','#1D8348'),
                 'salt': ('#7B241C','#FADBD8','#C0392B'),
                 'water': ('#1A5276','#D6EAF8','#5DADE2')}

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

    fig = Figure(figsize=(12, 12), facecolor='white')
    canvas = FigureCanvasAgg(fig)
    ax  = fig.add_subplot(111, aspect='equal')
    ax.axis('off')
    ax.set_xlim(_xlo, _xhi)
    ax.set_ylim(_ylo, _yhi)

    # Binding pocket background
    ax.add_patch(_mpatches.Circle((0, 0), 4.2, color='#EAF2FF', zorder=0, alpha=0.6))
    ax.add_patch(_mpatches.Circle((0, 0), 4.2, color='#AED6F1', fill=False,
                            linewidth=1.2, linestyle='--', zorder=0, alpha=0.4))
    ax.text(0, -3.7, 'Binding Pocket', ha='center', fontsize=7.5,
            color='#85929E', style='italic', zorder=1)

    # Interaction lines
    name2idx = {a['name']: i for i, a in enumerate(lig)}
    for c in contacts:
        rpos = res_2d[c['key']]
        li   = name2idx.get(c['lig_atom']['name'], 0) if c.get('lig_atom') else 0
        lap  = lig_2d[li]
        itype = c.get('itype', 'hbond' if c.get('is_hbond') else
                      'salt' if c.get('is_salt') else 'contact')
        lw, ls, col, show_dist = _ITYPE_STYLE.get(itype, _ITYPE_STYLE['contact'])
        ax.plot([lap[0], rpos[0]], [lap[1], rpos[1]],
                lw=lw, ls=ls, color=col, alpha=0.9, zorder=2,
                solid_capstyle='round')
        if show_dist:
            mx, my = (lap[0]+rpos[0])/2, (lap[1]+rpos[1])/2
            tc, fc, ec = _DIST_COL.get(itype, ('#555','#EEE','#999'))
            ax.text(mx, my, f"{c['dist']:.1f} Å",
                    fontsize=8.5, ha='center', va='center', fontweight='bold',
                    color=tc,
                    bbox=dict(fc=fc, ec=ec, alpha=0.88,
                              boxstyle='round,pad=0.22', linewidth=0.8),
                    zorder=8)

    # Ligand bonds
    for i, a in enumerate(lig):
        for j, b in enumerate(lig):
            if j <= i: continue
            if np.linalg.norm(a['pos'] - b['pos']) < 1.85:
                p1, p2 = lig_2d[i], lig_2d[j]
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]],
                        color='#2C3E50', lw=2.8, solid_capstyle='round',
                        zorder=4, alpha=0.85)

    # Ligand atoms
    for i, a in enumerate(lig):
        xy  = lig_2d[i]
        col = _IM_ELEM_COLORS.get(a['elem'], _IM_ELEM_COLORS['other'])
        r   = {'C':0.17,'N':0.19,'O':0.19,'F':0.17,'S':0.22}.get(a['elem'], 0.15)
        ax.add_patch(_mpatches.Circle(xy, r, color=col, zorder=5, ec='white', lw=1.4))
        if a['elem'] not in ('C',):
            ax.text(xy[0], xy[1], a['elem'],
                    ha='center', va='center', fontsize=6.5,
                    color='white', fontweight='bold', zorder=6)

    # Residue boxes
    bw, bh = 1.5, 0.68
    for c in contacts:
        rpos  = res_2d[c['key']]
        col   = _IM_RES_COLORS.get(c['resname'], '#717D7E')
        label = f"{c['resname']} {c['resnum']}"
        itype = c.get('itype', 'contact')
        _BADGE = {
            'hbond': 'H-bond', 'halogen': 'halogen bond', 'salt': 'salt bridge',
            'water': 'water bridge', 'pistack': 'π-stack', 'pication': 'π-cation',
            'hydrophobic': 'hydrophobic', 'contact': 'contact',
        }
        badge = _BADGE.get(itype, 'contact')
        ax.add_patch(_FancyBboxPatch(
            (rpos[0]-bw/2+0.04, rpos[1]-bh/2-0.04), bw, bh,
            boxstyle='round,pad=0.1', facecolor='#C0C0C0',
            alpha=0.22, zorder=5, linewidth=0))
        ax.add_patch(_FancyBboxPatch(
            (rpos[0]-bw/2, rpos[1]-bh/2), bw, bh,
            boxstyle='round,pad=0.1', facecolor=col,
            edgecolor='white', linewidth=1.6, alpha=0.95, zorder=6))
        ax.text(rpos[0], rpos[1]+0.10, label,
                ha='center', va='center', fontsize=9.5,
                fontweight='bold', color='white', zorder=7)
        ax.text(rpos[0], rpos[1]-0.18, badge,
                ha='center', va='center', fontsize=7,
                color='white', alpha=0.9, zorder=7)

    # Unified legend (bottom, all entries in one block)
    _leg = [
        _mpatches.Patch(color='#E67E22', label='H-bond'),
        _mpatches.Patch(color='#1D8348', label='Halogen bond'),
        _mpatches.Patch(color='#C0392B', label='Salt bridge'),
        _mpatches.Patch(color='#5DADE2', label='Water bridge'),
        _mpatches.Patch(color='#2471A3', label='π-stack'),
        _mpatches.Patch(color='#7D3C98', label='π-cation'),
        _mpatches.Patch(color='#BDC3C7', label='Contact'),
        _Line2D([],[],color='none', label=''),
        _mpatches.Patch(color='#C0392B', label='ASP/GLU'),
        _mpatches.Patch(color='#2471A3', label='ARG/LYS'),
        _mpatches.Patch(color='#1E8449', label='HIS'),
        _mpatches.Patch(color='#7D3C98', label='TRP/PHE'),
        _mpatches.Patch(color='#BA4A00', label='TYR'),
        _mpatches.Patch(color='#117A65', label='SER/THR/ASN/GLN'),
        _mpatches.Patch(color='#626567', label='Hydrophobic'),
        _Line2D([],[],color='none', label=''),
    ]
    # Append atom entries inline so everything sits in one box
    for elem, ec in [('C','#2C3E50'),('N','#1A5276'),('O','#A93226'),('F','#1D8348')]:
        _leg.append(_mpatches.Patch(color=ec, label=f'Lig {elem}'))

    # Dynamic legend: anchor just below the lowest content point (circle bottom
    # or lowest residue box, whichever is further down), with a small gap.
    # Uses data-space coordinates so placement adapts to content density.
    _y_res_bottom = (min(pos[1] - bh / 2 for pos in res_2d.values())
                     if res_2d else -4.2)
    _y_legend_top = min(_y_res_bottom, -4.2) - 0.35
    _cx_data      = (_xlo + _xhi) / 2
    ax.legend(handles=_leg, loc='upper center', ncol=4,
              bbox_to_anchor=(_cx_data, _y_legend_top),
              bbox_transform=ax.transData,
              fontsize=8.5, frameon=True, framealpha=0.95,
              edgecolor='#CCCCCC',
              title=f'{"PLIP interactions" if mode=="plip" else "Interactions"}  |  Residue type  |  Ligand atoms',
              title_fontsize=8.5, columnspacing=0.8, handlelength=1.1,
              borderpad=0.6)

    try:
        canvas.print_figure(str(out_png), dpi=CFG.VIS_FIGURE_DPI, bbox_inches='tight',
                            pad_inches=0.04, facecolor='white', edgecolor='none')
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

        lig, pro = _im_parse_pdb(pdb_path)
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
            _im_render_diagram(lig_2d, lig, res_2d, contacts, out_png, mode='plip')
            rendered = True

        return "Success" if rendered else "Failed"
    except Exception as _e:
        _fig_log(f"  [!] PLIP render failed for {base_name}: {type(_e).__name__}: {_e}")
        return "Failed"
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)


# -------------------------------------------------------------------------------
# Step 6.45: Matplotlib interaction diagram (pure Python — no PyMOL required)
# -------------------------------------------------------------------------------

_IM_CONTACT_DIST  = CFG.CATALYTIC_DIST_CUTOFF
_IM_HBOND_DIST    = CFG.THRESHOLD_HB_DIST_MAX
_IM_HBOND_DONORS  = {'N','NH1','NH2','NE','NE2','ND1','ND2','NZ','OG','OG1','OH','NE1'}
_IM_HBOND_ACC     = {'O','OD1','OD2','OE1','OE2','OG','OG1','OH','NE2','ND1',
                     'O1','O2','O3','O4','O5','O6'}
_IM_RES_COLORS = {
    'ASP':'#C0392B','GLU':'#C0392B',
    'ARG':'#2471A3','LYS':'#2471A3',
    'HIS':'#1E8449','TRP':'#7D3C98','TYR':'#BA4A00',
    'SER':'#117A65','THR':'#117A65','ASN':'#117A65','GLN':'#117A65',
    'PHE':'#6C3483',
    'LEU':'#626567','ILE':'#626567','VAL':'#626567',
    'ALA':'#626567','GLY':'#626567','PRO':'#626567',
    'MET':'#7E5109','CYS':'#7E5109',
}
_IM_ELEM_COLORS = {
    'C':'#2C3E50','N':'#1A5276','O':'#A93226',
    'F':'#1D8348','S':'#D4AC0D','other':'#717D7E',
}
_IM_HYDROPHOBIC = {'LEU','ILE','VAL','PHE','TRP','PRO','MET','ALA','GLY','CYS'}
_IM_AA3 = {
    'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL',
    'UNK','MSE','SEC','PYL','HIE','HID','HIP','ASH','GLH','LYN',
    'CYX','CYM',
}


def _im_parse_pdb(path):
    """Return (lig_atoms, pro_residues) from a prepared PDB file.

    Schrödinger prep places non-standard backbone residues as HETATM in chain L.
    Guard: exclude AA/UNK residues from ligand atoms even when chain == 'L'.
    """
    parser = _PDBParser(QUIET=True)
    struct = parser.get_structure("X", str(path))
    model  = struct[0]
    lig, pro = [], {}
    for chain in model:
        for res in chain:
            hetflag = res.id[0].strip()
            resname = res.resname.strip()
            heavy   = [a for a in res if a.element and a.element.strip() not in ('H', '')]
            is_lig_chain = chain.id == 'L' and resname not in _IM_AA3
            is_hetatm    = bool(hetflag) and hetflag != 'W' and resname not in _IM_AA3
            if is_lig_chain or is_hetatm:
                for a in heavy:
                    lig.append({'name': a.name.strip(),
                                'elem': a.element.strip().upper(),
                                'pos':  a.coord.copy()})
            elif not hetflag or resname in _IM_AA3:
                if not heavy: continue
                key = (chain.id, res.id[1], resname)
                atom_dict = {a.name.strip(): a.coord.copy() for a in heavy}
                pro[key] = {'resname': resname, 'resnum': res.id[1],
                            'chain': chain.id, 'atoms': atom_dict,
                            'center': np.mean([a.coord for a in heavy], axis=0)}
    return lig, pro


def _im_find_contacts(lig, pro):
    lig_pos = np.array([a['pos'] for a in lig])
    contacts = []
    for key, res in pro.items():
        prot_pos = np.array(list(res['atoms'].values()))
        prot_nms = list(res['atoms'].keys())
        diffs = lig_pos[:, None, :] - prot_pos[None, :, :]
        dists = np.linalg.norm(diffs, axis=-1)
        mind  = dists.min()
        if mind > _IM_CONTACT_DIST:
            continue
        li, pi = np.unravel_index(dists.argmin(), dists.shape)
        la, pname = lig[li], prot_nms[pi]
        is_hbond = (mind <= _IM_HBOND_DIST and
                    la['elem'] in {'O', 'N'} and
                    (pname in _IM_HBOND_DONORS or pname in _IM_HBOND_ACC))
        charged  = pname[:2] in {'NH', 'NZ', 'NE', 'ND', 'OD', 'OE'}
        is_salt  = (not is_hbond and mind <= 4.5 and
                    la['elem'] in {'O', 'N'} and charged)
        contacts.append({'key': key, 'resname': res['resname'], 'resnum': res['resnum'],
                         'chain': res['chain'], 'dist': mind,
                         'is_hbond': is_hbond, 'is_salt': is_salt,
                         'lig_atom': la, 'prot_atom': pname,
                         'center': res['center']})
    return sorted(contacts, key=lambda x: x['dist'])


def _im_project(lig, contacts):
    """SVD project onto ligand principal plane; place residues on fixed ring."""
    lig_c3 = np.mean([a['pos'] for a in lig], axis=0)
    L = np.array([a['pos'] - lig_c3 for a in lig])
    if len(L) >= 2:
        _, _, Vt = np.linalg.svd(L, full_matrices=False)
        u1, u2 = Vt[0], Vt[1]
    else:
        u1, u2 = np.array([1.,0.,0.]), np.array([0.,1.,0.])
    lig_2d = np.array([[np.dot(a['pos']-lig_c3, u1),
                        np.dot(a['pos']-lig_c3, u2)] for a in lig])
    span = np.max(np.linalg.norm(lig_2d, axis=1)) if len(lig_2d) else 1.0
    scale = 3.0 / max(span, 0.5)
    lig_2d *= scale
    ring_r = 5.5
    res_2d = {}
    for c in contacts:
        v = c['center'] - lig_c3
        x2 = np.dot(v, u1) * scale
        y2 = np.dot(v, u2) * scale
        d  = np.sqrt(x2**2 + y2**2)
        if d < 0.1: d = ring_r
        res_2d[c['key']] = np.array([ring_r * x2/d, ring_r * y2/d])
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
        lig, pro = _im_parse_pdb(pdb_path)
        if not lig:
            return "Failed"
        contacts = _im_find_contacts(lig, pro)
        if not contacts:
            return "Failed"
        # Add itype field so shared renderer works
        for c in contacts:
            if c['is_hbond']:
                c['itype'] = 'hbond'
            elif c['is_salt']:
                c['itype'] = 'salt'
            elif c['resname'] in _IM_HYDROPHOBIC:
                c['itype'] = 'hydrophobic'
            else:
                c['itype'] = 'contact'
        lig_2d, res_2d = _im_project(lig, contacts)
        lig_2d = _im_separate_atoms(lig_2d)
        _im_render_diagram(lig_2d, lig, res_2d, contacts, out_png, mode='distance')
        return "Success"
    except Exception as _e:
        _fig_log(f"  [!] InteractionMap failed for {base_name}: {_e}")
        return "Failed"


# -------------------------------------------------------------------------------
# Step 6.5: Software manager
# -------------------------------------------------------------------------------
class SoftwareManager:
    """Handles verification and automated installation of visual tools."""

    def __init__(self):
        self.status = {}

    def check_all(self):
        _fig_log(f"\n{ConsoleColours.BOLD}>>> Checking Visualisation Engine Availability{ConsoleColours.ENDC}")

        # 1. PyMOL — prefer the binary in the running conda env over the system one
        # (system /usr/bin/pymol on Ubuntu 24.04 uses python3-pymol 2.5 which
        # calls `from imp import find_module`; imp was removed in Python 3.12)
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
            import plip
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
                _fig_log(f"  • PyMOL    : https://pymol.org/")
            if not self.status.get("PLIP"):
                _fig_log(f"  • PLIP     : pip install plip")
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
            _parts = pdb_path.stem.split('_')
            for _i, _p in enumerate(_parts[:-1]):
                _candidate = f"{_p}_{_parts[_i + 1]}"
                if _candidate in smi_map:
                    res_name = _parts[_i + 1]   # e.g. 'TFA'
                    break
        try:
            with open(pdb_path, 'r') as f:
                for line in f:
                    if (line.startswith("HETATM") or line.startswith("ATOM")) and len(line) > 21 and line[21] == 'L':
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
# Step 6.6: Phase 2 entry point
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

    # Collect all target directories
    dirs = sorted(ext_dir.glob("1_*_Raw_Complexes")) + sorted(ext_dir.glob("2_*_Prepared_Complexes"))
    ctrl_root = next((ext_dir / n for n in ["4_Controls", "Controls"] if (ext_dir / n).exists()), None)
    if ctrl_root:
        if (ctrl_root / "1_Raw").exists():
            dirs.append(ctrl_root / "1_Raw")
        if (ctrl_root / "2_Prepared").exists():
            dirs.append(ctrl_root / "2_Prepared")

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
    workers   = max(1, multiprocessing.cpu_count() - 2)

    _fig_log(f"\n{ConsoleColours.BOLD}Multi-Engine Figure Generation Pipeline (Concurrency: {workers}){ConsoleColours.ENDC}")

    # Build separate task queues:
    # - subprocess_tasks: PyMOL + PLIP (spawn external processes — thread-safe)
    # - imap_tasks: InteractionMap (matplotlib Agg — NOT thread-safe, run serially)
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
            subprocess_tasks.append(('PyMOL', _run_pymol, task_base + (lig_name, lig_num, has_f, sw_status, C)))
            subprocess_tasks.append(('PLIP',  _run_plip,  task_base + (sw_status, C)))
            imap_tasks.append((pdb.stem, pdb, fig_root))

    _plip_tasks  = [t for t in subprocess_tasks if t[0] == 'PLIP']
    _pymol_tasks = [t for t in subprocess_tasks if t[0] == 'PyMOL']
    # PyMOL workers capped at 2: ray-tracing is CPU-bound; high concurrency causes
    # CPU saturation, massively inflated render times, and timeout-related orphan processes.
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
        results[cname]['InteractionMap'] = status
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


if __name__ == "__main__":
    import time as _time
    _t0 = _time.perf_counter()
    main()
    _utils_mod.print_elapsed(_t0, "06_Top-N_Extraction_FAcDs.py")
