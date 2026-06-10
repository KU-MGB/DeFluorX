#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 08  |  MD Thermodynamics & QM/MM Engine
===============================================================================

Molecular dynamics trajectory analysis: near-attack conformation (NAC) geometry,
8-residue Dream Team catalytic machinery tracking, WaterMap thermodynamic
integration, Desmond EAF ligand dynamics, and QSite automation for QM/MM
SN2 reaction-coordinate scans.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 June 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 08_MD_Thermodynamics_QMMM_Engine_FAcDs.py
  Role          : Trajectory analysis engine; terminal computational step before
                  QM/MM (outputs ideal frame + QSite .inp files).
  Imports from  : 00_02_Project_Config_FAcDs.py  (CFG — all thresholds + tier metadata)
                  00_03_Project_Utils_FAcDs.py   (ConsoleColours, geometric utilities)
  Reads         : <Run>/7_Physics_Validation/MolecularDynamics/*Rank_N*/*-out.cms
                                                               /*_trj/   (any dir name containing Rank_N)
                                                               /*.eaf
                  <Run>/7_Physics_Validation/WaterMaps/*Rank_N*.csv  (Maestro WM export)
                  <Run>/7_Physics_Validation/WaterMaps/*Rank_N*/*_wm.maegz
                  <Run>/1_Boltz2_Production/7_Boltz2_FAcDs_Ranked_*.csv
                  <Run>/1_Boltz2_Production/6_Boltz2_FAcDs_Master_*.csv
  Writes        : <Run>/8_MD_Thermodynamics_Results/Rank_N_<Name>/
                    - <Name>_NAC_Data.csv          (per-frame geometry + DT)
                    - <Name>_NAC_Dashboard.png     (2-panel figure)
                    - <Name>_Ideal_Final.maegz      (best frame for QSite)
                    - <Name>_QSite_SN2.inp         (QM/MM scan input)
                  <Run>/8_MD_Thermodynamics_Results/08_MD_Master_Ranking.csv
  Upstream      : 07_SID_Post_Processing_FAcDs.py → produces *_SID-out.eaf consumed here
                  06_Top-N_Extraction_FAcDs.py   → provides ranked structures & IDs
                  02_Production_FAcDs.py         → master CSV with alignment maps
                  06_Top-N_Extraction_FAcDs.py   → visual reporting of hits (figures merged into Step 06)
  Downstream    : None (terminal step; QSite .inp feeds Schrödinger QSite/Jaguar)
────────────────────────────────────────────────────────────────────────────────

# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Schrödinger Python: MUST be executed using `run python3` within a
     Schrödinger environment (requires `schrodinger.analysis.topology`).
  2. Trajectory Compatibility: Strictly expects Desmond CMS/trj format; will
     not process GROMACS or AMBER trajectories without prior conversion.
  3. Memory & I/O: Reading large trajectories (1000+ frames) at stride 1 is
     extremely I/O intensive; recommend high-speed NVMe or local scratch disk.
  4. WaterMap Export: Relies on manual Maestro export of WaterMap CSV/MAE files
     matching the Rank_N naming convention.
───────────────────────────────────────────────────────────────────────────────

Usage:
    python 08_MD_Thermodynamics_QMMM_Engine_FAcDs.py Boltz-2_Run_20260309T085406Z
    python 08_MD_Thermodynamics_QMMM_Engine_FAcDs.py Boltz-2_Run_20260309T085406Z --stride 5 --ranks 3 --lig PFAS
    python 08_MD_Thermodynamics_QMMM_Engine_FAcDs.py Boltz-2_Run_20260309T085406Z --nuc 85 --base 250 --acid 112

Arguments:
    run_dir           Positional. Boltz-2 run folder name or prefix (e.g.
                      Boltz-2_Run_20260309T085406Z). Glob-expanded to match timestamp
                      suffix automatically. Resolves into 7_Physics_Validation.
    --dir    DIR      Alternative to positional (legacy). Same resolution logic.
    --lig    RESNAME  Ligand residue name in the CMS system file.  Default: LIG
    --stride N        Analyse every N-th trajectory frame (1 = all frames).
                      Higher values trade accuracy for speed.     Default: 1
    --ranks  N        Number of top-ranked MD jobs to process (Rank_1 … Rank_N).
                      Default: auto-detected from directories in MolecularDynamics/
    --nuc    RESNUM   Fallback nucleophile residue number if alignment map absent.
                      Default: 110  (FAcD canonical Asp110)
    --base   RESNUM   Fallback catalytic base residue number.
                      Default: 277  (FAcD canonical His277)
    --acid   RESNUM   Fallback catalytic acid residue number.
                      Default: 134  (FAcD canonical Asp134)
    --csv    PATH     Path to 02_Production_FAcDs.py master CSV for triad mapping and
                      alignment map. Auto-detected from sibling Boltz-2_Run_*
                      directories if omitted.

── Engine features ────────────────────────────────────────────────────────────
  1. Performance striding: configurable frame-sampling interval.
  2. 8-Residue Dream Team tracking: per-frame distances for all catalytic
     machinery (Nuc, Clamp1, Clamp2, Acid, StabH, StabW, StabY, Base),
     mapped from 02_Production_FAcDs.py master CSV alignment map.
  3. WaterMap CSV integration: Maestro-exported thermodynamic statistics
     (dG, dH, -TdS, occupancy, H-bond counts) added to master CSV.
  4. WaterMap spatial scoring: per-frame dG-weighted water blockade using
     site coordinates from *_wm.maegz.
  5. Desmond EAF integration: ligand surface area (MSA) and radius of
     gyration (RG) extracted from pl_interact_survey EAF files.
  6. Water blockade: dG-weighted count of waters obstructing the SN2 runway.
  7. Walden pre-organisation: improper dihedral check for TS flattening (×1.1).
  8. QSite automation: M06-2X/6-31+G(d,p) .inp generation for coordinate scan.
  9. 3D Smart-Lock: geometry-biased triad & fluorine-cradle detection.
 10. Rich progress bars and colour-coded PASS/FAIL NAC reporting.
 11. Master aggregation: 08_MD_Master_Ranking.csv.
───────────────────────────────────────────────────────────────────────────────

Scientific references
─────────────────────
  Boltz-2 structure prediction:
    Passaro, S. et al. (2025) bioRxiv 2025.06.14.659707.
    DOI: https://doi.org/10.1101/2025.06.14.659707
  ColabFold MSA server:
    Mirdita, M. et al. (2022) Nature Methods 19:679–682.
    DOI: https://doi.org/10.1038/s41592-022-01488-1
  FAcD crystal structure & catalytic mechanism (PDB 3R3U):
    Chan, P.W.Y., Yakunin, A.F., Edwards, E.A. & Pai, E.F. (2011) JACS 133:7461–7468.
    DOI: https://doi.org/10.1021/ja200277d
  DEHA4 experimental defluorination validation (Delftia acidovorans D4B):
    Farajollahi, S. et al. (2024) ACS Omega 9(26):28546–28555.
    DOI: https://doi.org/10.1021/acsomega.4c02517
  Dream Team 8-residue catalytic machinery (triad distances):
    Holmquist, M. (2000) Curr Protein Pept Sci 1:209–235.
    DOI: https://doi.org/10.2174/1389203003381405
    Verschueren, K.H.G. et al. (1993) Crystallographic analysis of the catalytic mechanism of
      haloalkane dehalogenase. Nature 363:693–698.
    DOI: https://doi.org/10.1038/363693a0
  NAC geometry (Lightstone–Bruice criteria):
    Lightstone, F.C. & Bruice, T.C. (1996) JACS 118:2595–2605.
    DOI: https://doi.org/10.1021/ja952589l
    Bruice, T.C. (2002) Acc Chem Res 35:139–148.
    DOI: https://doi.org/10.1021/ar0001665
    Hur, S. & Bruice, T.C. (2003) PNAS 100:12015–12020.
    DOI: https://doi.org/10.1073/pnas.1534873100
  Bürgi–Dunitz angle (auxiliary carbonyl-addition metric):
    Bürgi, H.B., Dunitz, J.D. & Shefter, E. (1973) JACS 95:5065–5067.
    DOI: https://doi.org/10.1021/ja00796a058
    Bürgi, H.B., Dunitz, J.D., Lehn, J.M. & Wipff, G. (1974) Tetrahedron 30:1563–1572.
    DOI: https://doi.org/10.1016/S0040-4020(01)90678-7
  WaterMap thermodynamic hydration-site analysis:
    Abel, R., Young, T., Farid, R., Berne, B.J. & Friesner, R.A. (2008) JACS 130:2817–2831.
    DOI: https://doi.org/10.1021/ja0771033
  Desmond MD engine:
    Bowers, K.J. et al. (2006) Scalable algorithms for molecular dynamics simulations on
      commodity clusters. SC 06: Proc. ACM/IEEE Conf. Supercomputing.
    DOI: https://doi.org/10.1109/SC.2006.54
  QSite QM/MM level of theory (M06-2X/6-31+G(d,p)):
    Zhao, Y. & Truhlar, D.G. (2008) Theor Chem Acc 120:215–241.
    DOI: https://doi.org/10.1007/s00214-007-0310-x
    Rosta, E., Klähn, M. & Warshel, A. (2006) J Phys Chem B 110:2934–2941.
    DOI: https://doi.org/10.1021/jp057109j
    Murphy, R.B. et al. (2000) J Comput Chem 21:1442–1457.
    DOI: https://doi.org/10.1002/1096-987X(200012)21:16<1442::AID-JCC3>3.0.CO;2-O
"""

import sys
import os

# ===============================================================================
# SCHRÖDINGER BOOTSTRAP
# ===============================================================================
# Auto-sets SCHRODINGER=/opt/schrodinger if the env var is absent.
# When invoked with plain `python`, re-invokes transparently via
# $SCHRODINGER/run so the Schrödinger Python interpreter is used.
# Both forms are equivalent:
#   python 08_MD_Thermodynamics_QMMM_Engine_FAcDs.py Boltz-2_Run_20260309T085406Z
#   $SCHRODINGER/run 08_MD_Thermodynamics_QMMM_Engine_FAcDs.py Boltz-2_Run_20260309T085406Z
import subprocess as _sp

if "SCHRODINGER" not in os.environ:
    os.environ["SCHRODINGER"] = "/opt/schrodinger"

try:
    from schrodinger.application.desmond.packages import traj, topo
    from schrodinger import structure
    from schrodinger.structutils import measure
except ImportError:
    _run_exe = os.path.join(os.environ["SCHRODINGER"], "run")
    if os.path.isfile(_run_exe):
        sys.exit(_sp.call([_run_exe, os.path.abspath(__file__)] + sys.argv[1:]))
    print(f"CRITICAL: Schrödinger suite not found at {os.environ['SCHRODINGER']}.\n"
          f"Set the SCHRODINGER environment variable to your installation path.")
    sys.exit(1)
# ===============================================================================

import re
import time
import argparse
import logging
import warnings
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import numpy as np
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import matplotlib.patheffects as pe
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

try:
    from rich.progress import (Progress, SpinnerColumn, BarColumn,
                               TextColumn, TimeElapsedColumn, MofNCompleteColumn)
    from rich.console import Console as _RichConsole
    from rich.table import Table as _RichTable
    _RICH_AVAILABLE = True
    _rcon = _RichConsole()
except ImportError:
    _RICH_AVAILABLE = False
    _rcon = None

import importlib.util as _ilu


def _load_module(name: str, path: Path):
    """Load a Python file as a module regardless of its filename."""
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_REPO_DIR  = Path(__file__).resolve().parent
_cfg_mod   = _load_module("ProjectConfig", _REPO_DIR / "00_02_Project_Config_FAcDs.py")
_utils_mod = _load_module("ProjectUtils",  _REPO_DIR / "00_03_Project_Utils_FAcDs.py")
CFG        = _cfg_mod.CFG()

# ConsoleColours sourced from 00_03_Project_Utils (single canonical definition).
ConsoleColours  = _utils_mod.ConsoleColours
SEPARATOR_HEAVY = _utils_mod.SEPARATOR_HEAVY
SEPARATOR_LIGHT = _utils_mod.SEPARATOR_LIGHT

# Console helpers and logging setup sourced from 00_03_Project_Utils.
_console_title          = _utils_mod.console_title
_console_info           = _utils_mod.console_info
_console_sep            = _utils_mod.console_separator
_setup_logging          = _utils_mod.setup_logging
get_mic_vector          = _utils_mod.get_mic_vector
_mic_dists_2d           = _utils_mod.mic_dists_2d
_mic_vecs_1d            = _utils_mod.mic_vecs_1d
clean_spines            = _utils_mod.clean_spines
_calc_improper_dihedral = _utils_mod.calculate_improper_dihedral
_ensure_box_3x3             = _utils_mod._ensure_box_3x3
find_nucleophile_od_fallback = _utils_mod.find_nucleophile_od_fallback

PLOT_LOCK = threading.Lock()


# ===============================================================================
# SECTION 1: GLOBAL CONSTANTS & CONFIGURATION
# ===============================================================================
# All thresholds sourced from 00_02_Project_Config_FAcDs.py (CFG).
# Fallback literals are numerically identical — activate only when CFG is
# unavailable (e.g., standalone testing outside the repository).

# ── § 1.1  NAC geometry thresholds (SN2 attack) ──────────────────────────────
THRESHOLD_STRICT_NAC_DIST   = CFG.NAC_DIST_STRICT    # 3.2 Å — strict nucleophile–C
THRESHOLD_STRICT_NAC_ANGLE  = CFG.NAC_ANGLE_STRICT   # 155°  — strict O–C–F angle
THRESHOLD_RELAXED_NAC_DIST  = CFG.NAC_DIST_RELAXED   # 3.8 Å — relaxed nucleophile–C
THRESHOLD_RELAXED_NAC_ANGLE = CFG.NAC_ANGLE_RELAXED  # 145°  — relaxed O–C–F angle
POCKET_RESIDENCY_DIST       = CFG.POCKET_RESIDENCY_DIST  # 8.0 Å — pocket-bound cutoff

# ── § 1.2  Catalytic triad integrity thresholds (MD-calibrated) ───────────────
# Crystal thresholds (4.5 Å NB / 7.0 Å BA) + 2.0 Å thermal-fluctuation buffer
# for 300 K Desmond simulations; see CFG §6 for derivation.
THRESHOLD_TRIAD_NB          = CFG.THRESHOLD_TRIAD_NB_MD   # 6.5 Å — Nuc–Base (crystal = 4.5 Å)
THRESHOLD_TRIAD_BA          = CFG.THRESHOLD_TRIAD_BA_MD   # 9.0 Å — Base–Acid (crystal = 7.0 Å)

# ── § 1.3  Physical chemistry scoring parameters ──────────────────────────────
_WALDEN_IMPROPER_MAX        = CFG.WALDEN_IMPROPER_MAX  # 15.0° — improper dihedral for TS geometry
_SOLVENT_RESTYPES           = CFG.SOLVENT_RESTYPES     # residue names for water blockade detection
_SCORE_W_DIST               = CFG.SCORE_DIST_WEIGHT    # 100.0 — per Å below relaxed distance
_SCORE_W_ANGLE              = CFG.SCORE_ANGLE_WEIGHT   #   5.0 — per ° above relaxed angle
_SCORE_W_BLOCK              = CFG.SCORE_BLOCKADE_WEIGHT #  50.0 — per water blockade unit
_SCORE_W_WM                 = CFG.SCORE_WATERMAP_WEIGHT #  10.0 — per kcal mol⁻¹ WaterMap dG unit

# ── § 1.4  8-Residue Dream Team reference mapping (PDB 3R3U / DEHA4) ─────────
# Reference residue numbers in canonical FAcD; alignment map from master CSV
# translates these to enzyme-specific sequential numbers per job.
DREAM_TEAM_REF = CFG.DREAM_TEAM_REFS

SEPARATOR = "-" * 80
logger    = None  # Initialised in main()


# ===============================================================================
# SECTION 2: CONSOLE WRAPPER FUNCTIONS
# ===============================================================================

def console_title(msg: str) -> None:
    _console_title(msg, logger)


def console_info(msg: str) -> None:
    _console_info(msg, logger)


def console_separator(heavy: bool = True) -> None:
    _console_sep(logger, heavy=heavy)


def _print_labeled(label: str, ansi_col: str, msg: str,
                   log_label: str | None = None, level: str = "info") -> None:
    """Print a prefixed status line with ANSI colour.
    Colour passes through pipeline tee to the terminal; file logger stays clean."""
    _log = log_label or label
    if _rcon:
        _rcon.print(f"  [bold]{label}[/bold] {msg}", markup=False)
    else:
        print(f"  {ansi_col}{label}\033[0m {msg}", flush=True)
    if logger:
        getattr(logger, level)(f"{_log} | {msg}")


def console_nac_pass(msg: str) -> None:
    _print_labeled("NAC PASS",  "\033[92m", msg)


def console_nac_geom_pass(msg: str) -> None:
    _print_labeled("GEOM PASS", "\033[93m", msg, log_label="GEOM PASS")


def console_nac_fail(msg: str) -> None:
    _print_labeled("NAC FAIL",  "\033[91m", msg)


def console_qmm_ready(msg: str) -> None:
    _print_labeled("QM/MM READY", "\033[94m", msg)


def console_scientific_alert(msg: str) -> None:
    _print_labeled("SCIENTIFIC ALERT", "\033[95m", msg, level="warning")


def console_watermap_warning(msg: str) -> None:
    _print_labeled("WATERMAP WARNING", "\033[93m", msg, level="warning")


# ===============================================================================
# SECTION 3: GEOMETRY & MATHEMATICAL UTILITIES
# ===============================================================================

# NOTE: This function is NOT equivalent to 00_03_Project_Utils.calculate_min_distance.
# It uses the Schrödinger frame API (frame.pos(idx)) rather than numpy arrays.
# The API divergence is intentional — required for Schrödinger/Maestro integration.
def calculate_min_distance(frame, indices_A: list, indices_B: list) -> float:
    """Minimum PBC-corrected distance between two atom index sets."""
    if not indices_A or not indices_B:
        return np.nan
    box = frame.box if hasattr(frame, 'box') else None
    min_d = float('inf')
    for a in indices_A:
        pos_a = frame.pos(a)
        for b in indices_B:
            d = np.linalg.norm(get_mic_vector(pos_a, frame.pos(b), box))
            if d < min_d:
                min_d = d
    return min_d


def extract_hybrid_smart_system(cms_model, tr, lig_resname: str,
                                mapped_nuc: int, mapped_base: int, mapped_acid: int):
    """
    3D geometry-based auto-detection of catalytic triad and fluorine cradle.

    Mapped residue numbers act as soft bias hints (a few-Å bonus within ±15
    residues — CFG.SMART_LOCK_BIAS_DIST), not hard constraints: a clearly closer
    geometric candidate can still win over the sequence-aligned hint. Hard safety
    nets remain SMART_LOCK_NUC_MAX_DIST and the Nuc–Base frame-0 sanity check.

    Returns (idx_nuc, idx_base, idx_acid, idx_cradle, cf_pairs) as atom index
    lists, or (None, None, None, None, None) on failure.
    """
    # 1. Ligand C-F pairs
    lig_atoms = cms_model.select_atom(f"(res.ptype '{lig_resname}')")
    if not lig_atoms:
        console_info(f"    {ConsoleColours.FAIL}[!] Ligand '{lig_resname}' not found.{ConsoleColours.ENDC}")
        return None, None, None, None, None

    cf_pairs = []
    for i in lig_atoms:
        atom = cms_model.atom[i]
        if atom.atomic_number == 9:
            for bond in atom.bond:
                if bond.atom2.atomic_number == 6:
                    cf_pairs.append((bond.atom2.index, atom.index))
    if not cf_pairs:
        console_info(f"    {ConsoleColours.FAIL}[!] No C-F bonds found in ligand.{ConsoleColours.ENDC}")
        return None, None, None, None, None

    # 2. Index all protein sidechain atoms by residue
    backbone_atoms = {'N', 'C', 'CA', 'O', 'H', 'HA'}
    res_dict = {}
    for atom in cms_model.atom:
        if atom.pdbname.strip() in backbone_atoms:
            continue
        res_key = f"{atom.chain}_{atom.resnum}"
        if res_key not in res_dict:
            res_dict[res_key] = {
                'ptype': atom.getResidue().pdbres.strip(),
                'O_idx': [], 'N_idx': [], 'heavy_idx': [],
                'resnum': atom.resnum, 'chain': atom.chain
            }
        if atom.atomic_number == 8:  res_dict[res_key]['O_idx'].append(atom.index)
        if atom.atomic_number == 7:  res_dict[res_key]['N_idx'].append(atom.index)
        if atom.atomic_number != 1:  res_dict[res_key]['heavy_idx'].append(atom.index)

    frame_0 = tr[0]
    lig_c_idxs = [c for c, f in cf_pairs]

    # 3. Nucleophile: closest ASP/GLU/SER oxygen to ligand C, biased by mapped_nuc
    best_nuc_key = None
    min_eff      = float('inf')
    actual_dist  = float('inf')
    for res_key, data in res_dict.items():
        if not data['O_idx'] or data['ptype'] not in {'ASP', 'GLU', 'ASH', 'GLH'}:
            continue
        d     = calculate_min_distance(frame_0, lig_c_idxs, data['O_idx'])
        bonus = CFG.SMART_LOCK_BIAS_DIST if (mapped_nuc and abs(data['resnum'] - mapped_nuc) <= CFG.SMART_LOCK_RESNUM_WINDOW) else 0.0
        if (d + bonus) < min_eff:
            min_eff = d + bonus; best_nuc_key = res_key; actual_dist = d

    if not best_nuc_key or actual_dist > CFG.SMART_LOCK_NUC_MAX_DIST:
        # ── Oδ Orientation Fallback ────────────────────────────────────────────
        # Primary geometry search exhausted.  Delegate to shared utility
        # find_nucleophile_od_fallback() (00_03_Project_Utils_FAcDs.py).
        # The function accepts plain NumPy arrays only; extract positions here
        # before calling so CMS atom-group objects never enter the utility.
        # Threshold aligns with CFG.NAC_ANGLE_RELAXED (BRAIN.md §7).
        _box       = frame_0.box if hasattr(frame_0, 'box') else None
        c_idx, f_idx = cf_pairs[0]
        c_pos_np   = np.array(frame_0.pos(c_idx), dtype=float)
        f_pos_np   = np.array(frame_0.pos(f_idx), dtype=float)
        # Build candidate list: (resnum, od1_np, od2_np) for each ASP/ASH residue.
        _asp_cands = []
        for rk, rd in res_dict.items():
            if rd['ptype'] not in {'ASP', 'ASH'}:
                continue
            od_idxs = [i for i in rd['O_idx']
                       if cms_model.atom[i].pdbname.strip() in ('OD1', 'OD2')]
            if not od_idxs:
                continue
            od1_np = np.array(frame_0.pos(od_idxs[0]), dtype=float)
            od2_np = np.array(frame_0.pos(od_idxs[1]), dtype=float) if len(od_idxs) > 1 else None
            _asp_cands.append((rk, od1_np, od2_np))
        fb_result = find_nucleophile_od_fallback(
            _asp_cands, c_pos_np, f_pos_np,
            CFG.SMART_LOCK_OD_FALLBACK_DIST, CFG.SMART_LOCK_OD_FALLBACK_ANGLE
        )
        if fb_result is not None:
            best_nuc_key, _fb_od_pos = fb_result
            actual_dist = float(np.linalg.norm(_fb_od_pos - c_pos_np))
            console_info(
                f"    {ConsoleColours.WARNING}↳ Oδ Fallback: Nuc → "
                f"{res_dict[best_nuc_key]['ptype']} {res_dict[best_nuc_key]['resnum']} "
                f"[Chain {res_dict[best_nuc_key]['chain']}]  "
                f"d={actual_dist:.2f} Å  "
                f"(primary Smart-Lock exhausted){ConsoleColours.ENDC}"
            )
        else:
            console_info(
                f"    {ConsoleColours.FAIL}[!] Nucleophile Not Found (fallback exhausted): "
                f"no ASP Oδ within {CFG.SMART_LOCK_OD_FALLBACK_DIST:.1f} Å "
                f"with backside geometry ≥ {CFG.SMART_LOCK_OD_FALLBACK_ANGLE:.0f}°.{ConsoleColours.ENDC}"
            )
            return None, None, None, None, None

    idx_nuc    = res_dict[best_nuc_key]['O_idx']
    nuc_chain  = res_dict[best_nuc_key]['chain']
    nuc_resnum = res_dict[best_nuc_key]['resnum']
    console_info(f"    {ConsoleColours.OKBLUE}↳ 3D Smart-Lock: Nuc  → "
                 f"{res_dict[best_nuc_key]['ptype']} {nuc_resnum} "
                 f"[Chain {nuc_chain}]{ConsoleColours.ENDC}")

    # 4. Base: closest HIS nitrogen to nucleophile, biased by mapped_base
    idx_base = []
    best_base_key = None; min_eff = float('inf')
    for res_key, data in res_dict.items():
        if not data['N_idx'] or data['ptype'] not in {'HIS', 'HIP', 'HIE', 'HID'}:
            continue
        d     = calculate_min_distance(frame_0, idx_nuc, data['N_idx'])
        bonus = CFG.SMART_LOCK_BIAS_DIST if (mapped_base and abs(data['resnum'] - mapped_base) <= CFG.SMART_LOCK_RESNUM_WINDOW) else 0.0
        bonus += CFG.SMART_LOCK_CHAIN_BIAS if data['chain'] == nuc_chain else 0.0
        if (d + bonus) < min_eff:
            min_eff = d + bonus; best_base_key = res_key
    if best_base_key:
        idx_base = res_dict[best_base_key]['N_idx']
        console_info(f"    {ConsoleColours.OKBLUE}↳ 3D Smart-Lock: Base → "
                     f"{res_dict[best_base_key]['ptype']} {res_dict[best_base_key]['resnum']} "
                     f"[Chain {res_dict[best_base_key]['chain']}]{ConsoleColours.ENDC}")

    # 5. Acid: closest ASP/GLU oxygen to base, biased by mapped_acid
    idx_acid = []
    if idx_base:
        best_acid_key = None; min_eff = float('inf')
        for res_key, data in res_dict.items():
            if (res_key == best_nuc_key or not data['O_idx']
                    or data['ptype'] not in {'ASP', 'GLU', 'ASH', 'GLH'}):
                continue
            d     = calculate_min_distance(frame_0, idx_base, data['O_idx'])
            bonus = CFG.SMART_LOCK_BIAS_DIST if (mapped_acid and abs(data['resnum'] - mapped_acid) <= CFG.SMART_LOCK_RESNUM_WINDOW) else 0.0
            bonus += CFG.SMART_LOCK_CHAIN_BIAS if data['chain'] == nuc_chain else 0.0
            if (d + bonus) < min_eff:
                min_eff = d + bonus; best_acid_key = res_key
        if best_acid_key:
            idx_acid = res_dict[best_acid_key]['O_idx']
            console_info(f"    {ConsoleColours.OKBLUE}↳ 3D Smart-Lock: Acid → "
                         f"{res_dict[best_acid_key]['ptype']} {res_dict[best_acid_key]['resnum']} "
                         f"[Chain {res_dict[best_acid_key]['chain']}]{ConsoleColours.ENDC}")

    # 6. Fluorine cradle: TRP/TYR heavy atoms within fluorine cradle radius of nucleophile
    idx_cradle = []
    _fcr = CFG.F_CRADLE_RADIUS
    for res_key, data in res_dict.items():
        if data['ptype'] in {'TRP', 'TYR'}:
            if calculate_min_distance(frame_0, idx_nuc, data['heavy_idx']) <= _fcr:
                idx_cradle.extend(data['heavy_idx'])

    console_info(f"    {ConsoleColours.OKBLUE}↳ 3D Smart-Lock: Fluorine Cradle → "
                 f"{len(idx_cradle)} TRP/TYR heavy atoms.{ConsoleColours.ENDC}")
    console_info(f"    {ConsoleColours.OKBLUE}↳ Polyfluorinated Engine: "
                 f"{len(cf_pairs)} C-F bonds tracked.{ConsoleColours.ENDC}")

    return idx_nuc, idx_base, idx_acid, idx_cradle, cf_pairs


# ===============================================================================
# SECTION 4: DATA HANDLING & LOADING FUNCTIONS
# ===============================================================================

def parse_mapping(map_str: str) -> dict:
    """Parses 'ASP110:ASP112 | HIS277:HIS280' alignment string into {ref_num: tgt_num}."""
    aln_dict = {}
    if not isinstance(map_str, str) or map_str.lower() == 'nan':
        return aln_dict
    for p in map_str.split(" | "):
        try:
            m = re.match(r"([A-Z]{3})(\d+):([A-Z]{3}|GAP)(\d+|GAP)", p)
            if m:
                _, ref_num, _, tgt_num = m.groups()
                if tgt_num != "GAP":
                    aln_dict[int(ref_num)] = int(tgt_num)
        except Exception:
            continue
    return aln_dict


def load_watermap_csv(wm_csv_path: Path) -> list:
    """
    Load WaterMap thermodynamic sites from a Maestro 'Analyse WaterMap'
    CSV export (Rank_N.csv).

    Expected columns: Site, Occupancy, dH, -TdS, dG, #HB(WW), #HB(PW), #HB(LW).
    Returns list of dicts without spatial coordinates — use for statistics only.
    """
    if not wm_csv_path.exists():
        console_info(f"    {ConsoleColours.WARNING}[!] WaterMap CSV not found: {wm_csv_path.name}{ConsoleColours.ENDC}")
        return []
    sites = []
    try:
        df = pd.read_csv(wm_csv_path)
        # Flexible column detection (strip units from names)
        def _col_base(c: str) -> str:
            return re.sub(r'\s*\([^)]*\)', '', c).strip()
        stripped = {_col_base(c): c for c in df.columns}

        def _find(keys):
            for k in keys:
                if k in stripped: return stripped[k]
            for k in keys:
                for raw_col in df.columns:
                    if k.lower() in raw_col.lower():
                        return raw_col
            return None

        site_col  = _find(['Site', '#', 'Num', 'Number'])
        dg_col    = _find(['dG', 'deltaG', 'Free Energy'])
        dh_col    = _find(['dH', 'deltaH', 'Enthalpy'])
        tds_col   = _find(['-TdS', 'TdS', 'Entropy'])
        occ_col   = _find(['Occupancy', 'Occ'])
        hbpw_col  = _find(['#HB(PW)', 'HB(PW)', 'HB_PW', 'Protein'])
        hblw_col  = _find(['#HB(LW)', 'HB(LW)', 'HB_LW', 'Ligand'])

        if not dg_col:
            console_info(f"    {ConsoleColours.WARNING}[!] WaterMap CSV: dG column not found in {wm_csv_path.name}{ConsoleColours.ENDC}")
            return []

        for i, row in df.iterrows():
            try:
                sites.append({
                    'site':     int(row[site_col])   if site_col and pd.notna(row[site_col])  else i,
                    'dG':       float(row[dg_col])   if pd.notna(row[dg_col])                 else 0.0,
                    'dH':       float(row[dh_col])   if dh_col  and pd.notna(row[dh_col])     else 0.0,
                    'TdS':      float(row[tds_col])  if tds_col and pd.notna(row[tds_col])    else 0.0,
                    'occupancy':float(row[occ_col])  if occ_col and pd.notna(row[occ_col])    else 0.0,
                    'hb_pw':    float(row[hbpw_col]) if hbpw_col and pd.notna(row[hbpw_col]) else 0.0,
                    'hb_lw':    float(row[hblw_col]) if hblw_col and pd.notna(row[hblw_col]) else 0.0,
                })
            except (ValueError, TypeError):
                continue

        console_info(f"    {ConsoleColours.OKGREEN}✔ WaterMap CSV: {len(sites)} thermodynamic sites loaded.{ConsoleColours.ENDC}")
    except Exception as e:
        console_info(f"    {ConsoleColours.WARNING}[!] WaterMap CSV load failed: {e}{ConsoleColours.ENDC}")
    return sites


def compute_wm_csv_stats(wm_data: list) -> dict:
    """
    Summary statistics from a Maestro WaterMap CSV export.

    The CSV is already filtered to the active-site region by Maestro, so no
    spatial cutoff is applied here. Stable sites (dG < 0) represent ordered
    waters that must be displaced — a catalytic penalty. Unstable sites
    (dG > 0) are energetically favourable for ligand entry.
    """
    empty = {
        'WM_N_Sites': 0, 'WM_N_Stable': 0, 'WM_N_Unstable': 0,
        'WM_Mean_dG': np.nan, 'WM_Min_dG': np.nan, 'WM_Max_dG': np.nan,
        'WM_Mean_dH': np.nan, 'WM_Mean_TdS': np.nan,
        'WM_Mean_HB_PW': np.nan, 'WM_Mean_HB_LW': np.nan,
    }
    if not wm_data:
        return empty
    dgs  = [s['dG']  for s in wm_data]
    dhs  = [s['dH']  for s in wm_data]
    tdss = [s['TdS'] for s in wm_data]
    hbpw = [s['hb_pw'] for s in wm_data]
    hblw = [s['hb_lw'] for s in wm_data]
    return {
        'WM_N_Sites':    len(wm_data),
        'WM_N_Stable':   sum(1 for d in dgs if d < 0),
        'WM_N_Unstable': sum(1 for d in dgs if d >= 0),
        'WM_Mean_dG':    round(float(np.mean(dgs)),  3),
        'WM_Min_dG':     round(float(np.min(dgs)),   3),
        'WM_Max_dG':     round(float(np.max(dgs)),   3),
        'WM_Mean_dH':    round(float(np.mean(dhs)),  3),
        'WM_Mean_TdS':   round(float(np.mean(tdss)), 3),
        'WM_Mean_HB_PW': round(float(np.mean(hbpw)), 3),
        'WM_Mean_HB_LW': round(float(np.mean(hblw)), 3),
    }


def load_watermap_sites(wm_path: Path) -> list:
    """
    Load WaterMap site positions and dG values from a .maegz file.
    Used for per-frame spatial scoring (blockade and dG weighting).
    """
    sites = []
    if not wm_path.exists():
        console_info(f"    {ConsoleColours.WARNING}[!] WaterMap .maegz not found: {wm_path.name}{ConsoleColours.ENDC}")
        return sites
    try:
        with structure.StructureReader(str(wm_path)) as reader:
            st = next(reader)
            for atom in st.atom:
                sites.append({
                    'pos': np.array(atom.xyz),
                    'dG':  atom.property.get('r_watermap_deltaG', 0.0),
                    'num': atom.property.get('i_watermap_site_num', 0)
                })
        console_info(f"    {ConsoleColours.OKGREEN}✔ WaterMap .maegz: {len(sites)} spatial sites loaded.{ConsoleColours.ENDC}")
    except Exception as e:
        console_info(f"    {ConsoleColours.WARNING}[!] WaterMap .maegz load failed: {e}{ConsoleColours.ENDC}")
    return sites


def extract_8residue_indices(cms_model, dream_mapped: dict) -> dict:
    """
    Returns {role: [sidechain_heavy_atom_indices]} for all Dream Team
    catalytic residues.

    dream_mapped: {role_name: target_resnum} built from alignment map +
    DREAM_TEAM_REF. Roles with no mapped residue number return empty list.

    Builds a resnum→atoms lookup once, then iterates 8 times — efficient
    for full solvated MD systems with tens of thousands of atoms.
    """
    backbone = frozenset({'N', 'C', 'CA', 'O', 'H', 'HA', 'OXT', 'H1', 'H2', 'H3'})
    # Pre-build residue number → [(pdbname, atomic_number, index)]
    res_atom_map: dict = {}
    for atom in cms_model.atom:
        rn = atom.resnum
        if rn not in res_atom_map:
            res_atom_map[rn] = []
        res_atom_map[rn].append((atom.pdbname.strip(), atom.atomic_number, atom.index))

    result = {role: [] for role in DREAM_TEAM_REF}
    for role, resnum in dream_mapped.items():
        if resnum is None or resnum not in res_atom_map:
            continue
        indices = [idx for pdb, an, idx in res_atom_map[resnum]
                   if pdb not in backbone and an != 1]
        result[role] = indices
        if indices:
            console_info(f"    {ConsoleColours.OKBLUE}↳ Dream Team [{role:8s}]: "
                         f"res {resnum:4d} → {len(indices)} sidechain heavy atoms{ConsoleColours.ENDC}")
        else:
            console_info(f"    {ConsoleColours.WARNING}↳ Dream Team [{role:8s}]: "
                         f"res {resnum} not found in CMS (may have been renumbered or removed by PrepWizard — check step 05 output).{ConsoleColours.ENDC}")
    return result


def find_eaf_file(job_folder: Path):
    """
    Locate protein-ligand EAF file with priority:
      P1: *-out*pl*.eaf  (Rank_1: explicit pl_interact_survey)
      P2: *SID-out.eaf   (Rank_2/5: simulation-ID output)
      P3: *out*.eaf      (any with 'out' in name)
      P4: largest *.eaf  (fallback — Rank_3/4 naming)
    Returns Path or None.
    """
    for pattern in ("*-out*pl*.eaf", "*SID-out.eaf", "*out*.eaf"):
        candidates = sorted(job_folder.glob(pattern))
        if candidates:
            return candidates[0]
    candidates = sorted(job_folder.glob("*.eaf"),
                        key=lambda p: p.stat().st_size, reverse=True)
    return candidates[0] if candidates else None


def load_eaf_scalar_series(eaf_path, keyword: str) -> np.ndarray:
    """
    Parse a per-frame scalar Result array from a Desmond EAF ASCII file.

    EAF format: {keyword = { ... Result = [v1 v2 ... vN] ... }}
    Only blocks whose Result array contains purely numeric values are matched
    (excludes the ProtLigInter contact matrix which contains strings).

    Returns np.ndarray of float32 (empty on failure).
    """
    if not eaf_path or not Path(eaf_path).exists():
        return np.array([], dtype=np.float32)
    try:
        with open(eaf_path) as f:
            content = f.read()
        # Find the keyword block — search from its opening brace
        kw_start = content.find(f'{{{keyword} = {{')
        if kw_start == -1:
            kw_start = content.find(f'{{{keyword}={{')
        if kw_start == -1:
            return np.array([], dtype=np.float32)
        # From keyword start, find Result = [numeric-only array]
        segment = content[kw_start: kw_start + 3_000_000]  # 3 MB window
        m = re.search(r'Result\s*=\s*\[([\d\.\-\s+eE]+)\]', segment)
        if m:
            vals = [float(x) for x in m.group(1).split()]
            return np.array(vals, dtype=np.float32)
    except Exception as e:
        if logger:
            logger.warning(f"EAF {keyword} load failed ({eaf_path}): {e}")
    return np.array([], dtype=np.float32)


def format_job_label(job_name: str, rank: int) -> str:
    """Formats verbose job names to 'Rx_ProteinName_LigandName'."""
    clean = job_name.replace("desmond_md_job_", "").replace("_Prepared", "").replace("_CONTROL", "")
    clean = re.sub(r"^Rank_\d+_", "", clean)
    parts = [p for p in clean.split('_') if p]
    if len(parts) >= 3:
        prot = parts[1]; lig = parts[-1]
    elif len(parts) == 2:
        prot = parts[0]; lig = parts[1]
    else:
        prot = clean[:10]; lig = "LIG"
    return f"R{int(rank)}_{prot}_{lig}"


def load_triad_mapping(csv_path: Path) -> dict:
    """Loads catalytic triad residue numbers and static metrics from master CSV."""
    if not csv_path or not csv_path.exists():
        return {}
    mapping = {}
    try:
        df = pd.read_csv(csv_path, dtype=str)
        triad_col = next(
            (c for c in ["Active_Site_Triad_Map", "Mapped_to_Control_Cat_Triad"] if c in df.columns), None)
        if not triad_col:
            return {}

        def get_col(keywords):
            for k in keywords:
                for col in df.columns:
                    if k.lower() in col.lower():
                        return col
            return None

        dist_col  = get_col(['nucleophile_ligand_distance', 'nac_dist', 'distance_nuc', 'dist_nuc', 'nuc_dist', 'nuc', 'distance', 'dist'])
        angle_col = get_col(['attack_angle', 'nac_angle', 'angle', 'ang'])
        cat_col   = get_col(['category', 'tier', 'grade'])
        conf_col  = get_col(['boltz_confidence', 'binding_probability', 'confidence', 'prob'])
        fluor_col = get_col(['fluorous', 'contacts', 'interaction'])
        hbond_col = get_col(['count_hydrogen_bond', 'hydrogen_bond', 'hbond', 'h_bond'])
        heavy_col = get_col(['ligand_heavy_atoms', 'heavy_atoms', 'heavy'])
        aln_col   = get_col(['Full_Sequence_Alignment_Map', 'alignment_map', 'sequence_map'])

        def safe_float(val):
            if pd.isna(val) or str(val).strip().lower() in ['nan', 'none', '']:
                return np.nan
            try:
                m = re.search(r"[-+]?\d*\.\d+|[-+]?\d+", str(val))
                return float(m.group()) if m else np.nan
            except Exception:
                return np.nan

        for _, row in df.iterrows():
            job_idx = str(row.get('job_index', '')).strip().zfill(7)
            if not job_idx or job_idx == '0000000':
                m = re.search(r"(\d{7})", str(row.get('job_name', '')))
                if m: job_idx = m.group(1)

            triad_str = str(row.get(triad_col, ''))
            if triad_str and triad_str.lower() != 'nan':
                nuc_m  = re.search(r"Nuc:[a-zA-Z]+(\d+)",  triad_str)
                base_m = re.search(r"Base:[a-zA-Z]+(\d+)", triad_str)
                acid_m = re.search(r"Acid:[a-zA-Z]+(\d+)", triad_str)
                if job_idx:
                    mapping[job_idx] = {
                        'nuc':  int(nuc_m.group(1))  if nuc_m  else None,
                        'base': int(base_m.group(1)) if base_m else None,
                        'acid': int(acid_m.group(1)) if acid_m else None,
                        'alignment_map': str(row[aln_col]) if aln_col and pd.notna(row.get(aln_col)) else '',
                        'static_data': {
                            'Static_Distance':    safe_float(row[dist_col])  if dist_col  else np.nan,
                            'Static_Angle':       safe_float(row[angle_col]) if angle_col else np.nan,
                            'Category':           str(row[cat_col]) if cat_col and pd.notna(row[cat_col]) else 'Unknown',
                            'Boltz_Confidence':   safe_float(row[conf_col])  if conf_col  else np.nan,
                            'Fluorous_Contacts':  safe_float(row[fluor_col]) if fluor_col else np.nan,
                            'Static_HBonds':      safe_float(row[hbond_col]) if hbond_col else np.nan,
                            'Ligand_Heavy_Atoms': safe_float(row[heavy_col]) if heavy_col else np.nan,
                        }
                    }
        return mapping
    except Exception as e:
        if logger: logger.error(f"Error loading triad mapping: {e}")
        return {}


def find_latest_master_csv(base_dir: Path):
    """
    Finds the most recent 6_Boltz2_FAcDs_Master_*.csv (or any *_Master*.csv fallback)
    in the run dir's 1_Boltz2_Production folder.
    base_dir is expected to be 7_Physics_Validation; parent is the Boltz-2_Run_* dir.
    """
    run_dir   = base_dir.parent  # Boltz-2_Run_*
    prod_dir  = run_dir / "1_Boltz2_Production"
    all_csvs  = []
    if prod_dir.exists():
        # Priority 1: canonical PFAS-27 master naming (6_Boltz2_FAcDs_Master_*.csv)
        all_csvs.extend(prod_dir.glob("6_Boltz2_FAcDs_Master_*.csv"))
        # Priority 2: any *_Master*.csv in the same dir
        if not all_csvs:
            all_csvs.extend(prod_dir.glob("*_Master*.csv"))
    # Fallback: search all sibling Boltz-2_Run_* dirs (handles multi-run setups)
    if not all_csvs:
        for sibling in run_dir.parent.glob("Boltz-2_Run_*"):
            if not sibling.is_dir(): continue
            sp = sibling / "1_Boltz2_Production"
            if sp.exists():
                all_csvs.extend(sp.glob("6_Boltz2_FAcDs_Master_*.csv"))
                all_csvs.extend(sp.glob("*_Master*.csv"))
    return sorted(all_csvs, key=lambda x: x.stat().st_mtime)[-1] if all_csvs else None


# ===============================================================================
# SECTION 5: VISUALISATION ENGINES
# ===============================================================================

def generate_individual_dashboard(df: pd.DataFrame, job_name: str,
                                  output_path: Path, stats: dict) -> None:
    """2-panel per-job dashboard: SN2 scatter and dual-trace anchoring time series."""
    with PLOT_LOCK:
        sns.set_theme(style="whitegrid", context="paper")
        plt.rcParams.update({'font.family': 'sans-serif'})

        fig = plt.figure(figsize=(15, 6.5))
        gs  = gridspec.GridSpec(1, 2, width_ratios=[1, 1.2], wspace=0.15)

        # Panel 1: SN2 scatter with catalytic zones
        ax1 = fig.add_subplot(gs[0])
        ax1.add_patch(plt.Rectangle(
            (0, THRESHOLD_RELAXED_NAC_ANGLE), THRESHOLD_RELAXED_NAC_DIST,
            180 - THRESHOLD_RELAXED_NAC_ANGLE, color='#74C476', alpha=0.3, zorder=0))
        ax1.add_patch(plt.Rectangle(
            (0, THRESHOLD_STRICT_NAC_ANGLE), THRESHOLD_STRICT_NAC_DIST,
            180 - THRESHOLD_STRICT_NAC_ANGLE, color='#006D2C', alpha=0.4, zorder=0))

        sc = ax1.scatter(df["NAC_Distance_A"], df["NAC_Angle_Deg"],
                         c=df["Frame"], cmap="viridis", s=20, alpha=0.7, edgecolor='none', zorder=2)
        if len(df.dropna(subset=["NAC_Distance_A", "NAC_Angle_Deg"])) > 10:
            sns.kdeplot(data=df, x="NAC_Distance_A", y="NAC_Angle_Deg", ax=ax1,
                        levels=5, color="#111111", linewidths=1.0, alpha=0.5, zorder=3,
                        warn_singular=False)

        ax1.axvline(THRESHOLD_RELAXED_NAC_DIST, color='#D55E00', linestyle='--', linewidth=2.5)
        ax1.axhline(THRESHOLD_RELAXED_NAC_ANGLE, color='#0072B2', linestyle='--', linewidth=2.5)
        ax1.set_title("Thermodynamic S_N2 Reaction Trajectory", fontsize=12, fontweight='bold', pad=10)
        ax1.set_xlabel("Nucleophile – Ligand Distance (Å)", fontweight='bold', fontsize=10)
        ax1.set_ylabel("Attack Angle: O–C–F (°)",           fontweight='bold', fontsize=10)
        ax1.set_ylim(0, 180); ax1.set_xlim(left=0)
        
        # Legend moved to bottom right and contains zones
        ax1.legend(handles=[
            Line2D([0], [0], color='#D55E00', linestyle='--', lw=2.5,
                   label=f'Relaxed Dist < {THRESHOLD_RELAXED_NAC_DIST}Å'),
            Line2D([0], [0], color='#0072B2', linestyle='--', lw=2.5,
                   label=f'Relaxed Angle > {THRESHOLD_RELAXED_NAC_ANGLE}°'),
            Patch(facecolor='#74C476', alpha=0.3, label='Relaxed S_N2 Zone'),
            Patch(facecolor='#006D2C', alpha=0.4, label='Strict S_N2 Zone'),
        ], loc='lower right', frameon=True, framealpha=0.95,
           edgecolor='#E2E8F0', fancybox=True, fontsize=9)
           
        cbar = plt.colorbar(sc, ax=ax1, pad=0.02)
        cbar.set_label("Simulation Frame (Time)", rotation=270, labelpad=15,
                       fontweight='bold', fontsize=10)
        clean_spines(ax1)

        # Panel 2: Lock & Key dual-trace
        ax2 = fig.add_subplot(gs[1])
        window = min(25, max(1, len(df) // 10))
        df = df.copy()
        df['Warhead_Smooth'] = df["NAC_Distance_A"].rolling(window=window, min_periods=1).mean()

        ax2.plot(df["Frame"], df["NAC_Distance_A"],   color='#D55E00', linewidth=1.0, alpha=0.15)
        ax2.plot(df["Frame"], df['Warhead_Smooth'],   color='#D55E00', linewidth=2.5, alpha=0.95,
                 label='Warhead Anchor (Nuc – LigC)')
        ax2.axhline(4.5, color='#D55E00', linestyle=':', linewidth=1.5, alpha=0.7)

        if "Tail_Cradle_Dist_A" in df.columns and not df["Tail_Cradle_Dist_A"].isna().all():
            df['Tail_Smooth'] = df["Tail_Cradle_Dist_A"].rolling(window=window, min_periods=1).mean()
            ax2.plot(df["Frame"], df["Tail_Cradle_Dist_A"], color='#0072B2', linewidth=1.0, alpha=0.15)
            ax2.plot(df["Frame"], df['Tail_Smooth'],        color='#0072B2', linewidth=2.5, alpha=0.95,
                     label='Tail Anchor (Cradle – LigF)')
            ax2.axhline(6.0, color='#0072B2', linestyle=':', linewidth=1.5, alpha=0.7)

        data_max = df["NAC_Distance_A"].max() if not df["NAC_Distance_A"].isna().all() else 12.0
        ax2.set_ylim(1.5, max(12.0, data_max * 1.4))
        ax2.set_title("Lock & Key: Dynamic Active Site Anchoring", fontsize=12, fontweight='bold', pad=10)
        ax2.set_xlabel("Simulation Frame",     fontweight='bold', fontsize=10)
        ax2.set_ylabel("Interaction Distance (Å)", fontweight='bold', fontsize=10)
        ax2.legend(loc='upper right', frameon=True, framealpha=1.0,
                   edgecolor='#E2E8F0', fancybox=True, fontsize=9)
        clean_spines(ax2)

        summary_text = (
            f"THERMODYNAMIC OUTCOME\n"
            f"-------------------------\n"
            f"Pocket Retention    : {stats.get('Pocket_Retention_Pct', 0.0):.1f}% of simulation\n"
            f"Relaxed Viability   : {stats.get('Catalytic_Viability_Pct', 0.0):.1f}% (While bound)\n"
            f"Strict Viability    : {stats.get('Strict_Viability_Pct', 0.0):.1f}% (While bound)\n"
            f"MD Min Distance     : {stats.get('MD_Min_NAC_Dist_A', 0.0):.2f} Å\n"
            f"MD Avg Attack Angle : {stats.get('MD_Avg_NAC_Angle_Deg', 0.0):.1f}°\n"
            f"Triad Integrity     : {stats.get('Triad_Integrity_Pct', 0.0):.1f}%\n"
            f"WM Stable Sites     : {stats.get('WM_N_Stable', 'N/A')}\n"
            f"WM Mean dG          : {stats.get('WM_Mean_dG', float('nan')):.2f} kcal/mol"
        )
        ax2.text(0.03, 0.96, summary_text, transform=ax2.transAxes,
                 fontsize=9, fontfamily='monospace', va='top', ha='left', zorder=10,
                 bbox=dict(facecolor='white', edgecolor='#CBD5E1',
                           boxstyle='round,pad=0.6', alpha=1.0))

        # Title removed as requested by the user, but metadata is written to the log
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

    if logger:
        logger.info(
            f"    [Dashboard Meta] Generated dashboard for {job_name} at {output_path.name}.\n"
            f"      - Panel 1: S_N2 Reaction Trajectory scatter plot of Attack Angle (O-C-F) against Nucleophile-Ligand Distance.\n"
            f"        Points represent simulation frames colored from start (purple/dark) to end (yellow/light).\n"
            f"        Shaded green regions define the Catalytic Zones (Relaxed: distance <= {THRESHOLD_RELAXED_NAC_DIST}A, angle >= {THRESHOLD_RELAXED_NAC_ANGLE}°; "
            f"Strict: distance <= {THRESHOLD_STRICT_NAC_DIST}A, angle >= {THRESHOLD_STRICT_NAC_ANGLE}°).\n"
            f"      - Panel 2: Lock & Key active site anchoring time series. Shows the nucleophile-carbon warhead distance (orange) "
            f"and cradle-fluorine tail anchor distance (blue) over simulation frames to assess binding durability."
        )


def generate_global_comparative_dashboard(out_dir: Path, df_master: pd.DataFrame) -> None:
    """3-panel: violin (dist), violin (angle), scatter landscape with catalytic zones."""
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update({'font.family': 'sans-serif'})

    all_data = []
    for _, row in df_master.iloc[::-1].iterrows():
        csv_path = next(out_dir.glob(f"**/{row['Job_Name']}_NAC_Data.csv"), None)
        if csv_path and csv_path.exists():
            df_job = pd.read_csv(csv_path)
            df_job['Job']  = format_job_label(row['Job_Name'], row['Dynamic_Rank'])
            df_job['Rank'] = row['Dynamic_Rank']
            all_data.append(df_job)

    if not all_data:
        return
    combined_df   = pd.concat(all_data, ignore_index=True)
    
    # Consistent colour map based on ascending Dynamic_Rank
    sorted_df = df_master.sort_values('Dynamic_Rank', ascending=True)
    sorted_labels = [format_job_label(row['Job_Name'], row['Dynamic_Rank']) for _, row in sorted_df.iterrows()]
    job_colour_map = dict(zip(sorted_labels, sns.color_palette("tab20", n_colors=len(sorted_labels))))

    # Vectorised piecewise distance transformation function
    def transform_distance(x):
        return np.where(x <= 4.0, x, 4.0 + np.log(np.maximum(x - 3.0, 1e-9)))

    combined_df = combined_df.copy()
    combined_df["NAC_Distance_A_Transformed"] = transform_distance(combined_df["NAC_Distance_A"])

    fig = plt.figure(figsize=(22, 10))
    gs  = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.6], wspace=0.15)

    # Panel 1: Nucleophile Distance Distribution violin plot
    ax1 = fig.add_subplot(gs[0])
    sns.violinplot(data=combined_df, y="Job", x="NAC_Distance_A_Transformed", ax=ax1,
                   palette=job_colour_map, inner="quartile", linewidth=1.2)
    ax1.axvspan(0, transform_distance(THRESHOLD_RELAXED_NAC_DIST), color='#009E73', alpha=0.15, zorder=0)
    ax1.axvline(transform_distance(THRESHOLD_RELAXED_NAC_DIST), color='#D55E00', linestyle='--', linewidth=2)
    ax1.set_title("Nucleophile Distance Distribution\n(Distribution over MD Frames)", fontweight='bold')
    ax1.set_xlabel("Distance (Å) [non-linear scale]"); ax1.set_ylabel("")
    
    dist_ticks = [0, 1, 2, 3, 4, 5, 10, 15, 20, 30, 40, 50]
    ax1.set_xticks([transform_distance(t) for t in dist_ticks])
    ax1.set_xticklabels([str(t) for t in dist_ticks], rotation=90)
    ax1.set_xlim(left=0)
    
    fig.canvas.draw()
    for lbl in ax1.get_yticklabels():
        lbl.set_color(job_colour_map.get(lbl.get_text(), '#333333'))
        lbl.set_fontweight('bold'); lbl.set_fontsize(11)
    clean_spines(ax1)

    # Panel 2: Attack Angle Distribution violin plot
    ax2 = fig.add_subplot(gs[1], sharey=ax1)
    sns.violinplot(data=combined_df, y="Job", x="NAC_Angle_Deg", ax=ax2,
                   palette=job_colour_map, inner="quartile", linewidth=1.2)
    ax2.axvspan(THRESHOLD_RELAXED_NAC_ANGLE, 180, color='#009E73', alpha=0.15, zorder=0)
    ax2.axvline(THRESHOLD_RELAXED_NAC_ANGLE, color='#0072B2', linestyle='--', linewidth=2)
    ax2.set_title("Attack Angle Distribution\n(Distribution over MD Frames)", fontweight='bold')
    ax2.set_xlabel("Angle (°)"); ax2.set_ylabel("")
    ax2.set_xlim(0, 180) # strictly physical bounds
    ax2.tick_params(labelleft=False)
    ax2.tick_params(axis='x', labelrotation=90)
    clean_spines(ax2)

    # Panel 3: Global Catalytic Landscape scatter plot
    ax3 = fig.add_subplot(gs[2])
    sns.scatterplot(data=combined_df, x="NAC_Distance_A_Transformed", y="NAC_Angle_Deg",
                    hue="Job", palette=job_colour_map, s=35, alpha=0.5,
                    edgecolor='none', ax=ax3, legend=False)
    
    # Draw zones using transformed distance widths
    relaxed_w = transform_distance(THRESHOLD_RELAXED_NAC_DIST)
    strict_w  = transform_distance(THRESHOLD_STRICT_NAC_DIST)
    
    ax3.add_patch(plt.Rectangle(
        (0, THRESHOLD_RELAXED_NAC_ANGLE), relaxed_w,
        180 - THRESHOLD_RELAXED_NAC_ANGLE, color='#74C476', alpha=0.3, zorder=0))
    ax3.add_patch(plt.Rectangle(
        (0, THRESHOLD_STRICT_NAC_ANGLE), strict_w,
        180 - THRESHOLD_STRICT_NAC_ANGLE, color='#006D2C', alpha=0.4, zorder=0))
        
    ax3.axvline(transform_distance(THRESHOLD_RELAXED_NAC_DIST), color='#D55E00', linestyle='--', linewidth=2)
    ax3.axhline(THRESHOLD_RELAXED_NAC_ANGLE, color='#0072B2', linestyle='--', linewidth=2)
    ax3.set_title("Global Catalytic Landscape (colours match Y-axis labels)", fontweight='bold')
    ax3.set_xlabel("Nucleophile – Ligand Distance (Å) [non-linear scale]")
    ax3.set_ylabel("Attack Angle: O–C–F (°)")
    
    ax3.set_xticks([transform_distance(t) for t in dist_ticks])
    ax3.set_xticklabels([str(t) for t in dist_ticks], rotation=90)
    ax3.set_xlim(0, transform_distance(50))
    ax3.set_ylim(0, 180)
    
    # Legend at bottom left containing both lines and zones
    ax3.legend(handles=[
        Line2D([0], [0], color='#D55E00', linestyle='--', lw=2,
               label=f'Distance < {THRESHOLD_RELAXED_NAC_DIST}Å'),
        Line2D([0], [0], color='#0072B2', linestyle='--', lw=2,
               label=f'Angle > {THRESHOLD_RELAXED_NAC_ANGLE}°'),
        Patch(facecolor='#74C476', alpha=0.3, label='Relaxed S_N2 Zone'),
        Patch(facecolor='#006D2C', alpha=0.4, label='Strict S_N2 Zone'),
    ], loc='lower left', frameon=True, framealpha=0.95,
       edgecolor='#E2E8F0', fancybox=True)
    clean_spines(ax3)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        plt.tight_layout()
    out_path = out_dir / "09_MD_Comparative_Analysis.png"
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    console_info(f"    Comparative Dashboard Saved : {out_path.resolve()}")


def generate_viability_bar_chart(out_dir: Path, df_master: pd.DataFrame) -> None:
    """Catalytic viability summary bar chart: parallel nested bars inside a master track.

    Sub-bars show percentage of total simulation time:
      - Pocket Retention (faint opacity of candidate colour)
      - Triad Integrity (medium opacity of candidate colour)
      - Relaxed Catalysis (high opacity of candidate colour)
      - Strict Catalysis (solid candidate colour)

    Right annotation panel shows per-candidate: WaterMap mean ΔG, minimum NAC
    distance, and Dream Team residues matched.
    """
    sns.set_theme(style="whitegrid", context="paper")
    plt.rcParams.update({'font.family': 'sans-serif'})
    from matplotlib.colors import to_rgba

    df_plot = df_master.sort_values('Dynamic_Rank', ascending=True).copy()

    # Consistent colour map based on ascending Dynamic_Rank (tab20)
    sorted_labels = [format_job_label(r['Job_Name'], r['Dynamic_Rank']) for _, r in df_plot.iterrows()]
    job_colour_map = dict(zip(sorted_labels, sns.color_palette("tab20", n_colors=len(sorted_labels))))

    n_rows  = len(df_plot)
    row_h   = 0.85
    fig_h   = max(4.0, n_rows * row_h + 1.8)
    fig, (ax, ax_ann) = plt.subplots(
        1, 2, figsize=(15, fig_h),
        gridspec_kw={'width_ratios': [3, 1], 'wspace': 0.04})

    bar_h       = 0.72
    y_positions = list(range(n_rows))

    for i, (_, row) in enumerate(df_plot.iterrows()):
        y = y_positions[i]
        job_label = format_job_label(row['Job_Name'], row['Dynamic_Rank'])
        rank_color = job_colour_map.get(job_label, '#475569')

        pocket = float(row.get('Pocket_Retention_Pct', 0.0) or 0.0)
        triad_o  = float(row.get('Triad_Integrity_Pct',     0.0) or 0.0)
        viab_r   = float(row.get('Catalytic_Viability_Pct', 0.0) or 0.0)
        viab_s   = float(row.get('Strict_Viability_Pct',    0.0) or 0.0)

        # Convert pocket-bound percentages to absolute percentage of total simulation time
        triad_total   = triad_o * pocket / 100.0
        relaxed_total = viab_r * pocket / 100.0
        strict_total  = viab_s * pocket / 100.0

        # Background master container representing 100% of simulation
        bg_face = to_rgba(rank_color, alpha=0.05)
        bg_edge = to_rgba(rank_color, alpha=0.60)
        ax.barh(y, 100.0, height=bar_h, color=bg_face, edgecolor=bg_edge, linewidth=1.2, zorder=1)

        # Plot the 4 parallel sub-bars inside the track with progressive opacities of the rank color
        sub_bars = [
            (pocket,        to_rgba(rank_color, alpha=0.35), -0.21, '*'),
            (triad_total,   to_rgba(rank_color, alpha=0.55), -0.07, '**'),
            (relaxed_total, to_rgba(rank_color, alpha=0.75),  0.07, '***'),
            (strict_total,  to_rgba(rank_color, alpha=1.0),   0.21, '****')
        ]
        
        sub_bar_h = 0.13
        for val, col, offset, star in sub_bars:
            if val > 1e-4:
                ax.barh(y + offset, val, height=sub_bar_h, color=col, edgecolor='none', zorder=2)
                # Value label (in solid black)
                if val >= 8.0:
                    ax.text(val - 1.0, y + offset, f'{val:.1f}%',
                            ha='right', va='center', fontsize=6.5, fontweight='bold',
                            color='black', zorder=4)
                else:
                    ax.text(val + 0.5, y + offset, f'{val:.1f}%',
                            ha='left', va='center', fontsize=6.5, fontweight='bold',
                            color='black', zorder=4)
            else:
                # Value label for 0% (in solid black)
                ax.text(0.5, y + offset, "0.0%",
                        ha='left', va='center', fontsize=6.5, fontweight='bold',
                        color='black', zorder=4)
            
            # Place the category identifier star(s) after the 100% mark (e.g. at x = 103.0) in black colour
            ax.text(103.0, y + offset, star, ha='left', va='center',
                    fontsize=7, fontweight='bold', color='black', zorder=4)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [format_job_label(r['Job_Name'], r['Dynamic_Rank']) for _, r in df_plot.iterrows()],
        fontsize=10.5, fontweight='bold')
    ax.set_xlim(0, 112)
    ax.set_ylim(-0.65, n_rows - 0.35)
    ax.invert_yaxis()
    ax.set_xlabel("Percentage of Simulation Time (%)", fontweight='bold', fontsize=11)
    ax.axvline(100, color='#94A3B8', linestyle=':', linewidth=1.0, alpha=0.6, zorder=1)
    clean_spines(ax)

    # Set y-axis tick label colors to match job colors
    fig.canvas.draw()
    for lbl in ax.get_yticklabels():
        lbl.set_color(job_colour_map.get(lbl.get_text(), '#333333'))

    # ── Annotation panel ──────────────────────────────────────────────────────
    ax_ann.set_xlim(0, 1)
    ax_ann.set_ylim(-0.65, n_rows - 0.35)
    ax_ann.invert_yaxis()
    ax_ann.set_yticks([])
    ax_ann.set_xticks([])
    for sp in ax_ann.spines.values():
        sp.set_visible(False)
    ax_ann.axvline(0.0, color='#CBD5E1', linewidth=0.8)

    hdr_y = -0.45
    ax_ann.text(0.20, hdr_y, 'WM ΔG\n(kcal/mol)', ha='center', va='center',
                 fontsize=8, fontweight='bold', color='#334155')
    ax_ann.text(0.55, hdr_y, 'min d_NAC\n(Å)',     ha='center', va='center',
                 fontsize=8, fontweight='bold', color='#334155')
    ax_ann.text(0.87, hdr_y, 'DT\n(n)',            ha='center', va='center',
                 fontsize=8, fontweight='bold', color='#334155')

    for i, (_, row) in enumerate(df_plot.iterrows()):
        y    = y_positions[i]
        wm   = row.get('WM_Mean_dG',        float('nan'))
        dist = row.get('MD_Min_NAC_Dist_A',  float('nan'))
        dt   = row.get('Dream_Team_Mapped',  0)

        try:    wm_str   = f'{float(wm):.2f}'
        except: wm_str   = 'N/A'
        try:    dist_str = f'{float(dist):.2f}'
        except: dist_str = 'N/A'
        try:    dt_str   = str(int(float(dt or 0)))
        except: dt_str   = '–'

        dist_col = ('#15803D' if (dist_str != 'N/A' and float(dist) < 3.8)
                    else ('#EA580C' if (dist_str != 'N/A' and float(dist) < 5.0)
                    else '#DC2626'))

        ax_ann.text(0.20, y, wm_str,   ha='center', va='center', fontsize=9,
                    color='#334155', fontweight='bold')
        ax_ann.text(0.55, y, dist_str, ha='center', va='center', fontsize=9,
                    color=dist_col,  fontweight='bold')
        ax_ann.text(0.87, y, dt_str,   ha='center', va='center', fontsize=9,
                    color='#6D28D9', fontweight='bold')

    # ── Legend ────────────────────────────────────────────────────────────────
    from matplotlib.patches import Patch
    c_neutral = '#475569'
    legend_handles = [
        Patch(facecolor=to_rgba(c_neutral, alpha=0.35), edgecolor='none', label='Pocket Retention (*)'),
        Patch(facecolor=to_rgba(c_neutral, alpha=0.55), edgecolor='none', label='Triad Integrity (Total) (**)'),
        Patch(facecolor=to_rgba(c_neutral, alpha=0.75), edgecolor='none', label='Relaxed Catalysis (Total) (***)'),
        Patch(facecolor=to_rgba(c_neutral, alpha=1.0),  edgecolor='none', label='Strict Catalysis (Total) (****)'),
    ]
    ax.legend(handles=legend_handles, loc='lower left', bbox_to_anchor=(0.0, 1.02),
              frameon=True, framealpha=0.95, edgecolor='#CBD5E1', fontsize=7.5, ncol=4, columnspacing=0.8)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        plt.tight_layout()
    out_path = out_dir / "10_MD_Viability_Summary.png"
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    console_info(f"    Viability Bar Chart Saved   : {out_path.resolve()}")


# ===============================================================================
# SECTION 6: PHYSICAL CHEMISTRY MODULES
# ===============================================================================

def _blockade_vec(nuc_pos: np.ndarray, lig_c_pos: np.ndarray,
                  sol_positions: np.ndarray, box, wm_sites,
                  block_r: float, match_r: float) -> float:
    """Vectorized water blockade — replaces check_water_blockade when positions
    are already loaded from the pre-load cache.  Processes all solvent atoms
    simultaneously with numpy instead of looping in Python."""
    runway_vec = get_mic_vector(lig_c_pos, nuc_pos, box)
    runway_len = float(np.linalg.norm(runway_vec))
    if runway_len < 1e-3 or sol_positions.shape[0] == 0:
        return 0.0
    unit_r = runway_vec / runway_len

    # PBC-corrected nuc→solvent vectors for all waters at once
    vecs = sol_positions - nuc_pos[None, :]                     # (n_sol, 3)
    if box is not None:
        try:
            b3    = _ensure_box_3x3(box)
            inv_b = np.linalg.inv(b3)
            frac  = vecs @ inv_b
            frac -= np.round(frac)
            vecs  = frac @ b3
        except Exception:
            pass

    projs = vecs @ unit_r                                        # (n_sol,)
    mask  = (projs > 0.5) & (projs < runway_len - 0.5)
    if not np.any(mask):
        return 0.0

    perp_lens = np.linalg.norm(vecs[mask] - np.outer(projs[mask], unit_r), axis=1)
    blocking  = perp_lens < block_r
    if not np.any(blocking):
        return 0.0

    if not wm_sites:
        return float(np.sum(blocking))

    blocking_pos = sol_positions[mask][blocking]                # (n_block, 3)
    total = 0.0
    for sp in blocking_pos:
        best_d  = float('inf')
        best_dg = 0.0
        for site in wm_sites:
            d = np.linalg.norm(get_mic_vector(site['pos'], sp, box))
            if d < best_d:
                best_d = d; best_dg = site['dG']
        weight = (1.0 / (1.0 + np.exp(best_dg)) * 2.0) if best_d < match_r else 1.0
        total += weight
    return total


def check_water_blockade(frame, nuc_pos, lig_c_pos, box, sol_indices,
                         wm_sites=None) -> float:
    """Checks if explicit waters obstruct the SN2 runway; weights by WaterMap dG."""
    runway_vec = get_mic_vector(lig_c_pos, nuc_pos, box)
    runway_len = np.linalg.norm(runway_vec)
    if runway_len < 1e-3:
        return 0.0
    unit_runway = runway_vec / runway_len
    total = 0.0
    _block_r = CFG.WATERMAP_BLOCKADE_RADIUS
    _match_r = CFG.WATERMAP_MATCH_RADIUS
    for s_idx in sol_indices:
        sol_pos = frame.pos(s_idx)
        vec_ns  = get_mic_vector(sol_pos, nuc_pos, box)
        proj    = np.dot(vec_ns, unit_runway)
        if 0.5 < proj < (runway_len - 0.5):
            perp = vec_ns - proj * unit_runway
            if np.linalg.norm(perp) < _block_r:
                weight = 1.0
                if wm_sites:
                    min_d = float('inf')
                    best_dg = 0.0
                    for site in wm_sites:
                        d = np.linalg.norm(get_mic_vector(site['pos'], sol_pos, box))
                        if d < min_d:
                            min_d = d
                            best_dg = site['dG']
                    if min_d < _match_r:
                        weight = 1.0 / (1.0 + np.exp(best_dg)) * 2.0
                total += weight
    return total


def calculate_walden_bonus(frame, c_idx, f_idx, other_subst_indices, box) -> float:
    """Returns 1.1 if transition-state flattening (|improper dihedral| < threshold), else 1.0."""
    if len(other_subst_indices) < 2 or _calc_improper_dihedral is None:
        return 1.0
    try:
        angle = _calc_improper_dihedral(
            frame.pos(c_idx), frame.pos(f_idx),
            frame.pos(other_subst_indices[0]), frame.pos(other_subst_indices[1]),
            box,
        )
        if abs(angle) < _WALDEN_IMPROPER_MAX:
            return 1.1
    except Exception:
        pass
    return 1.0


# ===============================================================================
# SECTION 7: QSITE AUTOMATION
# ===============================================================================

def generate_qsite_inputs(maegz_path: Path, job_name: str,
                          nuc_num, stab_f_num, lig_c_idx, nuc_o_idx,
                          base_num=None, acid_num=None,
                          lig_charge: int = None) -> Path:
    """
    Writes a QSite coordinate-scan .inp for M06-2X/6-31+G(d,p) SN2 QM/MM.

    QM region includes the full catalytic triad + fluoride stabiliser:
      LIG   — substrate (bond-making/breaking atoms)
      Nuc   — Asp110 sidechain (nucleophilic oxygen)
      Base  — His277 sidechain (proton-transfer base; must be QM for relay step)
      Acid  — Asp134 sidechain (proton-transfer acid; must be QM for relay step)
      StabH — His155 sidechain (fluoride stabiliser)

    Excluding Base (His277) and Acid (Asp134) from QM would force the proton-
    transfer half of the SN2 mechanism onto the MM force field, making full-
    mechanism reaction-coordinate scans chemically meaningless.
    """
    inp_path   = maegz_path.parent / f"{job_name}_QSite_SN2.inp"
    qm_regions = ['  QM_REGION asl="res.ptype LIG"']
    if nuc_num:
        qm_regions.append(f'  QM_REGION asl="res.num {nuc_num}  AND sidechain"')
    if base_num:
        qm_regions.append(f'  QM_REGION asl="res.num {base_num} AND sidechain"')
    if acid_num:
        qm_regions.append(f'  QM_REGION asl="res.num {acid_num} AND sidechain"')
    if stab_f_num:
        qm_regions.append(f'  QM_REGION asl="res.num {stab_f_num} AND sidechain"')

    # Resolve total QM region charge by summing the ligand charge and catalytic sidechain charges.
    # Asp110 (Nuc) is deprotonated (-1). His277 (Base) and His155 (StabH) are neutral (0).
    # Asp134 (Acid) is protonated in the resting state (0).
    if lig_charge is None:
        jn_upper = job_name.upper()
        lig_charge = next(
            (v for k, v in CFG.LIGAND_QM_CHARGES.items() if k.upper() in jn_upper),
            CFG.QSITE_CHARGE,
        )
    qm_charge = lig_charge
    if nuc_num:
        qm_charge -= 1  # Asp110 is deprotonated (-1)
    # His277, His155, and protonated Asp134 are neutral (0), contributing no additional charge.

    # Coordinate scan: Nu–C distance from NAC start → product.
    # Start = CFG.QSITE_SCAN_START (3.5 Å, the strict NAC distance).
    # Step  = CFG.QSITE_SCAN_STEP  (−0.1 Å per point).
    # Steps = CFG.QSITE_SCAN_NSTEPS (23 points → 3.5 → 1.3 Å; 3.5 + −0.1×22 = 1.3).
    _scan_start  = CFG.QSITE_SCAN_START
    _scan_step   = CFG.QSITE_SCAN_STEP
    _scan_nsteps = CFG.QSITE_SCAN_NSTEPS
    content = (
        "# QSite Coordinate Scan: SN2 Dehalogenation\n"
        "# Level of Theory: M06-2X / 6-31+G(d,p)\n"
        "BEGIN_QSITE_INPUT\n"
        f"  STRUCTURE_FILE {maegz_path.name}\n"
        "  METHOD dft\n"
        f"  FUNCTIONAL {CFG.QSITE_FUNCTIONAL}\n"
        f"  BASIS {CFG.QSITE_BASIS_SET}\n"
        f"  CHARGE {qm_charge}\n"
        f"  MULTIPLICITY {CFG.QSITE_MULT}\n"
        + "\n".join(qm_regions) + "\n"
        "  JOB_TYPE scan\n"
        f"  SCAN_COORD distance {nuc_o_idx} {lig_c_idx} "
        f"{_scan_start} {_scan_step} {_scan_nsteps}\n"
        # Implicit solvation (screening-grade — see CFG.QSITE_SOLVATION note;
        # use explicit-solvent QM/MM-FEP for publication-grade barriers).
        f"  SOLVATION_METHOD {CFG.QSITE_SOLVATION}\n"
        "  MAX_ITER 100\n"
        "END_QSITE_INPUT\n"
    )
    with open(inp_path, "w") as f:
        f.write(content)
    return inp_path


# ===============================================================================
# SECTION 8: CORE ANALYSIS ENGINE
# ===============================================================================

# Number of threads used for parallel DTR frame pre-loading.
# Each thread opens its own read_traj instance (thread-safe).
# For HDD: keep at 2-4 (seek contention); for SSD: use cpu_count().
_N_PRELOAD_WORKERS = min(4, max(1, (os.cpu_count() or 4) - 2))

# Solvent sphere radius (Å) around the nucleophile for water blockade.
# Only solvent within this radius at frame 0 is tracked — reduces per-frame
# work from ~10,000 solvent atoms to ~100-200 while covering the SN2 runway.
_SOL_SPHERE_RADIUS = 20.0


def _load_frame_chunk(chunk: list, trj_path: str, atom_indices: list,
                      t_start: float, sim_span: float, n_tr: int) -> list:
    """Worker: load a contiguous chunk of frames from a fresh trajectory reader.
    Each thread has its own reader to avoid shared file-handle races.
    Returns list of (fi, f_idx, pos_array, box, frame_time)."""
    local_tr = traj.read_traj(trj_path)
    results  = []
    for fi, f_idx in chunk:
        frame = local_tr[f_idx]
        # Batch position read: one C-extension call per frame
        pos = frame.pos(atom_indices)                  # (n_atoms, 3) ndarray
        box = frame.box if hasattr(frame, 'box') else None
        t   = frame.time if hasattr(frame, 'time') else (t_start + f_idx * sim_span / n_tr)
        results.append((fi, f_idx, pos, box, t))
    return results


def _eaf_at(series: np.ndarray, frame_t: float, t_start: float, eaf_dt: float) -> float:
    """Return EAF scalar value at the MD frame time, interpolated by index."""
    if len(series) == 0:
        return np.nan
    idx = max(0, min(int((frame_t - t_start) / eaf_dt), len(series) - 1))
    return float(series[idx])


def process_single_job(rank: int, work_dir: Path, df_ranked: pd.DataFrame,
                       master_out_dir: Path, lig_resname: str, stride: int,
                       triad_override: dict = None,
                       fallback_nuc:  int = DREAM_TEAM_REF.get('Nuc', 110),
                       fallback_base: int = DREAM_TEAM_REF.get('Base', 277),
                       fallback_acid: int = DREAM_TEAM_REF.get('Acid', 134)):
    """Orchestrates the full analysis pipeline for one MD trajectory."""
    try:
        row = df_ranked[df_ranked['Scientific_Rank'] == rank].iloc[0]
    except (IndexError, KeyError):
        console_info(f"    {ConsoleColours.WARNING}[!] No ranked entry for Rank {rank}. Skipping.{ConsoleColours.ENDC}")
        return None

    job_name   = row['job_name']
    # Flexible discovery: try every directory pattern Desmond/pipeline may produce
    _md_root   = work_dir / "MolecularDynamics"
    _candidates = [
        _md_root / f"desmond_md_job_Rank_{rank}",
        _md_root / f"Results_MD_Simulation_Rank_{rank}",
        _md_root / f"desmond_md_Rank_{rank}",
        _md_root / f"md_job_Rank_{rank}",
    ]
    # Also glob for any directory containing "Rank_{rank}" (covers custom names)
    _candidates += sorted(_md_root.glob(f"*Rank_{rank}*")) if _md_root.exists() else []
    job_folder = next((p for p in _candidates if p.is_dir()), None)
    if job_folder is None:
        console_info(f"    {ConsoleColours.FAIL}[!] Job folder missing for Rank {rank} in {_md_root}{ConsoleColours.ENDC}")
        return None

    _dir_label  = re.sub(r'_\d{1,3}_[A-Za-z][A-Za-z0-9]+$', '', job_name)
    job_out_dir = master_out_dir / f"Rank_{rank}_{_dir_label}"
    job_out_dir.mkdir(parents=True, exist_ok=True)
    console_info(f"Processing Rank {rank} [Stride={stride}]: {ConsoleColours.OKBLUE}{job_name}{ConsoleColours.ENDC}")

    # ── Load trajectory ────────────────────────────────────────────────────────
    # First try flat lookup (single-stage jobs: *-out.cms and *_trj at folder root).
    # Multi-stage restart jobs store each stage in a *_N-out.tgz archive; extract
    # on demand, then scan recursively for the final stage CMS and all trj segments.
    import tarfile as _tarfile

    def _auto_extract_tgz(folder: Path) -> None:
        """Extract *.tgz archives in folder only when their content is absent."""
        for _tgz in sorted(folder.glob("*-out.tgz")):
            try:
                with _tarfile.open(str(_tgz)) as _tf:
                    _first_entry = _tf.getnames()[0] if _tf.getnames() else ""
                    _first_path  = folder / _first_entry.split('/')[0]
                    if _first_path.exists():
                        continue
                    console_info(f"      [+] Extracting {_tgz.name} ...")
                    _tf.extractall(path=str(folder))
            except Exception as _te:
                console_info(f"      [!] TGZ extract failed ({_tgz.name}): {_te}")

    _cms_flat = list(job_folder.glob("*-out.cms"))
    _trj_flat = list(job_folder.glob("*_trj"))

    if not _cms_flat or not _trj_flat:
        _auto_extract_tgz(job_folder)
        # Recursively collect all *-out.cms / *_trj from extracted stage subdirs
        def _seg_key(p: Path, suffix: str) -> int:
            _m = re.search(r'_(\d+)' + re.escape(suffix), p.name)
            return int(_m.group(1)) if _m else 0
        _cms_flat = sorted(job_folder.glob("**/*-out.cms"),
                           key=lambda p: _seg_key(p, '-out.cms'))
        _trj_flat = sorted(job_folder.glob("**/*_trj"),
                           key=lambda p: _seg_key(p, '_trj'))

    if not _cms_flat or not _trj_flat:
        console_info(f"    {ConsoleColours.FAIL}[!] Missing .cms or _trj in {job_folder.name}{ConsoleColours.ENDC}")
        return None

    cms_path = _cms_flat[-1]   # final stage CMS (highest segment number)
    msys_model, cms_model = topo.read_cms(str(cms_path))

    if len(_trj_flat) == 1:
        tr = traj.read_traj(str(_trj_flat[0]))
    else:
        # Multi-segment restart: concatenate all trajectory segments in order
        _seg_frames: list = []
        for _seg_trj in _trj_flat:
            _seg_frames.extend(traj.read_traj(str(_seg_trj)))
        tr = _seg_frames
        console_info(f"      [+] Multi-segment trajectory: {len(_trj_flat)} segments → "
                     f"{len(tr):,} total frames")

    # ── WaterMap spatial sites (maegz) — flexible naming discovery ───────────
    _wm_root = work_dir / "WaterMaps"
    # Find any subdirectory containing "Rank_{rank}" (any capitalisation/separator)
    _wm_dir_candidates = (
        sorted(_wm_root.glob(f"*[Rr]ank*{rank}*")) if _wm_root.exists() else []
    )
    wm_maegz = None
    for _wmd in _wm_dir_candidates:
        if not _wmd.is_dir():
            continue
        # Accept *_wm.maegz or *wm*.maegz anywhere in the dir tree (depth ≤ 2)
        for _pat in ("*_wm.maegz", "*wm*.maegz", "*.maegz"):
            _hits = list(_wmd.glob(_pat)) + list(_wmd.glob(f"**/{_pat}"))
            if _hits:
                wm_maegz = _hits[0]; break
        if wm_maegz:
            break
    wm_sites = load_watermap_sites(wm_maegz)

    # ── WaterMap CSV statistics — flexible naming discovery ───────────────────
    wm_csv_path = None
    if _wm_root.exists():
        _csv_hits = (
            sorted(_wm_root.glob(f"*[Rr]ank*{rank}*.csv")) +
            sorted(_wm_root.glob(f"*[Rr]ank_{rank}.csv"))
        )
        if _csv_hits:
            wm_csv_path = _csv_hits[0]
    wm_csv_data = load_watermap_csv(wm_csv_path)

    # ── Desmond EAF: ligand dynamics (surface area + radius of gyration) ──────
    # Requires systemd-oomd to be masked before running (prevents OOM kill near end):
    #   sudo systemctl stop systemd-oomd && sudo systemctl mask systemd-oomd.socket
    # After job completes, restore with:
    #   sudo systemctl unmask systemd-oomd.socket && sudo systemctl start systemd-oomd
    eaf_path = find_eaf_file(job_folder)
    eaf_msa  = load_eaf_scalar_series(eaf_path, "Molecular_Surface_Area")
    eaf_rg   = load_eaf_scalar_series(eaf_path, "Rad_Gyration")
    # Estimate EAF time step: EAF samples at ~1 ps; use trajectory time range
    n_tr     = len(tr)
    t_start  = tr[0].time  if hasattr(tr[0], 'time') else 0.0
    t_end    = tr[-1].time if hasattr(tr[-1], 'time') else float(n_tr)
    sim_span = max(t_end - t_start, 1.0)
    eaf_dt   = sim_span / len(eaf_msa) if len(eaf_msa) > 0 else 1.0

    # ── Resolve residue number hints: alignment map > triad_override > CLI ─────
    aln_dict     = parse_mapping(row.get('Full_Sequence_Alignment_Map', ''))
    to           = triad_override or {}
    nuc_hint  = aln_dict.get(DREAM_TEAM_REF.get('Nuc', 110))  or to.get('nuc')  or fallback_nuc
    base_hint = aln_dict.get(DREAM_TEAM_REF.get('Base', 277)) or to.get('base') or fallback_base
    acid_hint = aln_dict.get(DREAM_TEAM_REF.get('Acid', 134)) or to.get('acid') or fallback_acid
    stab_f_num   = aln_dict.get(DREAM_TEAM_REF.get('Stab_H', 155))

    # ── 3D geometry-based Smart-Lock ───────────────────────────────────────────
    idx_nuc, idx_base, idx_acid, idx_cradle, cf_pairs = extract_hybrid_smart_system(
        cms_model, tr, lig_resname, nuc_hint, base_hint, acid_hint)
    if not idx_nuc or not cf_pairs:
        console_info(f"    {ConsoleColours.FAIL}[!] Smart-Lock failed. Skipping.{ConsoleColours.ENDC}")
        return None

    # ── Nucleophile centroid at frame 0 (reference for WM statistics) ──────────
    _f0         = tr[0]

    warhead_c = list(set(c for c, f in cf_pairs))
    lig_f     = list(set(f for c, f in cf_pairs))
    nuc_num   = cms_model.atom[idx_nuc[0]].resnum
    # Base and Acid residue numbers from Smart-Lock — needed for QM region expansion.
    base_num  = cms_model.atom[idx_base[0]].resnum if idx_base else None
    acid_num  = cms_model.atom[idx_acid[0]].resnum if idx_acid else None

    # ── Mapping sanity check: Nuc–Base distance in frame 0 (GEM-BUG-1 Fix 2) ───
    # PrepWizard residue renumbering can silently mis-map the triad; catch it here
    # before any per-frame analysis runs.
    if idx_nuc and idx_base:
        nuc_pos_frame0  = np.mean([_f0.pos(i) for i in idx_nuc],  axis=0)
        base_pos_frame0 = np.mean([_f0.pos(i) for i in idx_base], axis=0)
        nuc_base_frame0_dist = np.linalg.norm(nuc_pos_frame0 - base_pos_frame0)
        if nuc_base_frame0_dist > CFG.SMART_LOCK_MAPPING_SANITY_DIST:
            console_info(
                f"  [!] MAPPING ERROR: Nuc–Base distance in frame 0 = "
                f"{nuc_base_frame0_dist:.1f} Å — exceeds sanity limit "
                f"({CFG.SMART_LOCK_MAPPING_SANITY_DIST:.0f} Å). "
                f"PrepWizard residue renumbering suspected. "
                f"Check index_offset.json for this job."
            )
            return None

    # ── 8-Residue Dream Team index extraction ──────────────────────────────────
    # Smart-Lock resnums serve as fallbacks when the alignment map lacks an entry,
    # ensuring Nuc/Base/Acid Dream Team columns are populated even for distant homologs.
    _sl_fallback = {'Nuc': nuc_num, 'Base': base_num, 'Acid': acid_num}
    dream_mapped = {
        role: (_v if (_v := aln_dict.get(ref_num)) is not None else _sl_fallback.get(role))
        for role, ref_num in DREAM_TEAM_REF.items()
    }
    console_info("  Mapping Dream Team catalytic machinery:")
    dt_indices = extract_8residue_indices(cms_model, dream_mapped)
    n_mapped   = sum(1 for v in dt_indices.values() if v)
    console_info(f"    {ConsoleColours.OKGREEN}✔ {n_mapped}/{len(DREAM_TEAM_REF)} Dream Team residues mapped.{ConsoleColours.ENDC}")

    # Filter out ptypes with special characters (e.g. '/') that break ASL parsing.
    _safe_restypes = [r for r in _SOLVENT_RESTYPES if r.isalnum() or '_' in r]
    _sol_asl       = " OR ".join(f"res.ptype {r}" for r in _safe_restypes)
    sol_indices    = list(cms_model.select_atom(f"({_sol_asl}) AND a.el O"))

    # ===============================================================================
    # Trajectory pre-loading: batch-read all strided frames into RAM
    # Replaces per-frame frame.pos() calls with fast numpy array indexing.
    # ===============================================================================
    _t_pre = time.time()

    # 1. Filter solvent to a sphere around the active site (95%+ reduction).
    _nuc_cen_0 = (np.mean([tr[0].pos(i) for i in idx_nuc], axis=0)
                  if idx_nuc else None)
    if _nuc_cen_0 is not None and sol_indices:
        _sol_use = [s for s in sol_indices
                    if np.linalg.norm(np.array(tr[0].pos(s)) - _nuc_cen_0)
                       <= _SOL_SPHERE_RADIUS]
    else:
        _sol_use = []

    # 2. Walden substituents for all warhead carbons (needed for dihedral).
    _walden_subs = sorted(set(
        b.atom2.index for c_idx in warhead_c
        for b in cms_model.atom[c_idx].bond
        if b.atom2.atomic_number not in (1, 9)   # exclude H and F
    ))

    # 3. Build the master atom list to pre-load.
    _preload_set = (
        list(idx_nuc) + list(warhead_c) + list(lig_f) +
        list(idx_base or []) + list(idx_acid or []) +
        list(idx_cradle or []) + _walden_subs +
        [a for idxs in dt_indices.values() for a in idxs] +
        _sol_use
    )
    _preload_atoms = sorted(set(_preload_set))
    _a2l           = {a: i for i, a in enumerate(_preload_atoms)}   # atom_idx → local
    _preload_list  = list(_preload_atoms)                            # for frame.pos()
    _n_pre         = len(_preload_atoms)

    # 4. Pre-compute local index arrays (avoids per-frame dict lookups).
    _li_nuc    = [_a2l[a] for a in idx_nuc]
    _li_c      = [_a2l[a] for a in warhead_c]
    _li_f      = [_a2l[a] for a in lig_f]
    _li_base   = [_a2l[a] for a in idx_base]   if idx_base   else []
    _li_acid   = [_a2l[a] for a in idx_acid]   if idx_acid   else []
    _li_cradle = [_a2l[a] for a in idx_cradle] if idx_cradle else []
    _li_sol    = [_a2l[a] for a in _sol_use]
    _li_dt     = {role: [_a2l[a] for a in idxs]
                  for role, idxs in dt_indices.items() if idxs}
    # Map warhead C atom index → local indices of its bonded F atoms
    _li_cf     = {_a2l[c]: [_a2l[f] for (cc, f) in cf_pairs if cc == c]
                  for c in warhead_c}
    # Walden: local index of each C and its substituent pair
    _li_wsubs  = {_a2l[c]: [_a2l[s] for s in _walden_subs
                             if any(b.atom2.index == s
                                    for b in cms_model.atom[c].bond)]
                  for c in warhead_c}

    # 5. Chunked frame streaming — load / analyse / free _CHUNK_SIZE frames at a time.
    #    Peak RAM is O(_CHUNK_SIZE × _n_pre × 3 × 8 bytes) regardless of trajectory
    #    length; stride=1 over 100 k frames is safe within 64 GB.
    _CHUNK_SIZE  = 5_000           # ~60 MB per chunk at 500 atoms
    _stride_list = list(range(0, len(tr), stride))
    _nf          = len(_stride_list)
    _n_chunks    = (_nf + _CHUNK_SIZE - 1) // _CHUNK_SIZE

    console_info(
        f"  Chunked streaming: {_nf} frames × {_n_pre} atoms "
        f"({len(_sol_use)} solvent) | chunk={_CHUNK_SIZE} | {_n_chunks} chunks | "
        f"setup {time.time()-_t_pre:.1f}s"
    )

    _t_loop = time.time()

    def _iter_frames():
        """Yield (fi, f_idx, pos, box, frame_t) one frame at a time.

        Loads _CHUNK_SIZE frames into a temporary array, yields them all, then
        frees the array before moving to the next chunk.  This keeps peak memory
        bounded to a single chunk rather than the full trajectory.
        """
        for _ci in range(_n_chunks):
            _cs = _ci * _CHUNK_SIZE
            _ce = min(_cs + _CHUNK_SIZE, _nf)
            _cz = _ce - _cs
            _cpos = np.empty((_cz, _n_pre, 3), dtype=np.float64)
            _cbox = [None] * _cz
            _ctim = [0.0]  * _cz
            for _li in range(_cz):
                _fx    = _stride_list[_cs + _li]
                _frame = tr[_fx]
                _cpos[_li] = _frame.pos(_preload_list)
                _cbox[_li] = _frame.box  if hasattr(_frame, 'box')  else None
                _ctim[_li] = (_frame.time if hasattr(_frame, 'time')
                               else t_start + _fx * sim_span / n_tr)
            for _li in range(_cz):
                yield _cs + _li, _stride_list[_cs + _li], _cpos[_li], _cbox[_li], _ctim[_li]
            del _cpos, _cbox, _ctim
            _elapsed_c = time.time() - _t_loop
            _done_c    = _ce
            _rate_c    = _done_c / _elapsed_c if _elapsed_c > 1e-6 else 0.0
            _eta_c     = int((_nf - _done_c) / _rate_c) if _rate_c > 0.0 else 0
            console_info(
                f"  Chunk {_ci + 1}/{_n_chunks}: {_done_c}/{_nf} frames "
                f"({_rate_c:.0f} fr/s, ETA ~{_eta_c}s)"
            )

    def _r3(v):
        """Round to 3 dp; return NaN for None or non-finite float."""
        try:
            f = float(v)
            return round(f, 3) if not np.isnan(f) else np.nan
        except (TypeError, ValueError):
            return np.nan

    results = []
    best_score      = -1e9
    ideal_frame_idx = -1
    ideal_geom      = {}
    n_pocket = n_relaxed = n_strict = n_triad = n_nac = 0
    _wm_radius = CFG.WATERMAP_SITE_RADIUS

    # ===============================================================================
    # Per-frame analysis loop (chunked streaming — _CHUNK_SIZE frames at a time)
    # ===============================================================================
    for fi, f_idx, _p, box, frame_t in _iter_frames():

        # ── Position arrays from cache (numpy indexing — no API calls) ──────────
        _pos_nuc_f = _p[_li_nuc]  if _li_nuc   else np.empty((0, 3))
        _pos_c_f   = _p[_li_c]    if _li_c     else np.empty((0, 3))
        _pos_f_f   = _p[_li_f]    if _li_f     else np.empty((0, 3))

        # ── Warhead: closest nucleophile O to C-F carbon ──────────────────────
        if _pos_nuc_f.size and _pos_c_f.size:
            _nc_d    = _mic_dists_2d(_pos_nuc_f, _pos_c_f, box)
            _flat    = int(np.argmin(_nc_d))
            _nl, _cl = divmod(_flat, len(warhead_c))
            min_nuc_dist = float(_nc_d[_nl, _cl])
            best_nuc_idx   = idx_nuc[_nl]
            best_ca_idx   = warhead_c[_cl]
        else:
            min_nuc_dist = float('inf')
            best_nuc_idx = best_ca_idx = None

        if min_nuc_dist <= POCKET_RESIDENCY_DIST:
            n_pocket += 1

        # ── O–C–F attack angle (backside-attack geometry) ─────────────────────
        max_ang = 0.0; best_f_idx = None
        if best_nuc_idx is not None and _pos_f_f.size:
            _li_c_best  = _a2l[best_ca_idx]
            _pos_c_best = _p[_li_c_best]
            _bf_li      = _li_cf.get(_li_c_best, [])         # local F indices for this C
            _bonded_f   = [f for c, f in cf_pairs if c == best_ca_idx]
            if _bf_li and _bonded_f:
                _v_cf_b  = np.array([get_mic_vector(_p[fli], _pos_c_best, box)
                                     for fli in _bf_li])
                _cf_lens = np.linalg.norm(_v_cf_b, axis=1)
                _bond_mask = _cf_lens > 1e-6
                if np.any(_bond_mask):
                    _v_cn = get_mic_vector(_p[_a2l[best_nuc_idx]], _pos_c_best, box)
                    _n_cn = np.linalg.norm(_v_cn)
                    if _n_cn > 1e-6:
                        _vf_v = _v_cf_b[_bond_mask]
                        _lf_v = _cf_lens[_bond_mask]
                        _dots = np.dot(_vf_v, _v_cn) / (_lf_v * _n_cn)
                        _angs = np.degrees(np.arccos(np.clip(_dots, -1.0, 1.0)))
                        _bi   = int(np.argmax(_angs))
                        max_ang        = float(_angs[_bi])
                        best_f_idx = _bonded_f[int(np.where(_bond_mask)[0][_bi])]

        # ── Tail anchor: fluorine to halide cradle ────────────────────────────
        dist_tail = np.nan
        if _li_cradle and _pos_f_f.size:
            dist_tail = float(np.min(_mic_dists_2d(_pos_f_f, _p[_li_cradle], box)))

        # ── Catalytic triad integrity ─────────────────────────────────────────
        nb_dist = ba_dist = float('inf')
        _pos_base_f = _p[_li_base] if _li_base else np.empty((0, 3))
        if _pos_nuc_f.size and _pos_base_f.size:
            nb_dist = float(np.min(_mic_dists_2d(_pos_nuc_f, _pos_base_f, box)))
        if _pos_base_f.size and _li_acid:
            ba_dist = float(np.min(_mic_dists_2d(_pos_base_f, _p[_li_acid], box)))
        triad_ok = (nb_dist <= THRESHOLD_TRIAD_NB and ba_dist <= THRESHOLD_TRIAD_BA)
        if triad_ok and min_nuc_dist <= POCKET_RESIDENCY_DIST:
            n_triad += 1

        if min_nuc_dist <= THRESHOLD_RELAXED_NAC_DIST and max_ang >= THRESHOLD_RELAXED_NAC_ANGLE:
            n_nac += 1

        # ── Dream Team per-frame distances to warhead carbon ─────────────────
        dt_dists: dict = {}
        for role, li_role in _li_dt.items():
            if li_role and _pos_c_f.size:
                dt_dists[role] = float(np.min(_mic_dists_2d(_p[li_role], _pos_c_f, box)))
            else:
                dt_dists[role] = np.nan

        # Proton relay: Nuc–Base and Base–Acid residue–residue distances
        dt_nuc_base = dt_base_acid = np.nan
        _li_dt_n = _li_dt.get('Nuc');  _li_dt_b = _li_dt.get('Base')
        _li_dt_a = _li_dt.get('Acid')
        if _li_dt_n and _li_dt_b:
            dt_nuc_base  = float(np.min(_mic_dists_2d(_p[_li_dt_n], _p[_li_dt_b], box)))
        if _li_dt_b and _li_dt_a:
            dt_base_acid = float(np.min(_mic_dists_2d(_p[_li_dt_b], _p[_li_dt_a], box)))

        # ── EAF ligand dynamics (interpolated to frame time) ──────────────────
        frame_msa = _eaf_at(eaf_msa, frame_t, t_start, eaf_dt)
        frame_rg  = _eaf_at(eaf_rg,  frame_t, t_start, eaf_dt)

        # ── Water blockade (vectorized over pre-loaded solvent positions) ──────
        _pos_sol_f = _p[_li_sol] if _li_sol else np.empty((0, 3))
        blockades  = (
            _blockade_vec(_p[_a2l[best_nuc_idx]], _p[_a2l[best_ca_idx]],
                          _pos_sol_f, box, wm_sites,
                          CFG.WATERMAP_BLOCKADE_RADIUS, CFG.WATERMAP_MATCH_RADIUS)
            if best_nuc_idx else 0.0
        )

        # ── Walden TS flattening (improper dihedral, inline with cached pos) ──
        walden_flat = False
        walden = 1.0
        if best_f_idx and best_ca_idx:
            _wsub_li = _li_wsubs.get(_a2l[best_ca_idx], [])
            if len(_wsub_li) >= 2 and _calc_improper_dihedral is not None:
                try:
                    angle = _calc_improper_dihedral(
                        _p[_a2l[best_ca_idx]], _p[_a2l[best_f_idx]],
                        _p[_wsub_li[0]], _p[_wsub_li[1]], box)
                    if abs(angle) < _WALDEN_IMPROPER_MAX:
                        walden = 1.1; walden_flat = True
                except Exception:
                    pass

        # ── WaterMap per-frame thermodynamic contribution ─────────────────────
        wm_score = 0.0
        if best_ca_idx:
            _pos_c_wm = _p[_a2l[best_ca_idx]]
            for site in wm_sites:
                d = np.linalg.norm(get_mic_vector(site['pos'], _pos_c_wm, box))
                if d < _wm_radius:
                    wm_score += (_wm_radius - d) * site['dG']

        # ── Composite frame score (QM/MM frame selection only) ────────────────
        # NOTE: This score is intentionally linear and unbounded — it is used
        # exclusively to rank MD frames for QSite input selection, NOT for
        # comparison with the Step 02 sigmoid-based ranking (which is bounded
        # 0–1 over structures). The two scores are mathematically incompatible
        # and must not be cross-compared. The Walden 1.1× bonus is a geometric
        # heuristic (TS D₃ₕ pre-organisation), not an energy contribution.
        _productive = (np.clip(THRESHOLD_RELAXED_NAC_DIST - min_nuc_dist, 0, None) * _SCORE_W_DIST
                       + np.clip(max_ang - THRESHOLD_RELAXED_NAC_ANGLE, 0, None) * _SCORE_W_ANGLE)
        _penalty    = blockades * _SCORE_W_BLOCK + max(0.0, -wm_score) * _SCORE_W_WM
        score       = _productive * walden - _penalty

        if triad_ok and min_nuc_dist <= THRESHOLD_RELAXED_NAC_DIST and max_ang >= THRESHOLD_RELAXED_NAC_ANGLE:
            n_relaxed += 1
        if triad_ok and min_nuc_dist <= THRESHOLD_STRICT_NAC_DIST and max_ang >= THRESHOLD_STRICT_NAC_ANGLE:
            n_strict += 1
        if score > best_score and _productive > 0.0:
            best_score = score; ideal_frame_idx = f_idx
            ideal_geom = {'nuc_o': best_nuc_idx, 'lig_c': best_ca_idx}

        results.append({
            # Core NAC geometry
            "Frame":               f_idx,
            "NAC_Distance_A":      _r3(min_nuc_dist),
            "Tail_Cradle_Dist_A":  _r3(dist_tail),
            "NAC_Angle_Deg":       round(max_ang, 2),
            "Consensus":           round(score, 2),
            "Triad_NB":            _r3(nb_dist),
            "Triad_BA":            _r3(ba_dist),
            "Walden_TS_Flat":      int(walden_flat),
            "Pocket_Bound":        int(min_nuc_dist <= POCKET_RESIDENCY_DIST),
            "NAC_Geom_Pass":       int(min_nuc_dist <= THRESHOLD_RELAXED_NAC_DIST
                                       and max_ang >= THRESHOLD_RELAXED_NAC_ANGLE),
            "NAC_Strict_Pass":     int(min_nuc_dist <= THRESHOLD_STRICT_NAC_DIST
                                       and max_ang >= THRESHOLD_STRICT_NAC_ANGLE),
            "Triad_Intact":        int(triad_ok),
            "Catalytic_Pass":      int(triad_ok
                                       and min_nuc_dist <= THRESHOLD_RELAXED_NAC_DIST
                                       and max_ang >= THRESHOLD_RELAXED_NAC_ANGLE),
            # Dream Team distances to warhead carbon
            "DT_Nuc_LigC_A":    _r3(dt_dists.get('Nuc',    np.nan)),
            "DT_Clamp1_LigC_A": _r3(dt_dists.get('Clamp1', np.nan)),
            "DT_Clamp2_LigC_A": _r3(dt_dists.get('Clamp2', np.nan)),
            "DT_Acid_LigC_A":   _r3(dt_dists.get('Acid',   np.nan)),
            "DT_StabH_LigC_A":  _r3(dt_dists.get('Stab_H', np.nan)),
            "DT_StabW_LigC_A":  _r3(dt_dists.get('Stab_W', np.nan)),
            "DT_StabY_LigC_A":  _r3(dt_dists.get('Stab_Y', np.nan)),
            # Proton relay distances (residue–residue, not to ligand)
            "DT_Nuc_Base_A":    _r3(dt_nuc_base),
            "DT_Base_Acid_A":   _r3(dt_base_acid),
            # EAF ligand dynamics
            "EAF_MSA":          _r3(frame_msa),
            "EAF_RG":           _r3(frame_rg),
        })

    # ===============================================================================
    # Post-loop statistics
    # ===============================================================================
    total = len(results)
    # Viability denominators use n_pocket (bound frames), not total frames.
    # Catalytic_Viability_Pct = fraction of BOUND frames that are catalytically
    # competent — the meaningful metric for a pre-reactive ensemble.
    bound_rows = [r for r in results if r.get('Pocket_Bound', 0) == 1]

    def _mean_col(col: str, rows=None) -> float:
        src  = rows if rows is not None else results
        vals = [r[col] for r in src if col in r and r[col] == r[col]]  # NaN check
        return round(float(np.mean(vals)), 3) if vals else np.nan

    def _min_col(col: str, rows=None) -> float:
        src  = rows if rows is not None else results
        vals = [r[col] for r in src if col in r and r[col] == r[col]]  # NaN check
        return round(float(np.min(vals)), 3) if vals else np.nan

    stats = {
        "Job_Name":                 job_name,
        "Scientific_Rank":          rank,
        "Total_Frames":             total,
        "Frames_In_Pocket":         n_pocket,
        "Frames_NAC_Geom_Only":     n_nac,
        "Frames_Triad_Intact":      n_triad,
        "Frames_Relaxed_Catalysis": n_relaxed,
        "Frames_Strict_Catalysis":  n_strict,
        "Pocket_Retention_Pct":     round((n_pocket  / total)    * 100, 2) if total    else 0.0,
        "Catalytic_Viability_Pct":  round((n_relaxed / n_pocket) * 100, 2) if n_pocket else 0.0,
        "Strict_Viability_Pct":     round((n_strict  / n_pocket) * 100, 2) if n_pocket else 0.0,
        "Triad_Integrity_Pct":      round((n_triad   / n_pocket) * 100, 2) if n_pocket else 0.0,
        "NAC_Geom_Only_Pct":        round((n_nac     / n_pocket) * 100, 2) if n_pocket else 0.0,
        "MD_Avg_NAC_Dist_A":        _mean_col("NAC_Distance_A"),
        "MD_Min_NAC_Dist_A":        _min_col("NAC_Distance_A"),
        "MD_Avg_NAC_Angle_Deg":     _mean_col("NAC_Angle_Deg"),
        "MD_Best_Score":            round(best_score, 2),
        # Dream Team mean distances (pocket-bound frames only)
        "DT_Nuc_LigC_Mean_A":       _mean_col("DT_Nuc_LigC_A",    bound_rows),
        "DT_Clamp1_Mean_A":         _mean_col("DT_Clamp1_LigC_A", bound_rows),
        "DT_Clamp2_Mean_A":         _mean_col("DT_Clamp2_LigC_A", bound_rows),
        "DT_Acid_Mean_A":           _mean_col("DT_Acid_LigC_A",   bound_rows),
        "DT_StabH_Mean_A":          _mean_col("DT_StabH_LigC_A",  bound_rows),
        "DT_StabW_Mean_A":          _mean_col("DT_StabW_LigC_A",  bound_rows),
        "DT_StabY_Mean_A":          _mean_col("DT_StabY_LigC_A",  bound_rows),
        "DT_Nuc_Base_Mean_A":       _mean_col("DT_Nuc_Base_A",    bound_rows),
        "DT_Base_Acid_Mean_A":      _mean_col("DT_Base_Acid_A",   bound_rows),
        "Mean_Triad_NB_A":          _mean_col("Triad_NB"),
        "Mean_Triad_BA_A":          _mean_col("Triad_BA"),
        "Dream_Team_Mapped":        n_mapped,
        # WaterMap CSV thermodynamic statistics
        **compute_wm_csv_stats(wm_csv_data),
        # EAF ligand dynamics summary
        "EAF_MSA_Mean":     round(float(np.mean(eaf_msa)),  2) if len(eaf_msa) > 0 else np.nan,
        "EAF_MSA_Min":      round(float(np.min(eaf_msa)),   2) if len(eaf_msa) > 0 else np.nan,
        "EAF_RG_Mean":      round(float(np.mean(eaf_rg)),   3) if len(eaf_rg)  > 0 else np.nan,
        "EAF_N_Frames":     int(len(eaf_msa)),
        "EAF_File":         eaf_path.name if eaf_path else "Not found",
    }
    for k in ['degrader_tier', 'Dist_Nucleophile_ASP110', 'SN2_Attack_Angle',
              'Boltz_Model_Confidence']:
        if k in row: stats[k] = row[k]

    _n_mapped_log = stats.get('Dream_Team_Mapped', 0)

    # ── Output: per-frame CSV with rolling EAF smoothing ──────────────────────
    df_res = pd.DataFrame(results)
    if 'EAF_MSA' in df_res.columns and df_res['EAF_MSA'].notna().any():
        df_res['EAF_MSA_Smooth'] = df_res['EAF_MSA'].rolling(window=50, min_periods=1).mean()
    df_res.to_csv(job_out_dir / f"{job_name}_NAC_Data.csv", index=False)

    generate_individual_dashboard(
        df_res, job_name, job_out_dir / f"{job_name}_NAC_Dashboard.png", stats)

    if ideal_frame_idx != -1:
        topo.update_cms(cms_model, tr[ideal_frame_idx])
        mae_path = job_out_dir / f"{job_name}_Ideal_Final.maegz"
        cms_model.write(str(mae_path))
        inp_path = generate_qsite_inputs(mae_path, job_name, nuc_num, stab_f_num,
                                         ideal_geom['lig_c'], ideal_geom['nuc_o'],
                                         base_num=base_num, acid_num=acid_num)
        console_qmm_ready(f"Best frame: {ideal_frame_idx} | Score: {best_score:.2f} | "
                          f"QSite: {mae_path.name}")
        # ── Post-write sanity check: verify the artifacts THIS step writes ─────
        # (QSite itself is an external downstream Schrödinger/Jaguar run; this
        #  script only produces the .maegz best-frame structure and the .inp scan
        #  definition that feed it. Do NOT check for a QSite *output* here.)
        for _artifact in (mae_path, inp_path):
            if not _artifact.exists() or _artifact.stat().st_size == 0:
                console_info(
                    f"  [!] QSite input missing or empty for {job_name}: {_artifact.name}"
                )

    return stats


# ===============================================================================
# SECTION 9: MAIN EXECUTION
# ===============================================================================

def _resolve_work_dir(raw: str) -> Path:
    """
    Resolves a run dir argument to the 7_Physics_Validation path.

    Handles three input forms:
      1. Exact path to 7_Physics_Validation (already correct).
      2. Exact Boltz-2_Run_* dir (appends 7_Physics_Validation).
      3. Partial prefix such as 'Boltz-2_Run_20260309T085406Z' — glob-expanded to
         match the full timestamped directory on disk.
    """
    p = Path(raw)
    if not p.is_absolute():
        p = Path.cwd() / p

    # Case 1: already points into Physics_Validation
    if "Physics_Validation" in p.parts or p.name == "7_Physics_Validation":
        p.mkdir(parents=True, exist_ok=True)
        return p.resolve()

    # Case 2/3: p is or looks like a Boltz-2_Run_* dir (possibly a prefix)
    if p.exists() and p.is_dir():
        pv = p / "7_Physics_Validation"
        pv.mkdir(parents=True, exist_ok=True)
        return pv.resolve()

    # Partial prefix — glob parent for matching timestamped directories
    parent = p.parent
    stem   = p.name
    matches = sorted(parent.glob(f"{stem}*"),
                     key=lambda x: x.stat().st_mtime)
    if matches:
        pv = matches[-1] / "7_Physics_Validation"
        pv.mkdir(parents=True, exist_ok=True)
        return pv.resolve()

    pv = p / "7_Physics_Validation"
    pv.mkdir(parents=True, exist_ok=True)
    return pv.resolve()


def main():
    parser = argparse.ArgumentParser(
        description="PFAS-27 MD Thermodynamics & QM/MM Engine")
    parser.add_argument("run_dir",  nargs="?", default=None,
                        help="Boltz-2 run folder name or prefix (e.g. Boltz-2_Run_20260309T085406Z)")
    parser.add_argument("--dir",    default=None,
                        help="Alternative to positional run_dir (legacy --dir flag)")
    parser.add_argument("--lig",    default="LIG", help="Ligand residue name")
    parser.add_argument("--stride", type=int, default=1,
                        help="Frame stride (1 = all frames, default).")
    parser.add_argument("--ranks",  type=int, default=None,
                        help="Number of ranked jobs (default: auto-detect from MolecularDynamics dirs)")
    parser.add_argument("--nuc",    type=int, default=DREAM_TEAM_REF.get('Nuc', 110),  help="Fallback nucleophile resnum")
    parser.add_argument("--base",   type=int, default=DREAM_TEAM_REF.get('Base', 277), help="Fallback base resnum")
    parser.add_argument("--acid",   type=int, default=DREAM_TEAM_REF.get('Acid', 134), help="Fallback acid resnum")
    parser.add_argument("--csv",    default=None,
                        help="Path to master CSV (auto-detected if omitted)")
    args = parser.parse_args()

    raw_dir  = args.run_dir or args.dir or "."
    work_dir = _resolve_work_dir(raw_dir)

    # Auto-detect all available MD rank indices — scan MD, WaterMaps, and
    # 8_MD_Thermodynamics_Results so that any previously processed rank is included.
    _auto_rank_list: list[int] = []
    if args.ranks is None:
        _found_ranks: set[int] = set()
        for _scan_root in [
            work_dir / "MolecularDynamics",
            work_dir / "WaterMaps",
            work_dir.parent / "8_MD_Thermodynamics_Results",
        ]:
            if _scan_root.exists():
                for _d in _scan_root.iterdir():
                    if _d.is_dir():
                        _m = re.search(r'[Rr]ank[_\s]?(\d+)', _d.name)
                        if _m:
                            _found_ranks.add(int(_m.group(1)))
        if _found_ranks:
            _auto_rank_list = sorted(_found_ranks)
            args.ranks = max(_found_ranks)
        else:
            args.ranks = 5

    master_out_dir = work_dir.parent / "8_MD_Thermodynamics_Results"
    master_out_dir.mkdir(parents=True, exist_ok=True)

    global logger
    logger = (_setup_logging(master_out_dir / "11_MD_Thermodynamics_Engine.log",
                             "08_md_thermo_engine")
              if _setup_logging else None)


    _utils_mod.print_script_banner(
        "08_MD_Thermodynamics_QMMM_Engine_FAcDs.py",
        "MD Thermodynamics  ·  QM/MM Frame Extraction  ·  NAC Validation",
    )
    console_info(f"Run Directory    : {work_dir.parent}")
    console_info(f"Physics Validation : {work_dir}")
    console_info(f"MD Simulations   : {work_dir / 'MolecularDynamics'}")
    console_info(f"WaterMaps        : {work_dir / 'WaterMaps'}")
    console_info(f"Output           : {master_out_dir}")
    console_info(f"Ligand Resname   : {args.lig}")
    console_info(f"Frame Stride     : {args.stride} (requested){' — all frames' if args.stride == 1 else f' — 1-in-{args.stride} sampled'}")
    console_separator()

    # ── Ranked CSV (7_Boltz2_FAcDs_Ranked_*.csv or any *_Ranked*.csv) ──────────
    prod_dir     = work_dir.parent / "1_Boltz2_Production"
    ranked_csvs  = (list(prod_dir.glob("7_Boltz2_FAcDs_Ranked_*.csv")) or
                    list(prod_dir.glob("*_Ranked*.csv")))
    if not ranked_csvs:
        console_info(f"{ConsoleColours.FAIL}Error: No ranked CSV found in {prod_dir}{ConsoleColours.ENDC}")
        sys.exit(1)
    ranked_csv_path = sorted(ranked_csvs, key=lambda x: x.stat().st_mtime)[-1]
    console_info(f"Ranked CSV       : {ranked_csv_path.name}")
    df_ranked = pd.read_csv(ranked_csv_path, low_memory=False)

    # ── Ranked CSV for triad mapping ──────────────────────────────────────────
    triad_csv = Path(args.csv) if args.csv else ranked_csv_path
    triad_map = {}
    if triad_csv and triad_csv.exists():
        console_info(f"Ranked CSV       : {triad_csv.name}")
        triad_map = load_triad_mapping(triad_csv)
        console_info(f"Triad Entries    : {len(triad_map)} jobs mapped")
    else:
        console_info(f"{ConsoleColours.WARNING}Ranked CSV       : Not found — using alignment map / fallback defaults.{ConsoleColours.ENDC}")

    # Build rank → triad_override lookup
    rank_to_override = {}
    for _, csv_row in df_ranked.iterrows():
        m = re.search(r"(\d{7})", str(csv_row.get('job_name', '')))
        if m:
            job_idx = m.group(1).zfill(7)
            if job_idx in triad_map:
                rank_val = int(csv_row.get('Scientific_Rank', 0))
                rank_to_override[rank_val] = triad_map[job_idx]

    console_separator()

    master_stats   = []
    _stats_lock    = threading.Lock()
    _rank_list     = _auto_rank_list if _auto_rank_list else list(range(1, args.ranks + 1))
    _n_workers     = min(len(_rank_list), max(1, (os.cpu_count() or 4) - 2))

    console_info(f"  Effective Stride : {args.stride}{' (all frames)' if args.stride == 1 else f' (1-in-{args.stride})'}")
    console_separator()

    def _process_one_rank(r):
        """Process a single rank; called from worker thread.

        Returns (r, _log_lines, res) so the main thread can print all ranks
        in sorted order — preventing interleaved output when ranks run in
        parallel.  The header line is always the first entry in _log_lines.
        """
        _log_lines: list = []
        _log_lines.append(f"\n{SEPARATOR_LIGHT}")
        _log_lines.append(
            f"  Rank {r} / {len(_rank_list)}  —  MD Trajectory Analysis  [Rank_{r}]"
        )

        triad_override = rank_to_override.get(r)
        if triad_override:
            _log_lines.append(
                f"  [Rank {r}] Triad CSV linked: "
                f"Nuc[{triad_override.get('nuc')}] | "
                f"Base[{triad_override.get('base')}] | "
                f"Acid[{triad_override.get('acid')}]"
            )
        res = process_single_job(r, work_dir, df_ranked, master_out_dir,
                                 args.lig, args.stride,
                                 triad_override=triad_override,
                                 fallback_nuc=args.nuc,
                                 fallback_base=args.base,
                                 fallback_acid=args.acid)
        if res:
            if triad_override and 'static_data' in triad_override:
                for k, v in triad_override['static_data'].items():
                    if k not in res:
                        res[k] = v
            if "CONTROL" in res.get('Job_Name', '').upper():
                res['Category'] = 'Control'
            with _stats_lock:
                master_stats.append(res)
            # ── Per-rank result summary box ────────────────────────────────
            _pocket_pct    = res.get('Pocket_Retention_Pct', 0.0)
            _n_pocket      = res.get('Frames_In_Pocket', 0)
            _total_f       = res.get('Total_Frames', 0)
            _viab_pct      = res.get('Catalytic_Viability_Pct', 0.0)
            _n_relaxed     = res.get('Frames_Relaxed_Catalysis', 0)
            _strict_pct    = res.get('Strict_Viability_Pct', 0.0)
            _n_strict      = res.get('Frames_Strict_Catalysis', 0)
            _nac_geom_pct  = res.get('NAC_Geom_Only_Pct', 0.0)
            _n_nac_geom    = res.get('Frames_NAC_Geom_Only', 0)
            _triad_pct     = res.get('Triad_Integrity_Pct', 0.0)
            _n_triad_f     = res.get('Frames_Triad_Intact', 0)
            _dt_mapped     = res.get('Dream_Team_Mapped', 0)
            _wm_n          = res.get('WM_N_Sites', 0)
            _wm_stable     = res.get('WM_N_Stable', 0)
            _wm_dg         = res.get('WM_Mean_dG', float('nan'))
            _best_score    = res.get('MD_Best_Score', 0.0)
            _nb_mean       = res.get('Mean_Triad_NB_A', float('nan'))
            _ba_mean       = res.get('Mean_Triad_BA_A', float('nan'))
            _job_lbl       = res.get('Job_Name', f'Rank {r}')
            try:
                _wm_dg_str = f"{_wm_dg:.2f} kcal/mol" if _wm_dg == _wm_dg else "—"
            except Exception:
                _wm_dg_str = "—"
            _nb_str = f"{_nb_mean:.2f} Å" if _nb_mean == _nb_mean else "n/a"
            _ba_str = f"{_ba_mean:.2f} Å" if _ba_mean == _ba_mean else "n/a"

            # Three-tier outcome: full catalytic / geometry-only / fail
            if _viab_pct >= CFG.VIABILITY_PASS_THRESHOLD:
                _nac_status = "NAC PASS"
                _nac_col    = ConsoleColours.OKGREEN
            elif _nac_geom_pct >= CFG.VIABILITY_PASS_THRESHOLD:
                _nac_status = "GEOM PASS  (triad weak — verify His distance)"
                _nac_col    = ConsoleColours.WARNING
            else:
                _nac_status = "NAC FAIL"
                _nac_col    = ConsoleColours.FAIL

            _srows = [
                ("Pocket Retention",           f"{_pocket_pct:.2f}%   ({_n_pocket:,} / {_total_f:,} frames)"),
                ("Triad Integrity",            f"{_triad_pct:.2f}%   ({_n_triad_f:,} / {_n_pocket:,})  mean Nuc–Base: {_nb_str}"),
                ("NAC Geom Only  (no triad)",  f"{_nac_geom_pct:.2f}%   ({_n_nac_geom:,} / {_n_pocket:,})"),
                ("Relaxed Viability (triad+)", f"{_viab_pct:.2f}%   ({_n_relaxed:,} / {_n_pocket:,})"),
                ("Strict  Viability (triad+)", f"{_strict_pct:.2f}%   ({_n_strict:,} / {_n_pocket:,})"),
                ("Dream Team",                 f"{_dt_mapped} / {len(DREAM_TEAM_REF)} residues mapped"),
                ("WaterMap",                   f"{_wm_n} sites  |  {_wm_stable} stable  |  {_wm_dg_str}"),
                ("Best Frame Score",           f"{_best_score:.2f}" if _best_score > -1e8 else "—"),
                ("NAC Outcome",                _nac_status),
            ]
            # Column widths — cap total box at 124 chars to stay within 140-col terminals.
            # Key column is fixed; value column truncated with ellipsis if needed.
            _MAX_BOX = 124         # total including "  " indent + 2 border chars
            _sk = max(len(row[0]) for row in _srows)
            _sv_raw = max(len(row[1]) for row in _srows)
            _sv = min(_sv_raw, _MAX_BOX - _sk - 10)  # 10 = 2×│ + 4×space each side
            # Truncate long values so they fit
            def _trunc(s, w):
                return (s[:w - 1] + '…') if len(s) > w else s
            _title = f"  Rank {r}  ·  {_job_lbl}"[:_sk + _sv + 7]
            _log_lines.append(f"  ┌{'─' * (_sk + _sv + 9)}┐")
            _log_lines.append(f"  │ {_title:<{_sk + _sv + 7}} │")
            _log_lines.append(f"  ├{'─'*(_sk+4)}┬{'─'*(_sv+4)}┤")
            for _k, _v in _srows:
                _vt = _trunc(_v, _sv)
                if _k == "NAC Outcome":
                    # Colour applied in terminal; log strips ANSI for file output
                    _vd = f"{_nac_col}{_vt}{ConsoleColours.ENDC}{' ' * (_sv - len(_vt))}"
                else:
                    _vd = f"{_vt:<{_sv}}"
                _log_lines.append(f"  │  {_k:<{_sk}}  │  {_vd}  │")
            _log_lines.append(f"  └{'─'*(_sk+4)}┴{'─'*(_sv+4)}┘")

            # Diagnostic note when triad gates all frames
            if _n_nac_geom > 0 and _n_relaxed == 0:
                _log_lines.append(
                    f"  ↳ Triad absent in all NAC frames — "
                    f"mean Nuc–Base: {_nb_str} (threshold ≤ {CFG.THRESHOLD_TRIAD_NB_MD:.1f} Å)  |  "
                    f"mean Base–Acid: {_ba_str} (threshold ≤ {CFG.THRESHOLD_TRIAD_BA_MD:.1f} Å)."
                )
                _log_lines.append(
                    f"    SN2 geometry present in {_n_nac_geom:,} frame(s). "
                    f"Consider checking His residue placement."
                )
        return r, _log_lines, res

    console_info(f"Parallel workers : {_n_workers} (of {len(_rank_list)} ranks)")
    print(flush=True)

    _failed_ranks: list  = []
    _all_results:  list  = []
    with ThreadPoolExecutor(max_workers=_n_workers,
                            thread_name_prefix="Rank") as _pool:
        _futures = {_pool.submit(_process_one_rank, r): r for r in _rank_list}
        for _fut in as_completed(_futures):
            try:
                _rank_num, _rank_log, _rank_res = _fut.result()
                _all_results.append((_rank_num, _rank_log, _rank_res))
            except Exception as _exc:
                _r = _futures[_fut]
                _failed_ranks.append(_r)
                console_info(f"{ConsoleColours.FAIL}[!] Rank {_r} failed: {_exc}{ConsoleColours.ENDC}")
                import traceback as _tb
                _tb.print_exc()

    # Print all rank output blocks in ascending rank order to prevent
    # interleaved lines from parallel workers.
    for _rank_num, _rank_log, _rank_res in sorted(_all_results, key=lambda x: x[0]):
        for _line in _rank_log:
            print(_line, flush=True)

    if _failed_ranks:
        console_info(f"{ConsoleColours.WARNING}[!] {len(_failed_ranks)} rank(s) raised exceptions: "
                     f"{sorted(_failed_ranks)}{ConsoleColours.ENDC}")

    if not master_stats:
        console_info(f"{ConsoleColours.FAIL}[CRITICAL] No MD data collected — "
                     f"all ranks failed or no trajectory data found.{ConsoleColours.ENDC}")
        sys.exit(1)

    if master_stats:
        console_separator()
        console_title("Finalising Global MD Rankings")

        df_master = pd.DataFrame(master_stats)
        df_master = df_master.sort_values(by="Catalytic_Viability_Pct", ascending=False)
        df_master.insert(0, "Dynamic_Rank", range(1, len(df_master) + 1))

        master_csv_path = master_out_dir / "08_MD_Master_Ranking.csv"
        df_master.to_csv(master_csv_path, index=False)
        console_info(f"Total Simulations Validated : {len(df_master)}")
        console_info(f"Master Ranking Sheet Saved  : {master_csv_path.resolve()}")
        console_separator()

        if _rcon:
            tbl = _RichTable(title="MD Thermodynamics Ranking Summary",
                             show_header=True, header_style="bold bright_blue",
                             box=None, show_lines=True)
            for col in ("Rank", "Job", "Pocket%", "Viability%", "Strict%",
                        "Triad%", "DT Mapped", "WM Stable", "WM dG", "Avg Dist Å"):
                tbl.add_column(col, justify="right" if col not in ("Job",) else "left")
            for _, row in df_master.iterrows():
                viab   = row.get("Catalytic_Viability_Pct", 0.0)
                colour = "green" if viab >= CFG.VIABILITY_HIGH_THRESHOLD else ("yellow" if viab >= CFG.VIABILITY_PASS_THRESHOLD else "red")
                tbl.add_row(
                    str(int(row["Dynamic_Rank"])),
                    format_job_label(row["Job_Name"], row["Dynamic_Rank"]),
                    f"{row.get('Pocket_Retention_Pct', 0):.1f}",
                    f"[{colour}]{viab:.1f}[/{colour}]",
                    f"{row.get('Strict_Viability_Pct', 0):.1f}",
                    f"{row.get('Triad_Integrity_Pct', 0):.1f}",
                    str(row.get('Dream_Team_Mapped', '-')),
                    str(row.get('WM_N_Stable', '-')),
                    f"{row.get('WM_Mean_dG', float('nan')):.2f}",
                    f"{row.get('MD_Avg_NAC_Dist_A', 0):.2f}",
                )
            _rcon.print(tbl)

        for label, fn in [
            ("Global Comparative Dashboard (3-panel)...", generate_global_comparative_dashboard),
            ("MD Viability & Retention Bar Chart...", generate_viability_bar_chart),
        ]:
            console_info(f"  ✔ Generating {label}")
            try:
                fn(master_out_dir, df_master)
            except Exception as e:
                console_info(f"  [!] {label.split('(')[0].strip()} failed: {e}")


if __name__ == "__main__":
    import time as _time
    _t0 = _time.perf_counter()
    main()
    _utils_mod.print_elapsed(_t0, "08_MD_Thermodynamics_QMMM_Engine_FAcDs.py")
