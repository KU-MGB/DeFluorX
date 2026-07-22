#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 07  |  MD + QM/MM Defluorination Engine
===============================================================================

Terminal computational step: turns the Desmond MD trajectories into a concrete,
ranked verdict on whether each candidate DEFLUORINATES — not merely binds.

It combines four evidence streams into the master ranking + figures:
  1. Kinetic pre-organisation — near-attack-conformation (NAC) geometry, the
     8-residue Dream Team catalytic machinery, WaterMap hydration, Desmond EAF
     ligand dynamics, and the CONTINUOUS strict-NAC dwell time (ns).
  2. Quantum barrier — QSite QM/MM SN2 relaxed scans over the top pre-organised
     frames, PARSED into ΔE‡ / ΔE_rxn and a departing-fluoride charge (→ −1 = F⁻);
     the reaction-profile figure is the direct proof of C–F cleavage.
  3. Reactive-pose thermodynamics — Prime MM-GBSA (Step 06) conditioned on the
     strict-NAC frames, with an energy-component decomposition + machinery
     engagement (classical support, not bond-breaking).
  4. Verdict — Is_Defluorinating gate + Defluor_Propensity (P(strict-NAC)·
     exp(−ΔE‡/RT)), and the whole-story Defluorination Landscape figure.

All thresholds, gate cut-offs, and figure colours come from CFG (SSOT).

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 20 July 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 07_MD_QMMM_Defluorination_FAcDs.py
  Role          : Trajectory analysis engine; terminal computational step before
                  QM/MM (outputs ideal frame + QSite .inp files).
  Imports from  : 00_01_Project_Config_FAcDs.py  (CFG — all thresholds + tier metadata)
                  00_02_Project_Utils_FAcDs.py   (ConsoleColours, geometric utilities)
  Reads         : <Run>/6_Physics_Validation/05_MD_Simulations/desmond_md_job_R_N/*-out.cms
                                                               /*_trj/   (dir carrying the _R_N rank token; *Rank_N* also matched)
                                                               /*.eaf
                  <Run>/6_Physics_Validation/03_WaterMaps/watermap_R_N.csv  (Step-06 WaterMap export)
                  <Run>/6_Physics_Validation/03_WaterMaps/watermap_R_N/*_wm.maegz
                  <Run>/1_Boltz2_Production/6_Boltz2_FAcDs_Ranked_*.csv
                  <Run>/1_Boltz2_Production/5_Boltz2_FAcDs_Master_*.csv
  Writes        : <Run>/7_MD_Thermodynamics_Results/Rank_N_<Name>/
                    - <Name>_NAC_Data.csv          (per-frame geometry + DT)
                    - <Name>_NAC_Dashboard.png     (2-panel figure)
                    - <Name>_Ideal_Final.maegz      (best frame for QSite)
                    - <Name>_QSite_SN2/            (primary QM/MM scan) + _QSite_SN2_f<frame>/ (extra ensemble frames)
                    - <Name>_QSite_Reaction_Profile.png   (PES vs reaction coordinate + departing-F charge → the C–F-cleavage proof)
                    - <Name>_MMGBSA_NAC_Decomposition.png (ΔG components: whole trajectory vs the reactive pose)
                    - <Name>_Machinery_Engagement.png     (per-residue distance to the warhead C + contact occupancy)
                  <Run>/7_MD_Thermodynamics_Results/01_MD_Master_Ranking.csv
                    (adds NAC dwell in ns, parsed QM/MM ΔE‡ / ΔE_rxn, departing-F
                     charge, NAC-conditioned MM-GBSA + component decomposition, and
                     the Defluor_Propensity / Is_Defluorinating verdict)
                  <Run>/7_MD_Thermodynamics_Results/05_Defluorination_Landscape.png
                    (whole-story figure: persistence × QM/MM barrier × binding)
                  <Run>/7_MD_Thermodynamics_Results/06_MMGBSA_Decomposition_AllRanks.png
                  <Run>/7_MD_Thermodynamics_Results/07_Machinery_Engagement_AllRanks.png
                    (the same two reactive-pose figures, merged across candidates)
  Upstream      : 06_Physics_Validation_FAcDs.py → runs WaterMap · System Builder · MD · SID · MM-GBSA;
                                                    produces the MD trajectories, WaterMap CSVs,
                                                    *_SID-out.eaf + Prime MM-GBSA summary consumed here
                  05_TopN_and_PDB_Preparation_FAcDs.py → provides ranked structures & IDs
                  02_Production_FAcDs.py         → master CSV with alignment maps
  Downstream    : None (terminal step; QSite .inp feeds Schrödinger QSite/Jaguar)
────────────────────────────────────────────────────────────────────────────────

# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Schrödinger Python: MUST be executed using `run python3` within a
     Schrödinger environment (requires `schrodinger.analysis.topology`).
  2. Trajectory Compatibility: Strictly expects Desmond CMS/trj format; will
     not process GROMACS or AMBER trajectories without prior conversion.
  3. Memory & I/O: Reading large trajectories (1000+ frames) at stride 1 is
     extremely I/O intensive; recommend high-speed NVMe or local scratch disk.
  4. WaterMap Input: Reads the WaterMap CSV/MAE that Step 06 produces under
     03_WaterMaps/ (watermap_R_N.csv + watermap_R_N/*_wm.maegz).
───────────────────────────────────────────────────────────────────────────────

Usage:
    python 07_MD_QMMM_Defluorination_FAcDs.py Boltz-2_Run_20260309T085406Z
    python 07_MD_QMMM_Defluorination_FAcDs.py Boltz-2_Run_20260309T085406Z --stride 5 --ranks 3 --lig PFAS
    python 07_MD_QMMM_Defluorination_FAcDs.py Boltz-2_Run_20260309T085406Z --nuc 85 --base 250 --acid 112

Arguments:
    run_dir           Positional. Boltz-2 run folder name or prefix (e.g.
                      Boltz-2_Run_20260309T085406Z). Glob-expanded to match timestamp
                      suffix automatically. Resolves into 6_Physics_Validation.
    --dir    DIR      Alternative to positional (legacy). Same resolution logic.
    --lig    RESNAME  Ligand residue name in the CMS system file.  Default: LIG
    --stride N        Analyse every N-th trajectory frame (1 = all frames).
                      Higher values trade accuracy for speed.     Default: 1
    --ranks  N        Number of top-ranked MD jobs to process (Rank_1 … Rank_N).
                      Default: auto-detected from directories in 05_MD_Simulations/
    --nuc    RESNUM   Fallback nucleophile residue number if alignment map absent.
                      Default: 110  (FAcD canonical Asp110)
    --base   RESNUM   Fallback catalytic base residue number.
                      Default: 280  (FAcD canonical His280, 3R3U numbering)
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
  8. QSite automation: B3LYP/6-31+G(d,p) QM/MM SN2 relaxed-scan generation +
     execution over the top CFG.QSITE_N_FRAMES pre-organised NAC frames (ensemble
     barrier, not a single-frame lower bound), with the scan output PARSED back
     into ΔE‡ (min/mean/σ) and ΔE_rxn.
  9. 3D Smart-Lock: geometry-biased triad & fluorine-cradle detection.
 10. Rich progress bars and colour-coded PASS/FAIL NAC reporting.
 11. Master aggregation: 01_MD_Master_Ranking.csv.
 12. NAC persistence: longest/mean continuous strict-NAC dwell converted to ns
     (real "time in position", not a frame-count fraction).
 13. NAC-conditioned MM-GBSA: ΔG_bind over the strict-NAC frames vs the global
     mean, PLUS an energy-component decomposition (Coulomb / vdW / Covalent strain
     / …) → *_MMGBSA_NAC_Decomposition.png, and the catalytic machinery's
     engagement — every Dream-Team residue's median distance to the warhead carbon
     with the share of reactive frames in which the contact holds, judged against
     the CFG criterion bands → *_Machinery_Engagement.png. Both are also written
     merged across candidates (13_* and 14_*).
 14. QM/MM reaction profile: the PES along the SN2 coordinate with ΔE‡ / ΔE_rxn
     and the departing-fluoride Mulliken charge (→ −1 = F⁻) — the direct proof of
     C–F cleavage → *_QSite_Reaction_Profile.png (cols F_Charge_Reactant/Product/Delta).
 15. Defluorination verdict: Is_Defluorinating gate + Defluor_Propensity
     ( P(strict-NAC)·exp(−ΔE‡/RT) ) — the concrete turnover claim, not affinity.
 16. Defluorination landscape figure (12_*): persistence × QM/MM barrier × binding.
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
  QSite QM/MM level of theory (B3LYP/6-31+G(d,p)):
    Becke, A.D. (1993) J Chem Phys 98:5648–5652. DOI: https://doi.org/10.1063/1.464913
    Lee, C., Yang, W. & Parr, R.G. (1988) Phys Rev B 37:785–789.
    DOI: https://doi.org/10.1103/PhysRevB.37.785
    (B3LYP is required: QSite frozen-orbital QM/MM cuts reject meta-GGA hybrids
     such as M06-2X and dispersion-corrected variants such as B3LYP-D3.)
    Rosta, E., Klähn, M. & Warshel, A. (2006) J Phys Chem B 110:2934–2941.
    DOI: https://doi.org/10.1021/jp057109j
    Murphy, R.B. et al. (2000) J Comput Chem 21:1442–1457.
    DOI: https://doi.org/10.1002/1096-987X(200012)21:16<1442::AID-JCC3>3.0.CO;2-O
"""

import sys
import os
import shutil

# =============================================================================
# SCHRÖDINGER BOOTSTRAP
# =============================================================================
# Auto-sets SCHRODINGER=/opt/schrodinger if the env var is absent.
# When invoked with plain `python`, re-invokes transparently via
# $SCHRODINGER/run so the Schrödinger Python interpreter is used.
# Both forms are equivalent:
#   python 07_MD_QMMM_Defluorination_FAcDs.py Boltz-2_Run_20260309T085406Z
#   $SCHRODINGER/run 07_MD_QMMM_Defluorination_FAcDs.py Boltz-2_Run_20260309T085406Z
import subprocess as _sp

if "SCHRODINGER" not in os.environ:
    os.environ["SCHRODINGER"] = "/opt/schrodinger"

try:
    from schrodinger.application.desmond.packages import traj, topo
    from schrodinger import structure
except ImportError:
    _run_exe = os.path.join(os.environ["SCHRODINGER"], "run")
    if os.path.isfile(_run_exe):
        sys.exit(_sp.call([_run_exe, os.path.abspath(__file__)] + sys.argv[1:]))
    print(f"CRITICAL: Schrödinger suite not found at {os.environ['SCHRODINGER']}.\n"
          f"Set the SCHRODINGER environment variable to your installation path.")
    sys.exit(1)
# =============================================================================

import re
import time
import argparse
import warnings
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

# CPU usage cap (total cores − 2; mirrors CFG.PREP_CPU_RESERVE). Reserve 2 cores
# for OS/desktop stability by limiting the thread-pool maths libraries (BLAS /
# MKL / OpenMP / NumExpr). Must precede numpy import to take effect;
# setdefault() preserves any value exported by the caller or pipeline runner.
_CPU_CAP = str(max(1, (os.cpu_count() or 4) - 2))
for _tv in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_tv, _CPU_CAP)

import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator

try:
    from rich.console import Console as _RichConsole
    from rich.table import Table as _RichTable
    _RICH_AVAILABLE = True
    _rcon = _RichConsole()
except ImportError:
    _RICH_AVAILABLE = False
    _rcon = None

import importlib.util as _ilu


class LazyTrajectory:
    """Memory-efficient lazy trajectory wrapper that concatenates multiple segment paths
    without loading all frames or segment contents in memory at once.
    """
    def __init__(self, trj_paths):
        self.trj_paths = list(trj_paths)
        self._readers = [traj.read_traj(str(p)) for p in self.trj_paths]
        self._lengths = [len(r) for r in self._readers]
        self._total_frames = sum(self._lengths)

    def __len__(self):
        return self._total_frames

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            raise NotImplementedError("Slices not supported")
        curr = idx
        if curr < 0:
            curr += self._total_frames
        if curr < 0 or curr >= self._total_frames:
            raise IndexError(f"Index {idx} out of range")
        for r, l in zip(self._readers, self._lengths):
            if curr < l:
                return r[curr]
            curr -= l
        raise IndexError(f"Index {idx} out of range")


def _load_module(name: str, path: Path):
    """Load a Python file as a module regardless of its filename."""
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_REPO_DIR  = Path(__file__).resolve().parent
_cfg_mod   = _load_module("ProjectConfig", _REPO_DIR / "00_01_Project_Config_FAcDs.py")
_utils_mod = _load_module("ProjectUtils",  _REPO_DIR / "00_02_Project_Utils_FAcDs.py")
CFG        = _cfg_mod.CFG()

# The 3R3U × FA positive control gets one distinct label + colour across every 07 figure, so it
# never blends into the candidate fluoroacetate cohort.
_CTRL_LABEL  = "3R3U-FA"
_CTRL_COLOUR = CFG.VIS_ACCENT["control"]

# ConsoleColours sourced from 00_02_Project_Utils (single canonical definition).
apply_figure_style = _utils_mod.apply_figure_style
write_json_atomic  = _utils_mod.write_json_atomic
auto_label_colour  = _utils_mod.auto_label_colour
ConsoleColours  = _utils_mod.ConsoleColours
SEPARATOR_HEAVY = _utils_mod.SEPARATOR_HEAVY
SEPARATOR_LIGHT = _utils_mod.SEPARATOR_LIGHT

# Console helpers and logging setup sourced from 00_02_Project_Utils.
_console_title          = _utils_mod.console_title
_console_info           = _utils_mod.console_info
_console_sep            = _utils_mod.console_separator
_setup_logging          = _utils_mod.setup_logging
get_mic_vector          = _utils_mod.get_mic_vector
_mic_dists_2d           = _utils_mod.mic_dists_2d
clean_spines            = _utils_mod.clean_spines
_calc_improper_dihedral = _utils_mod.calculate_improper_dihedral
_ensure_box_3x3             = _utils_mod._ensure_box_3x3
"""
THE NUCLEOPHILE IS AN ASPARTATE. CFG SAYS SO; THIS SCRIPT NOW ASKS IT.

Both call sites hardcoded {'ASP', 'GLU', 'ASH', 'GLH'} — a set that admits GLUTAMATE, which
CFG.ROLE_EXPECTED_RESIDUES has never allowed and whose own comment forbids by name ("no Glu drift"). An
unmapped glutamate could therefore usurp the mapped aspartate through the soft distance bias and be
reported as the attacking residue.

ASH is retained because it is the same aspartate under a force-field name for its protonated form — a
question of NOMENCLATURE, not of identity. Whether that aspartate is actually deprotonated (and so able
to attack) is decided in Step 05 by enforce_catalytic_protonation, which strips the acidic proton. The
recognition set names the residue; the protonation policy makes it a nucleophile.
"""
_NUCLEOPHILE_RESIDUES = frozenset(CFG.ROLE_EXPECTED_RESIDUES["Nucleophile"])
_ACID_RESIDUES         = frozenset(CFG.ROLE_EXPECTED_RESIDUES["Acid_Catalyst"])
_BASE_RESIDUES         = frozenset(CFG.ROLE_EXPECTED_RESIDUES["Base_Catalyst"])

find_nucleophile_od_fallback = _utils_mod.find_nucleophile_od_fallback

PLOT_LOCK = threading.Lock()


# =============================================================================
# SECTION 1: GLOBAL CONSTANTS & CONFIGURATION
# =============================================================================
# All thresholds sourced from 00_01_Project_Config_FAcDs.py (CFG).
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


# =============================================================================
# SECTION 2: CONSOLE WRAPPER FUNCTIONS
# =============================================================================

import threading
thread_logger = threading.local()

_strip_ansi = getattr(_utils_mod, '_strip_ansi', lambda x: x)


def _sigmoid07(x: float, k: float, x0: float) -> float:
    """The graded NAC term used to pick the attacking oxygen (CFG §4.1).

    Identical in form to Step 02's sigmoid, so both stages score a near-attack conformation on the
    same scale and agree on which aspartate oxygen is attacking. A frame whose distance and angle
    are judged by one rule in the ranking and another in the trajectory analysis is not comparable
    with itself.
    """
    return float(1.0 / (1.0 + np.exp(-np.clip(k * (x - x0), -60.0, 60.0))))


def console_title(msg: str) -> None:
    if hasattr(thread_logger, 'lines'):
        thread_logger.lines.append(f"  [Rank {thread_logger.rank}] {msg}")
        if logger:
            logger.info(f"Rank {thread_logger.rank} | {_strip_ansi(msg)}")
    else:
        _console_title(msg, logger)


def console_info(msg: str) -> None:
    if hasattr(thread_logger, 'lines'):
        for line in str(msg).split('\n'):
            thread_logger.lines.append(f"  [Rank {thread_logger.rank}] {line}")
        if logger:
            logger.info(f"Rank {thread_logger.rank} | {_strip_ansi(msg)}")
    else:
        _console_info(msg, logger)


def _mask_oomd_at_start():
    """Prompt for the (optional) sudo password up front and mask systemd-oomd, so the memory-heavy
    QSite QM/MM phase near the end of the run is not OOM-killed. Prompting at the START lets the user
    walk away — the run stays unattended. Returns a restore callable (registered with atexit). If sudo
    is unavailable or skipped, the run proceeds unprotected (NORMAL mode). No-op when not on a TTY."""
    if not sys.stdin.isatty():
        return lambda: None
    console_info("OPTIONAL — protect this run from the Linux out-of-memory killer.")
    console_info("  QSite QM/MM holds large systems in memory for hours; systemd-oomd can kill it. "
                 "Masking needs root.")
    console_info("  Enter your sudo password to mask systemd-oomd, or press Enter / Ctrl-D to skip.")
    try:
        with open("/dev/tty") as _tty:
            _ok = _sp.run(["sudo", "-v"], stdin=_tty).returncode == 0
    except Exception:
        _ok = _sp.run(["sudo", "-v"]).returncode == 0 if True else False
    if not _ok:
        console_info("  [NORMAL MODE] sudo unavailable/skipped — systemd-oomd NOT masked.")
        return lambda: None
    _stop = threading.Event()

    def _keepalive():
        while not _stop.wait(60):
            _sp.run(["sudo", "-vn"], capture_output=True)
    threading.Thread(target=_keepalive, daemon=True).start()
    console_info("  Masking systemd-oomd — restored automatically when Step 07 exits.")
    _sp.run(["sudo", "systemctl", "stop", "systemd-oomd"], capture_output=True)
    _sp.run(["sudo", "systemctl", "mask", "systemd-oomd.socket"], capture_output=True)

    def _restore():
        _stop.set()
        console_info("Restoring systemd-oomd services...")
        _sp.run(["sudo", "systemctl", "unmask", "systemd-oomd.socket"], capture_output=True)
        _sp.run(["sudo", "systemctl", "start", "systemd-oomd"], capture_output=True)
    return _restore


def console_separator(heavy: bool = True) -> None:
    if hasattr(thread_logger, 'lines'):
        sep = SEPARATOR_HEAVY if heavy else SEPARATOR_LIGHT
        thread_logger.lines.append(f"  [Rank {thread_logger.rank}] {sep}")
    else:
        _console_sep(logger, heavy=heavy)


def _print_labeled(label: str, ansi_col: str, msg: str,
                   log_label: str | None = None, level: str = "info") -> None:
    """Print a prefixed status line with ANSI colour.
    Colour passes through pipeline tee to the terminal; file logger stays clean."""
    _log = log_label or label
    if hasattr(thread_logger, 'lines'):
        thread_logger.lines.append(f"  [Rank {thread_logger.rank}] {ansi_col}{label}\033[0m {msg}")
    else:
        if _rcon:
            _rcon.print(f"  [bold]{label}[/bold] {msg}", markup=False)
        else:
            print(f"  {ansi_col}{label}\033[0m {msg}", flush=True)
    if logger:
        r_prefix = f"Rank {thread_logger.rank} | " if hasattr(thread_logger, 'rank') else ""
        getattr(logger, level)(f"{r_prefix}{_log} | {msg}")


def console_qmm_ready(msg: str) -> None:
    _print_labeled("QM/MM READY", "\033[94m", msg)


# =============================================================================
# SECTION 3: GEOMETRY & MATHEMATICAL UTILITIES
# =============================================================================

# NOTE: This function is NOT equivalent to 00_02_Project_Utils.calculate_min_distance.
# It uses the Schrödinger frame API (frame.pos(idx)) rather than numpy arrays.
# The API divergence is intentional — required for Schrödinger/Maestro integration.
def calculate_min_distance(frame, indices_A: list, indices_B: list) -> float:
    """Minimum PBC-corrected distance between two atom index sets (vectorised).

    The O(N×M) pairwise MIC distances are computed in one numpy call via
    _mic_dists_2d instead of a Python double loop — the per-atom frame.pos()
    gather is unavoidable, but the distance maths is fully vectorised.
    """
    if not indices_A or not indices_B:
        return np.nan
    box = frame.box if hasattr(frame, 'box') else None
    pos_a = np.asarray([frame.pos(a) for a in indices_A], dtype=float)   # (nA, 3)
    pos_b = np.asarray([frame.pos(b) for b in indices_B], dtype=float)   # (nB, 3)
    return float(_mic_dists_2d(pos_a, pos_b, box).min())


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

    # Restrict the warhead to the scissile α-carbon — the carbon adjacent to the ligand
    # carboxylate head — mirroring Step 02's reactive-centre gating. Only C–F bonds ON
    # that α-carbon are candidate scissile bonds, so an internal CF2/CF3 of a polyfluoro
    # decoy (e.g. PFOA) is not modelled as the reaction centre. TFA's α-CF3 is retained
    # (its scissile carbon IS the α-carbon). Falls back to all C–F bonds, logged, when no
    # carboxylate/α-carbon can be identified (e.g. a non-carboxylate chemotype).
    def _carboxylate_alpha_carbon():
        """
        The α-carbon is the FLUORINATED carbon next to a ligand carboxylate. Every carboxylate is
        examined, not just the first one found: a dicarboxylic ligand has two heads, and returning
        the first α-carbon in atom order can lock the reaction centre onto the unfluorinated or
        solvent-exposed end. Candidates are ranked by how many C–F bonds they carry, so the
        defluorination site is the one the enzyme could actually act on; a carboxylate whose
        α-carbon bears no fluorine cannot be the scissile centre and is skipped.
        """
        _lig_idx = {int(i) for i in lig_atoms}
        _cf_carbons = {c for c, _f in cf_pairs}
        _cands = []
        for _i in lig_atoms:
            _a = cms_model.atom[_i]
            if _a.atomic_number != 6:
                continue
            _o_neigh = [b.atom2 for b in _a.bond if b.atom2.atomic_number == 8]
            if len(_o_neigh) < 2:                        # not a carboxylate carbon (C bonded to ≥2 O)
                continue
            for b in _a.bond:
                _nb = b.atom2
                if _nb.atomic_number == 6 and _nb.index in _lig_idx and _nb.index in _cf_carbons:
                    _cands.append(_nb.index)
        if not _cands:
            return None
        return max(_cands, key=lambda _c: sum(1 for _cc, _f in cf_pairs if _cc == _c))

    _alpha_idx = _carboxylate_alpha_carbon()
    if _alpha_idx is not None:
        _alpha_cf = [(c, f) for (c, f) in cf_pairs if c == _alpha_idx]
        if _alpha_cf:
            cf_pairs = _alpha_cf
        else:
            console_info("    [i] α-carbon carries no C–F bond; retaining all ligand C–F bonds.")
    else:
        console_info("    [i] No ligand carboxylate/α-carbon identified; retaining all ligand C–F bonds.")

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

    # 3. Nucleophile: closest ASP (aspartate Oδ) to ligand C, biased by mapped_nuc
    best_nuc_key = None
    min_eff      = float('inf')
    actual_dist  = float('inf')
    for res_key, data in res_dict.items():
        if not data['O_idx'] or data['ptype'] not in _NUCLEOPHILE_RESIDUES:
            continue
        d     = calculate_min_distance(frame_0, lig_c_idxs, data['O_idx'])
        bonus = CFG.SMART_LOCK_BIAS_DIST if (mapped_nuc and abs(data['resnum'] - mapped_nuc) <= CFG.SMART_LOCK_RESNUM_WINDOW) else 0.0
        if (d + bonus) < min_eff:
            min_eff = d + bonus; best_nuc_key = res_key; actual_dist = d

    if not best_nuc_key or actual_dist > CFG.SMART_LOCK_NUC_MAX_DIST:
        # ── Oδ Orientation Fallback ────────────────────────────────────────────
        # Primary geometry search exhausted.  Delegate to shared utility
        # find_nucleophile_od_fallback() (00_02_Project_Utils_FAcDs.py).
        # The function accepts plain NumPy arrays only; extract positions here
        # before calling so CMS atom-group objects never enter the utility.
        # Threshold aligns with CFG.NAC_ANGLE_RELAXED (BRAIN.md §7).
        c_idx, f_idx = cf_pairs[0]
        c_pos_np   = np.array(frame_0.pos(c_idx), dtype=float)
        f_pos_np   = np.array(frame_0.pos(f_idx), dtype=float)
        # Build candidate list: (resnum, od1_np, od2_np) for each ASP/ASH residue.
        _asp_cands = []
        for rk, rd in res_dict.items():
            if rd['ptype'] not in _ACID_RESIDUES:
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
        # CFG's set also carries the AMBER/CHARMM histidine names (HSD/HSE/HSP). The hardcoded set here
        # listed only the OPLS ones, so a base named HSD/HSE/HSP would have been invisible — the Smart-Lock
        # would have reported "no base found" on a structure that has one.
        if not data['N_idx'] or data['ptype'] not in _BASE_RESIDUES:
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

    # 5. Acid: closest ASP (aspartate Oδ) to base, biased by mapped_acid
    idx_acid = []
    if idx_base:
        best_acid_key = None; min_eff = float('inf')
        for res_key, data in res_dict.items():
            if (res_key == best_nuc_key or not data['O_idx']
                    or data['ptype'] not in _ACID_RESIDUES):
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

    # 6. Fluorine cradle: TRP/TYR/HIS heavy atoms within fluorine cradle radius of
    #    nucleophile. HIS is the third fluoride-stabilising H-bond donor (His155 in
    #    native FAcD). The catalytic base His (best_base_key) is excluded so it is
    #    not double-counted as both the general base and a cradle stabiliser.
    idx_cradle = []
    _fcr = CFG.F_CRADLE_RADIUS
    _cradle_types = {'TRP', 'TYR', 'HIS', 'HIP', 'HIE', 'HID'}
    for res_key, data in res_dict.items():
        if res_key == best_base_key:
            continue
        if data['ptype'] in _cradle_types:
            if calculate_min_distance(frame_0, idx_nuc, data['heavy_idx']) <= _fcr:
                idx_cradle.extend(data['heavy_idx'])

    console_info(f"    {ConsoleColours.OKBLUE}↳ 3D Smart-Lock: Fluorine Cradle → "
                 f"{len(idx_cradle)} TRP/TYR/HIS heavy atoms.{ConsoleColours.ENDC}")
    console_info(f"    {ConsoleColours.OKBLUE}↳ Polyfluorinated Engine: "
                 f"{len(cf_pairs)} C-F bonds tracked.{ConsoleColours.ENDC}")

    return idx_nuc, idx_base, idx_acid, idx_cradle, cf_pairs


# =============================================================================
# SECTION 4: DATA HANDLING & LOADING FUNCTIONS
# =============================================================================

def parse_mapping(map_str: str) -> dict:
    """Parses 'ASP110:ASP112 | HIS280:HIS282' alignment string into {ref_num: tgt_num} (3R3U reference numbering)."""
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


def identity_from_cms(folder: "Path", df_ranked) -> "str | None":
    """The MD job's true model job_name, read from the token embedded in its
    -out.cms, matched against df_ranked['job_name'].

    The MD folders (desmond_md_job_R_1, R_2, R_11 …) are whichever cases had
    MD-ready trajectories, so a folder's R-number need not equal its ranked
    position. The prepared model name is baked into every -out.cms (e.g.
    '0032129_1190_A0A2U3PT06_9BRAD_26_Fluoroacetate'), so read it and match the
    ranked row on identity rather than trusting R_N == Scientific_Rank. Returns
    None when no -out.cms / recognisable token is present (caller then falls back
    to the Scientific_Rank match)."""
    if 'job_name' not in getattr(df_ranked, 'columns', []):
        return None
    _cms = next(iter(sorted(folder.glob("*-out.cms"))), None) \
        or next(iter(sorted(folder.glob("**/*-out.cms"))), None)
    if _cms is None:
        return None
    try:
        _txt = _cms.read_text(errors="ignore")
    except Exception:
        return None
    _jn   = {str(x) for x in df_ranked['job_name'].dropna()}
    _toks = set(re.findall(r'\d{7}_\w+', _txt))
    _hit  = next((t for t in _toks if t in _jn), None)          # exact match
    if _hit:
        return _hit
    for t in sorted(_toks, key=len, reverse=True):              # tolerate a suffix
        for jn in _jn:
            if jn and (t.startswith(jn) or jn.startswith(t)):
                return jn
    return None


def mapped_resnum(row, col: str) -> "int | None":
    """Residue number from a ranked-sheet Mapped_* cell (e.g. 'ASP110' → 110, 'HIS280' → 280).

    These columns carry this homolog's alignment-against-control position for each
    catalytic role and differ substantially between jobs, so the QM region must be
    built from them per job rather than from reference numbering or a geometric
    guess. Returns None for empty/NaN/GAP/unmapped cells so the caller can fall
    back (alignment map → Smart-Lock geometry)."""
    try:
        _v = row.get(col)
    except AttributeError:
        _v = None
    if _v is None:
        return None
    _s = str(_v).strip()
    if not _s or _s.lower() in ("nan", "gap", "none", "unmapped", "-"):
        return None
    m = re.search(r"(\d+)", _s)
    return int(m.group(1)) if m else None


def load_watermap_csv(wm_csv_path: Path) -> list:
    """
    Load WaterMap thermodynamic sites from a Maestro 'Analyse WaterMap'
    CSV export (Rank_N.csv).

    Expected columns: Site, Occupancy, dH, -TdS, dG, #HB(WW), #HB(PW), #HB(LW).
    Returns list of dicts without spatial coordinates — use for statistics only.
    """
    if wm_csv_path is None or not wm_csv_path.exists():
        _name = wm_csv_path.name if wm_csv_path is not None else "(none found)"
        console_info(f"    {ConsoleColours.WARNING}[!] WaterMap CSV not found: {_name}{ConsoleColours.ENDC}")
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
    WaterMap hydration sites and their ΔG, from the WaterMap OUTPUT .maegz.

    A site is an atom carrying `r_watermap_deltaG`, and only such atoms are read. Accepting every
    atom in the file would take a WaterMap INPUT structure (`*_gpu-in.maegz`) — the protein itself —
    as thousands of ΔG = 0 "hydration sites". A file with no ΔG-bearing atom is not a water map and
    yields nothing.
    """
    sites = []
    if wm_path is None or not wm_path.exists():
        _name = wm_path.name if wm_path is not None else "(none found)"
        console_info(f"    {ConsoleColours.WARNING}[!] WaterMap .maegz not found: {_name}{ConsoleColours.ENDC}")
        return sites
    try:
        with structure.StructureReader(str(wm_path)) as reader:
            st = next(reader)
            for atom in st.atom:
                _dg = atom.property.get('r_watermap_deltaG')
                if _dg is None:          # not a hydration site — protein, ligand, bulk solvent
                    continue
                sites.append({
                    'pos': np.array(atom.xyz),
                    'dG':  float(_dg),
                    'num': atom.property.get('i_watermap_site_num', 0)
                })
        if not sites:
            console_info(f"    {ConsoleColours.WARNING}[!] {wm_path.name} carries no "
                         f"r_watermap_deltaG atoms — not a WaterMap result; no sites used."
                         f"{ConsoleColours.ENDC}")
        else:
            console_info(f"    {ConsoleColours.OKGREEN}✔ WaterMap: {len(sites)} hydration sites "
                         f"loaded from {wm_path.name}.{ConsoleColours.ENDC}")
    except Exception as e:
        console_info(f"    {ConsoleColours.WARNING}[!] WaterMap .maegz load failed: {e}{ConsoleColours.ENDC}")
    return sites


def check_md_equilibration(md_dir: Path, job_name: str) -> dict:
    """
    Verify that the NPT ensemble actually equilibrated before any frame is treated as a sample.

    Every statistic this script reports — NAC occupancy, dwell time, MM-GBSA — is an equilibrium
    average, and an equilibrium average over a system that is still relaxing is not an average of
    anything. Desmond's relaxation protocol usually settles the box within the first nanosecond, but
    'usually' is not a measurement: a badly packed box, a clashing prepared structure or a failed
    barostat all show up here as a volume that keeps drifting, and nowhere else.

    Desmond writes the thermodynamic stream to <job>.ene, whose header names the columns:
        0:time (ps)  1:E  2:E_p  3:E_k  4:E_c  5:E_x  6:E_f  7:P (bar)  8:V (A^3)  9:T (K)

    The box volume is the slow coordinate — it is what the barostat is still working on long after
    the temperature has settled — so equilibration is declared from V, and T is checked separately as
    a thermostat sanity test. Equilibration time = the first point after which the volume stays within
    MD_EQUIL_V_TOL_PCT of the production-window mean. The residual drift is then measured as a linear
    slope over that window: a box that is still shrinking or swelling has not equilibrated, however
    tight its instantaneous fluctuations look.

    Returns a dict of the measured quantities plus MD_Equilibrated; on any parse failure it returns
    MD_Equilibrated = None (unknown), never a silent True.
    """
    _out = {"MD_Equilibrated": None, "MD_Equil_Time_ps": np.nan,
            "MD_Mean_T_K": np.nan, "MD_SD_T_K": np.nan,
            "MD_Mean_V_A3": np.nan, "MD_V_Drift_Pct_per_ns": np.nan,
            "MD_Mean_P_bar": np.nan}
    try:
        _ene = md_dir / f"{job_name}.ene"
        if not _ene.exists() or _ene.stat().st_size == 0:
            console_info(f"    [!] No {_ene.name} — NPT equilibration cannot be verified; "
                         f"frame statistics are reported without it.")
            return _out

        _rows = []
        with _ene.open() as _fh:
            for _ln in _fh:
                if _ln.startswith("#"):
                    continue
                _p = _ln.split()
                if len(_p) < 10:
                    continue
                try:
                    _rows.append((float(_p[0]), float(_p[7]), float(_p[8]), float(_p[9])))
                except ValueError:
                    continue
        if len(_rows) < 100:
            return _out

        _a = np.asarray(_rows, float)                      # time, P, V, T
        _t, _P, _V, _T = _a[:, 0], _a[:, 1], _a[:, 2], _a[:, 3]

        # Production window: everything after the first MD_EQUIL_SKIP_FRAC of the run. Its mean is
        # the reference the equilibration time is measured against.
        _i0 = int(len(_t) * float(CFG.MD_EQUIL_SKIP_FRAC))
        _v_ref = float(_V[_i0:].mean())
        _tol = _v_ref * float(CFG.MD_EQUIL_V_TOL_PCT) / 100.0

        """
        Equilibration is judged on BLOCK MEANS of the volume, not on instantaneous values. The
        instantaneous box volume of an equilibrated NPT system still spikes past any sane tolerance a
        few times in a million steps — measured on this project's own trajectory, 0.02 % of points
        exceed 1 % of the mean while the box is demonstrably settled (residual drift −0.0001 %/ns).
        A test that demands every later point stay inside the tolerance therefore never passes until
        the final steps, and would discard an entire equilibrated trajectory. Averaging into blocks
        removes the fluctuation and leaves the drift, which is the thing being asked about.
        """
        _nb = min(int(CFG.MD_EQUIL_BLOCKS), max(2, len(_t) // 2))
        _bs = len(_t) // _nb
        _tb = _t[:_nb * _bs].reshape(_nb, _bs).mean(axis=1)
        _vb = _V[:_nb * _bs].reshape(_nb, _bs).mean(axis=1)

        """
        A block is settled when the REST of the run is settled — but `.all()` demands that every later
        block, without exception, sits inside the tolerance. One anomalous block near the end of an
        otherwise perfectly equilibrated run then invalidates everything before it, and equilibration is
        declared at the last block: a 1 µs trajectory collapses to a handful of production frames
        because of a single spike. The test is therefore on the FRACTION of later blocks that are
        settled, which is the question actually being asked — is the box still moving? — and is not
        hostage to one outlier.
        """
        _within = np.abs(_vb - _v_ref) <= _tol
        _need = float(CFG.MD_EQUIL_BLOCK_FRAC_MIN)
        _equil_b = next((_i for _i in range(_nb) if _within[_i:].mean() >= _need), _nb - 1)
        """
        The block test decided WHEN the box settled but never decided WHETHER it did. `_equil_b` fell back
        to the last block when no window met the fraction — and the verdict below then ignored that
        entirely and asked only whether the residual DRIFT was small. A box that oscillates around a
        stable mean has near-zero drift and no settled window at all: it passed. The fraction of settled
        blocks over the production window is therefore carried into the verdict as a first-class term.
        """
        _blocks_ok = bool(_within[_equil_b:].mean() >= _need)
        # The block's leading edge, not its midpoint: frames from the start of the settled block on
        # are samples.
        _equil_t = float(_t[_equil_b * _bs])

        # Residual volume drift over the production window, as % of the mean per nanosecond.
        _tw, _vw = _t[_i0:], _V[_i0:]
        _slope = float(np.polyfit(_tw, _vw, 1)[0]) if len(_tw) > 2 else 0.0   # A^3 per ps
        _drift = (_slope * 1000.0) / _v_ref * 100.0                            # % per ns

        _mean_T = float(_T[_i0:].mean())
        _ok_T = abs(_mean_T - float(CFG.MD_EQUIL_TARGET_T)) <= float(CFG.MD_EQUIL_T_TOL_K)
        _ok_V = abs(_drift) <= float(CFG.MD_EQUIL_V_DRIFT_MAX_PCT_NS)

        _out.update({
            "MD_Equilibrated": bool(_ok_T and _ok_V and _blocks_ok),
            "MD_Equil_Blocks_Settled_Frac": round(float(_within[_equil_b:].mean()), 3),
            "MD_Equil_Time_ps": round(_equil_t, 1),
            "MD_Mean_T_K": round(_mean_T, 2),
            "MD_SD_T_K": round(float(_T[_i0:].std()), 2),
            "MD_Mean_V_A3": round(_v_ref, 1),
            "MD_V_Drift_Pct_per_ns": round(_drift, 4),
            "MD_Mean_P_bar": round(float(_P[_i0:].mean()), 2),
        })

        if _out["MD_Equilibrated"]:
            console_info(f"    [i] NPT equilibrated at {_equil_t / 1000.0:.2f} ns — "
                         f"T {_mean_T:.1f}±{_out['MD_SD_T_K']:.1f} K, "
                         f"V drift {_drift:+.3f} %/ns.")
        else:
            _why = []
            if not _ok_T:
                _why.append(f"T {_mean_T:.1f} K deviates from {CFG.MD_EQUIL_TARGET_T:g} K "
                            f"by more than {CFG.MD_EQUIL_T_TOL_K:g} K")
            if not _ok_V:
                _why.append(f"box volume still drifting at {_drift:+.3f} %/ns "
                            f"(limit {CFG.MD_EQUIL_V_DRIFT_MAX_PCT_NS:g})")
            console_info(f"    {ConsoleColours.WARNING}[!] NPT NOT equilibrated: "
                         f"{'; '.join(_why)}. Frame averages from this trajectory are not "
                         f"equilibrium averages.{ConsoleColours.ENDC}")
        return _out
    except Exception as _e:                                 # noqa: BLE001
        console_info(f"    [!] Equilibration check failed ({_e}) — reported as unknown.")
        return _out


def load_watermap_reference_ca(wm_path: Path) -> dict:
    """The Cα coordinates of the structure the WaterMap sites were computed IN, keyed by residue.

    The sites are static coordinates in the WaterMap input's frame. The MD protein diffuses and
    tumbles through the box — measured on this project's own trajectory, the Cα centroid moves
    30-45 Å over 1 µs while the fold stays rigid (2-3 Å RMSD once superimposed). Comparing an MD
    coordinate with a static site coordinate therefore compares two unrelated frames, and every
    site match is meaningless without first superimposing the two.

    Returns {resnum: xyz}; the caller fits the frame onto these and moves the sites with it.
    """
    ref = {}
    _in = None
    for _pat in ("*_gpu-in.maegz", "*-in.maegz", "*_wm.maegz"):
        _hits = sorted(wm_path.parent.glob(_pat)) if wm_path else []
        if _hits:
            _in = _hits[0]
            break
    if _in is None:
        return ref
    try:
        with structure.StructureReader(str(_in)) as reader:
            st = next(reader)
            for a in st.atom:
                if a.pdbname.strip() == "CA":
                    ref[int(a.resnum)] = np.array(a.xyz)
    except Exception as e:
        console_info(f"    [!] WaterMap reference frame unreadable ({e}) — sites cannot be aligned.")
    return ref


def unwrap_ca_trace(ca_xyz, box):
    """Make a Cα trace whole across the periodic boundary before it is superimposed.

    Raw Desmond frames are WRAPPED: a protein straddling a box face has part of its backbone
    re-imaged to the far side, so its Cartesian trace is split in two even though the fold is
    perfectly intact. Kabsch cannot superimpose a split trace onto an intact reference — the
    RMSD explodes, and the fold check then rejects a good frame as a denatured one. The failure
    is silent and geometry-dependent: it fires only for the frames in which the protein happens
    to sit on a boundary.

    Consecutive Cα are ~3.8 Å apart, far below any half-box length, so the trace can be rebuilt
    without ambiguity by walking it and placing each atom at the nearest periodic image of its
    predecessor. This never moves an already-whole trace, so it is a no-op on frames that do not
    straddle the boundary.
    """
    _xyz = np.asarray(ca_xyz, dtype=float)
    if box is None or len(_xyz) < 2:
        return _xyz
    out = _xyz.copy()
    for i in range(1, len(out)):
        out[i] = out[i - 1] + get_mic_vector(out[i], out[i - 1], box)
    return out


def kabsch_transform(ref_xyz: np.ndarray, frame_xyz: np.ndarray):
    """The rigid transform that carries `ref_xyz` onto `frame_xyz` (rotation, then translation),
    with the post-superposition Cα RMSD.

    Standard Kabsch superposition, with the reflection guard: without it a degenerate SVD can
    return an improper rotation (a mirror image of the protein).

    The RMSD is returned because the transform is only meaningful while the fold is the same fold.
    A rigid transform will happily map any two point sets onto each other, so a frame in which the
    protein has unfolded, or in which the atom correspondence has broken, still yields a rotation —
    and every WaterMap site carried through it lands somewhere arbitrary. The caller thresholds on
    the RMSD (CFG.MD_FOLD_RMSD_MAX) and discards such frames rather than measuring against them.
    """
    _rc, _fc = ref_xyz.mean(axis=0), frame_xyz.mean(axis=0)
    _a, _b = ref_xyz - _rc, frame_xyz - _fc
    _U, _S, _Vt = np.linalg.svd(_a.T @ _b)
    _d = np.sign(np.linalg.det(_Vt.T @ _U.T))
    _R = _Vt.T @ np.diag([1.0, 1.0, _d]) @ _U.T
    _t = _fc - _R @ _rc
    _rmsd = float(np.sqrt((((ref_xyz @ _R.T + _t) - frame_xyz) ** 2).sum(axis=1).mean()))
    return _R, _t, _rmsd


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
    """Formats verbose job names to 'Rx_ProteinName_LigandName'. The 3R3U × FA positive control is
    labelled distinctly (Rx_3R3U-FA) so it never reads as a candidate fluoroacetate."""
    if str(job_name).startswith(str(getattr(CFG, "CONTROL_JOB_PREFIX", "0000000"))) or "3R3U" in str(job_name).upper():
        return f"R{int(rank)}_{_CTRL_LABEL}"
    clean = job_name.replace("desmond_md_job_", "").replace("_Prepared", "").replace("_CONTROL", "")
    clean = re.sub(r"^R(?:ank)?_\d+_", "", clean)
    parts = [p for p in clean.split('_') if p]
    if len(parts) >= 5:
        prot = "_".join(parts[2:-2])
        lig = parts[-1]
    elif len(parts) >= 3:
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
            # Short keywords (<=4 chars) match only as a whole underscore/space token,
            # so 'ang' cannot hit 'range', 'dist' cannot hit 'distribution', and 'prob'
            # cannot hit 'problem'. Longer keywords keep substring matching.
            import re as _re
            for k in keywords:
                kl = k.lower()
                for col in df.columns:
                    cl = col.lower()
                    if len(kl) <= 4:
                        if kl in _re.split(r'[^a-z0-9]+', cl):
                            return col
                    elif kl in cl:
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


# =============================================================================
# SECTION 5: VISUALISATION ENGINES
# =============================================================================

# -----------------------------------------------------------------------------
# SECTION 5a: PER-RANK & COMPARATIVE DASHBOARDS
# The standard MD-validation figures — per-rank NAC dashboard, the global
# comparative dashboard, and the viability/retention bar chart.
# -----------------------------------------------------------------------------
def generate_individual_dashboard(df: pd.DataFrame, job_name: str,
                                  output_path: Path, stats: dict) -> None:
    """2-panel per-job dashboard: SN2 scatter and dual-trace anchoring time series."""
    with PLOT_LOCK:
        sns.set_theme(style="whitegrid", context="paper")
        apply_figure_style(CFG)

        _C = CFG.DEFLUOR_FIG_COLOUR
        fig = plt.figure(figsize=(15, 6.5))
        gs  = gridspec.GridSpec(1, 2, width_ratios=[1, 1.2], wspace=0.15)

        # Panel 1: SN2 scatter with catalytic zones
        ax1 = fig.add_subplot(gs[0])
        ax1.add_patch(plt.Rectangle(
            (0, THRESHOLD_RELAXED_NAC_ANGLE), THRESHOLD_RELAXED_NAC_DIST,
            180 - THRESHOLD_RELAXED_NAC_ANGLE, color=_C['zone_relaxed'], alpha=0.3, zorder=0))
        ax1.add_patch(plt.Rectangle(
            (0, THRESHOLD_STRICT_NAC_ANGLE), THRESHOLD_STRICT_NAC_DIST,
            180 - THRESHOLD_STRICT_NAC_ANGLE, color=_C['zone_strict'], alpha=0.4, zorder=0))

        sc = ax1.scatter(df["NAC_Distance_A"], df["NAC_Angle_Deg"],
                         c=df["Frame"], cmap="viridis", s=20, alpha=0.7, edgecolor='none', zorder=2)
        if len(df.dropna(subset=["NAC_Distance_A", "NAC_Angle_Deg"])) > 10:
            sns.kdeplot(data=df, x="NAC_Distance_A", y="NAC_Angle_Deg", ax=ax1,
                        levels=5, color=_C['kde'], linewidths=1.0, alpha=0.5, zorder=3,
                        warn_singular=False)

        ax1.axvline(THRESHOLD_RELAXED_NAC_DIST, color=_C['warhead'], linestyle='--', linewidth=2.5)
        ax1.axhline(THRESHOLD_RELAXED_NAC_ANGLE, color=_C['tail'], linestyle='--', linewidth=2.5)
        ax1.set_xlabel("Nucleophile – ligand distance (Å)  ·  S$_N$2 reaction trajectory")
        ax1.set_ylabel("Attack Angle: O–C–F (°)")
        ax1.set_ylim(0, 180); ax1.set_xlim(left=0)
        
        # Legend moved to bottom right and contains zones
        ax1.legend(handles=[
            Line2D([0], [0], color=_C['warhead'], linestyle='--', lw=2.5,
                   label=f'Relaxed Dist < {THRESHOLD_RELAXED_NAC_DIST}Å'),
            Line2D([0], [0], color=_C['tail'], linestyle='--', lw=2.5,
                   label=f'Relaxed Angle > {THRESHOLD_RELAXED_NAC_ANGLE}°'),
            Patch(facecolor=_C['zone_relaxed'], alpha=0.3, label='Relaxed S_N2 Zone'),
            Patch(facecolor=_C['zone_strict'], alpha=0.4, label='Strict S_N2 Zone'),
        ], loc='lower right', frameon=True,
           edgecolor=_C['legend_edge'], fancybox=True, fontsize=CFG.VIS_FONT_LEGEND)
           
        cbar = plt.colorbar(sc, ax=ax1, pad=0.02)
        cbar.set_label("Simulation Frame (Time)", rotation=270, labelpad=15)
        clean_spines(ax1)

        # Panel 2: Lock & Key dual-trace
        ax2 = fig.add_subplot(gs[1])
        window = min(25, max(1, len(df) // 10))
        df = df.copy()
        df['Warhead_Smooth'] = df["NAC_Distance_A"].rolling(window=window, min_periods=1).mean()

        ax2.plot(df["Frame"], df["NAC_Distance_A"],   color=_C['warhead'], linewidth=1.0, alpha=0.15)
        ax2.plot(df["Frame"], df['Warhead_Smooth'],   color=_C['warhead'], linewidth=2.5, alpha=0.95,
                 label='Warhead Anchor (Nuc – LigC)')
        ax2.axhline(CFG.NAC_DIST_RELAXED, color=_C['warhead'], linestyle=':', linewidth=1.5, alpha=0.7)

        if "Tail_Cradle_Dist_A" in df.columns and not df["Tail_Cradle_Dist_A"].isna().all():
            df['Tail_Smooth'] = df["Tail_Cradle_Dist_A"].rolling(window=window, min_periods=1).mean()
            ax2.plot(df["Frame"], df["Tail_Cradle_Dist_A"], color=_C['tail'], linewidth=1.0, alpha=0.15)
            ax2.plot(df["Frame"], df['Tail_Smooth'],        color=_C['tail'], linewidth=2.5, alpha=0.95,
                     label='Tail Anchor (Cradle – LigF)')
            ax2.axhline(CFG.MECH_CRADLE_RADIUS, color=_C['tail'], linestyle=':', linewidth=1.5, alpha=0.7)

        data_max = df["NAC_Distance_A"].max() if not df["NAC_Distance_A"].isna().all() else 12.0
        ax2.set_ylim(1.5, max(12.0, data_max * 1.4))
        ax2.set_xlabel("Simulation frame")
        ax2.set_ylabel("Active-site anchoring — interaction distance (Å)")
        ax2.legend(loc='upper right', frameon=True,
                   edgecolor=_C['legend_edge'], fancybox=True, fontsize=CFG.VIS_FONT_LEGEND)
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
                 fontsize=CFG.VIS_FONT_ANNOT, va='top', ha='left', zorder=10,
                 bbox=dict(facecolor='white', edgecolor=CFG.VIS_MD_PALETTE["border"],
                           boxstyle='round,pad=0.6', alpha=1.0))

        # No on-figure title; the descriptive metadata is written to the log instead
        plt.savefig(output_path, dpi=int(CFG.VIS_FIGURE_DPI), bbox_inches='tight')
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
    apply_figure_style(CFG)
    _C = CFG.DEFLUOR_FIG_COLOUR

    all_data = []
    for _, row in df_master.iloc[::-1].iterrows():
        csv_path = next(iter(sorted(out_dir.glob(f"**/{row['Job_Name']}_NAC_Data.csv"))), None)
        if csv_path and csv_path.exists():
            df_job = pd.read_csv(csv_path)
            df_job['Job']  = format_job_label(row['Job_Name'], row['Scientific_Rank'])
            df_job['Rank'] = row['Scientific_Rank']
            all_data.append(df_job)

    if not all_data:
        return
    combined_df   = pd.concat(all_data, ignore_index=True)
    
    # Consistent colour map based on ascending Scientific_Rank
    sorted_df = df_master.sort_values('Scientific_Rank', ascending=True)
    sorted_labels = [format_job_label(row['Job_Name'], row['Scientific_Rank']) for _, row in sorted_df.iterrows()]
    job_colour_map = dict(zip(sorted_labels, sns.color_palette("tab20", n_colors=len(sorted_labels))))
    for _lab in list(job_colour_map):          # 3R3U-FA control gets the one distinct control colour
        if _CTRL_LABEL in str(_lab):
            job_colour_map[_lab] = _CTRL_COLOUR

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
                   order=sorted_labels, palette=job_colour_map, inner="quartile", linewidth=1.2)
    ax1.axvspan(0, transform_distance(THRESHOLD_RELAXED_NAC_DIST), color=CFG.VIS_ACCENT["green"], alpha=0.15, zorder=0)
    ax1.axvline(transform_distance(THRESHOLD_RELAXED_NAC_DIST), color=CFG.VIS_ACCENT["vermillion"], linestyle='--', linewidth=2)
    ax1.set_xlabel("Nucleophile–ligand distance (Å) over MD frames [non-linear scale]"); ax1.set_ylabel("")
    
    dist_ticks = [0, 1, 2, 3, 4, 5, 10, 15, 20, 30, 40, 50]
    ax1.set_xticks([transform_distance(t) for t in dist_ticks])
    ax1.set_xticklabels([str(t) for t in dist_ticks], rotation=90)
    ax1.set_xlim(left=0)
    
    fig.canvas.draw()
    for lbl in ax1.get_yticklabels():
        lbl.set_color(job_colour_map.get(lbl.get_text(), CFG.VIS_MD_PALETTE["text_default"]))
        lbl.set_fontweight('bold'); lbl.set_fontsize(11)
    clean_spines(ax1)

    # Panel 2: Attack Angle Distribution violin plot
    ax2 = fig.add_subplot(gs[1], sharey=ax1)
    sns.violinplot(data=combined_df, y="Job", x="NAC_Angle_Deg", ax=ax2,
                   order=sorted_labels, palette=job_colour_map, inner="quartile", linewidth=1.2)
    ax2.axvspan(THRESHOLD_RELAXED_NAC_ANGLE, 180, color=CFG.VIS_ACCENT["green"], alpha=0.15, zorder=0)
    ax2.axvline(THRESHOLD_RELAXED_NAC_ANGLE, color=CFG.VIS_ACCENT["blue"], linestyle='--', linewidth=2)
    ax2.set_xlabel("S$_N$2 attack angle (°) over MD frames"); ax2.set_ylabel("")
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
        180 - THRESHOLD_RELAXED_NAC_ANGLE, color=_C['zone_relaxed'], alpha=0.3, zorder=0))
    ax3.add_patch(plt.Rectangle(
        (0, THRESHOLD_STRICT_NAC_ANGLE), strict_w,
        180 - THRESHOLD_STRICT_NAC_ANGLE, color=_C['zone_strict'], alpha=0.4, zorder=0))
        
    ax3.axvline(transform_distance(THRESHOLD_RELAXED_NAC_DIST), color=CFG.VIS_ACCENT["vermillion"], linestyle='--', linewidth=2)
    ax3.axhline(THRESHOLD_RELAXED_NAC_ANGLE, color=CFG.VIS_ACCENT["blue"], linestyle='--', linewidth=2)
    ax3.set_xlabel("Nucleophile–ligand distance (Å) [non-linear scale]  ·  global catalytic landscape (point colour = job)")
    ax3.set_ylabel("S$_N$2 attack angle O–C–F (°)")
    
    ax3.set_xticks([transform_distance(t) for t in dist_ticks])
    ax3.set_xticklabels([str(t) for t in dist_ticks], rotation=90)
    ax3.set_xlim(0, transform_distance(50))
    ax3.set_ylim(0, 180)
    
    # Legend at bottom left containing both lines and zones
    ax3.legend(handles=[
        Line2D([0], [0], color=CFG.VIS_ACCENT["vermillion"], linestyle='--', lw=2,
               label=f'Distance < {THRESHOLD_RELAXED_NAC_DIST}Å'),
        Line2D([0], [0], color=CFG.VIS_ACCENT["blue"], linestyle='--', lw=2,
               label=f'Angle > {THRESHOLD_RELAXED_NAC_ANGLE}°'),
        Patch(facecolor=_C['zone_relaxed'], alpha=0.3, label='Relaxed S_N2 Zone'),
        Patch(facecolor=_C['zone_strict'], alpha=0.4, label='Strict S_N2 Zone'),
    ], loc='lower left', frameon=True,
       edgecolor=CFG.VIS_MD_PALETTE["border_light"], fancybox=True)
    clean_spines(ax3)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        plt.tight_layout()
    out_path = out_dir / "02_MD_Comparative_Analysis.png"
    plt.savefig(out_path, dpi=int(getattr(CFG, "VIS_FIGURE_DPI", 300)), bbox_inches='tight')
    plt.close(fig)
    console_info(f"    Comparative Dashboard Saved : {out_path.resolve()}")


def generate_comparative_residue_engagement(out_dir: Path, df_master: pd.DataFrame) -> None:
    """
    One comparative view of catalytic-machinery engagement across every SN2 case.

    Each row is a candidate (ligand + Scientific_Rank), each column a catalytic
    residue role (nucleophile, acid, the His stabiliser and the Trp/Tyr fluoride
    cradle). The cell is that residue's mean distance to the warhead carbon over
    the strict-NAC frames (DT_<role>_NAC_Mean_A) — the same engagement metric the
    per-case *_MMGBSA_NAC_Decomposition figure shows, here pooled so one can read
    off, at a glance, which residue closes in (green, short distance) and which
    stays disengaged (red, long distance) in each case. Lower = more engaged.
    """
    _roles = [
        ("Nucleophile\n(Asp)",   "DT_Nuc_NAC_Mean_A"),
        ("Acid\n(Asp)",          "DT_Acid_NAC_Mean_A"),
        ("Stabiliser\n(His)",    "DT_StabH_NAC_Mean_A"),
        ("Cradle\n(Trp)",        "DT_StabW_NAC_Mean_A"),
        ("Cradle\n(Tyr)",        "DT_StabY_NAC_Mean_A"),
    ]
    _cols  = [(lbl, col) for lbl, col in _roles if col in df_master.columns]
    if not _cols or df_master.empty:
        console_info("    [!] Comparative residue engagement skipped — no DT_*_NAC_Mean_A columns.")
        return

    _sorted = df_master.sort_values("Scientific_Rank", ascending=True)
    _labels = [format_job_label(r["Job_Name"], r["Scientific_Rank"]) for _, r in _sorted.iterrows()]
    _matrix = _sorted[[c for _, c in _cols]].apply(pd.to_numeric, errors="coerce")
    _matrix.index   = _labels
    _matrix.columns = [lbl for lbl, _ in _cols]

    """
    The colour scale runs over the CRITERION BANDS, not over an arbitrary distance window: green at
    the reactive contact (NAC_DIST_STRICT — the geometry the mechanism requires) through to red at
    the electrostatic limit (THRESHOLD_SALT_BRIDGE — beyond which the residue is not in contact at
    all). The engagement figure judges the same distances against the same four cut-offs, so a cell
    and a bar now mean the same thing. Anything past the outer band saturates red: how far beyond
    'not in contact' a residue sits carries no further meaning.
    """
    _contact = float(CFG.NAC_DIST_STRICT)
    _outer   = float(CFG.THRESHOLD_SALT_BRIDGE)

    _h = max(3.2, 0.42 * len(_labels) + 1.6)
    _w = max(6.0, 1.5 * len(_cols) + 2.5)
    fig, ax = plt.subplots(figsize=(_w, _h))
    _cmap = plt.get_cmap(CFG.ENGAGE_HEATMAP_CMAP).copy()
    _cmap.set_bad(color=CFG.ENGAGE_HEATMAP_NAN)   # missing residue → grey
    sns.heatmap(
        _matrix, ax=ax, cmap=_cmap, vmin=_contact, vmax=_outer,
        annot=True, fmt=".1f", annot_kws={"fontsize": CFG.VIS_FONT_ANNOT},
        linewidths=0.6, linecolor="white",
        cbar_kws={"label": f"Mean distance to warhead C in strict-NAC frames (Å)\n"
                           f"{_contact:g} = reactive contact · {_outer:g} = electrostatic limit"},
    )
    ax.set_xlabel("Catalytic residue role  ·  green = engaged, red = out of contact")
    ax.set_ylabel("SN2 case (ligand · Scientific_Rank)")
    ax.tick_params(axis="x", labelrotation=0)
    ax.tick_params(axis="y", labelrotation=0)
    plt.setp(ax.get_yticklabels(), fontsize=8)

    out_path = out_dir / "04_Comparative_Residue_Engagement.png"
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        plt.savefig(out_path, dpi=int(getattr(CFG, "VIS_FIGURE_DPI", 300)), bbox_inches="tight")
    plt.close(fig)
    console_info(f"    Comparative Residue Engagement Saved : {out_path.resolve()}")


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
    apply_figure_style(CFG)
    from matplotlib.colors import to_rgba

    df_plot = df_master.sort_values('Scientific_Rank', ascending=True).copy()

    # Consistent colour map based on ascending Scientific_Rank (tab20)
    sorted_labels = [format_job_label(r['Job_Name'], r['Scientific_Rank']) for _, r in df_plot.iterrows()]
    job_colour_map = dict(zip(sorted_labels, sns.color_palette("tab20", n_colors=len(sorted_labels))))
    for _lab in list(job_colour_map):          # 3R3U-FA control gets the one distinct control colour
        if _CTRL_LABEL in str(_lab):
            job_colour_map[_lab] = _CTRL_COLOUR

    n_rows  = len(df_plot)
    row_h   = 0.85
    fig_h   = max(4.0, n_rows * row_h + 1.8)
    fig, (ax, ax_ann) = plt.subplots(
        1, 2, figsize=(15, fig_h),
        gridspec_kw={'width_ratios': [2.8, 1], 'wspace': 0.01})

    bar_h       = 0.72
    y_positions = list(range(n_rows))

    for i, (_, row) in enumerate(df_plot.iterrows()):
        y = y_positions[i]
        job_label = format_job_label(row['Job_Name'], row['Scientific_Rank'])
        rank_color = job_colour_map.get(job_label, CFG.VIS_MD_PALETTE["rank_default"])

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

        # Plot the 4 parallel sub-bars with distinct functional colours
        sub_bars = [
            (pocket,        CFG.VIS_MD_PALETTE["slate"], -0.21),
            (triad_total,   CFG.VIS_MD_PALETTE["triad"], -0.07),
            (relaxed_total, CFG.VIS_MD_PALETTE["relaxed"],  0.07),
            (strict_total,  CFG.VIS_MD_PALETTE["strict"],   0.21)
        ]
        
        sub_bar_h = 0.13
        for val, col, offset in sub_bars:
            if val > 1e-4:
                ax.barh(y + offset, val, height=sub_bar_h, color=col, edgecolor='none', zorder=2)
                # Value label (in solid black)
                if val >= 8.0:
                    ax.text(val - 1.0, y + offset, f'{val:.1f}%',
                            ha='right', va='center', fontsize=7.0, fontweight='bold',
                            color='black', zorder=4)
                else:
                    ax.text(val + 0.5, y + offset, f'{val:.1f}%',
                            ha='left', va='center', fontsize=7.0, fontweight='bold',
                            color='black', zorder=4)
            else:
                # Value label for 0% (in solid black)
                ax.text(0.5, y + offset, "0.0%",
                        ha='left', va='center', fontsize=7.0, fontweight='bold',
                        color='black', zorder=4)

    ax.set_yticks(y_positions)
    ax.set_yticklabels(
        [format_job_label(r['Job_Name'], r['Scientific_Rank']) for _, r in df_plot.iterrows()],
        fontsize=10.5, fontweight='bold')
    ax.set_xlim(0, 105)
    ax.set_ylim(-0.65, n_rows - 0.35)
    ax.invert_yaxis()
    ax.set_xlabel("Percentage of Simulation Time (%)")
    ax.axvline(100, color=CFG.VIS_MD_PALETTE["slate"], linestyle=':', linewidth=1.0, alpha=0.6, zorder=1)
    clean_spines(ax)

    # Set y-axis tick label colours to match job colours
    fig.canvas.draw()
    for lbl in ax.get_yticklabels():
        lbl.set_color(job_colour_map.get(lbl.get_text(), CFG.VIS_MD_PALETTE["text_default"]))

    # ── Annotation panel ──────────────────────────────────────────────────────
    ax_ann.set_xlim(0, 1)
    ax_ann.set_ylim(-0.65, n_rows - 0.35)
    ax_ann.invert_yaxis()
    ax_ann.set_yticks([])
    ax_ann.set_xticks([])
    for sp in ax_ann.spines.values():
        sp.set_visible(False)
    ax_ann.axvline(0.0, color=CFG.VIS_MD_PALETTE["border"], linewidth=0.8)

    hdr_y = -0.45
    ax_ann.text(0.10, hdr_y, 'WM ΔG\n(kcal/mol)', ha='center', va='center',
                 fontsize=7.0, fontweight='bold', color=CFG.VIS_MD_PALETTE["text_dark"])
    ax_ann.text(0.30, hdr_y, 'WM_N\n(stable)',     ha='center', va='center',
                 fontsize=7.0, fontweight='bold', color=CFG.VIS_MD_PALETTE["text_dark"])
    ax_ann.text(0.50, hdr_y, 'min d_NAC\n(Å)',     ha='center', va='center',
                 fontsize=7.0, fontweight='bold', color=CFG.VIS_MD_PALETTE["text_dark"])
    ax_ann.text(0.72, hdr_y, 'avg a_NAC\n(°)',     ha='center', va='center',
                 fontsize=7.0, fontweight='bold', color=CFG.VIS_MD_PALETTE["text_dark"])
    ax_ann.text(0.92, hdr_y, 'DT\n(n)',            ha='center', va='center',
                 fontsize=7.0, fontweight='bold', color=CFG.VIS_MD_PALETTE["text_dark"])

    for i, (_, row) in enumerate(df_plot.iterrows()):
        y    = y_positions[i]
        wm   = row.get('WM_Mean_dG',        float('nan'))
        wmn  = row.get('WM_N_Stable',       float('nan'))
        dist = row.get('MD_Min_NAC_Dist_A',  float('nan'))
        ang  = row.get('MD_Avg_NAC_Angle_Deg', float('nan'))
        dt   = row.get('Dream_Team_Mapped',  0)

        try:    wm_str   = f'{float(wm):.2f}'
        except (ValueError, TypeError): wm_str   = 'N/A'
        try:    wmn_str  = str(int(float(wmn)))
        except (ValueError, TypeError): wmn_str  = 'N/A'
        try:    dist_str = f'{float(dist):.2f}'
        except (ValueError, TypeError): dist_str = 'N/A'
        try:    ang_str  = f'{float(ang):.1f}°'
        except (ValueError, TypeError): ang_str  = 'N/A'
        try:    dt_str   = str(int(float(dt or 0)))
        except (ValueError, TypeError): dt_str   = '–'

        dist_col = (CFG.VIS_MD_PALETTE["good"] if (dist_str != 'N/A' and float(dist) < CFG.NAC_DIST_RELAXED)
                    else (CFG.VIS_MD_PALETTE["marginal"] if (dist_str != 'N/A' and float(dist) < CFG.NAC_DIST_MARGINAL_MAX)
                    else CFG.VIS_MD_PALETTE["fail"]))

        ang_col = (CFG.VIS_MD_PALETTE["good"] if (ang_str != 'N/A' and float(ang) > CFG.NAC_ANGLE_RELAXED)
                   else (CFG.VIS_MD_PALETTE["marginal"] if (ang_str != 'N/A' and float(ang) > CFG.SN2_ANGLE_MARGINAL_MIN)
                   else CFG.VIS_MD_PALETTE["fail"]))

        ax_ann.text(0.10, y, wm_str,   ha='center', va='center', fontsize=8,
                    color=CFG.VIS_MD_PALETTE["text_dark"], fontweight='bold')
        ax_ann.text(0.30, y, wmn_str,  ha='center', va='center', fontsize=8,
                    color=CFG.VIS_MD_PALETTE["accent_blue"], fontweight='bold')
        ax_ann.text(0.50, y, dist_str, ha='center', va='center', fontsize=8,
                    color=dist_col,  fontweight='bold')
        ax_ann.text(0.72, y, ang_str,  ha='center', va='center', fontsize=8,
                    color=ang_col,   fontweight='bold')
        ax_ann.text(0.92, y, dt_str,   ha='center', va='center', fontsize=8,
                    color=CFG.VIS_MD_PALETTE["accent_purple"], fontweight='bold')


    # ── Legend ────────────────────────────────────────────────────────────────
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor=CFG.VIS_MD_PALETTE["slate"], edgecolor='none', label='Pocket Retention'),
        Patch(facecolor=CFG.VIS_MD_PALETTE["triad"], edgecolor='none', label='Triad Integrity (Total)'),
        Patch(facecolor=CFG.VIS_MD_PALETTE["relaxed"], edgecolor='none', label='Relaxed Catalysis (Total)'),
        Patch(facecolor=CFG.VIS_MD_PALETTE["strict"], edgecolor='none', label='Strict Catalysis (Total)'),
    ]
    ax.legend(handles=legend_handles, loc='lower left', bbox_to_anchor=(0.0, 1.02),
              frameon=True,  edgecolor=CFG.VIS_MD_PALETTE["slate"],  ncol=4)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        plt.tight_layout()
    out_path = out_dir / "03_MD_Viability_Summary.png"
    plt.savefig(out_path, dpi=int(getattr(CFG, "VIS_FIGURE_DPI", 300)), bbox_inches='tight')
    plt.close(fig)
    console_info(f"    Viability Bar Chart Saved   : {out_path.resolve()}")


# -----------------------------------------------------------------------------
# SECTION 5b: DEFLUORINATION FIGURES (verdict landscape + reactive-state decomp)
# Colours/thresholds from CFG.DEFLUOR_FIG_COLOUR / DEFLUOR_ENGAGE_* (SSOT).
# -----------------------------------------------------------------------------
def generate_defluorination_landscape(out_dir: Path, df_master: pd.DataFrame) -> None:
    """The whole-story figure. Every candidate is placed by catalytic PERSISTENCE
    (x — longest continuous strict-NAC dwell, ns) against its QM/MM SN2 BARRIER
    (y — ΔE‡ kcal/mol, lower = more accessible transition state). Marker size
    encodes NAC-conditioned MM-GBSA binding strength (a stable reactive pose) and
    colour encodes the fused defluorination propensity. The shaded green quadrant
    is the competence gate (dwell ≥ Y and ΔE‡ ≤ Z): points inside it, brightest,
    are predicted to defluorinate; tight binders outside it are the "bind-but-do-
    not-react" decoys. One figure separates catalysis from mere affinity."""
    try:
        d = df_master.copy()
        _lab = "Job_Name" if "Job_Name" in d.columns else d.columns[0]
        d["_label"] = [format_job_label(r.get(_lab, "?"), r.get("Scientific_Rank", i + 1))
                       for i, (_, r) in enumerate(d.iterrows())]
        _has_bar = "QSite_Barrier_kcal" in d.columns and pd.to_numeric(
            d["QSite_Barrier_kcal"], errors="coerce").notna().any()
        x = pd.to_numeric(d.get("NAC_Dwell_Max_ns"), errors="coerce")
        if _has_bar:
            y = pd.to_numeric(d["QSite_Barrier_kcal"], errors="coerce")
            ylab = "QM/MM S$_N$2 barrier  ΔE‡  (kcal/mol) — lower = more reactive"
        else:
            y = pd.to_numeric(d.get("Strict_Viability_Pct"), errors="coerce")
            ylab = "Strict-NAC viability (%) — QM/MM barrier pending"
        m = x.notna() & y.notna()
        if int(m.sum()) == 0:
            console_info("    [!] Defluorination landscape skipped — no dwell/barrier data yet.")
            return
        d, x, y = d[m].reset_index(drop=True), x[m].reset_index(drop=True), y[m].reset_index(drop=True)

        _dg = None
        for _c in ("MMGBSA_dG_NAC_Mean_kcal", "MMGBSA_dG_Global_Mean_kcal", "MMGBSA_dG_ArithMean_kcal"):
            if _c in d.columns and pd.to_numeric(d[_c], errors="coerce").notna().any():
                _dg = pd.to_numeric(d[_c], errors="coerce").abs(); break
        _sz = ((60 + 260 * (_dg / _dg.max())).fillna(110)
               if _dg is not None and _dg.max() and _dg.max() > 0 else pd.Series(130, index=d.index))
        _col = pd.to_numeric(d.get("Defluor_Propensity_Norm"), errors="coerce")

        Y = float(getattr(CFG, "DEFLUOR_DWELL_MIN_NS", 1.0))
        Z = float(getattr(CFG, "DEFLUOR_BARRIER_MAX_KCAL", 22.0))
        fig, ax = plt.subplots(figsize=(11, 7.5))
        # Floors keep the axes non-degenerate when the data is sparse/all-zero (e.g. a single
        # rank or a strided test) — a singular xlim/ylim otherwise warns and collapses the frame.
        _xhi = max(float(x.max()) * 1.12, Y * 1.5, 1.0)
        _ylo = min(float(y.min()) * 0.9, 0.0)
        _yhi = max(float(y.max()) * 1.12, (Z * 1.25 if _has_bar else float(y.max()) * 1.12), _ylo + 1.0)
        _C = CFG.DEFLUOR_FIG_COLOUR
        ax.set_xlim(0, _xhi); ax.set_ylim(_ylo, _yhi)
        if _has_bar:
            ax.add_patch(plt.Rectangle((Y, _ylo), _xhi - Y, Z - _ylo,
                                       color=_C["gate"], alpha=0.09, zorder=0))
            ax.axhline(Z, color=_C["gate_line"], ls="--", lw=1.2, zorder=1)
            ax.axvline(Y, color=_C["gate_line"], ls="--", lw=1.2, zorder=1)
            ax.text(_xhi * 0.98, _ylo + (Z - _ylo) * 0.5,
                    f"defluorination-competent\n(dwell ≥ {Y:g} ns, ΔE‡ ≤ {Z:g})",
                    ha="right", va="center", fontsize=8, color=_C["gate_text"], style="italic")
        sc = ax.scatter(x, y, s=_sz, c=(_col if _col.notna().any() else _C["scatter"]),
                        cmap="viridis", vmin=0, vmax=1, edgecolor=_C["edge"],
                        linewidth=0.8, alpha=0.92, zorder=5)
        for xi, yi, lab in zip(x, y, d["_label"]):
            ax.annotate(str(lab), (xi, yi), fontsize=7, xytext=(4, 4),
                        textcoords="offset points", zorder=6)
        if _col.notna().any():
            cb = fig.colorbar(sc, ax=ax, pad=0.02)
            cb.set_label("Defluorination propensity  (log-scaled, best = 1)", fontsize=9)
        ax.set_xlabel("Catalytic persistence — longest continuous strict-NAC dwell (ns)")
        ax.set_ylabel(ylab)
        clean_spines(ax)
        # Disclosure: the MD this dwell is measured on runs under the Step-06 ligand positional
        # restraint (CFG.MD_RESTRAIN_LIGAND), so persistence reflects NAC geometry SUSTAINED under
        # restraint — the ligand is held near its pose by design, not free to escape. The QM/MM ΔE‡
        # (y-axis / colour) is the unrestrained arbiter of turnover.
        if getattr(CFG, "MD_RESTRAIN_LIGAND", False):
            fig.text(0.5, 0.005,
                     "NAC dwell measured under the Step-06 ligand positional restraint — persistence is "
                     "restraint-sustained, not spontaneous; the QM/MM ΔE‡ is the unrestrained turnover arbiter.",
                     ha="center", va="bottom", fontsize=7.5, color="0.45", wrap=True)
        out_path = out_dir / "05_Defluorination_Landscape.png"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            plt.savefig(out_path, dpi=int(getattr(CFG, "VIS_FIGURE_DPI", 300)), bbox_inches="tight")
        plt.close(fig)
        console_info(f"    Defluorination Landscape Saved : {out_path.resolve()}")
    except Exception as _e:
        console_info(f"    [!] Defluorination landscape failed ({_e}).")


# -----------------------------------------------------------------------------
# SECTION 5c: REACTIVE-POSE FIGURES (MM-GBSA decomposition · machinery engagement)
#
# Both figures ask the same question — what CHANGES when the ligand reaches the
# reactive geometry — and both are built from the per-frame tables the run already
# writes: <job>_NAC_Data.csv (geometry, one row per frame) and the frame-stamped
# MM-GBSA CSV from Step 06. Each is produced per candidate and once merged across
# candidates. Every colour, cut-off and threshold comes from CFG.
# -----------------------------------------------------------------------------
"""
The catalytic machinery, in the order it acts: the nucleophile attacks the warhead carbon, the
acid/base pair runs the proton chemistry, the clamp holds the carboxylate, the cradle stabilises the
departing fluoride. Each residue's per-frame distance column is paired with its CFG role-group
colour and with the ranked CSV's alignment column that names the residue in THIS homolog.
"""
_ENGAGE_ROLES = [
    ("DT_Nuc_LigC_A",    "Nucleophile",         "Nuc",        "Mapped_Nucleophile"),
    ("DT_Base_LigC_A",   "Acid/base catalysis", "Base",       "Mapped_Base"),
    ("DT_Acid_LigC_A",   "Acid/base catalysis", "Acid",       "Mapped_Acid"),
    ("DT_Clamp1_LigC_A", "Carboxylate clamp",   "Clamp 1",    "Mapped_Clamp1"),
    ("DT_Clamp2_LigC_A", "Carboxylate clamp",   "Clamp 2",    "Mapped_Clamp2"),
    ("DT_StabH_LigC_A",  "Fluoride pocket",     "Cradle His", "Mapped_Stabiliser_H"),
    ("DT_StabW_LigC_A",  "Fluoride pocket",     "Cradle Trp", "Mapped_Stabiliser_W"),
    ("DT_StabY_LigC_A",  "Fluoride pocket",     "Cradle Tyr", "Mapped_Stabiliser_Y"),
]
_MMGBSA_TERMS = ("Coulomb", "vdW", "Solv_GB", "Lipo", "Hbond", "Packing", "Covalent", "SelfCont")


def _engage_zones() -> list:
    """The criteria that define engagement — every one of them a CFG constant."""
    return [
        (float(CFG.NAC_DIST_STRICT),       "reactive contact"),
        (float(CFG.THRESHOLD_HB_DIST_MAX), "H-bond range"),
        (float(CFG.NAC_DIST_RELAXED),      "relaxed NAC"),
        (float(CFG.THRESHOLD_SALT_BRIDGE), "electrostatic range"),
    ]


def _darken(colour, f: float = 0.55) -> tuple:
    """A darker shade of a series colour, so a mark drawn over its own bars stays readable."""
    return tuple(c * f for c in mcolors.to_rgb(colour))


def _master_container(ax, x: float, colour, width: float = 0.9,
                      lo: float = 0.0, hi: float = 1.0) -> None:
    """A master container: a faint, colour-outlined bar drawn BEHIND a group of child bars.

    It carries no value — it exists to make the group read as one object and to give the group a
    colour identity the eye can follow across the panel.
    """
    ax.bar(x, hi - lo, width, bottom=lo, color=mcolors.to_rgba(colour, 0.10),
           edgecolor=colour, linewidth=1.4, zorder=2)


def _reactive_frames(nac: pd.DataFrame) -> "tuple[np.ndarray, str]":
    """The reactive frames, and which criterion produced them.

    Strict NAC is the definition the mechanism rests on, but a candidate can have almost none of
    them (a pose that reaches the reactive geometry only a handful of times in a microsecond).
    Falling back to the geometric NAC keeps the figure honest — the label says which criterion was
    used, so a thin ensemble can never be mistaken for a rich one.
    """
    strict = nac.loc[nac["NAC_Strict_Pass"] == 1, "Frame"].to_numpy(dtype=int) \
        if "NAC_Strict_Pass" in nac.columns else np.array([], dtype=int)
    if len(strict) >= 30:
        return strict, "strict NAC"
    geom = nac.loc[nac["NAC_Geom_Pass"] == 1, "Frame"].to_numpy(dtype=int) \
        if "NAC_Geom_Pass" in nac.columns else np.array([], dtype=int)
    return geom, "geometric NAC"


def _load_reactive_pose_data(out_dir: Path) -> list:
    """Gather, per candidate: the per-frame NAC table, the frame-stamped MM-GBSA table, and the
    alignment map that names each catalytic residue in that homolog."""
    run_dir = out_dir.parent
    md = run_dir / "6_Physics_Validation" / "05_MD_Simulations"
    out = []
    for d in sorted(out_dir.glob("Rank_*")):
        m = re.match(r"Rank_(\d+)_", d.name)
        if not m:
            continue
        rank = int(m.group(1))
        nac_csv = next(iter(sorted(d.glob(f"*{CFG.SUFFIX_NAC_DATA}"))), None)
        _mgdir = md / f"desmond_md_job_R_{rank}"
        mg_csv = _mgdir / f"desmond_md_job_R_{rank}{CFG.SUFFIX_MMGBSA_CSV}"
        if not mg_csv.is_file():        # tolerant fallback, same as the engine's discovery ladder
            mg_csv = next(iter(sorted(_mgdir.glob("*mmgbsa*.csv"))), mg_csv)
        if nac_csv is None:
            continue
        job = nac_csv.name.replace(CFG.SUFFIX_NAC_DATA, "")
        _full = re.sub(r"^\d+_", "", job.split("_")[-1]) if "_" in job else job
        # The 3R3U × FA positive control shares the ligand 'FA' with the candidate fluoroacetate;
        # label it distinctly (3R3U-FA) so every downstream figure names + colours it as the control.
        _is_ctrl = job.startswith(str(getattr(CFG, "CONTROL_JOB_PREFIX", "0000000"))) or "3R3U" in job.upper()
        lig = _CTRL_LABEL if _is_ctrl else CFG.VIS_LIGAND_SHORT.get(_full.lower(), _full)
        out.append({
            "rank": rank, "job": job, "ligand": lig, "is_control": _is_ctrl, "mapped": {},
            "nac": pd.read_csv(nac_csv),
            "mmgbsa": pd.read_csv(mg_csv) if mg_csv.is_file() else pd.DataFrame(),
            "dir": d,
        })

    _prod = run_dir / "1_Boltz2_Production"
    ranked = (_utils_mod.latest_by_mtime(_prod.glob(CFG.GLOB_RANKED_CSV))
              or _utils_mod.latest_by_mtime(_prod.glob("*Ranked*.csv")))
    if ranked is not None:
        rk = pd.read_csv(ranked, low_memory=False)
        for e in out:
            row = rk[rk["job_name"] == e["job"]]
            if not row.empty:
                r0 = row.iloc[0]
                e["mapped"] = {col: (str(r0[col]) if col in row.columns and pd.notna(r0[col]) else "")
                               for _, _, _, col in _ENGAGE_ROLES}
    return sorted(out, key=lambda r: r["rank"])


def _mmgbsa_components(mg: pd.DataFrame, frames) -> dict:
    """Median of each ΔG component, over all scored frames or over the reactive subset.

    The join is on the FRAME NUMBER, never on row position: the MM-GBSA CSV is written with a
    stride, so its row i is not frame i.
    """
    df = mg if frames is None else mg[mg["Frame"].isin(set(int(f) for f in frames))]
    out = {}
    for t in _MMGBSA_TERMS:
        col = f"r_psp_MMGBSA_dG_Bind_{t}"
        if col in df.columns and len(df):
            s = pd.to_numeric(df[col], errors="coerce").dropna()
            if not s.empty:
                out[t] = float(s.median())
    return out


def plot_mmgbsa_decomposition(out_dir: Path, ranks: list, merged: bool) -> None:
    """Energy components: the whole trajectory against the reactive (NAC) pose.

    The question is not 'what holds the ligand' but 'what CHANGES when the ligand reaches the
    reactive geometry', so the two bars per component are the same quantity over two ensembles —
    the faint bar is every scored frame, the solid bar only the reactive ones. A favourable Coulomb
    shift on reaching the NAC is electrostatic pre-organisation for the SN2.

    Each component carries a master container: symbolic, encoding nothing quantitative, it simply
    makes the pair read as one object. On a single-candidate panel the bars take their COMPONENT's
    colour (matching container and tick label); merged, the colour must separate the CANDIDATES.
    """
    _dpi = int(CFG.VIS_FIGURE_DPI)
    _ink = CFG.MMGBSA_INK
    _pal = list(CFG.MMGBSA_RANK_PALETTE)
    _f_leg = float(CFG.VIS_FONT_LEGEND)
    _min_kcal = float(CFG.DEFLUOR_COMPONENT_MIN_KCAL)
    cmap = plt.get_cmap("tab10")
    _scored = [r for r in ranks if not r["mmgbsa"].empty and "Frame" in r["mmgbsa"].columns]
    if not _scored:
        console_info("    [!] MM-GBSA decomposition skipped — no frame-stamped MM-GBSA CSV.")
        return

    for entry in ([None] if merged else _scored):
        rr = _scored if merged else [entry]
        fig, ax = plt.subplots(figsize=(14 if merged else 12, 6.6))
        """
        Drop the components that carry no signal — but decide that from the DATA, and only when the
        term is below CFG.DEFLUOR_COMPONENT_MIN_KCAL in every candidate AND both ensembles. A term
        that is zero for two candidates and non-zero for the third is a difference BETWEEN them and
        must stay. The omitted terms are named under the panel with their largest magnitude: an
        empty box is noise, but a silently deleted term is a lie.
        """
        keep, dropped = [], []
        for t in _MMGBSA_TERMS:
            mx = 0.0
            for r in _scored:                      # judge against every candidate, even here
                fr, _ = _reactive_frames(r["nac"])
                for vals in (_mmgbsa_components(r["mmgbsa"], None),
                             _mmgbsa_components(r["mmgbsa"], fr)):
                    v = vals.get(t)
                    if v is not None and v == v:
                        mx = max(mx, abs(v))
            (keep if mx >= _min_kcal else dropped).append((t, mx))
        terms = [t for t, _ in keep]
        if not terms:
            plt.close(fig)
            continue
        width = 0.78 / (len(rr) * 2)
        hdl = []

        for gi, r in enumerate(rr):
            fr, crit = _reactive_frames(r["nac"])
            allc = _mmgbsa_components(r["mmgbsa"], None)
            nacc = _mmgbsa_components(r["mmgbsa"], fr)
            n_scored = len(set(int(f) for f in fr) &
                           set(pd.to_numeric(r["mmgbsa"]["Frame"], errors="coerce")
                               .dropna().astype(int)))
            col = _CTRL_COLOUR if r.get("is_control") else _pal[gi % len(_pal)]
            bar_cols = ([cmap(i % 10) for i in range(len(terms))] if not merged else col)
            for si, (vals, alpha) in enumerate(((allc, 0.40), (nacc, 1.0))):
                off = (gi * 2 + si) * width - 0.39 + width / 2
                ax.bar(np.arange(len(terms)) + off,
                       [vals.get(t, np.nan) for t in terms], width=width * 0.92,
                       color=bar_cols, alpha=alpha, edgecolor=_ink["outline"], linewidth=0.6,
                       zorder=3)
            """
            One legend row per candidate, not two. The faint/solid pair means the same thing for
            every candidate, so it is stated ONCE as a neutral shade key. On a single-candidate
            panel the bars carry the component colours, so a coloured swatch would claim a meaning
            it does not have and the candidate is named by the text alone.
            """
            hdl.append(Patch(facecolor=col if merged else "none",
                             edgecolor=_ink["outline"] if merged else "none",
                             label=f"#{r['rank']} {r['ligand']} · n={n_scored:,}"))

        hdl += [Patch(facecolor=_ink["muted"], alpha=0.40, edgecolor=_ink["outline"],
                      label="whole traj."),
                Patch(facecolor=_ink["muted"], alpha=1.0, edgecolor=_ink["outline"],
                      label="reactive pose")]

        # The container is drawn AFTER the children so it can span the data, at a lower zorder so
        # the child bars sit inside it.
        _lo, _hi = ax.get_ylim()
        _pad = 0.06 * (_hi - _lo)
        for i, t in enumerate(terms):
            _master_container(ax, i, cmap(i % 10), width=0.94,
                              lo=_lo + _pad * 0.2, hi=_hi - _pad * 0.2)
        ax.set_ylim(_lo, _hi)

        ax.axhline(0.0, color=_ink["outline"], linewidth=0.9, zorder=2.5)
        ax.set_xticks(np.arange(len(terms)))
        ax.set_xticklabels([t.replace("_", " ") for t in terms], rotation=25, ha="right")
        for lbl, i in zip(ax.get_xticklabels(), range(len(terms))):
            lbl.set_color(cmap(i % 10))
            lbl.set_fontweight("bold")
        ax.set_ylabel("ΔG component (kcal/mol)")
        hdl.append(Patch(facecolor="none", edgecolor="none", label="negative = favours binding"))
        """
        The omitted terms are NAMED, but in the legend rather than under the axis: a term dropped for
        being numerically dead is a footnote and belongs at footnote size. Nothing is hidden — an
        empty box is noise, a silently deleted term is a lie.
        """
        if dropped:
            hdl.append(Patch(facecolor="none", edgecolor="none",
                             label="omitted (|ΔG| < %.1f): %s"
                                   % (_min_kcal, ", ".join(t.replace("_", " ")
                                                           for t, _ in dropped))))
        # Tick every DEFLUOR_COMPONENT_TICK_KCAL: matplotlib's default step is coarser than most of
        # the components themselves — an H-bond term of −2.7 cannot be read off a 20 kcal/mol grid.
        _tick = float(CFG.DEFLUOR_COMPONENT_TICK_KCAL)
        ax.yaxis.set_major_locator(MultipleLocator(_tick))
        ax.yaxis.set_minor_locator(MultipleLocator(_tick / 2))
        ax.grid(alpha=CFG.VIS_GRID_ALPHA, linewidth=CFG.VIS_GRID_LINEWIDTH, axis="y")
        ax.grid(alpha=CFG.VIS_GRID_ALPHA * 0.5, linewidth=CFG.VIS_GRID_LINEWIDTH * 0.7,
                axis="y", which="minor")
        ax.set_axisbelow(True)
        ax.legend(handles=hdl, loc="upper right", fontsize=_f_leg, frameon=True)

        out_path = (out_dir / "06_MMGBSA_Decomposition_AllRanks.png" if merged
                    else rr[0]["dir"] / f"{rr[0]['job']}_MMGBSA_NAC_Decomposition.png")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            plt.savefig(out_path, dpi=_dpi, bbox_inches="tight")
        plt.close(fig)
        console_info(f"    MM-GBSA decomposition saved : {out_path.name}")


def plot_machinery_engagement(out_dir: Path, ranks: list, merged: bool) -> None:
    """Which catalytic residue is actually in reach of the warhead carbon, judged by the criteria
    that define the interaction.

    A distance alone says nothing; a distance against the CFG criterion that defines the
    interaction says everything, so the cut-offs are drawn as reference bands. Each bar is the
    MEDIAN over the reactive frames with the interquartile range as the whisker, and each carries
    the share of those frames in which the contact held — the median cannot say how PERSISTENT a
    contact is, and two residues with the same median can differ many-fold in occupancy.
    """
    _dpi = int(CFG.VIS_FIGURE_DPI)
    _ink = CFG.MMGBSA_INK
    _pal = list(CFG.MMGBSA_RANK_PALETTE)
    _f_leg, _f_ann = float(CFG.VIS_FONT_LEGEND), float(CFG.VIS_FONT_ANNOT)
    _roles = CFG.ACTIVE_SITE_ROLE_GROUP_COLOUR
    _bands = list(CFG.DEFLUOR_BAND_COLOURS)
    _occ_cut = float(CFG.THRESHOLD_SALT_BRIDGE)
    _zones = _engage_zones()
    _ranks = [r for r in ranks if len(_reactive_frames(r["nac"])[0])]
    if not _ranks:
        console_info("    [!] Machinery engagement skipped — no reactive frames.")
        return

    for entry in ([None] if merged else _ranks):
        rr = _ranks if merged else [entry]
        fig, ax = plt.subplots(figsize=(13.5 if merged else 11, 6.6))
        width = 0.8 / len(rr)
        rank_hdl = []

        """
        The ceiling comes from the DATA — the tallest upper quartile across every candidate and
        residue — and is fixed BEFORE anything is drawn, since the containers and the bands are
        sized against it. Every candidate shares one ceiling, so the panels stay comparable.
        """
        _q3max = 0.0
        for _r in _ranks:
            _frx, _ = _reactive_frames(_r["nac"])
            _sx = _r["nac"][_r["nac"]["Frame"].isin(set(int(f) for f in _frx))]
            for _c, _, _, _ in _ENGAGE_ROLES:
                if _c in _sx:
                    _v = pd.to_numeric(_sx[_c], errors="coerce").dropna()
                    if not _v.empty:
                        _q3max = max(_q3max, float(_v.quantile(.75)))
        _ytop = max(float(CFG.DEFLUOR_ENGAGE_Y_MIN_TOP), float(np.ceil(_q3max + 0.6)))

        """
        The criteria as reference bands, drawn behind everything. Each band is tagged with a LETTER
        in the margin OUTSIDE the axes — written inside, the tag lands on whichever bar happens to
        reach that height — and spelled out in the legend. Each letter takes its own band's colour,
        darkened to be legible.
        """
        prev = 0.0
        band_hdl = []
        _blend = ax.get_yaxis_transform()          # x in axes fraction, y in data units
        for i, (cut, lab) in enumerate(_zones):
            tag = chr(ord("A") + i)
            c = _bands[i % len(_bands)]
            ax.axhspan(prev, cut, color=c, alpha=float(CFG.DEFLUOR_BAND_ALPHA), zorder=0)
            ax.axhline(cut, color=_darken(c, 0.75), linestyle="--", linewidth=1.0, alpha=0.85,
                       zorder=1)
            ax.annotate(tag, xy=(1.0, (prev + cut) / 2), xycoords=_blend, xytext=(7, 0),
                        textcoords="offset points", fontsize=_f_ann + 2.0, fontweight="bold",
                        color=_darken(c, 0.62), va="center", ha="left", zorder=6,
                        annotation_clip=False)
            band_hdl.append(Patch(facecolor=c, alpha=0.45, edgecolor=_darken(c, 0.75),
                                  label=f"{tag} · {lab} ({cut:g} Å)"))
            prev = cut

        for gi, r in enumerate(rr):
            fr, crit = _reactive_frames(r["nac"])
            sub = r["nac"][r["nac"]["Frame"].isin(set(int(f) for f in fr))]
            """
            On a single-candidate panel the bar takes its ROLE's colour — the same colour as the
            container around it and the tick label beneath it, so the three read as one object.
            Merged, the colour must separate the CANDIDATES, since the role is already given by the
            container they share.
            """
            col = _CTRL_COLOUR if r.get("is_control") else _pal[gi % len(_pal)]
            bar_cols = ([_roles.get(role, col) for _, role, _, _ in _ENGAGE_ROLES]
                        if not merged else col)
            meds, occs, q1s, q3s = [], [], [], []
            for c, _, _, _ in _ENGAGE_ROLES:
                s = (pd.to_numeric(sub[c], errors="coerce").dropna()
                     if c in sub else pd.Series(dtype=float))
                meds.append(float(s.median()) if not s.empty else np.nan)
                q1s.append(float(s.quantile(.25)) if not s.empty else np.nan)
                q3s.append(float(s.quantile(.75)) if not s.empty else np.nan)
                occs.append(100.0 * float((s <= _occ_cut).mean()) if not s.empty else np.nan)
            off = gi * width - 0.4 + width / 2
            err = np.array([[m - a for m, a in zip(meds, q1s)],
                            [b - m for m, b in zip(meds, q3s)]])
            _lab = f"#{r['rank']} {r['ligand']} · {crit.replace('geometric', 'geom.')} · n={len(fr):,}"
            ax.bar(np.arange(len(_ENGAGE_ROLES)) + off, meds, width=width * 0.9,
                   color=bar_cols, alpha=0.9, edgecolor=_ink["outline"], linewidth=0.6,
                   yerr=err, capsize=3, error_kw=dict(ecolor=_ink["muted"], lw=0.9), zorder=3,
                   label=(_lab if merged else None))
            if not merged:
                rank_hdl = [Patch(facecolor="none", edgecolor="none", label=_lab)]

            for i, m in enumerate(meds):
                if m != m:
                    continue
                """
                Two different quantities, kept apart so they cannot be read as one: the MEDIAN
                carries its unit and sits above the whisker cap (on the bar top a long whisker
                crosses it), the OCCUPANCY sits inside the bar head, where the bar's own colour says
                whose it is. The residue name holds the foot of the bar, so the head is free.
                """
                _y = max(m, q3s[i] if q3s[i] == q3s[i] else m)
                _c = _ink["dark"] if not merged else _darken(col, 0.55)
                ax.text(i + off, min(_y + 0.10, _ytop - 0.25), f"{m:.1f} Å", ha="center",
                        va="bottom", fontsize=_f_ann, fontweight="bold", color=_c, zorder=10,
                        path_effects=[pe.withStroke(linewidth=2.4, foreground=_ink["light"])])
                _bar_c = bar_cols[i] if isinstance(bar_cols, list) else bar_cols
                _on_bar = auto_label_colour(CFG, _bar_c)
                ax.text(i + off, m - 0.13, f"{occs[i]:.0f}%", rotation=90, ha="center", va="top",
                        fontsize=_f_ann - 0.5 if merged else _f_ann, fontweight="bold",
                        color=_on_bar, zorder=10)

            # Every bar names the residue it measures, written up the inside of the bar: the residue
            # behind a role is per-homolog (ASP110 here, ASP109 there), so it belongs on the bar and
            # not on a tick shared between candidates.
            mp = r.get("mapped", {})
            for i, (_, _, _, mcol) in enumerate(_ENGAGE_ROLES):
                resn = mp.get(mcol, "")
                if resn and meds[i] == meds[i]:
                    _bar_c = bar_cols[i] if isinstance(bar_cols, list) else bar_cols
                    ax.text(i + off, 0.12, resn, rotation=90, ha="center", va="bottom",
                            fontsize=_f_ann - 0.7 if merged else _f_ann - 0.3, fontweight="bold",
                            color=auto_label_colour(CFG, _bar_c), zorder=8)

        # Master containers group the residues by CATALYTIC ROLE (colours from CFG). The outline is
        # kept a hair inside the axis; flush to the limit it merges with the spine and reads as if
        # the box had burst through it.
        _ylo, _yhi = 0.0, _ytop * 0.99
        start = 0
        for i in range(len(_ENGAGE_ROLES) + 1):
            if i == len(_ENGAGE_ROLES) or _ENGAGE_ROLES[i][1] != _ENGAGE_ROLES[start][1]:
                c = _roles.get(_ENGAGE_ROLES[start][1], _ink["muted"])
                x0, x1 = start - 0.47, (i - 1) + 0.47
                _master_container(ax, (x0 + x1) / 2, c, width=(x1 - x0), lo=_ylo, hi=_yhi)
                start = i
        ax.set_ylim(0, _ytop)

        ax.set_xticks(np.arange(len(_ENGAGE_ROLES)))
        ax.set_xticklabels([lbl for _, _, lbl, _ in _ENGAGE_ROLES])
        for lb, (_, role, _, _) in zip(ax.get_xticklabels(), _ENGAGE_ROLES):
            lb.set_color(_roles.get(role, _ink["muted"]))
            lb.set_fontweight("bold")
        ax.set_ylabel("Distance to warhead C (Å)")
        _hdl_stats = Patch(facecolor="none", edgecolor="none",
                           label="bars = median · whiskers = IQR")
        ax.yaxis.set_major_locator(MultipleLocator(1.0))
        ax.yaxis.set_minor_locator(MultipleLocator(0.5))
        ax.grid(alpha=CFG.VIS_GRID_ALPHA, linewidth=CFG.VIS_GRID_LINEWIDTH, axis="y")
        ax.set_axisbelow(True)

        _hdl = (ax.get_legend_handles_labels()[0] if merged else rank_hdl)
        _hdl += [_hdl_stats,
                 Patch(facecolor="none", edgecolor="none", label="Å above bar = median distance"),
                 Patch(facecolor="none", edgecolor="none",
                       label=f"% in bar = frames within {_occ_cut:g} Å")] + band_hdl
        ax.legend(handles=_hdl, loc="upper right", fontsize=_f_leg, frameon=True)

        out_path = (out_dir / "07_Machinery_Engagement_AllRanks.png" if merged
                    else rr[0]["dir"] / f"{rr[0]['job']}_Machinery_Engagement.png")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            plt.savefig(out_path, dpi=_dpi, bbox_inches="tight")
        plt.close(fig)
        console_info(f"    Machinery engagement saved : {out_path.name}")


def generate_reactive_pose_figures(out_dir: Path, df_master: pd.DataFrame) -> None:
    """Both reactive-pose figures, per candidate and merged, from the per-frame tables on disk."""
    apply_figure_style(CFG)
    ranks = _load_reactive_pose_data(out_dir)
    if not ranks:
        console_info("    [!] Reactive-pose figures skipped — no per-frame NAC tables found.")
        return
    for _merged in (False, True):
        plot_mmgbsa_decomposition(out_dir, ranks, merged=_merged)
        plot_machinery_engagement(out_dir, ranks, merged=_merged)


# =============================================================================
# SECTION 6: PHYSICAL CHEMISTRY MODULES
# =============================================================================

_MIC_WARNED = False   # one-shot guard: warn (not silently pass) if a blockade MIC correction fails

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
        except Exception as _e:
            # Discard this frame rather than continue with uncorrected (non-PBC)
            # vectors, which would corrupt the blockade metric near box edges. The
            # NaN score never wins `score > best_score`, so the frame is dropped
            # from QM/MM selection; if the whole box is malformed every frame is
            # dropped and the rank skips gracefully (ideal_frame_idx stays -1).
            global _MIC_WARNED
            if not _MIC_WARNED:
                console_info(f"    {ConsoleColours.WARNING}[!] Water-blockade MIC correction "
                             f"failed ({_e}); discarding affected frame(s) — check the "
                             f"trajectory box tensor.{ConsoleColours.ENDC}")
                _MIC_WARNED = True
            return float("nan")

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
        """
        The weight is a Boltzmann-like occupancy of the WaterMap site, so the exponent must be
        DIMENSIONLESS. exp(dG) with dG in kcal/mol is not: it implicitly divides by 1 kcal/mol rather than
        by RT, which at 300 K would mean an effective temperature near 503 K, and it silently rescales with
        any change of energy unit.

        dG/RT is dimensionless, and RT is the same RT the rest of the pipeline uses (CFG, 300 K). The
        exponent is clipped because a strongly stabilised site can otherwise overflow exp().
        """
        _rt_wm = float(CFG.GAS_CONSTANT_KCAL) * float(CFG.MMGBSA_TEMPERATURE_K)
        _z_wm = float(np.clip(best_dg / _rt_wm, -60.0, 60.0))
        weight = (1.0 / (1.0 + np.exp(_z_wm)) * 2.0) if best_d < match_r else 1.0
        total += weight
    return total




# =============================================================================
# SECTION 7: QSITE AUTOMATION
# =============================================================================

# Solvation-droplet radius (Å) for the QSite .mae: only solvent within this
# distance of the ligand is retained, trimming the full periodic box to a
# tractable local MM region. Defined here because write_qsite_droplet() takes it
# as a default argument, evaluated when the function is defined below.
_QSITE_DROPLET_RADIUS = float(CFG.QSITE_DROPLET_RADIUS)


# -----------------------------------------------------------------------------
# SECTION 7a: QM/MM INPUT PREPARATION & EXECUTION
# Trims the periodic box to a solvation droplet, writes the QSite/Jaguar SN2
# relaxed-scan .in, and launches the QM/MM job (idempotent per output folder).
# -----------------------------------------------------------------------------
def write_qsite_droplet(cms_model, path: Path, lig_resname: str,
                        radius: float = _QSITE_DROPLET_RADIUS) -> str:
    """Write an uncompressed QSite .mae trimmed to a solvation droplet: the full
    protein + ligand + only solvent molecules with an atom within `radius` Å of
    any ligand atom. Trims the full periodic water box to a tractable local MM
    region for QM/MM. Falls back to writing the full structure on any error, so
    QSite always receives a valid input. Returns a short status string for logs.
    """
    try:
        st = cms_model.fsys_ct.copy()
        lig_xyz = np.array([a.xyz for a in st.atom
                            if a.pdbres.strip() == lig_resname])
        if lig_xyz.size == 0:
            st.write(str(path))
            return "full (no ligand atoms found — not trimmed)"
        '''
        Counter-ions are trimmed on the SAME radius as the water (CFG.COUNTERION_RESTYPES). Trimming
        only the water would leave every Na+/Cl- in the box behind with its hydration shell deleted:
        a bare point charge in the frozen MM shell, screened by vacuum instead of by bulk water, and
        polarising the QM Hamiltonian from wherever it happened to diffuse. Ions inside the droplet
        keep their water and are kept with it.

        Distances go through the minimum-image convention. In the normal path this changes nothing —
        center_cms has already wrapped every atom into a primary cell centred on the ligand, so each
        atom lies within half a box length of it and the minimum image IS the atom itself. That is a
        property of the centring, not of this function, and a droplet cut on raw Cartesian distances
        would silently lose a solvation lobe the moment the centring did not hold. The convention is
        applied here so the trim is correct on its own terms.
        '''
        _trim_res = ({r.strip().upper() for r in CFG.SOLVENT_RESTYPES} |
                     {r.strip().upper() for r in CFG.COUNTERION_RESTYPES})
        _box = getattr(cms_model, "box", None)
        _mol_atoms = defaultdict(list)
        for a in st.atom:
            if a.pdbres.strip().upper() in _trim_res:
                _mol_atoms[a.molecule_number].append(a.index)
        _del, _n_ion_cut = [], 0
        for _mol, _aidxs in _mol_atoms.items():
            _coords = np.array([st.atom[i].xyz for i in _aidxs])
            _dmin = float(np.min(_mic_dists_2d(_coords, lig_xyz, _box)))
            if _dmin > radius:
                _del.extend(_aidxs)
                if st.atom[_aidxs[0]].pdbres.strip().upper() not in {
                        r.strip().upper() for r in CFG.SOLVENT_RESTYPES}:
                    _n_ion_cut += 1
        if _del:
            st.deleteAtoms(_del)
        """
        The droplet's BOUNDARY. Cutting the periodic box to a finite ball of water leaves a free
        surface — nothing holds the outermost molecules, and they relax into vacuum during the
        relaxed scan, distorting the electrostatics the QM region sits in.

        QSite takes the boundary from a per-atom property, `i_i_constraint`, and Jaguar reports back
        how many atoms it accepted as frozen and as constrained (verified against a live QM/MM job:
        the values map 0 = free, 1 = FROZEN, 2 = CONSTRAINED). Three zones around the ligand:

            ≤ QSITE_FREE_RADIUS     free       — the reaction centre and its first solvation shell
            ≤ QSITE_BUFFER_RADIUS   restrained — can respond to the reaction, cannot drift
            beyond                  frozen     — the droplet surface, held rigid

        The protein is never frozen inside the free/buffer radii, so the catalytic machinery keeps
        every degree of freedom the chemistry needs.
        """
        _free_r = float(CFG.QSITE_FREE_RADIUS)
        _buf_r = float(CFG.QSITE_BUFFER_RADIUS)
        _n_free = _n_con = _n_frz = 0
        _all_xyz = np.array([a.xyz for a in st.atom], dtype=float)
        _d_all = _mic_dists_2d(_all_xyz, lig_xyz, _box).min(axis=1)
        for a, _d in zip(st.atom, _d_all):
            _d = float(_d)
            if _d <= _free_r:
                _v = 0
            elif _d <= _buf_r:
                _v = 2
            else:
                _v = 1
            a.property['i_i_constraint'] = _v
            _n_free += _v == 0
            _n_con += _v == 2
            _n_frz += _v == 1
        st.write(str(path))
        return (f"droplet r={radius:.1f} Å (removed {len(_del)} solvent/ion atoms, "
                f"{_n_ion_cut} counter-ion(s); boundary: "
                f"{_n_free} free ≤{_free_r:g} Å, {_n_con} restrained ≤{_buf_r:g} Å, "
                f"{_n_frz} frozen beyond)")
    except Exception as exc:
        cms_model.fsys_ct.write(str(path))
        return f"full (droplet trim failed: {exc})"


def generate_qsite_inputs(mae_path: Path, job_name: str,
                          nuc_num, stab_f_num, lig_c_idx, nuc_o_idx,
                          base_num=None, acid_num=None,
                          lig_charge: int = None, lig_resname: str = "LIG",
                          cradle_nums=None) -> Path:
    """
    Write a valid QSite/Jaguar QM/MM relaxed-scan .in for the SN2
    dehalogenation reaction coordinate (Nu_O···C_lig distance scan).

    The format follows the genuine Jaguar/QSite input specification — verified
    against `schrodinger.application.qsite.input.QSiteInput` and the Jaguar
    `scan-relaxed` example:
      MAEFILE      — the structure (uncompressed .mae, same folder)
      &gen         — DFT functional / basis / QM charge / QM-MM mode / relaxed opt
      &qmregion    — ligand molecule (QM) + catalytic-sidechain QM/MM cuts
      &zvar/&coord — relaxed scan of the Nu_O–C_lig distance (start → product)

    QM region = LIG substrate + Asp110 (Nuc) + His280 (Base) + Asp134 (Acid)
    + His155 (StabH) sidechains + the fluoride-cradle aromatic/H-bond donors
    (TRP/TYR/HIS, `cradle_nums`) + up to `_QM_WATER_MAX` coordinating waters
    within `_QM_WATER_RADIUS` of the reaction centre. Excluding Base/Acid from QM
    would push the proton-transfer half of the mechanism onto the MM force field;
    excluding the cradle donors / first-shell waters would leave the F⁻ leaving
    group under-stabilised and bias the QM/MM barrier upward.

    NOTE: no implicit-solvation keyword is written. The extracted frame retains
    its explicit TIP3P water box in the MM region, so adding an implicit model
    (an implicit `SOLVATION_METHOD sgb`) would double-count solvation. Strip the box
    and re-add `isolv` only for an implicit-solvent QM/MM variant.
    """
    inp_path = mae_path.parent / f"{job_name}_QSite_SN2.in"

    # ── Read structure (needed for charge, ligand molid, and cut resolution) ──
    # A truncated or empty .mae yields no structure: raise the reason rather than a bare
    # StopIteration from inside a generator, and close the handle either way.
    with structure.StructureReader(str(mae_path)) as _r:
        st = next(iter(_r), None)
    if st is None:
        raise ValueError(f"QSite input structure is empty or unreadable: {mae_path}")
    lig_mol = next(
        (m.number for m in st.molecule
         if any(a.pdbres.strip() == lig_resname for a in m.atom)), None)
    if lig_mol is None:
        # The ligand IS the QM region's reactive core and its charge is always added to molchg; if it
        # is not found, omitting it while still counting its charge gives Jaguar an electron-count
        # mismatch that aborts the scan. Fail loudly instead.
        raise ValueError(f"QSite QM region: ligand '{lig_resname}' not found in {mae_path.name} — "
                         "cannot build a charge-consistent QM region.")

    # Each QM/MM cut row must name the protein molecule id (QSite matches the
    # residue within that molecule). Resolve it via the residue's backbone Cα,
    # which also disambiguates the catalytic residue from any water that happens
    # to share the same residue number. Returns (molid, chain) or None.
    def _resolve_cut(resnum):
        """
        Resolve a residue's QM/MM frozen-orbital boundary as a list of (QM atom, MM atom) bonds.

        The standard sidechain cut is Cα–Cβ: the sidechain is QM, the backbone stays MM. Two residue
        types cannot be cut that way, and both are handled with their own boundary rather than being
        dropped from the QM region — an alignment-mapped catalytic role CAN land on either, and a
        silently missing cradle residue changes the barrier being computed.

          GLYCINE has no Cβ. Its only 'sidechain' is the second α-hydrogen, so the boundary has to be
          taken across the backbone instead: a DOUBLE cut at N–Cα and C–Cα leaves Cα and its hydrogens
          in the QM region with the peptide N and the carbonyl C in MM.

          PROLINE has a Cβ, so the glycine guard never sees it, but its sidechain closes back onto the
          backbone nitrogen. A lone Cα–Cβ cut severs the pyrrolidine ring and leaves the Cδ–N bond
          crossing the boundary uncapped. The ring is closed with a SECOND cut at Cδ–N, so the whole
          C₃ bridge is QM and both of its attachments to the backbone are proper frozen-orbital cuts.

        Matching is on the residue's Cα, which also disambiguates the catalytic residue from any water
        sharing its residue number. A blank chain id is accepted — a single-chain MD frame carries one,
        and requiring a chain would drop every catalytic residue from the QM region.
        """
        a = next((a for a in st.atom
                  if a.resnum == resnum and a.pdbname.strip() == "CA"), None)
        if a is None:
            return None
        _rn = next((b.pdbres.strip() for b in st.atom if b.resnum == resnum), "?")
        _names = {b.pdbname.strip() for b in st.atom if b.resnum == resnum
                  and b.molecule_number == a.molecule_number and b.chain == a.chain}

        if _rn.upper().startswith("PRO"):
            if {"CB", "CD", "N"} <= _names:
                return (a.molecule_number, a.chain, [("CB", "CA"), ("CD", "N")])
            console_info(f"    [!] Residue {_rn}{resnum} is a proline with an incomplete ring "
                         f"(missing {sorted({'CB', 'CD', 'N'} - _names)}) — dropped from the QM region.")
            return None

        if "CB" not in _names:                      # glycine (or any Cβ-less residue)
            if {"N", "C"} <= _names:
                console_info(f"    [i] Residue {_rn}{resnum} has no Cβ — QM boundary taken across the "
                             f"backbone instead (double cut N–Cα, C–Cα).")
                return (a.molecule_number, a.chain, [("CA", "N"), ("CA", "C")])
            console_info(f"    [!] Residue {_rn}{resnum} has neither Cβ nor a complete backbone — "
                         f"dropped from the QM region.")
            return None

        return (a.molecule_number, a.chain, [("CB", "CA")])

    # ── Resolve every catalytic cut once → (resnum, molid, chain). The residue
    #    identities (nuc/base/acid/stab) are alignment- and Smart-Lock-derived
    #    upstream and passed in, never hardcoded, so the QM region adapts to each
    #    homolog's own numbering. Reused for both the charge sum and the cut table.
    _resolved_cuts = []
    _seen = set()
    for rn in (nuc_num, base_num, acid_num, stab_f_num, *(cradle_nums or [])):
        if not rn or rn in _seen:        # de-dup: Base/Acid/cradle may share a residue
            continue
        _seen.add(rn)
        _r = _resolve_cut(rn)
        if _r is not None:
            _resolved_cuts.append((rn, _r[0], _r[1], _r[2]))   # resnum, molid, chain, [(qm, mm), …]

    # Cap the QM region size. The list is priority-ordered (nucleophile, base, acid,
    # stab-F, then cradle), so truncating keeps the catalytic core and drops the
    # furthest cradle residues. A very large QM region (e.g. 17 residues) inflates the
    # electron count and makes a molchg / electron-parity mismatch likely — which
    # Jaguar rejects with 'incorrect molecular charge' and then skips every scan point.
    _max_qm = int(getattr(CFG, "QSITE_MAX_QM_RESIDUES", 0))
    if _max_qm > 0 and len(_resolved_cuts) > _max_qm:
        _dropped = len(_resolved_cuts) - _max_qm
        _resolved_cuts = _resolved_cuts[:_max_qm]
        console_info(f"    [QSite] QM region capped to {_max_qm} residues "
                     f"(dropped {_dropped} furthest cradle residue(s)) to keep the "
                     f"charge/electron count tractable.")

    # ── QM region net charge — derived from the structure, not assumed ─────
    # QM atoms = full ligand + each catalytic residue's sidechain beyond the
    # Cα–Cβ cut (CA/backbone remain MM). Instead of assuming protonation states
    # (a fixed `lig − 1` would count only the nucleophile and presume a
    # protonated acid and neutral His), sum the *actual* formal charges those
    # atoms carry in the prepared structure. This keeps `molchg` correct for
    # whatever PrepWizard assigned (deprotonated Asp −1, neutral His, protonated
    # catalytic acid 0, HIP +1, …) and for any residue numbering, so QSite/Jaguar
    # does not abort on a QM charge/electron-count mismatch.
    # Include terminal backbone atoms (OXT + N-terminal H1/H2/H3): if a catalytic /
    # QM residue sits at a chain terminus, OXT would otherwise be counted as a
    # sidechain atom and add a spurious −1 to qm_charge, while the QM/MM cut at CB–CA
    # leaves OXT in the MM region → a charge/electron-count mismatch that aborts QSite.
    _backbone = {'N', 'C', 'CA', 'O', 'H', 'HA', 'OXT', 'H1', 'H2', 'H3'}   # Smart-Lock convention (L429) + termini

    def _sidechain_formal_charge(molid, chain, resnum, cut_pairs):
        """
        The QM atoms of this residue are the ones on the QM side of its own boundary, so the charge
        sum has to follow the cut that was actually made — not a fixed 'everything but the backbone'
        rule. For a glycine cut across the backbone the QM side is Cα and its hydrogens; for every
        other residue it is the sidechain beyond Cβ (proline included, whose C₃ bridge is all
        sidechain-named). Both come out formally neutral, but deriving it from the cut keeps the
        molchg correct if a future boundary ever encloses a charged backbone atom.
        """
        _qm_side = {q for q, _m in cut_pairs}
        if _qm_side == {"CA"}:                      # backbone double cut (glycine): QM = Cα + its H
            return sum(
                a.formal_charge for a in st.atom
                if a.resnum == resnum and a.chain == chain
                and a.molecule_number == molid
                and a.pdbname.strip() in ("CA", "HA", "HA2", "HA3")
            )
        return sum(
            a.formal_charge for a in st.atom
            if a.resnum == resnum and a.chain == chain
            and a.molecule_number == molid
            and a.pdbname.strip() not in _backbone
        )

    # Ligand charge stays authoritative from CFG (curated per PFAS); only the
    # protein-sidechain contribution is read dynamically from the structure.
    if lig_charge is None:
        jn_upper = job_name.upper()
        lig_charge = next(
            (v for k, v in CFG.LIGAND_QM_CHARGES.items() if k.upper() in jn_upper),
            CFG.QSITE_CHARGE,
        )
    _sidechain_charge = sum(
        _sidechain_formal_charge(molid, chain, rn, cut_pairs)
        for rn, molid, chain, cut_pairs in _resolved_cuts
    )
    qm_charge = int(round(lig_charge + _sidechain_charge))

    # ── &qmregion QM/MM cut table ──────────────────────────────────────────
    # One row per boundary bond: Cα–Cβ for a standard sidechain, two rows for a glycine
    # (backbone N–Cα and C–Cα) or a proline (Cα–Cβ plus the ring-closing Cδ–N).
    _cuts = [
        f"  {molid:>5}  {chain:>4}  {rn:>6}  {_qm:>7}  {_mm:>7}"
        for rn, molid, chain, cut_pairs in _resolved_cuts
        for _qm, _mm in cut_pairs
    ]

    # ── Coordinating catalytic waters → whole-molecule QM ──────────────────────
    # First-shell waters within _QM_WATER_RADIUS of the reaction centre (scissile
    # C, leaving F, nucleophile O) stabilise the departing fluoride; add the
    # nearest few as full QM molecules so the leaving-group energy is not left to
    # the MM force field. Self-contained on the structure's coordinates.
    _qm_water_mols = []
    try:
        _centres = []
        for _ai in (nuc_o_idx, lig_c_idx):
            if _ai and 1 <= int(_ai) <= st.atom_total:
                _centres.append(np.array(st.atom[int(_ai)].xyz))
        for a in st.atom:
            if a.pdbres.strip() == lig_resname and (a.element or "").strip() == "F":
                _centres.append(np.array(a.xyz))
        _water_res = {r.strip().upper() for r in CFG.SOLVENT_RESTYPES}
        if _centres:
            _cen_arr = np.asarray(_centres)
            _cand = []
            for a in st.atom:
                if (a.element or "").strip() == "O" and a.pdbres.strip().upper() in _water_res:
                    _d = float(np.min(np.linalg.norm(_cen_arr - np.array(a.xyz), axis=1)))
                    if _d <= _QM_WATER_RADIUS:
                        _cand.append((_d, a.molecule_number))
            for _d, _mol in sorted(_cand):
                if _mol != lig_mol and _mol not in _qm_water_mols:
                    _qm_water_mols.append(_mol)
                if len(_qm_water_mols) >= _QM_WATER_MAX:
                    break
    except Exception as _wexc:
        # An empty QM water shell CHANGES the physics (the fluoride's first solvation shell goes
        # classical), so it is never a silent fallback.
        console_info(f"    [!] QM water selection failed ({_wexc}) — the QM region will carry no "
                     f"explicit waters.")
        _qm_water_mols = []

    # ── Relaxed scan: Nu_O–C_lig distance, NAC start → product ────────────
    _start = CFG.QSITE_SCAN_START
    _end   = CFG.QSITE_SCAN_START + CFG.QSITE_SCAN_STEP * (CFG.QSITE_SCAN_NSTEPS - 1)

    # SCF robustness for the QM/MM relaxed scan. The diffuse 6-31+G** basis on a large QM
    # region causes near-linear-dependence ("small singular value") and DIIS blow-ups that
    # abort scan points (observed: ~half the points fatal). Use a non-diffuse basis for the
    # scan geometry (diffuse adds little to a RELATIVE barrier), a level shift + raised
    # iteration cap to force convergence, and nofail so one hard point can't kill the scan.
    _scan_basis = CFG.QSITE_SCAN_BASIS.replace("(d,p)", "**").replace("(d)", "*")
    _gen = [
        f"basis={_scan_basis}",
        "igeopt=1",                       # relaxed (constrained-optimised) scan
        f"molchg={qm_charge}",            # QM-region net charge (authoritative — see per-step patch in run_qsite)
        f"dftname={CFG.QSITE_FUNCTIONAL}",
        "mmqm=1",                         # enable QM/MM
        "impversion=huge",
        f"vshift={CFG.QSITE_SCF_VSHIFT:g}",    # SCF level shift (stabilises convergence)
        f"maxit={int(CFG.QSITE_SCF_MAXIT)}",   # max SCF iterations
        f"iacc={int(CFG.QSITE_SCF_IACC)}",     # SCF accuracy grid (1 = robust/fast)
        "nofail=1",                                            # a non-converged point is skipped, not fatal
        "mulken=1",                                            # print the Mulliken population analysis:
                                                               # without it Jaguar writes NO charge table
                                                               # and the departing-fluoride charge — the
                                                               # electronic proof of C–F cleavage — cannot
                                                               # be parsed at all (verified 13 July 2026;
                                                               # 'mulliken' and 'ipop' are rejected)
    ]
    if CFG.QSITE_MULT != 1:
        _gen.append(f"multip={CFG.QSITE_MULT}")

    """
    &mmkey — the MM half of the QM/MM Hamiltonian.

    The classical region must not be scored on a different force field from the trajectory it came
    from. The Desmond system was built and propagated under OPLS4, so the MM environment around the
    QM region should be OPLS4 too — but QSite does not offer it. Schrodinger's own QSite driver
    exposes exactly two (-ffield choices=['OPLS_2005','OPLS3e']) and emits the string keyword
    qsite_ff='opls3e' for the newer one. OPLS_2005 is the DEFAULT, so an unset &mmkey scores an OPLS4
    trajectory on a force field two generations older.

CFG.QSITE_MM_FF is therefore EMPTY, so the &mmkey is empty and Impact falls back to its default
    OPLS_2005: opls3e / S-OPLS cannot be used because it aborts every QSite job at Impact line 22 on a
    frozen-orbital cut (00_01 §QSITE_MM_FF). This is a DECLARED LIMITATION, not a fix — OPLS_2005 is a
    generation below the trajectory's OPLS4 — but it is the only MM force field QSite runs with these
    cuts. The QM region, where the bond breaks, is unaffected; the mismatch is in the classical
    environment around it.
    """
    _mmkey = f"&mmkey\n{CFG.QSITE_MM_FF}\n&\n" if getattr(CFG, "QSITE_MM_FF", "") else "&mmkey\n&\n"
    if getattr(CFG, "QSITE_MM_FF", ""):
        console_info(f"    [i] QSite &mmkey MM force field: {CFG.QSITE_MM_FF} — the closest QSite offers to "
                     f"the trajectory's OPLS4 (QSite supports only OPLS_2005 / OPLS3e).")
    else:
        console_info("    [!] QSite &mmkey carries no force-field flag — the MM region falls back to "
                     "OPLS_2005 while the trajectory was propagated under OPLS4.")

    content = (
        f"MAEFILE: {mae_path.name}\n"
        "&gen\n" + "\n".join(_gen) + "\n&\n"
        + _mmkey
        + "&qmregion\n"
        " molid chain  resnum   qmatom   mmatom\n"
        + ("\n".join(_cuts) + "\n" if _cuts else "")
        + " molid theory\n"
        + (f"     {lig_mol}     qm\n" if lig_mol is not None else "")   # molid 0 is falsy but valid
        + "".join(f"     {wm}     qm\n" for wm in _qm_water_mols)
        + "&\n"
        "&zvar\n"
        f"r = {_start:g} to {_end:g} in {CFG.QSITE_SCAN_NSTEPS}\n"
        "&\n"
        "&coord\n"
        f" {nuc_o_idx} {lig_c_idx} # r\n"
        "&\n"
    )
    with open(inp_path, "w") as f:
        f.write(content)
    return inp_path


def run_qsite(qsite_dir: Path, inp_path: Path, job_name: str, rank: int) -> bool:
    """
    Launch the QSite/Jaguar executable on a generated .in, writing all output
    inside `qsite_dir`. Returns True if a job was launched, False if QSite is
    unavailable.

    Run policy mirrors the PDB-preparation cache: the caller skips the whole
    step when `qsite_dir` already exists, so this only runs for fresh folders.
    `subprocess` is given an explicit `cwd`, which is per-process and therefore
    thread-safe (unlike `os.chdir`) under the worker pool.
    """
    qsite_exe = os.path.join(os.environ.get("SCHRODINGER", ""), "qsite")
    if not os.path.isfile(qsite_exe):
        console_info(f"    {ConsoleColours.WARNING}[Rank {rank}] qsite executable not found "
                     f"at {qsite_exe} — skipping launch (input written).{ConsoleColours.ENDC}")
        return False

    jobname = f"{job_name}_QSite_SN2"
    cmd = [qsite_exe, "-WAIT", "-PARALLEL", str(_QSITE_PROCS),
           "-jobname", jobname, inp_path.name]
    print(f"  [Rank {rank}] Launching QSite ({_QSITE_PROCS} proc): {inp_path.name}", flush=True)

    # Live progress: QM/MM relaxed scans run for many minutes per frame, during
    # which a blocking call looks frozen. Launch non-blocking (still -WAIT, so the
    # driver process stays alive until the job finishes) and emit a heartbeat every
    # CFG.QSITE_PROGRESS_INTERVAL_SEC. Progress is read from the Jaguar log: count
    # completed scan points against CFG.QSITE_SCAN_NSTEPS for a k/N bar, falling
    # back to elapsed-time only when no markers are present yet.
    _interval = max(5, int(getattr(CFG, "QSITE_PROGRESS_INTERVAL_SEC", 20)))
    _total    = max(1, int(getattr(CFG, "QSITE_SCAN_NSTEPS", 1)))
    _log_path = qsite_dir / f"{jobname}.log"

    def _scan_done():
        try:
            txt = _log_path.read_text(errors="ignore")
        except OSError:
            return None
        for _m in ("Scan point", "scan point", "Geometry optimization has converged"):
            _n = txt.count(_m)
            if _n:
                return _n
        return None

    try:
        proc = _sp.Popen(cmd, cwd=str(qsite_dir),
                         stdout=_sp.DEVNULL, stderr=_sp.STDOUT)
    except Exception as e:
        console_info(f"    {ConsoleColours.FAIL}[Rank {rank}] QSite launch failed: {e}{ConsoleColours.ENDC}")
        return False

    _t0 = time.time()
    while proc.poll() is None:
        time.sleep(_interval)
        _elapsed = time.time() - _t0
        _done = _scan_done()
        if _done is not None:
            _pct  = min(100, int(100 * _done / _total))
            _fill = _pct // 5
            _bar  = "#" * _fill + "-" * (20 - _fill)
            # Linear ETA from the scan points already converged; shown only once at least one point is
            # done and the scan is not finished, so an idle start or a completed scan prints no estimate.
            _eta  = (f" | ETA ~{_elapsed * (_total - _done) / _done / 60:.1f} min"
                     if 0 < _done < _total else "")
            print(f"  [Rank {rank}] QSite {jobname}: [{_bar}] {_done}/{_total} pts "
                  f"({_pct}%) | {_elapsed / 60:.1f} min{_eta}", flush=True)
        else:
            print(f"  [Rank {rank}] QSite {jobname}: running… {_elapsed / 60:.1f} min elapsed",
                  flush=True)

    _rc = proc.returncode
    _elapsed = time.time() - _t0
    if _rc == 0:
        print(f"  [Rank {rank}] QSite {jobname} finished in {_elapsed / 60:.1f} min (rc=0).", flush=True)
    else:
        console_info(f"    {ConsoleColours.WARNING}[Rank {rank}] QSite {jobname} exited rc={_rc} "
                     f"after {_elapsed / 60:.1f} min — check {_log_path.name}.{ConsoleColours.ENDC}")
    # Self-check: surface the common QSite failure where the QM-region charge /
    # electron count is inconsistent (Jaguar 'incorrect molecular charge', odd
    # electrons) and it silently skips every scan point → an empty/NaN barrier.
    # Report it loudly so it is never mistaken for a converged result.
    try:
        _reason = _qsite_scan_failure_reason(qsite_dir, jobname)
        if _reason:
            console_info(f"    {ConsoleColours.FAIL}[Rank {rank}] QSite produced no valid "
                         f"reaction coordinate: {_reason}{ConsoleColours.ENDC}")
    except Exception:
        pass
    return True


def _qsite_scan_failure_reason(qsite_dir: Path, job_name: str) -> "str | None":
    """Read the QSite/Jaguar .out and return a human-readable reason when the scan
    failed to produce a valid PES (charge/electron mismatch, skipped points), else
    None. Used to convert Jaguar's silent 'Skipping to next scan point' into a clear
    diagnosis rather than an empty barrier."""
    try:
        _o = next((p for p in ([qsite_dir / f"{job_name}.out"] + sorted(qsite_dir.glob("*.out")))
                   if p.exists() and p.stat().st_size > 0), None)
        if _o is None:
            return None
        _t = _o.read_text(errors="ignore")

        # Impact-level abort: Impact dies during parameter assignment, BEFORE any SCF or
        # scan point, so there is no "Skipping" line to catch further down — the failure is
        # otherwise invisible and the QSite launcher still exits rc=0. The S-OPLS case is the
        # one this pipeline meets: OPLS3e/OPLS4 (S-OPLS) is rejected on a frozen-orbital-cut
        # QM region (QM bonded to MM, not H-capped), which is why CFG.QSITE_MM_FF selects
        # OPLS_2005. Surface the exact Impact error so a dead job is never read as a barrier.
        _die = re.search(r"%IMPACT-E \(die\):.*?\n(?:\s*%IMPACT-E:.*\n?)*", _t)
        if _die or "Impact exited with an error" in _t:
            _msg = re.sub(r"\s+", " ", _die.group(0)).strip() if _die else "Impact exited with an error"
            if re.search(r"S-OPLS forcefield is only available", _t):
                return ("Impact rejected the MM force field: S-OPLS (OPLS3e/OPLS4) is not allowed "
                        "on a frozen-orbital-cut QM region. Set CFG.QSITE_MM_FF empty (OPLS_2005). "
                        f"[{_msg}]")
            return f"Impact aborted before any SCF step — no barrier possible. [{_msg}]"

        _skips = _t.count("Skipping to next scan point")
        if re.search(r"incorrect molecular charge|Odd number of electrons", _t, re.I) or _skips:
            _m = re.search(r"Molecular charge:\s*(-?\d+)", _t)
            _got = f" (Jaguar reads charge {_m.group(1)})" if _m else ""
            return (f"QM-region charge/electron mismatch{_got} — {_skips} scan point(s) "
                    f"skipped. Reduce CFG.QSITE_MAX_QM_RESIDUES or check molchg/protonation.")

        """
        Charge-consistency check even when the SCF converges. QSite derives the QM
        charge from the QM atoms' partial charges at the QM/MM boundary and can
        override the `molchg` written in &gen (e.g. the intended −2 becomes a more
        negative value when each Cα–Cβ cut leaks ~−1 of backbone charge into the QM
        count). The SCF then runs on the wrong electron count and the barrier is
        for the wrong charge state — an error that never raises. Compare the
        `molchg` requested in the .in with the net charge Jaguar actually used and
        warn on any mismatch so it is caught during the run, not after.
        """
        _inp = next((p for p in ([qsite_dir / f"{job_name}_QSite_SN2.in"] + sorted(qsite_dir.glob("*_QSite_SN2.in")))
                     if p.exists()), None)
        _mreq = re.search(r"molchg\s*=\s*(-?\d+)", _inp.read_text(errors="ignore")) if _inp else None
        _mrun = re.search(r"net molecular charge:\s*(-?\d+)", _t)
        if _mreq and _mrun and int(_mreq.group(1)) != int(_mrun.group(1)):
            return (f"QM charge OVERRIDDEN — requested molchg={_mreq.group(1)} but Jaguar "
                    f"ran net charge {_mrun.group(1)} (QM/MM boundary re-derivation). The "
                    f"barrier is for the wrong charge state; verify QM-region protonation / "
                    f"cut boundaries before trusting ΔE‡.")
    except Exception:
        return None
    return None


# -----------------------------------------------------------------------------
# SECTION 7b: QM/MM SCAN PARSING → BARRIER, REACTION PROFILE, FLUORIDE CHARGE
# Turns the Jaguar/QSite relaxed-scan .out into ΔE‡ / ΔE_rxn, the PES, and the
# departing-F charge. Best-effort regex parsing (Jaguar output varies by version);
# every failure path returns NaN/empty and is non-fatal.
# -----------------------------------------------------------------------------
def _extract_scan_energies(text: str) -> "list[float]":
    """Per-scan-point energy series, returned in KCAL/MOL, from a Jaguar/QSite
    relaxed-scan .out.

    QSite prints each converged scan point's QM/MM total as
    ``Total Energy of the system...... -X.XXXXXE+03 kcal/mol`` — already kcal/mol,
    so this layout is read as-is. Hartree-based Jaguar layouts (scan-summary table,
    SCFE) are supported as fallbacks and converted with HARTREE_TO_KCAL. Returns []
    if fewer than two points parse.

    CAVEAT: a scan whose .out contains 'Skipping to next scan point' has NOT fully
    converged (points were dropped); the returned series is then short and any
    barrier from it is only a partial estimate — check QSite_NScan against
    CFG.QSITE_SCAN_NSTEPS before trusting ΔE‡."""
    _h2k = float(getattr(CFG, "HARTREE_TO_KCAL", 627.509474))
    # Restrict to the scan region: any 'Total Energy' printed BEFORE the first
    # 'Geometry scan coordinates' block is the pre-scan initial-structure reference,
    # not a scan point, and would corrupt the reactant baseline (a huge false barrier).
    _scan0 = re.search(r"Geometry scan coordinates", text, re.IGNORECASE)
    _body = text[_scan0.start():] if _scan0 else text
    # (1) PRIMARY — QSite QM/MM per-point total energy, ALREADY in kcal/mol.
    _qmmm = re.findall(
        r"Total Energy of the system\.*\s*(-?\d[\d.]*(?:[eE][+-]?\d+)?)\s*kcal", _body, re.IGNORECASE)
    if len(_qmmm) >= 2:
        return [float(x) for x in _qmmm]
    # (2) Fallback — Jaguar geometry-scan summary table (hartree → kcal).
    # Search _body (post geometry-scan header), not raw text, so a pre-scan baseline energy
    # cannot be captured as a scan point.
    _tbl = re.findall(r"^\s*\d+\s+[-\d.]+\s+(-\d+\.\d{4,})\s*$", _body, re.MULTILINE)
    if len(_tbl) >= 2:
        return [float(x) * _h2k for x in _tbl]
    # (3) Fallback — Jaguar SCFE converged energies (hartree → kcal); _body only, same reason.
    _scfe = re.findall(r"SCFE:.*?(-\d+\.\d{4,})", _body)
    if len(_scfe) >= 2:
        return [float(x) * _h2k for x in _scfe]
    return []


def _extract_scan_coordinates(text: str) -> "list[float]":
    """The constrained scan value actually used at each point, read from the Jaguar output.

    Jaguar echoes the active constraint ("  r = 3.5#") before EVERY geometry-optimisation step, not
    once per scan point, so consecutive duplicates are collapsed — leaving one value per point (the
    constraint is monotonic along the scan). Reading it is the only way to keep the PES aligned when
    nofail=1 drops a non-converged point.
    """
    _raw = [float(_m) for _m in re.findall(r"^\s*r\s*=\s*(-?\d+(?:\.\d+)?)#", text, re.M)]
    _out = []
    for _v in _raw:
        if not _out or _v != _out[-1]:
            _out.append(_v)
    return _out


def _extract_fluoride_charge_series(text: str) -> "list[float]":
    """The departing fluorine's Mulliken charge at each scan point.

    Jaguar prints the population analysis as a label row over a charge row:

        Atom       C1           F2           H3
        Charge    0.35999     -0.37483     -0.33090

    so the charge is taken from the atoms LABELLED F. Reading 'the most negative float in the block'
    instead would return a carboxylate oxygen (~-0.7, squarely inside any fluoride window) and report
    it as the fluoride.

    ONE fluorine is followed across the whole scan, identified by its Jaguar atom label. The departing
    F is resolved at the PRODUCT end, where it has become fluoride (→ ~-0.9) and the spectators remain
    near -0.3, and that same label is then read back at every scan point. Taking the most negative F
    independently at each point would let a spectator fluorine stand in for the leaving group at the
    reactant end — where all the fluorines are still near-degenerate — so the reported reactant and
    product charges could describe two different atoms.
    """
    series = []
    _blocks = []
    for _blk in re.split(r"(?i)Atomic charges from Mulliken population analysis", text)[1:]:
        """
        The block is read to its own end — the 'sum of atomic charges' line Jaguar prints after the
        table — not to a fixed character budget. A QM region of eight residues plus waters wraps the
        charge table over many label/charge row pairs, and a fixed cut can fall in the middle of it
        and silently drop the fluorine.
        """
        _end = re.search(r"(?i)sum of atomic charges", _blk)
        _lines = _blk[:_end.start() if _end else len(_blk)].splitlines()
        _labels, _charges = [], []
        for _k, _ln in enumerate(_lines):
            if _ln.strip().startswith("Atom") and _k + 1 < len(_lines):
                _chg = _lines[_k + 1]
                if _chg.strip().startswith("Charge"):
                    _labels += _ln.split()[1:]
                    _charges += _chg.split()[1:]
        _f = {}
        for _lab, _c in zip(_labels, _charges):
            if _lab[:1].upper() == "F" and _lab[1:].isdigit():
                try:
                    _f[_lab] = float(_c)
                except ValueError:
                    continue
        if _f:
            _blocks.append(_f)

    if not _blocks:
        return series
    # The leaving fluorine, fixed at the product end and then followed backwards through the scan.
    _leaving = min(_blocks[-1], key=_blocks[-1].get)
    for _blk_f in _blocks:
        series.append(_blk_f.get(_leaving, min(_blk_f.values())))
    return series

def parse_qsite_profile(qsite_dir: Path, job_name: str) -> dict:
    """Parse the QM/MM SN2 relaxed scan into a full reaction profile: the barrier
    and reaction energy (kcal/mol), the per-point PES, the reconstructed reaction
    coordinate (Nu_O···C distance, Å), and a best-effort departing-fluoride charge.
    The scan runs NAC (reactant, r≈3.5 Å) → product (r≈1.3 Å); ΔE‡ = E_max−E_react
    proves C–F cleavage is surmountable, ΔE_rxn = E_product−E_react whether it is
    downhill, and the fluoride charge → ~−0.9 whether F actually leaves as F⁻."""
    _empty = {"QSite_Barrier_kcal": np.nan, "QSite_dErxn_kcal": np.nan, "QSite_NScan": 0,
              "coord": [], "energy_kcal": [], "f_charge": [],
              "F_Charge_Reactant": np.nan, "F_Charge_Product": np.nan, "F_Charge_Delta": np.nan}
    try:
        _outs = ([qsite_dir / f"{job_name}_QSite_SN2.out"] + sorted(qsite_dir.glob("*.out")))
        text = ""
        for _o in _outs:
            if _o.exists() and _o.stat().st_size > 0:
                text = _o.read_text(errors="ignore"); break
        if not text:
            return _empty
        e = _extract_scan_energies(text)   # already in kcal/mol
        if len(e) < 3:
            return _empty
        """
        The reactant is the LOWEST point before the barrier top, not simply the first scan point.
        The scan starts from a constrained geometry that relaxes as the optimisation proceeds, so
        e[0] is generally not the reactant minimum, and measuring from it reports the activation
        energy of whichever geometry the scan happened to start from.

        The transition state must lie AFTER the reactant minimum, and it must be a real maximum.

        When the global maximum IS the first scan point, the profile is monotonically downhill from
        an unrelaxed start: there is no barrier anywhere on the sampled coordinate. The honest report
        is then NO BARRIER — not 0.0 kcal/mol. A reported 0.0 would read downstream as 'this reaction
        is barrierless', the strongest possible claim, when what actually happened is that the scan
        never resolved a transition state. `Is_Defluorinating` keys on the barrier, so a fabricated
        0.0 would promote precisely the jobs whose scans failed.

        The reactant minimum is likewise searched only BEFORE the maximum. For a strongly exothermic
        profile the global minimum is the product, and measuring the barrier down from the product
        would report the reverse barrier.
        """
        _imax = e.index(max(e))
        if _imax == 0:
            console_info("    [!] QSite scan is monotonically downhill from the first point — no "
                         "transition state on the sampled coordinate. Barrier reported as NaN, "
                         "not 0.0: the scan did not resolve a TS.")
            return _empty
        """
        A transition state is a MAXIMUM the reaction passes THROUGH: the energy must come back down on
        the far side. If the highest point is the last point sampled, the scan is still climbing when
        it ends — a steric wall, a scan window that stops short of the saddle, or a coordinate that
        simply does not cross one. Reporting e[-1] − e[reactant] as a barrier there invents a ΔE‡ for a
        reaction the scan never showed happening, and an invented barrier is worse than no barrier: it
        would be carried forward as a turnover number.
        """
        if _imax >= len(e) - 1:
            console_info("    [!] QSite scan is still climbing at the last point — the maximum is an "
                         "endpoint, not a saddle the reaction passes through. No TS resolved (a steric "
                         "wall or a scan window that stops short of it). Barrier reported as NaN.")
            return _empty
        _react = min(e[:_imax])          # reactant well: strictly before the TS, never the product
        _energy_kcal = [round(x - _react, 3) for x in e]   # relative to the reactant minimum
        """
        The reaction coordinate is READ from the output, never rebuilt by index: nofail=1 means a
        non-converged point is skipped, so the i-th energy is not necessarily the i-th grid value
        and an index-built coordinate silently shifts the whole PES.
        """
        _coord = _extract_scan_coordinates(text)
        if len(_coord) != len(e):
            _start = float(CFG.QSITE_SCAN_START)
            _step = float(CFG.QSITE_SCAN_STEP)
            _coord = [round(_start + _step * i, 3) for i in range(len(e))]
            console_info(f"    [!] Scan coordinate not parseable from the QSite output "
                         f"({len(e)} energies) — falling back to the CFG grid; a skipped scan point "
                         f"would misalign the PES.")
        _fq = _extract_fluoride_charge_series(text)
        """
        The reactant fluoride charge is read at the SAME scan point the barrier is measured from
        (the pre-TS minimum), not at the first point of the scan. Quoting the charge of a geometry
        that is not the reactant, alongside a barrier that is measured from the reactant, describes
        two different states as one.
        """
        _i_react = e.index(_react) if _react in e else 0
        _fq_react = _fq[_i_react] if len(_fq) > _i_react else (_fq[0] if _fq else np.nan)
        _fq_prod = _fq[-1] if _fq else np.nan
        return {
            # The barrier is the post-reactant maximum minus the reactant minimum — e[_imax], not
            # max(e), so a downhill-from-the-start profile cannot report a 0.0 kcal/mol barrier.
            "QSite_Barrier_kcal": round(e[_imax] - _react, 2),
            "QSite_dErxn_kcal":   round(e[-1] - _react, 2),
            "QSite_NScan":        len(e),
            "coord":              _coord,
            "energy_kcal":        _energy_kcal,
            "f_charge":           _fq,
            "F_Charge_Reactant":  round(_fq_react, 3) if _fq_react == _fq_react else np.nan,
            "F_Charge_Product":   round(_fq_prod, 3) if _fq_prod == _fq_prod else np.nan,
            "F_Charge_Delta":     (round(_fq_prod - _fq_react, 3)
                                   if (_fq_react == _fq_react and _fq_prod == _fq_prod) else np.nan),
        }
    except Exception:
        return _empty


def parse_qsite_barrier(qsite_dir: Path, job_name: str) -> dict:
    """Scalar barrier/reaction-energy/fluoride subset of parse_qsite_profile, for
    the multi-frame ensemble aggregation. Returns NaNs when nothing parseable."""
    _p = parse_qsite_profile(qsite_dir, job_name)
    return {k: _p[k] for k in ("QSite_Barrier_kcal", "QSite_dErxn_kcal", "QSite_NScan",
                               "F_Charge_Reactant", "F_Charge_Product", "F_Charge_Delta")}


def plot_qsite_reaction_profile(out_path: Path, job_name: str, rank, profile: dict) -> None:
    """The direct 'did it defluorinate' figure: the QM/MM potential-energy surface
    along the SN2 reaction coordinate (Nu_O···C compression), with the activation
    barrier ΔE‡ and reaction energy ΔE_rxn marked, and — where parseable — the
    departing-fluoride Mulliken charge dropping toward −1 (F becoming free F⁻).

    Called from process_single_job, which runs inside the ThreadPoolExecutor, so the
    figure is built under PLOT_LOCK. Pyplot's figure registry is global state: two ranks
    plotting concurrently can interleave into one another's figure, and a savefig on the
    CURRENT figure can then write a different candidate's PES under this candidate's name.
    Nothing raises — the file is written and looks plausible — so the lock is the only
    thing standing between the QM/MM phase and a mislabelled reaction profile. The save
    goes through the figure object rather than the pyplot global for the same reason.
    """
    try:
        x = profile.get("coord") or []
        y = profile.get("energy_kcal") or []
        if len(x) < 3 or len(y) < 3:
            return
        n = min(len(x), len(y))
        x, y = x[:n], y[:n]
        _C = CFG.DEFLUOR_FIG_COLOUR
        with PLOT_LOCK:
            fig, ax = plt.subplots(figsize=(8.5, 5.5))
            try:
                ax.plot(x, y, "-o", color=_C["pes"], lw=2, ms=4, zorder=3, label="QM/MM PES")
                _imax = int(np.argmax(y))
                ax.scatter([x[_imax]], [y[_imax]], s=120, color=_C["ts"], zorder=5, label="transition state")
                ax.scatter([x[-1]], [y[-1]], s=90, color=_C["product"], zorder=5, label="product")
                ax.annotate(f"ΔE‡ = {profile.get('QSite_Barrier_kcal', float('nan')):.1f} kcal/mol",
                            (x[_imax], y[_imax]), xytext=(6, 8), textcoords="offset points",
                            fontsize=9, color=_C["ts"], fontweight="bold")
                ax.annotate(f"ΔE_rxn = {profile.get('QSite_dErxn_kcal', float('nan')):.1f}",
                            (x[-1], y[-1]), xytext=(6, -12), textcoords="offset points",
                            fontsize=9, color=_C["product"], fontweight="bold")
                ax.set_xlabel("Reaction coordinate — Nu(O)···C distance (Å), reactant → product",
                              fontweight="bold")
                ax.set_ylabel("Relative QM/MM energy (kcal/mol)")
                ax.invert_xaxis()   # NAC (large r) on the left → product (small r) on the right
                clean_spines(ax)
                _fq = profile.get("f_charge") or []
                if len(_fq) >= 3:
                    _fc = _C["f_charge"]
                    ax2 = ax.twinx()
                    _xf = [x[min(int(i * (n - 1) / (len(_fq) - 1)), n - 1)] for i in range(len(_fq))]
                    ax2.plot(_xf, _fq, "--s", color=_fc, lw=1.4, ms=3, alpha=0.85, label="departing-F charge")
                    ax2.set_ylabel("Mulliken charge on departing F (→ −1 = fluoride)", color=_fc)
                    ax2.tick_params(axis="y", labelcolor=_fc)
                # C–F cleavage verdict (make the "did the bond break" answer unmistakable).
                _fq_final = (_fq[-1] if len(_fq) >= 1 else None)
                _scf = profile.get("scissile_f_label", "scissile C–F")
                if _fq_final is not None and _fq_final <= CFG.QSITE_F_CHARGE_CLEAVED:
                    ax.text(0.5, 0.02, f"C–F CLEAVED — {_scf} F → {_fq_final:+.2f} e (free fluoride)",
                            transform=ax.transAxes, ha="center", va="bottom", fontsize=10, fontweight="bold",
                            color=_C["product"],
                            bbox=dict(boxstyle="round,pad=0.3", fc=_C["cleaved_bg"], ec=_C["product"], alpha=0.95))
                elif _fq_final is not None:
                    ax.text(0.5, 0.02, f"C–F intact — {_scf} F charge {_fq_final:+.2f} e (bond not broken)",
                            transform=ax.transAxes, ha="center", va="bottom", fontsize=10, fontweight="bold",
                            color=_C["ts"],
                            bbox=dict(boxstyle="round,pad=0.3", fc=_C["intact_bg"], ec=_C["ts"], alpha=0.95))
                ax.legend(loc="upper left")
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", UserWarning)
                    fig.savefig(out_path, dpi=int(CFG.VIS_FIGURE_DPI), bbox_inches="tight")
            finally:
                plt.close(fig)
        console_info(f"    QSite reaction profile saved : {out_path.name}")
    except Exception as _e:
        console_info(f"    [!] QSite reaction-profile plot failed ({_e}).")


# =============================================================================
# SECTION 8: CORE ANALYSIS ENGINE
# =============================================================================

# Number of threads used for parallel DTR frame pre-loading.
# Each thread opens its own read_traj instance (thread-safe).
# For HDD: keep at 2-4 (seek contention); for SSD: use cpu_count().
_N_PRELOAD_WORKERS = min(4, CFG.GLOBAL_MAX_WORKERS)

# Solvent sphere radius (Å) around the nucleophile for water blockade.
# Solvent that enters this radius in ANY sampled frame is tracked — reduces
# per-frame work from ~10,000 solvent atoms to ~100-200 while covering the SN2
# runway across the whole trajectory (not just frame 0).
_SOL_SPHERE_RADIUS = getattr(CFG, "SOLVENT_SPHERE_RADIUS", 20.0)

# Number of evenly-spaced frames sampled to build the blockade solvent superset.
# The frame-0-only sphere missed waters that diffuse into the runway later,
# biasing the blockade metric low on long trajectories; the union over these
# samples removes that bias while keeping a fixed pre-load set.
_SOL_SPHERE_SAMPLE_FRAMES = max(1, int(CFG.SOLVENT_SPHERE_SAMPLE_FRAMES))

# QSite execution settings. Initialised from CFG and overridden per-run by main()
# from the CLI. Kept as module globals (CFG is a frozen dataclass and cannot be
# mutated). Read-only inside the worker threads.
_QSITE_RUN = CFG.QSITE_RUN
_QSITE_PROCS = CFG.QSITE_PROCS

# QSite/Jaguar QM/MM engine (Impact `main1h`) is effectively SINGLE-THREADED for these small
# frozen-cut QM regions — measured at ~100% of ONE core regardless of `-PARALLEL N`. Serialising
# QSite therefore wasted cores-1 of the budget. Instead we run many QM/MM scans CONCURRENTLY, each
# on its own core, bounded by this semaphore (set in main() to min(cores-2, RAM cap)). Per-job
# procs is forced to 1 (the -PARALLEL flag does nothing here but spawn idle helpers). The prep that
# mutates the shared cms_model stays sequential; only the independent Jaguar subprocesses run in
# parallel. Total concurrent scans is capped by both this semaphore and the number of QM/MM jobs
# actually available (ranks x QSITE_N_FRAMES).
_QSITE_SEM = threading.Semaphore(1)   # reassigned in main() to the concurrency cap


def _qsite_concurrency() -> int:
    """Concurrent single-threaded QSite scans to run: min(cores-2, RAM budget). Each Jaguar QM
    job needs ~1.5 GB; keep 60% of MemAvailable for QSite so the MD/analysis keep headroom."""
    cores = max(1, (os.cpu_count() or 4) - 2)
    ram_cap = cores
    try:
        with open("/proc/meminfo") as _f:
            _avail_kb = next(int(_l.split()[1]) for _l in _f if _l.startswith("MemAvailable"))
        ram_cap = max(1, int((_avail_kb / 1024.0 / 1024.0) * 0.60 / 1.5))
    except Exception:
        pass
    return max(1, min(cores, ram_cap))
# Serialises multi-line progress prints from parallel rank threads so lines never interleave.
_PROGRESS_LOCK = threading.Lock()
# True only when a single rank is processed (e.g. the --ranks 1 test): the per-frame counter
# then updates in place with a carriage return; with multiple parallel ranks it appends
# throttled lines instead (in-place \r from concurrent threads would garble).
_PROGRESS_CR = False

# QM-region coordinating waters: solvent O within this radius (Å) of the
# scissile carbon / leaving fluorine / nucleophile oxygen enters the QM region
# as a whole molecule (F⁻ leaving-group stabilisation). Capped to keep the QM
# electron count tractable.
_QM_WATER_RADIUS = float(CFG.QSITE_QM_WATER_RADIUS)
_QM_WATER_MAX    = int(CFG.QSITE_QM_WATER_MAX)


def _eaf_at(series: np.ndarray, frame_t: float, t_start: float, eaf_dt: float) -> float:
    """Return EAF scalar value at the MD frame time, interpolated by index."""
    if len(series) == 0:
        return np.nan
    idx = max(0, min(int((frame_t - t_start) / eaf_dt), len(series) - 1))
    return float(series[idx])


# -----------------------------------------------------------------------------
# SECTION 8a: NAC PERSISTENCE (continuous strict-NAC dwell → nanoseconds)
# The real "time in position": longest/mean uninterrupted strict-NAC run, not the
# frame-count fraction a flickering ligand can inflate.
# -----------------------------------------------------------------------------
def _nac_dwell_stats(flags: "list[int]", ns_per_frame: float) -> dict:
    """Continuous-residence statistics for a per-frame strict-NAC boolean series.

    Distinguishes genuine catalytic pre-organisation (the ligand SITTING in the
    reactive geometry) from mere frame-count fraction (which a ligand flickering
    in and out can inflate). Returns the longest and mean continuous run and the
    total, each converted to nanoseconds via ns_per_frame (= total simulated ns /
    analysed frames), plus the number of distinct NAC episodes.
    """
    runs, cur = [], 0
    for f in flags:
        if f:
            cur += 1
        elif cur:
            runs.append(cur); cur = 0
    if cur:
        runs.append(cur)
    _max = max(runs) if runs else 0
    _mean = float(np.mean(runs)) if runs else 0.0
    _tot = int(sum(runs))
    return {
        "NAC_Dwell_Max_ns":   round(_max  * ns_per_frame, 3),
        "NAC_Dwell_Mean_ns":  round(_mean * ns_per_frame, 3),
        "NAC_Total_ns":       round(_tot  * ns_per_frame, 3),
        "NAC_Episodes":       len(runs),
    }


# -----------------------------------------------------------------------------
# SECTION 8b: PER-JOB TRAJECTORY ANALYSIS ENGINE
# NAC geometry + Dream-Team tracking + WaterMap/EAF/MM-GBSA integration + QM/MM
# frame selection, producing the per-rank NAC_Data.csv, dashboards, and stats row.
# -----------------------------------------------------------------------------
def process_single_job(rank: int, work_dir: Path, df_ranked: pd.DataFrame,
                       master_out_dir: Path, lig_resname: str, stride: int,
                       triad_override: dict = None,
                       fallback_nuc:  int = DREAM_TEAM_REF.get('Nuc', 110),
                       fallback_base: int = DREAM_TEAM_REF.get('Base', 280),
                       fallback_acid: int = DREAM_TEAM_REF.get('Acid', 134)):
    """Orchestrates the full analysis pipeline for one MD trajectory."""

    # ── Locate the MD job folder for this R-number first ───────────────────────
    # Flexible discovery: try every directory pattern Desmond/pipeline may produce.
    _md_root   = work_dir / "05_MD_Simulations"
    _candidates = [
        _md_root / f"desmond_md_job_R_{rank}",
        _md_root / f"desmond_md_job_Rank_{rank}",
        _md_root / f"Results_MD_Simulation_Rank_{rank}",
        _md_root / f"desmond_md_Rank_{rank}",
        _md_root / f"md_job_Rank_{rank}",
    ]
    # Also glob for any directory carrying the rank token: an exact '_R_{rank}'
    # suffix (avoids matching R_1 against R_10/R_11) or a '*Rank_{rank}*' name.
    if _md_root.exists():
        _candidates += sorted(_md_root.glob(f"*_R_{rank}"))
        _candidates += sorted(_md_root.glob(f"*Rank_{rank}*"))
    job_folder = next((p for p in _candidates if p.is_dir()), None)
    if job_folder is None:
        console_info(f"    {ConsoleColours.FAIL}[!] Job folder missing for Rank {rank} in {_md_root}{ConsoleColours.ENDC}")
        return None

    # ── Select the ranked row by the case this MD job ACTUALLY contains ─────────
    # The MD folders are whichever candidates had MD-ready trajectories, so the
    # folder's R-number need not equal its Scientific_Rank. Match the ranked row
    # on the model identity embedded in the -out.cms; fall back to
    # Scientific_Rank == rank only when the identity cannot be read, and warn on
    # any disagreement so a mis-mapped case is never analysed silently.
    _cms_id = identity_from_cms(job_folder, df_ranked)
    row = None
    if _cms_id is not None:
        _match = df_ranked[df_ranked['job_name'] == _cms_id]
        if not _match.empty:
            row = _match.iloc[0]
            _row_rank = row.get('Scientific_Rank')
            if _row_rank is not None and int(_row_rank) != int(rank):
                console_info(f"    {ConsoleColours.WARNING}[!] Folder R_{rank} contains "
                             f"'{_cms_id}' (Scientific_Rank {_row_rank}) — mapping by cms "
                             f"identity, not the folder number.{ConsoleColours.ENDC}")
    if row is None:
        try:
            row = df_ranked[df_ranked['Scientific_Rank'] == rank].iloc[0]
            if _cms_id is not None:
                console_info(f"    {ConsoleColours.WARNING}[!] cms identity '{_cms_id}' not "
                             f"found in ranked sheet — fell back to Scientific_Rank == {rank}."
                             f"{ConsoleColours.ENDC}")
        except (IndexError, KeyError):
            console_info(f"    {ConsoleColours.WARNING}[!] No ranked entry for Rank {rank}. Skipping.{ConsoleColours.ENDC}")
            return None

    job_name   = row['job_name']
    print(f"  [Rank {rank}] Initialising analysis for: {ConsoleColours.OKBLUE}{job_name}{ConsoleColours.ENDC}", flush=True)

    _dir_label  = re.sub(r'_\d{1,3}_[A-Za-z][A-Za-z0-9]+$', '', job_name)
    job_out_dir = master_out_dir / f"Rank_{rank}_{_dir_label}"
    job_out_dir.mkdir(parents=True, exist_ok=True)
    console_info(f"Processing Rank {rank} [Stride={stride}]: {ConsoleColours.OKBLUE}{job_name}{ConsoleColours.ENDC}")

    # ── Load trajectory ────────────────────────────────────────────────────────
    # First try flat lookup (single-stage jobs: *-out.cms and *_trj at folder root).
    """
    The flat globs are ordered by SEGMENT NUMBER because cms_path is taken as _cms_flat[-1] — the final
    stage of the run. Filesystem order is arbitrary on Linux, so an unsorted flat glob would hand '[-1]'
    whichever segment the directory happened to list last and call it the final one.

    No tgz-extract / recursive fallback: multisim leaves the PRODUCTION stage flat at the job-dir root
    ({job}_trj + {job}-out.cms) and archives ONLY the earlier equilibration stages as *_N-out.tgz.
    Extracting those and concatenating them as production (as an earlier version did) analyses a
    Brownie/NVT/NPT relax as if it were production, with no root .ene so every relax frame counts as
    post-equilibration. Step 06 now refuses that on the build side; 07 mirrors it — a missing flat
    {job}_trj means the MD did not finish, so skip the rank.
    """
    def _seg_key(p: Path, suffix: str) -> int:
        _m = re.search(r'_(\d+)' + re.escape(suffix), p.name)
        return int(_m.group(1)) if _m else 0

    _cms_flat = sorted(job_folder.glob("*-out.cms"), key=lambda p: _seg_key(p, '-out.cms'))
    _trj_flat = sorted(job_folder.glob("*_trj"), key=lambda p: _seg_key(p, '_trj'))

    if not _cms_flat or not _trj_flat:
        console_info(f"    {ConsoleColours.FAIL}[!] No flat {job_folder.name}/*_trj — the production MD "
                     f"did not finish (only equilibration stages are archived); skipping rank.{ConsoleColours.ENDC}")
        return None

    cms_path = _cms_flat[-1]   # final stage CMS (highest segment number)
    msys_model, cms_model = topo.read_cms(str(cms_path))

    tr = LazyTrajectory(_trj_flat)
    print(f"  [Rank {rank}] Loading trajectory: {len(_trj_flat)} segment(s) | {len(tr):,} total frames...", flush=True)

    # ── WaterMap spatial sites (maegz) — flexible naming discovery ───────────
    _wm_root = work_dir / "03_WaterMaps"
    # Find any subdirectory carrying the rank token: 'watermap_R_N'
    # (exact '_R_{rank}' suffix so R_1 ≠ R_10) or a '*Rank_N*' name.
    _wm_dir_candidates = (
        sorted(_wm_root.glob(f"*_R_{rank}")) + sorted(_wm_root.glob(f"*[Rr]ank*{rank}*"))
        if _wm_root.exists() else []
    )
    wm_maegz = None
    for _wmd in _wm_dir_candidates:
        if not _wmd.is_dir():
            continue
        # Accept *_wm.maegz or *wm*.maegz anywhere in the dir tree (depth ≤ 2)
        for _pat in ("*_wm.maegz", "*wm*.maegz", "*.maegz"):
            _hits = sorted(_wmd.glob(_pat)) + sorted(_wmd.glob(f"**/{_pat}"))
            if _hits:
                wm_maegz = _hits[0]; break
        if wm_maegz:
            break
    wm_sites = load_watermap_sites(wm_maegz)

    # ── WaterMap CSV statistics — flexible naming discovery ───────────────────
    wm_csv_path = None
    if _wm_root.exists():
        _csv_hits = (
            sorted(_wm_root.glob(f"*_R_{rank}.csv")) +
            sorted(_wm_root.glob(f"*_R_{rank}_*.csv")) +
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
    _has_time = hasattr(tr[0], 'time') and hasattr(tr[-1], 'time')
    if not _has_time:
        # Without frame timestamps the fallback below is a frame INDEX, not picoseconds — the
        # equilibration cut and the NAC dwell would then be in frames, not ns. Desmond trajectories
        # always carry .time, so this is a loud warning rather than a silent unit switch.
        console_info("    [!] Trajectory frames carry no .time — equilibration/dwell fall back to FRAME "
                     "indices (not ps); check the trajectory if these times look wrong.")
    t_start  = tr[0].time  if _has_time else 0.0
    t_end    = tr[-1].time if _has_time else float(n_tr)
    sim_span = max(t_end - t_start, 1.0)
    eaf_dt   = sim_span / len(eaf_msa) if len(eaf_msa) > 0 else 1.0

    # ── Resolve residue number hints: alignment map > triad_override > CLI ─────
    aln_dict     = parse_mapping(row.get('Full_Sequence_Alignment_Map', ''))
    to           = triad_override or {}
    nuc_hint  = aln_dict.get(DREAM_TEAM_REF.get('Nuc', 110))  or to.get('nuc')  or fallback_nuc
    base_hint = aln_dict.get(DREAM_TEAM_REF.get('Base', 280)) or to.get('base') or fallback_base
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

    # sorted(): these become ARRAY INDEX orders downstream, and a set has none.
    warhead_c = sorted(set(c for c, f in cf_pairs))
    lig_f     = sorted(set(f for c, f in cf_pairs))
    nuc_num   = cms_model.atom[idx_nuc[0]].resnum
    # Base and Acid residue numbers from Smart-Lock — needed for QM region expansion.
    base_num  = cms_model.atom[idx_base[0]].resnum if idx_base else None
    acid_num  = cms_model.atom[idx_acid[0]].resnum if idx_acid else None
    # Fluoride-cradle residue numbers (TRP/TYR/HIS donors) → QM region cuts.
    cradle_nums = sorted({cms_model.atom[i].resnum for i in (idx_cradle or [])})

    # ── Mapping sanity check: Nuc–Base distance in frame 0 ─────────────────────
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

    """
    QSite QM-region residues are mapped from the ranked sheet's explicit Mapped_*
    columns for THIS job, every time. Each homolog's catalytic positions differ
    substantially from the FAcD control (the sheet stores that alignment), so the
    QM cuts must follow the per-job mapping rather than reference numbering or the
    3D geometric detector — the latter can latch onto nearby TRP/TYR/HIS that are
    not the true cradle (e.g. 38/44/47/68), dropping the real Trp/Tyr cradle and
    leaving the departing F under-stabilised. Resolution order per role:
    Mapped_* column → parsed Full_Sequence_Alignment_Map (dream_mapped) →
    Smart-Lock geometry. cradle = Mapped_Stabiliser_W + Mapped_Stabiliser_Y.
    """
    _qm_nuc  = mapped_resnum(row, 'Mapped_Nucleophile')  or dream_mapped.get('Nuc')    or nuc_num
    _qm_base = mapped_resnum(row, 'Mapped_Base')         or dream_mapped.get('Base')   or base_num
    _qm_acid = mapped_resnum(row, 'Mapped_Acid')         or dream_mapped.get('Acid')   or acid_num
    _qm_stab = (mapped_resnum(row, 'Mapped_Stabiliser_H') or dream_mapped.get('Stab_H')
                or stab_f_num
                or (cms_model.atom[dt_indices['Stab_H'][0]].resnum if dt_indices.get('Stab_H') else None))
    _map_cradle = [mapped_resnum(row, 'Mapped_Stabiliser_W'),
                   mapped_resnum(row, 'Mapped_Stabiliser_Y')]
    _qm_cradle  = sorted({r for r in _map_cradle if r})
    if not _qm_cradle:    # fall back to alignment map, then the geometric scan
        _qm_cradle = sorted({dream_mapped[r] for r in ('Stab_W', 'Stab_Y') if dream_mapped.get(r)}) \
                     or cradle_nums
    # The per-job mapped residues (_qm_*) are the authoritative QM-region set passed
    # to QSite below; the geometry-derived nuc/base/acid (nuc_num/base_num/acid_num)
    # stay in force for the per-frame distance analysis and are left untouched.
    console_info(f"    [QSite] QM region mapped from ranked sheet — "
                 f"Nuc {_qm_nuc}, Base {_qm_base}, Acid {_qm_acid}, StabH {_qm_stab}, "
                 f"cradle {_qm_cradle}")

    # Filter out ptypes with special characters (e.g. '/') that break ASL parsing.
    _safe_restypes = [r for r in _SOLVENT_RESTYPES if r.isalnum() or '_' in r]
    _sol_asl       = " OR ".join(f"res.ptype {r}" for r in _safe_restypes)
    sol_indices    = list(cms_model.select_atom(f"({_sol_asl}) AND a.el O"))

    """
    Counter-ion capping of the reactive centre.

    A PFAS carboxylate pairs strongly with Na⁺, and the FAcD active site is an anion trap: the Asp
    nucleophile, the Asp acid and the substrate carboxylate all carry negative charge. System Builder's
    ion-exclusion region keeps counter-ions out of the site at BUILD time, but nothing stops one
    diffusing in during the simulation and coordinating the nucleophile's Oδ or the ligand head, where
    it screens the very charge that drives the SN2 attack and can sit for tens of nanoseconds.

    Such frames are geometrically indistinguishable from productive NAC frames — the distance and the
    angle can both look ideal — so they must be FLAGGED, not silently averaged into the NAC statistics.
    A frame is capped when a cation sits within CFG.CATION_CAP_DIST of either the nucleophile Oδ or a
    ligand carboxylate oxygen (inner-sphere coordination; Na⁺–O carboxylate contact is ~2.4 Å).
    """
    cation_indices = list(cms_model.select_atom("a.el Na OR a.el K OR a.el Mg OR a.el Ca"))
    try:
        lig_head_o = list(cms_model.select_atom(f'(res.ptype "{lig_resname}") AND a.el O'))
    except Exception:
        lig_head_o = []

    # ===============================================================================
    # Trajectory pre-loading: batch-read all strided frames into RAM
    # Replaces per-frame frame.pos() calls with fast numpy array indexing.
    # ===============================================================================
    _t_pre = time.time()

    # 1. Filter solvent to a sphere around the active site (95%+ reduction).
    #    UNION across evenly-spaced sampled frames (centred on the nucleophile in
    #    each) rather than freezing the sphere at frame 0, so any water that ever
    #    obstructs the SN2 runway is tracked — the frozen sphere biased the
    #    blockade metric low on long trajectories.
    if idx_nuc and sol_indices:
        _n_fr       = len(tr)
        _sample_idx = (np.unique(np.linspace(0, _n_fr - 1,
                                             min(_SOL_SPHERE_SAMPLE_FRAMES, _n_fr)).astype(int))
                       if _n_fr > 1 else np.array([0]))
        _sol_arr = np.asarray(sol_indices)
        _keep    = np.zeros(len(_sol_arr), dtype=bool)
        for _fi in _sample_idx:
            _fr   = tr[int(_fi)]
            _cen  = np.mean([_fr.pos(i) for i in idx_nuc], axis=0)
            _spos = np.array([_fr.pos(s) for s in sol_indices])
            _keep |= (np.linalg.norm(_spos - _cen, axis=1) <= _SOL_SPHERE_RADIUS)
        _sol_use = _sol_arr[_keep].tolist()
    else:
        _sol_use = []

    # 2. Walden substituents for all warhead carbons (needed for dihedral).
    # For SN2 TS flattening, the 3 non-leaving substituents on the sp3 carbon
    # become coplanar. We must include all bonded atoms (including H and F).
    _walden_subs = sorted(set(
        b.atom2.index for c_idx in warhead_c
        for b in cms_model.atom[c_idx].bond
    ))

    # 3. Build the master atom list to pre-load.
    _preload_set = (
        list(idx_nuc) + list(warhead_c) + list(lig_f) +
        list(idx_base or []) + list(idx_acid or []) +
        list(idx_cradle or []) + _walden_subs +
        [a for idxs in dt_indices.values() for a in idxs] +
        _sol_use + list(cation_indices) + list(lig_head_o)
    )
    """
    The Cα atoms are preloaded ALONGSIDE the reaction-centre atoms, because the WaterMap sites are
    static coordinates in the WaterMap input's frame while the MD protein diffuses and tumbles
    (measured here: the Cα centroid moves 30-45 Å over 1 µs). Every frame is therefore superimposed
    onto that reference before a site distance is taken; the sites are moved WITH the protein, so
    the distances stay in the frame's own coordinate system and the minimum-image convention
    remains valid.
    """
    _wm_ref_ca = load_watermap_reference_ca(wm_maegz) if wm_sites else {}
    _ca_atoms  = ([a.index for a in cms_model.atom
                   if a.pdbname.strip() == "CA" and int(a.resnum) in _wm_ref_ca]
                  if _wm_ref_ca else [])
    _preload_set = list(_preload_set) + _ca_atoms
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
    _li_ca     = [_a2l[a] for a in _ca_atoms]
    _li_cat    = [_a2l[a] for a in cation_indices]
    _li_ligo   = [_a2l[a] for a in lig_head_o]
    _ref_ca_xyz = (np.array([_wm_ref_ca[int(cms_model.atom[a].resnum)] for a in _ca_atoms])
                   if _ca_atoms else np.empty((0, 3)))
    _wm_pos_ref = (np.array([site['pos'] for site in wm_sites])
                   if (wm_sites and len(_ca_atoms) >= 3) else np.empty((0, 3)))
    if wm_sites and len(_ca_atoms) < 3:
        console_info(f"    {ConsoleColours.WARNING}[!] WaterMap sites cannot be aligned to the "
                     f"trajectory (no matching Cα reference) — the protein tumbles, so site "
                     f"distances would be meaningless. WaterMap scoring is disabled for this rank."
                     f"{ConsoleColours.ENDC}")
        wm_sites = []
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

    print(f"  [Rank {rank}] Starting chunked streaming: {_nf} frames | setup completed in {time.time()-_t_pre:.1f}s", flush=True)
    console_info(
        f"  Chunked streaming: {_nf} frames × {_n_pre} atoms "
        f"({len(_sol_use)} solvent) | chunk={_CHUNK_SIZE} | {_n_chunks} chunks | "
        f"setup {time.time()-_t_pre:.1f}s"
    )

    """
    NPT equilibration is verified BEFORE any frame is treated as a sample. Frames recorded before
    the box settled are still written to the per-frame CSV (flagged Post_Equilibration = 0) but are
    excluded from the NAC counters, the dwell statistics and the QM/MM frame pool: an equilibrium
    average taken over a relaxing system is not an equilibrium average, and a still-shrinking box
    produces perfectly well-formed NAC numbers that mean nothing.
    """
    _equil = check_md_equilibration(job_folder, job_folder.name)
    _equil_t_ps = float(_equil.get("MD_Equil_Time_ps") or 0.0)
    if not np.isfinite(_equil_t_ps):
        _equil_t_ps = 0.0
    n_pre_equil = 0

    _t_loop = time.time()

    def _iter_frames():
        """Yield (fi, f_idx, pos, box, frame_t) one frame at a time.

        Loads _CHUNK_SIZE frames using a fresh local LazyTrajectory wrapper,
        yields them all, and then deletes the wrapper to free cached coordinate memory.
        """
        for _ci in range(_n_chunks):
            _cs = _ci * _CHUNK_SIZE
            _ce = min(_cs + _CHUNK_SIZE, _nf)
            _cz = _ce - _cs
            _local_tr = LazyTrajectory(_trj_flat)
            _cpos = np.empty((_cz, _n_pre, 3), dtype=np.float64)
            _cbox = [None] * _cz
            _ctim = [0.0]  * _cz
            for _li in range(_cz):
                _fx    = _stride_list[_cs + _li]
                _frame = _local_tr[_fx]
                _cpos[_li] = _frame.pos(_preload_list)
                _cbox[_li] = _frame.box  if hasattr(_frame, 'box')  else None
                _ctim[_li] = (_frame.time if hasattr(_frame, 'time')
                               else t_start + _fx * sim_span / n_tr)
            for _li in range(_cz):
                yield _cs + _li, _stride_list[_cs + _li], _cpos[_li], _cbox[_li], _ctim[_li]
            del _cpos, _cbox, _ctim, _local_tr
            _elapsed_c = time.time() - _t_loop
            _done_c    = _ce
            _rate_c    = _done_c / _elapsed_c if _elapsed_c > 1e-6 else 0.0
            _eta_c     = int((_nf - _done_c) / _rate_c) if _rate_c > 0.0 else 0
            print(f"  [Rank {rank}] Processing: {_done_c}/{_nf} frames ({_rate_c:.0f} fr/s, ETA ~{_eta_c}s)", flush=True)
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
    _qm_candidates: list = []   # (score, f_idx, nuc_o_idx, lig_c_idx) for every productive NAC frame
    n_pocket = n_relaxed = n_strict = n_triad = n_nac = 0
    n_fold_reject = 0        # frames whose fold no longer superimposes on the WaterMap reference
    n_cation_capped = 0      # frames with a counter-ion coordinating the nucleophile or the ligand head
    _wm_radius = CFG.WATERMAP_SITE_RADIUS

    # ===============================================================================
    # Per-frame analysis loop (chunked streaming — _CHUNK_SIZE frames at a time)
    # ===============================================================================
    # Live progress cadence: ~100 updates across the trajectory (never every frame).
    _prog_step = max(1, _nf // 100)
    for _fi, f_idx, _p, box, frame_t in _iter_frames():
        # Live per-rank frame counter — "Processing: Rank_N: k/total frames (pct%)".
        # In place (\r) for a single rank; throttled appended lines when ranks run parallel.
        if (_fi % _prog_step == 0) or (_fi + 1 == _nf):
            _pct_fr = int(100 * (_fi + 1) / _nf) if _nf else 100
            _msg_fr = f"  Processing: Rank_{rank}: {_fi + 1:,}/{_nf:,} frames ({_pct_fr}%)"
            with _PROGRESS_LOCK:
                if _PROGRESS_CR and sys.stdout.isatty():
                    print(f"\r{_msg_fr}   ", end="", flush=True)
                    if _fi + 1 == _nf:
                        print(flush=True)
                else:
                    print(_msg_fr, flush=True)

        """
        Frames recorded before the box equilibrated are not samples of the equilibrium ensemble, so
        they are excluded from every counter below and from the QM/MM frame pool. They are still
        written to the per-frame CSV with Post_Equilibration = 0, so the exclusion is visible rather
        than silent.
        """
        post_equil = bool(float(frame_t) >= _equil_t_ps)
        if not post_equil:
            n_pre_equil += 1

        # ── Position arrays from cache (numpy indexing — no API calls) ──────────
        _pos_nuc_f = _p[_li_nuc]  if _li_nuc   else np.empty((0, 3))
        _pos_c_f   = _p[_li_c]    if _li_c     else np.empty((0, 3))
        _pos_f_f   = _p[_li_f]    if _li_f     else np.empty((0, 3))

        """
        ── Which oxygen attacks which carbon ─────────────────────────────────────────────────────
        FAcD's nucleophile is the aspartate carboxylate: it attacks the α-carbon and displaces the
        fluoride. Its two oxygens are resonance-equivalent, so the one that REACTS is the one lined
        up for backside attack, not the one that happens to be nearer — they sit ~2.2 Å apart and can
        point in different directions, so choosing by distance can report a side-on approach for a
        frame whose other oxygen is properly anti-periplanar.

        The pair is therefore chosen by the backside O–C–F angle it produces (distance breaks ties).

        The NAC gates are then evaluated on THAT pair: `nac_nuc_dist` is the approach distance of the
        oxygen that supplied the angle. Gating on the minimum over both oxygens while measuring the
        angle on the other one describes a nucleophile that does not exist — one oxygen lending its
        reach, its partner lending its trajectory — and lets a frame pass the NAC criterion on an
        oxygen that is not attacking. `min_nuc_dist` is retained as the true closest approach over
        both oxygens, which is what pocket residency means, and is reported alongside.
        """
        if _pos_nuc_f.size and _pos_c_f.size:
            _nc_d    = _mic_dists_2d(_pos_nuc_f, _pos_c_f, box)
            min_nuc_dist = float(_nc_d.min())          # closest approach, over BOTH oxygens
            nac_nuc_dist = min_nuc_dist                # approach of the ATTACKING oxygen (set below)
            _best_key = None
            best_nuc_idx = best_ca_idx = None
            for _ni, _n_atom in enumerate(idx_nuc):
                for _ci, _c_atom in enumerate(warhead_c):
                    _bf = _li_cf.get(_a2l[_c_atom], [])
                    if not _bf:
                        continue
                    _cpos = _p[_a2l[_c_atom]]
                    _v_cn = get_mic_vector(_p[_a2l[_n_atom]], _cpos, box)
                    _n_cn = np.linalg.norm(_v_cn)
                    if _n_cn <= 1e-6:
                        continue
                    _v_cf = np.array([get_mic_vector(_p[_fl], _cpos, box) for _fl in _bf])
                    _l_cf = np.linalg.norm(_v_cf, axis=1)
                    _ok = _l_cf > 1e-6
                    if not np.any(_ok):
                        continue
                    _a = np.degrees(np.arccos(np.clip(
                        np.dot(_v_cf[_ok], _v_cn) / (_l_cf[_ok] * _n_cn), -1.0, 1.0)))
                    """
                    The attacking oxygen satisfies BOTH near-attack conditions at once — short
                    approach and linear trajectory — because the SN2 needs them on the same atom.
                    Scored jointly with the same graded NAC terms Step 02 uses, so the two stages
                    agree on which oxygen is attacking. Ranking on angle alone selects a
                    well-aligned oxygen that can sit far out of reach, and the frame is then
                    credited with a trajectory no nucleophile could travel.
                    """
                    _ang_j = float(_a.max())
                    _d_j = float(_nc_d[_ni, _ci])
                    _nac_j = (_sigmoid07(_d_j, CFG.SOFT_K_NUC, CFG.NAC_DIST_STRICT)
                              * _sigmoid07(_ang_j, CFG.SOFT_K_ANG, CFG.NAC_ANGLE_STRICT))
                    _key = (_nac_j, _ang_j, -_d_j)
                    if _best_key is None or _key > _best_key:
                        _best_key, best_nuc_idx, best_ca_idx = _key, _n_atom, _c_atom
                        nac_nuc_dist = _d_j
            if best_nuc_idx is None:                   # no bonded F resolvable — fall back
                _flat = int(np.argmin(_nc_d))
                _nl, _cl = divmod(_flat, len(warhead_c))
                best_nuc_idx, best_ca_idx = idx_nuc[_nl], warhead_c[_cl]
                nac_nuc_dist = min_nuc_dist
        else:
            min_nuc_dist = nac_nuc_dist = float('inf')
            best_nuc_idx = best_ca_idx = None

        if post_equil and min_nuc_dist <= POCKET_RESIDENCY_DIST:
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
        if post_equil and triad_ok and min_nuc_dist <= POCKET_RESIDENCY_DIST:
            n_triad += 1

        # ── Counter-ion capping of the reactive centre ────────────────────────
        # A Na⁺ coordinating the nucleophile Oδ or the ligand carboxylate screens the charge the
        # SN2 depends on, while leaving the NAC distance and angle looking perfectly productive.
        cat_nuc_d = cat_lig_d = np.inf
        if _li_cat:
            _pos_cat = _p[_li_cat]
            if _pos_nuc_f.size:
                cat_nuc_d = float(np.min(_mic_dists_2d(_pos_cat, _pos_nuc_f, box)))
            if _li_ligo:
                cat_lig_d = float(np.min(_mic_dists_2d(_pos_cat, _p[_li_ligo], box)))
        cation_capped = bool(min(cat_nuc_d, cat_lig_d) <= CFG.CATION_CAP_DIST)
        if post_equil and cation_capped:
            n_cation_capped += 1

        if post_equil and nac_nuc_dist <= THRESHOLD_RELAXED_NAC_DIST and max_ang >= THRESHOLD_RELAXED_NAC_ANGLE:
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

        """
        The WaterMap sites are carried into THIS frame before any distance is measured: the frame's
        Cα atoms are superimposed on the reference the sites were computed in, and the same rigid
        transform is applied to the site coordinates. Without it the comparison is between two
        unrelated coordinate frames — the protein has diffused tens of ångströms — and every site
        match is noise.
        """
        _wm_frame = wm_sites
        _fold_rmsd = np.nan
        if len(_wm_pos_ref):
            _R, _t, _fold_rmsd = kabsch_transform(_ref_ca_xyz, unwrap_ca_trace(_p[_li_ca], box))
            if _fold_rmsd > float(CFG.MD_FOLD_RMSD_MAX):
                """
                The fold in this frame no longer superimposes on the WaterMap reference. Carrying the
                hydration sites through a transform fitted to a mismatched fold places them at
                arbitrary positions, and every site-based blockade term computed from them is noise.
                The frame keeps its geometry-only metrics; its WaterMap contribution is withheld.
                """
                _wm_frame = []
                n_fold_reject += 1
            else:
                _wm_moved = _wm_pos_ref @ _R.T + _t
                _wm_frame = [{'pos': _wm_moved[_k], 'dG': wm_sites[_k]['dG'], 'num': wm_sites[_k]['num']}
                             for _k in range(len(wm_sites))]

        # ── Water blockade (vectorized over pre-loaded solvent positions) ──────
        _pos_sol_f = _p[_li_sol] if _li_sol else np.empty((0, 3))
        blockades  = (
            _blockade_vec(_p[_a2l[best_nuc_idx]], _p[_a2l[best_ca_idx]],
                          _pos_sol_f, box, _wm_frame,
                          CFG.WATERMAP_BLOCKADE_RADIUS, CFG.WATERMAP_MATCH_RADIUS)
            if best_nuc_idx else 0.0
        )

        # ── Walden TS flattening (improper dihedral, inline with cached pos) ──
        walden_flat = False
        walden = 1.0
        if best_f_idx and best_ca_idx:
            # The 3 equatorial substituents are all atoms bonded to Cα, except the leaving F
            _wsub_li = [li for li in _li_wsubs.get(_a2l[best_ca_idx], []) if li != _a2l[best_f_idx]]
            if len(_wsub_li) >= 3 and _calc_improper_dihedral is not None:
                try:
                    angle = _calc_improper_dihedral(
                        _p[_a2l[best_ca_idx]], _p[_wsub_li[0]],
                        _p[_wsub_li[1]], _p[_wsub_li[2]], box)
                    if abs(angle) < _WALDEN_IMPROPER_MAX:
                        walden = 1.1; walden_flat = True
                except Exception:
                    pass

        # ── WaterMap per-frame thermodynamic contribution ─────────────────────
        wm_score = 0.0
        if best_ca_idx:
            _pos_c_wm = _p[_a2l[best_ca_idx]]
            for site in _wm_frame:
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
        _productive = (np.clip(THRESHOLD_RELAXED_NAC_DIST - nac_nuc_dist, 0, None) * _SCORE_W_DIST
                       + np.clip(max_ang - THRESHOLD_RELAXED_NAC_ANGLE, 0, None) * _SCORE_W_ANGLE)
        _penalty    = blockades * _SCORE_W_BLOCK + max(0.0, -wm_score) * _SCORE_W_WM
        score       = _productive * walden - _penalty

        if post_equil and triad_ok and nac_nuc_dist <= THRESHOLD_RELAXED_NAC_DIST and max_ang >= THRESHOLD_RELAXED_NAC_ANGLE:
            n_relaxed += 1
        if post_equil and triad_ok and nac_nuc_dist <= THRESHOLD_STRICT_NAC_DIST and max_ang >= THRESHOLD_STRICT_NAC_ANGLE:
            n_strict += 1
        if cation_capped or not post_equil:
            """
            Disqualified as a QM/MM starting structure.

            A cation-capped frame would hand the QM region a Na⁺ screening the very nucleophile it is
            meant to model, and the computed barrier would describe that ion pair rather than the
            enzyme. A pre-equilibration frame is a snapshot of a box that is still relaxing, not of
            the equilibrium ensemble. Both are still reported, flagged, and (for capping) still
            contribute geometry to the NAC statistics; neither can be the structure a barrier is
            measured on.
            """
            _productive = 0.0
        if _productive > 0.0:
            # Candidate pool for multi-frame QM/MM: every pre-organised NAC frame,
            # ranked by score. The top CFG.QSITE_N_FRAMES give an ensemble barrier
            # (min/mean/σ) rather than a single best-frame lower bound.
            _qm_candidates.append((score, f_idx, best_nuc_idx, best_ca_idx))
        if score > best_score and _productive > 0.0:
            # Single most pre-organised (highest-score) NAC frame — kept as the
            # primary QM/MM frame (its barrier is the lower bound of the ensemble).
            best_score = score; ideal_frame_idx = f_idx

        results.append({
            # Core NAC geometry
            "Frame":               f_idx,
            # NAC_Distance_A is the approach of the ATTACKING oxygen (the one carrying NAC_Angle_Deg);
            # Nuc_Min_Dist_A is the closest approach over both Oδ, which is what pocket residency means.
            "NAC_Distance_A":      _r3(nac_nuc_dist),
            "Nuc_Min_Dist_A":      _r3(min_nuc_dist),
            "Tail_Cradle_Dist_A":  _r3(dist_tail),
            "NAC_Angle_Deg":       round(max_ang, 2),
            "Consensus":           round(score, 2),
            "Triad_NB":            _r3(nb_dist),
            "Triad_BA":            _r3(ba_dist),
            "Walden_TS_Flat":      int(walden_flat),
            "Fold_RMSD_A":         _r3(_fold_rmsd),
            # Counter-ion capping: a cation inside CFG.CATION_CAP_DIST of the nucleophile Oδ or the
            # ligand carboxylate. Such a frame's NAC geometry is real but its electrostatics are not.
            "Cation_Nuc_Dist_A":   _r3(cat_nuc_d if np.isfinite(cat_nuc_d) else np.nan),
            "Cation_LigO_Dist_A":  _r3(cat_lig_d if np.isfinite(cat_lig_d) else np.nan),
            "Cation_Capped":       int(cation_capped),
            # 0 = recorded before the box equilibrated: excluded from every counter and from the
            # QM/MM frame pool, but kept in the CSV so the exclusion is auditable.
            "Post_Equilibration":  int(post_equil),
            "Pocket_Bound":        int(min_nuc_dist <= POCKET_RESIDENCY_DIST),
            "NAC_Geom_Pass":       int(nac_nuc_dist <= THRESHOLD_RELAXED_NAC_DIST
                                       and max_ang >= THRESHOLD_RELAXED_NAC_ANGLE),
            "NAC_Strict_Pass":     int(nac_nuc_dist <= THRESHOLD_STRICT_NAC_DIST
                                       and max_ang >= THRESHOLD_STRICT_NAC_ANGLE),
            "Triad_Intact":        int(triad_ok),
            "Catalytic_Pass":      int(triad_ok
                                       and nac_nuc_dist <= THRESHOLD_RELAXED_NAC_DIST
                                       and max_ang >= THRESHOLD_RELAXED_NAC_ANGLE),
            # Dream Team distances to warhead carbon
            "DT_Nuc_LigC_A":    _r3(dt_dists.get('Nuc',    np.nan)),
            "DT_Base_LigC_A":   _r3(dt_dists.get('Base',   np.nan)),
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
    """
    Every average below is taken over the SAMPLED frames — those recorded after the box equilibrated.
    The pre-equilibration frames stay in the per-frame CSV (flagged) but must not enter an equilibrium
    average, and they must not sit in a percentage denominator either: counting them would deflate
    every occupancy by the length of the relaxation phase.
    """
    sampled = [r for r in results if r.get('Post_Equilibration', 1) == 1]
    total = len(sampled)
    total_frames_read = len(results)

    """
    THE BOX CAN SETTLE WHILE THE PROTEIN DOES NOT.

    check_md_equilibration() reads the .ene stream, so it sees the box volume and the temperature and
    nothing else. A variant that is unfolding — or whose active site is being prised open to swallow a
    bulky PFAS tail — reaches a stable volume and a stable temperature and passes a barostat-only gate.
    That is precisely the system the gate exists to reject, and it is the one it cannot see.

    The backbone Ca RMSD is already measured per frame (Kabsch, against the starting structure), so the
    structural half of the verdict is taken here, where those frames exist: over the sampled window the
    mean Ca RMSD must be within CFG.MD_EQUIL_CA_RMSD_MAX_A. Equilibration is thermodynamic AND
    structural, and a run that fails either one is not equilibrated.
    """
    _ca_vals = [r["Fold_RMSD_A"] for r in sampled
                if r.get("Fold_RMSD_A") is not None and r["Fold_RMSD_A"] == r["Fold_RMSD_A"]]
    _ca_mean = float(np.mean(_ca_vals)) if _ca_vals else float("nan")
    _ca_ok = bool(np.isfinite(_ca_mean) and _ca_mean <= float(CFG.MD_EQUIL_CA_RMSD_MAX_A))

    """
    AN UNVERIFIABLE GATE IS NOT A PASSED GATE.

    A test written `if MD_Equilibrated and _ca_vals and not _ca_ok` fails open: when _ca_vals is EMPTY —
    no WaterMap reference, so no Kabsch superposition and every Fold_RMSD_A is NaN — the `and _ca_vals`
    short-circuits and the structural half of the verdict is skipped in silence. A trajectory whose fold
    is never checked would then be reported as equilibrated, which is the one outcome this check exists to
    prevent.

    Missing evidence is now its own state: MD_Equilibrated becomes None (UNKNOWN), never True. None is
    already this function's contract for 'could not determine' — it is what it returns on a parse failure
    — so every consumer that guards on `is True` treats it as not-equilibrated, and nothing silently
    inherits a pass it never earned.
    """
    if _equil.get("MD_Equilibrated") is True:
        if not _ca_vals:
            console_info(f"    [!] {job_name}: box and thermostat settled, but the backbone could NOT be "
                         f"checked (no Ca RMSD — WaterMap reference absent). Equilibration is UNKNOWN, "
                         f"not confirmed.")
            _equil["MD_Equilibrated"] = None
        elif not _ca_ok:
            console_info(f"    [!] {job_name}: box and thermostat settled, but the BACKBONE did not — "
                         f"mean Ca RMSD {_ca_mean:.2f} A > {float(CFG.MD_EQUIL_CA_RMSD_MAX_A):.2f} A. "
                         f"The fold is still moving; not equilibrated.")
            _equil["MD_Equilibrated"] = False
    _equil["MD_Equil_CA_RMSD_A"] = round(_ca_mean, 3) if np.isfinite(_ca_mean) else np.nan

    # Viability denominators use n_pocket (bound frames), not total frames.
    # Catalytic_Viability_Pct = fraction of BOUND frames that are catalytically
    # competent — the meaningful metric for a pre-reactive ensemble.
    bound_rows = [r for r in sampled if r.get('Pocket_Bound', 0) == 1]

    def _mean_col(col: str, rows=None) -> float:
        src  = rows if rows is not None else sampled
        vals = [r[col] for r in src if col in r and r[col] == r[col]]  # NaN check
        return round(float(np.mean(vals)), 3) if vals else np.nan

    def _min_col(col: str, rows=None) -> float:
        src  = rows if rows is not None else sampled
        vals = [r[col] for r in src if col in r and r[col] == r[col]]  # NaN check
        return round(float(np.min(vals)), 3) if vals else np.nan

    stats = {
        "Job_Name":                 job_name,
        "Scientific_Rank":          rank,
        # Total_Frames is the SAMPLED count (post-equilibration) — the denominator of every
        # percentage below. Frames_Read is what the trajectory actually contained.
        "Total_Frames":             total,
        "Frames_Read":              total_frames_read,
        "Frames_Pre_Equilibration": n_pre_equil,
        **_equil,
        "Frames_In_Pocket":         n_pocket,
        "Frames_NAC_Geom_Only":     n_nac,
        "Frames_Triad_Intact":      n_triad,
        "Frames_Relaxed_Catalysis": n_relaxed,
        "Frames_Strict_Catalysis":  n_strict,
        # Frames whose fold no longer superimposes on the WaterMap reference (Cα RMSD >
        # CFG.MD_FOLD_RMSD_MAX): their hydration term is withheld, their geometry is kept.
        "Frames_Fold_Rejected":     n_fold_reject,
        # Frames with a counter-ion coordinating the nucleophile Oδ or the ligand carboxylate. Their
        # NAC geometry can look ideal while the catalytic charge is screened, so they are barred from
        # QM/MM frame selection and reported here. A high percentage means the ion-exclusion region
        # used at system-build time was too small, or the salt concentration is parking Na⁺ in the site.
        "Frames_Cation_Capped":     n_cation_capped,
        "Cation_Capped_Pct":        round((n_cation_capped / total) * 100, 2) if total else 0.0,
        "Pocket_Retention_Pct":     round((n_pocket  / total)    * 100, 2) if total    else 0.0,
        "Catalytic_Viability_Pct":  round((n_relaxed / n_pocket) * 100, 2) if n_pocket else 0.0,
        "Strict_Viability_Pct":     round((n_strict  / n_pocket) * 100, 2) if n_pocket else 0.0,
        "Triad_Integrity_Pct":      round((n_triad   / n_pocket) * 100, 2) if n_pocket else 0.0,
        "NAC_Geom_Only_Pct":        round((n_nac     / n_pocket) * 100, 2) if n_pocket else 0.0,
        # Continuous strict-NAC residence in ns (real "time in position", not frame
        # fraction). ns_per_frame = total simulated ns / ALL analysed frames (uniform stride):
        # sim_span covers the whole trajectory, so its denominator must be every analysed frame
        # (total_frames_read), not the post-equilibration subset — otherwise the full span is
        # divided by fewer frames and every dwell is scaled up.
        # ns_per_frame is NaN when the trajectory carries no .time: sim_span is then a FRAME COUNT, not
        # ps, so a dwell built from it would invent ~1 ps/frame (10 ns from 10,000 frames). NaN makes
        # _nac_dwell_stats yield NaN dwells and _verdict withholds — as the stride>1 leg already does —
        # rather than reporting a fabricated residence time.
        **_nac_dwell_stats([r.get("NAC_Strict_Pass", 0) for r in sampled],
                           (((sim_span / 1000.0) / total_frames_read)
                            if (_has_time and total_frames_read) else float("nan"))),
        # The stride the dwell was measured at. A run of consecutive ANALYSED frames is only
        # evidence of continuous residence when every frame was analysed: at stride > 1 the
        # ligand may leave and re-enter the reactive geometry between two samples and the whole
        # interval still counts as one unbroken dwell, biasing it upward in units of stride × dt.
        # The verdict gates on this, so the stride travels with it.
        "MD_Analysis_Stride":       int(stride),
        "Sim_Total_ns":             round(sim_span / 1000.0, 2),
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
    for k in [CFG.COL_TIER, CFG.COL_NUC_DIST, CFG.COL_SN2, CFG.COL_CONF]:
        if k in row: stats[k] = row[k]


    # ── NAC-conditioned MM-GBSA: total ΔG, energy-component decomposition, and
    #    catalytic-machinery engagement, all restricted to the strict-NAC frames ──
    # Shows WHERE the binding energy comes from in the reactive pose (Coulomb =
    # electrostatic pre-organisation for the SN2, Covalent = ligand strain toward
    # the TS, etc.) and whether the nucleophile + fluoride cradle are geometrically
    # engaged. Best-effort: the per-frame Prime MM-GBSA CSV (Step 06) is row-aligned
    # to trajectory frames; component columns are matched by name.
    _nac_fr = [r["Frame"] for r in results if r.get("NAC_Strict_Pass", 0)]
    _decomp = {}   # component → (nac_mean, global_mean) for the decomposition plot
    try:
        _mmg = (sorted(job_folder.glob("*prime*mmgbsa*.csv"))
                or sorted(job_folder.glob("*mmgbsa*.csv")))
        if _mmg:
            _mdf = pd.read_csv(_mmg[0])
            """
            Map trajectory frame → MM-GBSA row by the frame index Step 06 stamps on the CSV
            ('Frame'). Never by row POSITION: Prime scores every CFG.MMGBSA_STEP_SIZE-th frame,
            so row i is frame i·step, and positional lookup would silently attribute one frame's
            energy to another. Only a CSV with no Frame column (an old every-frame run) falls
            back to position, where row i genuinely is frame i.
            """
            _row_of_frame = None
            if "Frame" in _mdf.columns:
                _fr_idx = pd.to_numeric(_mdf["Frame"], errors="coerce")
                _row_of_frame = {int(f): i for i, f in enumerate(_fr_idx) if f == f}
            elif len(_mdf) < len(results):
                """
                No Frame column AND fewer rows than trajectory frames ⇒ the CSV was strided but
                never stamped. Row i is NOT frame i, so positional lookup would attribute the
                wrong frame's energy. Refuse it: a missing number beats a confidently wrong one.
                Re-run Step 06 (it stamps the frame index) to recover this metric.
                """
                console_info(f"    [!] MM-GBSA CSV has {len(_mdf)} rows for {len(results)} frames "
                             f"and no 'Frame' column — it was strided but not frame-stamped. "
                             f"Skipping NAC-conditioned MM-GBSA rather than mis-aligning frames; "
                             f"re-run Step 06 to stamp it.")
                raise ValueError("MM-GBSA CSV lacks frame indices")
            _nac_scored = ([f for f in _nac_fr if f in _row_of_frame] if _row_of_frame is not None
                           else [f for f in _nac_fr if 0 <= f < len(_mdf)])
            if _nac_fr and not _nac_scored:
                console_info(f"    [!] None of the {len(_nac_fr)} strict-NAC frames were scored by "
                             f"MM-GBSA — NAC-conditioned ΔG unavailable (lower CFG.MMGBSA_STEP_SIZE).")
            elif _row_of_frame is not None and len(_nac_scored) < len(_nac_fr):
                console_info(f"    [i] NAC-conditioned MM-GBSA uses {len(_nac_scored)} of "
                             f"{len(_nac_fr)} strict-NAC frames (the rest fall between the "
                             f"MM-GBSA stride's sampled frames).")

            def _nac_vs_global(col):
                """Mean over the scored strict-NAC frames, the global mean, and the NAC sample's
                size and dispersion. n and SD are returned rather than discarded: they are what
                decides whether the difference between the two means can be claimed at all."""
                _s = pd.to_numeric(_mdf[col], errors="coerce")
                _g = float(_s.mean())
                _v = [float(_s.iloc[_row_of_frame[f] if _row_of_frame is not None else f])
                      for f in _nac_scored]
                _v = [x for x in _v if x == x]
                _sd = float(np.std(_v, ddof=1)) if len(_v) >= 2 else np.nan
                return ((float(np.mean(_v)) if _v else np.nan),
                        (_g if _g == _g else np.nan), len(_v), _sd)

            _dgc = getattr(CFG, "MMGBSA_DG_COLUMN", "r_psp_MMGBSA_dG_Bind")
            if _dgc not in _mdf.columns:
                _c = [c for c in _mdf.columns if re.search(r"dg.?bind", c, re.I)]
                _dgc = _c[0] if _c else None
            if _dgc is not None:
                _nac, _glob, _n_nac, _sd_nac = _nac_vs_global(_dgc)
                stats["MMGBSA_dG_Global_Mean_kcal"] = round(_glob, 2) if _glob == _glob else np.nan
                stats["MMGBSA_dG_NAC_Mean_kcal"]    = round(_nac, 2) if _nac == _nac else np.nan
                stats["MMGBSA_dG_NAC_SD_kcal"]      = round(_sd_nac, 2) if _sd_nac == _sd_nac else np.nan
                stats["MMGBSA_NAC_Frames_Scored"]   = int(_n_nac)
                '''
                The penalty is the difference of two means, and it is only reportable if the NAC mean
                is. Prime's per-frame ΔG_bind scatter runs to several kcal/mol, so below
                CFG.MMGBSA_NAC_MIN_FRAMES the difference is dominated by the sampling noise of the
                smaller sample. The mean and SD above are kept — they are the evidence — but the
                penalty itself is withheld rather than printed to two decimals from one frame.
                '''
                if _n_nac >= int(CFG.MMGBSA_NAC_MIN_FRAMES) and _nac == _nac and _glob == _glob:
                    stats["MMGBSA_NAC_Penalty_kcal"] = round(_nac - _glob, 2)
                else:
                    stats["MMGBSA_NAC_Penalty_kcal"] = np.nan
                    if 0 < _n_nac < int(CFG.MMGBSA_NAC_MIN_FRAMES):
                        console_info(f"    [!] NAC-conditioned MM-GBSA has only {_n_nac} scored "
                                     f"strict-NAC frame(s) (< CFG.MMGBSA_NAC_MIN_FRAMES = "
                                     f"{int(CFG.MMGBSA_NAC_MIN_FRAMES)}) — mean reported, "
                                     f"penalty withheld.")
            # Energy-component decomposition (the "which forces" breakdown).
            _components = {"Coulomb": r"coulomb", "vdW": r"vdw|van.?der.?waals",
                           "Covalent": r"covalent", "H-bond": r"h.?bond",
                           "Lipo": r"lipo", "Packing": r"packing",
                           "SolvGB": r"solv.?gb|gb\b|solvation", "SelfCont": r"self.?cont"}
            for _name, _pat in _components.items():
                # Require a dG_Bind DELTA column (not the absolute Complex/Receptor/
                # Ligand energies, which are thousands of kcal/mol and would swamp the
                # decomposition). e.g. r_psp_MMGBSA_dG_Bind_Coulomb, not *_Complex_Coulomb.
                _hit = [c for c in _mdf.columns
                        if re.search(r"dg.?bind", c, re.I) and re.search(_pat, c, re.I)
                        and not re.search(r"complex|receptor|ligand", c, re.I)]
                if _hit:
                    _nac, _glob = _nac_vs_global(_hit[0])
                    _decomp[_name] = (_nac, _glob)
                    stats[f"MMGBSA_{_name}_NAC_Mean_kcal"] = round(_nac, 2) if _nac == _nac else np.nan
    except Exception as _e:
        console_info(f"    [!] NAC-conditioned MM-GBSA decomposition skipped ({_e}).")

    # Catalytic-machinery engagement: mean distance of the nucleophile + fluoride
    # cradle + clamps to the warhead carbon over the strict-NAC frames (did the
    # machinery actually close in during the reactive windows?).
    _dt_nac = {}
    _dt_cols = {"Nuc (Asp)": "DT_Nuc_LigC_A", "Base": "DT_Base_LigC_A",
                "Acid": "DT_Acid_LigC_A", "Clamp1": "DT_Clamp1_LigC_A",
                "Clamp2": "DT_Clamp2_LigC_A", "Cradle-His": "DT_StabH_LigC_A",
                "Cradle-Trp": "DT_StabW_LigC_A", "Cradle-Tyr": "DT_StabY_LigC_A"}
    _nac_set = set(_nac_fr)
    for _lbl, _col in _dt_cols.items():
        _vv = [r[_col] for r in results
               if r.get("Frame") in _nac_set and _col in r and r[_col] == r[_col]]
        if _vv:
            _dt_nac[_lbl] = float(np.mean(_vv))
            stats[f"{_col.replace('_LigC_A','')}_NAC_Mean_A"] = round(float(np.mean(_vv)), 2)

    # ── Output: per-frame CSV with rolling EAF smoothing ──────────────────────
    df_res = pd.DataFrame(results)
    if 'EAF_MSA' in df_res.columns and df_res['EAF_MSA'].notna().any():
        df_res['EAF_MSA_Smooth'] = df_res['EAF_MSA'].rolling(window=50, min_periods=1).mean()
    _utils_mod.atomic_write_csv(df_res, job_out_dir / f"{job_name}_NAC_Data.csv")

    """
    The per-job statistics, written beside the per-frame table. Everything the dashboard prints and
    the master ranking carries — pocket retention, viability, triad integrity, NAC dwell, WaterMap
    terms, the Dream-Team distances — is computed once, in the frame loop, and until now existed
    only in memory: the master CSV is written at the END of the whole run, so a figure could not be
    redrawn, nor a number checked, without re-reading a 100,000-frame trajectory. NumPy scalars are
    cast so the file is plain JSON.
    """
    def _plain(v):
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (np.floating,)):
            return None if np.isnan(v) else float(v)
        if isinstance(v, float) and np.isnan(v):
            return None
        return v

    _stats_path = job_out_dir / f"{job_name}_MD_Stats.json"
    write_json_atomic(_stats_path, {k: _plain(v) for k, v in stats.items()})
    console_info(f"    Per-job statistics saved : {_stats_path.name}")

    print(f"  [Rank {rank}] Generating dashboard chart...", flush=True)
    generate_individual_dashboard(
        df_res, job_name, job_out_dir / f"{job_name}_NAC_Dashboard.png", stats)

    if ideal_frame_idx != -1:
        # Top-N pre-organised NAC frames for QM/MM: the ensemble barrier (min/mean/σ)
        # over several reactive frames is defensible where a single best-frame value
        # is only a lower bound. Frames are the highest-scoring distinct NAC frames.
        _n_qm = max(1, int(getattr(CFG, "QSITE_N_FRAMES", 1)))
        _seen, _sel = set(), []
        for _cand in sorted(_qm_candidates, key=lambda t: t[0], reverse=True):
            if _cand[1] in _seen:
                continue
            _seen.add(_cand[1]); _sel.append(_cand)
            if len(_sel) >= _n_qm:
                break

        def _prep_qsite_frame(_fi: int, _geom: dict, _folder: Path, _primary: bool):
            """Snap → PBC-repair → droplet → .in for ONE frame. Returns the .in path if the frame
            still needs a Jaguar run, or None if the folder is already present (cached → parse-only).

            This step MUTATES the shared cms_model (update/make_whole/center), so it must run
            sequentially across frames — only the QM/MM run itself is parallelised (below). Repairs
            periodic wrapping and re-centres the box on the LIGAND so a surface active site sits at
            the box middle and the droplet trim keeps its first-shell waters (vacuum artefact →
            SCF divergence otherwise).
            """
            topo.update_cms(cms_model, tr[_fi])
            topo.make_whole_cms(msys_model, cms_model)
            _gids = topo.asl2gids(cms_model, f"res.ptype {lig_resname}")
            topo.center_cms(msys_model, _gids, cms_model)
            if _primary:
                cms_model.fsys_ct.write(str(job_out_dir / f"{job_name}_Ideal_Final.maegz"))
            if _folder.exists():
                console_info(f"    [Rank {rank}] QSite folder exists — will parse: {_folder.name}")
                return None
            _folder.mkdir(parents=True, exist_ok=True)
            # QSite/Jaguar reads uncompressed .mae; trim the periodic box to a
            # local solvation droplet so the QM/MM MM region stays tractable.
            _mae = _folder / f"{job_name}_Ideal_Final.mae"
            _ds = write_qsite_droplet(cms_model, _mae, lig_resname)
            print(f"  [Rank {rank}] QSite .mae solvent (frame {_fi}): {_ds}", flush=True)
            _inp = generate_qsite_inputs(
                _mae, job_name, _qm_nuc, _qm_stab,
                _geom['lig_c'], _geom['nuc_o'],
                base_num=_qm_base, acid_num=_qm_acid, lig_resname=lig_resname,
                cradle_nums=_qm_cradle)
            for _a in (_mae, _inp):
                if not _a.exists() or _a.stat().st_size == 0:
                    console_info(f"    [!] QSite input missing/empty for {job_name}: {_a.name}")
            if _primary:
                console_qmm_ready(f"Best frame: {_fi} | Score: {best_score:.2f} | "
                                  f"QSite input: {_inp.name}")
            return _inp

        print(f"  [Rank {rank}] QM/MM: {len(_sel)} frame(s) (best {ideal_frame_idx}, "
              f"score {best_score:.2f}) — ensemble SN2 barrier.", flush=True)
        _folds = [(job_out_dir / f"{job_name}_QSite_SN2") if _k == 0
                  else (job_out_dir / f"{job_name}_QSite_SN2_f{_cand[1]}")
                  for _k, _cand in enumerate(_sel)]
        # PHASE 1 — sequential prep (mutates cms_model); collect the frames that still need a run.
        _runjobs = []
        for _k, _cand in enumerate(_sel):
            _inp = _prep_qsite_frame(_cand[1], {'nuc_o': _cand[2], 'lig_c': _cand[3]}, _folds[_k], _k == 0)
            if _inp is not None and _QSITE_RUN:
                _runjobs.append((_folds[_k], _inp))
        # PHASE 2 — run the QM/MM scans CONCURRENTLY (each Jaguar QM engine is single-threaded),
        # globally bounded by _QSITE_SEM so parallel ranks × frames never exceed the core/RAM cap.
        if _QSITE_RUN and _runjobs:
            def _run_one(_job):
                _fold_r, _inp_r = _job
                with _QSITE_SEM:
                    console_info(f"    [Rank {rank}] QSite scan (concurrent) — {_fold_r.name}")
                    run_qsite(_fold_r, _inp_r, job_name, rank)
            with ThreadPoolExecutor(max_workers=len(_runjobs), thread_name_prefix=f"QSiteR{rank}") as _qpool:
                list(_qpool.map(_run_one, _runjobs))
        elif not _QSITE_RUN:
            for _fold_r in _folds:
                console_info(f"    [Rank {rank}] QSite run disabled (--no-run-qsite) — inputs in {_fold_r.name}")
        # PHASE 3 — parse every selected frame (freshly run or cached), sequentially.
        _barriers, _derxns = [], []
        for _k, _cand in enumerate(_sel):
            _fold = _folds[_k]
            _res = parse_qsite_barrier(_fold, job_name)
            _b = _res.get("QSite_Barrier_kcal")
            if _b == _b:   # not NaN → a barrier was parsed
                _barriers.append(_b); _derxns.append(_res["QSite_dErxn_kcal"])
            if _k == 0:
                # Primary frame: full reaction profile (PES + departing-fluoride
                # charge) → the direct "did it defluorinate" figure + F-charge cols.
                _prof = parse_qsite_profile(_fold, job_name)
                plot_qsite_reaction_profile(
                    job_out_dir / f"{job_name}_QSite_Reaction_Profile.png", job_name, rank, _prof)
                for _fk in ("F_Charge_Reactant", "F_Charge_Product", "F_Charge_Delta"):
                    stats[_fk] = _prof.get(_fk, np.nan)
        if _barriers:
            '''
            The reported barrier is the RATE-WEIGHTED ENSEMBLE barrier,

                ΔE‡_ens = −RT · ln ⟨ exp(−ΔE‡ᵢ / RT) ⟩ ,

            not the minimum over the scored frames. A minimum is an extreme-value statistic, not a
            property of the ensemble: its downward bias grows with the number of frames that happened
            to parse, so a candidate with four failed scans and one lucky low barrier would outrank a
            candidate with five consistent ones. The bias compounds with the frame SELECTION that
            precedes this, which already favours the most TS-like geometry — taking a minimum
            afterwards maximises the same quantity twice and reports the result as a barrier.

            The exponential average is the right correction because it is what the RATE actually
            averages: the observable is ⟨k⟩ ∝ ⟨exp(−ΔE‡/RT)⟩, and inverting that gives the effective
            barrier the ensemble would exhibit. It is well behaved at both limits — for frames of
            equal barrier it returns that barrier exactly, and when one frame lies far below the rest
            it returns approximately min + RT·ln(N), i.e. it re-applies precisely the penalty that
            the best-of-N search removed. ΔE_rxn is averaged with the SAME Boltzmann weights, so the
            two halves of the verdict describe one ensemble rather than one frame each.

            The minimum is still reported, as QSite_Barrier_Min_kcal, and the number of frames that
            were ATTEMPTED is recorded next to the number that were SCORED — a silent parse failure
            is otherwise indistinguishable from a frame that was never run.
            '''
            _RT_q  = float(CFG.GAS_CONSTANT_KCAL) * float(CFG.MMGBSA_TEMPERATURE_K)
            _b_arr = np.asarray(_barriers, dtype=float)
            _d_arr = np.asarray(_derxns,   dtype=float)
            # Shift by the minimum before exponentiating: exp(-ΔE/RT) underflows for ΔE ≳ 200 kcal,
            # and the shift cancels exactly in the log, so this is algebra, not an approximation.
            _b_min = float(np.min(_b_arr))
            _w     = np.exp(-(_b_arr - _b_min) / _RT_q)
            _b_ens = _b_min - _RT_q * float(np.log(np.mean(_w)))
            _wsum  = float(np.sum(_w))
            _d_ens = float(np.sum(_w * _d_arr) / _wsum) if _wsum > 0 else float(np.mean(_d_arr))

            stats["QSite_Barrier_kcal"]        = round(_b_ens, 2)   # rate-weighted ensemble ΔE‡
            stats["QSite_Barrier_Min_kcal"]    = round(_b_min, 2)   # most accessible single frame
            stats["QSite_Barrier_Mean_kcal"]   = round(float(np.mean(_b_arr)), 2)
            stats["QSite_Barrier_SD_kcal"]     = (round(float(np.std(_b_arr, ddof=1)), 2)
                                                  if len(_b_arr) > 1 else np.nan)
            stats["QSite_dErxn_kcal"]          = round(_d_ens, 2)
            stats["QSite_N_Frames_Scored"]     = int(len(_b_arr))
            stats["QSite_N_Frames_Attempted"]  = int(len(_sel))

            _sd_txt = ("n/a" if len(_b_arr) < 2
                       else f"{stats['QSite_Barrier_SD_kcal']:.1f}")
            print(f"  [Rank {rank}] QM/MM ΔE‡(ensemble) = {_b_ens:.1f} kcal/mol "
                  f"(min {_b_min:.1f}, mean {stats['QSite_Barrier_Mean_kcal']:.1f} ± {_sd_txt}, "
                  f"scored {len(_b_arr)}/{len(_sel)} frames) | "
                  f"ΔE_rxn = {stats['QSite_dErxn_kcal']:.1f} kcal/mol", flush=True)
            if len(_b_arr) < len(_sel):
                console_info(f"    [!] {len(_sel) - len(_b_arr)} of {len(_sel)} QM/MM frame(s) did not "
                             f"yield a barrier — the ensemble average is over the {len(_b_arr)} that did.")

    print(f"  [Rank {rank}] Completed analysis successfully.", flush=True)
    return stats


# =============================================================================
# SECTION 9: MAIN EXECUTION
# =============================================================================

def _resolve_work_dir(raw: str) -> Path:
    """
    Resolves a run dir argument to the 6_Physics_Validation path.

    Handles three input forms:
      1. Exact path to 6_Physics_Validation (already correct).
      2. Exact Boltz-2_Run_* dir (appends 6_Physics_Validation).
      3. Partial prefix such as 'Boltz-2_Run_20260309T085406Z' — glob-expanded to
         match the full timestamped directory on disk.
    """
    p = Path(raw)
    if not p.is_absolute():
        p = Path.cwd() / p

    # Case 1: already points into Physics_Validation
    if "Physics_Validation" in p.parts or p.name == "6_Physics_Validation":
        p.mkdir(parents=True, exist_ok=True)
        return p.resolve()

    # Case 2/3: p is or looks like a Boltz-2_Run_* dir (possibly a prefix)
    if p.exists() and p.is_dir():
        pv = p / "6_Physics_Validation"
        pv.mkdir(parents=True, exist_ok=True)
        return pv.resolve()

    # Partial prefix — glob parent for matching timestamped directories
    parent = p.parent
    stem   = p.name
    matches = sorted(parent.glob(f"{stem}*"),
                     key=lambda x: (x.stat().st_mtime, x.name))   # name breaks an mtime tie deterministically
    if matches:
        pv = matches[-1] / "6_Physics_Validation"
        pv.mkdir(parents=True, exist_ok=True)
        return pv.resolve()

    pv = p / "6_Physics_Validation"
    pv.mkdir(parents=True, exist_ok=True)
    return pv.resolve()


def main():
    # Collapse stacked separator rules: a caller prints a rule, a helper prints its own,
    # and the log grows triple bars with nothing between them.
    _utils_mod.install_console_rule_filter()
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
                        help="Number of ranked jobs (default: auto-detect from 05_MD_Simulations dirs)")
    parser.add_argument("--nuc",    type=int, default=DREAM_TEAM_REF.get('Nuc', 110),  help="Fallback nucleophile resnum")
    parser.add_argument("--base",   type=int, default=DREAM_TEAM_REF.get('Base', 280), help="Fallback base resnum")
    parser.add_argument("--acid",   type=int, default=DREAM_TEAM_REF.get('Acid', 134), help="Fallback acid resnum")
    parser.add_argument("--csv",    default=None,
                        help="Path to master CSV (auto-detected if omitted)")
    parser.add_argument("--workers", type=int, default=None,
                        help="Number of parallel workers (default: auto-detect based on CPU cores)")
    parser.add_argument("--no-run-qsite", action="store_true",
                        help="Only write QSite .in/.mae inputs; do not launch the QSite executable.")
    parser.add_argument("--qsite-procs", type=int, default=None,
                        help="CPUs per QSite job, qsite -PARALLEL (default 1 — the QM engine is single-threaded; scans run concurrently instead).")
    args = parser.parse_args()

    # QSite execution policy: CFG default, overridable per-run from the CLI.
    global _QSITE_RUN, _QSITE_PROCS, _PROGRESS_CR, _QSITE_SEM
    _QSITE_RUN = CFG.QSITE_RUN and not args.no_run_qsite
    # QSite QM/MM (Impact main1h) is single-threaded here, so -PARALLEL >1 only spawns idle
    # helpers: default 1 proc per job, and instead run many scans concurrently (below).
    _QSITE_PROCS = args.qsite_procs if args.qsite_procs is not None else 1
    _QSITE_SEM = threading.Semaphore(_qsite_concurrency())
    console_info(f"QSite concurrency: up to {_qsite_concurrency()} single-threaded QM/MM scans in "
                 f"parallel (cores-2 vs RAM cap), {_QSITE_PROCS} proc/job")

    # Optional sudo up front so the run is fully unattended (systemd-oomd masked for the QSite phase).
    import atexit as _atexit
    _atexit.register(_mask_oomd_at_start())

    raw_dir  = args.run_dir or args.dir or "."
    work_dir = _resolve_work_dir(raw_dir)


    # Auto-detect all available MD rank indices — scan MD, WaterMaps, and
    # 7_MD_Thermodynamics_Results so that every rank already processed is included.
    # Always scan the ranks that actually exist on disk — the MD cohort is SPARSE (whole-library
    # Scientific_Rank, e.g. {1, 2, 8}), so `--ranks N` must not mean the literal 1..N (that silently
    # skips R_8 when N=3). It means the N lowest-numbered ranks that were actually run.
    _auto_rank_list: list[int] = []
    _found_ranks: set[int] = set()
    for _scan_root in [
        work_dir / "05_MD_Simulations",
        work_dir / "03_WaterMaps",
        work_dir.parent / "7_MD_Thermodynamics_Results",
    ]:
        if _scan_root.exists():
            for _d in sorted(_scan_root.iterdir()):
                if _d.is_dir():
                    _m = re.search(r'(?:_R_|[Rr]ank[_\s]?)(\d+)', _d.name)
                    if _m:
                        _found_ranks.add(int(_m.group(1)))
    _avail = sorted(_found_ranks)
    if args.ranks is None:
        _auto_rank_list = _avail
        args.ranks = max(_found_ranks) if _found_ranks else 5
    else:
        _auto_rank_list = _avail[:args.ranks] if _avail else list(range(1, args.ranks + 1))

    master_out_dir = work_dir.parent / "7_MD_Thermodynamics_Results"
    master_out_dir.mkdir(parents=True, exist_ok=True)

    global logger
    logger = (_setup_logging(master_out_dir / "00_MD_Thermodynamics.log",
                             "md_thermo_engine")
              if _setup_logging else None)


    _utils_mod.print_script_banner(
        "07_MD_QMMM_Defluorination_FAcDs.py",
        "MD Thermodynamics  ·  QM/MM Frame Extraction  ·  NAC Validation",
    )
    console_info(f"Run Directory    : {work_dir.parent}")
    console_info(f"Physics Validation : {work_dir}")
    console_info(f"MD Simulations   : {work_dir / '05_MD_Simulations'}")
    console_info(f"WaterMaps        : {work_dir / '03_WaterMaps'}")
    console_info(f"Output           : {master_out_dir}")
    console_info(f"Ligand Resname   : {args.lig}")
    console_info(f"Frame Stride     : {args.stride} (requested){' — all frames' if args.stride == 1 else f' — 1-in-{args.stride} sampled'}")
    console_separator()

    # ── Ranked CSV (6_Boltz2_FAcDs_Ranked_*.csv or any *_Ranked*.csv) ──────────
    prod_dir     = work_dir.parent / "1_Boltz2_Production"
    ranked_csvs  = (sorted(prod_dir.glob(CFG.GLOB_RANKED_CSV)) or
                    sorted(prod_dir.glob("*_Ranked*.csv")))
    if not ranked_csvs:
        console_info(f"{ConsoleColours.FAIL}Error: No ranked CSV found in {prod_dir}{ConsoleColours.ENDC}")
        sys.exit(1)
    # Newest wins; the name breaks a tie, so two files written in the same second cannot swap places
    # between runs.
    ranked_csv_path = sorted(ranked_csvs, key=lambda x: (x.stat().st_mtime, x.name))[-1]
    console_info(f"Ranked CSV       : {ranked_csv_path.name}")
    df_ranked = pd.read_csv(ranked_csv_path, low_memory=False)

    # ── Ranked CSV for triad mapping ──────────────────────────────────────────
    triad_csv = Path(args.csv) if args.csv else ranked_csv_path
    triad_map = {}
    if triad_csv and triad_csv.exists():
        if args.csv and Path(args.csv).resolve() != ranked_csv_path.resolve():
            console_info(f"Triad CSV        : {triad_csv.name}")
        triad_map = load_triad_mapping(triad_csv)
        console_info(f"Triad Entries    : {len(triad_map)} jobs mapped")
    else:
        console_info(f"{ConsoleColours.WARNING}Triad CSV        : Not found — using alignment map / fallback defaults.{ConsoleColours.ENDC}")

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
    # Frame/SN2 analysis runs parallel across ranks at 100% of total_cores-2 (QSite is
    # serialised separately, so it cannot oversubscribe). Explicit --workers overrides.
    _cores_budget  = max(1, (os.cpu_count() or 4) - 2)
    _n_workers     = args.workers if args.workers is not None else max(1, min(len(_rank_list), _cores_budget))
    # Single rank → in-place \r frame counter; multiple parallel ranks → throttled append lines.
    _PROGRESS_CR   = (len(_rank_list) == 1)

    console_info(f"  Effective Stride : {args.stride}{' (all frames)' if args.stride == 1 else f' (1-in-{args.stride})'}")
    console_separator()

    def _process_one_rank(r):
        thread_logger.rank = r
        thread_logger.lines = []

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
        if hasattr(thread_logger, 'lines'):
            _log_lines.extend(thread_logger.lines)
            del thread_logger.lines
        if hasattr(thread_logger, 'rank'):
            del thread_logger.rank
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

    console_info(f"Parallel workers : {_n_workers} (of {len(_rank_list)} ranks) — "
                 f"frame/SN2 analysis parallel @ cores-2={_cores_budget}; QSite QM/MM serialised (1 at a time)")
    print(flush=True)

    _failed_ranks: list  = []
    _all_results:  list  = []
    with ThreadPoolExecutor(max_workers=_n_workers,
                            thread_name_prefix="Rank") as _pool:
        _futures = {_pool.submit(_process_one_rank, r): r for r in _rank_list}
        _completed_count = 0
        for _fut in as_completed(_futures):
            try:
                _rank_num, _rank_log, _rank_res = _fut.result()
                _all_results.append((_rank_num, _rank_log, _rank_res))
                _completed_count += 1
                if sys.stdout.isatty():
                    print(f"\r  [FAcDs Pipeline] Progress: {_completed_count}/{len(_rank_list)} ranks completed...", end="", flush=True)
                else:
                    print(f"  [FAcDs Pipeline] Progress: {_completed_count}/{len(_rank_list)} ranks completed.", flush=True)
            except Exception as _exc:
                _r = _futures[_fut]
                _failed_ranks.append(_r)
                _completed_count += 1
                console_info(f"{ConsoleColours.FAIL}[!] Rank {_r} failed: {_exc}{ConsoleColours.ENDC}")
                import traceback as _tb
                _tb.print_exc()
        if sys.stdout.isatty():
            print("\n", flush=True)

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
        # Stable (mergesort) sort with a deterministic secondary key so equal viability
        # values yield a reproducible Dynamic_Rank across runs (default quicksort is
        # unstable and would shuffle ties non-deterministically).
        _tiebreak = next((c for c in ("Scientific_Rank", "Job", "Job_Name", "Label",
                                      "Protein", "Ligand") if c in df_master.columns), None)
        _sort_cols = ["Catalytic_Viability_Pct"] + ([_tiebreak] if _tiebreak else [])
        _ascending = [False] + ([True] if _tiebreak else [])
        df_master = df_master.sort_values(by=_sort_cols, ascending=_ascending,
                                          kind="mergesort").reset_index(drop=True)
        df_master.insert(0, "Dynamic_Rank", range(1, len(df_master) + 1))

        # Merge Prime MM-GBSA ΔG_bind (Step 06) onto the master by Scientific_Rank so the
        # terminal ranking couples non-covalent binding thermodynamics with the QSite
        # reaction barrier. Both keys are coerced to numeric before the join. A missing or
        # partial summary (Step 06 not run) is non-fatal — the ΔG columns are left absent.
        _mmgbsa_csv = (work_dir / "05_MD_Simulations"
                       / getattr(CFG, "MMGBSA_OUTPUT_SUBDIR", "Prime_MMGBSA")
                       / CFG.FILE_MMGBSA_SUMMARY)
        if _mmgbsa_csv.exists() and "Scientific_Rank" in df_master.columns:
            try:
                _mdf = pd.read_csv(_mmgbsa_csv)
                _keep = [c for c in _mdf.columns if c.startswith("MMGBSA_dG")]
                if "Scientific_Rank" in _mdf.columns and _keep:
                    _mdf = _mdf[["Scientific_Rank"] + _keep].copy()
                    _mdf["Scientific_Rank"] = pd.to_numeric(_mdf["Scientific_Rank"], errors="coerce")
                    df_master["Scientific_Rank"] = pd.to_numeric(df_master["Scientific_Rank"], errors="coerce")
                    df_master = df_master.merge(_mdf, on="Scientific_Rank", how="left")
                    console_info(f"Merged Prime MM-GBSA ΔG_bind for {int(_mdf['Scientific_Rank'].notna().sum())} rank(s).")
            except Exception as _e:
                console_info(f"{ConsoleColours.WARNING}MM-GBSA merge skipped ({_e}).{ConsoleColours.ENDC}")
        elif not _mmgbsa_csv.exists():
            console_info("MM-GBSA summary not found (Step 06 not run) — master written without ΔG_bind columns.")

        """
        ── Defluorination verdict + propensity ───────────────────────────────────────────────────────

            Defluor_Propensity = P(strict-NAC) · exp(−ΔE‡ / RT)

        WHAT THIS IS, AND WHAT IT IS NOT. It is a monotonic RANKING PROXY. It is not a rate, and it must
        not be quoted as one.

        The exponent carries ΔE‡ — the ELECTRONIC barrier QSite returns from the scan. A rate constant
        needs the GIBBS barrier: Eyring's k = (k_B·T/h)·exp(−ΔG‡/RT), and ΔG‡ = ΔE‡ + ZPE + thermal
        corrections − TΔS‡. This pipeline runs no frequency calculation, so it has none of those terms.
        For an enzymatic SN2 the activation entropy is not a rounding error, and the prefactor is absent
        entirely.

        What survives is the ORDERING: for two candidates treated identically, a lower ΔE‡ and a higher
        NAC persistence give a higher propensity, and that is the only claim made of it.
        Defluor_Propensity_Norm rescales it onto 0–1 for readability. The rescaling is min-max in
        LOG space, because the raw quantity is a Boltzmann factor spanning many orders of magnitude;
        Defluor_Propensity_Log10 carries the value it is derived from.

        Is_Defluorinating is the boolean gate: strict persistence AND real dwell AND a surmountable
        barrier AND a non-uphill SN2 product. Thresholds are CFG (SSOT). A missing barrier yields
        "Barrier pending" — the claim is withheld, never converted into a false positive.
        """
        _RT = float(CFG.GAS_CONSTANT_KCAL) * float(CFG.MMGBSA_TEMPERATURE_K)
        def _p_strict(r) -> float:
            return max(0.0, float(r.get("Strict_Viability_Pct", 0.0) or 0.0)) / 100.0
        _prop = []
        for _, _r in df_master.iterrows():
            _bar = _r.get("QSite_Barrier_kcal", np.nan)
            _prop.append(_p_strict(_r) * float(np.exp(-float(_bar) / _RT)) if _bar == _bar else np.nan)
        df_master["Defluor_Propensity"] = _prop

        '''
        The normalisation is done in LOG space, because the quantity being normalised is a Boltzmann
        factor. exp(-ΔE‡/RT) with RT = 0.596 kcal/mol moves about seven orders of magnitude for every
        10 kcal/mol of barrier, so a linear p / max(p) rounds every candidate but the leader to
        0.0000 — an 18 kcal/mol barrier already prints as zero against a 12 kcal/mol one — and the
        landscape figure that uses this as its colour channel becomes a single bright point on a
        uniformly dark field. Ordering is identical either way; log scaling is what makes the
        SPACING between candidates visible, which is the column's only purpose.

        Logs are taken analytically rather than from the exponentiated value: log10(p) is
        log10(p_strict) − ΔE‡ / (RT · ln 10), which stays finite where exp(-ΔE‡/RT) would underflow
        to exactly zero and take log10 to −inf. A candidate with no strict-NAC population has a
        genuinely zero propensity and is floored at the bottom of the scale rather than dropped.
        '''
        _ln10 = float(np.log(10.0))
        _logp = []
        for _, _r in df_master.iterrows():
            _bar = _r.get("QSite_Barrier_kcal", np.nan)
            _ps  = _p_strict(_r)
            if _bar != _bar:
                _logp.append(np.nan)                       # no barrier → no claim
            elif _ps <= 0.0:
                _logp.append(-np.inf)                      # no strict-NAC population → floor
            else:
                _logp.append(float(np.log10(_ps)) - float(_bar) / (_RT * _ln10))
        _fin = [v for v in _logp if v == v and np.isfinite(v)]
        _lo, _hi = (min(_fin), max(_fin)) if _fin else (np.nan, np.nan)
        _span = (_hi - _lo) if (_fin and _hi > _lo) else 0.0
        df_master["Defluor_Propensity_Log10"] = [
            round(v, 3) if (v == v and np.isfinite(v)) else np.nan for v in _logp]
        df_master["Defluor_Propensity_Norm"] = [
            np.nan if v != v else
            (0.0 if not np.isfinite(v) else
             (1.0 if _span <= 0.0 else round((v - _lo) / _span, 4)))
            for v in _logp]

        def _verdict(r):
            """
            The competence claim (the verdict returns 'Defluorination-competent', not an observed
            turnover). Every leg must be MEASURED and must PASS: a quantity that could not
            be computed withholds the verdict, it does not satisfy it. A missing ΔE_rxn is not a
            downhill ΔE_rxn — the thermodynamic leg asserts the SN2 product is not uphill, and an
            absent number is no evidence that it isn't. The barrier and the reaction energy are
            treated identically for that reason; both come from the same QSite parse, and if that
            parse gave only one of them the surviving number cannot carry the other's claim.
            """
            _sv  = float(r.get("Strict_Viability_Pct", 0) or 0)
            _dw  = float(r.get("NAC_Dwell_Max_ns", 0) or 0)
            _bar = r.get("QSite_Barrier_kcal", np.nan)
            _der = r.get("QSite_dErxn_kcal", np.nan)
            _std = int(r.get("MD_Analysis_Stride", 1) or 1)
            if _std > 1:
                # The dwell leg cannot be evaluated on sub-sampled frames, and a defluorination
                # verdict must not be able to change with an I/O performance flag.
                return 0, f"Dwell unmeasurable at stride {_std} — re-run at stride 1"
            if _bar != _bar:
                return 0, "Barrier pending"
            if _der != _der:
                return 0, "Reaction energy pending"
            # Defluorination is, by definition, the C–F bond breaking: the product fluoride must have
            # actually delocalised to free-fluoride (Mulliken charge ≤ QSITE_F_CHARGE_CLEAVED). A low
            # barrier and downhill ΔE_rxn are necessary but not sufficient — a scan can look favourable
            # without releasing the fluoride. Same parse as the barrier, so treat a missing value as
            # pending, not as a fail.
            _fqp = r.get("F_Charge_Product", np.nan)
            if _fqp != _fqp:
                return 0, "Fluoride charge pending"
            if float(_fqp) > CFG.QSITE_F_CHARGE_CLEAVED:
                return 0, "C–F not cleaved (fluoride not released)"
            _ok = (_sv >= CFG.DEFLUOR_STRICT_VIABILITY_MIN_PCT
                   and _dw >= CFG.DEFLUOR_DWELL_MIN_NS
                   and float(_bar) <= CFG.DEFLUOR_BARRIER_MAX_KCAL
                   and float(_der) <= CFG.DEFLUOR_DERXN_MAX_KCAL)
            return (1, "Defluorination-competent") if _ok else (0, "Binds, not competent")
        _v = [_verdict(r) for _, r in df_master.iterrows()]
        df_master["Is_Defluorinating"] = [x[0] for x in _v]
        df_master["Defluor_Verdict"]   = [x[1] for x in _v]
        # Turnover ranking (highest propensity = rank 1; unscored ranks sort last).
        df_master["Defluor_Rank"] = (
            df_master["Defluor_Propensity"].rank(ascending=False, method="min", na_option="bottom").astype("Int64"))
        _ncomp = int(sum(x[0] for x in _v))
        console_info(f"Defluorination-competent candidates : {_ncomp} / {len(df_master)} "
                     f"(gate: strict≥{CFG.DEFLUOR_STRICT_VIABILITY_MIN_PCT}%, "
                     f"dwell≥{CFG.DEFLUOR_DWELL_MIN_NS} ns, ΔE‡≤{CFG.DEFLUOR_BARRIER_MAX_KCAL}, "
                     f"ΔE_rxn≤{CFG.DEFLUOR_DERXN_MAX_KCAL} kcal/mol)")

        master_csv_path = master_out_dir / "01_MD_Master_Ranking.csv"
        _utils_mod.atomic_write_csv(df_master, master_csv_path)
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
                    format_job_label(row["Job_Name"], row["Scientific_Rank"]),
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
            ("Reactive-Pose Figures (MM-GBSA decomposition · machinery engagement)...",
             generate_reactive_pose_figures),
            ("Comparative Residue Engagement (all SN2 cases)...", generate_comparative_residue_engagement),
            ("MD Viability & Retention Bar Chart...", generate_viability_bar_chart),
            ("Defluorination Landscape (persistence × barrier × binding)...", generate_defluorination_landscape),
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
    _utils_mod.print_elapsed(_t0, "07_MD_QMMM_Defluorination_FAcDs.py")
