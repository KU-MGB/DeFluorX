#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 03  |  Validation, Ranking & Visualisation
===============================================================================
The definitive "Judge": merges physics-based structural validation with
confidence metrics from Boltz-2, performs Pareto optimisation, rescues
hidden-gem candidates, and generates a publication figure suite organised
folder-by-folder (seven numbered folders: Ramachandran, dataset & alignment
overview, AI prediction quality, catalytic geometry & mechanism, ligand
interactions & chemical space, PFAS scope & synthesis, and pocket-fit / multi-model
diagnostics).

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 05 July 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 03_Validation_Figures_FAcDs.py
  Role          : Scoring, ranking, and visual reporting of Boltz-2 predictions.
  Imports from  : 00_02_Project_Config_FAcDs.py  (CFG — tier colours, vis params)
                  00_03_Project_Utils_FAcDs.py   (ConsoleColours, setup_logging,
                                                 console_info, console_separator)
  Reads         : <Run>/1_Boltz2_Production/*_Ranked_*.csv  (falls back to *_Master_*.csv)
                  (Master CSV written by 02_Production_FAcDs.py; latest file selected)
  Writes        : <Run>/3_Validation_Figures/  (figures grouped folder-by-folder)
                    01_Ramachandran/                        Ramachandran_*.png
                    02_Dataset_and_Alignment_Overview/      01_*.png … 04_*.png
                    03_AI_Confidence_Quality/               01_*.png … 03_*.png
                    04_Catalytic_Geometry_and_Mechanism/    01_*.png … 10_*.png
                    05_Ligand_Interactions_and_Chemical_Space/ 01_*.png … 08_*.png
                    06_PFAS_Scope_and_Synthesis/            01_*.png … 13_*.png
                    07_Diagnostic_and_MultiModel_Trends/    01_*.png … 09_*.png  (pocket-fit + consensus + competence)
                  <Run>/3_Validation_Figures/03_Final_Validated_Master.csv
                  <Run>/3_Validation_Figures/04_ACTION_Rescue_Hidden_Gems.csv
                  <Run>/3_Validation_Figures/05_Figure_Descriptions.txt
                  <Run>/3_Validation_Figures/06_Analysis_Log.txt
  Upstream      : 02_Production_FAcDs.py → writes the master ranked CSV (incl. the pocket-fit
                  columns active_site_volume, ligand_volume, pocket_occupancy, fit_ratio,
                  ligand_fits) consumed here
  Downstream    : 04_Dendrogram_FAcDs.py  → reads 03_Final_Validated_Master.csv
───────────────────────────────────────────────────────────────────────────────

── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Memory overhead: Generating 24+ high-resolution figures simultaneously can
     spike RAM usage; recommended 32GB+ for large (>10k row) datasets.
  2. CSV Schema Sensitivity: Relies on the exact column naming convention from
     Step 02; custom CSV modifications will break the scoring logic.
  3. Pareto Convergence: Non-dominated sorting complexity scales O(N log N);
     runs with extremely high objective-conflict may take several minutes.
  4. Visualisation Headless: Matplotlib MUST use the 'Agg' backend; the script
     forces this, but environment-specific Tkinter conflicts may still occur.
───────────────────────────────────────────────────────────────────────────────

Usage:
    -------------------------------------------------------------------------------
    1. Activate Environment
       conda activate PFAS

    2. Run Validation (Auto-detects latest run if argument omitted)
       python 03_Validation_Figures_FAcDs.py Boltz-2_Run_20260309T085406Z
    -------------------------------------------------------------------------------

-------------------------------------------------------------------------------
Purpose:
    The definitive "Judge" of the Boltz-2 pipeline. Merges advanced data-driven
    validation with comprehensive visual reporting.

    1. Audits Physics (Geometry) vs. AI (Confidence).
    2. Performs Pareto Optimisation (Non-dominated sorting) to find optimal trade-offs.
       * NOTE: Uses Fast Sort-and-Sweep algorithm for large datasets (>50k rows).
    3. Rescues "Hidden Gems" that AI missed but Physics loves.
    4. Generates a complete suite of high-resolution scientific figures.

    Phylogenetic analysis is handled by the downstream 04_Dendrogram_FAcDs.py script.

-------------------------------------------------------------------------------
Outputs (Saved in <Run_Folder>/3_Validation_Figures/):
    [Data — written to the folder root]
    • 03_Final_Validated_Master.csv         <-- THE FINAL DATASET
    • 04_ACTION_Rescue_Hidden_Gems.csv      <-- MANUAL REVIEW LIST
    • 05_Figure_Descriptions.txt            <-- Per-figure description log
    • 06_Analysis_Log.txt                   <-- Detailed Execution Log

    [Figures — grouped folder-by-folder; each folder numbered in narrative order]
    ── 01_Ramachandran/ ── backbone-geometry validation of the controls
    • Ramachandran_3R3U_Crystal.png · Ramachandran_{DeHa4,3R3U}_Control.png (+ _vs_Crystal)

    ── 02_Dataset_and_Alignment_Overview/ ── dataset + sequence-alignment overview
    • 01_Active_Site_Residue_Mapping_Coverage.png  <-- Per-residue alignment coverage (+ MISSED.csv)
    • 02_Tier_Distribution.png                     <-- Catalytic tier counts + model-selection pie
    • 03_Alignment_Grades.png                      <-- Sequence-identity grades + KDE by tier
    • 04_Tier_Grade_Distribution.png               <-- Tier × alignment-grade cross-tabulation

    ── 03_AI_Confidence_Quality/ ── Boltz-2 confidence metrics
    • 01_AI_Quality_Assessment.png                 <-- Boltz confidence boxes + pTM/ipTM panel
    • 02_AI_Quality_Space.png                      <-- {CFG.TIER_TOP} pTM/ipTM density + thumbnails
    • 03_pTM_vs_ipTM_by_Tier.png                   <-- pTM vs ipTM 2-D scatter by tier

    ── 04_Catalytic_Geometry_and_Mechanism/ ── structural + mechanistic geometry
    • 01_ActiveSite_RMSD_by_Tier.png               <-- Active-site RMSD vs reference control (median bar)
    • 02_Mech_State_CrossTab.png                   <-- Halide stabilisation × carboxylate clamp heatmap
    • 03_Mechanistic_Score_by_Tier.png             <-- Mechanistic score mean ± CI (median annotated)
    • 04_SN2_Angle_by_Tier.png                     <-- SN2 attack-angle ECDF by tier
    • 05_Mechanism_Geometry_Scatter.png            <-- SN2 angle vs nucleophile-distance scatter
    • 06_Mechanistic_Quality_Space.png             <-- {CFG.TIER_TOP} mechanistic density + thumbnails
    • 07_Feature_Correlations.png                  <-- Spearman ρ feature correlation heatmap
    • 08_Tier_Quality_DotPlot.png                  <-- Multi-metric tier-quality Cleveland dot plot
    • 09_Mechanistic_Fingerprint.png               <-- Per-tier catalytic feature profile (radar / spider; config §5 components)
    • 10_TwoCriteria_Tier_Logic.png                <-- 3-panel: Criterion A gates B (a) · B-ECDF separates tiers (b) · SN2 dead-end BDE×occlusion gate (c)

    ── 05_Ligand_Interactions_and_Chemical_Space/ ── interaction profile + chemical space
    • 01_Molecular_Interaction_Profile.png         <-- Bond-type profile + fluorine engagement
    • 02_Interaction_Quality_Space.png             <-- {CFG.TIER_TOP} interaction density + thumbnails
    • 03_Fluorine_Engagement_by_Tier.png           <-- Fluorine engagement ratio box + trend line
    • 04_Catalytic_Quality_vs_Inhibition.png       <-- Soft-catalytic-score violin + active-site contact density
    • 05_ActiveSite_Contact_Density_by_Tier.png    <-- Active-site contact density box + strip
    • 06_Chemical_Space_Map.png                    <-- UMAP chemical-space manifold by tier
    • 07_Chemical_Space_Landscape.png              <-- {CFG.TIER_TOP} KDE density + structure thumbnails
    • 08_Binding_Energetics.png                    <-- Binding-probability violin by tier

    ── 06_PFAS_Scope_and_Synthesis/ ── multi-metric synthesis + publication assembly
    • 01_Radar_TopHits.png                         <-- Radar: top-5 hits vs worst-tier baseline
    • 02_Radar_TierReps.png                        <-- Radar: one representative per tier
    • 03_Tier_Success_Rates.png                    <-- Tier % success: substrate vs inhibitor geometry
    • 04_Conf_SN2_Landscape.png                    <-- Per-tier medians in confidence × SN2-angle space
    • 05_Conflict_Composition.png                  <-- Conflict category × tier stacked bar
    • 06_Hidden_Gems_DeepDive.png                  <-- Physics-good / AI-missed slope chart
    • 07_Category_Overlap_Euler.png                <-- Multi-set overlap Euler diagram
    • 08_Top25_Multitarget_Proteins.png            <-- Top-25 multitarget proteins, stacked by tier
    • 09_TopTier_Protein_PFAS_Breakdown.png        <-- Top-tier proteins: PFAS degraded (colour=ligand)
    • 10_Sankey_Workflow.png                       <-- Sankey: nuc-dist → SN2 angle → final tier
    • 11_PFAS_Size_Hexbin_Landscape.png            <-- Chain-length hexbin + tier scatter + rolling median
    • 12_PFAS_Size_Composition_Merged.png          <-- Per bin: outcome + degrader-tier stacked bars
    • 13_PFAS_Carbon_Confidence.png             <-- Catalytic competence vs MW per carbon group

    ── 07_Diagnostic_and_MultiModel_Trends/ ── pocket-fit + multi-model consensus diagnostics
    • 01_Pocket_vs_Ligand_Volume.png               <-- Cavity vs ligand volume, y=x steric fit boundary
    • 02_Pocket_Occupancy_by_Carbon_Number.png     <-- Occupancy violin+box per PFAS carbon number, median trend
    • 03_Occupancy_vs_Competence.png               <-- Competence vs occupancy + binned-median trend
    • 04_Ligand_Fit_Rate_by_Ligand.png             <-- % complexes passing the steric fit test, per ligand
    • 05_MultiModel_Consensus_by_Tier.png          <-- Multi-model degrader consensus, mean ± CI by tier
    • 06_Confidence_vs_Consensus.png               <-- Confidence vs cross-model consensus hexbin + trend
    • 07_Quality_and_Competence_Diagnostics.png    <-- Confidence×competence + mech-score×penalty scatter (tier)
    • 08_Size_Preference_Containment.png           <-- Effective-mech distribution + means + hit-rate + pocket containment vs ligand size
    • 09_Reactive_Engagement.png                   <-- Reactive-C→catalytic-residue distance + properly-positioned fraction (vs hit-rate) by carbon number

-------------------------------------------------------------------------------
Scientific References:
    1. Data handling & numerics:
       - Harris, C.R. et al. (2020) NumPy. Nature 585:357–362. DOI: https://doi.org/10.1038/s41586-020-2649-2
       - McKinney, W. (2010) Data Structures for Statistical Computing in Python. Proc 9th Python in Science Conf 56–61. DOI: https://doi.org/10.25080/Majora-92bf1922-00a
       - Virtanen, P. et al. (2020) SciPy 1.0. Nature Methods 17:261–272. DOI: https://doi.org/10.1038/s41592-019-0686-2
    2. Plotting:
       - Hunter, J.D. (2007) Matplotlib. Comput Sci Eng 9:90–95. DOI: https://doi.org/10.1109/MCSE.2007.55
       - Waskom, M.L. (2021) seaborn. J Open Source Softw 6:3021. DOI: https://doi.org/10.21105/joss.03021
    3. Dimensionality reduction & feature scaling:
       - McInnes, L., Healy, J. & Melville, J. (2018) UMAP. J Open Source Softw 3:861. DOI: https://doi.org/10.21105/joss.00861
       - Pedregosa, F. et al. (2011) scikit-learn: Machine Learning in Python. J Mach Learn Res 12:2825–2830.
    4. Image composition (Pillow / PIL fork):
       - Clark, A. (2015) Pillow. https://python-pillow.org
    5. Molecular rendering (optional structure figures):
       - The PyMOL Molecular Graphics System, Schrödinger, LLC. https://pymol.org
    6. All scientific thresholds/criteria plotted here are defined in
       00_02_Project_Config_FAcDs.py — see that module's Scientific References
       for the underlying primary literature (NAC, Maestro criteria, mech score, etc.).
===============================================================================
"""

# =============================================================================
# SECTION 1: SYSTEM CONFIGURATION & IMPORTS
# =============================================================================

# -------------------------------------------------------------------------------
# Step 1.1: Standard Library Imports
# -------------------------------------------------------------------------------
import sys
import shutil
import argparse
import warnings
from pathlib import Path
import contextlib
from datetime import datetime

# -------------------------------------------------------------------------------
# Step 1.2: Scientific Stack Imports
# -------------------------------------------------------------------------------
"""
CPU usage cap (total cores − 2; mirrors CFG.PREP_CPU_RESERVE). Reserve 2 cores
for OS/desktop stability by limiting the thread-pool maths libraries (BLAS /
MKL / OpenMP / NumExpr). Must precede numpy/scipy import to take effect;
setdefault() preserves any value exported by the caller or pipeline runner.
"""
import os as _os
_CPU_CAP = str(max(1, (_os.cpu_count() or 4) - 2))
for _tv in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_tv, _CPU_CAP)

import numpy as np
import pandas as pd
import seaborn as sns
import umap
import matplotlib
matplotlib.use("Agg")   # non-interactive backend — renders on headless HPC/GPU nodes (no DISPLAY)
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, gaussian_kde
from sklearn.preprocessing import MinMaxScaler, RobustScaler, StandardScaler
from sklearn.decomposition import PCA

# Bioinformatics Imports
from matplotlib.patches import ConnectionPatch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from PIL import Image, ImageDraw, ImageFont


# -------------------------------------------------------------------------------
# Step 1.3: Pipeline modules (00_02) via importlib
# -------------------------------------------------------------------------------
import importlib.util as _ilu

def _load_module(name: str, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod  = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_cfg_mod   = _load_module("ProjectConfig", Path(__file__).resolve().parent / "00_02_Project_Config_FAcDs.py")
_utils_mod = _load_module("ProjectUtils",  Path(__file__).resolve().parent / "00_03_Project_Utils_FAcDs.py")
CFG        = _cfg_mod.CFG()

ConsoleColours  = _utils_mod.ConsoleColours
SEPARATOR_HEAVY = _utils_mod.SEPARATOR_HEAVY
SEPARATOR_LIGHT = _utils_mod.SEPARATOR_LIGHT
_setup_logging  = _utils_mod.setup_logging
_console_title  = _utils_mod.console_title
_console_info   = _utils_mod.console_info
_console_sep    = _utils_mod.console_separator

# -------------------------------------------------------------------------------
# Step 1.4: Global Configuration
# -------------------------------------------------------------------------------

# Suppress warnings for clean output
warnings.filterwarnings("ignore")

# Weighted Scoring (Used for "Ensemble Score" heuristic fallback)
WEIGHTS = CFG.VIS_FALLBACK_WEIGHTS

"""
Tier palette — sourced from CFG.TIER_COLOUR (Okabe-Ito colourblind-safe).
Local overrides: CFG.TIER_DECOY uses a softer grey for unlabelled entries;
"Error" is a 03-specific indicator for analytics failures (not in CFG).
"""
TIER_PALETTE = {
    **CFG.TIER_COLOUR,                 # decoy colour comes from CFG.TIER_COLOUR (single source)
    "Error": "#FF6B6B",                # 03-specific analytics-failure indicator (not in CFG)
}

# Specific order for tiers to ensure logical plotting (Tier_1A to Tier_5_Decoy)
TIER_ORDER_LOGIC = list(CFG.TIER_ORDER)

CONFLICT_PALETTE = dict(CFG.CONFLICT_COLOUR)   # sourced from CFG § 8.7

# MD-ready cohort — the "super-best" MD-selected complexes, overlaid as gold stars on
# tier plots (they may sit inside Tier_1A). Single source for the marker + the subset.
_MD_STAR_KW = dict(marker="*", s=300, facecolor="#FFD400", edgecolor="black", linewidths=1.1, zorder=9)

def _md_ready_df(df: pd.DataFrame) -> pd.DataFrame:
    """Rows flagged MD_Selected — the super-best cohort taken to MD."""
    if "MD_Selected" not in df.columns:
        return df.iloc[0:0]
    _m = df["MD_Selected"].astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes"])
    return df[_m]

def _md_ready_stars_cat(ax, df: pd.DataFrame, tier_order, valcol: str):
    """Overlay the MD-ready cohort as gold stars on a categorical (per-tier) axis where
    x = tier index and y = valcol; same-tier stars are x-jittered so they never overlap.
    Returns a Line2D legend handle (to append to the figure's handle list), or None."""
    _md = _md_ready_df(df)
    _pos = {t: i for i, t in enumerate(tier_order)}
    _drawn = 0
    for _t, _grp in _md.groupby("degrader_tier"):
        if _t not in _pos:
            continue
        _vals = pd.to_numeric(_grp[valcol], errors="coerce").dropna().values
        _xs = np.linspace(-0.18, 0.18, len(_vals)) if len(_vals) > 1 else [0.0]
        for _dx, _v in zip(_xs, _vals):
            ax.scatter([_pos[_t] + _dx], [_v], **_MD_STAR_KW)
            _drawn += 1
    if not _drawn:
        return None
    return Line2D([0], [0], marker="*", ls="", markerfacecolor="#FFD400",
                  markeredgecolor="black", markersize=12, label=f"MD-ready (n={_drawn})")


def _auto_label_colour(bg, threshold: float = 0.5) -> str:
    """Single source for every auto-contrast label colour in this module.

    Returns CFG.VIS_BAR_LABEL_COLOURS_ON_LIGHT[0] (near-black) on light fills and
    CFG.VIS_BAR_LABEL_COLOURS_ON_DARK[0] (white) on dark fills, chosen from the WCAG
    relative luminance of `bg` (a hex string or any Matplotlib colour). The single
    luminance helper shared by every figure — one coefficient set and cutoff.
    """
    import matplotlib.colors as _mc
    try:
        r, g, b = _mc.to_rgb(bg)
    except Exception:
        return CFG.VIS_BAR_LABEL_COLOURS_ON_LIGHT[0]
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return (CFG.VIS_BAR_LABEL_COLOURS_ON_LIGHT[0] if lum > threshold
            else CFG.VIS_BAR_LABEL_COLOURS_ON_DARK[0])

# Modern Plotting Theme (Publication Quality)
sns.set_theme(style="whitegrid", context="paper", font_scale=1.3)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], # Nature preference
    "axes.grid": True,
    "grid.color": "#F0F0F0",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": CFG.VIS_FIGURE_DPI,
    "savefig.bbox": "tight"
})

SEPARATOR = "-" * 80


# =============================================================================
# SECTION 2: LOGGING & UTILITIES
# =============================================================================

logger = None

def console_title(msg: str) -> None:
    _console_title(msg, logger)

def console_info(msg: str) -> None:
    _console_info(msg, logger)

def console_separator() -> None:
    _console_sep(logger, heavy=True)

def _aux_dir(out_dir: Path) -> Path:
    """Subfolder for the non-figure data deliverables (CSV/TXT), keeping the
    3_Validation_Figures root clutter-free — only figure folders + this one."""
    d = out_dir / "00_Analysis_Data"
    d.mkdir(parents=True, exist_ok=True)
    return d


ReportManager = _utils_mod.ReportManager   # shared logger (00_03)


def _make_reporter(out_dir: Path):
    """Construct the shared ReportManager with this step's log path/header/logger."""
    return ReportManager(
        _aux_dir(out_dir) / "06_Analysis_Log.txt",
        "BOLTZ-2 VALIDATION & PHYLOGENY REPORT ",
        separator=SEPARATOR_LIGHT, rule_width=80, log_fn=console_info)

# Alignment grade is provided by _utils_mod.get_alignment_grade(val, CFG).


# =============================================================================
# SECTION 3: ANALYSIS ALGORITHMS (PARETO OPTIMISATION)
# =============================================================================

# -------------------------------------------------------------------------------
# Step 3.1: Pareto Optimisation
# -------------------------------------------------------------------------------
def calculate_pareto_fronts(df: pd.DataFrame, objectives: list, maximize: list) -> pd.Series:
    """ Identifies non-dominated sorting fronts (Pareto Frontiers). """
    data = df[objectives].copy().values
    n_points = data.shape[0]

    weighted_data = data.copy()
    for i, is_max in enumerate(maximize):
        if is_max: weighted_data[:, i] = -weighted_data[:, i]

    pareto_ranks = np.full(n_points, 999, dtype=int)
    population_ids = np.arange(n_points)

    MAX_RANKS = CFG.VIS_MAX_RANKS_DISPLAY
    current_rank = 1

    while len(population_ids) > 0 and current_rank <= MAX_RANKS:
        subset = weighted_data[population_ids]

        if subset.shape[1] == 2:
            sorted_args = np.lexsort((subset[:, 1], subset[:, 0]))
            sorted_ids = population_ids[sorted_args]
            sorted_data = subset[sorted_args]

            min_y_seen = float("inf")
            current_front_indices = []
            _last_xy = (None, None)   # last point kept on this front (for exact-duplicate ties)

            for i in range(len(sorted_data)):
                x_val = sorted_data[i, 0]
                y_val = sorted_data[i, 1]
                # Strictly-better points join the front; exact duplicates of the last
                # kept point share its rank (consistent with the M>2 dominance check).
                if y_val < min_y_seen or (x_val == _last_xy[0] and y_val == _last_xy[1]):
                    if y_val < min_y_seen:
                        min_y_seen = y_val
                    current_front_indices.append(sorted_ids[i])
                    _last_xy = (x_val, y_val)

            pareto_ranks[current_front_indices] = current_rank
            mask_keep = np.isin(population_ids, current_front_indices, invert=True)
            population_ids = population_ids[mask_keep]

        else:
            """
            Vectorised dominance check — no Python loop.
            boe[i,j,m]: subset[j,m] <= subset[i,m]  (j better-or-equal to i on obj m)
            b[i,j,m]:   subset[j,m] <  subset[i,m]  (j strictly better on obj m)
            j dominates i iff all-objectives boe AND any-objective b.
            """
            _boe      = subset[None, :, :] <= subset[:, None, :]   # (N, N, M)
            _b        = subset[None, :, :] <  subset[:, None, :]   # (N, N, M)
            _j_dom_i  = np.all(_boe, axis=2) & np.any(_b, axis=2) # (N, N)
            np.fill_diagonal(_j_dom_i, False)
            is_dominated = np.any(_j_dom_i, axis=1)                # (N,)

            local_front_indices = np.where(~is_dominated)[0]
            global_front_indices = population_ids[local_front_indices]
            pareto_ranks[global_front_indices] = current_rank
            population_ids = np.delete(population_ids, local_front_indices)

        current_rank += 1

    return pd.Series(pareto_ranks, index=df.index)

# =============================================================================
# SECTION 4: DATA PROCESSING & FIGURES (PIPELINE)
# =============================================================================

def load_and_prep_data(prod_dir: Path, reporter: ReportManager) -> tuple[pd.DataFrame, list[str]]:
    """
    Prefer ranked CSV (has Scientific_Rank + Ranking_Score_Calc used by downstream steps).
    Fall back to master CSV if ranked not yet generated.
    """
    candidates = sorted(list(prod_dir.glob("*_Ranked_*.csv")))
    if not candidates: raise FileNotFoundError("No Ranked CSV found in Production folder.")

    target = candidates[-1]
    reporter.log(f"Input Data Source: {target.resolve()}")
    df = pd.read_csv(target, low_memory=False)

    if "job_name" in df.columns:
        df = df.drop_duplicates(subset=["job_name"], keep="last")

    tier_map = CFG.TIER_RANK
    df["tier_numeric"] = df["degrader_tier"].map(tier_map).fillna(0)

    col_map = {
        "binding_likelihood_computed": "Binding_Probability",
        "Binding_Probability_Score": "Binding_Probability",
        "custom_affinity_score": "Chemical_Affinity_Score",
        "confidence_score": "Boltz_Model_Confidence",
        "interaction_density": "Interaction_Density_Norm",
        "sn2_attack_angle": "SN2_Attack_Angle",
        "Active_Site_RMSD_to_Control": "Active_Site_RMSD",
    }
    df.rename(columns=col_map, inplace=True)

    # tier_numeric is the target label — it is deliberately NOT an input feature,
    # else PCA/UMAP would align their principal axes with the tier and manufacture
    # a circular, guaranteed tier separation (target leakage). It is used only to
    # orient the PC1 sign after projection and to colour points post hoc.
    candidate_features = [
        "Boltz_Model_Confidence", "iptm", "mean_plddt",
        "Interaction_Density_Norm", "Chemical_Affinity_Score", "Binding_Probability",
        "count_hydrogen_bond", "count_salt_bridge"
    ]

    valid_features = []
    for f in candidate_features:
        if f in df.columns:
            df[f] = pd.to_numeric(df[f], errors="coerce").fillna(0)
            if df[f].std() > 0.001: valid_features.append(f)

    reporter.log(f"Loaded {len(df)} complexes. Valid Features for Analysis: {len(valid_features)}")
    return df, valid_features

def _median_ci95(vals, n_boot: int = 1000, seed: int = 99):
    """Bootstrap 95% confidence interval of the MEDIAN (percentile method).

    Returns (lo, hi). Used for error bands drawn around a plotted median: an
    SEM·t interval is only valid around a mean and misrepresents the median on
    skewed distributions. Falls back to a degenerate (m, m) for n < 2.
    """
    v = pd.to_numeric(pd.Series(vals), errors="coerce").dropna().to_numpy(dtype=float)
    if v.size < 2:
        m = float(v[0]) if v.size else 0.0
        return m, m
    rng = np.random.default_rng(seed)
    meds = np.median(rng.choice(v, size=(n_boot, v.size), replace=True), axis=1)
    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))


def perform_advanced_ranking(df: pd.DataFrame, features: list[str], out_dir: Path, reporter: ReportManager):
    """Multi-objective ranking of complexes.

    Scales the feature matrix (robust, with a StandardScaler fallback on zero-IQR
    columns), derives a PCA PC1 score (sign-aligned to tier_numeric), a CFG-weighted
    composite score, Pareto fronts (confidence vs binding probability) and a UMAP
    embedding. Adds the Score/Ensemble/Pareto/UMAP columns to df in place and
    returns it.
    """
    reporter.section("Step 1/5 — Multi-Objective Ranking (PCA & Pareto)  [analysis, no figure folder]")

    x = df[features].dropna()
    '''
    Guard the scaler: a feature whose inter-quartile range collapses to 0
    (heavily skewed or constant tail) makes RobustScaler divide by zero and
    propagate inf/NaN into PCA and UMAP. Fall back to StandardScaler when any
    column has zero IQR (it divides by standard deviation instead).
    '''
    _iqr = x.quantile(0.75) - x.quantile(0.25)
    if x.shape[1] > 0 and (_iqr <= 0).any():
        x_scaled = StandardScaler().fit_transform(x)
    else:
        x_scaled = RobustScaler().fit_transform(x)

    # Default in case PCA is skipped (len(x) < 2)
    score_pca = np.full(len(x), 50.0)

    # PCA
    if len(x) >= 2:
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(x_scaled)

        # Orient PC1 by correlating with an unsupervised INPUT feature that tracks catalytic
        # quality (mechanistic_score → competence_score → first feature), NOT the supervised
        # tier label. PCA sign is mathematically arbitrary; fixing it against an input feature
        # keeps the "high PC1 = better" reading for the exploratory map without rotating the
        # latent space to the human-assigned category (no target leakage).
        check_col = next((c for c in ("mechanistic_score", "competence_score") if c in x.columns), features[0])
        corr, _ = spearmanr(pcs[:, 0], df.loc[x.index, check_col])

        score_pca = MinMaxScaler(feature_range=(0, 100)).fit_transform(pcs[:, 0].reshape(-1, 1)).flatten()
        if corr < 0:
            score_pca = 100.0 - score_pca
            reporter.log(
                f"  -> PC1 direction: Spearman ρ(PC1, {check_col}) = {corr:.4f} — NEGATIVE."
                f"  PC1 scores mathematically inverted so high PC1 = higher catalytic quality."
                f"  [Methods note: PC1 sign chosen by correlation with the input feature "
                f"  {check_col} (unsupervised); inversion applied when ρ < 0.]"
            )
        else:
            reporter.log(
                f"  -> PC1 direction: Spearman ρ(PC1, {check_col}) = {corr:.4f} — POSITIVE."
                f"  PC1 preserved (high PC1 = higher catalytic quality, no inversion needed)."
            )

        df.loc[x.index, "Score_PCA_Data"] = score_pca
    else:
        reporter.log("  ! Skipping PCA: Insufficient data points (< 2).")
        df.loc[x.index, "Score_PCA_Data"] = 50.0

    # Weighted Score
    score_weighted = np.zeros(len(x))
    total_w = 0.0
    for col, w in WEIGHTS.items():
        if col in df.columns:
            _cv = pd.to_numeric(df.loc[x.index, col], errors="coerce").to_numpy(dtype=float)
            # Percentile-robust min-max (1st-99th): a single extreme value does not
            # compress the rest of the column into a narrow band (unlike a plain
            # outlier-sensitive MinMaxScaler), while keeping the [0,1] scale.
            _lo, _hi = np.nanpercentile(_cv, 1), np.nanpercentile(_cv, 99)
            norm_col = (np.clip((_cv - _lo) / (_hi - _lo), 0.0, 1.0)
                        if _hi > _lo else np.zeros_like(_cv))
            score_weighted += norm_col * w
            total_w += w
    if total_w > 0: score_weighted = (score_weighted / total_w) * 100
    df.loc[x.index, "Score_Weighted"] = score_weighted

    # Pareto — a pure confidence-vs-binding front optimises for tight-binding high-confidence
    # poses and promotes catalytically dead competitive inhibitors as "optimal degraders". The
    # front must therefore see a catalytic (NAC-viability) axis. Rather than add a third
    # objective (which forces the O(N²) dense-dominance path → tens of GB on the full frame),
    # gate the CANDIDATES on catalytic viability first — a mechanistic-score floor (falling back
    # to the is_degrader flag) — then run the fast 2-objective non-dominated sort on the viable
    # subset only. Non-viable poses are held off the front (rank 999), so no inhibitor surfaces
    # as an optimal degrader while the O(N log N) scalability is preserved.
    if "Boltz_Model_Confidence" in df.columns and "Binding_Probability" in df.columns:
        _mech_obj = next((c for c in ("mechanistic_score_effective", "mechanistic_score") if c in df.columns), None)
        if _mech_obj:
            _viable = pd.to_numeric(df[_mech_obj], errors="coerce").fillna(0.0) >= CFG.TIER_MECH_MIN["Tier_2A"]
        elif "is_degrader" in df.columns:
            _viable = df["is_degrader"].astype(bool)
        else:
            _viable = pd.Series(True, index=df.index)
        reporter.log(
            f"  -> Calculating Pareto Frontiers (Confidence vs Binding Probability) over "
            f"{int(_viable.sum())} catalytically-viable poses"
            + (f" ({_mech_obj} ≥ {CFG.TIER_MECH_MIN['Tier_2A']})..." if _mech_obj else " (is_degrader)...")
        )
        df["Pareto_Rank"] = 999
        if _viable.any():
            _pr = calculate_pareto_fronts(df.loc[_viable], ["Boltz_Model_Confidence", "Binding_Probability"], [True, True])
            df.loc[_viable, "Pareto_Rank"] = _pr
    else:
        df["Pareto_Rank"] = 0

    # Ensemble_Score is exploratory only (PCA/weighted blend for colouring and
    # chemical-space maps). The authoritative rank downstream is the ranked CSV's
    # Scientific_Rank / competence_score from Step 02; Ensemble_Score is never used
    # as a sort or rank key.
    df.loc[x.index, "Ensemble_Score"] = (0.5 * score_pca) + (0.5 * score_weighted)
    df["Ensemble_Score"] = df["Ensemble_Score"].fillna(0)

    # UMAP
    if len(x) >= 5:
        try:
            reporter.log("  -> Computing UMAP Manifold for Chemical Space Map...")
            reducer = umap.UMAP(n_neighbors=min(15, len(x)-1), min_dist=0.1, random_state=42)
            umap_map = reducer.fit_transform(x_scaled)
            df.loc[x.index, "UMAP_X"] = umap_map[:, 0]
            df.loc[x.index, "UMAP_Y"] = umap_map[:, 1]
        except Exception as e:
            reporter.log(f"  ! UMAP Calculation Skipped: {e}")
    else:
        reporter.log("  ! Skipping UMAP: Insufficient data points (< 5).")

    if len(x) >= 2:
        pd.DataFrame(pca.components_.T, columns=["PC1", "PC2"], index=features).to_csv(_aux_dir(out_dir) / "02_PCA_Loadings.csv", index=True)

    return df

def analyse_conflicts(df: pd.DataFrame, out_dir: Path, reporter: ReportManager):
    """Flag conflict/opportunity cases and write the validated master CSV.

    Classifies each complex (e.g. hidden gems = high mechanistic merit but
    low confidence; decoys = high confidence but poor mechanism), logs the
    per-class counts, and saves the final validated master table to out_dir.
    """
    reporter.section("Step 2/5 — Conflict & Opportunity Analysis  [analysis, no figure folder]")

    def classify(row):
        tier = row.get("degrader_tier", CFG.TIER_DECOY)
        conf = row.get("Boltz_Model_Confidence", 0.0)

        _high_quality = [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]
        _hi = CFG.CONFLICT_CONF_HIGH
        _lo = CFG.CONFLICT_CONF_LOW
        if tier in _high_quality and conf >= _hi: return "Consensus High"
        if tier in [CFG.TIER_POOR, CFG.TIER_DECOY, "Error"] and conf < _lo: return "Consensus Low"
        if tier in _high_quality and conf < _hi: return "Hidden Gem"
        if tier in [CFG.TIER_POOR, CFG.TIER_DECOY] and conf >= _hi: return "Decoy"
        return "Ambiguous"

    df["Conflict_Category"] = df.apply(classify, axis=1)

    gems = df[df["Conflict_Category"] == "Hidden Gem"].sort_values("Pareto_Rank")
    gems.to_csv(_aux_dir(out_dir) / "04_ACTION_Rescue_Hidden_Gems.csv", index=False)

    reporter.log(f"Hidden Gems (Rescue Target) : {len(gems)}")
    reporter.log(f"Decoys (Potential Artifacts): {len(df[df['Conflict_Category'] == 'Decoy'])}")

    return df

# =============================================================================
# SECTION 4B: Tier_1A Landscape Companion Figures (the *b variants: 05b, 11b, 13b,
#             14b, 17b, 18b, 19b, 20b, 24b, 26b)
# =============================================================================

"""
Companion figures overlaying Tier_1A structural thumbnails onto the same
chemical spaces as their parent figures, showing where the Tier_1A hits
sit relative to the full 58 k-complex dataset.
"""

# ── Design constants ──────────────────────────────────────────────────────────
_TT_ENTRY_COLS = ["#E74C3C", "#17A589", "#27AE60", "#8E44AD", "#E67E22"]
_TT_STAR_FILL  = "#2ECC71"
_TT_STAR_S     = CFG.VIS_TT_STAR_SIZE
_TT_TIER_S     = CFG.VIS_TIER_SIZES
_TT_TIER_A     = CFG.VIS_TIER_ALPHAS
_TT_SHRINK_B   = CFG.VIS_TT_SHRINK_BORDER
_TT_IMG_PX     = CFG.VIS_TT_IMG_PX
_TT_TEXT_PX    = CFG.VIS_TT_TEXT_PX
_TT_BORDER_PX  = CFG.VIS_TT_BORDER_PX
_TT_PAD        = CFG.VIS_TT_PAD
_TT_FONT_SIZE  = CFG.VIS_TT_FONT_SIZE
# Resolve DejaVuSans portably: prefer the matplotlib-bundled copy (cross-platform),
# fall back to the common Linux path, else None (PIL default font is used downstream).
def _resolve_dejavu_font() -> str | None:
    try:
        from matplotlib import font_manager as _fm
        return _fm.findfont("DejaVu Sans", fallback_to_default=True)
    except Exception:
        _p = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
        return _p if _os.path.exists(_p) else None

_TT_FONT_PATH  = _resolve_dejavu_font()
_TT_AX_W       = CFG.VIS_TT_AX_WIDTH
_TT_LEFT_X0    = 0.058
_TT_RIGHT_X0   = 0.808
_TT_INSET_Y    = [0.60, 0.22, 0.69, 0.43, 0.17]
_TT_MARGINS    = dict(left=0.22, right=0.74, top=0.84, bottom=0.16)


def _tt_find_cif(pred_jobs: Path, job_name: str):
    jd = pred_jobs / job_name
    if not jd.exists():
        return None
    bc = jd / "Best_Complex"
    cifs = list(bc.glob("*.cif")) if bc.exists() else list(jd.rglob("*.cif"))
    return cifs[0] if cifs else None


def _tt_render_one(cif_path: Path, out_png: Path, width=800, height=800) -> bool:
    try:
        import pymol
        pymol.finish_launching(["pymol", "-cq"])
        from pymol import cmd
        cmd.reinitialize()
        cmd.feedback("disable", "all", "everything")
        cmd.load(str(cif_path), "mol")
        cmd.select("prot",   "mol and chain A")
        cmd.select("lig",    "mol and chain L")
        cmd.select("pocket", "byres (lig around 5) and chain A")
        cmd.hide("everything", "all")
        cmd.show("cartoon", "prot")
        cmd.color("grey70", "prot")
        cmd.show("surface", "prot")
        cmd.color("grey70", "prot")
        cmd.set("transparency",         0.30, "prot")
        cmd.set("cartoon_transparency", 0.15, "prot")
        cmd.show("sticks", "pocket")
        cmd.hide("sticks", "pocket and name N+CA+C+O")
        cmd.color("grey80", "pocket and elem C")
        cmd.set("stick_radius", 0.12, "pocket")
        cmd.show("sticks", "lig")
        cmd.color("magenta",  "lig and elem C")
        cmd.color("red",      "lig and elem O")
        cmd.color("cyan",     "lig and elem N")
        cmd.color("limegreen","lig and elem F")
        cmd.color("sulfur",   "lig and elem S")
        cmd.set("stick_radius", 0.20, "lig")
        cmd.distance("hbonds", "lig", "pocket", CFG.THRESHOLD_HB_DIST_MAX, mode=2)
        cmd.set("dash_color", "yellow", "hbonds")
        cmd.set("dash_width", 2.5)
        cmd.set("dash_gap",   0.30)
        cmd.hide("labels", "hbonds")
        if cmd.count_atoms("pocket and name CA") > 0:
            cmd.label("pocket and name CA", '"  %s%s" % (resn, resi)')
        cmd.set("label_color",   "white")
        cmd.set("label_size",    -0.45)
        cmd.set("label_font_id", 7)
        cmd.orient("pocket")
        cmd.zoom("pocket", buffer=CFG.VIS_TT_ZOOM_BUFFER)
        cmd.bg_color("black")
        cmd.viewport(width, height)
        cmd.set("ray_shadows",           "off")
        cmd.set("ambient",               0.6)
        cmd.set("ray_opaque_background", 1)
        cmd.png(str(out_png), width=width, height=height, ray=1, quiet=1)
        return out_png.exists()
    except Exception:
        return False


def _tt_render_all(pa: pd.DataFrame, pred_jobs: Path, thumb_dir: Path, reporter) -> list:
    thumb_dir.mkdir(exist_ok=True)
    imgs = []
    # Render at most CFG.VIS_MAX_THUMBNAILS thumbnails (top-N rows of the
    # rank-sorted Pareto frame) so a populous Tier_1A cannot flood the panel.
    for local_i, (_, row) in enumerate(pa.head(CFG.VIS_MAX_THUMBNAILS).iterrows()):
        png = thumb_dir / f"pa_{local_i}.png"
        if not png.exists():
            cif = _tt_find_cif(pred_jobs, str(row.get("job_name", "")))
            if cif:
                _tt_render_one(cif, png)
        imgs.append(np.array(Image.open(png).convert("RGBA")) if png.exists() else None)
    return imgs


def _tt_make_composite(img_arr, label_str: str, colour_hex: str) -> np.ndarray:
    PAD = _TT_PAD
    br, bg, bb = (int(colour_hex[i:i+2], 16) for i in (1, 3, 5))
    W, H_s = _TT_IMG_PX, _TT_TEXT_PX
    src = Image.fromarray(img_arr).convert("RGBA")
    inner = _TT_IMG_PX - 2 * PAD
    src_small = src.resize((inner, inner), Image.LANCZOS)
    canvas = Image.new("RGBA", (W, W + H_s), (15, 15, 25, 255))
    canvas.paste(src_small, (PAD, PAD))
    draw = ImageDraw.Draw(canvas)
    draw.line([(  _TT_BORDER_PX, W), (W - _TT_BORDER_PX, W)],
              fill=(br, bg, bb, 200), width=4)
    usable_w = W - 2 * _TT_BORDER_PX - 8
    lines = label_str.split("\n")
    font = ImageFont.load_default()
    for fs in range(_TT_FONT_SIZE, 14, -2):
        try:
            f = ImageFont.truetype(_TT_FONT_PATH, fs)
        except Exception:
            break
        if max(draw.textbbox((0, 0), l, font=f)[2] for l in lines) <= usable_w:
            font = f; break
    y_cur = W + 18
    for line in lines:
        bb_box = draw.textbbox((0, 0), line, font=font)
        tw, lh = bb_box[2] - bb_box[0], bb_box[3] - bb_box[1]
        draw.text(((W - tw) // 2, y_cur), line, font=font, fill=(220, 220, 230, 255))
        y_cur += lh + 10
    draw.rectangle([_TT_BORDER_PX // 2, _TT_BORDER_PX // 2,
                    W - _TT_BORDER_PX // 2 - 1, W + H_s - _TT_BORDER_PX // 2 - 1],
                   outline=(br, bg, bb, 255), width=_TT_BORDER_PX)
    return np.array(canvas)


def _tt_draw_scatter(ax, df):
    for tier in [t for t in reversed(CFG.TIER_ORDER) if t != CFG.TIER_TOP]:
        sub = df[df["degrader_tier"] == tier]
        if sub.empty: continue
        ax.scatter(sub["_X"], sub["_Y"],
                   c=TIER_PALETTE.get(tier, "#BBB"),
                   s=_TT_TIER_S.get(tier, 4), alpha=_TT_TIER_A.get(tier, 0.2),
                   linewidths=0, rasterized=True, zorder=2)


def _tt_draw_stars(ax, pa):
    for local_i, (_, row) in enumerate(pa.iterrows()):
        ec   = _TT_ENTRY_COLS[local_i % len(_TT_ENTRY_COLS)]
        fill = "#F1C40F" if local_i == 4 else _TT_STAR_FILL
        sz   = int(_TT_STAR_S * 0.82) if local_i == 4 else _TT_STAR_S
        ax.scatter(row["_X"], row["_Y"], c=fill, s=sz, marker="*",
                   edgecolors=ec, linewidths=1.2, zorder=10)


def _tt_draw_thumbnails(fig, ax, pa, imgs):
    fw, fh = fig.get_size_inches()
    # Show the top CFG.VIS_MAX_THUMBNAILS Tier_1A entries — distribute evenly
    # between left and right columns.
    n_max = min(len(pa), len(imgs), CFG.VIS_MAX_THUMBNAILS)
    _ax_w = 0.085 if n_max > 5 else _TT_AX_W
    ax_h = _ax_w * (fw / fh) * (_TT_IMG_PX + _TT_TEXT_PX) / _TT_IMG_PX
    """
    Y layout: split into left (even indices) and right (odd indices) columns. Thumbnails
    are bottom-anchored and stacked with a fixed gap so boxes never touch; if the natural
    height would overflow the column, the whole thumbnail is shrunk (aspect preserved).
    """
    _y_margin, _y_top = 0.08, 0.90
    _y_gap   = 0.022                                  # minimum vertical gap between stacked thumbnails
    _span    = _y_top - _y_margin
    # Right column fills first, capped at 5 per side. A cohort of ≤5 sits entirely on the
    # right (the empty left margin is then cropped by bbox_inches="tight", narrowing the
    # figure); a cohort >5 is split evenly across both columns (6 → 3+3, 7 → 4+3, 10 → 5+5).
    if n_max <= 5:
        n_right, n_left = n_max, 0
    else:
        n_right = (n_max + 1) // 2
        n_left  = n_max - n_right
    _n_col   = max(n_left, n_right, 1)
    _max_h   = (_span - (_n_col - 1) * _y_gap) / _n_col
    if ax_h > _max_h:                                 # shrink to fit (keep aspect) → no overlap
        _ax_w *= _max_h / ax_h
        ax_h = _max_h
    def _ycols(n):
        if n == 0:
            return []
        _block = n * ax_h + (n - 1) * _y_gap
        _start = _y_margin + max(0.0, (_span - _block) / 2.0)   # centre the column vertically
        return [_start + i * (ax_h + _y_gap) for i in range(n)]
    _y_left  = _ycols(n_left)
    _y_right = _ycols(n_right)
    for local_i, (_, row) in enumerate(pa.head(n_max).iterrows()):
        colour = _TT_ENTRY_COLS[local_i % len(_TT_ENTRY_COLS)]
        if local_i < n_right:                        # right column fills first
            side, col_idx, x0 = "right", local_i, _TT_RIGHT_X0
            y0 = _y_right[col_idx]
        else:                                        # overflow to the left column
            side, col_idx, x0 = "left", local_i - n_right, _TT_LEFT_X0
            y0 = _y_left[col_idx]
        sx, sy = row["_X"], row["_Y"]
        img    = imgs[local_i]
        if img is None: continue
        prot = str(row.get("Protein_Name", "")).replace("_Control","").split("_")[0][:11]
        lig  = str(row.get("Ligand_Name", ""))
        lig  = lig[:13] + "…" if len(lig) > 13 else lig
        sn2  = row.get("SN2_Attack_Angle",        float("nan"))
        conf = row.get("Boltz_Model_Confidence",   float("nan"))
        ipt  = row.get("iptm",                     float("nan"))
        idn  = row.get("Interaction_Density_Norm", float("nan"))
        label = (f"{prot}  ·  {lig}\n"
                 f"SN2: {sn2:.1f}°   Conf: {conf:.3f}\n"
                 f"ipTM: {ipt:.3f}  IDens: {idn:.2f}")
        comp  = _tt_make_composite(img, label, colour)
        ax_in = fig.add_axes([x0, y0, _ax_w, ax_h])
        h_c, w_c = comp.shape[:2]
        ax_in.imshow(comp, aspect="auto", interpolation="lanczos",
                     extent=[0, w_c, h_c, 0])
        ax_in.set_xlim(0, w_c); ax_in.set_ylim(h_c, 0)
        ax_in.set_axis_off()
        axy = y0 + ax_h / 2
        axx = (x0 + _ax_w) if side == "left" else x0
        fig.add_artist(ConnectionPatch(
            xyA=(axx, axy), coordsA="figure fraction",
            xyB=(sx, sy),   coordsB="data",
            axesA=None, axesB=ax,
            arrowstyle="-|>", color=colour, lw=2.5,
            mutation_scale=22, shrinkA=4, shrinkB=_TT_SHRINK_B, zorder=11,
        ))


def _tt_legend_handles(df):
    n_tt = len(df[df["degrader_tier"] == CFG.TIER_TOP])
    # Symbol-type header entries so the reader understands both glyphs
    h = [
        Line2D([0],[0], marker="o", color="w",
               markerfacecolor="#888888", markeredgecolor="none",
               markersize=7, alpha=0.55, label="Individual complex  (scatter dot)"),
        Line2D([0],[0], marker="*", color="w",
               markerfacecolor=_TT_STAR_FILL, markeredgecolor=_TT_ENTRY_COLS[0],
               markersize=14, label=f"{CFG.TIER_TOP} ★ highlighted  (n={n_tt})"),
        Line2D([0],[0], marker="*", color="w",
               markerfacecolor="#F1C40F", markeredgecolor=_TT_ENTRY_COLS[4],
               markersize=11, label=f"{CFG.TIER_TOP} ★ gold (TFA duplicate)"),
    ]
    for t in [t for t in CFG.TIER_ORDER if t != CFG.TIER_TOP]:
        sub = df[df["degrader_tier"] == t]
        if sub.empty:
            continue
        h.append(Line2D([0],[0], marker="o", color="w",
                        markerfacecolor=TIER_PALETTE.get(t, "#BBB"),
                        markersize=7, alpha=0.85, label=f"{t}  (n={len(sub):,})"))
    return h


def _tt_add_legend(fig, handles):
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.02), ncol=4,
               fontsize=8.5, framealpha=0.95,
               title=f"Degrader tier  (★ = {CFG.TIER_TOP} highlighted; ● = scatter background)",
               title_fontsize=9.0,
               borderpad=0.8, labelspacing=0.6, handletextpad=0.5)


def _tt_new_fig():
    fig, ax = plt.subplots(figsize=(18, 12))
    fig.subplots_adjust(**_TT_MARGINS)
    ax.set_facecolor("#FAFAFA")
    return fig, ax


def _tt_style(ax, xlabel, ylabel, title=None):
    # `title` accepted for call-site compatibility but intentionally not rendered
    # (figures carry no titles).
    ax.set_xlabel(xlabel, fontsize=11, labelpad=6)
    ax.set_ylabel(ylabel, fontsize=11, labelpad=6)
    ax.tick_params(labelsize=9)
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.18, color="#AAAAAA", linewidth=0.5)


# ── Figure 18b — UMAP Chemical Space + PA Landscape ──────────────────────────
def _fig_18b_tt_landscape(df, pa, imgs, out_dir: Path, reporter):
    try:
        dv = df.dropna(subset=["UMAP_X","UMAP_Y"]).copy()
        dv["_X"] = dv["UMAP_X"]; dv["_Y"] = dv["UMAP_Y"]
        pax = pa.copy(); pax["_X"] = pa["UMAP_X"]; pax["_Y"] = pa["UMAP_Y"]
        xy  = np.vstack([dv["_X"].values, dv["_Y"].values])
        kde = gaussian_kde(xy, bw_method=0.18)
        x0,x1 = dv["_X"].min()-1, dv["_X"].max()+1
        y0,y1 = dv["_Y"].min()-1, dv["_Y"].max()+1
        xi,yi = np.linspace(x0,x1,220), np.linspace(y0,y1,220)
        Xi,Yi = np.meshgrid(xi,yi)
        Zi = kde(np.vstack([Xi.ravel(),Yi.ravel()])).reshape(Xi.shape)
        # The density floor (1e-6) caps the relative free energy at
        # F_max = -0.596·ln(1e-6) ≈ 8.23 kcal/mol; tie vmax to it so the full
        # colourbar range is populated instead of leaving the upper ~45% empty.
        _dens_floor = 1e-6
        _F_MAX = float(-0.596 * np.log(_dens_floor))
        F  = np.clip(-0.596*np.log(np.clip(Zi/Zi.max(), _dens_floor, None)), 0, _F_MAX)
        fig, ax = _tt_new_fig()
        cf = ax.contourf(Xi,Yi,F, levels=45, cmap="Blues_r", vmin=0, vmax=_F_MAX, alpha=0.90)
        ax.contour(Xi,Yi,F, levels=18, colors="white", linewidths=0.30, alpha=0.38)
        cb = fig.colorbar(cf, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("KDE Population Density  (darker = more complexes)", fontsize=9.5)
        cb.set_ticks([0, 2, 4, 6, 8])
        _tt_draw_scatter(ax, dv); _tt_draw_stars(ax, pax)
        _tt_draw_thumbnails(fig, ax, pax, imgs)
        _tt_style(ax, "UMAP 1", "UMAP 2", "")
        _tt_add_legend(fig, _tt_legend_handles(dv))
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        out = out_dir / "Figure_18b_TT_Chemical_Space_Landscape.png"
        fig.savefig(out, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight", facecolor="white")
        plt.close(fig)
    except Exception as e:
        reporter.log(f"  ! Figure 18b skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig


# ── Figure 13b — SN2 Angle × Confidence + PA Landscape ───────────────────────
def _fig_13b_tt_mechanistic(df, pa, imgs, out_dir: Path, reporter):
    try:
        dv = df.dropna(subset=["SN2_Attack_Angle","Boltz_Model_Confidence"]).copy()
        dv["_X"] = pd.to_numeric(dv["SN2_Attack_Angle"],      errors="coerce")
        dv["_Y"] = pd.to_numeric(dv["Boltz_Model_Confidence"], errors="coerce")
        dv = dv.dropna(subset=["_X","_Y"])
        pax = pa.copy()
        pax["_X"] = pd.to_numeric(pa["SN2_Attack_Angle"],      errors="coerce")
        pax["_Y"] = pd.to_numeric(pa["Boltz_Model_Confidence"], errors="coerce")
        ylo = dv["_Y"].quantile(0.01) - 0.005
        yhi = dv["_Y"].max() + 0.005
        fig, ax = _tt_new_fig()
        hb = ax.hexbin(dv["_X"], dv["_Y"], gridsize=60, cmap="YlOrBr",
                       mincnt=1, linewidths=0.15, alpha=0.85, zorder=1,
                       extent=[dv["_X"].min(), dv["_X"].max(), ylo, yhi])
        cb = fig.colorbar(hb, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("Complex count per bin", fontsize=9.5)
        for tlbl, tcol in [(CFG.TIER_TOP, _TT_STAR_FILL),
                           (CFG.TIER_ORDER[1], TIER_PALETTE[CFG.TIER_ORDER[1]]),
                           (CFG.TIER_ORDER[2], TIER_PALETTE[CFG.TIER_ORDER[2]])]:
            angle = CFG.TIER_ANGLE_MIN[tlbl]   # single source: tier angle threshold from CFG
            ax.axvline(angle, color=tcol, ls="--", lw=1.6, alpha=0.85, zorder=4)
            ax.text(angle+0.3, yhi-0.002, f"≥{angle}°\n{tlbl}",
                    fontsize=5.5, color=tcol, fontweight="bold", va="top")
        ax.axhline(dv["_Y"].median(), color="#555", ls=":", lw=0.9, alpha=0.55, zorder=3)
        _tt_draw_scatter(ax, dv); _tt_draw_stars(ax, pax)
        _tt_draw_thumbnails(fig, ax, pax, imgs)
        _tt_style(ax,
                  "SN2 Attack Angle (°) — higher = near-ideal nucleophilic trajectory",
                  "Boltz Model Confidence (AI structural quality, higher = better)",
                  "")
        _tt_add_legend(fig, _tt_legend_handles(dv))
        ax.set_ylim(ylo, yhi)
        out = out_dir / "Figure_13b_TT_Mechanistic_Quality_Space.png"
        fig.savefig(out, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight", facecolor="white")
        plt.close(fig)
    except Exception as e:
        reporter.log(f"  ! Figure 13b skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig


# ── Figure 05b — Confidence × ipTM + PA Landscape ────────────────────────────
def _fig_05b_tt_ai_quality(df, pa, imgs, out_dir: Path, reporter):
    try:
        dv = df.dropna(subset=["Boltz_Model_Confidence","iptm"]).copy()
        dv["_X"] = dv["Boltz_Model_Confidence"]; dv["_Y"] = dv["iptm"]
        pax = pa.copy()
        pax["_X"] = pa["Boltz_Model_Confidence"]; pax["_Y"] = pa["iptm"]
        xlo,xhi = dv["_X"].quantile(0.01)-0.01, dv["_X"].max()+0.005
        ylo,yhi = dv["_Y"].quantile(0.01)-0.01, dv["_Y"].max()+0.005
        fig, ax = _tt_new_fig()
        hb = ax.hexbin(dv["_X"], dv["_Y"], gridsize=60, cmap="YlOrBr",
                       mincnt=1, linewidths=0.15, alpha=0.85, zorder=1,
                       extent=[xlo,xhi,ylo,yhi])
        cb = fig.colorbar(hb, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("Complex count per bin", fontsize=9.5)
        for v, fn in [(dv["_X"].median(), ax.axvline),(dv["_Y"].median(), ax.axhline)]:
            fn(v, color="#555", ls="--", lw=0.9, alpha=0.55, zorder=3)
        _tt_draw_scatter(ax, dv); _tt_draw_stars(ax, pax)
        _tt_draw_thumbnails(fig, ax, pax, imgs)
        _tt_style(ax,
                  "Boltz Model Confidence (higher = better folding quality)",
                  "ipTM — Interface Predicted TM-score",
                  "")
        _tt_add_legend(fig, _tt_legend_handles(dv))
        ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
        out = out_dir / "Figure_05b_TT_AI_Quality_Space.png"
        fig.savefig(out, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight", facecolor="white")
        plt.close(fig)
    except Exception as e:
        reporter.log(f"  ! Figure 05b skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig


# ── Figure 14b — Interaction Density × Count + PA Landscape ──────────────────
def _fig_14b_tt_interactions(df, pa, imgs, out_dir: Path, reporter):
    try:
        dv = df.dropna(subset=["Interaction_Density_Norm","num_interactions"]).copy()
        dv["_X"] = pd.to_numeric(dv["Interaction_Density_Norm"], errors="coerce")
        dv["_Y"] = pd.to_numeric(dv["num_interactions"],          errors="coerce")
        dv = dv.dropna(subset=["_X","_Y"])
        pax = pa.copy()
        pax["_X"] = pd.to_numeric(pa["Interaction_Density_Norm"], errors="coerce")
        pax["_Y"] = pd.to_numeric(pa["num_interactions"],          errors="coerce")
        # np.nanmax guards an empty Pareto frame (pax[...].max() → NaN → NaN axis limits → set_xlim failure).
        xhi = float(np.nanmax([dv["_X"].quantile(0.99)+0.5, (pax["_X"].max() if len(pax) else np.nan)+0.5]))
        yhi = float(np.nanmax([dv["_Y"].quantile(0.99)+2,   (pax["_Y"].max() if len(pax) else np.nan)+2]))
        dv  = dv[(dv["_X"] <= xhi) & (dv["_Y"] <= yhi)]
        fig, ax = _tt_new_fig()
        hb = ax.hexbin(dv["_X"], dv["_Y"], gridsize=60, cmap="YlOrBr",
                       mincnt=1, linewidths=0.15, alpha=0.85, zorder=1)
        cb = fig.colorbar(hb, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("Complex count per bin", fontsize=9.5)
        xmed, ymed = dv["_X"].median(), dv["_Y"].median()
        ax.axvline(xmed, color="#555", ls="--", lw=0.9, alpha=0.55, zorder=3)
        ax.axhline(ymed, color="#555", ls="--", lw=0.9, alpha=0.55, zorder=3)
        ax.text(xmed+0.05, 0.5, f"median {xmed:.2f}", fontsize=7.5, color="#555")
        ax.text(0.1, ymed+0.3, f"median {ymed:.0f}",  fontsize=7.5, color="#555")
        _tt_draw_scatter(ax, dv); _tt_draw_stars(ax, pax)
        _tt_draw_thumbnails(fig, ax, pax, imgs)
        _tt_style(ax,
                  "Interaction Density Norm (interactions per ligand heavy atom)",
                  "Total Interactions (num_interactions, all contact types)",
                  "")
        _tt_add_legend(fig, _tt_legend_handles(dv))
        ax.set_xlim(0, xhi); ax.set_ylim(0, yhi)
        out = out_dir / "Figure_14b_TT_Interaction_Quality_Space.png"
        fig.savefig(out, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight", facecolor="white")
        plt.close(fig)
    except Exception as e:
        reporter.log(f"  ! Figure 14b skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig


_FIG_MAPPING = {
        # ── 02_Dataset_and_Alignment_Overview ──
        "Figure_01_Active_Site_Residue_Mapping_Coverage.png": "02_Dataset_and_Alignment_Overview/01_Active_Site_Residue_Mapping_Coverage.png",
        "Figure_02_Tier_Distribution.png": "02_Dataset_and_Alignment_Overview/02_Tier_Distribution.png",
        "Figure_03_Alignment_Grades.png": "02_Dataset_and_Alignment_Overview/03_Alignment_Grades.png",
        "Figure_04_Tier_Grade_Distribution.png": "02_Dataset_and_Alignment_Overview/04_Tier_Grade_Distribution.png",
        # ── 03_AI_Confidence_Quality ──
        "Figure_05a_AI_Quality_Assessment.png": "03_AI_Confidence_Quality/01_AI_Quality_Assessment.png",
        "Figure_05b_TT_AI_Quality_Space.png": "03_AI_Confidence_Quality/02_AI_Quality_Space.png",
        "Figure_06_pTM_vs_ipTM_by_Tier.png": "03_AI_Confidence_Quality/03_pTM_vs_ipTM_by_Tier.png",
        # ── 04_Catalytic_Geometry_and_Mechanism ──
        "Figure_07_ActiveSite_RMSD_by_Tier.png": "04_Catalytic_Geometry_and_Mechanism/01_ActiveSite_RMSD_by_Tier.png",
        "Figure_11b_Mechanistic_Fingerprint.png": "04_Catalytic_Geometry_and_Mechanism/09_Mechanistic_Fingerprint.png",
        "Figure_08_Feature_Correlations.png": "04_Catalytic_Geometry_and_Mechanism/07_Feature_Correlations.png",
        "Figure_09_Tier_Quality_DotPlot.png": "04_Catalytic_Geometry_and_Mechanism/08_Tier_Quality_DotPlot.png",
        "Figure_10_Mech_State_CrossTab.png": "04_Catalytic_Geometry_and_Mechanism/02_Mech_State_CrossTab.png",
        "Figure_11_Mechanistic_Score_by_Tier.png": "04_Catalytic_Geometry_and_Mechanism/03_Mechanistic_Score_by_Tier.png",
        "Figure_12_SN2_Angle_by_Tier.png": "04_Catalytic_Geometry_and_Mechanism/04_SN2_Angle_by_Tier.png",
        "Figure_13a_Mechanism_Geometry_Scatter.png": "04_Catalytic_Geometry_and_Mechanism/05_Mechanism_Geometry_Scatter.png",
        "Figure_13b_TT_Mechanistic_Quality_Space.png": "04_Catalytic_Geometry_and_Mechanism/06_Mechanistic_Quality_Space.png",
        # ── 05_Ligand_Interactions_and_Chemical_Space ──
        "Figure_14a_Molecular_Interaction_Profile.png": "05_Ligand_Interactions_and_Chemical_Space/01_Molecular_Interaction_Profile.png",
        "Figure_17b_Binding_Energetics.png": "05_Ligand_Interactions_and_Chemical_Space/08_Binding_Energetics.png",
        "Figure_14b_TT_Interaction_Quality_Space.png": "05_Ligand_Interactions_and_Chemical_Space/02_Interaction_Quality_Space.png",
        "Figure_15_Fluorine_Engagement_by_Tier.png": "05_Ligand_Interactions_and_Chemical_Space/03_Fluorine_Engagement_by_Tier.png",
        "Figure_16_Catalytic_Quality_vs_Inhibition.png": "05_Ligand_Interactions_and_Chemical_Space/04_Catalytic_Quality_vs_Inhibition.png",
        "Figure_17_ActiveSite_Contact_Density_by_Tier.png": "05_Ligand_Interactions_and_Chemical_Space/05_ActiveSite_Contact_Density_by_Tier.png",
        "Figure_18a_Chemical_Space_Map.png": "05_Ligand_Interactions_and_Chemical_Space/06_Chemical_Space_Map.png",
        "Figure_18b_TT_Chemical_Space_Landscape.png": "05_Ligand_Interactions_and_Chemical_Space/07_Chemical_Space_Landscape.png",
        # ── 06_PFAS_Scope_and_Synthesis ──
        "Figure_19a_Radar_TopHits.png": "06_PFAS_Scope_and_Synthesis/01_Radar_TopHits.png",
        "Figure_19b_Radar_TierReps.png": "06_PFAS_Scope_and_Synthesis/02_Radar_TierReps.png",
        "Figure_20a_Tier_Success_Rates.png": "06_PFAS_Scope_and_Synthesis/03_Tier_Success_Rates.png",
        "Figure_20b_Conf_SN2_Landscape.png": "06_PFAS_Scope_and_Synthesis/04_Conf_SN2_Landscape.png",
        "Figure_21_Conflict_Composition.png": "06_PFAS_Scope_and_Synthesis/05_Conflict_Composition.png",
        "Figure_22_Hidden_Gems_DeepDive.png": "06_PFAS_Scope_and_Synthesis/06_Hidden_Gems_DeepDive.png",
        "Figure_23_Category_Overlap_Euler.png": "06_PFAS_Scope_and_Synthesis/07_Category_Overlap_Euler.png",
        "Figure_24a_Top25_Multitarget_Proteins.png": "06_PFAS_Scope_and_Synthesis/08_Top25_Multitarget_Proteins.png",
        "Figure_24b_TopTier_Protein_PFAS_Breakdown.png": "06_PFAS_Scope_and_Synthesis/09_TopTier_Protein_PFAS_Breakdown.png",
        "Figure_25_Sankey_Workflow.png": "06_PFAS_Scope_and_Synthesis/10_Sankey_Workflow.png",
        "Figure_26a_PFAS_Size_Hexbin_Landscape.png": "06_PFAS_Scope_and_Synthesis/11_PFAS_Size_Hexbin_Landscape.png",
        "Figure_26b_PFAS_Size_Composition_Merged.png": "06_PFAS_Scope_and_Synthesis/12_PFAS_Size_Composition_Merged.png",
        "Figure_26c_PFAS_Carbon_Confidence.png": "06_PFAS_Scope_and_Synthesis/13_PFAS_Carbon_Confidence.png",
        # ── Two-criteria tier logic + dead-end feasibility (single 3-panel figure,
        #    filed in the catalytic folder) ──
        "Figure_30_TwoCriteria_Tier_Logic.png": "04_Catalytic_Geometry_and_Mechanism/10_TwoCriteria_Tier_Logic.png",
}


@contextlib.contextmanager
def _redirect_savefig(out_dir: Path, reporter):
    """Route each Figure_NN_*.png save to its numbered folder (per _FIG_MAPPING)
    for the duration of the block, restoring the original savefig on exit — even
    on exception, so a failure mid-suite can never leak the patched savefig into
    later routines or the rest of the process. Target folders are created on entry.
    """
    import matplotlib
    for _rel in _FIG_MAPPING.values():
        (out_dir / _rel).parent.mkdir(parents=True, exist_ok=True)
    _orig_plt = plt.savefig
    _orig_fig = matplotlib.figure.Figure.savefig
    def _redirect(fname):
        name = fname.name if isinstance(fname, Path) else Path(fname).name
        if name in _FIG_MAPPING:
            reporter.log(f'  \u2714 Saved: {_FIG_MAPPING[name]}')
            return out_dir / _FIG_MAPPING[name]
        return fname
    plt.savefig = lambda fname, *a, **k: _orig_plt(_redirect(fname), *a, **k)
    matplotlib.figure.Figure.savefig = (
        lambda self, fname, *a, **k: _orig_fig(self, _redirect(fname), *a, **k))
    try:
        yield
    finally:
        plt.savefig = _orig_plt
        matplotlib.figure.Figure.savefig = _orig_fig

def generate_comprehensive_figures(df: pd.DataFrame, features: list[str], out_dir: Path, reporter: ReportManager):
    """Render the full figure suite. The _redirect_savefig context manager routes
    plt.savefig / Figure.savefig to the numbered folders and always restores them
    on exit, so an exception mid-suite can never leak the patched savefig into
    later routines or the rest of the process.
    """
    with _redirect_savefig(out_dir, reporter):
        return _generate_comprehensive_figures_impl(df, features, out_dir, reporter)


def _generate_comprehensive_figures_impl(df: pd.DataFrame, features: list[str], out_dir: Path, reporter: ReportManager):
    """Render the full publication figure suite (folders 02–07).

    Builds the dataset/confidence/geometry/interaction/scope figures (02–06), the
    diagnostic trends (07) and the two-criteria active-site/feasibility figures
    (08). Each "Figure_NN_*.png" save is routed to its numbered folder by the
    _redirect_savefig context manager wrapping this call (see
    generate_comprehensive_figures). Side-effecting (writes PNGs); returns None.
    """
    reporter.section("Step 4/5 — Main Validation Figure Suite (Publication Quality)  [writes folders 02–06]")

    existing_tiers = [t for t in TIER_ORDER_LOGIC if t in df["degrader_tier"].unique()]

    # ── Global publication-quality rcParams ───────────────────────────────────
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.titlesize": 13, "axes.titleweight": "bold",
        "axes.labelsize": 11,
        "xtick.labelsize": 9, "ytick.labelsize": 9,
        "legend.fontsize": 8.5, "legend.framealpha": 0.92,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": "#EBEBEB", "grid.linewidth": 0.6,
    })

    def _text_color(hex_bg: str, threshold: float = 0.5) -> str:
        """Auto-contrast label colour for hex_bg (delegates to _auto_label_colour)."""
        return _auto_label_colour(hex_bg, threshold)

    reporter.section("  Folder 02_Dataset_and_Alignment_Overview — coverage, tiers, alignment grades")
    # --- Figure 01: Active-site residue mapping coverage across all variants ---
    """
    For each of the 8 canonical FAcD catalytic residues, the number of protein variants
    in which the residue was resolved during the positional alignment (Alignment_Stats.csv,
    written per protein by 02_Production). A residue counts as mapped when its key in
    active_site_mapping resolves to a target index. A companion CSV lists, per residue,
    the sequence IDs that failed to map it.
    """
    import json as _json01
    from collections import Counter as _Counter01
    _aln_csv = out_dir.parent / "1_Boltz2_Production" / "3_Sequence_Reference_Data" / "Active_Site_Alignments" / CFG.FILE_ALIGNMENT_STATS
    if _aln_csv.exists():
        _aa3 = {"ASP": "Asp", "ARG": "Arg", "HIS": "His", "TRP": "Trp", "TYR": "Tyr"}
        _ref = CFG.REF_ACTIVE_SITE_MAP
        # Order residues by mechanistic role (CFG.ACTIVE_SITE_ROLE_ORDER).
        _keys = [k for k in CFG.ACTIVE_SITE_ROLE_ORDER if k in _ref]
        _labels01 = [f"{_aa3.get(_ref[k]['res'], _ref[k]['res'].title())}{_ref[k]['id']}\n({k})" for k in _keys]

        _adf = pd.read_csv(_aln_csv)
        if "protein" in _adf.columns:
            # Exclude the control jobs (their protein names end in "_Control") so the
            # coverage count and the companion missed-CSV reflect only the query set.
            _adf = _adf[~_adf["protein"].astype(str).str.endswith("_Control")]
            _adf = _adf.drop_duplicates("protein")
        _n01 = len(_adf)
        _cnt01 = _Counter01()
        _missed01 = {k: [] for k in _keys}
        _all_mapped01 = 0
        for _prot, _s in zip(_adf.get("protein", range(_n01)), _adf["active_site_mapping"].fillna("")):
            try:
                _m01 = _json01.loads(_s)
            except Exception:
                continue
            _found01 = 0
            for k in _keys:
                if _m01.get(k) is not None:
                    _cnt01[k] += 1; _found01 += 1
                else:
                    _missed01[k].append(str(_prot))
            if _found01 == len(_keys):
                _all_mapped01 += 1

        _counts01 = [_cnt01[k] for k in _keys]
        _pcts01 = [100.0 * c / _n01 if _n01 else 0.0 for c in _counts01]

        _col_labels01 = {k: f"{_aa3.get(_ref[k]['res'], _ref[k]['res'].title())}{_ref[k]['id']}_{k}" for k in _keys}
        _maxlen01 = max((len(v) for v in _missed01.values()), default=0)
        pd.DataFrame({
            _col_labels01[k]: _missed01[k] + [""] * (_maxlen01 - len(_missed01[k])) for k in _keys
        }).to_csv(_aux_dir(out_dir) / "01_Active_Site_Residue_Mapping_Missed.csv", index=False)

        fig, ax = plt.subplots(figsize=(11, 6))
        # Colour each bar by its functional role group (CFG single source) so residues
        # of the same mechanistic role share a colour.
        _groups01 = [CFG.ACTIVE_SITE_ROLE_GROUP.get(k, "") for k in _keys]
        _colours01 = [CFG.ACTIVE_SITE_ROLE_GROUP_COLOUR.get(g, "#999999") for g in _groups01]
        ax.bar(range(len(_keys)), _counts01, color=_colours01, edgecolor="white", linewidth=0.8, zorder=2)
        ax.axhline(_n01, color="grey", linestyle="--", linewidth=1, zorder=1)
        '''
        Data labels INSIDE each bar: vertical (90°) "count | pct%" anchored at the
        bar's middle-top and extending downward. Each datum — count, "|" separator and
        percentage — is drawn as its own segment in a distinct CFG hue; the triple is
        selected by bar-fill luminance (light segments on dark fills, dark on light) and
        carries a contrast halo so every segment stays legible against any role colour.
        '''
        _lbl_off01 = (_n01 * 0.01) if _n01 else 0.01
        fig.canvas.draw()                       # realise the layout so text extents are valid
        _renderer01 = fig.canvas.get_renderer()
        # Pixels-per-data-unit on the y-axis, for converting rendered text extents to data
        # heights (vertical labels stack along y).
        _ppd01 = (ax.transData.transform((0, 1.0))[1]
                  - ax.transData.transform((0, 0.0))[1]) or 1.0
        _gap01 = 5.0 / _ppd01                    # 5-px inter-segment gap, in data units
        for i, (c, p) in enumerate(zip(_counts01, _pcts01)):
            _on_dark01 = _text_color(_colours01[i]) == "white"
            _seg_cols01 = (CFG.VIS_BAR_LABEL_COLOURS_ON_DARK if _on_dark01
                           else CFG.VIS_BAR_LABEL_COLOURS_ON_LIGHT)
            _halo01 = "black" if _on_dark01 else "white"
            _segs01 = [f"{c:,}", " | ", f"{p:.1f}%"]
            # Measure each segment's rendered length (data units) so the segments can be
            # centred end-to-end with no overlap; the 90°-rotated label reads bottom-to-top.
            _hts01 = []
            for _s01 in _segs01:
                _probe01 = ax.text(0, 0, _s01, rotation=90, ha="center", va="center",
                                   fontsize=8, fontweight="bold")
                _probe01.draw(_renderer01)
                _hts01.append(_probe01.get_window_extent(_renderer01).height / _ppd01)
                _probe01.remove()
            _stack01 = sum(_hts01) + _gap01 * (len(_hts01) - 1)
            _acc01 = max(c - _lbl_off01 - _stack01, _lbl_off01)
            for _s01, _sc01, _h01 in zip(_segs01, _seg_cols01, _hts01):
                ax.text(i, _acc01 + _h01 / 2.0, _s01,
                        rotation=90, ha="center", va="center",
                        fontsize=8, fontweight="bold", color=_sc01, zorder=6,
                        path_effects=[pe.withStroke(linewidth=1.4, foreground=_halo01)])
                _acc01 += _h01 + _gap01
        ax.set_xticks(range(len(_keys)))
        ax.set_xticklabels(_labels01, fontsize=9)
        ax.set_ylabel("Variants with residue mapped")
        ax.set_ylim(0, _n01 * 1.08 if _n01 else 1)
        # Legend: one entry per role group (in mechanistic order) + the all-variants line,
        # laid out as a single row in the top-left.
        from matplotlib.patches import Patch as _Patch01
        _seen01 = []
        for g in _groups01:
            if g not in _seen01:
                _seen01.append(g)
        _handles01 = [_Patch01(facecolor=CFG.ACTIVE_SITE_ROLE_GROUP_COLOUR.get(g, "#999999"),
                               edgecolor="white", label=g) for g in _seen01]
        _handles01.append(Line2D([0], [0], color="grey", linestyle="--", label=f"All variants (n={_n01:,})"))
        ax.legend(handles=_handles01, loc="upper left", ncol=len(_handles01),
                  fontsize=8, frameon=True, columnspacing=1.0, handletextpad=0.5)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_01_Active_Site_Residue_Mapping_Coverage.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()
    else:
        reporter.log(f"  ! Skipped Figure 01: Alignment_Stats.csv not found at {_aln_csv}")

    # --- Figure 02: Tier Distribution + Model Selection (pie inset) ---
    if "degrader_tier" in df.columns:
        total_complexes = len(df)
        model_col = "best_model_name" if "best_model_name" in df.columns else None
        counts6 = [len(df[df["degrader_tier"] == t]) for t in existing_tiers]

        # Adaptive figure height: scale with the tallest bar so the pie doesn't float
        _max_cnt6 = max(counts6) if counts6 else 1
        _fig_h6 = max(5.5, min(10.0, 4.5 + _max_cnt6 / 12000))

        fig, ax = plt.subplots(figsize=(11, _fig_h6))
        bar_x6 = np.arange(len(existing_tiers))

        ax.bar(bar_x6, counts6,
               color=[TIER_PALETTE.get(t, "#999") for t in existing_tiers],
               edgecolor="white", linewidth=0.8, width=0.65, zorder=2)
        for xi, cnt in zip(bar_x6, counts6):
            pct = cnt / total_complexes * 100 if total_complexes else 0
            _pct_str6 = (f"{pct:.3f}" if pct < 0.1 else
                         f"{pct:.2f}" if pct < 1.0 else
                         f"{pct:.1f}")
            ax.text(xi, cnt + _max_cnt6 * 0.004,
                    f"{cnt:,}\n({_pct_str6}%)",
                    ha="center", va="bottom", fontsize=9, fontweight="bold", zorder=5)

        # Y-axis ceiling: round up to the nearest 5000, then add 30% headroom for labels
        import math as _math6
        _y_ceil6 = _math6.ceil(_max_cnt6 / 5000) * 5000
        ax.set_ylim(0, _y_ceil6 * 1.28)
        ax.set_xticks(bar_x6)
        ax.set_xticklabels(existing_tiers, rotation=35, ha="right", fontsize=10)
        for tick, tier in zip(ax.get_xticklabels(), existing_tiers):
            tick.set_color(TIER_PALETTE.get(tier, "black"))
            tick.set_fontweight("bold")
        ax.set_xlabel("Degrader Tier", fontsize=11)
        ax.set_ylabel("Number of complexes", fontsize=11)
        ax.yaxis.grid(True, linestyle=":", alpha=0.3, zorder=0)
        ax.set_axisbelow(True)

        # Pie inset: lower-left corner over the empty space above the bars. Model
        # name and percentage read directly inside each wedge (contrast-coloured
        # against the slice fill), so no separate boxed legend is needed.
        if model_col:
            model_counts = df[model_col].value_counts().sort_index()
            _pie_colours = list(CFG.VIS_PIE_MODEL_COLOURS)[:len(model_counts)]
            _model_labels6 = []
            for raw_key in model_counts.index:
                raw_name = str(raw_key)
                if raw_name.lower().startswith("model_"):
                    raw_name = "M" + raw_name.split("_", 1)[1]
                elif raw_name.lower().startswith("model"):
                    raw_name = "M" + raw_name[5:].lstrip("_")
                _model_labels6.append(raw_name)
            ax_pie = ax.inset_axes([0.02, 0.28, 0.40, 0.56])
            wedges6, _, autotexts = ax_pie.pie(
                model_counts.values,
                colors=_pie_colours,
                autopct="%1.0f%%",
                startangle=90,
                pctdistance=0.62,
                labeldistance=None,
                wedgeprops=dict(edgecolor="white", linewidth=1.6))
            for _w6, _lab6, _col6, _at6 in zip(wedges6, _model_labels6, _pie_colours, autotexts):
                _tc6 = _text_color(_col6)
                _at6.set_fontsize(7.5)
                _at6.set_color(_tc6)
                _at6.set_fontweight("bold")
                _ang6 = np.deg2rad((_w6.theta1 + _w6.theta2) / 2.0)
                ax_pie.text(0.82 * np.cos(_ang6), 0.82 * np.sin(_ang6), _lab6,
                            ha="center", va="center", fontsize=7.5,
                            fontweight="bold", color=_tc6)
            ax_pie.set_frame_on(False)

        plt.tight_layout()
        plt.savefig(out_dir / "Figure_02_Tier_Distribution.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # --- Figure 03: Sequence Identity — Histogram (grade bands) overlaid with tier KDE lines ---
    """
    Single chart: histogram bars show the overall grade distribution (how many proteins
    per identity band); KDE lines per tier overlay on the same x-axis, revealing whether
    different catalytic tiers cluster at different identity levels.
    """
    if "identity_pct" in df.columns:
        # Use existing single-letter Alignment_Grade from CSV; only compute if missing
        if "Alignment_Grade" not in df.columns or df["Alignment_Grade"].isnull().all():
            def _simple_grade(v):
                try:
                    val = float(v)
                except (ValueError, TypeError):
                    return "I"
                return _utils_mod.get_alignment_grade(val, CFG)
            df["Alignment_Grade"] = df["identity_pct"].apply(_simple_grade)
        # Ensure single-letter format (convert "Grade A (...)" → 'A' if needed)
        _ag = df["Alignment_Grade"].astype(str)
        if _ag.str.startswith("Grade ").any():
            df["Alignment_Grade"] = _ag.str.extract(r"Grade\s+([A-I])", expand=False).fillna("I")
        id_plot = df.dropna(subset=["identity_pct"]).copy()
        id_plot["identity_pct"] = pd.to_numeric(id_plot["identity_pct"], errors="coerce").clip(0, 100)
        id_plot = id_plot.dropna(subset=["identity_pct"])

        fig, ax = plt.subplots(figsize=(13, 7))
        ax2 = ax.twinx()   # secondary y-axis for KDE density

        # Grade band background shading — sourced from CFG § 8.8
        grade_bands = list(CFG.GRADE_BANDS)
        # Per-protein grade counts (one row per unique protein)
        _prot_dedup02 = id_plot.drop_duplicates(subset=["Protein_Name"]) if "Protein_Name" in id_plot.columns else id_plot
        _n_total_prot02 = len(_prot_dedup02)
        for lo, hi, col, lbl in grade_bands:
            ax.axvspan(lo, hi, alpha=0.20, color=col, zorder=0)
            ax.text((lo + hi) / 2, 0.972, f"Grade {lbl}", ha="center", va="top",
                    fontsize=7.5, color=col, fontweight="bold",
                    transform=ax.get_xaxis_transform(),
                    bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none",
                              alpha=0.70), zorder=11)
            _cnt_g02 = int((_prot_dedup02["Alignment_Grade"] == lbl).sum())
            ax.text((lo + hi) / 2, 0.952, f"{_cnt_g02}/{_n_total_prot02}", ha="center", va="top",
                    fontsize=7.0, color=col,
                    transform=ax.get_xaxis_transform(), zorder=11)
        ax.set_axisbelow(True)

        # Primary: stacked histogram of proteins by tier (all tiers stacked per bin)
        bin_edges = np.arange(0, 101, 5)  # 5-pct-wide bins
        tier_hist_data = []
        tier_hist_labels = []
        tier_hist_colors = []
        for tier in existing_tiers:
            sub_t = id_plot[id_plot["degrader_tier"] == tier]
            if sub_t.empty:
                continue
            tier_hist_data.append(sub_t["identity_pct"].values)
            tier_hist_labels.append(f"{tier}  (n={len(sub_t):,})")
            tier_hist_colors.append(TIER_PALETTE.get(tier, "#999"))
        if tier_hist_data:
            ax.hist(tier_hist_data, bins=bin_edges, stacked=True,
                    color=tier_hist_colors, label=tier_hist_labels,
                    edgecolor="white", linewidth=0.4, alpha=0.75, zorder=2)

        """
        Secondary: KDE density line per tier — ridge (joyplot) approach with vertical offsets
        Each tier is staggered upward so curves do not overlap (ridge/joyplot style).
        """
        from scipy.stats import gaussian_kde as _gkde
        _kde_tiers_valid = [t for t in existing_tiers
                            if len(id_plot[id_plot["degrader_tier"] == t]) >= 10]
        _ridge_step = 0.008   # vertical offset per tier index
        _tt_kde_med_xy = None   # captured KDE-median dot position of the top tier (on ax2)
        for _tier_idx, tier in enumerate(_kde_tiers_valid):
            sub_t = id_plot[id_plot["degrader_tier"] == tier]
            kde_x = np.linspace(0, 100, 400)
            try:
                kde_y = _gkde(sub_t["identity_pct"])(kde_x)
            except (np.linalg.LinAlgError, ValueError):
                # Zero-variance tier (all sequences share one identity) → singular
                # covariance; skip this ridge rather than crash the whole figure.
                continue
            _tier_offset = _tier_idx * _ridge_step
            kde_y_shifted = kde_y + _tier_offset
            _tcol = TIER_PALETTE.get(tier, "#999")
            ax2.plot(kde_x, kde_y_shifted, color=_tcol,
                     linewidth=2.5, alpha=0.75, zorder=5 + _tier_idx,
                     linestyle="-")
            # Fill under each shifted KDE curve for ridge visual separation
            ax2.fill_between(kde_x, _tier_offset, kde_y_shifted,
                             color=_tcol, alpha=0.15, zorder=4 + _tier_idx)
            # Median tick mark on the KDE line
            med_x = float(sub_t["identity_pct"].median())
            med_y = float(np.atleast_1d(_gkde(sub_t["identity_pct"])(np.array([med_x])))[0])
            ax2.scatter([med_x], [med_y + _tier_offset], color=_tcol,
                        s=55, zorder=6 + _tier_idx, edgecolors="black", linewidths=0.8)
            if tier == CFG.TIER_TOP:
                _tt_kde_med_xy = (med_x, med_y + _tier_offset)

        # Anchor KDE axis at zero; show full density range with 10% headroom
        _kde_ymax = ax2.get_ylim()[1]
        ax2.set_ylim(0, _kde_ymax * 1.10)

        # Grade threshold vertical lines
        for pct_v, _lbl_v, col_v in [(90, "Grade A boundary", "#1B7837"),
                                      (60, "Grade D boundary", "#E69F00"),
                                      (40, "Grade F boundary", "#D55E00")]:
            ax.axvline(x=pct_v, color=col_v, linestyle="--", alpha=0.65,
                       linewidth=1.3, zorder=4)

        ax.set_xlim(0, 100)
        ax.set_xticks(range(0, 101, 10))
        ax.tick_params(axis="x", labelsize=9)
        # Left axis (count): steel-blue tick labels matching the grey count gridlines
        ax.tick_params(axis="y", labelsize=9, labelcolor="#2C6FAC")
        ax.yaxis.label.set_color("#2C6FAC")
        # Right axis (KDE density): amber tick labels matching the amber density gridlines
        ax2.tick_params(axis="y", labelsize=9, labelcolor="#A06000")
        ax2.yaxis.label.set_color("#A06000")
        ax.set_xlabel("Sequence Identity to Reference DeHa4 Control  (%)", fontsize=11)
        ax.set_ylabel("Number of Protein–Ligand Complexes  (stacked by tier)", fontsize=11)
        ax2.set_ylabel("KDE Density  (probability density)", fontsize=9, rotation=270, labelpad=14)

        """
        Dual-colour grid system:
        Left-axis (count) → steel-blue lines to match left tick labels
        Right-axis (KDE density) → amber lines to match right tick labels
        Grid strictly behind everything: draw grids at zorder=0, bars at zorder=2+
        ax.set_axisbelow only affects the primary axis; force zorder on both axes
        """
        ax.set_axisbelow(True)
        ax2.set_axisbelow(True)
        # Raise ax2 (KDE lines) above ax (bars) so KDE lines render on top
        ax2.set_zorder(ax.get_zorder() + 1)
        ax2.patch.set_visible(False)   # keep ax background visible
        ax.yaxis.grid(True, color="#2C6FAC", linewidth=0.55, linestyle="-", alpha=0.25, zorder=0)
        ax.xaxis.grid(False)
        ax2.yaxis.grid(True, color="#A06000", linewidth=0.55, linestyle="--", alpha=0.20, zorder=0)

        _ymax02 = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 200
        ax.annotate("Tiers converge\nnear 60% identity",
                    xy=(60, _ymax02 * 0.45),
                    xytext=(70, _ymax02 * 0.60),
                    fontsize=7.0, color="#444", style="italic",
                    arrowprops=dict(arrowstyle="->", color="#888", lw=0.9))
        # Tier_1A visibility annotation — arrow points to the Tier_1A KDE-curve median
        # (a visible dot), not the near-invisible count bar. Drawn on the KDE axis and
        # lifted above the ridge lines with an opaque box so it never hides behind them.
        _tt_sub02 = id_plot[id_plot["degrader_tier"] == CFG.TIER_TOP]
        if len(_tt_sub02) > 0 and _tt_kde_med_xy is not None:
            _tt_med_x = float(_tt_sub02["identity_pct"].median())
            ax2.annotate(
                f"{CFG.TIER_TOP}  (n={len(_tt_sub02):,})\nmedian identity {_tt_med_x:.0f}%",
                xy=_tt_kde_med_xy, xycoords="data",
                xytext=(0.55, 0.46), textcoords=ax2.transAxes,
                fontsize=7.5, ha="left", va="center",
                color=TIER_PALETTE.get(CFG.TIER_TOP, "#2E7D52"),
                fontweight="bold",
                arrowprops=dict(arrowstyle="->", lw=1.4,
                                color=TIER_PALETTE.get(CFG.TIER_TOP, "#2E7D52")),
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec=TIER_PALETTE.get(CFG.TIER_TOP, "#2E7D52"),
                          alpha=1.0, linewidth=1.1),
                zorder=30,
            )
        # Single unified legend: compound handles (bar patch + KDE line with median dot)
        from matplotlib.lines import Line2D as _L7
        _bar_hdls_02, _bar_lbls_02 = ax.get_legend_handles_labels()
        _unified_handles = []
        _unified_labels  = []
        for _bh, _bl in zip(_bar_hdls_02, _bar_lbls_02):
            # Extract tier name from label like f"Tier_1A  (n=8)"
            _tname = _bl.split("  ")[0] if "  " in _bl else _bl
            _tcol  = TIER_PALETTE.get(_tname, "#999")
            # Compound: bar patch (facecolor) + dot marker
            _unified_handles.append(_bh)
            _kdot = _L7([0],[0], marker="o", color=_tcol, linewidth=1.5,
                        markerfacecolor=_tcol, markeredgecolor="black",
                        markeredgewidth=0.5, markersize=5)
            _unified_handles.append(_kdot)
            _unified_labels.append(_bl)
            _unified_labels.append("")  # blank label for the dot entry pairs
        """
        Merge pairs: interleave bar and dot into the same rows (ncol=2 pairs)
        Simpler: build single-entry handles with dual marker using HandlerTuple
        """
        try:
            from matplotlib.legend_handler import HandlerTuple as _HT02
            _pair_handles = []
            _pair_labels  = []
            for _bh, _bl in zip(_bar_hdls_02, _bar_lbls_02):
                _tname = _bl.split("  ")[0] if "  " in _bl else _bl
                _tcol  = TIER_PALETTE.get(_tname, "#999")
                _kdot  = _L7([0],[0], marker="o", color=_tcol, linewidth=1.5,
                             markerfacecolor=_tcol, markeredgecolor="black",
                             markeredgewidth=0.5, markersize=5)
                _pair_handles.append((_bh, _kdot))
                _pair_labels.append(_tname)
            _leg7 = ax.legend(_pair_handles, _pair_labels,
                              loc="upper left", ncol=3, fontsize=8.5,
                              framealpha=0.92, fancybox=True,
                              handler_map={tuple: _HT02(ndivide=None, pad=0.2)},
                              title="Tier (bars = counts, KDE lines ● median)",
                              title_fontsize=8.0,
                              bbox_to_anchor=(0.01, 0.93))
        except Exception:
            _leg7 = ax.legend(_bar_hdls_02, _bar_lbls_02,
                              loc="upper left", ncol=3, fontsize=8.5,
                              framealpha=0.92, fancybox=True,
                              title="Tier (bars = counts, KDE lines ● median)",
                              title_fontsize=8.0,
                              bbox_to_anchor=(0.01, 0.93))
        _leg7.set_zorder(20)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_03_Alignment_Grades.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # --- Figure 04: Tier × Alignment Grade — Stacked 100% Horizontal Bar Chart ---
    """
    Each tier is a horizontal bar whose width = 100%.  The bar is divided into grade
    bands (I worst → A best) using the same green-to-brown colour ramp as Fig 07,
    making cross-figure reading intuitive.  Absolute counts are annotated inside wide
    segments; narrow segments get small external labels.
    """
    if "identity_pct" in df.columns and "degrader_tier" in df.columns:
        try:
            f11_df = df.dropna(subset=["identity_pct"]).copy()
            f11_df["identity_pct"] = pd.to_numeric(f11_df["identity_pct"], errors="coerce").clip(0, 100)
            f11_df = f11_df.dropna(subset=["identity_pct"])
            grade_cuts = list(CFG.ALIGN_GRADE_BINS)          # single source (CFG)
            _glab = CFG.ALIGN_GRADE_LABELS
            grade_lbls = []
            for _i, _lb in enumerate(_glab):
                _lo, _hi = grade_cuts[_i], grade_cuts[_i + 1]
                if _i == 0:               grade_lbls.append(f"{_lb} (<{_hi}%)")
                elif _i == len(_glab) - 1: grade_lbls.append(f"{_lb} (≥{_lo}%)")
                else:                      grade_lbls.append(f"{_lb} ({_lo}–{_hi}%)")
            # Grade colours — sourced from CFG § 8.8
            grade_colors = dict(CFG.GRADE_COLOUR_FULL)
            f11_df["Grade"] = pd.cut(f11_df["identity_pct"], bins=grade_cuts,
                                     labels=grade_lbls, right=False)
            ct11 = pd.crosstab(f11_df["degrader_tier"], f11_df["Grade"])
            ct11 = ct11.reindex(index=[t for t in existing_tiers if t in ct11.index],
                                columns=grade_lbls, fill_value=0)
            ct11_pct = ct11.div(ct11.sum(axis=1), axis=0).fillna(0) * 100

            from scipy.stats import chi2_contingency as _chi2_cont
            try:
                _chi2_f11, _p_f11, _dof_f11, _ = _chi2_cont(ct11.values)
                _chi2_str = (f"χ²({_dof_f11}) = {_chi2_f11:.1f},  "
                             f"p {'< 0.001' if _p_f11 < 0.001 else f'= {_p_f11:.3f}'}")
            except Exception:
                _chi2_str = ""

            tiers_f11 = list(ct11.index)
            n_tiers11 = len(tiers_f11)
            fig, ax = plt.subplots(figsize=(13, max(5.5, n_tiers11 * 0.85 + 2.5)))
            ax.set_axisbelow(True)

            bar_y = np.arange(n_tiers11)
            left11 = np.zeros(n_tiers11)
            # Track cumulative position for top-bar grade labels
            _grade_label_pos11 = {}   # grade → x-midpoint in top bar
            for grade in grade_lbls:
                vals = ct11_pct[grade].values if grade in ct11_pct.columns else np.zeros(n_tiers11)
                raw_n = ct11[grade].values if grade in ct11.columns else np.zeros(n_tiers11, dtype=int)
                ax.barh(bar_y, vals, left=left11, height=0.62,
                        color=grade_colors.get(grade, "#999"),
                        edgecolor="white", linewidth=0.6, zorder=2)
                # Track midpoint for top bar (last tier = highest y)
                _top_v = float(vals[-1]) if len(vals) > 0 else 0.0
                _top_l = float(left11[-1]) if len(left11) > 0 else 0.0
                if _top_v >= 1.5:
                    _grade_label_pos11[grade] = _top_l + _top_v / 2
                for yi, vi, li, ni in zip(bar_y, vals, left11, raw_n):
                    if vi < 2.5:
                        continue
                    mid_x = li + vi / 2
                    # Show both % and n for all visible segments
                    if vi >= 6:
                        txt = f"{vi:.0f}%\n(n={int(ni):,})"
                        fs = 6.5
                    else:
                        txt = f"{vi:.0f}%"
                        fs = 6
                    txt_col = "white" if grade_colors.get(grade, "#999")[1:] < "AAAAAA" else "#111"
                    ax.text(mid_x, yi, txt, ha="center", va="center",
                            fontsize=fs, color=txt_col, fontweight="bold", zorder=3,
                            linespacing=1.1)
                left11 = left11 + vals

            # ── Grade bands identical to Fig02 — fixed identity-% positions ──
            """
            grade_bands = [(lo, hi, col, lbl), ...] with same positions as Fig02.
            Both figures share a 0–100 x-axis so the visual layout matches exactly.
            """
            _g11_lbl_to_full = {g.split(" ")[0]: g for g in grade_lbls}
            _g11_total       = ct11.sum(axis=0)
            _n_total11       = int(ct11.values.sum())
            for lo, hi, col, lbl in grade_bands:
                # Background shading — same alpha/style as Fig02
                ax.axvspan(lo, hi, alpha=0.20, color=col, zorder=0)
                # Dashed vertical divider at left boundary (skip x=0)
                if lo > 0:
                    ax.axvline(lo, color=col, lw=1.0, ls="--", alpha=0.55, zorder=3)
                # Grade header label at midpoint — same style as Fig02
                _mid = (lo + hi) / 2
                ax.text(_mid, 1.015, f"Grade {lbl}",
                        ha="center", va="bottom", fontsize=7.5, color=col,
                        fontweight="bold", transform=ax.get_xaxis_transform(),
                        clip_on=False,
                        bbox=dict(boxstyle="round,pad=0.08", fc="white",
                                  ec="none", alpha=0.70), zorder=11)
                _full_lbl = _g11_lbl_to_full.get(lbl, "")
                _gn = int(_g11_total.get(_full_lbl, 0))
                ax.text(_mid, 1.001, f"{_gn:,}/{_n_total11:,}",
                        ha="center", va="bottom", fontsize=7.0, color=col,
                        transform=ax.get_xaxis_transform(), clip_on=False, zorder=11)

            # Tier row labels
            _tier_n11 = ct11.sum(axis=1)
            ax.set_yticks(bar_y)
            ax.set_yticklabels([f"{t}  (n={_tier_n11.loc[t]:,})" for t in tiers_f11], fontsize=9.5)
            for tick, tier in zip(ax.get_yticklabels(), tiers_f11):
                tick.set_color(TIER_PALETTE.get(tier, "black"))
                tick.set_fontweight("bold")

            # X-axis ticks every 10% — identical to Fig02
            ax.set_xlim(0, 100)
            ax.set_xticks(range(0, 101, 10))
            ax.set_xticklabels([f"{v}%" for v in range(0, 101, 10)],
                               fontsize=8, rotation=45, ha="right")
            ax.set_xlabel("Proportion of complexes in each sequence identity grade  (%)",
                          fontsize=11)
            ax.set_ylabel("Degrader Tier", fontsize=11)
            ax.xaxis.grid(True, color="#DCDCDC", linewidth=0.65, zorder=0)
            ax.set_axisbelow(True)

            if _chi2_str:
                ax.text(0.99, 0.01, _chi2_str, transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=8.5, style="italic",
                        color="#333",
                        bbox=dict(boxstyle="round,pad=0.25", fc="white",
                                  ec="#ccc", alpha=0.88, linewidth=0.6))
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_04_Tier_Grade_Distribution.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Figure 04 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # ── Companion thumbnail setup (shared by Figs 04b, 12b, 13b, 17b) ─────────
    # Source the PyMOL structure thumbnails from the MD-ready cohort (MD_Selected,
    # SECTION 18) so only the ~9-10 complexes actually taken to MD are rendered —
    # not all Tier_1A. Falls back to the top tier if the flag is absent (older CSV).
    _md_col = getattr(CFG, "MD_SELECTED_COL", "MD_Selected")
    if _md_col in df.columns:
        _sel = df[_md_col].astype(str).str.lower().isin(["true", "1", "1.0"])
        _pa = df[_sel].reset_index(drop=True)
    elif "degrader_tier" in df.columns:
        _pa = df[df["degrader_tier"] == CFG.TIER_TOP].reset_index(drop=True)
    else:
        _pa = pd.DataFrame()
    # Rank-sort so the capped thumbnail set (CFG.VIS_MAX_THUMBNAILS) takes the best
    # entries; Scientific_Rank ascending = best first, else Boltz confidence.
    if not _pa.empty:
        if "Scientific_Rank" in _pa.columns:
            _pa = _pa.sort_values("Scientific_Rank", ascending=True).reset_index(drop=True)
        elif "Boltz_Model_Confidence" in _pa.columns:
            _pa = _pa.sort_values("Boltz_Model_Confidence", ascending=False).reset_index(drop=True)
    _imgs = []
    try:
        if not _pa.empty:
            _pred_jobs = out_dir.parent / "1_Boltz2_Production" / "4_Prediction_Jobs"
            _thumb_dir = out_dir / "_tt_thumbnails"
            _thumb_dir.mkdir(parents=True, exist_ok=True)
            _imgs = _tt_render_all(_pa, _pred_jobs, _thumb_dir, reporter)
    except Exception:
        pass
    _tt_has_imgs = not _pa.empty and sum(i is not None for i in _imgs) > 0

    reporter.section("  Folder 03_AI_Confidence_Quality — Boltz-2 confidence + pTM/ipTM")
    # --- Figure 05: AI Quality — Boltz Confidence (boxes) + pTM & ipTM (lines) overlaid ---
    """
    Single chart: box plots show per-tier distribution of Boltz Model Confidence (primary
    y-axis). Overlaid connected dot-lines show per-tier median pTM and ipTM on the same
    scale (both 0–1), making agreement and divergence between the three AI metrics visible.
    """
    _has_conf  = "Boltz_Model_Confidence" in df.columns and "degrader_tier" in df.columns
    _has_ptm   = "ptm" in df.columns and "iptm" in df.columns and "degrader_tier" in df.columns
    if _has_conf or _has_ptm:
        try:
            fig, ax = plt.subplots(figsize=(13, 7))

            # Quality zone bands — edges/colours from CFG (single source)
            _cbh, _cba, _cbc = CFG.CONF_BAND_HIGH, CFG.CONF_BAND_ACCEPTABLE, CFG.CONF_BAND_COLOURS
            ax.axhspan(_cbh, 1.01, alpha=0.10, color=_cbc["high"], zorder=0)
            ax.axhspan(_cba, _cbh, alpha=0.09, color=_cbc["acceptable"], zorder=0)
            ax.axhspan(0.00, _cba, alpha=0.07, color=_cbc["below"], zorder=0)
            ax.axhline(y=_cbh, color=_cbc["high"], linestyle="--", linewidth=1.3, alpha=0.80, zorder=1)
            ax.axhline(y=_cba, color=_cbc["acceptable"], linestyle=":", linewidth=1.3, alpha=0.80, zorder=1)
            """
            Zone labels — placed in DATA y-coordinates (get_yaxis_transform) AFTER ylim is set,
            so they land at the correct zone midpoints regardless of the dynamic y-axis range.
            (Labels are defined here as closures and drawn below after ax.set_ylim is called.)
            """

            tier_x = {tier: i for i, tier in enumerate(existing_tiers)}

            # Boltz Confidence — box plots (primary, fills the background)
            if _has_conf:
                sns.boxplot(data=df, x="degrader_tier", y="Boltz_Model_Confidence",
                            order=existing_tiers, palette=TIER_PALETTE, linewidth=1.2,
                            flierprops=dict(marker="o", markersize=2, alpha=0.25),
                            width=0.55, zorder=2, ax=ax)
                # Median text above each box
                tier_order_f10 = list(existing_tiers)
                medians_f10 = []
                for i, tier in enumerate(tier_order_f10):
                    med = df.loc[df["degrader_tier"] == tier, "Boltz_Model_Confidence"].median()
                    medians_f10.append(float(med) if not np.isnan(med) else float("nan"))
                    if not np.isnan(med):
                        ax.text(i, med + 0.003, f"{med:.3f}", ha="center", va="bottom",
                                fontsize=8, color="black", fontweight="bold", zorder=7)
                # Tier_5_Decoy overconfidence annotation
                _decoy_vals = [v for t, v in zip(tier_order_f10, medians_f10) if t == CFG.TIER_DECOY]
                if _decoy_vals and _decoy_vals[0] > 0.95:
                    _di = tier_order_f10.index(CFG.TIER_DECOY)
                    ax.annotate(f"{CFG.TIER_DECOY}\noverconfidence",
                                xy=(_di, _decoy_vals[0]),
                                xytext=(_di + 0.4, _decoy_vals[0] - 0.02),
                                fontsize=7.5, color="#B22222", style="italic",
                                arrowprops=dict(arrowstyle="->", color="#B22222", lw=0.9),
                                bbox=dict(boxstyle="round,pad=0.2", fc="#FFF5F5",
                                          ec="#B22222", alpha=0.9, lw=0.8))

            # pTM — connected median line with CI band
            if _has_ptm:
                ptm_xs, ptm_meds, ptm_lo, ptm_hi = [], [], [], []
                iptm_xs, iptm_meds, iptm_lo, iptm_hi = [], [], [], []
                for tier in existing_tiers:
                    sub = df.loc[df["degrader_tier"] == tier].dropna(subset=["ptm", "iptm"])
                    if len(sub) < 3:
                        continue
                    x_pos = tier_x[tier]
                    for metric, xs_l, meds_l, lo_l, hi_l, col in [
                        ("ptm", ptm_xs, ptm_meds, ptm_lo, ptm_hi, "#0072B2"),
                        ("iptm", iptm_xs, iptm_meds, iptm_lo, iptm_hi, "#CC79A7"),
                    ]:
                        vals = sub[metric].dropna()
                        med_v = float(vals.median())
                        lo_ci, hi_ci = _median_ci95(vals)
                        xs_l.append(x_pos)
                        meds_l.append(med_v)
                        lo_l.append(lo_ci)
                        hi_l.append(hi_ci)

                for xs_l, meds_l, lo_l, hi_l, col, lbl, mk in [
                    (ptm_xs, ptm_meds, ptm_lo, ptm_hi, "#0072B2", "pTM  (fold confidence)", "s"),
                    (iptm_xs, iptm_meds, iptm_lo, iptm_hi, "#CC79A7", "ipTM  (interface confidence)", "^"),
                ]:
                    if not xs_l:
                        continue
                    ax.plot(xs_l, meds_l, color=col, linewidth=2.0, zorder=5,
                            marker=mk, markersize=8, markeredgecolor="black",
                            markeredgewidth=0.8, label=lbl)
                    ax.fill_between(xs_l, lo_l, hi_l, color=col, alpha=0.18, zorder=4)

            # Y-axis zoom: start just below data minimum
            all_vals_f10 = []
            if _has_conf:
                all_vals_f10.extend(df["Boltz_Model_Confidence"].dropna().tolist())
            if _has_ptm:
                all_vals_f10.extend(df["ptm"].dropna().tolist())
                all_vals_f10.extend(df["iptm"].dropna().tolist())
            if all_vals_f10:
                y_lo_f10 = max(0.0, float(np.percentile(all_vals_f10, 1)) - 0.02)
                ax.set_ylim(y_lo_f10, 1.01)
            else:
                y_lo_f10 = 0.0
            # Broken-axis indicator at y_lo_f10 (axis is truncated — not starting at 0)
            if y_lo_f10 > 0.01:
                _bax_kw = dict(transform=ax.transAxes, color="#333", clip_on=False,
                               linewidth=2.2)
                _d = 0.022   # diagonal tick half-length in axes coords
                _bax_y = 0.0  # bottom of axes
                # Left y-axis break marks
                ax.plot((-_d, +_d), (_bax_y - _d, _bax_y + _d), **_bax_kw)
                ax.plot((-_d * 2.2, +_d * 0.2), (_bax_y - _d, _bax_y + _d), **_bax_kw)
                # Right side break marks
                ax.plot((1 - _d, 1 + _d), (_bax_y - _d, _bax_y + _d), **_bax_kw)
                ax.plot((1 - _d * 2.2, 1 + _d * 0.2), (_bax_y - _d, _bax_y + _d), **_bax_kw)
                ax.text(0.01, 0.015, f"↕ axis break (starts at {CFG.CONF_BAND_ACCEPTABLE:.2f})",
                        transform=ax.transAxes,
                        fontsize=7.5, color="#555", ha="left", va="bottom",
                        style="italic")
            """
            Zone labels in DATA coordinates — each label at the midpoint of its
            VISIBLE zone slice (zones clipped to the actual ylim after set_ylim).
            get_yaxis_transform(): x=axes fraction [0,1], y=data coordinates
            """
            _yt10 = ax.get_yaxis_transform()
            _y_top10 = 1.01
            """
            High confidence zone: 0.90–1.01 — fixed at top-right via transAxes so it is
            always visible regardless of ylim and never lands on Tier_1A box bodies.
            """
            _hc_vis_lo10 = max(0.90, y_lo_f10)
            if _hc_vis_lo10 < _y_top10 - 0.005:
                ax.text(0.01, 0.98, f"High confidence (≥{CFG.CONF_BAND_HIGH:.2f})", color="#007A50",
                        fontsize=8, ha="left", va="top", fontweight="bold",
                        style="italic", transform=ax.transAxes, zorder=6,
                        bbox=dict(boxstyle="round,pad=0.15", fc="#D6F5EB", ec=CFG.CONF_BAND_COLOURS["high"],
                                  alpha=0.85, linewidth=0.6))
            # Acceptable zone: 0.80–0.90
            _acc_vis_lo10 = max(0.80, y_lo_f10)
            _acc_vis_hi10 = min(0.90, _y_top10)
            if _acc_vis_hi10 > _acc_vis_lo10 + 0.005:
                _acc_mid10 = (_acc_vis_lo10 + _acc_vis_hi10) / 2
                ax.text(0.01, _acc_mid10, f"Acceptable ({CFG.CONF_BAND_ACCEPTABLE:.2f}–{CFG.CONF_BAND_HIGH:.2f})", color="#8A6000",
                        fontsize=8, ha="left", va="center", fontweight="bold",
                        style="italic", transform=_yt10, zorder=6,
                        bbox=dict(boxstyle="round,pad=0.15", fc="#FFF3CC", ec=CFG.CONF_BAND_COLOURS["acceptable"],
                                  alpha=0.85, linewidth=0.6))
            # Below threshold zone: y_lo–0.80
            _bel_vis_hi10 = min(0.80, _y_top10)
            if _bel_vis_hi10 > y_lo_f10 + 0.005:
                _bel_mid10 = (y_lo_f10 + _bel_vis_hi10) / 2
                ax.text(0.01, _bel_mid10, f"Below threshold (<{CFG.CONF_BAND_ACCEPTABLE:.2f})", color="#A03000",
                        fontsize=8, ha="left", va="center", fontweight="bold",
                        style="italic", transform=_yt10, zorder=6,
                        bbox=dict(boxstyle="round,pad=0.15", fc="#FDECEA", ec=CFG.CONF_BAND_COLOURS["below"],
                                  alpha=0.85, linewidth=0.6))

            ax.set_xticks(range(len(existing_tiers)))
            ax.set_xticklabels(existing_tiers, rotation=40, ha="right", fontsize=9)
            for _tick10x, _tier10x in zip(ax.get_xticklabels(), existing_tiers):
                _tick10x.set_color(TIER_PALETTE.get(_tier10x, "black"))
                _tick10x.set_fontweight("bold")
            ax.tick_params(axis="y", labelsize=9, labelcolor="#2C6FAC")
            ax.yaxis.label.set_color("#2C6FAC")
            ax.yaxis.grid(True, color="#2C6FAC", linewidth=0.5, linestyle="-", alpha=0.20, zorder=0)
            ax.set_axisbelow(True)
            ax.set_xlabel("Degrader Tier", fontsize=11)
            ax.set_ylabel("AI Quality Score  (Boltz Confidence / pTM / ipTM;  0–1 scale)", fontsize=11)

            # Legend: single horizontal row, top-right, above data
            from matplotlib.patches import Patch as _Patch10
            conf_patch = _Patch10(facecolor="#AAAAAA", edgecolor="black", linewidth=0.8,
                                  label="Boltz Confidence  (box)")
            handles_f10, labels_f10 = ax.get_legend_handles_labels()
            _all_h10 = [conf_patch] + handles_f10
            _all_l10 = ["Boltz Confidence  (box)"] + labels_f10
            _leg10 = ax.legend(_all_h10, _all_l10,
                               loc="upper right", fontsize=8,
                               framealpha=0.92, fancybox=True,
                               ncol=len(_all_h10))   # all in one row
            _leg10.set_zorder(20)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_05a_AI_Quality_Assessment.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Figure 05 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # Figure 05b: PA companion (Confidence × ipTM landscape)
    if _tt_has_imgs:
        _fig_05b_tt_ai_quality(df, _pa, _imgs, out_dir, reporter)

    # --- Figure 06: pTM vs ipTM — 2-D Scatter (REDESIGN) ---
    """
    Each complex is a point: X = pTM (fold confidence), Y = ipTM (interface confidence).
    Colour = degrader tier. The identity diagonal (pTM = ipTM) divides the space:
      ABOVE diagonal → interface more confident than fold = preferred state
      BELOW diagonal → fold more confident than interface = potential decoy risk
    Per-tier median diamonds summarise cluster positions.
    This is fundamentally more informative than paired box plots.
    """
    if "ptm" in df.columns and "iptm" in df.columns and "degrader_tier" in df.columns:
        try:
            f18_df = df.dropna(subset=["ptm", "iptm"]).copy()
            f18_df["ptm"]  = pd.to_numeric(f18_df["ptm"],  errors="coerce")
            f18_df["iptm"] = pd.to_numeric(f18_df["iptm"], errors="coerce")
            f18_df = f18_df.dropna(subset=["ptm", "iptm"])

            all_vals18 = np.concatenate([f18_df["ptm"].values, f18_df["iptm"].values])
            all_vals18 = all_vals18[~np.isnan(all_vals18)]
            _lo18 = max(0.0, float(np.percentile(all_vals18, 0.5)) - 0.01) if len(all_vals18) else 0.0
            _hi18 = min(1.01, float(np.percentile(all_vals18, 99.5)) + 0.01) if len(all_vals18) else 1.01

            from matplotlib.gridspec import GridSpec as _GS18
            fig = plt.figure(figsize=(13, 11))
            _gs18 = _GS18(2, 2, figure=fig,
                          width_ratios=[4, 1], height_ratios=[1, 4],
                          hspace=0.04, wspace=0.04)
            ax       = fig.add_subplot(_gs18[1, 0])  # main scatter
            ax_top   = fig.add_subplot(_gs18[0, 0], sharex=ax)  # top marginal
            ax_right = fig.add_subplot(_gs18[1, 1], sharey=ax)  # right marginal
            fig.add_subplot(_gs18[0, 1]).set_visible(False)      # empty corner
            ax.set_axisbelow(True)

            # ── Three coloured zone backgrounds ────────────────────────────────
            """
            Zone 1 (above diagonal): ipTM > pTM — preferred state, interface well modelled
            Zone 2 (diagonal band ±0.03): pTM ≈ ipTM — balanced confidence
            Zone 3 (below diagonal): pTM > ipTM — fold stronger than interface (decoy risk)
            """
            _diag_pts = np.linspace(_lo18, _hi18, 300)
            from matplotlib.patches import Polygon as _Poly18
            _band = 0.03
            # Zone 1: above diagonal (ipTM > pTM)
            _verts_up = ([(x, x + _band) for x in _diag_pts] +
                         [(_hi18, _hi18), (_lo18, _hi18)])
            ax.add_patch(_Poly18(_verts_up, closed=True,
                                 facecolor="#D4EFDF", alpha=0.45, zorder=0, linewidth=0))
            # Zone 3: below diagonal (pTM > ipTM)
            _verts_dn = ([(x, x - _band) for x in _diag_pts] +
                         [(_hi18, _lo18), (_lo18, _lo18)])
            ax.add_patch(_Poly18(_verts_dn, closed=True,
                                 facecolor="#FADBD8", alpha=0.45, zorder=0, linewidth=0))
            # Zone 2 (diagonal band) is implicitly the white strip between the two coloured zones

            # Identity diagonal line
            ax.plot(_diag_pts, _diag_pts, color="#555", linestyle="--", linewidth=1.3,
                    alpha=0.60, zorder=2)
            # Band boundary dashes
            ax.plot(_diag_pts, _diag_pts + _band, color="#555", linestyle=":",
                    linewidth=0.7, alpha=0.40, zorder=2)
            ax.plot(_diag_pts, _diag_pts - _band, color="#555", linestyle=":",
                    linewidth=0.7, alpha=0.40, zorder=2)

            """
            Zone text — placed BELOW the upper-left legend to avoid overlap.
            Legend with ncol=2 is ~2 rows tall; empirically offset to 0.60 axes-y.
            """
            ax.text(0.04, 0.60, "ipTM > pTM\nInterface well-modelled",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=8.5, color="#1A7A4A", fontweight="bold", style="italic",
                    bbox=dict(boxstyle="round,pad=0.22", fc="#D4EFDF", ec="#1A7A4A",
                              alpha=0.80, linewidth=0.7))
            ax.text(0.97, 0.97, "pTM = ipTM\n(±0.03 band)",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=7.5, color="#555", style="italic")
            ax.text(0.96, 0.10, "pTM > ipTM\nFold > interface\n(potential decoy)",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=8.5, color="#A93226", fontweight="bold", style="italic",
                    bbox=dict(boxstyle="round,pad=0.22", fc="#FADBD8", ec="#A93226",
                              alpha=0.80, linewidth=0.7))

            """
            Rendering strategy:
            Tier_3 tier (n~34k) → hexbin background density (YlOrBr)
            All other tiers   → scatter with tier colour, drawn on top in order
            This keeps tier colours distinguishable while showing the density of the bulk.
            """
            for tier in sorted(existing_tiers, key=lambda t: t == CFG.TIER_TOP):
                sub18 = f18_df[f18_df["degrader_tier"] == tier]
                if sub18.empty:
                    continue
                ptm_vals  = sub18["ptm"].values
                iptm_vals = sub18["iptm"].values
                col18     = TIER_PALETTE.get(tier, "#999")
                if tier == CFG.TIER_ORDER[4]:
                    ax.hexbin(ptm_vals, iptm_vals, gridsize=45, cmap="YlOrBr",
                              mincnt=1, alpha=0.65, zorder=1)
                else:
                    _is_pa   = tier == CFG.TIER_TOP
                    _is_top  = tier in (CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2])
                    ax.scatter(ptm_vals, iptm_vals,
                               c=col18,
                               alpha=0.92 if _is_pa else (0.55 if _is_top else 0.20),
                               s=110 if _is_pa else (40 if _is_top else 10),
                               edgecolors="black" if _is_pa else "none",
                               linewidths=1.2 if _is_pa else 0,
                               marker="*" if _is_pa else "o",
                               zorder=10 if _is_pa else (5 if _is_top else 3),
                               label=None)

            """
            Per-tier median diamonds — NO inline text annotation;
            coordinates are embedded in the legend label instead
            """
            from matplotlib.lines import Line2D as _L18
            _leg18 = []
            # Header entries: show both symbol types so the reader knows what each glyph means
            _leg18.append(_L18([0], [0], marker="o", color="w",
                               markerfacecolor="#888888", markeredgecolor="none",
                               markersize=7, alpha=0.55,
                               label="Individual complex  (hexbin/scatter)"))
            _leg18.append(_L18([0], [0], marker="D", color="w",
                               markerfacecolor="#888888", markeredgecolor="black",
                               markersize=9, markeredgewidth=1.0,
                               label="Tier median  (◆ diamond)"))
            for tier in existing_tiers:
                sub18 = f18_df[f18_df["degrader_tier"] == tier]
                if len(sub18) < 2:
                    continue
                mx18, my18 = float(sub18["ptm"].median()), float(sub18["iptm"].median())
                col18 = TIER_PALETTE.get(tier, "#999")
                ax.scatter([mx18], [my18], marker="D", s=160,
                           color=col18, edgecolors="black", linewidths=1.3, zorder=8)
                # Tier_1A callout arrow on main scatter
                if tier == CFG.TIER_TOP:
                    ax.annotate(
                        f"{CFG.TIER_TOP}\n(n={len(sub18)}, ★)",
                        xy=(mx18, my18),
                        xytext=(mx18 - 0.035, my18 + 0.015),
                        fontsize=8, color=col18, fontweight="bold", zorder=12,
                        arrowprops=dict(arrowstyle="->", color=col18, lw=1.2,
                                        shrinkA=6, shrinkB=4),
                        bbox=dict(boxstyle="round,pad=0.25", fc="white",
                                  ec=col18, alpha=0.92, linewidth=1.2))
                # Median coords live in the legend label — no floating text cluttering the plot
                _leg18.append(_L18([0], [0], marker="D", color="w",
                                   markerfacecolor=col18, markeredgecolor="black",
                                   markersize=9, markeredgewidth=1.0,
                                   label=f"{tier}  ◆({mx18:.3f}, {my18:.3f})  n={len(sub18):,}"))

            # Marginal distributions — KDE lines (smoother than histograms)
            from scipy.stats import gaussian_kde as _gkde18
            existing_tiers18 = [t for t in existing_tiers if t in f18_df["degrader_tier"].values]
            _kde_x18 = np.linspace(_lo18, _hi18, 300)
            for tier in existing_tiers18:
                sub18_m = f18_df[f18_df["degrader_tier"] == tier]
                _min_kde = 1 if tier == CFG.TIER_TOP else 10
                if len(sub18_m) < _min_kde:
                    continue
                _col18m = TIER_PALETTE.get(tier, "#999")
                _is_pa18 = tier == CFG.TIER_TOP
                _lw18    = 2.8 if _is_pa18 else (2.0 if tier == CFG.TIER_ORDER[1] else 1.3)
                _alp18   = 0.95 if _is_pa18 else 0.80
                _zo18    = 10  if _is_pa18 else 5
                try:
                    _kde_ptm  = _gkde18(sub18_m["ptm"].values)(_kde_x18)
                    _kde_iptm = _gkde18(sub18_m["iptm"].values)(_kde_x18)
                    ax_top.plot(_kde_x18, _kde_ptm,  color=_col18m, lw=_lw18,
                                alpha=_alp18, zorder=_zo18)
                    ax_top.fill_between(_kde_x18, _kde_ptm,
                                        alpha=0.25 if _is_pa18 else 0.10,
                                        color=_col18m, zorder=_zo18 - 1)
                    ax_right.plot(_kde_iptm, _kde_x18, color=_col18m, lw=_lw18,
                                  alpha=_alp18, zorder=_zo18)
                    ax_right.fill_betweenx(_kde_x18, _kde_iptm,
                                           alpha=0.25 if _is_pa18 else 0.10,
                                           color=_col18m, zorder=_zo18 - 1)
                except Exception:
                    pass
            # Annotate Tier_1A KDE peak on both marginal axes
            _tt_col_f05 = TIER_PALETTE.get(CFG.TIER_TOP, "#009E73")
            _tt_sub_f05 = f18_df[f18_df["degrader_tier"] == CFG.TIER_TOP]
            if len(_tt_sub_f05) >= 2:
                try:
                    _kde_tt_ptm  = _gkde18(_tt_sub_f05["ptm"].values)(_kde_x18)
                    _kde_tt_iptm = _gkde18(_tt_sub_f05["iptm"].values)(_kde_x18)
                    _peak_ptm_x  = float(_kde_x18[np.argmax(_kde_tt_ptm)])
                    _peak_ptm_y  = float(np.max(_kde_tt_ptm))
                    _peak_iptm_x = float(_kde_x18[np.argmax(_kde_tt_iptm)])
                    _peak_iptm_y = float(np.max(_kde_tt_iptm))
                    # Top KDE: annotation placed BELOW the peak (away from legend box)
                    ax_top.annotate(
                        f"{CFG.TIER_TOP} ▲",
                        xy=(_peak_ptm_x, _peak_ptm_y),
                        xytext=(_peak_ptm_x + 0.008, _peak_ptm_y * 0.40),
                        fontsize=7, color=_tt_col_f05, fontweight="bold",
                        ha="left", zorder=15,
                        arrowprops=dict(arrowstyle="->", color=_tt_col_f05, lw=1.0,
                                        shrinkA=2, shrinkB=2),
                        bbox=dict(boxstyle="round,pad=0.14", fc="white",
                                  ec=_tt_col_f05, alpha=0.95, linewidth=0.8))
                    # Right KDE: annotation to the left of peak where density is low
                    ax_right.annotate(
                        f"{CFG.TIER_TOP} ▲",
                        xy=(_peak_iptm_y, _peak_iptm_x),
                        xytext=(_peak_iptm_y * 0.45, _peak_iptm_x + 0.008),
                        fontsize=7, color=_tt_col_f05, fontweight="bold",
                        ha="right", zorder=15,
                        arrowprops=dict(arrowstyle="->", color=_tt_col_f05, lw=1.0,
                                        shrinkA=2, shrinkB=2),
                        bbox=dict(boxstyle="round,pad=0.14", fc="white",
                                  ec=_tt_col_f05, alpha=0.95, linewidth=0.8))
                except Exception:
                    pass
            plt.setp(ax_top.get_xticklabels(), visible=False)
            plt.setp(ax_right.get_yticklabels(), visible=False)
            # Marginal density panels: top = pTM (shared x), right = ipTM (shared y).
            # No panel titles — orientation follows the shared 2-D axes.
            ax_top.set_ylabel("pTM density", fontsize=7, labelpad=2)
            ax_right.set_xlabel("ipTM density", fontsize=7, labelpad=2)
            ax_top.set_axisbelow(True)
            ax_right.set_axisbelow(True)
            ax_top.tick_params(axis="y", labelsize=6)
            ax_right.tick_params(axis="x", labelsize=6)
            for _sp in list(ax_top.spines.values()) + list(ax_right.spines.values()):
                _sp.set_visible(False)

            ax.set_xlim(_lo18, _hi18)
            ax.set_ylim(_lo18, _hi18)
            ax.set_xlabel("pTM — Global Fold Confidence  (0–1)", fontsize=11)
            ax.set_ylabel("ipTM — Interface Confidence  (0–1)", fontsize=11)
            ax.tick_params(axis="both", labelsize=9)
            ax.grid(True, color="#EBEBEB", linewidth=0.5, alpha=0.6, zorder=1)
            # Legend: placed in the top KDE marginal (left side) to keep scatter uncluttered
            ax_top.legend(handles=_leg18, loc="upper left", ncol=3, fontsize=7.0,
                          framealpha=0.92, fancybox=True,
                          title="Tier  (◆ = median pTM, ipTM)", title_fontsize=7.0,
                          bbox_to_anchor=(0.0, 1.0))
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_06_pTM_vs_ipTM_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Figure 06 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    reporter.section("  Folder 04_Catalytic_Geometry_and_Mechanism — structural + mechanistic geometry")
    # --- Figure 07: Active Site RMSD vs. Reference Control by Tier ---
    # Column alias detection: try multiple plausible names before skipping
    _rmsd17_aliases = [
        "Active_Site_RMSD_to_Control",
        "Active_Site_RMSD", "active_site_rmsd", "ActiveSite_RMSD", "RMSD_ActiveSite",
        "rmsd_active_site", "binding_site_rmsd", "Binding_Site_RMSD",
        "RMSD", "rmsd", "structural_rmsd", "Structural_RMSD",
        "site_rmsd", "Site_RMSD", "AS_RMSD", "as_rmsd",
        "reference_rmsd", "Reference_RMSD", "RMSD_to_reference",
    ]
    # Broad partial-match fallback: any column with 'rmsd' in the name
    _rmsd17_col = next((c for c in _rmsd17_aliases if c in df.columns), None)
    if _rmsd17_col is None:
        _rmsd17_col = next((c for c in df.columns if "rmsd" in c.lower()), None)
    if _rmsd17_col is None:
        _rmsd_candidates = [c for c in df.columns if "rmsd" in c.lower() or "active_site" in c.lower()]
        reporter.log(f"  ! Figure 07 skipped: no RMSD column found. "
                     f"All column names tried: {_rmsd17_aliases}. "
                     f"Partial matches in CSV (containing 'rmsd'/'active_site'): {_rmsd_candidates}. "
                     f"Full column list: {list(df.columns)}")
    if _rmsd17_col and "degrader_tier" in df.columns:
        try:
            if _rmsd17_col != "Active_Site_RMSD":
                df = df.rename(columns={_rmsd17_col: "Active_Site_RMSD"})
            f17_df = df.dropna(subset=["Active_Site_RMSD"]).copy()
            f17_df["Active_Site_RMSD"] = pd.to_numeric(f17_df["Active_Site_RMSD"], errors="coerce")
            f17_df = f17_df.dropna(subset=["Active_Site_RMSD"])
            """
            Show ALL data — no percentile clip so small-n tiers (e.g. Tier_1A n=5)
            retain every value.  y_ceil caps the display at 4 Å, but is extended if any
            small-n tier has a value beyond that so its dots always land inside the axes.
            """
            f17_plot = f17_df.copy()
            valid_tiers_f17 = [t for t in existing_tiers if t in f17_plot["degrader_tier"].values]
            _sn_max17 = 0.0
            for _t17s in valid_tiers_f17:
                _sv17 = f17_plot.loc[f17_plot["degrader_tier"] == _t17s, "Active_Site_RMSD"]
                if 0 < len(_sv17) <= 30:
                    _sn_max17 = max(_sn_max17, float(_sv17.max()))
            _raw_max = float(f17_df["Active_Site_RMSD"].max())
            _re, _ra, _rd = CFG.RMSD_BAND_EXCELLENT, CFG.RMSD_BAND_ACCEPTABLE, CFG.RMSD_BAND_DIVERGED
            _violin_cap = _rd   # violin body capped at the diverged band (CFG)
            _raw_pct98  = float(f17_df["Active_Site_RMSD"].quantile(0.98)) if len(f17_df) > 0 else _rd
            y_ceil      = max(_raw_pct98 * 1.08, _rd + 0.6)  # show 98th percentile + headroom
            # Restrict the violin/strip data to in-cap values by DROPPING (not
            # clamping) rows above _violin_cap. Clamping piled a point-mass at the
            # cap that the KDE smoothed into a false density bulge, and made the
            # per-point strips render outliers twice (once clamped at the cap, once
            # at their true value in the outlier scatter below). Values above the
            # cap are shown once, at their true coordinates, by that outlier scatter.
            f17_plot = f17_plot[f17_plot["Active_Site_RMSD"] <= _violin_cap].copy()
            # Compute % of data in each zone (across all tiers together, using unclipped values)
            total_n = len(f17_df)
            pct_green  = (f17_df["Active_Site_RMSD"] <= _re).sum() / total_n * 100
            pct_yellow = ((f17_df["Active_Site_RMSD"] > _re) & (f17_df["Active_Site_RMSD"] <= _ra)).sum() / total_n * 100
            pct_red    = (f17_df["Active_Site_RMSD"] > _ra).sum() / total_n * 100
            fig, ax = plt.subplots(figsize=(12, 7))
            # Coloured quality zones — colours from CFG (single source); RMSD edges 1.0/2.0 Å
            _cbc = CFG.CONF_BAND_COLOURS
            ax.axhspan(0.0, _re, alpha=0.13, color=_cbc["high"], zorder=0)   # green  = Excellent
            ax.axhspan(_re, _ra, alpha=0.13, color=_cbc["acceptable"], zorder=0)   # amber  = Acceptable
            ax.axhspan(_ra, y_ceil + 0.1, alpha=0.13, color=_cbc["below"], zorder=0)  # red = Diverged
            # Zone boundary lines
            ax.axhline(y=_re, color=_cbc["high"], linestyle="--", alpha=0.7, linewidth=1.1)
            ax.axhline(y=_ra, color=_cbc["below"], linestyle="--", alpha=0.7, linewidth=1.1)
            # Violin plot — uses clipped data so bodies remain within the continuous axis
            sns.violinplot(data=f17_plot, x="degrader_tier", y="Active_Site_RMSD",
                           order=valid_tiers_f17, palette=TIER_PALETTE,
                           inner="quartile", cut=0, linewidth=1.1, ax=ax, zorder=3)
            # Connected median trend across tiers (markers joined by a line; value
            # labels boxed just below each node) — mirrors the Fig 16 penalty line.
            _med_xs17, _med_ys17, _med_raw17 = [], [], []
            for _im17, _tier17m in enumerate(valid_tiers_f17):
                _med17m = f17_df.loc[f17_df["degrader_tier"] == _tier17m,
                                     "Active_Site_RMSD"].median()
                if np.isnan(_med17m):
                    continue
                _med_xs17.append(_im17)
                _med_ys17.append(min(float(_med17m), _violin_cap))
                _med_raw17.append(float(_med17m))
            if _med_xs17:
                ax.plot(_med_xs17, _med_ys17, color="#CC0000", linewidth=2.2,
                        marker="D", markersize=9, markeredgecolor="black",
                        markeredgewidth=0.9, zorder=8, label="Median RMSD (per tier)")
                for _xi17, _yi17, _rv17 in zip(_med_xs17, _med_ys17, _med_raw17):
                    ax.text(_xi17, _yi17 - 0.13, f"{_rv17:.2f}", ha="center", va="top",
                            fontsize=7.5, color="#CC0000", fontweight="bold", zorder=9,
                            bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                      ec="#CC0000", alpha=0.9, linewidth=0.6))
            # Very sparse strip for large-n tiers only
            _large_n_tiers17 = [t for t in valid_tiers_f17
                                 if len(f17_plot[f17_plot["degrader_tier"] == t]) > 30]
            if _large_n_tiers17:
                _f17_large = f17_plot[f17_plot["degrader_tier"].isin(_large_n_tiers17)]
                sns.stripplot(data=_f17_large, x="degrader_tier", y="Active_Site_RMSD",
                              order=valid_tiers_f17, color="black", alpha=0.04,
                              size=1.5, jitter=True, ax=ax, zorder=2)
            """
            For tiers with very few data points (≤30), draw individual dots prominently
            so every value is visible — avoids them being lost inside the violin body
            """
            for _i17, _tier17 in enumerate(valid_tiers_f17):
                _sub17 = f17_plot.loc[f17_plot["degrader_tier"] == _tier17, "Active_Site_RMSD"].dropna()
                if len(_sub17) <= 30:
                    _col17 = TIER_PALETTE.get(_tier17, "#999")
                    _jit17 = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(_sub17))
                    ax.scatter(_i17 + _jit17, _sub17.values,
                               color=_col17, edgecolors="black", s=55,
                               linewidths=0.8, alpha=0.88, zorder=7)
            # Outlier scatter: actual RMSD values above 3 Å (unclipped)
            for _ti, tier in enumerate(valid_tiers_f17):
                sub_f6 = f17_df[f17_df["degrader_tier"] == tier]
                _out = sub_f6[sub_f6["Active_Site_RMSD"] > _violin_cap]["Active_Site_RMSD"].values
                if len(_out) == 0:
                    continue
                _jit_f6 = np.random.default_rng(17 + _ti).uniform(-0.22, 0.22, size=len(_out))
                ax.scatter(_ti + _jit_f6, _out,
                           color=TIER_PALETTE.get(tier, "#999"), s=38, alpha=0.65,
                           marker="o", edgecolors="white", linewidths=0.4, zorder=4)
                if len(_out) >= 1:
                    _lbl_col = TIER_PALETTE.get(tier, "#999")
                    ax.text(_ti, y_ceil + 0.04, f"+{len(_out)}",
                            ha="center", va="bottom", fontsize=9.5,
                            color=_lbl_col, fontweight="bold",
                            zorder=30, clip_on=False,
                            bbox=dict(boxstyle="round,pad=0.18", fc="white",
                                      ec=_lbl_col, alpha=1.0, linewidth=1.0))

            # Capped axis, fixed ticks, and dashed break marker at 3 Å
            ax.axhline(3.0, color="#888", linestyle=":", linewidth=0.9, alpha=0.7)
            ax.set_ylim(0, y_ceil)
            _y_ticks_f6 = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
            _y_extra_f6 = [t for t in np.arange(3.5, y_ceil, 0.5) if t <= y_ceil - 0.1]
            _y_ticks_f6 += _y_extra_f6
            ax.set_yticks(_y_ticks_f6)
            ax.set_yticklabels([f"{t:.1f}" for t in _y_ticks_f6])
            ax.tick_params(axis="y", labelsize=7)
            ax.set_xlabel("Degrader Tier", fontsize=11)
            ax.set_ylabel("Active Site RMSD (Å)  — lower = better", fontsize=10)
            ax.text(0.99, 0.99, " ",
                    ha="right", va="top", fontsize=7, color="#888", style="italic",
                    transform=ax.transAxes)
            # Build two-line x-tick labels: tier name + n / median stats (using unclipped values)
            _xtick_labels17 = []
            for _tier17 in valid_tiers_f17:
                _sub17 = f17_df.loc[f17_df["degrader_tier"] == _tier17,
                                    "Active_Site_RMSD"].dropna()
                if len(_sub17) > 0:
                    _med17 = float(_sub17.median())
                    _xtick_labels17.append(
                        f"{_tier17}\nn={len(_sub17):,}  med={_med17:.2f}Å"
                    )
                else:
                    _xtick_labels17.append(_tier17)
            ax.set_xticks(range(len(valid_tiers_f17)))
            ax.set_xticklabels(_xtick_labels17, rotation=40, ha="right", fontsize=7.5)
            for _tick17, _tier17 in zip(ax.get_xticklabels(), valid_tiers_f17):
                _tick17.set_color(TIER_PALETTE.get(_tier17, "black"))
                _tick17.set_fontweight("bold")

            # Per-zone count labels (90°-rotated) inside each RMSD zone band
            _re, _ra, _rd = CFG.RMSD_BAND_EXCELLENT, CFG.RMSD_BAND_ACCEPTABLE, CFG.RMSD_BAND_DIVERGED
            _cbc = CFG.CONF_BAND_COLOURS
            _zone_count_defs = [
                (0.0, _re,    _cbc["high"],       f"<{_re:.0f} Å"),
                (_re, _ra,    _cbc["acceptable"], f"{_re:.0f}–{_ra:.0f} Å"),
                (_ra, _rd,    _cbc["below"],      f"{_ra:.0f}–{_rd:.0f} Å"),
                (_rd, np.inf, "#990000",          f">{_rd:.0f} Å"),
            ]
            for _zlo, _zhi, _zcol, _zlbl in _zone_count_defs:
                _zmask = (f17_df["Active_Site_RMSD"] >= _zlo) & (f17_df["Active_Site_RMSD"] < _zhi)
                _zn = int(_zmask.sum())
                _zy = min(_zhi, y_ceil) * 0.98
                ax.text(len(valid_tiers_f17) - 0.3, _zy,
                        f"n={_zn:,}  {_zlbl}",
                        ha="center", va="top", fontsize=6.5, color=_zcol, fontweight="bold",
                        rotation=90, zorder=25, clip_on=False,
                        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=_zcol,
                                  alpha=1.0, linewidth=0.5))

            # Zone band legend (top-left) for the shaded quality bands
            from matplotlib.patches import Patch as _P6z
            _zone_legend = [
                _P6z(facecolor="#D4EDDA", alpha=0.55, edgecolor="none", label=f"Excellent  (< 1.0 Å)  {pct_green:.0f}%"),
                _P6z(facecolor="#FFF3CD", alpha=0.55, edgecolor="none", label=f"Acceptable  (1.0–2.0 Å)  {pct_yellow:.0f}%"),
                _P6z(facecolor="#F8D7DA", alpha=0.55, edgecolor="none", label=f"Diverged  (> 2.0 Å)  {pct_red:.0f}%"),
            ]
            _tier_handles, _tier_labels = ax.get_legend_handles_labels()
            _combined_handles = _zone_legend + _tier_handles
            _combined_labels  = [h.get_label() for h in _zone_legend] + _tier_labels
            _mdh17 = _md_ready_stars_cat(ax, f17_df, valid_tiers_f17, "Active_Site_RMSD")
            if _mdh17 is not None:
                _combined_handles = _combined_handles + [_mdh17]
                _combined_labels  = _combined_labels + [_mdh17.get_label()]
            ax.legend(_combined_handles, _combined_labels,
                      loc="upper left", fontsize=8, framealpha=0.92)
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_07_ActiveSite_RMSD_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Figure 07 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 08: Feature Correlation Matrix (Spearman ρ) ---
    _corr_col_labels = {
        "Boltz_Model_Confidence":      "Conf",
        "iptm":                        "ipTM",
        "Binding_Probability":         "BindP",
        "Interaction_Density_Norm":    "IntDen",
        "Chemical_Affinity_Score":     "ChemAff",
        "tier_numeric":                "Tier",
        "Pareto_Rank":                 "Pareto",
        "SN2_Attack_Angle":            "SN2°",
        "identity_pct":                "SeqID%",
        "SN2_Trajectory_Deviation_A":  "TrajDev",
        "soft_catalytic_score":        "SoftCat",
    }
    # Map tier_numeric alias — try every plausible name written by different pipeline versions
    for _alias in ["tier_val", "tier_value", "Tier_Score", "tier_score",
                   "Tier_Numeric", "degrader_tier_numeric", "Degrader_Tier_Numeric"]:
        if _alias in df.columns and "tier_numeric" not in df.columns:
            df["tier_numeric"] = pd.to_numeric(df[_alias], errors="coerce")
    # Last-resort: derive numeric rank from the categorical degrader_tier string
    if "tier_numeric" not in df.columns and "degrader_tier" in df.columns:
        _tier_rank_map = {t: i for i, t in enumerate(TIER_ORDER_LOGIC)}
        df["tier_numeric"] = df["degrader_tier"].map(_tier_rank_map)
    corr_cols = [c for c in _corr_col_labels if c in df.columns]
    if len(corr_cols) >= 2:
        from scipy import stats as _sp_stats

        """
        Pairwise-complete Spearman (pandas corr with min_periods) — no median
        imputation. Median imputation inflates tied ranks and biases ρ for features
        with systematic missingness (features observed only for a subset of complexes); pairwise
        deletion uses each pair's jointly observed rows instead.
        """
        sub_f3 = df[corr_cols].apply(lambda s: pd.to_numeric(s, errors="coerce"))
        _min_n   = 20   # min jointly-observed rows for a correlation / p-value (shared floor)
        corr_raw = sub_f3.corr(method="spearman", min_periods=_min_n)
        _n_f3 = int(len(sub_f3))   # total complexes available (pairwise n ≤ this)

        # Per-pair p-values on jointly observed rows (nan_policy='omit'). The p-value
        # matrix honours the SAME _min_n sample-size floor as corr_raw: pairs below it
        # stay NaN so Benjamini–Hochberg excludes them rather than counting them as
        # p=1.0 (which would inflate the hypothesis count m and over-correct the real
        # tests). p_mat is therefore NaN-initialised, not ones-initialised.
        p_mat = np.full((len(corr_cols), len(corr_cols)), np.nan)
        for _ii, _c1 in enumerate(corr_cols):
            for _jj, _c2 in enumerate(corr_cols):
                if _ii != _jj:
                    _joint_n = int((sub_f3[_c1].notna() & sub_f3[_c2].notna()).sum())
                    if _joint_n < _min_n:
                        continue   # underpowered → leave NaN (excluded from BH)
                    try:
                        _, _pv = _sp_stats.spearmanr(sub_f3[_c1], sub_f3[_c2],
                                                     nan_policy="omit")
                        p_mat[_ii, _jj] = _pv if np.isfinite(_pv) else np.nan
                    except Exception:
                        p_mat[_ii, _jj] = np.nan
        p_df_f3 = pd.DataFrame(p_mat, index=corr_raw.index, columns=corr_raw.columns)

        # Benjamini-Hochberg FDR on the UNIQUE upper-triangular pairs only
        # (N(N-1)/2 tests, not N(N-1)); the adjusted matrix is then mirrored.
        _n_feat_f7 = len(corr_cols)
        '''
        Only correct over ACTUALLY-TESTED pairs: pairs masked out by the
        min_periods sample-size threshold are NaN in p_mat and must be excluded
        from the Benjamini–Hochberg hypothesis count, else m is inflated and the
        adjusted p-values become overly conservative. Masked pairs keep p_adj = 1.
        '''
        _pv_flat_f7, _pv_idx_f7 = [], []
        for _ii7 in range(_n_feat_f7):
            for _jj7 in range(_ii7 + 1, _n_feat_f7):
                _pv7 = p_mat[_ii7, _jj7]
                if _pv7 is None or (isinstance(_pv7, float) and np.isnan(_pv7)):
                    continue
                _pv_flat_f7.append(_pv7)
                _pv_idx_f7.append((_ii7, _jj7))
        _pv_arr = np.array(_pv_flat_f7)
        _sort_idx7 = np.argsort(_pv_arr)
        _m7 = len(_pv_arr)
        _adj7 = np.ones(_m7)
        _prev_adj = 1.0
        for _rank7 in range(_m7 - 1, -1, -1):
            _si7 = _sort_idx7[_rank7]
            _adj7[_si7] = min(_prev_adj, _pv_arr[_si7] * _m7 / (_rank7 + 1))
            _prev_adj = _adj7[_si7]
        _adj7 = np.clip(_adj7, 0, 1)
        p_adj_mat_f7 = np.ones((_n_feat_f7, _n_feat_f7))
        for _k7, (_ii7, _jj7) in enumerate(_pv_idx_f7):
            p_adj_mat_f7[_ii7, _jj7] = _adj7[_k7]
            p_adj_mat_f7[_jj7, _ii7] = _adj7[_k7]   # mirror to keep matrix symmetric
        p_adj_df_f7 = pd.DataFrame(p_adj_mat_f7, index=corr_raw.index, columns=corr_raw.columns)

        """
        Curated feature ORDER grouped by family (replaces the opaque clustering
        reorder). The degradation outcome (Tier) goes first, then each feature
        family — this makes the two regimes self-evident and lets us draw
        labelled separators and spotlight the Tier column below.
        """
        _fam_groups_f7 = [
            ("Outcome",              ["Tier"]),
            ("Catalytic geometry",   ["SN2°", "SoftCat", "TrajDev"]),
            ("Affinity / seq / rank", ["ChemAff", "SeqID%", "Pareto"]),
            ("Binding & confidence", ["BindP", "IntDen", "Conf", "ipTM"]),
        ]
        _cur_lbls_f7 = [_corr_col_labels[c] for c in corr_cols]
        _present_f7  = set(_cur_lbls_f7)
        _display_names = [lbl for _, lbls in _fam_groups_f7 for lbl in lbls
                          if lbl in _present_f7]
        _display_names += [l for l in _cur_lbls_f7 if l not in _display_names]
        _ord = [_cur_lbls_f7.index(lbl) for lbl in _display_names]
        corr_raw  = corr_raw.iloc[_ord, _ord]
        p_df_f3   = p_df_f3.iloc[_ord, _ord]
        p_adj_df_f7 = p_adj_df_f7.iloc[_ord, _ord]
        corr_raw.index   = _display_names
        corr_raw.columns = _display_names
        p_df_f3.index    = _display_names
        p_df_f3.columns  = _display_names
        p_adj_df_f7.index   = _display_names
        p_adj_df_f7.columns = _display_names

        # Build annotation: ρ value + significance stars
        def _sig_stars(p):
            return "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))

        annot_f3 = pd.DataFrame("", index=corr_raw.index, columns=corr_raw.columns)
        for _ri in corr_raw.index:
            for _ci in corr_raw.columns:
                _rho  = corr_raw.loc[_ri, _ci]
                _pv   = p_adj_df_f7.loc[_ri, _ci]
                if pd.isna(_rho):
                    annot_f3.loc[_ri, _ci] = "—"
                    continue
                _star = _sig_stars(_pv)
                annot_f3.loc[_ri, _ci] = f"{_rho:.2f}\n{_star}" if _star else f"{_rho:.2f}"

        # k=1 excludes diagonal from mask so self-correlation (1.00) is shown
        mask_f3 = np.triu(np.ones_like(corr_raw, dtype=bool), k=1)
        for _dn in corr_raw.index:
            annot_f3.loc[_dn, _dn] = "1.00"
        """
        Dendrogram removed — features are reordered by the curated family grouping
        (_fam_groups_f7 / _ord above), not by clustering. Single-axes heatmap only.
        """
        _, ax_f3 = plt.subplots(figsize=(12, 11))
        sns.heatmap(corr_raw, mask=mask_f3, annot=annot_f3, fmt="",
                    cmap="coolwarm", vmin=-1, vmax=1, center=0,
                    square=True, linewidths=0.5,
                    annot_kws={"size": 8.5},
                    cbar_kws={"label": "Spearman ρ  (−1 = perfect negative, +1 = perfect positive)",
                              "shrink": 0.75},
                    ax=ax_f3)
        # Significance / method key — on top, above the matrix.
        ax_f3.text(0.5, 1.03,
                   f"★ p<0.05   ★★ p<0.01   ★★★ p<0.001  (BH FDR, unique pairs)   |   Features grouped by family   |   pairwise-complete Spearman, n ≤ {_n_f3:,}",
                   transform=ax_f3.transAxes, ha="center", va="bottom",
                   fontsize=8, color="#555", style="italic",
                   bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#ccc",
                             alpha=0.88, linewidth=0.7))
        # Labels are now short (≤7 chars) — no truncation needed
        ax_f3.set_xticklabels([l.get_text() for l in ax_f3.get_xticklabels()],
                              fontsize=9, rotation=45, ha="right")
        ax_f3.set_yticklabels([l.get_text() for l in ax_f3.get_yticklabels()],
                              fontsize=9, rotation=90, va="center", ha="right")

        # Colour significance stars: white on dark cells (|ρ|>0.5), dark on light cells
        for _txt_art in ax_f3.texts:
            _raw = _txt_art.get_text()
            _lines = _raw.split("\n")
            if len(_lines) == 2 and _lines[1] and all(c == "*" for c in _lines[1]):
                try:
                    _rho_val = float(_lines[0])
                    _txt_art.set_color("white" if abs(_rho_val) > 0.5 else "#1a1a1a")
                except ValueError:
                    pass

        # ── Self-explanatory overlays ─────────────────────────────────────────
        from matplotlib.patches import Rectangle as _Rect7
        _N7 = len(_display_names)
        _fam_lookup_f7 = {l: g for g, ls in _fam_groups_f7 for l in ls}
        _fams_seq_f7   = [_fam_lookup_f7.get(l, "Other") for l in _display_names]

        """
        1) Family separator lines at each group boundary — clipped to the lower
        (data) triangle so they never run through the empty upper-right area
        that holds the key + take-away note.
        """
        for _i7 in range(1, _N7):
            if _fams_seq_f7[_i7] != _fams_seq_f7[_i7 - 1]:
                ax_f3.plot([0, _i7], [_i7, _i7], color="#222", linewidth=1.7, zorder=5)
                ax_f3.plot([_i7, _i7], [_i7, _N7], color="#222", linewidth=1.7, zorder=5)

        # 2) Family bracket labels under the matrix (one per group).
        for _gname, _glbls in _fam_groups_f7:
            _gi = [i for i, l in enumerate(_display_names) if l in _glbls]
            if not _gi:
                continue
            ax_f3.annotate(_gname, xy=((min(_gi) + max(_gi) + 1) / 2.0, _N7 + 1.15),
                           xycoords="data", ha="center", va="top", fontsize=8.5,
                           fontweight="bold", color="#333", annotation_clip=False)

        """
        3) Spotlight the Tier column — every feature's correlation with the
        degradation outcome is the question the figure answers.
        """
        if "Tier" in _display_names:
            _ti7 = _display_names.index("Tier")
            ax_f3.add_patch(_Rect7((_ti7, _ti7), 1, _N7 - _ti7, fill=False,
                            edgecolor="#111111", linewidth=2.4, zorder=6))

        """
        4) Take-away note — in the empty upper-right triangle (separators above
        are clipped to the data triangle, so the note sits clear of them).
        """
        ax_f3.text(
            0.605, 0.90,
            "• Degradation Tier tracks CATALYTIC GEOMETRY\n"
            "   (SN2°↑, SoftCat↑, TrajDev↓;  |ρ| ≈ 0.7–0.8)\n"
            "• Tier is ~independent of BINDING & CONFIDENCE\n"
            "   (Conf, ipTM, BindP, IntDen;  |ρ| < 0.1)\n"
            "→ binding well ≠ degrading well",
            transform=ax_f3.transAxes, ha="left", va="top", fontsize=9.8,
            color="#1a1a1a", linespacing=1.4,
            bbox=dict(boxstyle="round,pad=0.6", fc="#FFF8E7", ec="#C9A227",
                      alpha=0.95, linewidth=1.2))

        plt.tight_layout()
        plt.savefig(out_dir / "Figure_08_Feature_Correlations.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # --- Figure 09: Tier Quality Summary — Cleveland Dot Plot (multi-metric) ---
    """
    Each metric occupies a horizontal row; coloured dots show each tier's normalised
    mean score.  A thin grey range line connects the min-to-max dot for that metric,
    immediately revealing which metrics discriminate tiers the most.
    Raw mean values annotate each dot for quantitative readability.
    """
    _f08_metric_map = {
        "Boltz_Model_Confidence": "AI Confidence", "mechanistic_score": "Mech. Score",
        "Binding_Probability": "Binding Prob.", "identity_pct": "Seq. Identity (%)",
        "SN2_Attack_Angle": "SN2 Angle (°)", "Active_Site_RMSD": "RMSD (Å, inv.)",
        "SN2_Trajectory_Deviation_A": "SN2 Traj. Dev. (Å)",
        "soft_catalytic_score": "Soft Catalytic Score",
    }
    _f08_invert = {"Active_Site_RMSD", "SN2_Trajectory_Deviation_A"}
    _f08_cols = [c for c in _f08_metric_map if c in df.columns and "degrader_tier" in df.columns]
    if len(_f08_cols) >= 2:
        try:
            _diag08_parts = []
            for _t08d in existing_tiers:
                _s08d = df[df["degrader_tier"] == _t08d]
                _v08d = pd.to_numeric(
                    _s08d["mechanistic_score"] if "mechanistic_score" in _s08d.columns
                    else pd.Series(dtype=float), errors="coerce").dropna()
                if len(_v08d):
                    _diag08_parts.append(f"{_t08d}:{float(_v08d.mean()):.2f}(n={len(_v08d)})")
                else:
                    _diag08_parts.append(f"{_t08d}:NA")
            pass  # Fig08 diag suppressed
            _f08_rows = []
            for tier in existing_tiers:
                sub = df[df["degrader_tier"] == tier]
                if sub.empty:
                    continue
                row = {"Tier": tier}
                for col in _f08_cols:
                    vals = pd.to_numeric(sub[col], errors="coerce").dropna()
                    # Mean point estimate — matches the per-dot annotation and the
                    # mean bootstrap CI (_bsci) drawn for the same tier/metric.
                    row[_f08_metric_map[col]] = float(vals.mean()) if len(vals) else np.nan
                _f08_rows.append(row)
            if _f08_rows:
                _f08_df   = pd.DataFrame(_f08_rows).set_index("Tier")
                _f08_norm = _f08_df.copy()
                for col in _f08_norm.columns:
                    orig_key = next((k for k, v in _f08_metric_map.items() if v == col), None)
                    mn, mx = float(_f08_norm[col].min()), float(_f08_norm[col].max())
                    _f08_norm[col] = (_f08_norm[col] - mn) / (mx - mn) if mx > mn else 0.5
                    if orig_key in _f08_invert:
                        _f08_norm[col] = 1 - _f08_norm[col]

                metrics_f08 = list(_f08_norm.columns)
                tiers_f08   = list(_f08_norm.index)
                n_met, n_tier = len(metrics_f08), len(tiers_f08)

                # ── Pre-pass: collect tied-metric flags (needed by both CI loop and drawing) ──
                _any_tied08     = False
                _tied_metrics08 = []   # [(metric_name, raw_value), ...]
                _met_tied_map08 = {}   # metric -> bool
                for _metric_pre in metrics_f08:
                    _nvs_pre = [float(_f08_norm.loc[t, _metric_pre]) for t in tiers_f08]
                    _nvs_pre = [v for v in _nvs_pre if not np.isnan(v)]
                    _is_tied_pre = (len(_nvs_pre) > 1 and
                                    (max(_nvs_pre) - min(_nvs_pre)) < 0.02)
                    _met_tied_map08[_metric_pre] = _is_tied_pre
                    if _is_tied_pre:
                        _any_tied08 = True

                # Bootstrap 95% CI helper
                def _bsci(vals, n_boot=1000):
                    if len(vals) < 2:
                        return float(vals[0]) if len(vals) else (0.0, 0.0)
                    rng_b = np.random.default_rng(99)
                    # Vectorised: one (n_boot, n) resample draw, mean along axis 1.
                    _samples = rng_b.choice(np.asarray(vals), size=(n_boot, len(vals)), replace=True)
                    _meds = np.mean(_samples, axis=1)
                    return float(np.percentile(_meds, 2.5)), float(np.percentile(_meds, 97.5))

                # Precompute per-tier bootstrap CIs for each metric (raw scale)
                _f08_ci = {}  # (tier, metric) -> (lo_norm, hi_norm)
                for col in _f08_cols:
                    disp = _f08_metric_map[col]
                    mn_raw = float(_f08_df[disp].min())
                    mx_raw = float(_f08_df[disp].max())
                    _inv = col in _f08_invert
                    for tier in tiers_f08:
                        sub_ci = df[df["degrader_tier"] == tier]
                        raw_vals = pd.to_numeric(sub_ci[col], errors="coerce").dropna().values
                        if len(raw_vals) >= 2:
                            lo_r, hi_r = _bsci(raw_vals)
                            if mx_raw > mn_raw:
                                lo_n = (lo_r - mn_raw) / (mx_raw - mn_raw)
                                hi_n = (hi_r - mn_raw) / (mx_raw - mn_raw)
                            else:
                                lo_n, hi_n = 0.5, 0.5
                            if _inv:
                                lo_n, hi_n = 1 - hi_n, 1 - lo_n
                            '''
                            Clip to the (widened) display range, not the [0,1] data
                            range: a tier whose mean sits at 0 or 1 has a CI that
                            legitimately overhangs the ends, and clamping at 0/1 would
                            truncate that uncertainty at the axis wall.
                            '''
                            _f08_ci[(tier, disp)] = (float(np.clip(lo_n, -0.15, 1.15)),
                                                     float(np.clip(hi_n, -0.15, 1.15)))
                        else:
                            _f08_ci[(tier, disp)] = None

                # 3 rows × 1 col — pair metrics two per row
                _met_pairs = [(metrics_f08[i*2], metrics_f08[i*2+1] if i*2+1 < n_met else None)
                              for i in range(3)]
                fig, axes_f08 = plt.subplots(3, 1, figsize=(max(11, n_tier * 1.6 + 3), 9),
                                             sharex=False)
                fig.subplots_adjust(hspace=0.35)

                for _row_idx, ax08 in enumerate(axes_f08):
                    _row_mets = [m for m in _met_pairs[_row_idx] if m is not None]
                    ax08.set_axisbelow(True)
                    ax08.xaxis.grid(True, color="#EBEBEB", linewidth=0.7, linestyle="-", zorder=0)
                    ax08.yaxis.grid(False)

                    for _rmi, metric in enumerate(_row_mets):
                        _base_y = _rmi * 1.0   # row spacing within subplot
                        norm_vals  = _f08_norm[metric].values
                        valid_mask = ~np.isnan(norm_vals)
                        _is_tied08 = _met_tied_map08.get(metric, False)

                        _dot_data08 = []
                        for _ti, tier in enumerate(tiers_f08):
                            nv = float(_f08_norm.loc[tier, metric])
                            rv = float(_f08_df.loc[tier, metric])
                            if np.isnan(nv):
                                continue
                            _dot_data08.append((nv, rv, TIER_PALETTE.get(tier, "#999"), tier))
                        _dot_data08.sort(key=lambda d: d[0])
                        _nd08 = len(_dot_data08)
                        _all_nvs08 = [d[0] for d in _dot_data08]
                        _is_tied08 = bool(_nd08 > 1 and (max(_all_nvs08) - min(_all_nvs08)) < 0.02)

                        if _is_tied08:
                            _raw_tied08 = float(_dot_data08[0][1]) if _dot_data08 else float("nan")
                            if not any(m == metric for m, _ in _tied_metrics08):
                                _tied_metrics08.append((metric, _raw_tied08))

                        # Range line
                        if valid_mask.sum() >= 2:
                            if _is_tied08:
                                ax08.hlines(_base_y, 0.0, 1.0, colors="#D4A000",
                                            linewidth=1.0, linestyle="--", alpha=0.55, zorder=1)
                            else:
                                ax08.hlines(_base_y, float(np.nanmin(norm_vals)),
                                            float(np.nanmax(norm_vals)),
                                            colors="#CCCCCC", linewidth=2.5, zorder=1)

                        # Bootstrap CI bars drawn at base_y (before nudge known) — always visible
                        for nv, rv, col_f08, tier in _dot_data08:
                            _ci08 = _f08_ci.get((tier, metric))
                            if _ci08:
                                _clo = min(_ci08[0], nv - 0.013)
                                _chi = max(_ci08[1], nv + 0.013)
                                ax08.hlines(_base_y, _clo, _chi,
                                            colors=col_f08, linewidth=7.0, alpha=0.25,
                                            zorder=2)

                        # Y-nudge for overlapping dots
                        _y_nudge08 = [0.0] * _nd08
                        _gi08 = 0
                        while _gi08 < _nd08:
                            _gj08 = _gi08
                            while (_gj08 + 1 < _nd08 and
                                   abs(_dot_data08[_gj08+1][0] - _dot_data08[_gi08][0]) < 0.05):
                                _gj08 += 1
                            _gsize = _gj08 - _gi08 + 1
                            if _gsize > 1:
                                _offs = np.linspace(-0.08*(_gsize-1)/2, 0.08*(_gsize-1)/2, _gsize)
                                for _k, _off in enumerate(_offs):
                                    _y_nudge08[_gi08+_k] = float(_off)
                            _gi08 = _gj08 + 1

                        _dot_s08   = 110 if _is_tied08 else 190
                        _dot_alp08 = 0.38 if _is_tied08 else 1.0
                        for _di08, (nv, rv, col_f08, tier) in enumerate(_dot_data08):
                            ax08.scatter(nv, _base_y + _y_nudge08[_di08],
                                         color=col_f08, s=_dot_s08, zorder=4,
                                         alpha=_dot_alp08,
                                         edgecolors="black", linewidths=0.7)

                        # Value labels
                        if not _is_tied08:
                            _ylev08 = [0] * _nd08
                            for _j08 in range(1, _nd08):
                                if abs(_dot_data08[_j08][0] - _dot_data08[_j08-1][0]) < 0.10:
                                    _ylev08[_j08] = 1 - _ylev08[_j08-1]
                            _y_tops08 = [0.17, 0.30]
                            for _di08, (nv, rv, col_f08, tier) in enumerate(_dot_data08):
                                _nv_c = max(0.04, min(1.06, nv))
                                _y_pos08 = _base_y + _y_nudge08[_di08] + _y_tops08[_ylev08[_di08]]
                                ax08.text(_nv_c, _y_pos08, f"{rv:.2f}",
                                          ha="center", va="bottom", fontsize=6.5,
                                          color=col_f08, fontweight="bold", zorder=5,
                                          rotation=90,
                                          bbox=dict(boxstyle="round,pad=0.07", fc="white",
                                                    ec="none", alpha=0.75))

                    ax08.set_yticks([i * 1.0 for i in range(len(_row_mets))])
                    ax08.set_yticklabels(_row_mets, fontsize=10.5)
                    ax08.set_ylim(-0.5, len(_row_mets) - 0.5 + 0.55)
                    ax08.set_xlim(-0.18, 1.20)
                    ax08.spines[["top", "right"]].set_visible(False)

                for _ax08x in axes_f08:
                    _ax08x.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
                    _ax08x.set_xticklabels(["0\n(worst)", "0.25", "0.50", "0.75", "1.0\n(best)"],
                                           fontsize=8.5)
                axes_f08[-1].set_xlabel(
                    "Normalised Mean Score per Tier  (0 = worst, 1 = best within metric)",
                    fontsize=11)

                from matplotlib.lines import Line2D as _L08
                _leg08 = [_L08([0],[0], marker="o", color="w",
                               markerfacecolor=TIER_PALETTE.get(t,"#999"),
                               markeredgecolor="black", markeredgewidth=0.8,
                               markersize=10, label=t) for t in tiers_f08]
                _leg08.append(_L08([0],[0], color="none", linewidth=0,
                                   label="★ RMSD metrics: inverted — lower RMSD = better score"))
                _leg08.append(_L08([0],[0], color="#999", linewidth=5.5, alpha=0.3,
                                   label="95% bootstrap CI on mean (1 000 resamples)"))
                if _any_tied08:
                    _leg08.append(_L08([0],[0], color="#D4A000", linewidth=1.5,
                                       linestyle="--",
                                       label="⚡ Amber dashed = metric tied across all tiers"))
                    for _tm_name, _tm_val in _tied_metrics08:
                        _leg08.append(_L08([0],[0], color="none", linewidth=0,
                                           label=f"   {_tm_name}: all tiers = {_tm_val:.2f}"))
                axes_f08[0].legend(handles=_leg08, loc="upper left",
                                   bbox_to_anchor=(0.0, 1.0),
                                   ncol=max(1, (len(_leg08) + 1) // 2),
                                   fontsize=7.5, framealpha=0.92, fancybox=True)
                plt.tight_layout()
                plt.savefig(out_dir / "Figure_09_Tier_Quality_DotPlot.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
                plt.close()
        except Exception as e:
            reporter.log(f"  ! Figure 09 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 10: Mechanistic State Cross-Tab (Halide Stabilisation × Carboxylate Clamp) ---
    # Column alias detection for both mechanistic columns
    for _hs_alias in [
        "Has_Halide_Stabilisation",
        "Halide_Stabilisation", "halide_stabilisation", "halide_stabilization",
        "HalideStabilisation", "halide_shield", "Halide_Shield",
        "Aromatic_Shield", "aromatic_shield", "halide_interaction", "Halide_Interaction",
        "has_halide_shield", "Has_Halide_Shield", "halide_clamp", "Halide_Clamp",
    ]:
        if _hs_alias in df.columns and "Halide_Stabilisation" not in df.columns:
            df["Halide_Stabilisation"] = df[_hs_alias]
    # Partial-match fallback: any column with 'halide' or 'stabilisa'/'shield'
    if "Halide_Stabilisation" not in df.columns:
        _hs_fb = next((c for c in df.columns
                       if any(k in c.lower() for k in ("halide", "stabilisa", "shield"))), None)
        if _hs_fb:
            df["Halide_Stabilisation"] = df[_hs_fb]

    for _cc_alias in [
        "Carboxylate_Clamp", "carboxylate_clamp", "CarboxylateClamp",
        "carboxylate_clamping", "Carboxylate_Clamping",
        "Asp_Clamp", "asp_clamp", "ASP_Clamp", "asp_clamping",
        "carboxylate_contact", "Carboxylate_Contact",
        "has_carboxylate_clamp", "Has_Carboxylate_Clamp",
    ]:
        if _cc_alias in df.columns and "Carboxylate_Clamp" not in df.columns:
            df["Carboxylate_Clamp"] = df[_cc_alias]
    # Partial-match fallback: any column with 'clamp' or 'carboxyl'
    if "Carboxylate_Clamp" not in df.columns:
        _cc_fb = next((c for c in df.columns
                       if any(k in c.lower() for k in ("clamp", "carboxyl", "asp_clamp"))), None)
        if _cc_fb:
            df["Carboxylate_Clamp"] = df[_cc_fb]
    if "Halide_Stabilisation" not in df.columns or "Carboxylate_Clamp" not in df.columns:
        _mech_missing = []
        if "Halide_Stabilisation" not in df.columns:
            _mech_missing.append("Halide_Stabilisation")
        if "Carboxylate_Clamp" not in df.columns:
            _mech_missing.append("Carboxylate_Clamp")
        _mech_candidates = [c for c in df.columns
                            if any(k in c.lower() for k in ("halide", "clamp", "carboxyl", "stabil"))]
        reporter.log(f"  ! Figure 10 skipped: columns not found: {_mech_missing}. "
                     f"Partial matches in CSV: {_mech_candidates}. "
                     f"Full column list: {list(df.columns)}")
    if "Halide_Stabilisation" in df.columns and "Carboxylate_Clamp" in df.columns and "degrader_tier" in df.columns:
        try:
            mech_df = df.copy()
            mech_df["HS"] = mech_df["Halide_Stabilisation"].map(
                lambda x: "Stabilised" if str(x).strip().lower() in ("true", "1", "yes") else "Unstabilised")
            mech_df["CC"] = mech_df["Carboxylate_Clamp"].map(
                lambda x: "Clamped" if str(x).strip().lower() in ("true", "1", "yes") else "Unclamped")
            mech_df["Mech_State"] = mech_df["HS"] + "\n" + mech_df["CC"]

            state_order = ["Stabilised\nClamped", "Stabilised\nUnclamped",
                           "Unstabilised\nClamped", "Unstabilised\nUnclamped"]
            cross = pd.crosstab(mech_df["degrader_tier"], mech_df["Mech_State"])
            cross = cross.reindex(index=[t for t in existing_tiers if t in cross.index],
                                  columns=[s for s in state_order if s in cross.columns], fill_value=0)

            # Build annotation: "N\n(X%)" where % is row-wise (per tier)
            row_totals = cross.sum(axis=1)
            annot_labels = cross.copy().astype(object)
            for row_idx in cross.index:
                for col_idx in cross.columns:
                    n = cross.loc[row_idx, col_idx]
                    pct = (n / row_totals[row_idx] * 100) if row_totals[row_idx] > 0 else 0
                    annot_labels.loc[row_idx, col_idx] = f"{n:,}\n({pct:.1f}%)"

            # Normalise rows to percentage for the heatmap display
            _cross_raw = cross.copy()
            cross_pct = cross.div(cross.sum(axis=1), axis=0).fillna(0) * 100

            # Chi-square on raw counts
            from scipy.stats import chi2_contingency as _chi2_f09
            try:
                _chi2_09, _p_09, _dof_09, _ = _chi2_f09(_cross_raw.values)
                _chi2_stat = (f"χ²({_dof_09}) = {_chi2_09:.1f},  "
                              f"p {'< 0.001' if _p_09 < 0.001 else f'= {_p_09:.3f}'}")
            except Exception:
                _chi2_stat = ""

            _, ax19 = plt.subplots(figsize=(11, 7))
            sns.heatmap(cross_pct, annot=annot_labels, fmt="", cmap="YlOrRd", linewidths=0.5,
                        annot_kws={"size": 10, "weight": "bold"},
                        cbar_kws={"label": "% of tier"}, ax=ax19)
            if _chi2_stat:
                ax19.text(0.99, 0.02, _chi2_stat, transform=ax19.transAxes,
                          ha="right", va="bottom", fontsize=8.5, style="italic",
                          color="#333",
                          bbox=dict(boxstyle="round,pad=0.25", fc="white",
                                    ec="#ccc", alpha=0.95, linewidth=0.6))
            ax19.set_xlabel("Mechanistic State", fontsize=11)
            ax19.set_ylabel("Degrader Tier", fontsize=11)
            ax19.set_xticklabels(ax19.get_xticklabels(), rotation=15, ha="right", fontsize=9)
            # Colour y-tick labels by tier (same palette as all other figures)
            for tick in ax19.get_yticklabels():
                tier_name = tick.get_text()
                tick.set_color(TIER_PALETTE.get(tier_name, "black"))
                tick.set_fontweight("bold")
                tick.set_rotation(0)
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_10_Mech_State_CrossTab.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Figure 10 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 11: Mechanistic Score — Mean±CI dot plot ---
    if "mechanistic_score" in df.columns and "degrader_tier" in df.columns:
        from scipy import stats as _scipy_stats
        f13_data = df.dropna(subset=["mechanistic_score"])
        fig, ax = plt.subplots(figsize=(13, 8))
        _mfs, _mfm = CFG.MECH_FP_BAND_STRONG, CFG.MECH_FP_BAND_MODERATE   # CFG single source
        _cbc = CFG.CONF_BAND_COLOURS
        ax.axhspan(_mfs, 1.01, alpha=0.10, color=_cbc["high"], zorder=0, label="_nolegend_")
        ax.axhspan(_mfm, _mfs, alpha=0.09, color=_cbc["acceptable"], zorder=0, label="_nolegend_")
        ax.axhspan(0.00, _mfm, alpha=0.08, color=_cbc["below"], zorder=0, label="_nolegend_")
        ax.text(0.99, 0.942, f"Strong zone  (≥{_mfs:.2f})", color="#007A50",
                fontsize=8, ha="right", va="center", fontweight="bold",
                style="italic", transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.15", fc="#D6F5EB", ec=CFG.CONF_BAND_COLOURS["high"],
                          alpha=0.85, linewidth=0.6))
        ax.text(0.99, 0.798, f"Moderate zone  ({_mfm:.2f}–{_mfs:.2f})", color="#8A6000",
                fontsize=8, ha="right", va="center", fontweight="bold",
                style="italic", transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.15", fc="#FFF3CC", ec=CFG.CONF_BAND_COLOURS["acceptable"],
                          alpha=0.85, linewidth=0.6))
        ax.text(0.99, 0.435, f"Weak zone  (<{_mfm:.2f})", color="#A03000",
                fontsize=8, ha="right", va="center", fontweight="bold",
                style="italic", transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.15", fc="#FDECEA", ec=CFG.CONF_BAND_COLOURS["below"],
                          alpha=0.85, linewidth=0.6))
        ax.axhline(y=_mfs, color=_cbc["high"], linestyle="--", linewidth=1.2, alpha=0.7)
        ax.axhline(y=_mfm, color=_cbc["below"], linestyle=":", linewidth=1.0, alpha=0.7)
        sns.stripplot(data=f13_data, x="degrader_tier", y="mechanistic_score",
                      order=existing_tiers, palette=TIER_PALETTE, alpha=0.15, size=3.0,
                      jitter=0.28, ax=ax, zorder=1)
        _f10_stats = {}
        for _t10 in existing_tiers:
            _v10 = f13_data.loc[f13_data["degrader_tier"] == _t10,
                                 "mechanistic_score"].dropna()
            if len(_v10) >= 2:
                _f10_stats[_t10] = {
                    "mean": float(_v10.mean()),
                    "ci95": float(_scipy_stats.sem(_v10) * _scipy_stats.t.ppf(0.975, len(_v10) - 1)),
                    "med":  float(_v10.median()),
                    "n":    len(_v10),
                    "vals": _v10,
                }
        for i, tier in enumerate(existing_tiers):
            if tier not in _f10_stats:
                continue
            _st = _f10_stats[tier]
            mean_v, ci95, med_v, vals = _st["mean"], _st["ci95"], _st["med"], _st["vals"]
            ax.plot([i, i], [mean_v - ci95, mean_v + ci95], color="black", linewidth=2.2, zorder=4)
            ax.scatter([i], [mean_v], color=TIER_PALETTE.get(tier, "#999"), s=180,
                       zorder=5, edgecolors="black", linewidths=1.2)
            pct_strong   = (vals >= _mfs).sum() / len(vals) * 100
            pct_moderate = ((vals >= _mfm) & (vals < _mfs)).sum() / len(vals) * 100
            pct_weak     = (vals < _mfm).sum() / len(vals) * 100
            if pct_strong > 0:
                ax.text(i, (_mfs + 1.0) / 2, f"{pct_strong:.0f}%", ha="center", va="center",
                        fontsize=7, color="#007A50", fontweight="bold", zorder=10,
                        bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["high"],
                                  alpha=0.75, linewidth=0.5))
            if pct_moderate > 0:
                ax.text(i, (_mfm + _mfs) / 2, f"{pct_moderate:.0f}%", ha="center", va="center",
                        fontsize=7, color="#8A6000", fontweight="bold", zorder=10,
                        bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["acceptable"],
                                  alpha=0.75, linewidth=0.5))
            if pct_weak > 0:
                ax.text(i, _mfm / 2, f"{pct_weak:.0f}%", ha="center", va="center",
                        fontsize=7, color="#A03000", fontweight="bold", zorder=10,
                        bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["below"],
                                  alpha=0.75, linewidth=0.5))
        ax.set_ylim(-0.02, 1.06)
        ax.set_yticks(sorted({0.0, 0.2, 0.4, 0.6, 0.8, 1.0, round(_mfm, 2), round(_mfs, 2)}))
        ax.set_xlabel("Degrader Tier", fontsize=11)
        ax.set_ylabel("Mechanistic Score  (anchor set + graded SN2 angle; config §5.1)", fontsize=10)
        ax.set_xticks(range(len(existing_tiers)))
        _xtlbl10 = []
        for _t10 in existing_tiers:
            if _t10 in _f10_stats:
                _s = _f10_stats[_t10]
                _lbl = f"{_t10}\nμ={_s['mean']:.3f}  ±{_s['ci95']:.3f}"
                if abs(_s["med"] - _s["mean"]) > 0.005:
                    _lbl += f"\nmed={_s['med']:.3f}"
                _lbl += f"\nn={_s['n']:,}"
            else:
                _lbl = _t10
            _xtlbl10.append(_lbl)
        ax.set_xticklabels(_xtlbl10, rotation=40, ha="right", fontsize=7.5)
        for _tick13x, _tier13x in zip(ax.get_xticklabels(), existing_tiers):
            _tick13x.set_color(TIER_PALETTE.get(_tier13x, "black"))
            _tick13x.set_fontweight("bold")
        ax.yaxis.grid(True, color="#DCDCDC", linewidth=0.6, alpha=0.70, zorder=0)
        ax.set_axisbelow(True)
        _means10 = [_f10_stats[t]["mean"] for t in existing_tiers if t in _f10_stats]
        if _means10 and (max(_means10) - min(_means10)) < 0.05:
            _mean10_all = float(np.mean(_means10))
            ax.axhline(y=_mean10_all, color="#D4A000", linewidth=1.5,
                       linestyle="--", alpha=0.65, zorder=2)
            ax.text(0.99, 0.52,
                    f"⚡ All tier means within {max(_means10)-min(_means10):.3f} "
                    f"— minimal inter-tier separation\n"
                    f"(mech-score cluster near {_mean10_all:.3f}; see also Fig 08 tied-metric row)",
                    transform=ax.transAxes, ha="right", va="center", fontsize=8,
                    color="#9A6B00", fontweight="bold", style="italic",
                    bbox=dict(boxstyle="round,pad=0.15", fc="#FFF8E7", ec="#D4A000",
                              alpha=0.90, linewidth=0.8), zorder=8)
        from matplotlib.lines import Line2D as _L10
        _lh10 = [
            _L10([0],[0], marker="o", color="none", markerfacecolor="#888888",
                 markeredgecolor="none", markersize=5, alpha=0.35,
                 label="Individual complex"),
            _L10([0],[0], marker="D", color="none", markerfacecolor="#888888",
                 markeredgecolor="black", markersize=8, markeredgewidth=1.2,
                 label="Tier mean (◆)"),
            _L10([0],[0], color="black", linewidth=2.0, label="95% confidence interval"),
        ]
        _mdh10 = _md_ready_stars_cat(ax, df, existing_tiers, "mechanistic_score")
        if _mdh10 is not None:
            _lh10.append(_mdh10)
        ax.legend(handles=_lh10, loc="lower left", fontsize=8.5,
                  framealpha=0.92, fancybox=True,
                  ncol=4).set_zorder(20)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_11_Mechanistic_Score_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # --- Figure 11b: Mechanistic Fingerprint — per-tier catalytic feature profile ---
    """
    Catalytic "fingerprint" per degrader tier rendered as a radar / spider chart: each
    spoke is a mechanistic feature, each filled polygon a tier, and the radius its tier-mean
    engagement (0–1). Boolean requirements → fraction engaged; the SN2 attack angle is
    normalised (angle / 180°, ideal back-side = 1.0); continuous scores → clipped mean.
    Built from the FAcDs mechanistic-score components (config §5), it shows which catalytic
    requirements each tier satisfies and how the profile degrades down the tiers —
    complementing the scalar Figure 11 (mech-score mean ± CI).
    """
    _fp_feats = [
        ("Carboxylate\nhead",          "head_is_carboxylate",      "bool"),
        ("Halide\nstabilisation",      "Has_Halide_Stabilisation", "bool"),
        ("Carboxylate\nclamp",         "Has_Carboxylate_Clamp",    "bool"),
        ("SN2 alignment\n(angle/180°)", "sn2_alignment_score",     "angle"),
        ("Mechanistic\nscore",         "mechanistic_score",        "score"),
    ]
    _fp_avail = [(lbl, col, kind) for (lbl, col, kind) in _fp_feats if col in df.columns]
    if _fp_avail and "degrader_tier" in df.columns:
        fig = None
        try:
            _fp_tiers = [t for t in existing_tiers if (df["degrader_tier"] == t).any()]
            _fp_rows, _fp_n = [], []
            for _ft in _fp_tiers:
                _fp_sub = df[df["degrader_tier"] == _ft]
                _fp_n.append(len(_fp_sub))
                _fp_row = []
                for _lbl, _col, _kind in _fp_avail:
                    if _kind == "bool":
                        _bs = _fp_sub[_col].astype(str).str.strip().str.lower().isin(
                            ["true", "1", "1.0", "yes"])
                        _fp_row.append(float(_bs.mean()) if len(_bs) else np.nan)
                    elif _kind == "angle":
                        _vs = (pd.to_numeric(_fp_sub[_col], errors="coerce") / 180.0).clip(0, 1)
                        _fp_row.append(float(_vs.mean()) if _vs.notna().any() else np.nan)
                    else:
                        _vs = pd.to_numeric(_fp_sub[_col], errors="coerce").clip(0, 1)
                        _fp_row.append(float(_vs.mean()) if _vs.notna().any() else np.nan)
                _fp_rows.append(_fp_row)
            _fpm = pd.DataFrame(_fp_rows, index=_fp_tiers,
                                columns=[l for l, _, _ in _fp_avail])
            '''
            Radar (spider) chart in the same visual style as 01_Radar_TopHits:
            one spoke per mechanistic feature, one filled polygon per tier, engagement
            plotted 0–1 directly so the red rings read as absolute engagement. Series
            capped at CFG.VIS_RADAR_MAX_HITS to keep the chart readable.
            '''
            _fp_tiers = _fp_tiers[:CFG.VIS_RADAR_MAX_HITS]
            _fp_labels = [l for l, _, _ in _fp_avail]
            _N_fp = len(_fp_labels)
            _fp_angles = [n / float(_N_fp) * 2 * np.pi for n in range(_N_fp)] + [0]
            fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
            # Conventional radar orientation: first spoke at top (12 o'clock), spokes laid
            # out clockwise so the polygon reads like a dial.
            ax.set_theta_offset(np.pi / 2)
            ax.set_theta_direction(-1)
            _fp_rings = [0.2, 0.4, 0.6, 0.8, 1.0]
            _fp_ring_ang = np.linspace(0, 2 * np.pi, 300)
            for _rl in _fp_rings:
                ax.plot(_fp_ring_ang, [_rl] * 300, color="#CC2222", linewidth=0.55,
                        alpha=0.45, zorder=1, linestyle="-", solid_capstyle="round")
            for _ti, _t in enumerate(_fp_tiers):
                _yv = _fpm.loc[_t].values.tolist()
                _vals = _yv + [_yv[0]]
                _c = TIER_PALETTE.get(_t, "#999999")
                ax.plot(_fp_angles, _vals, linewidth=2.2, linestyle="solid",
                        label=f"{_t}  (n={_fp_n[_ti]:,})", color=_c, zorder=3)
                ax.fill(_fp_angles, _vals, color=_c, alpha=0.08, zorder=2)
                ax.scatter(_fp_angles[:-1], _vals[:-1], color=_c, s=55, zorder=5,
                           edgecolors="white", linewidths=0.8)
            ax.set_xticks(_fp_angles[:-1])
            ax.set_xticklabels([])
            '''
            Spoke labels drawn horizontally hugging the outer ring (same proximity as
            01_Radar_TopHits); horizontal alignment follows each spoke's displayed
            direction (clockwise from top) so text sits beside its spoke without
            overrunning the plot or colliding with the radial tick ring.
            '''
            _lbl_r = 1.10 * 1.04          # hug the outer ring (ylim top = 1.10)
            for _ang_s, _lbl_s in zip(_fp_angles[:-1], _fp_labels):
                _disp = np.pi / 2 - _ang_s          # displayed angle (theta dir = -1, offset = +90°)
                _dx, _dy = np.cos(_disp), np.sin(_disp)
                _ha_s = "left" if _dx > 0.10 else "right" if _dx < -0.10 else "center"
                _va_s = "bottom" if _dy > 0.10 else "top" if _dy < -0.10 else "center"
                ax.text(_ang_s, _lbl_r, _lbl_s, ha=_ha_s, va=_va_s,
                        fontsize=8.5, fontweight="bold", linespacing=0.95)
            ax.set_yticks(_fp_rings)
            ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"],
                               color="#CC2222", size=7.5, fontweight="bold")
            # Radial tick labels parked in an empty wedge between spokes (no data line there).
            ax.set_rlabel_position(np.degrees(np.pi / float(_N_fp)))
            ax.set_ylim(0, 1.10)
            ax.yaxis.grid(False)
            ax.xaxis.grid(True, color="#CCCCCC", linewidth=0.6, alpha=0.7)
            _leg11b = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.18),
                                ncol=min(3, len(_fp_tiers)), fontsize=8,
                                framealpha=0.88, fancybox=True)
            _leg11b.set_zorder(20)
            for _lt, _t in zip(_leg11b.get_texts(), _fp_tiers):
                _lt.set_color(TIER_PALETTE.get(_t, "black"))
            plt.tight_layout(rect=[0, 0.15, 1, 1])
            plt.savefig(out_dir / "Figure_11b_Mechanistic_Fingerprint.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        except Exception as e:
            reporter.log(f"  ! Figure 11b skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig
        finally:
            if fig is not None:
                plt.close(fig)

    # --- Figure 12: SN2 Attack Angle — ECDF by Tier ---
    if "SN2_Attack_Angle" in df.columns and "degrader_tier" in df.columns:
        plot_df = df[df["SN2_Attack_Angle"] > 0].copy()
        valid_tiers_f14 = [t for t in existing_tiers if t in plot_df["degrader_tier"].values]
        fig, ax = plt.subplots(figsize=(12, 7))
        # Shade tier quality zones + count labels inside each zone band.
        # Zone lower bounds derive from CFG.TIER_ANGLE_MIN (single source); the upper
        # bound of each zone is the next-higher tier's minimum (top zone capped at 180°).
        _zt11   = [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]
        _zcol11 = ["#009E73", "#56B4E9", "#0072B2", "#CC79A7"]
        _zlo11  = [CFG.TIER_ANGLE_MIN[t] for t in _zt11]
        _zhi11  = [180.0] + _zlo11[:-1]
        _zone_defs11 = [(_zlo11[i], _zhi11[i], _zcol11[i]) for i in range(len(_zt11))]
        for lo, hi, col in _zone_defs11:
            ax.axvspan(lo, hi, alpha=0.07, color=col)
            ax.text((lo + hi) / 2, 1.025, f"≥{lo}°", ha="center", va="bottom",
                    fontsize=7.5, color=col, fontweight="bold", transform=ax.get_xaxis_transform())
            _zn11 = int(((plot_df["SN2_Attack_Angle"] >= lo) & (plot_df["SN2_Attack_Angle"] < hi)).sum())
            ax.text((lo + hi) / 2, 0.10, f"n={_zn11:,}",
                    ha="center", va="center", fontsize=6.0, color=col, fontweight="bold",
                    rotation=90, zorder=10,
                    transform=ax.get_xaxis_transform(),
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.68))
        # ECDF per tier — median value labels are collected here and placed after the
        # loop in alternating vertical bands so clustered medians never overlap.
        _med_pts = []
        for tier in valid_tiers_f14:
            tier_angles = np.sort(plot_df.loc[plot_df["degrader_tier"] == tier, "SN2_Attack_Angle"].values)
            ecdf_y = np.arange(1, len(tier_angles) + 1) / len(tier_angles)
            col = TIER_PALETTE.get(tier, "#999")
            ax.step(tier_angles, ecdf_y, where="post", color=col, linewidth=2.5,
                    label=f"{tier}  (n={len(tier_angles):,})")
            _n_ecdf = len(tier_angles)
            _eps_dkw = np.sqrt(np.log(2 / 0.05) / (2 * _n_ecdf))
            _ecdf_y_lo = np.clip(ecdf_y - _eps_dkw, 0, 1)
            _ecdf_y_hi = np.clip(ecdf_y + _eps_dkw, 0, 1)
            ax.fill_between(tier_angles, _ecdf_y_lo, _ecdf_y_hi,
                            color=TIER_PALETTE.get(tier, "#999"), alpha=0.15, zorder=1)
            med = float(np.median(tier_angles))
            med_y = float(np.interp(med, tier_angles, ecdf_y))
            ax.scatter([med], [med_y], color=col, s=70, zorder=6, edgecolors="black", linewidths=0.8)
            _med_pts.append((med, med_y, col))
        # Median value labels: sort by angle and place in four alternating vertical
        # bands (offsets in points) so clustered medians (e.g. 158/162/168°) never
        # collide; each label is joined to its dot by a thin leader arrow.
        _band14 = [40, -40, 24, -24]
        for _k, (_mx, _my, _mc) in enumerate(sorted(_med_pts, key=lambda z: z[0])):
            _oy = _band14[_k % len(_band14)]
            ax.annotate(f"{_mx:.0f}°", (_mx, _my), xytext=(0, _oy), textcoords="offset points",
                        fontsize=8, color=_mc, fontweight="bold", ha="center",
                        va="bottom" if _oy > 0 else "top", zorder=7,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=_mc, alpha=0.85, linewidth=0.6),
                        arrowprops=dict(arrowstyle="->", color=_mc, lw=0.7, shrinkA=0, shrinkB=3))
        # MD-ready overlay — the super-best MD-selected complexes (incl. those inside
        # Tier_1A) marked as gold stars on their own tier's ECDF curve.
        _mdsel14 = plot_df[plot_df.get("MD_Selected", pd.Series(False, index=plot_df.index))
                           .astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes"])]
        for _, _mr in _mdsel14.iterrows():
            _ma = float(_mr["SN2_Attack_Angle"]); _mt = _mr["degrader_tier"]
            _ta = np.sort(plot_df.loc[plot_df["degrader_tier"] == _mt, "SN2_Attack_Angle"].values)
            if len(_ta) == 0:
                continue
            _myv = float(np.interp(_ma, _ta, np.arange(1, len(_ta) + 1) / len(_ta)))
            ax.scatter([_ma], [_myv], marker="*", s=300, facecolor="#FFD400",
                       edgecolor="black", linewidths=1.1, zorder=9)
        if len(_mdsel14):
            ax.scatter([], [], marker="*", s=180, facecolor="#FFD400", edgecolor="black",
                       linewidths=1.0, label=f"MD-ready hits (n={len(_mdsel14)})")

        # Vertical threshold lines
        # Threshold lines sourced from CFG so they always match the tier gates
        # (e.g. Tier_1A = 170°, not a stale hardcoded literal).
        _ecdf_tiers = [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]
        _ecdf_angles = [CFG.TIER_ANGLE_MIN[t] for t in _ecdf_tiers]
        for angle, tier_key in zip(_ecdf_angles, _ecdf_tiers):
            ax.axvline(x=angle, color=TIER_PALETTE.get(tier_key, "#999"), linestyle="--", alpha=0.65, linewidth=1.0)
        ax.set_xlim(30, 185)
        ax.set_ylim(0, 1.09)
        ax.set_xticks([30, 45, 90, 130] + sorted(set(_ecdf_angles)) + [180])
        ax.set_yticks([i/10 for i in range(0, 11)])
        ax.set_yticklabels([f"{i*10}%" for i in range(0, 11)])
        ax.set_xlabel("SN2 Attack Angle (°)  — 180° = ideal linear nucleophilic back-attack", fontsize=11)
        ax.set_ylabel("Cumulative Fraction of Complexes in Tier", fontsize=11)
        # upper left — ECDF lines fan right so top-left is always clear
        from matplotlib.lines import Line2D as _L11
        _hdr11 = _L11([0],[0], color="none", linewidth=0, label="Tier (●=median)")
        _h11, _l11 = ax.get_legend_handles_labels()
        _leg14 = ax.legend(handles=[_hdr11] + _h11, labels=["Tier (●=median)"] + _l11,
                           loc="upper left", bbox_to_anchor=(0.01, 0.99),
                           fontsize=6, framealpha=0.92, fancybox=True,
                           ncol=max(2, (len(_h11) + 1) // 2))
        _leg14.set_zorder(20)
        # Thin horizontal grid lines at every 10% ECDF level for easy reading
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color="#DCDCDC", linewidth=0.6, linestyle="-", alpha=0.85, zorder=0)
        ax.xaxis.grid(True, color="#EBEBEB", linewidth=0.5, linestyle=":", alpha=0.7, zorder=0)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_12_SN2_Angle_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # --- Figure 13a: SN2 Mechanism Geometry Scatter ---
    if "Dist_Nucleophile" in df.columns and "SN2_Attack_Angle" in df.columns:
        _n_raw9 = len(df)
        plot_df12 = df[
            (df["Dist_Nucleophile"] > 0) &
            (df["Dist_Nucleophile"] < 5.0) &   # plot display window (Å), not a scientific gate
            (df["SN2_Attack_Angle"] >= 0)
        ].copy()
        n_excluded12 = _n_raw9 - len(plot_df12)
        plot_df12["Stabilised"] = plot_df12.get("halide_stabilisation_score",
                                                 pd.Series(0.0, index=plot_df12.index)) > 0.5
        plot_df12["Status"] = plot_df12["Stabilised"].map({True: "Aromatic Shield", False: "Unstabilised"})

        fig, ax = plt.subplots(figsize=(12, 8))
        markers12 = {"Aromatic Shield": "o", "Unstabilised": "X"}
        plot_df12_sorted = plot_df12.copy()
        plot_df12_sorted["_tier_ord"] = plot_df12_sorted["degrader_tier"].map(
            {t: i for i, t in enumerate(reversed(existing_tiers))})
        plot_df12_sorted = plot_df12_sorted.sort_values("_tier_ord")
        sns.scatterplot(data=plot_df12_sorted, x="Dist_Nucleophile", y="SN2_Attack_Angle",
                        hue="degrader_tier", hue_order=existing_tiers, style="Status",
                        markers=markers12, palette=TIER_PALETTE, alpha=0.60, s=55, ax=ax)
        if ax.get_legend():
            ax.get_legend().remove()

        # KDE density overlay per tier (large-n tiers only)
        try:
            from scipy.stats import gaussian_kde as _gkde12a
            _xi12 = np.linspace(0.1, 4.9, 100)
            _yi12 = np.linspace(5, 180, 100)
            _XX12, _YY12 = np.meshgrid(_xi12, _yi12)
            _grid12 = np.vstack([_XX12.ravel(), _YY12.ravel()])
            for _zi12, tier in enumerate(existing_tiers):
                _td12 = plot_df12[plot_df12["degrader_tier"] == tier]
                if len(_td12) < 30:
                    continue
                _col12 = TIER_PALETTE.get(tier, "#999")
                try:
                    _kde12 = _gkde12a(np.vstack([
                        _td12["Dist_Nucleophile"].values,
                        _td12["SN2_Attack_Angle"].values
                    ]))
                    _ZZ12 = _kde12(_grid12).reshape(_XX12.shape)
                    _ZZ12 = _ZZ12 / _ZZ12.max()
                    ax.contour(_XX12, _YY12, _ZZ12, levels=[0.25, 0.50, 0.90],
                               colors=[_col12], linewidths=[0.5, 1.0, 1.6],
                               alpha=0.70, zorder=_zi12 + 2)
                except Exception:
                    pass
        except ImportError:
            pass

        """
        Thresholds sourced DIRECTLY from CFG (not hardcoded) so the guide lines
        exactly match the production tier gates in 02_Production.
        """
        _nuc_cfg12 = CFG.TIER_NUC_DIST
        _ang_cfg12 = CFG.TIER_ANGLE_MIN
        ax.axvspan(0, _nuc_cfg12[CFG.TIER_TOP], alpha=0.06, color="#009E73", label="_nolegend_")
        ax.axhspan(_ang_cfg12[CFG.TIER_TOP], 180, alpha=0.06, color="#009E73", label="_nolegend_")

        for _idx12, tier_key in enumerate([CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]):
            dist = _nuc_cfg12[tier_key]
            col = TIER_PALETTE.get(tier_key, "#999")
            ax.axvline(x=dist, color=col, linestyle="--", alpha=0.65, linewidth=1.2)
            ax.text(dist, 186 if _idx12 % 2 == 0 else 181, f"{tier_key}\n≤{dist}Å", color=col,
                    fontsize=6.5, va="bottom", ha="center", alpha=0.90,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.75, ec=col,
                              linewidth=0.5))

        for tier_key in [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]:
            angle = _ang_cfg12[tier_key]
            col = TIER_PALETTE.get(tier_key, "#999")
            ax.axhline(y=angle, color=col, linestyle=":", alpha=0.65, linewidth=1.2, zorder=1)
            ax.text(0.02, angle, f"≥{angle}° ({tier_key})", color=col,
                    fontsize=6.5, ha="left", va="bottom", alpha=0.90, zorder=6,
                    transform=ax.get_yaxis_transform(),
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                              alpha=0.70, linewidth=0))

        ax.set_xlim(0, 5.0)
        ax.set_ylim(0, 195)
        ax.set_yticks(range(0, 196, 20))
        _nuc_lbl = CFG.REF_ACTIVE_SITE_MAP["Nuc"]
        ax.set_xlabel(f'{_nuc_lbl["res"]}{_nuc_lbl["id"]} Nucleophile → Carbon Distance (Å)  — shorter = closer to reaction geometry',
                      fontsize=10)
        ax.set_ylabel("SN2 Attack Angle (°)  — 180° = perfect linear back-attack", fontsize=10)
        ax.text(0.01, 0.01,
                f"n = {len(plot_df12):,} complexes shown  |  {n_excluded12:,} excluded (dist = 0 or dist ≥ 5 Å)",
                transform=ax.transAxes, fontsize=7.5, color="#666", style="italic", va="bottom",
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#ccc", alpha=0.82))
        handles12, labels12 = ax.get_legend_handles_labels()
        clean_labels12, clean_handles12 = [], []
        for h12, l12 in zip(handles12, labels12):
            if l12 in ("degrader_tier", "Status"):
                continue
            clean_handles12.append(h12)
            clean_labels12.append(l12)
        ax.legend(clean_handles12, clean_labels12,
                  loc="upper left", bbox_to_anchor=(0.01, 0.99),
                  fontsize=7, framealpha=0.92, fancybox=True,
                  ncol=max(2, (len(clean_handles12) + 1) // 2),
                  handlelength=1.2, handletextpad=0.4)
        ax.set_axisbelow(True)
        ax.grid(color="#E0E0E0", linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_13a_Mechanism_Geometry_Scatter.png",
                    dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # Figure 13b: PA companion (mechanistic scatter overlay)
    if _tt_has_imgs:
        _fig_13b_tt_mechanistic(df, _pa, _imgs, out_dir, reporter)

    # --- Shared Data Prep: Interaction Bond Columns (used by Fig 13 + 14) ---
    int_cols_raw = {
        "hydrogen_bond":      ("count_hydrogen_bond",      "counts_hydrogen_bond"),
        "salt_bridge":        ("count_salt_bridge",        "counts_salt_bridge"),
        "halogen_contact":    ("count_halogen_contact",    "counts_halogen_contact"),
        "fluorine_polar":     ("count_fluorine_polar",     "counts_fluorine_polar"),
        "fluorous_hydrophobic": ("count_fluorous_hydrophobic", "counts_fluorous_hydrophobic"),
        "hydrophobic":        ("count_hydrophobic",        "counts_hydrophobic"),
    }
    int_label_map = {
        "hydrogen_bond": "H-Bond",
        "salt_bridge": "Salt Bridge",
        "halogen_contact": "Halogen",
        "fluorine_polar": "F-Polar",
        "fluorous_hydrophobic": "F-Hydrophobic",
        "hydrophobic": "Hydrophobic",
    }
    present_int_cols = {}
    for key, (alias1, alias2) in int_cols_raw.items():
        if alias1 in df.columns:
            present_int_cols[key] = alias1
        elif alias2 in df.columns:
            present_int_cols[key] = alias2

    # ===========================================================================
    # Two-criteria tier logic + dead-end feasibility — one 3-panel figure written
    # to the catalytic folder (04). Panel a: Criterion A (active-site integrity,
    # active_site_residues_correct) gates Criterion B (catalytic constellation).
    # Panel b: Criterion-B ECDF separates the tiers (per-tier mean diamonds).
    # Panel c: the SN2 dead-end chemistry gate (scissile C–F BDE × backside
    # occlusion) for the ligands that actually carry the penalty, plus the FA/DFA
    # controls. Column-guarded: an absent column skips the figure rather than failing.
    # ===========================================================================
    try:
        from matplotlib.lines import Line2D as _L2Dtc
        _need_tc = ["active_site_residues_correct", "catalytic_constellation_score",
                    "degrader_tier", "scissile_cf_bde", "sn2_backside_occlusion",
                    "sn2_dead_end", "Ligand_Name"]
        if not all(_c in df.columns for _c in _need_tc):
            reporter.log("  ! Two-criteria figure skipped: required columns absent")
        else:
            _BF = CFG.TIER_CONSTELLATION_MIN["Tier_1A"]
            _BDE_MAX = CFG.SCISSILE_CF_BDE_MAX
            _OCC_MAX = CFG.SN2_BACKSIDE_OCCL_MAX
            _tord = [t for t in TIER_ORDER_LOGIC if t in df["degrader_tier"].unique()]
            _figtc = plt.figure(figsize=(21, 6.4))
            _gstc = _figtc.add_gridspec(1, 3, width_ratios=[1.2, 1.05, 1.05], wspace=0.145)

            # ── panel a : Criterion A (integrity) gates Criterion B (constellation) ──
            _axa = _figtc.add_subplot(_gstc[0, 0])
            _A = pd.to_numeric(df["active_site_residues_correct"], errors="coerce") / 8.0
            _B = pd.to_numeric(df["catalytic_constellation_score"], errors="coerce")
            _sa = pd.DataFrame({"A": _A, "B": _B}).dropna()
            _abins = sorted(_sa["A"].round(3).unique())
            _adata = [_sa.loc[_sa["A"].round(3).eq(b), "B"].values for b in _abins]
            _ans = [len(v) for v in _adata]
            _parts = _axa.violinplot(_adata, positions=range(len(_abins)), widths=0.85, showextrema=False)
            for _pc, _col in zip(_parts["bodies"], plt.cm.RdYlGn(np.linspace(0.1, 0.9, len(_abins)))):
                _pc.set(facecolor=_col, alpha=0.7, edgecolor="#555", linewidth=0.5)
            _amed = [float(np.median(v)) if len(v) else np.nan for v in _adata]
            _axa.plot(range(len(_abins)), _amed, "-D", color="#C0392B", lw=2, ms=6, zorder=6, label="Median B per bin")
            _axa.axhline(_BF, ls="--", color="#C0392B", lw=1.6, label=f"Criterion-B floor for Tier_1A ({_BF:g})")
            for _i, _b in enumerate(_abins):
                _pct = 100.0 * (_sa.loc[_sa["A"].round(3).eq(_b), "B"] >= _BF).mean()
                _axa.text(_i, 1.03, f"{_pct:.0f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color="#2C6FAC")
            _axa.text(0.015, 0.965, f"Top % = share with B ≥ {_BF:g}", transform=_axa.transAxes,
                      fontsize=8, color="#2C6FAC", va="top", ha="left", fontweight="bold")
            _axa.set_xticks(range(len(_abins)))
            _axa.set_xticklabels([f"{b:.3g}\n(n={n:,})" for b, n in zip(_abins, _ans)], fontsize=8)
            _axa.set_ylim(0, 1.18)
            _axa.set_xlabel("Criterion A — active-site integrity\n(fraction of the 8 catalytic residues correctly placed)", fontsize=10)
            _axa.set_ylabel(f"Criterion B — catalytic constellation score\n(reactive-geometry match vs {CFG.REFERENCE_PDB_ID} crystal, 0–1)", fontsize=10)
            _axa.legend(loc="upper right", fontsize=7, ncol=2, framealpha=0.95, columnspacing=0.9, handletextpad=0.4)
            _axa.grid(True, color="#CCCCCC", linewidth=0.6, alpha=0.75, zorder=0); _axa.set_axisbelow(True)
            _axa.text(-0.13, 1.02, "a", transform=_axa.transAxes, fontsize=17, fontweight="bold")

            # ── panel b : Criterion-B ECDF by tier + arrowed per-tier mean values ──
            _axb = _figtc.add_subplot(_gstc[0, 1])
            _means = []
            for _t in _tord:
                _v = pd.to_numeric(df.loc[df["degrader_tier"].eq(_t), "catalytic_constellation_score"],
                                   errors="coerce").dropna().sort_values()
                if len(_v) < 5:
                    continue
                _axb.plot(_v.values, np.linspace(0, 1, len(_v)), color=TIER_PALETTE.get(_t, "#999"), lw=2, label=_t, zorder=3)
                _m = float(_v.mean()); _fr = float((_v <= _m).mean())
                _axb.scatter([_m], [_fr], s=90, color=TIER_PALETTE.get(_t, "#999"), edgecolor="black", lw=0.8, marker="D", zorder=6)
                _means.append((_t, _m, _fr))
            _placetc = {"Tier_5_Decoy": (0.20, 0.90), "Tier_4": (0.20, 0.63), "Tier_3": (0.50, 0.11),
                        "Tier_2B": (0.66, 0.04), "Tier_1A": (0.90, 0.28), "Tier_2A": (0.90, 0.42), "Tier_1B": (0.90, 0.56)}
            for _t, _m, _fr in _means:
                _lx, _ly = _placetc.get(_t, (min(_m + 0.10, 0.94), min(_fr + 0.06, 0.96)))
                _axb.annotate(f"{_t.replace('Tier_', 'T')} = {_m:.2f}", xy=(_m, _fr), xycoords="data",
                              xytext=(_lx, _ly), textcoords=_axb.transAxes, fontsize=8, fontweight="bold",
                              color=TIER_PALETTE.get(_t, "#999"), va="center", ha="center", zorder=8,
                              arrowprops=dict(arrowstyle="->", color=TIER_PALETTE.get(_t, "#999"), lw=1.0, alpha=0.85))
            # MD-ready overlay: super-best MD-selected complexes on their tier's B-ECDF curve.
            _mdb = _md_ready_df(df)
            for _, _mr in _mdb.iterrows():
                _bx = pd.to_numeric(pd.Series([_mr.get("catalytic_constellation_score")]), errors="coerce").iloc[0]
                _bt = _mr.get("degrader_tier")
                if pd.isna(_bx):
                    continue
                _bv = pd.to_numeric(df.loc[df["degrader_tier"].eq(_bt), "catalytic_constellation_score"],
                                    errors="coerce").dropna().sort_values().values
                if len(_bv) == 0:
                    continue
                _axb.scatter([_bx], [float(np.interp(_bx, _bv, np.linspace(0, 1, len(_bv))))], **_MD_STAR_KW)
            if len(_mdb):
                _axb.scatter([], [], marker="*", s=140, facecolor="#FFD400", edgecolor="black",
                             linewidths=0.9, label=f"MD-ready (n={len(_mdb)})")
            _axb.axvline(_BF, ls="--", color="#C0392B", lw=1.4, alpha=0.7)
            _axb.set_xlabel("Criterion B — catalytic constellation score", fontsize=10)
            _axb.set_ylabel("Cumulative fraction", fontsize=10)
            _axb.legend(loc="upper left", ncol=max(len(_means) + (1 if len(_mdb) else 0), 1),
                        fontsize=5.6, handlelength=0.9,
                        handletextpad=0.25, columnspacing=0.5, framealpha=0.92, borderpad=0.3)
            _axb.grid(True, color="#CCCCCC", linewidth=0.6, alpha=0.75, zorder=0); _axb.set_axisbelow(True)
            _axb.text(-0.13, 1.02, "b", transform=_axb.transAxes, fontsize=17, fontweight="bold")

            # ── panel c : SN2 dead-end chemistry gate (penalised ligands + FA/DFA) ──
            _axc = _figtc.add_subplot(_gstc[0, 2])
            _DEAD, _FEAS = "#D62728", "#1F77B4"
            _dd = df.dropna(subset=["scissile_cf_bde", "sn2_backside_occlusion"])
            _gg = _dd.groupby("Ligand_Name").agg(
                bde=("scissile_cf_bde", "median"), occ=("sn2_backside_occlusion", "median"),
                dead=("sn2_dead_end", lambda s: pd.to_numeric(s, errors="coerce").fillna(0).max())).reset_index()
            _gg["pen"] = (_gg["dead"] > 0) | (_gg["bde"] > _BDE_MAX) | (_gg["occ"] > _OCC_MAX)
            _keep = _gg[_gg["pen"] | _gg["Ligand_Name"].isin(["26_Fluoroacetate", "27_Difluoroacetate"])].copy()
            if len(_keep) >= 2:
                _keep["lig"] = _keep["Ligand_Name"].str.replace(r"^\d+_", "", regex=True)
                _keep["isdead"] = _keep["dead"] > 0
                _axc.axhspan(_BDE_MAX, _keep.bde.max() + 3, color=_DEAD, alpha=0.05)
                _axc.axvspan(_OCC_MAX, _keep.occ.max() + 0.4, color=_DEAD, alpha=0.05)
                _keep["_gx"] = _keep["occ"].round(1); _keep["_gy"] = _keep["bde"].round(0)
                for _gk, _grp in _keep.groupby(["_gx", "_gy"]):
                    _yo = 10.0
                    for _, _r in _grp.sort_values("lig").iterrows():
                        _cc = _DEAD if _r["isdead"] else _FEAS
                        _axc.scatter(_r["occ"], _r["bde"], s=85, color=_cc, edgecolor="w", lw=0.5, zorder=3)
                        _axc.annotate(_r["lig"], (_r["occ"], _r["bde"]), xytext=(0, _yo), textcoords="offset points",
                                      rotation=90, ha="center", va="bottom", fontsize=7, color=_cc, fontweight="bold")
                        _yo += len(_r["lig"]) * 4.4 + 7
                _axc.set_ylim(107.0, _keep.bde.max() + 12)
            _axc.axhline(_BDE_MAX, ls="--", color=_DEAD, lw=1.3); _axc.axvline(_OCC_MAX, ls="--", color=_DEAD, lw=1.3)
            _axc.set_xlabel("SN2 backside steric occlusion (Σ vdW, Å)", fontsize=10)
            _axc.set_ylabel("Scissile C–F bond-dissociation energy (kcal/mol)", fontsize=10)
            _axc.legend(handles=[_L2Dtc([], [], marker="o", ls="", color=_DEAD, label="SN2 dead-end (flagged)"),
                                 _L2Dtc([], [], marker="o", ls="", color=_FEAS, label="feasible α-attack (control)"),
                                 _L2Dtc([], [], ls="--", color=_DEAD, label=f"BDE×occ gate ({_BDE_MAX:g} · {_OCC_MAX:g})")],
                        loc="lower center", ncol=3, fontsize=6.3, framealpha=0.95, columnspacing=0.7,
                        handlelength=1.1, handletextpad=0.3, borderpad=0.3)
            _axc.grid(True, color="#CCCCCC", linewidth=0.6, alpha=0.75, zorder=0); _axc.set_axisbelow(True)
            _axc.text(-0.13, 1.02, "c", transform=_axc.transAxes, fontsize=17, fontweight="bold")

            plt.savefig(out_dir / "Figure_30_TwoCriteria_Tier_Logic.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close(_figtc)
    except Exception as e:
        reporter.log(f"  ! Two-criteria figure skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    reporter.section("  Folder 05_Ligand_Interactions_and_Chemical_Space — interaction profile + chemical space")
    # --- Figure 14: Molecular Interaction Profile — Stacked bars (bond types) + F-engagement line ---
    """
    Single chart, dual y-axis: stacked bars per tier show the total interaction count and
    its bond-type composition (left y-axis); overlaid connected dot-line shows the fluorine
    engagement ratio per tier on the right y-axis, revealing whether richer interaction
    profiles correlate with higher fluorine utilisation.
    """
    _has_int   = bool(present_int_cols) and "degrader_tier" in df.columns
    _has_f16   = ("interacting_fluorine_count" in df.columns and
                  "total_fluorine_count" in df.columns and "degrader_tier" in df.columns)
    if _has_int or _has_f16:
        try:
            # Pre-compute tier count for dynamic figure height
            n_tiers_f13 = len([t for t in existing_tiers if t in df["degrader_tier"].values]) if "degrader_tier" in df.columns else 6
            fig, (ax, ax_hm) = plt.subplots(1, 2, figsize=(20, max(5.5, n_tiers_f13 * 0.95 + 2.5)),
                                             gridspec_kw={"width_ratios": [3, 2], "wspace": 0.35})
            ax_r15 = ax.twinx()   # right y-axis for fluorine engagement line

            # Left axis: stacked bar chart — one bar per tier, stacked by bond type
            if _has_int:
                tier_int = df.groupby("degrader_tier")[[present_int_cols[k] for k in present_int_cols]].mean()
                tier_int = tier_int.reindex([t for t in existing_tiers if t in tier_int.index])
                tier_int.columns = [int_label_map.get(k, k) for k in present_int_cols]
                tier_int = tier_int.fillna(0)
                # Normalise to 100 % so each bar shows interaction-type composition
                _tier_int_abs     = tier_int.copy()
                _row_totals_abs   = _tier_int_abs.sum(axis=1)
                tier_int          = _tier_int_abs.div(_row_totals_abs, axis=0).fillna(0).mul(100)
                # Stacked bars: each bar = one tier; stack = bond type proportion
                bond_palette = {
                    "H-Bond": "#4C72B0", "Salt Bridge": "#DD8452",
                    "Halogen": "#55A868", "F-Polar": "#C44E52",
                    "F-Hydrophobic": "#8172B2", "Hydrophobic": "#937860",
                }
                bar_bottom = np.zeros(len(tier_int))
                bar_x = np.arange(len(tier_int))
                tier_labels = list(tier_int.index)
                for bond_type in list(int_label_map.values()):
                    if bond_type not in tier_int.columns:
                        continue
                    vals_pct = tier_int[bond_type].values           # percent for bar height
                    vals_abs = _tier_int_abs[bond_type].values      # absolute mean for annotation
                    ax.bar(bar_x, vals_pct, bottom=bar_bottom,
                           color=bond_palette.get(bond_type, "#999"),
                           edgecolor="white", linewidth=0.5, width=0.55,
                           label=bond_type, zorder=2)
                    # Annotate segment with absolute mean count; skip tiny segments
                    for xi, vi, va, bi in zip(bar_x, vals_pct, vals_abs, bar_bottom):
                        if vi < 3.0:
                            continue
                        mid_y = bi + vi / 2
                        _fs = 7 if vi >= 8.0 else 5.5
                        ax.text(xi, mid_y, f"{va:.1f}",
                                ha="center", va="center", fontsize=_fs,
                                color="white", fontweight="bold", zorder=3)
                    bar_bottom = bar_bottom + vals_pct
                # Total mean count label on top (absolute, not the 100% bar height)
                for xi, tot_abs in zip(bar_x, _row_totals_abs.values):
                    ax.text(xi, 102.0, f"Σ={tot_abs:.0f}",
                            ha="center", va="bottom", fontsize=8, fontweight="bold", color="black", zorder=4)
                ax.set_xticks(bar_x)
                ax.set_xticklabels(tier_labels, rotation=35, ha="right", fontsize=9)
                ax.set_ylim(0, 115)
                ax.set_yticks([0, 25, 50, 75, 100])
                ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=9)
                # Left axis: steel-blue tick labels matching the blue gridlines
                ax.tick_params(axis="y", labelsize=9, labelcolor="#2C6FAC")
                ax.yaxis.label.set_color("#2C6FAC")
                ax.set_xlabel("Degrader Tier", fontsize=11)
                ax.set_ylabel("Interaction type proportion  (%)  —  annotation = mean count", fontsize=11)

            # Right axis: Fluorine Engagement Ratio connected line per tier
            if _has_f16:
                f16_df = df.copy()
                f16_df["FER"] = np.where(
                    f16_df["total_fluorine_count"] > 0,
                    f16_df["interacting_fluorine_count"] / f16_df["total_fluorine_count"],
                    np.nan
                )
                f16_df = f16_df.dropna(subset=["FER"])
                fer_xs, fer_meds, fer_lo_ci, fer_hi_ci = [], [], [], []
                for i, tier in enumerate(existing_tiers):
                    if tier not in f16_df["degrader_tier"].values:
                        continue
                    vals = f16_df.loc[f16_df["degrader_tier"] == tier, "FER"].dropna()
                    if len(vals) < 2:
                        continue
                    med_v = float(vals.median())
                    lo_ci, hi_ci = _median_ci95(vals)
                    fer_xs.append(i)
                    fer_meds.append(med_v)
                    fer_lo_ci.append(max(0, lo_ci))
                    fer_hi_ci.append(min(1, hi_ci))
                if fer_xs:
                    ax_r15.plot(fer_xs, fer_meds, color="#E69F00", linewidth=2.6,
                                marker="o", markersize=9, markeredgecolor="black",
                                markeredgewidth=0.9, zorder=8,
                                label="Fluorine Engagement Ratio  (median ± 95% CI)")
                    ax_r15.fill_between(fer_xs, fer_lo_ci, fer_hi_ci,
                                        color="#E69F00", alpha=0.22, zorder=7)
                    # Filled zone backgrounds on right axis
                    ax_r15.axhspan(0.75, 1.20, alpha=0.15, color="#009E73", zorder=0)
                    ax_r15.axhspan(0.50, 0.75, alpha=0.13, color="#0072B2", zorder=0)
                    ax_r15.axhspan(0.00, 0.50, alpha=0.12, color="#D55E00", zorder=0)
                    ax_r15.axhline(y=0.5, color="#0072B2", linestyle="--",
                                   alpha=0.70, linewidth=1.3, zorder=6)
                    ax_r15.axhline(y=1.0, color="#009E73", linestyle=":",
                                   alpha=0.75, linewidth=1.3, zorder=6)
                    """
                    FER value labels: all placed in a fixed band above the right axis
                    (y = 1.08–1.20) to clear the stacked bars.  Cycle through 3 y-levels
                    and alternate x-offset ±0.15 so adjacent labels never touch.
                    """
                    for _fi15, (xi, yi) in enumerate(zip(fer_xs, fer_meds)):
                        _lbl_y15 = 1.08 + (_fi15 % 3) * 0.05   # 1.08 / 1.13 / 1.18
                        _x_off15 = 0.12 * (1 if _fi15 % 2 == 0 else -1)
                        ax_r15.annotate(f"{yi:.2f}",
                                        xy=(xi, yi), xycoords="data",
                                        xytext=(xi + _x_off15, _lbl_y15), textcoords="data",
                                        ha="center", va="bottom", fontsize=7.5,
                                        color="#A06000", fontweight="bold", zorder=10,
                                        arrowprops=dict(arrowstyle="->", color="#E69F00",
                                                        lw=0.5, mutation_scale=6,
                                                        shrinkA=0, shrinkB=2),
                                        bbox=dict(boxstyle="round,pad=0.09", fc="white",
                                                  ec=CFG.CONF_BAND_COLOURS["acceptable"], alpha=0.94, linewidth=0.5))
                ax_r15.set_ylim(0, 1.30)
                ax_r15.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
                ax_r15.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=9)
                ax_r15.set_ylabel("Fluorine Engagement Ratio  (interacting F / total F per ligand)",
                                  fontsize=11, rotation=270, labelpad=14, color="#A06000")
                ax_r15.tick_params(axis="y", labelcolor="#A06000", labelsize=9)
                """
                Rotated zone labels — left side of figure, written bottom-to-top
                ylim=(0, 1.30): Low midpoint=0.25→0.192; Mid midpoint=0.625→0.481; High midpoint=0.875→0.673
                """
                ax_r15.text(0.99, 0.192, "Low engagement (<50%)",
                            ha="center", va="center", rotation=90,
                            fontsize=7.5, color="#D55E00", style="italic", fontweight="bold",
                            transform=ax_r15.transAxes,
                            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["below"],
                                      alpha=0.80, linewidth=0.5))
                ax_r15.text(0.99, 0.481, "50% engagement",
                            ha="center", va="center", rotation=90,
                            fontsize=7.5, color="#0072B2", style="italic", fontweight="bold",
                            transform=ax_r15.transAxes,
                            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec="#0072B2",
                                      alpha=0.80, linewidth=0.5))
                ax_r15.text(0.99, 0.673, "Full engagement",
                            ha="center", va="center", rotation=90,
                            fontsize=7.5, color="#009E73", style="italic", fontweight="bold",
                            transform=ax_r15.transAxes,
                            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["high"],
                                      alpha=0.80, linewidth=0.5))

            """
            Dual-colour grid system: left = steel-blue (matches left tick labels),
            right = amber (matches right tick labels). Both BEHIND bars via set_axisbelow.
            """
            ax.set_axisbelow(True)
            ax.yaxis.grid(True, color="#2C6FAC", linewidth=0.55, linestyle="-",
                          alpha=0.20, zorder=0)
            ax.xaxis.grid(False)
            ax_r15.yaxis.grid(True, color="#A06000", linewidth=0.55, linestyle="--",
                              alpha=0.18, zorder=0)
            ax_r15.set_axisbelow(True)

            """
            Split legend: bond-type bars in lower-left (typically sparse area),
            FER line separately via ax_r15 annotation to avoid overlap with bars.
            """
            handles_bars, labels_bars = ax.get_legend_handles_labels()
            handles_line, labels_line = ax_r15.get_legend_handles_labels()
            _all_h15 = handles_bars + handles_line
            _all_l15 = labels_bars + labels_line
            # Single row across the top; if too many items, allow a second row (ncol = ceil/2)
            import math as _math15
            _ncol15 = len(_all_h15) if len(_all_h15) <= 8 else _math15.ceil(len(_all_h15) / 2)
            _leg15 = ax.legend(_all_h15, _all_l15,
                               loc="upper left", ncol=_ncol15,
                               fontsize=7.5, framealpha=0.93, fancybox=True,
                               borderpad=0.4, labelspacing=0.3,
                               handlelength=1.0, handletextpad=0.4,
                               columnspacing=0.8)
            _leg15.set_zorder(20)

            # Heatmap: tier × bond-type proportion (%) — only when interaction data exist
            if _has_int and "tier_int" in dir():
                import seaborn as _sns13
                _hm_data = tier_int.T  # bond types as rows, tiers as columns
                # tier_int is already normalised to 100% per bond type
                _sns13.heatmap(_hm_data, ax=ax_hm,
                               cmap="YlOrRd", vmin=0, vmax=100,
                               annot=True, fmt=".0f", annot_kws={"size": 8},
                               linewidths=0.4, linecolor="white",
                               cbar_kws={"label": "Proportion (%)", "shrink": 0.8})
                ax_hm.set_xlabel("Degrader Tier  —  Interaction Type Proportion (%)", fontsize=9.5)
                ax_hm.set_ylabel("Interaction Type", fontsize=10)
                ax_hm.tick_params(axis="x", labelrotation=45, labelsize=8.5)
                ax_hm.tick_params(axis="y", labelrotation=0,  labelsize=8)
            else:
                ax_hm.set_visible(False)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_14a_Molecular_Interaction_Profile.png",
                        dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Figure 14 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # Figure 14b: PA companion (interaction profile overlay)
    if _tt_has_imgs:
        _fig_14b_tt_interactions(df, _pa, _imgs, out_dir, reporter)

    # --- Figure 15: Fluorine Engagement Ratio by Tier (box + strip + median trend line) ---
    if "interacting_fluorine_count" in df.columns and "total_fluorine_count" in df.columns and "degrader_tier" in df.columns:
        try:
            f16_df = df.copy()
            f16_df["FER"] = np.where(
                f16_df["total_fluorine_count"] > 0,
                f16_df["interacting_fluorine_count"] / f16_df["total_fluorine_count"],
                np.nan)
            f16_df = f16_df.dropna(subset=["FER"])
            valid_t16 = [t for t in existing_tiers if t in f16_df["degrader_tier"].values]
            fig, ax = plt.subplots(figsize=(12, 7))

            ax.axhspan(0.75, 1.01, alpha=0.07, color="#009E73", zorder=0)
            ax.axhspan(0.50, 0.75, alpha=0.06, color="#E69F00", zorder=0)
            ax.axhspan(0.00, 0.50, alpha=0.06, color="#D55E00", zorder=0)
            ax.axhline(y=0.75, color="#009E73", linestyle="--", alpha=0.6, linewidth=1.1)
            ax.axhline(y=0.50, color="#E69F00", linestyle=":", alpha=0.6, linewidth=1.1)
            ax.set_axisbelow(True)
            # KDE violin + IQR box overlay (no boxplot — keeps figure clean)
            _skipped14 = []   # tiers with too few points for a KDE violin (kept visible via a note)
            for _ti14, tier14 in enumerate(valid_t16):
                _vals14 = f16_df.loc[f16_df["degrader_tier"] == tier14, "FER"].dropna().values
                if len(_vals14) < 4:
                    # Too few complexes for a stable KDE — draw the raw points so the
                    # tier is not silently omitted, and record it for the caption.
                    if len(_vals14):
                        ax.scatter([_ti14] * len(_vals14), _vals14, s=26, zorder=6,
                                   color=TIER_PALETTE.get(tier14, "#999"),
                                   edgecolors="black", linewidths=0.5)
                    _skipped14.append((tier14, len(_vals14)))
                    continue
                from scipy.stats import gaussian_kde as _gkde14
                try:
                    _kde14  = _gkde14(_vals14)
                    _y14    = np.linspace(np.min(_vals14), np.max(_vals14), 200)
                    _w14    = _kde14(_y14)
                    _w14   /= _w14.max()
                    _w14   *= 0.35   # half-width of violin
                    ax.fill_betweenx(_y14, _ti14 - _w14, _ti14 + _w14,
                                     alpha=0.55, color=TIER_PALETTE.get(tier14, "#999"))
                except (np.linalg.LinAlgError, ValueError):
                    # Degenerate (collinear) FER values give a singular KDE covariance;
                    # render the raw points as a jittered strip in place of the violin.
                    _jit14 = np.random.default_rng(0).uniform(-0.12, 0.12, size=len(_vals14))
                    ax.scatter(_ti14 + _jit14, _vals14, s=12, alpha=0.45,
                               color=TIER_PALETTE.get(tier14, "#999"),
                               edgecolors="none", zorder=2)
                # IQR box on top
                _q25, _q75 = np.percentile(_vals14, [25, 75])
                _med14 = np.median(_vals14)
                ax.vlines(_ti14, _q25, _q75, color="#333", linewidth=2.5, zorder=5)
                ax.scatter([_ti14], [_med14], color="white", s=30, zorder=6,
                           edgecolors="#333", linewidths=1.2)

            # Per-tier median statistics badge
            _med16_xs, _med16_ys = [], []
            for i, tier in enumerate(valid_t16):
                sub16 = f16_df.loc[f16_df["degrader_tier"] == tier, "FER"].dropna()
                if len(sub16) == 0:
                    continue
                med = float(sub16.median())
                n16 = len(sub16)
                _med16_xs.append(i)
                _med16_ys.append(med)
            # FER connecting line across tier medians
            if len(_med16_xs) >= 2:
                ax.plot(_med16_xs, _med16_ys, color="#E69F00", linewidth=2.2,
                        linestyle="--", alpha=0.85, zorder=4,
                        marker="D", markersize=6, markeredgecolor="#333",
                        markeredgewidth=0.7, label="FER trend (median)")
            for i, tier in enumerate(valid_t16):
                sub16 = f16_df.loc[f16_df["degrader_tier"] == tier, "FER"].dropna()
                if len(sub16) == 0:
                    continue
                med = float(sub16.median())
                n16 = len(sub16)
                ax.text(i, med + 0.022, f"med={med:.2f}\nn={n16:,}",
                        ha="center", va="bottom", fontsize=6.0, fontweight="bold",
                        zorder=7, color="#111",
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#ccc",
                                  linewidth=0.6, alpha=0.85))

            ax.set_ylim(0, 1.08)
            ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
            ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
            # Color y-tick labels to match zone backgrounds
            _ytick_colors16 = {
                "0%":   "#B84000",
                "25%":  "#B84000",
                "50%":  "#8A6000",
                "75%":  "#007A50",
                "100%": "#007A50",
            }
            for tick16 in ax.get_yticklabels():
                tick16.set_color(_ytick_colors16.get(tick16.get_text(), "#111"))
                tick16.set_fontweight("bold")
            """
            Start-y (axes fraction) = bottom border of each label's OWN zone, so
            every label sits on its category's lower boundary line (va='bottom',
            text rises into its zone). Low's zone bottom coincides with the axes
            floor, which is expected — High/Moderate now match that anchoring.
            """
            for (_zy14_start, _ztxt14, _zcol14, _zec14) in [
                (0.73, f"High engagement (≥{CFG.ENGAGEMENT_BAND_HIGH*100:.0f}%)", "#007A50", "#009E73"),
                (0.50, f"Moderate ({CFG.ENGAGEMENT_BAND_MODERATE*100:.0f}–{CFG.ENGAGEMENT_BAND_HIGH*100:.0f}%)", "#8A6000", "#E69F00"),
                (0.10, f"Low engagement (<{CFG.ENGAGEMENT_BAND_MODERATE*100:.0f}%)", "#B84000", "#D55E00"),
            ]:
                ax.text(0.985, _zy14_start, _ztxt14,
                        ha="center", va="bottom", clip_on=True,
                        fontsize=5.5, color=_zcol14, style="italic",
                        fontweight="bold", rotation=90, transform=ax.transAxes,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=_zec14,
                                  alpha=0.80, linewidth=0.5))
            ax.set_xlabel("Degrader Tier", fontsize=11)
            ax.set_ylabel("Fluorine Engagement Ratio  (interacting F / total F)", fontsize=11)
            ax.set_xticks(range(len(valid_t16)))
            ax.set_xticklabels(valid_t16, rotation=35, ha="right", fontsize=9)
            for tick16x, tier16x in zip(ax.get_xticklabels(), valid_t16):
                tick16x.set_color(TIER_PALETTE.get(tier16x, "black"))
                tick16x.set_fontweight("bold")

            # Legend — violin + IQR box description
            from matplotlib.patches import Patch as _P16
            from matplotlib.lines import Line2D as _L16
            _leg16_h = [
                _P16(facecolor="#888", alpha=0.55, label="KDE violin (per tier)"),
                _L16([0], [0], color="#333", linewidth=2.5, label="IQR (25–75%)"),
                _L16([0], [0], color="white", marker="o", markersize=6,
                     markeredgecolor="#333", markeredgewidth=1.2,
                     linestyle="None", label="Median"),
            ]
            _mdh16 = _md_ready_stars_cat(ax, f16_df, valid_t16, "FER")
            if _mdh16 is not None:
                _leg16_h.append(_mdh16)
            _leg16 = ax.legend(handles=_leg16_h,
                               loc="lower left", bbox_to_anchor=(0.01, 0.01),
                               fontsize=6, framealpha=0.92, fancybox=True, ncol=4,
                               labelspacing=0.2, handlelength=0.8, handletextpad=0.3)
            _leg16.set_zorder(20)

            if _skipped14:
                _sk_txt = ", ".join(f"{t} (n={n})" for t, n in _skipped14)
                ax.text(0.99, 0.015, f"n<4, shown as points (no KDE): {_sk_txt}",
                        transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=6.5, style="italic", color="#555", zorder=20)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_15_Fluorine_Engagement_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Figure 15 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 16: Catalytic quality vs active-site contact density per tier ---
    """
    Single chart, dual y-axis. The left axis shows the per-tier distribution of catalytic
    quality as violins — soft_catalytic_score (the continuous SN2-geometry composite),
    which genuinely varies across tiers. The right axis shows the per-tier mean
    active-site contact density (Interaction_Density_Norm) as a connected dot-line, pairing
    catalytic geometry with active-site engagement. Binding Probability is not plotted:
    it saturates near 0.99 across all tiers and carries no per-tier structure. Product
    inhibition is assessed downstream by the Step-06 MM-GBSA stage.
    """
    _qual_col  = next((c for c in ["soft_catalytic_score", "SN2_Attack_Angle", "mechanistic_score"] if c in df.columns), None)
    _qual_lbl  = {"soft_catalytic_score": "Soft catalytic score  (SN2-geometry composite)",
                  "SN2_Attack_Angle": "SN2 attack angle (°)",
                  "mechanistic_score": "Mechanistic score"}.get(_qual_col, str(_qual_col))
    _qual_fmt  = "{:.1f}" if _qual_col == "SN2_Attack_Angle" else "{:.3f}"
    _has_bind  = _qual_col is not None and "degrader_tier" in df.columns
    _has_dens  = "Interaction_Density_Norm" in df.columns and "degrader_tier" in df.columns
    if _has_bind or _has_dens:
        try:
            from scipy import stats as _sc_stats12
            fig, ax = plt.subplots(figsize=(13, 7))
            ax_r = ax.twinx()   # right y-axis for contact density line

            # Left axis: catalytic-quality violin per tier
            if _has_bind:
                f12_data = df.dropna(subset=[_qual_col])
                _qv12 = pd.to_numeric(f12_data[_qual_col], errors="coerce").dropna()
                _qlo12, _qhi12 = float(_qv12.quantile(0.01)), float(_qv12.quantile(0.99))
                _qpad12 = (_qhi12 - _qlo12) * 0.08 or 0.02
                sns.violinplot(data=f12_data, x="degrader_tier", y=_qual_col, order=existing_tiers,
                               palette=TIER_PALETTE, inner="box", linewidth=1.2, cut=0,
                               alpha=0.75, ax=ax)
                for i, tier in enumerate(existing_tiers):
                    med = f12_data.loc[f12_data["degrader_tier"] == tier, _qual_col].median()
                    if not np.isnan(med):
                        # Red median line across the violin width
                        ax.hlines(med, i - 0.22, i + 0.22, colors="#CC0000",
                                  linewidth=2.0, zorder=7, linestyle="-")
                        ax.text(i + 0.25, med, "  " + _qual_fmt.format(med), va="center", ha="left",
                                fontsize=7, color="#CC0000", fontweight="bold", zorder=8)
                # Top extends to 1.0 (0–1 composite) so the upper violin tails are
                # not clipped; floor keeps the zoomed lower bound.
                _qtop12 = min(1.0, max(_qhi12 + _qpad12, float(_qv12.max()) + _qpad12))
                if _qtop12 >= 0.9:
                    _qtop12 = 1.0
                ax.set_ylim(max(0.0, _qlo12 - _qpad12), _qtop12)
                ax.set_ylabel(_qual_lbl, fontsize=11)
            else:
                ax.set_visible(False)

            # Right axis: per-tier mean active-site contact density connected dot-line
            if _has_dens:
                f20_df = df.dropna(subset=["Interaction_Density_Norm"]).copy()
                valid_tiers_f20 = [t for t in existing_tiers if t in f20_df["degrader_tier"].values]
                pip_xs, pip_meds, pip_lo_ci, pip_hi_ci = [], [], [], []
                for tier in valid_tiers_f20:
                    x_pos = existing_tiers.index(tier) if tier in existing_tiers else 0
                    vals  = f20_df.loc[f20_df["degrader_tier"] == tier, "Interaction_Density_Norm"].dropna()
                    if len(vals) < 2:
                        continue
                    mean_v = float(vals.mean())
                    ci95  = float(_sc_stats12.sem(vals) * _sc_stats12.t.ppf(0.975, len(vals) - 1))
                    pip_xs.append(x_pos)
                    pip_meds.append(mean_v)
                    pip_lo_ci.append(mean_v - ci95)
                    pip_hi_ci.append(mean_v + ci95)
                if pip_xs:
                    ax_r.plot(pip_xs, pip_meds, color="#D55E00", linewidth=2.5,
                              marker="D", markersize=9, markeredgecolor="black",
                              markeredgewidth=0.9, zorder=7, label="Active-site contact density\n(mean ± 95 % CI)")
                    ax_r.fill_between(pip_xs, pip_lo_ci, pip_hi_ci, color="#D55E00",
                                      alpha=0.20, zorder=6)
                    # Mean annotation
                    for xi, yi in zip(pip_xs, pip_meds):
                        ax_r.text(xi, yi + (max(pip_meds) - min(pip_meds)) * 0.04,
                                  f"{yi:.2f}", ha="center", va="bottom", fontsize=8,
                                  color="#D55E00", fontweight="bold", zorder=8)
                ax_r.set_ylabel("Active-site contact density  (interactions per complex)",
                                fontsize=11, rotation=270, labelpad=14, color="#D55E00")
                ax_r.tick_params(axis="y", labelcolor="#D55E00", labelsize=9)

            ax.set_xticks(range(len(existing_tiers)))
            ax.set_xticklabels(existing_tiers, rotation=40, ha="right", fontsize=9)
            for _tick12x, _tier12x in zip(ax.get_xticklabels(), existing_tiers):
                _tick12x.set_color(TIER_PALETTE.get(_tier12x, "black"))
                _tick12x.set_fontweight("bold")
            # Left axis (binding prob.): steel-blue labels matching blue gridlines
            ax.tick_params(axis="y", labelsize=9, labelcolor="#2C6FAC")
            ax.yaxis.label.set_color("#2C6FAC")
            # Right axis (contact density): amber labels matching amber gridlines
            ax_r.tick_params(axis="y", labelsize=9, labelcolor="#A06000")
            ax_r.yaxis.label.set_color("#A06000")
            ax.set_xlabel("Degrader Tier", fontsize=11)

            # Dual-colour grids: left-axis = steel-blue; right-axis = amber
            ax.set_axisbelow(True)
            ax.yaxis.grid(True, color="#2C6FAC", linewidth=0.55, linestyle="-", alpha=0.20, zorder=0)
            ax.xaxis.grid(False)
            ax_r.yaxis.grid(True, color="#A06000", linewidth=0.55, linestyle="--", alpha=0.18, zorder=0)
            ax_r.set_axisbelow(True)

            # Legend — bottom-left (violins tend to be taller on right)
            if _has_dens and pip_xs:
                from matplotlib.lines import Line2D as _Line12
                _pip_handle = _Line12([0], [0], color="#D55E00", linewidth=2.5,
                                      marker="D", markersize=8,
                                      label="Active-site contact density  (mean ± 95% CI)")
                _leg12 = ax_r.legend(handles=[_pip_handle], loc="lower left",
                                     fontsize=8.5, framealpha=0.92, fancybox=True)
                _leg12.set_zorder(20)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_16_Catalytic_Quality_vs_Inhibition.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Figure 16 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 17: Active-site contact density by tier ---
    if "Interaction_Density_Norm" in df.columns and "degrader_tier" in df.columns:
        try:
            f20_df = df.dropna(subset=["Interaction_Density_Norm"]).copy()
            valid_t20 = [t for t in existing_tiers if t in f20_df["degrader_tier"].values]
            if not f20_df.empty and valid_t20:
                fig, ax = plt.subplots(figsize=(12, 6.5))
                _dens_all = pd.to_numeric(f20_df["Interaction_Density_Norm"],
                                          errors="coerce").dropna()
                # Clamp the visible range to p99 so outlier tails (cavity mis-detection)
                # do not compress the populated 0–p99 band into the lower axis.
                _p99_20 = float(_dens_all.quantile(0.99))
                _d_max  = float(min(max(1.0, _p99_20 * 1.08),
                                    float(_dens_all.max()) * 1.10 or 1.0))
                _n_hidden20 = int((_dens_all > _d_max).sum())
                # Data-driven low / acceptable / high zones from global tertiles.
                _z_lo20, _z_hi20 = (float(_dens_all.quantile(0.33)),
                                    float(_dens_all.quantile(0.66)))
                ax.axhspan(_z_hi20, _d_max, alpha=0.07, color="#009E73", zorder=0)
                ax.axhspan(_z_lo20, _z_hi20, alpha=0.06, color="#E69F00", zorder=0)
                ax.axhspan(0.0, _z_lo20, alpha=0.06, color="#D55E00", zorder=0)
                ax.axhline(_z_hi20, color="#009E73", linestyle="--", alpha=0.55, linewidth=1.0)
                ax.axhline(_z_lo20, color="#D55E00", linestyle=":", alpha=0.55, linewidth=1.0)
                sns.boxplot(data=f20_df, x="degrader_tier", y="Interaction_Density_Norm",
                            order=valid_t20, palette=TIER_PALETTE, linewidth=1.2,
                            showfliers=False, ax=ax)
                sns.stripplot(data=f20_df, x="degrader_tier", y="Interaction_Density_Norm",
                              order=valid_t20, color="black", alpha=0.10, size=2.0, jitter=True, ax=ax)
                # Per-tier mean trend line (diamond markers) over the median boxes.
                _mean20_xs, _mean20_ys = [], []
                for i, tier in enumerate(valid_t20):
                    sub20 = f20_df.loc[f20_df["degrader_tier"] == tier, "Interaction_Density_Norm"].dropna()
                    if len(sub20) == 0:
                        continue
                    _mean20_xs.append(i)
                    _mean20_ys.append(float(sub20.mean()))
                if len(_mean20_xs) >= 2:
                    ax.plot(_mean20_xs, _mean20_ys, color="#7A0177", linewidth=2.0,
                            linestyle="--", alpha=0.9, zorder=6, marker="D", markersize=6,
                            markeredgecolor="white", markeredgewidth=0.7,
                            label="Mean trend")
                for i, tier in enumerate(valid_t20):
                    sub20 = f20_df.loc[f20_df["degrader_tier"] == tier, "Interaction_Density_Norm"].dropna()
                    if len(sub20) == 0:
                        continue
                    med = float(sub20.median())
                    _bg_col20 = TIER_PALETTE.get(tier, "#666")
                    # Adaptive text: white on dark backgrounds, near-black on light ones
                    _txt_col20 = _text_color(_bg_col20)
                    ax.text(i, med, f"{med:.1f}", ha="center", va="center",
                            fontsize=8.5, fontweight="bold", zorder=7,
                            color=_txt_col20,
                            bbox=dict(boxstyle="round,pad=0.15", fc=_bg_col20,
                                      ec="white", linewidth=0.6, alpha=0.88))
                ax.set_xlabel("Degrader Tier", fontsize=11)
                ax.set_ylabel("Active-site contact density  (interactions per complex)", fontsize=11)
                ax.set_xticks(range(len(valid_t20)))
                ax.set_xticklabels(valid_t20, rotation=35, ha="right", fontsize=9)
                for _tick20x, _tier20x in zip(ax.get_xticklabels(), valid_t20):
                    _tick20x.set_color(TIER_PALETTE.get(_tier20x, "black"))
                    _tick20x.set_fontweight("bold")
                ax.set_axisbelow(True)
                ax.yaxis.grid(True, color="#2C6FAC", linewidth=0.55, linestyle="-",
                              alpha=0.22, zorder=0)
                ax.xaxis.grid(False)
                ax.tick_params(axis="y", labelcolor="#2C6FAC", labelsize=9)
                ax.yaxis.label.set_color("#2C6FAC")
                ax.set_ylim(bottom=0, top=_d_max)
                # Zone labels parked at the right margin.
                for _zy20, _zt20, _zc20 in [((_z_hi20 + _d_max) / 2, "High", "#1B7A4B"),
                                            ((_z_lo20 + _z_hi20) / 2, "Moderate", "#8A6000"),
                                            (_z_lo20 / 2, "Low", "#B84000")]:
                    ax.text(0.995, _zy20, _zt20, transform=ax.get_yaxis_transform(),
                            ha="right", va="center", fontsize=7.5, color=_zc20,
                            fontweight="bold", alpha=0.85, zorder=8)
                if _n_hidden20 > 0:
                    ax.text(0.01, 0.94,
                            f"{_n_hidden20:,} complexes > {_d_max:.1f} hidden (axis clamped at p99)",
                            transform=ax.transAxes, ha="left", va="top",
                            fontsize=7, color="#777", style="italic")
                ax.legend(loc="upper right", fontsize=7.5, framealpha=0.9)

                # Kruskal-Wallis significance note
                try:
                    from scipy.stats import kruskal as _kw20
                    _kw_groups = [f20_df.loc[f20_df["degrader_tier"] == t,
                                             "Interaction_Density_Norm"].dropna().values
                                  for t in valid_t20 if len(f20_df[f20_df["degrader_tier"]==t]) >= 3]
                    if len(_kw_groups) >= 2:
                        _, _kw_p = _kw20(*_kw_groups)
                        _pstr = ("p<0.001" if _kw_p < 0.001 else
                                 f"p={_kw_p:.3f}")
                        ax.text(0.01, 0.99, f"Kruskal–Wallis {_pstr}  (across tiers)",
                                transform=ax.transAxes, ha="left", va="top",
                                fontsize=8, color="#333", style="italic")
                except Exception:
                    pass

                plt.tight_layout()
                plt.savefig(out_dir / "Figure_17_ActiveSite_Contact_Density_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
                plt.close()
            else:
                reporter.log("  ! Figure 17 skipped: no valid tier groups for Interaction_Density_Norm.")
        except Exception as e:
            reporter.log(f"  ! Figure 17 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig
    else:
        reporter.log("  ! Figure 17 skipped: column 'Interaction_Density_Norm' or 'degrader_tier' absent.")

    # --- Figure 17b: Binding Energetics — Binding Probability (violin) per tier ---
    """
    Per-tier read-out of binding quality. The Binding_Probability distribution is shown
    as a violin (zoomed to where values cluster) with a red median bar per tier. Product
    inhibition is not plotted here; it is assessed downstream by the Step-06 MM-GBSA stage.
    """
    _be_bind = "Binding_Probability" if "Binding_Probability" in df.columns else None
    _be_has_bind = _be_bind is not None and "degrader_tier" in df.columns
    if _be_has_bind:
        fig = None
        try:
            fig, ax = plt.subplots(figsize=(13, 7))
            _bedata = df.dropna(subset=[_be_bind])
            _be_lo = max(0.0, float(_bedata[_be_bind].quantile(0.01)) - 0.003)
            # Headroom = 8% of the visible span so the upper violin tails are not
            # clipped at the frame (binding probability is tightly clustered).
            _be_max = float(_bedata[_be_bind].max())
            _be_hi = min(1.0, _be_max + max(0.0015, (_be_max - _be_lo) * 0.08))
            sns.violinplot(data=_bedata, x="degrader_tier", y=_be_bind, order=existing_tiers,
                           palette=TIER_PALETTE, inner="box", linewidth=1.2, cut=0,
                           alpha=0.75, ax=ax)
            _be_mxs, _be_mys = [], []
            for _i, _t in enumerate(existing_tiers):
                _m = _bedata.loc[_bedata["degrader_tier"] == _t, _be_bind].median()
                if not np.isnan(_m):
                    _be_mxs.append(_i); _be_mys.append(float(_m))
            # Single connected median trend (diamonds joined) — no per-violin bars.
            if len(_be_mxs) >= 2:
                ax.plot(_be_mxs, _be_mys, color="#CC0000", linewidth=2.2, linestyle="--",
                        alpha=0.9, marker="D", markersize=7, markeredgecolor="black",
                        markeredgewidth=0.7, zorder=8, label="Median trend")
                # Boxed value labels just below each node (Fig 03 style).
                _be_span = (_be_hi - _be_lo) if _be_hi > _be_lo else 1.0
                for _mx, _my in zip(_be_mxs, _be_mys):
                    ax.text(_mx, _my - _be_span * 0.03, f"{_my:.3f}", ha="center", va="top",
                            fontsize=6.5, fontweight="bold", color="#CC0000", zorder=9,
                            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#CCCCCC",
                                      linewidth=0.5, alpha=0.85))
            if _be_hi > _be_lo:
                ax.set_ylim(_be_lo, _be_hi)
            ax.set_ylabel("Binding probability  (sigmoid of SN2 geometry)", fontsize=11, color="#2C6FAC")
            ax.tick_params(axis="y", labelcolor="#2C6FAC", labelsize=9)
            ax.set_xticks(range(len(existing_tiers)))
            ax.set_xticklabels(existing_tiers, rotation=40, ha="right", fontsize=9)
            for _txk, _tt in zip(ax.get_xticklabels(), existing_tiers):
                _txk.set_color(TIER_PALETTE.get(_tt, "black")); _txk.set_fontweight("bold")
            ax.set_xlabel("Degrader tier", fontsize=11)
            ax.set_axisbelow(True)
            ax.yaxis.grid(True, color="#2C6FAC", linewidth=0.55, alpha=0.20, zorder=0)
            ax.xaxis.grid(False)
            from matplotlib.patches import Patch as _Patch17b
            from matplotlib.lines import Line2D as _Line17b
            _lh17b = [
                _Patch17b(facecolor="#9DB8D2", edgecolor="#2C6FAC",
                          label="Binding probability (per-tier violin)"),
                _Line17b([0], [0], color="#CC0000", linewidth=2.2,
                         label="Binding probability median"),
            ]
            # Single row, anchored top-left above the axes (keeps it off the violins).
            ax.legend(handles=_lh17b, loc="lower left", bbox_to_anchor=(0.0, 1.01),
                      ncol=len(_lh17b), fontsize=8, framealpha=0.92, fancybox=True,
                      columnspacing=1.2, handletextpad=0.5, borderaxespad=0.0).set_zorder(20)
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_17b_Binding_Energetics.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        except Exception as e:
            reporter.log(f"  ! Figure 17b skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig
        finally:
            if fig is not None:
                plt.close(fig)

    # --- Figure 18a: Chemical Space UMAP Manifold ---
    if "UMAP_X" in df.columns:
        # Degradability landscape: hexbin over the UMAP chemical space coloured by
        # the MEAN competence per bin, so chemical regions enriched for good degraders
        # light up. The MD-ready / elite hits are overlaid as labelled stars. This
        # replaces the former 58k-point tier scatter, which was an unreadable hairball.
        fig, ax = plt.subplots(figsize=(11, 7.5))
        _u = df.dropna(subset=["UMAP_X", "UMAP_Y"]).copy()
        _ux, _uy = _u["UMAP_X"].values, _u["UMAP_Y"].values
        _cmetric = next((c for c in ("competence_score", "mechanistic_score",
                                     "Boltz_Model_Confidence") if c in _u.columns), None)
        if _cmetric:
            _cvals = pd.to_numeric(_u[_cmetric], errors="coerce").fillna(0.0).values
            hb = ax.hexbin(_ux, _uy, C=_cvals, reduce_C_function=np.mean,
                           gridsize=45, cmap="viridis", mincnt=1, linewidths=0.15)
            cb = fig.colorbar(hb, ax=ax, fraction=0.035, pad=0.01)
            cb.set_label(f"Mean {_cmetric.replace('_', ' ')} per bin  (brighter = more degradable)",
                         fontsize=9)
        else:
            ax.hexbin(_ux, _uy, gridsize=45, cmap="Greys", mincnt=1)

        _md_col = getattr(CFG, "MD_SELECTED_COL", "MD_Selected")
        if _md_col in _u.columns and _u[_md_col].astype(str).str.lower().isin(["true", "1", "1.0"]).any():
            _elite = _u[_u[_md_col].astype(str).str.lower().isin(["true", "1", "1.0"])]
            _elite_lbl = f"MD-ready hits (n={len(_elite)})"
        else:
            _elite = _u[_u["degrader_tier"] == CFG.TIER_TOP]
            _elite_lbl = f"{CFG.TIER_TOP} hits (n={len(_elite)})"
        if not _elite.empty:
            ax.scatter(_elite["UMAP_X"], _elite["UMAP_Y"], marker="*", s=185,
                       color="#FFD400", edgecolors="black", linewidths=0.9,
                       zorder=8, label=_elite_lbl)
            if "Ligand_Name" in _elite.columns:
                _lig = _elite["Ligand_Name"].astype(str).str.replace(r"^\d+_", "", regex=True)
                for _xx, _yy, _lg in zip(_elite["UMAP_X"], _elite["UMAP_Y"], _lig):
                    ax.annotate(_lg, (_xx, _yy), textcoords="offset points", xytext=(5, 4),
                                fontsize=6.5, fontweight="bold", color="#1A1A1A", zorder=9)

        ax.set_xlabel("UMAP Dimension 1  (distances reflect chemical similarity)", fontsize=10)
        ax.set_ylabel("UMAP Dimension 2", fontsize=10)
        ax.grid(False)
        ax.legend(loc="lower left", fontsize=8, framealpha=0.92)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_18a_Chemical_Space_Map.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # Figure 18b: PA companion (chemical space landscape)
    if _tt_has_imgs:
        _fig_18b_tt_landscape(df, _pa, _imgs, out_dir, reporter)

    reporter.section("  Folder 06_PFAS_Scope_and_Synthesis — synthesis + publication assembly")
    # --- Figure 19: Candidate Radar: 18a (top-5 hits) + 18b (one per tier) ---
    _radar_labels = {
        "Boltz_Model_Confidence":   "AI Conf.",
        "iptm":                     "ipTM",
        "Interaction_Density_Norm": "Int.Den",
        "mean_plddt":               "pLDDT",
        "Binding_Probability":      "Bind.Prob",
        "SN2_Attack_Angle":         "SN2(°)",
        "mechanistic_score": "Mech",
        "Active_Site_RMSD":         "RMSD(Å)",
        "identity_pct":             "Seq.ID%",
        "SN2_Trajectory_Deviation_A": "Traj.Dev",
        "soft_catalytic_score":     "Soft.Cat",
    }

    def _draw_radar(rows_df, baseline_df, metrics, fig_path, subtitle,
                    label_mode="tier"):
        '''
        label_mode='rank'  → "Fluoroacetate (Rank 1)"
        label_mode=f'tier'  → "Fluoroacetate (Tier_1A)"
        '''
        if not metrics:
            return
        _all = pd.concat([rows_df[metrics], baseline_df[metrics]]) if not baseline_df.empty else rows_df[metrics]
        _scaler_vals = _all.copy()
        for c in metrics:
            mn, mx = float(_all[c].min()), float(_all[c].max())
            _scaler_vals[c] = (_all[c] - mn) / (mx - mn + 1e-9)
        N = len(metrics)
        angles = [n / float(N) * 2 * np.pi for n in range(N)] + [0]
        spoke_labels = [_radar_labels[m] for m in metrics]
        _, ax_r = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
        colours_r = ["#057759", "#0BF1E2", "#E69F00", "#CC79A7", "#0072B2",
                     "#56B4E9", "#F0E442", "#009E73", "#D55E00", "#CC79A7"]
        lname_col = next((c for c in ["Ligand_Name", "ligand"] if c in rows_df.columns), None)

        _ring_levels = [0.2, 0.4, 0.6, 0.8, 1.0]
        _ring_angles = np.linspace(0, 2 * np.pi, 300)
        for _rl in _ring_levels:
            ax_r.plot(_ring_angles, [_rl] * 300,
                      color="#CC2222", linewidth=0.55, alpha=0.45, zorder=1,
                      linestyle="-", solid_capstyle="round")

        import re as _re18
        for i in range(len(rows_df)):
            vals  = _scaler_vals.iloc[i].values.flatten().tolist() + [_scaler_vals.iloc[i].values[0]]
            _raw  = str(rows_df.iloc[i][lname_col]) if lname_col else f"Hit {i+1}"
            # Strip leading "26_" style number prefix
            _clean = _re18.sub(r"^\d+_", "", _raw)
            tier  = str(rows_df.iloc[i].get("degrader_tier", ""))
            if label_mode == "rank":
                _rank_val = int(rows_df.iloc[i].get("Scientific_Rank", i + 1))
                lname_lbl = f"{_clean}  (Rank {_rank_val})"
            else:
                lname_lbl = f"{_clean}  ({tier})"
            c = colours_r[i % len(colours_r)]
            ax_r.plot(angles, vals, linewidth=2.2, linestyle="solid",
                      label=lname_lbl, color=c, zorder=3)
            ax_r.fill(angles, vals, color=c, alpha=0.08, zorder=2)
            ax_r.scatter(angles[:-1], vals[:-1],
                         color=c, s=55, zorder=5, edgecolors="white", linewidths=0.8)

        if not baseline_df.empty:
            bvals = _scaler_vals.iloc[-1].values.flatten().tolist() + [_scaler_vals.iloc[-1].values[0]]
            ax_r.plot(angles, bvals, linewidth=2, linestyle="--", color="#999",
                      label=f"{CFG.TIER_DECOY} avg.  (reference)", zorder=3)
            ax_r.fill(angles, bvals, color="#999", alpha=0.05, zorder=2)
            ax_r.scatter(angles[:-1], bvals[:-1],
                         color="#999", s=45, zorder=5, edgecolors="white", linewidths=0.8)

        _radar_fs = 8
        ax_r.set_xticks(angles[:-1])
        ax_r.set_xticklabels([])   # clear default labels; draw manually below
        ax_r.tick_params(axis="x", pad=18)
        # Draw spoke labels facing outward from the centre
        for _ri_s, (_ang_s, _lbl_s) in enumerate(zip(angles[:-1], spoke_labels)):
            _ang_deg_s = np.degrees(_ang_s)
            if 85 < _ang_deg_s < 95 or 265 < _ang_deg_s < 275:
                _ha_s = "center"
            elif _ang_deg_s < 180:
                _ha_s = "left"
            else:
                _ha_s = "right"
            _rot_s = _ang_deg_s - 90
            if _ang_deg_s > 180:
                _rot_s += 180
            ax_r.text(_ang_s, ax_r.get_ylim()[1] * 1.15, _lbl_s,
                      ha=_ha_s, va="center", fontsize=_radar_fs,
                      rotation=_rot_s, rotation_mode="anchor")
        ax_r.set_yticks(_ring_levels)
        ax_r.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"],
                             color="#CC2222", size=7.5, fontweight="bold")
        ax_r.set_rlabel_position(45)
        ax_r.set_ylim(0, 1.10)
        ax_r.yaxis.grid(False)
        ax_r.xaxis.grid(True, color="#CCCCCC", linewidth=0.6, alpha=0.7)

        n_lines = len(rows_df) + (0 if baseline_df.empty else 1)
        ax_r.legend(loc="lower center", bbox_to_anchor=(0.5, -0.18),
                    ncol=min(3, n_lines), fontsize=_radar_fs,
                    framealpha=0.88, fancybox=True)
        plt.tight_layout(rect=[0, 0.15, 1, 1])
        plt.savefig(fig_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    try:
        metrics_r = [m for m in _radar_labels if m in df.columns]
        if metrics_r:
            worst_tier = existing_tiers[-1] if existing_tiers else None
            worst_avg  = (df[df["degrader_tier"] == worst_tier].mean(numeric_only=True).to_frame().T
                          if worst_tier else pd.DataFrame())

            # --- 18a: All Tier_1A complexes sorted by Scientific_Rank ---
            _pa18a = df[df["degrader_tier"] == CFG.TIER_TOP].copy()
            _sort_col18a = "Scientific_Rank" if "Scientific_Rank" in _pa18a.columns \
                           else "Boltz_Model_Confidence"
            _asc18a = _sort_col18a == "Scientific_Rank"
            _pa18a  = _pa18a.sort_values(_sort_col18a, ascending=_asc18a)
            if _pa18a.empty:
                _pa18a = df.sort_values("Boltz_Model_Confidence", ascending=False).head(8)
            # Cap plotted hits so the radar stays readable (top-N by ranking key).
            _pa18a = _pa18a.head(CFG.VIS_RADAR_MAX_HITS)
            _draw_radar(_pa18a, worst_avg, metrics_r,
                        out_dir / "Figure_19a_Radar_TopHits.png",
                        f"All {CFG.TIER_TOP} complexes vs. worst-tier baseline",
                        label_mode="rank")

            # --- 18b: One representative per tier (best Boltz confidence per tier) ---
            reps_b = []
            for t in existing_tiers:
                sub_t = df[df["degrader_tier"] == t]
                if not sub_t.empty:
                    reps_b.append(sub_t.sort_values("Boltz_Model_Confidence", ascending=False).iloc[[0]])
            rep_df = pd.concat(reps_b) if reps_b else df.head(len(existing_tiers))
            _draw_radar(rep_df, pd.DataFrame(), metrics_r,
                        out_dir / "Figure_19b_Radar_TierReps.png",
                        f"One best representative per tier ({CFG.TIER_TOP} → worst)",
                        label_mode="tier")
    except Exception as e:
        reporter.log(f"  ! Radar chart skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 20: Dual-Objective Tier Assessment — SN2 Geometry vs AI Confidence ---
    """
    Replaces the broken Pareto-rank approach (Binding_Probability was ~constant 0.99).
    New design: for each tier, compute % meeting SN2 ≥ 165° (geometry) and
    % meeting Boltz_Model_Confidence ≥ 0.75 (AI) — showing where each tier sits
    in the biologically meaningful substrate-vs-inhibitor space.
    """
    if "SN2_Attack_Angle" in df.columns and "Boltz_Model_Confidence" in df.columns:
        try:

            _d19 = df[["SN2_Attack_Angle", "Boltz_Model_Confidence", "degrader_tier"]].copy()
            _d19["SN2_Attack_Angle"]      = pd.to_numeric(_d19["SN2_Attack_Angle"],      errors="coerce")
            _d19["Boltz_Model_Confidence"] = pd.to_numeric(_d19["Boltz_Model_Confidence"], errors="coerce")
            _d19 = _d19.dropna()
            _tiers19 = [t for t in TIER_ORDER_LOGIC if t in _d19["degrader_tier"].unique()]

            # Per-tier metrics
            _rows19 = []
            for _t19 in _tiers19:
                _sub = _d19[_d19["degrader_tier"] == _t19]
                _n   = len(_sub)
                _pct_sn2  = float((_sub["SN2_Attack_Angle"] >= CFG.SUBSTRATE_ANGLE_MIN).mean() * 100)
                _pct_conf = float((_sub["Boltz_Model_Confidence"] >= CFG.SUBSTRATE_CONF_MIN).mean() * 100)
                _pct_both = float(((_sub["SN2_Attack_Angle"] >= CFG.SUBSTRATE_ANGLE_MIN) & (_sub["Boltz_Model_Confidence"] >= CFG.SUBSTRATE_CONF_MIN)).mean() * 100)
                _pct_inh  = float((_sub["SN2_Attack_Angle"] < CFG.INHIBITOR_ANGLE_MAX).mean() * 100)
                _med_sn2  = float(_sub["SN2_Attack_Angle"].median())
                _med_conf = float(_sub["Boltz_Model_Confidence"].median())
                _rows19.append({"tier": _t19, "n": _n,
                                "pct_substrate": _pct_both,
                                "pct_sn2_ok": _pct_sn2,
                                "pct_conf_ok": _pct_conf,
                                "pct_inhibitor": _pct_inh,
                                "med_sn2": _med_sn2,
                                "med_conf": _med_conf})
            _df19 = pd.DataFrame(_rows19)

            """
            Figure 20 split into two standalone panels: 19a (tier % success bars)
            and 19b (Conf × SN2 median scatter). Each saved as its own PNG.
            """
            fig19a, ax19L = plt.subplots(figsize=(7.5, 7))
            fig19b, ax19R = plt.subplots(figsize=(7.5, 7))

            # ── Panel 19a: tier % success (substrate vs inhibitor) ───────────
            _y19 = np.arange(len(_df19))
            _bar_h = 0.25
            """
            Master container bar per tier — grey track, outline coloured to match
            the tier's y-axis label colour. Shorter height leaves clear space
            between adjacent tier groups; outline 70% of previous thickness.
            """
            from matplotlib.colors import to_rgba as _to_rgba19
            _master_ec19 = [TIER_PALETTE.get(t, "#BDBDBD") for t in _df19["tier"].tolist()]
            """
            Faint tint of each tier's own colour as the track fill (alpha baked in
            so the coloured outline stays vivid) — cohesive, not flat grey.
            """
            _master_fc19 = [_to_rgba19(c, 0.13) for c in _master_ec19]
            ax19L.barh(_y19, 105, height=0.82, color=_master_fc19,
                       edgecolor=_master_ec19, linewidth=1.4, zorder=2)
            ax19L.barh(_y19 + _bar_h,  _df19["pct_substrate"].values,  height=_bar_h,
                       color="#009E73", alpha=0.85, label=f"Substrate (SN2≥{CFG.SUBSTRATE_ANGLE_MIN:.0f}° + Conf≥{CFG.SUBSTRATE_CONF_MIN:.2f})", zorder=3)
            ax19L.barh(_y19,            _df19["pct_sn2_ok"].values,     height=_bar_h,
                       color="#56B4E9", alpha=0.75, label=f"SN2≥{CFG.SUBSTRATE_ANGLE_MIN:.0f}° (any conf)", zorder=3)
            ax19L.barh(_y19 - _bar_h,  _df19["pct_inhibitor"].values,  height=_bar_h,
                       color="#D55E00", alpha=0.80, label=f"Potential inhibitor (SN2<{CFG.INHIBITOR_ANGLE_MAX:.0f}°)", zorder=3)

            """
            Value labels — inside (centred, white) when the bar is wide enough,
            else just outside the bar's right end (dark text).
            """
            for _i19, row19 in _df19.iterrows():
                _yi = int(_i19)
                for _pct_val, _yo in [(row19["pct_substrate"], _bar_h),
                                      (row19["pct_sn2_ok"], 0),
                                      (row19["pct_inhibitor"], -_bar_h)]:
                    if _pct_val < 1.0:
                        continue
                    if _pct_val >= 12:
                        ax19L.text(_pct_val / 2, _yi + _yo, f"{_pct_val:.0f}%",
                                   va="center", ha="center", fontsize=6,
                                   color="white", fontweight="bold", zorder=5)
                    else:
                        ax19L.text(_pct_val + 1.0, _yi + _yo, f"{_pct_val:.0f}%",
                                   va="center", ha="left", fontsize=6,
                                   color="#333", zorder=5)

            ax19L.set_yticks(_y19)
            ax19L.set_yticklabels(
                [f'{r["tier"]}  (n={r["n"]:,})' for _, r in _df19.iterrows()],
                fontsize=7
            )
            for tick, tier in zip(ax19L.get_yticklabels(), _df19["tier"].tolist()):
                tick.set_color(TIER_PALETTE.get(tier, "black"))
                tick.set_fontweight("bold")
            ax19L.set_xlim(0, 105)
            ax19L.axvline(25, color="#999", lw=0.7, ls=":", zorder=2.4)
            ax19L.axvline(50, color="#999", lw=0.7, ls=":", zorder=2.4)
            ax19L.axvline(75, color="#999", lw=0.7, ls=":", zorder=2.4)
            ax19L.set_xlabel("Percentage of complexes (%)", fontsize=9)
            ax19L.tick_params(axis="x", labelsize=8)
            # Legend: single row, 70% font (7.5→5.25), inside bottom-right
            ax19L.legend(loc="lower right", ncol=3, fontsize=5.25,
                         framealpha=0.92, fancybox=True, borderpad=0.4,
                         handlelength=1.0, handletextpad=0.3, columnspacing=0.7)
            ax19L.invert_yaxis()

            # ── Right panel: scatter of tier medians in Conf × SN2 space ─────
            # Background hexbin of all complexes (shows density landscape)
            _hb19 = ax19R.hexbin(_d19["Boltz_Model_Confidence"], _d19["SN2_Attack_Angle"],
                                  gridsize=40, cmap="Greys", mincnt=1, alpha=0.4, linewidths=0.2, zorder=1)

            # Substrate / inhibitor zone shading — bounds from CFG
            ax19R.axhspan(CFG.SUBSTRATE_ANGLE_MIN, 185, color="#009E73", alpha=0.10, zorder=0)
            ax19R.axhspan(0, CFG.INHIBITOR_ANGLE_MAX, color="#D55E00", alpha=0.08, zorder=0)
            ax19R.axhline(CFG.SUBSTRATE_ANGLE_MIN, color="#009E73", lw=1.2, ls="--", alpha=0.7, zorder=2)
            ax19R.axhline(CFG.INHIBITOR_ANGLE_MAX, color="#D55E00", lw=1.2, ls="--", alpha=0.7, zorder=2)
            ax19R.axvline(CFG.SUBSTRATE_CONF_MIN, color="#0072B2", lw=1.2, ls="--", alpha=0.7, zorder=2)

            # Tier median diamonds with IQR error bars — collect positions first
            _tier_pts19 = []   # (tier, cx, cy, colour)
            for _t19 in _tiers19:
                _sub = _d19[_d19["degrader_tier"] == _t19]
                if _sub.empty:
                    continue
                _cx = float(_sub["Boltz_Model_Confidence"].median())
                _cy = float(_sub["SN2_Attack_Angle"].median())
                _ex_lo = _cx - float(_sub["Boltz_Model_Confidence"].quantile(0.25))
                _ex_hi = float(_sub["Boltz_Model_Confidence"].quantile(0.75)) - _cx
                _ey_lo = _cy - float(_sub["SN2_Attack_Angle"].quantile(0.25))
                _ey_hi = float(_sub["SN2_Attack_Angle"].quantile(0.75)) - _cy
                _col19 = TIER_PALETTE.get(_t19, "#999")
                ax19R.errorbar(_cx, _cy,
                               xerr=[[_ex_lo], [_ex_hi]],
                               yerr=[[_ey_lo], [_ey_hi]],
                               fmt="D", color=_col19, ms=11,
                               ecolor=_col19, elinewidth=1.5, capsize=3,
                               markeredgecolor="black", markeredgewidth=0.8,
                               zorder=6, label=f"{_t19}  (med SN2={_cy:.0f}°, conf={_cx:.3f})")
                _tier_pts19.append((_t19, _cx, _cy, _col19))

            ax19R.set_xlim(0.72, 1.01)
            ax19R.set_ylim(-5, 185)

            """
            Tier name labels parked in the empty left column, evenly spread on y so
            none overlap; gentle curved arrows point from each label to its diamond.
            """
            _pts_sorted19 = sorted(_tier_pts19, key=lambda p: p[2], reverse=True)
            _n_lab19 = len(_pts_sorted19)
            _lab_ys19 = list(np.linspace(160, 25, _n_lab19)) if _n_lab19 > 1 else [90]
            for (_t19, _cx, _cy, _col19), _laby in zip(_pts_sorted19, _lab_ys19):
                ax19R.annotate(_t19, xy=(_cx, _cy),
                               xytext=(0.83, _laby),
                               ha="right", va="center",
                               fontsize=6.5, color=_col19, fontweight="bold",
                               arrowprops=dict(arrowstyle="-|>", color=_col19, lw=1.3,
                                               mutation_scale=15, shrinkA=3, shrinkB=5,
                                               connectionstyle="arc3,rad=0.12"),
                               bbox=dict(boxstyle="round,pad=0.18", fc="white",
                                         ec=_col19, alpha=0.92, linewidth=0.8),
                               zorder=8)

            ax19R.set_xlabel("Boltz Model Confidence", fontsize=9)
            ax19R.set_ylabel("SN2 Attack Angle (°)", fontsize=9)
            ax19R.tick_params(axis="both", labelsize=8)

            # Zone / threshold labels — all inside the axes, single line each
            ax19R.text(0.87, 181, "Substrate zone", color="#007A52", fontsize=7,
                       style="italic", ha="center", va="center")
            ax19R.text(0.87, 8, "Inhibitor zone", color="#C04000", fontsize=7,
                       style="italic", ha="center", va="center")
            ax19R.text(0.745, 95, f"Conf threshold = {CFG.SUBSTRATE_CONF_MIN:.2f}", color="#0072B2", fontsize=6.5,
                       ha="center", va="center", style="italic", rotation=90)

            # Legend: 3 columns, shrunk icons, fixed outside above the top-left corner
            ax19R.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02),
                         ncol=3, fontsize=5.25, framealpha=0.90, fancybox=True,
                         markerscale=0.45, handlelength=1.0, handletextpad=0.3,
                         columnspacing=0.8, labelspacing=0.3, borderpad=0.3)

            fig19a.savefig(out_dir / "Figure_20a_Tier_Success_Rates.png",
                           dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            fig19b.savefig(out_dir / "Figure_20b_Conf_SN2_Landscape.png",
                           dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close(fig19a)
            plt.close(fig19b)
        except Exception as e:
            reporter.log(f"  ! Figure 20 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 21: Conflict Composition — 100% Stacked Horizontal Bar per Tier ---
    """
    One bar per tier (best → worst on y-axis).  Left portion = favourable categories
    (Consensus High, Hidden Gem); right portion = concerning categories.
    A thin vertical divider at the "favourable total" position separates the two sides.
    Segment labels show % only for segments ≥ 5%.  Y-axis labels include tier n-count.
    """
    if "Conflict_Category" in df.columns and "degrader_tier" in df.columns:
        cat_colours_f2 = {
            "Consensus High": "#009E73", "Hidden Gem": "#CC79A7",
            CFG.TIER_DECOY: "#D55E00", "Consensus Low": "#777777", "Ambiguous": "#E69F00"}
        # Canonical plotting order: favourable → concerning, left → right
        cat_order_f2 = ["Consensus High", "Hidden Gem", "Ambiguous", "Consensus Low", "Decoy"]
        present_cats_f2 = [c for c in cat_order_f2 if c in df["Conflict_Category"].unique()]

        ct_f2 = pd.crosstab(df["degrader_tier"], df["Conflict_Category"])
        tier_rows_f2 = [t for t in existing_tiers if t in ct_f2.index]
        ct_f2 = ct_f2.reindex(index=tier_rows_f2, columns=present_cats_f2, fill_value=0)
        tier_n_f2   = ct_f2.sum(axis=1)
        # Grand total = all rows with a valid degrader_tier (not filtered by Conflict_Category)
        _total_rows_f2 = len(df)
        grand_total    = int(df["degrader_tier"].notna().sum())
        pct_f2         = ct_f2.div(tier_n_f2, axis=0).fillna(0) * 100
        # Detect control cases — generic mask first, then specific 3R3U / DeHa4
        _ctrl_mask = pd.Series(False, index=df.index)
        for _ctrl_col in ["is_control", "control", "Control", "is_Control"]:
            if _ctrl_col in df.columns:
                _ctrl_mask = df[_ctrl_col].astype(bool)
                break
        if not _ctrl_mask.any() and "job_name" in df.columns:
            _ctrl_mask = df["job_name"].str.lower().str.contains("control|ctrl|deha4_ref|3r3u", na=False)
        _n_controls = int(_ctrl_mask.sum())
        """
        Identify the two reference controls by their exact job-name patterns:
        3R3U  control: job_name contains '3r3u'  AND 'control'
        DeHa4 control: job_name contains 'deha4' AND 'control'
        E.g.: "0000000_4_3R3U_Control_26_Fluoroacetate"
        "0000000_02_DeHa4_Control_26_Fluoroacetate"
        """
        _ctrl_3R3U_mask  = pd.Series(False, index=df.index)
        _ctrl_deha4_mask = pd.Series(False, index=df.index)
        if "job_name" in df.columns:
            _jn_lower = df["job_name"].str.lower()
            _ctrl_3R3U_mask = (
                _jn_lower.str.contains("3r3u", na=False) &
                _jn_lower.str.contains("control|ctrl", na=False)
            )
            _ctrl_deha4_mask = (
                _jn_lower.str.contains("deha4", na=False) &
                _jn_lower.str.contains("control|ctrl", na=False)
            )
            # Fallback: if 'control' keyword absent, match on 3r3u / deha4 alone
            if not _ctrl_3R3U_mask.any():
                _ctrl_3R3U_mask = _jn_lower.str.contains("3r3u", na=False)
            if not _ctrl_deha4_mask.any():
                _ctrl_deha4_mask = _jn_lower.str.contains("deha4", na=False)
        _n_3R3U  = int(_ctrl_3R3U_mask.sum())
        _n_deha4 = int(_ctrl_deha4_mask.sum())

        # Y positions (top tier at top)
        tier_y_f2   = np.arange(len(tier_rows_f2) - 1, -1, -1)
        bar_h       = 0.60

        # Height: fit bars tightly — no whitespace padding below lowest tier
        _fh_f2 = max(4.0, len(tier_rows_f2) * 0.65 + 1.8)
        fig_f2, ax_f2 = plt.subplots(figsize=(13, _fh_f2))
        ax_f2.set_facecolor("#FAFAFA")

        left_f2 = np.zeros(len(tier_rows_f2))
        for cat in present_cats_f2:
            vals   = pct_f2[cat].values
            counts = ct_f2[cat].values
            ax_f2.barh(tier_y_f2, vals, left=left_f2, height=bar_h,
                       color=cat_colours_f2.get(cat, "#999"),
                       edgecolor="white", linewidth=0.7, zorder=2)
            for yi, vi, li, ni in zip(tier_y_f2, vals, left_f2, counts):
                if vi < 4:
                    continue
                txt_col = "white" if vi > 12 else "black"
                if vi >= 9:
                    lbl = f"{vi:.0f}%\n(n={int(ni):,})"
                    fs  = 7.0
                else:
                    lbl = f"{vi:.0f}%"
                    fs  = 7.5
                ax_f2.text(li + vi / 2, yi, lbl,
                           ha="center", va="center", fontsize=fs,
                           color=txt_col, fontweight="bold", zorder=4,
                           linespacing=1.1)
            left_f2 = left_f2 + vals

        # Vertical divider at boundary between favourable and concerning categories
        fav_cats_f2 = [c for c in ["Consensus High", "Hidden Gem"] if c in present_cats_f2]
        fav_pct_f2  = pct_f2[fav_cats_f2].sum(axis=1).values if fav_cats_f2 else None
        if fav_pct_f2 is not None:
            for yi, xv in zip(tier_y_f2, fav_pct_f2):
                ax_f2.plot([xv, xv], [yi - bar_h / 2, yi + bar_h / 2],
                           color="white", linewidth=2.2, zorder=5)

        # Reference line at 50%
        ax_f2.axvline(50, color="#555555", linewidth=1.1, linestyle="--", alpha=0.5, zorder=1)
        ax_f2.text(50.5, len(tier_rows_f2) - 0.5, "50%",
                   va="top", ha="left", fontsize=8, color="#555")

        # Y-axis labels: tier name + n-count in matching colour
        ax_f2.set_yticks(tier_y_f2)
        ax_f2.set_yticklabels(
            [f"{t}  (n={tier_n_f2.loc[t]:,})" for t in tier_rows_f2],
            fontsize=10)
        for tick, tier in zip(ax_f2.get_yticklabels(), tier_rows_f2):
            tick.set_color(TIER_PALETTE.get(tier, "black"))
            tick.set_fontweight("bold")

        ax_f2.set_xlim(0, 100)
        # Clamp y-range to exactly the number of tiers — prevents blank canvas below
        ax_f2.set_ylim(-0.55, len(tier_rows_f2) - 0.45)
        ax_f2.set_xlabel("Percentage of complexes within tier  (%)", fontsize=11)
        ax_f2.set_xticks([0, 25, 50, 75, 100])
        ax_f2.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
        ax_f2.set_axisbelow(True)
        ax_f2.xaxis.grid(True, alpha=0.30, linestyle=":", color="#888", zorder=0)
        #ax_f2.set_title("Conflict-Category Composition by Degrader Tier", fontsize=13, pad=10)

        # Top-right badge: total jobs | tiered count — lifted above the axes frame
        ax_f2.text(0.99, 1.06,
                   f"Total jobs: {_total_rows_f2:,}   |   Tiered: {grand_total:,}",
                   transform=ax_f2.transAxes, ha="right", va="bottom",
                   fontsize=8.5, color="#333", fontweight="bold",
                   bbox=dict(boxstyle="round,pad=0.18", fc="white", ec="#ccc",
                             alpha=0.88, linewidth=0.5))

        """
        Section labels — placed INSIDE the chart at the top of the bars so they do not
        collide with the "Total jobs | Tiered" badge that floats above the axes frame.
        get_xaxis_transform(): x = data coords (0–100 %), y = axes fraction [0,1].
        """
        if fav_cats_f2:
            mid_fav = np.mean(fav_pct_f2) / 2
            mid_con = np.mean(fav_pct_f2) + (100 - np.mean(fav_pct_f2)) / 2
            ax_f2.text(mid_fav, 1.015, "▶ Favourable",
                       ha="center", va="bottom", fontsize=9.5,
                       color="#009E73", fontweight="bold",
                       transform=ax_f2.get_xaxis_transform(),
                       bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=CFG.CONF_BAND_COLOURS["high"],
                                 alpha=0.80, linewidth=0.5))
            ax_f2.text(mid_con, 1.015, "Concerning ◀",
                       ha="center", va="bottom", fontsize=9.5,
                       color="#D55E00", fontweight="bold",
                       transform=ax_f2.get_xaxis_transform(),
                       bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=CFG.CONF_BAND_COLOURS["below"],
                                 alpha=0.80, linewidth=0.5))

        # Arrow annotations for small segments (< 4%) that couldn't fit inline text
        left_f2_tmp = np.zeros(len(tier_rows_f2))
        for cat in present_cats_f2:
            vals   = pct_f2[cat].values
            counts = ct_f2[cat].values
            for yi, vi, li, ni in zip(tier_y_f2, vals, left_f2_tmp, counts):
                if vi > 0 and vi < 4 and ni > 0:
                    mid_x = li + vi / 2
                    col_arr = cat_colours_f2.get(cat, "#999")
                    ax_f2.annotate(f"n={int(ni):,}",
                                   xy=(mid_x, yi),
                                   xytext=(mid_x, yi + 0.38),
                                   ha="center", va="bottom", fontsize=6.5,
                                   color=col_arr, fontweight="bold",
                                   arrowprops=dict(arrowstyle="->", color=col_arr,
                                                   lw=0.7, shrinkB=2),
                                   bbox=dict(boxstyle="round,pad=0.10", fc="white",
                                             ec=col_arr, alpha=0.85, linewidth=0.5),
                                   zorder=6)
            left_f2_tmp = left_f2_tmp + vals

        """
        Highlight each control complex individually — labelled protein+ligand
        (e.g. "DeHa4+TFA"). 6 controls: {3R3U, DeHa4} × {TFA, Fluoroacetate,
        Difluoroacetate}. Stacked vertically per tier row so they never overlap.
        """
        def _lig_short_f20(_idx):
            _ln = ""
            if "Ligand_Name" in df.columns:
                _ln = str(df.at[_idx, "Ligand_Name"])
            if (not _ln or _ln.lower() == "nan") and "job_name" in df.columns:
                _jn = str(df.at[_idx, "job_name"])
                _ln = _jn.split("Control_", 1)[1] if "Control_" in _jn else _jn
            # strip a leading numeric index like "25_TFA" → "TFA"
            if "_" in _ln and _ln.split("_", 1)[0].isdigit():
                _ln = _ln.split("_", 1)[1]
            return _ln.replace("_", " ").strip()

        _ctrl_specs_f20 = [(_ctrl_3R3U_mask, "3R3U", "#8B0000"),
                           (_ctrl_deha4_mask, "DeHa4", "#00408B")]
        _ctrl_by_tier_f20 = {}   # tier -> [(label, colour), ...]
        for _cmask, _cprot, _ccol in _ctrl_specs_f20:
            if not (_cmask.any() and "job_name" in df.columns):
                continue
            for _cidx in df.index[_cmask]:
                _ctier = df.at[_cidx, "degrader_tier"]
                if _ctier not in tier_rows_f2:
                    continue
                _ctrl_by_tier_f20.setdefault(_ctier, []).append(
                    (f"{_cprot}+{_lig_short_f20(_cidx)}", _ccol))
        for _ctier, _items in _ctrl_by_tier_f20.items():
            _cy = tier_y_f2[tier_rows_f2.index(_ctier)]
            _k  = len(_items)
            _offs = np.linspace(0.34, -0.34, _k) if _k > 1 else [0.0]
            for (_clbl, _ccol), _off in zip(_items, _offs):
                ax_f2.annotate(f"◀ {_clbl}", xy=(100, _cy),
                               xytext=(99, _cy + _off),
                               ha="right", va="center", fontsize=6.8,
                               color=_ccol, fontweight="bold",
                               arrowprops=dict(arrowstyle="->", color=_ccol, lw=0.9,
                                               shrinkA=1, shrinkB=2),
                               bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                         ec=_ccol, alpha=0.90, linewidth=0.5),
                               zorder=7)

        from matplotlib.patches import Patch as _P02
        _leg02 = [_P02(facecolor=cat_colours_f2.get(c, "#999"), edgecolor="white",
                       linewidth=0.5, label=c) for c in present_cats_f2]
        # Single-row legend below chart — ncol = number of categories so all fit on one line
        _leg02_obj = ax_f2.legend(handles=_leg02,
                                  loc="upper center", bbox_to_anchor=(0.5, -0.10),
                                  ncol=len(present_cats_f2),
                                  fontsize=8.5, framealpha=0.92, fancybox=True,
                                  title="Conflict Category", title_fontsize=8)
        _leg02_obj.set_zorder(20)

        fig_f2.tight_layout(rect=[0, 0.07, 1, 1])
        plt.savefig(out_dir / "Figure_21_Conflict_Composition.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # --- Figure 22: Hidden Gems — Slope Chart (SN2 Angle vs AI Confidence) ---
    """
    Each Hidden Gem complex is a lollipop row showing its absolute SN2 attack angle
    (physics metric, left dot) and Boltz Model Confidence (AI metric, right dot).
    A large gap means strong mechanistic geometry but undervalued by AI — the
    "hidden gem" story in a single visual stroke.
    """
    _phys_col21_abs = next((c for c in ["SN2_Attack_Angle", "Scientific_Rank", "Pareto_Rank"] if c in df.columns), None)
    _ai_col21_abs   = next((c for c in ["Boltz_Model_Confidence", "Ensemble_Data_Rank"] if c in df.columns), None)
    if "Conflict_Category" in df.columns and _phys_col21_abs and _ai_col21_abs:
        try:
            gems = df[df["Conflict_Category"] == "Hidden Gem"].copy()
            if not gems.empty:
                # Use absolute values for physics and AI metrics
                gems["_phys_pct21"] = pd.to_numeric(gems[_phys_col21_abs], errors="coerce")
                gems["_ai_pct21"]   = pd.to_numeric(gems[_ai_col21_abs],   errors="coerce")
                gems = gems.dropna(subset=["_phys_pct21", "_ai_pct21"])
                # Normalise to 0–100 scale for consistent axis display
                _phys_min21, _phys_max21 = float(gems["_phys_pct21"].min()), float(gems["_phys_pct21"].max())
                _ai_min21,   _ai_max21   = float(gems["_ai_pct21"].min()),   float(gems["_ai_pct21"].max())
                _phys_rng21 = _phys_max21 - _phys_min21 if _phys_max21 > _phys_min21 else 1.0
                _ai_rng21   = _ai_max21   - _ai_min21   if _ai_max21   > _ai_min21   else 1.0
                gems["_phys_pct21"] = (gems["_phys_pct21"] - _phys_min21) / _phys_rng21 * 100
                gems["_ai_pct21"]   = (gems["_ai_pct21"]   - _ai_min21)   / _ai_rng21   * 100
                gems["_gap21"]      = gems["_phys_pct21"] - gems["_ai_pct21"]

                lname_col21 = next((c for c in ["complex_id", "Ligand_Name", "ligand",
                                                "Protein_Name", "protein"]
                                    if c in gems.columns), None)

                # ── Horizontal dual-dot lollipop chart ──────────────────────────
                """
                Each Hidden Gem is a row.  Two coloured dots show its Physics rank
                (blue circle) and AI rank (orange diamond) on a shared x-axis (0–100%).
                A grey connecting line spans between the two dots; width ∝ gap.
                Sorted by AI rank (best at top) so the "most overlooked" gems sit
                at the bottom where they are easy to read.
                """
                from matplotlib.lines import Line2D as _L21
                gems_s21 = gems.sort_values("_ai_pct21", ascending=True)  # best at top after invert
                _ng21 = len(gems_s21)
                _y_pos21 = np.arange(_ng21)

                # Dynamic x-axis range so clustered low-percentile gems fill the axis
                _all_pct21 = (list(gems_s21["_phys_pct21"].values) +
                              list(gems_s21["_ai_pct21"].values))
                _x_data_max21 = max(_all_pct21) if _all_pct21 else 5.0
                _x_hi21 = float(np.clip(_x_data_max21 * 4, 20, 100))

                fig21, ax21 = plt.subplots(figsize=(12, max(4.5, _ng21 * 0.70 + 2.2)))
                ax21.set_facecolor("#FAFAFA")
                ax21.set_axisbelow(True)
                ax21.xaxis.grid(True, color="#DCDCDC", linewidth=0.7, zorder=0)

                # Background zone bands — tertile reference lines
                ax21.axvspan(0,  33, alpha=0.04, color="#D55E00", zorder=0)  # lower third
                ax21.axvspan(33, 67, alpha=0.03, color="#E69F00", zorder=0)  # middle
                ax21.axvspan(67, 100, alpha=0.04, color="#009E73", zorder=0)  # top third
                for _xv21 in [33, 67]:
                    ax21.axvline(_xv21, color="#CCCCCC", linewidth=0.8,
                                 linestyle="--", zorder=1)

                for _gi21, (_, row21) in enumerate(gems_s21.iterrows()):
                    _phys21 = float(row21["_phys_pct21"])
                    _ai21   = float(row21["_ai_pct21"])
                    _gap21  = float(row21["_gap21"])
                    _tier21 = (row21.get("degrader_tier", "")
                               if "degrader_tier" in gems.columns else "")
                    _col21  = TIER_PALETTE.get(_tier21, "#AA3377")
                    _y21    = _gi21

                    # Connecting line — width proportional to |gap|
                    _lw21 = np.clip(0.8 + abs(_gap21) / 25, 0.8, 3.5)
                    _lo_x21, _hi_x21 = sorted([_phys21, _ai21])
                    ax21.hlines(_y21, _lo_x21, _hi_x21,
                                colors="#CCCCCC", linewidth=_lw21,
                                alpha=0.75, zorder=2)

                    # Physics dot (blue circle, red outline)
                    ax21.scatter(_phys21, _y21, color="#0072B2", s=130, zorder=4,
                                 marker="o", edgecolors="red", linewidths=1.2)
                    ax21.text(_phys21, _y21 + 0.34, f"{_phys21:.0f}%",
                              va="bottom", ha="center", fontsize=5.5, color="#0072B2",
                              fontweight="bold", zorder=5)
                    # AI dot (golden diamond, no border)
                    ax21.scatter(_ai21, _y21, color="#E69F00", s=130, zorder=4,
                                 marker="D", edgecolors="none")
                    ax21.text(_ai21, _y21 - 0.34, f"{_ai21:.0f}%",
                              va="top", ha="center", fontsize=5.5, color="#E69F00",
                              fontweight="bold", zorder=5)

                    # Value labels just outside the right dot
                    _rightmost21 = max(_phys21, _ai21)
                    ax21.text(_rightmost21 + 1.2, _y21,
                              f"Δ={abs(_gap21):.1f}%",
                              va="center", ha="left", fontsize=7.5,
                              color="#555", fontweight="bold")

                # Y-axis: gem names (reversed so top row = best AI rank)
                _gem_labels21 = []
                for _gi_lbl21, (_, row21) in enumerate(gems_s21.iterrows()):
                    if lname_col21:
                        _gem_labels21.append(str(row21[lname_col21])[:28])
                    else:
                        _gem_labels21.append(f"Hidden Gem {_gi_lbl21 + 1}")
                ax21.set_yticks(_y_pos21)
                ax21.set_yticklabels(_gem_labels21, fontsize=8.5)
                for _tick21, (_, row21) in zip(ax21.get_yticklabels(), gems_s21.iterrows()):
                    _t21 = (row21.get("degrader_tier", "")
                            if "degrader_tier" in gems.columns else "")
                    _tick21.set_color(TIER_PALETTE.get(_t21, "#333"))
                    _tick21.set_fontweight("bold")

                ax21.set_xlim(-1, _x_hi21 + 2)
                _x21_phys_lbl = ("SN2 Attack Angle" if _phys_col21_abs == "SN2_Attack_Angle"
                                 else "Physics Rank")
                _x21_ai_lbl   = ("Boltz Confidence" if _ai_col21_abs == "Boltz_Model_Confidence"
                                 else "AI Rank")
                ax21.set_xlabel(
                    f"Normalised Score  ({_x21_phys_lbl} / {_x21_ai_lbl})  "
                    "(0 = lowest · 100 = highest)",
                    fontsize=10)
                _all_tick_vals21 = [0, 5, 10, 15, 20, 25, 33, 50, 67, 80, 90, 100]
                _tick_lbl_map21  = {0: "0%", 5: "5%", 10: "10%", 15: "15%",
                                    20: "20%", 25: "25%", 33: "33rd", 50: "50%",
                                    67: "67th", 80: "80%", 90: "90%", 100: "100%"}
                _vis_ticks21 = [t for t in _all_tick_vals21 if t <= _x_hi21 + 1]
                ax21.set_xticks(_vis_ticks21)
                ax21.set_xticklabels([_tick_lbl_map21[t] for t in _vis_ticks21], fontsize=8.5)
                ax21.set_ylim(-0.6, _ng21 - 0.4)

                # Zone labels at top — only render those within the visible x-range
                _zone_lbls21 = [
                    (16.5,  "Lower third", "#D55E00", 33),
                    (50.0,  "Middle third", "#8A6000", 67),
                    (83.5,  "Top third",   "#009E73", 100),
                ]
                for _zx21, _ztxt21, _zcol21, _zmax21 in _zone_lbls21:
                    if _zx21 <= _x_hi21 + 1:
                        ax21.text(_zx21, 0.97, _ztxt21, ha="center", va="top",
                                  fontsize=7.5, color=_zcol21, style="italic",
                                  transform=ax21.get_xaxis_transform())

                # Legend
                _lh21 = [
                    _L21([0],[0], marker="o", color="w", markerfacecolor="#0072B2",
                         markeredgecolor="red", markeredgewidth=1.2,
                         markersize=9, label=f"Physics score ({_x21_phys_lbl})"),
                    _L21([0],[0], marker="D", color="w", markerfacecolor="#E69F00",
                         markeredgecolor="none",
                         markersize=9, label=f"AI score ({_x21_ai_lbl})"),
                ]
                # Add tier entries if multiple tiers
                _tier_in_gems = [t for t in existing_tiers
                                 if "degrader_tier" in gems.columns and
                                 t in gems["degrader_tier"].values]
                _lh21 += [_L21([0],[0], color=TIER_PALETTE.get(t,"#999"), linewidth=3,
                                label=f'{t}  (n={len(gems[gems["degrader_tier"]==t]):,})')
                          for t in _tier_in_gems]
                ax21.legend(handles=_lh21, loc="lower right", fontsize=8,
                            framealpha=0.92, fancybox=True, ncol=2,
                            title="Marker type  |  Tier (line colour)",
                            title_fontsize=7.5)

                # Badge: summary stats
                n_gems21 = len(gems)
                med_gap21 = float(gems["_gap21"].median())
                ax21.text(0.99, 0.99,
                          f"Hidden Gems: {n_gems21:,}\nMedian |Δ|: {abs(med_gap21):.1f}%",
                          transform=ax21.transAxes, va="top", ha="right",
                          fontsize=8.5, color="#333",
                          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#ccc",
                                    alpha=0.88, linewidth=0.5))

                plt.tight_layout()
                df.drop(columns=["_phys_pct21", "_ai_pct21"], errors="ignore", inplace=True)
                plt.savefig(out_dir / "Figure_22_Hidden_Gems_DeepDive.png",
                            dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
                plt.close(fig21)
            else:
                reporter.log("  ! Figure 22 skipped: No Hidden Gems found in dataset.")
        except Exception as e:
            reporter.log(f"  ! Figure 22 skipped: {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 23: Category Overlap — Euler / Venn Diagram ---
    try:
        """
        Physics-Strong is a TRUE catalytic-geometry axis, independent of model confidence:
        the top tercile of soft_catalytic_score (continuous SN2-geometry composite),
        falling back to mechanistic_score or SN2 angle. Pareto_Rank is deliberately NOT
        used — it is built from Binding_Probability (near-constant ≈ 0.99) and confidence,
        so it would make the Physics circle a second copy of the AI circle.
        """
        _phys_col22 = next((c for c in ["soft_catalytic_score", "mechanistic_score", "SN2_Attack_Angle"] if c in df.columns), None)
        _ai_col22   = next((c for c in ["Boltz_Model_Confidence", "iptm", "ptm"] if c in df.columns), None)
        _tier_strong22 = {CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2]}

        if _phys_col22 and _ai_col22 and "degrader_tier" in df.columns:
            _n22 = len(df)
            _phys_vals22   = pd.to_numeric(df[_phys_col22], errors="coerce")
            _phys_thresh22 = float(_phys_vals22.quantile(0.67))   # top tercile of catalytic geometry (higher = better)
            _set_A = set(df.index[_phys_vals22 >= _phys_thresh22])
            _ai_thresh22   = float(pd.to_numeric(df[_ai_col22], errors="coerce").quantile(0.67))
            _set_B = set(df.index[pd.to_numeric(df[_ai_col22], errors="coerce") >= _ai_thresh22])
            _set_C = set(df.index[df["degrader_tier"].isin(_tier_strong22)])

            _ctrl_3R3U_idx22  = None
            _ctrl_deha4_idx22 = None
            if "job_name" in df.columns:
                _jn22 = df["job_name"].str.lower()
                _3r3u_m22  = (_jn22.str.contains("3r3u",  na=False) &
                              _jn22.str.contains("control|ctrl", na=False))
                _deha4_m22 = (_jn22.str.contains("deha4", na=False) &
                              _jn22.str.contains("control|ctrl", na=False))
                if _3r3u_m22.any():
                    _ctrl_3R3U_idx22  = df.index[_3r3u_m22][0]
                if _deha4_m22.any():
                    _ctrl_deha4_idx22 = df.index[_deha4_m22][0]

            _only_A   = len(_set_A - _set_B - _set_C)
            _only_B   = len(_set_B - _set_A - _set_C)
            _only_C   = len(_set_C - _set_A - _set_B)
            _AB_not_C = len((_set_A & _set_B) - _set_C)
            _AC_not_B = len((_set_A & _set_C) - _set_B)
            _BC_not_A = len((_set_B & _set_C) - _set_A)
            _ABC      = len(_set_A & _set_B & _set_C)

            from matplotlib.patches import Circle as _Circ22
            from matplotlib.patches import Patch as _Patch22
            from matplotlib.lines import Line2D as _L22

            _, ax22 = plt.subplots(figsize=(11, 9))
            ax22.set_aspect("equal")
            ax22.set_xlim(0.5, 9.5)
            ax22.set_ylim(0.4, 8.7)
            ax22.axis("off")

            _r22 = 2.65
            _cx_A, _cy_A = 3.5, 5.8
            _cx_B, _cy_B = 6.5, 5.8
            _cx_C, _cy_C = 5.0, 3.4
            for _cx22, _cy22, _col22 in [
                    (_cx_A, _cy_A, "#0072B2"),
                    (_cx_B, _cy_B, "#D55E00"),
                    (_cx_C, _cy_C, "#009E73")]:
                ax22.add_patch(_Circ22((_cx22, _cy22), _r22,
                                       facecolor=_col22, alpha=0.20,
                                       edgecolor=_col22, linewidth=2.2, zorder=2))

            _abc_cx, _abc_cy, _r_abc = 5.0, 5.10, 0.75
            ax22.add_patch(_Circ22((_abc_cx, _abc_cy), _r_abc,
                                   facecolor="#AAAAAA", alpha=0.15,
                                   edgecolor="#444", linewidth=2.0,
                                   linestyle=":", zorder=4))

            _set_PA22 = set(df.index[df["degrader_tier"] == CFG.TIER_TOP])
            _n_PA22   = len(_set_PA22)
            _tt_cx22, _tt_cy22, _r_PA22, _tt_col22 = _cx_C, _cy_C, 0.50, "#CC79A7"
            if _n_PA22 > 0:
                _fA = len(_set_PA22 & _set_A) / _n_PA22
                _fB = len(_set_PA22 & _set_B) / _n_PA22
                _fC = len(_set_PA22 & _set_C) / _n_PA22
                _ww = _fA + _fB + _fC
                if _ww > 0:
                    _tt_cx22 = (_fA * _cx_A + _fB * _cx_B + _fC * _cx_C) / _ww
                    _tt_cy22 = (_fA * _cy_A + _fB * _cy_B + _fC * _cy_C) / _ww
                ax22.add_patch(_Circ22((_tt_cx22, _tt_cy22), _r_PA22,
                                       facecolor=_tt_col22, alpha=0.28,
                                       edgecolor=_tt_col22, linewidth=2.5,
                                       linestyle="--", zorder=5))
                ax22.annotate(f"{CFG.TIER_TOP}\n(n={_n_PA22:,})",
                              xy=(_tt_cx22 - _r_PA22 * 0.7, _tt_cy22 + _r_PA22 * 0.7),
                              xytext=(_tt_cx22 - 2.2, _tt_cy22 + 1.4),
                              ha="center", va="bottom", fontsize=8,
                              color=_tt_col22, fontweight="bold", zorder=9,
                              bbox=dict(boxstyle="round,pad=0.20", fc="#F9E4F2",
                                        ec=_tt_col22, alpha=0.93, linewidth=1.3),
                              arrowprops=dict(arrowstyle="->", color=_tt_col22,
                                              lw=1.1, shrinkA=0, shrinkB=3))

            # Primary set labels — uniform size, bold, distinct filled-and-outlined box
            _main_lbl_kw22 = dict(ha="center", va="center", fontsize=11,
                                  fontweight="bold", linespacing=1.25, zorder=6)
            ax22.text(2.15, 7.25, "Physics-Strong", color="#0072B2", **_main_lbl_kw22,
                      bbox=dict(boxstyle="round,pad=0.40", fc="#D6EAF8", ec="#0072B2",
                                alpha=0.95, linewidth=1.8))
            ax22.text(7.45, 7.25, "AI-Strong", color="#D55E00", **_main_lbl_kw22,
                      bbox=dict(boxstyle="round,pad=0.40", fc="#FCE4D0", ec=CFG.CONF_BAND_COLOURS["below"],
                                alpha=0.95, linewidth=1.8))
            ax22.text(5.00, 1.55, "Tier-Strong", color="#009E73", **_main_lbl_kw22,
                      bbox=dict(boxstyle="round,pad=0.40", fc="#D4EFDF", ec=CFG.CONF_BAND_COLOURS["high"],
                                alpha=0.95, linewidth=1.8))

            _region_data22 = [
                (1.85, 5.80, _only_A,    "#0072B2", "Physics only"),
                (8.15, 5.80, _only_B,    "#D55E00", "AI only"),
                (5.00, 2.20, _only_C,    "#009E73", "Tier only"),
                (5.00, 6.90, _AB_not_C,  "#1A4080", "Physics ∩ AI"),
                (2.80, 4.40, _AC_not_B,  "#005840", "Physics ∩ Tier"),
                (7.20, 4.40, _BC_not_A,  "#803000", "AI ∩ Tier"),
                (5.00, 5.30, _ABC,       "#222",    "All three"),
            ]
            """
            Second-layer region labels — lighter square boxes, name normal / count bold,
            visually subordinate to the primary set labels above.
            """
            for _rx22, _ry22, _rn22, _rc22, _rl22 in _region_data22:
                _pct22 = _rn22 / _n22 * 100 if _n22 else 0
                ax22.text(_rx22, _ry22,
                          f"{_rl22}\n" + r"$\bf{" + f"{_rn22:,}" + r"}$" + f"  ({_pct22:.1f}%)",
                          ha="center", va="center", fontsize=7,
                          color=_rc22, linespacing=1.35,
                          bbox=dict(boxstyle="square,pad=0.32", fc="#FBFBFB", ec=_rc22,
                                    alpha=0.92, linewidth=0.9), zorder=5)

            _ctrl_items22 = [
                (_ctrl_3R3U_idx22,  "3R3U (control)",  "#8B0000"),
                (_ctrl_deha4_idx22, "DeHa4 (control)", "#00408B"),
            ]
            _ctrl_rpos22 = {
                (True,  True,  True):  (5.00, 5.10),
                (True,  True,  False): (5.00, 6.90),
                (True,  False, True):  (2.80, 4.40),
                (False, True,  True):  (7.20, 4.40),
                (True,  False, False): (1.85, 5.80),
                (False, True,  False): (8.15, 5.80),
                (False, False, True):  (5.00, 2.20),
                (False, False, False): (0.70, 0.70),
            }
            _ctrl_nudge22 = [(-0.45, +0.52), (+0.45, +0.52)]
            _ctrl_legend_h22 = []
            for _ci22, (_cidx22, _clbl22, _ccol22) in enumerate(_ctrl_items22):
                if _cidx22 is None:
                    continue
                _rk22 = (_cidx22 in _set_A, _cidx22 in _set_B, _cidx22 in _set_C)
                _bx22, _by22 = _ctrl_rpos22[_rk22]
                _nx22, _ny22 = _ctrl_nudge22[_ci22]
                _mx22, _my22 = _bx22 + _nx22, _by22 + _ny22
                ax22.scatter([_mx22], [_my22], marker="*", s=280,
                             color=_ccol22, edgecolors="black",
                             linewidths=0.6, zorder=16, clip_on=False)
                _ctrl_legend_h22.append(
                    _L22([0], [0], marker="*", color="w",
                         markerfacecolor=_ccol22, markeredgecolor="black",
                         markersize=11, markeredgewidth=0.6,
                         label=f"★  {_clbl22}"))

            _leg22 = [
                _Patch22(facecolor="#0072B2", alpha=0.45, edgecolor="#0072B2",
                         linewidth=1.2, label=f"Physics-Strong  (top-tercile {_phys_col22})"),
                _Patch22(facecolor="#D55E00", alpha=0.45, edgecolor="#D55E00",
                         linewidth=1.2, label="AI-Strong  (top-tercile confidence)"),
                _Patch22(facecolor="#009E73", alpha=0.45, edgecolor="#009E73",
                         linewidth=1.2, label=f"Tier-Strong  ({CFG.TIER_TOP}, {CFG.TIER_ORDER[1]}, {CFG.TIER_ORDER[2]})"),
            ]
            if _n_PA22 > 0:
                _leg22.append(_L22([0], [0], color="#CC79A7", linewidth=2.2,
                                   linestyle="--",
                                   label=f"{CFG.TIER_TOP} overlay  (n={_n_PA22:,})"))
            _leg22.append(_L22([0], [0], color="#444", linewidth=1.8,
                               linestyle=":", marker="o", markersize=12,
                               markerfacecolor="none", markeredgecolor="#444",
                               markeredgewidth=1.5,
                               label=f"All-three intersection  ({_ABC:,})"))
            _leg22.extend(_ctrl_legend_h22)
            # Legend pulled close under the diagram (footnote clutter removed)
            ax22.legend(handles=_leg22,
                        loc="upper center", bbox_to_anchor=(0.5, 0.02),
                        ncol=3, fontsize=8, framealpha=0.90, fancybox=True,
                        borderpad=0.5, columnspacing=1.2, handletextpad=0.5)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_23_Category_Overlap_Euler.png",
                        dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        else:
            reporter.log("  ! Figure 23 skipped: requires a catalytic-geometry column (soft_catalytic_score/mechanistic_score/SN2_Attack_Angle), a confidence column, and degrader_tier")
    except Exception as e:
        reporter.log(f"  ! Figure 23 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    _fig23_multitarget(df, out_dir, reporter)
    _fig24_sankey(df, out_dir, reporter)
    _fig25_pfas_size(df, out_dir, reporter)

    '''
    07_Diagnostic_and_MultiModel_Trends is the last on-disk folder, generated here
    after folder 06 so the folders are written in ascending order (01 → 07) and the
    run log reads folder-by-folder. Its filenames are not in the rename map, so the
    active savefig redirect passes them straight through to the diagnostic folder.
    '''
    _diag_dir = out_dir / "07_Diagnostic_and_MultiModel_Trends"
    _diag_dir.mkdir(parents=True, exist_ok=True)
    generate_additional_figures(df, _diag_dir, reporter)


    # savefig routing is installed/restored by the _redirect_savefig context
    # manager around this call, so no manual restore is needed here.


# =============================================================================
# SECTION 4C: Publication Assembly Figures (23–26: multitarget, top-tier breakdown,
#             Sankey workflow, PFAS-size composites)
# =============================================================================

def _fig23_multitarget(df: pd.DataFrame, out_dir: Path, reporter) -> None:
    """Figure 24: Top 25 multi-target proteins (stacked bar)."""
    try:
        import matplotlib.patches as _mp

        d = _utils_mod.standardise_dataframe_tiers(df.copy(), CFG)
        d["tier"] = d[CFG.COL_TIER].astype(str)

        tier_order = [t for t in CFG.TIER_ORDER if t in TIER_PALETTE]
        tier_colors = {t: TIER_PALETTE.get(t, "#999") for t in tier_order}
        d["tier_rank"] = d["tier"].map(CFG.TIER_RANK).fillna(0)
        d["prot"] = d[CFG.COL_PROT] if CFG.COL_PROT in d.columns else d.get("Protein_Name", d.get("protein", ""))
        d["lig"] = d[CFG.COL_LIG] if CFG.COL_LIG in d.columns else d.get("Ligand_Name", d.get("ligand", ""))
        d["Ensemble_Score"] = pd.to_numeric(
            d.get("Ensemble_Score", 0), errors="coerce"
        ).fillna(0)

        # Within-tier tiebreak on the authoritative Step-02 ranking (Scientific_Rank →
        # competence_score), NOT the figure-layer Ensemble_Score, so the representative
        # complex here matches the pipeline's chosen pose.
        if "Scientific_Rank" in d.columns:
            d["_best_key"] = pd.to_numeric(d["Scientific_Rank"], errors="coerce").fillna(1e9)
            _key_asc = True
        elif "competence_score" in d.columns:
            d["_best_key"] = pd.to_numeric(d["competence_score"], errors="coerce").fillna(0)
            _key_asc = False
        else:
            d["_best_key"] = pd.to_numeric(d.get("Ensemble_Score", 0), errors="coerce").fillna(0)
            _key_asc = False

        best_pairs = (
            d.sort_values(["prot", "lig", "tier_rank", "_best_key"],
                          ascending=[True, True, False, _key_asc])
             .groupby(["prot", "lig"], as_index=False)
             .first()
        )
        protein_summary = best_pairs.groupby("prot").agg(
            total_pfases=("lig", "nunique")
        ).reset_index()
        tier_counts_piv = (
            best_pairs.pivot_table(
                index="prot", columns="tier", values="lig",
                aggfunc="nunique", fill_value=0
            ).reindex(columns=tier_order, fill_value=0)
        )
        protein_summary = protein_summary.join(tier_counts_piv, on="prot").fillna(0)
        protein_summary = protein_summary.astype({t: int for t in tier_order})
        # Quality-weighted sort: sum of (tier_rank × count) then total as tiebreaker
        _tier_wt23 = {CFG.TIER_TOP: 6, CFG.TIER_ORDER[1]: 5, CFG.TIER_ORDER[2]: 4,
                      CFG.TIER_ORDER[3]: 3, CFG.TIER_ORDER[4]: 2, CFG.TIER_DECOY: 1}
        protein_summary["_qscore"] = sum(
            protein_summary[t] * _tier_wt23.get(t, 0) for t in tier_order
        )
        protein_summary = protein_summary.sort_values(
            ["_qscore", "total_pfases", CFG.TIER_TOP, CFG.TIER_ORDER[1]],
            ascending=False
        )
        top_proteins = protein_summary.head(25).set_index("prot")

        # ── Figure 24a: stacked bar — top-25 proteins ────────────────────────
        fig23a, ax = plt.subplots(figsize=(16, 9))
        """
        Decoys are NOT degraded — exclude the Tier-5 decoy segment so the bar
        length is the true count of PFAS actually degraded. This keeps the
        x-axis label honest and makes the breadth ranking visible, instead of
        every bar hitting the panel maximum (full panel = degraders + decoy).
        """
        bar_tiers = [t for t in tier_order if t != CFG.TIER_DECOY]
        bar_data = top_proteins[bar_tiers]
        bar_data.plot(kind="barh", stacked=True, ax=ax,
                      color=[tier_colors[t] for t in bar_tiers],
                      edgecolor="none", legend=False)
        ax.invert_yaxis()
        ax.set_xlabel(
            "Number of unique PFAS ligands degraded  "
            "(best tier per protein–ligand pair; Tier-5 decoys excluded)",
            fontsize=10)
        ax.set_ylabel("Protein  (top 25, ranked by quality-weighted degradation breadth)",
                      fontsize=10)
        ax.tick_params(axis="y", labelsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.xaxis.grid(True, color="#DDDDDD", linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        for i, prot in enumerate(bar_data.index):
            total  = int(bar_data.loc[prot].sum())
            pa_cnt = int(bar_data.loc[prot, CFG.TIER_TOP]) \
                     if CFG.TIER_TOP in bar_data.columns else 0
            _lbl = f"{total}" if pa_cnt == 0 else f"{total}  ★×{pa_cnt}"
            ax.text(total + 0.15, i, _lbl, va="center", fontsize=8.5,
                    color=tier_colors[CFG.TIER_TOP] if pa_cnt > 0 else "#444",
                    fontweight="bold" if pa_cnt > 0 else "normal")
        # Legend: only the tiers actually present in the plotted bars, single row.
        _present23  = [t for t in bar_tiers if int(bar_data[t].sum()) > 0]
        legend_h23  = [_mp.Patch(color=tier_colors[t], label=t.replace("_", " "))
                       for t in _present23]
        _any_top23  = int(bar_data[CFG.TIER_TOP].sum()) if CFG.TIER_TOP in bar_data.columns else 0
        _title23    = "Degrader tier" + (
            f"  (★×N = {CFG.TIER_TOP.replace('_', ' ')} count per protein)"
            if _any_top23 > 0 else "")
        fig23a.legend(handles=legend_h23, loc="lower center",
                      ncol=len(legend_h23), frameon=False, fontsize=9.5,
                      title=_title23, title_fontsize=9)
        # Figure title removed per request.
        plt.tight_layout(rect=[0, 0.07, 1, 1])
        plt.savefig(out_dir / "Figure_24a_Top25_Multitarget_Proteins.png",
                    dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close(fig23a)

        _fig23b_toptier_breakdown(best_pairs, tier_order, tier_colors, out_dir, reporter)

    except Exception as e:
        reporter.log(f"  ! Figure 24 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig


def _fig23b_toptier_breakdown(best_pairs, tier_order, tier_colors, out_dir, reporter) -> None:
    """Figure 24b: for every protein that reaches the TOP tier (CFG.TIER_TOP) on at
    least one PFAS, show the full set of PFAS it degrades as a stacked bar. Each
    segment = one PFAS ligand, coloured BY LIGAND; that ligand's tier is printed on
    the segment in short form (T1A, T1B, T2A, …). Reveals whether the elite
    proteins are also broad degraders across the panel."""
    try:
        import matplotlib as _mpl
        import matplotlib.patches as _mp

        top_tier = CFG.TIER_TOP
        bp = best_pairs.copy()
        bp["tier"] = bp["tier"].astype(str)

        top_prots = bp.loc[bp["tier"] == top_tier, "prot"].unique()
        if len(top_prots) == 0:
            reporter.log(f"  ! Figure 24b skipped: no proteins reached {top_tier}")
            return

        # Degrader pairs (drop decoy) for the top-tier proteins
        deg = bp[(bp["prot"].isin(top_prots)) & (bp["tier"] != CFG.TIER_DECOY)].copy()
        deg["tier_rank"] = deg["tier"].map(CFG.TIER_RANK).fillna(0)

        # Stable distinct ligand → colour map over the full PFAS panel
        lig_order = list(CFG.PFAS_ORDER) + sorted(
            l for l in deg["lig"].unique() if l not in CFG.PFAS_ORDER)
        _palette = list(_mpl.colormaps["tab20"].colors) + list(_mpl.colormaps["tab20b"].colors)
        lig_colour = {l: _palette[i % len(_palette)] for i, l in enumerate(lig_order)}

        def _tshort(t):
            return t.replace("Tier_", "T").replace("_Decoy", "")

        """
        Auto-contrast label colour: dark text on light blocks, white on dark
        (perceived luminance, ITU-R BT.601). Keeps every tier tag readable
        regardless of the ligand colour beneath it.
        """
        def _label_col(c):
            return _auto_label_colour(c)

        """
        Rank proteins: broadest PFAS coverage first (number degraded), then best
        tier quality (mean tier rank, higher = stronger) — best degrader on top,
        weakest at the bottom.
        """
        order = (deg.assign(_is_top=(deg["tier"] == top_tier).astype(int))
                    .groupby("prot")
                    .agg(_ntop=("_is_top", "sum"), _ntot=("lig", "nunique"),
                         _qual=("tier_rank", "mean"))
                    .sort_values(["_ntot", "_qual"], ascending=False))
        prots = list(order.index)[:25]

        fig, ax = plt.subplots(figsize=(15, max(4.0, 0.55 * len(prots) + 2)))
        # Nested heights: ligand block < tier-category frame.
        _BH, _CAT_H = 0.72, 0.82
        for yi, prot in enumerate(prots):
            sub = (deg[deg["prot"] == prot]
                   .sort_values("tier_rank", ascending=False)
                   .drop_duplicates("lig"))
            left = 0
            _tier_seq = []
            for _, r in sub.iterrows():
                col = lig_colour.get(r["lig"], "#999999")
                ax.barh(yi, 1.0, left=left, height=_BH, color=col, edgecolor="white",
                        linewidth=0.8, zorder=3)
                ax.text(left + 0.5, yi, _tshort(r["tier"]), ha="center", va="center",
                        fontsize=6.3, fontweight="bold", color=_label_col(col), zorder=5)
                _tier_seq.append(r["tier"])
                left += 1

            """
            Frames = each contiguous tier run, bordered in that tier's colour,
            so categories (e.g. 1×T1A, 6×T2A, …) read as distinct sub-groups.
            """
            _s = 0
            while _s < len(_tier_seq):
                _e = _s
                while _e + 1 < len(_tier_seq) and _tier_seq[_e + 1] == _tier_seq[_s]:
                    _e += 1
                _tc = tier_colors.get(_tier_seq[_s], "#666666")
                ax.add_patch(_mp.Rectangle((_s, yi - _CAT_H / 2), (_e - _s + 1), _CAT_H,
                             fill=False, edgecolor=_tc, linewidth=2.3, zorder=4))
                _s = _e + 1

            ntop = int(order.loc[prot, "_ntop"])
            ax.text(left + 0.25, yi, f"{int(left)}  ★×{ntop}", va="center",
                    fontsize=8, fontweight="bold",
                    color=tier_colors.get(top_tier, "#444444"))
        ax.set_xlim(-0.35, float(order["_ntot"].max()) + 3.0)

        ax.set_yticks(range(len(prots)))
        ax.set_yticklabels(prots, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Number of PFAS degraded  (each block = one PFAS; block label = its tier)",
                      fontsize=10)
        ax.set_ylabel(f"Top-tier proteins  (reach {_tshort(top_tier)} on ≥1 PFAS;  ★×N = {_tshort(top_tier)} count)",
                      fontsize=10)
        for _s in ("top", "right", "left"):
            ax.spines[_s].set_visible(False)
        ax.xaxis.grid(True, color="#DDDDDD", linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)

        shown = [l for l in lig_order
                 if l in set(deg[deg["prot"].isin(prots)]["lig"])]
        leg = [_mp.Patch(color=lig_colour[l], label=l) for l in shown]
        ax.legend(handles=leg, loc="center left", bbox_to_anchor=(1.01, 0.5),
                  fontsize=7, frameon=False, title="PFAS ligand", title_fontsize=8)
        plt.tight_layout()
        out = out_dir / "Figure_24b_TopTier_Protein_PFAS_Breakdown.png"
        plt.savefig(out, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        reporter.log(f"  ! Figure 24b skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig


def _fig24_sankey(df: pd.DataFrame, out_dir: Path, reporter) -> None:
    """
    Figure 25: Extended Sankey — multiple geometric constraints → final tier.
    Layout (up to 8 columns): All Complexes (thin) → Nuc Dist → Clamp ARG111 →
    Clamp ARG114 → Acid ASP134 → Base HIS277 → SN2 Angle → Final Tier.
    Columns absent from the data are silently dropped.
    """
    try:
        from matplotlib.path import Path as _MplPath
        from matplotlib.patches import PathPatch as _PathPatch
        import matplotlib.patches as _mp
        import textwrap as _tw24

        d = _utils_mod.standardise_dataframe_tiers(df.copy(), CFG)
        # Drop any stale computed columns so re-assignment never raises "already exists"
        for _stale24 in ["nuc_cat", "clamp_cat", "stab_cat", "triad_cat",
                         "ang_cat", "mech_cat", "tier_cat"]:
            if _stale24 in d.columns:
                d.drop(columns=[_stale24], inplace=True)

        """
        Binning — each column mirrors a REAL production tier gate, with cut
        points taken straight from CFG (not hand-picked). The five redundant
        per-residue ligand distances are replaced by the actual rule checks:
        nucleophile distance · carboxylate clamp · halide stabilisation ·
        catalytic-triad geometry · SN2 angle · mechanistic score.
        """
        def _num24(col, default=np.inf):
            return pd.to_numeric(d[col], errors="coerce") if col in d.columns \
                else pd.Series(default, index=d.index)

        # 1) Nucleophile distance — cuts = CFG.TIER_NUC_DIST gate values
        _nd = CFG.TIER_NUC_DIST
        _nuc_bins_24 = [-np.inf, _nd[CFG.TIER_TOP], _nd[CFG.TIER_ORDER[2]], _nd[CFG.TIER_ORDER[3]], _nd[CFG.TIER_ORDER[4]], np.inf]
        _nuc_lbls_24 = [f"≤{_nd[CFG.TIER_TOP]:.1f}Å",
                        f"{_nd[CFG.TIER_TOP]:.1f}–{_nd[CFG.TIER_ORDER[2]]:.1f}Å",
                        f"{_nd[CFG.TIER_ORDER[2]]:.1f}–{_nd[CFG.TIER_ORDER[3]]:.1f}Å",
                        f"{_nd[CFG.TIER_ORDER[3]]:.1f}–{_nd[CFG.TIER_ORDER[4]]:.1f}Å",
                        f">{_nd[CFG.TIER_ORDER[4]]:.1f}Å"]
        d["nuc_cat"] = (pd.cut(_num24("Dist_Nucleophile"),
                               bins=_nuc_bins_24, labels=_nuc_lbls_24)
                        .astype(str).fillna("Unknown")) if "Dist_Nucleophile" in d.columns else "Unknown"

        # 2) Carboxylate clamp (ARG111 OR ARG114 ≤5 Å) — production clamp_ok
        d["clamp_cat"] = np.where(_num24("carboxylate_clamp_integrity", 0).fillna(0) >= 0.5,
                                  "Clamp intact", "Clamp broken")
        # 3) Halide stabilisation (Trp/Tyr/polar sidechain ≤5.5 Å)
        d["stab_cat"] = np.where(_num24("halide_stabilisation_score", 0).fillna(0) >= 0.5,
                                 "Stabilised", "Unstabilised")
        # 4) Catalytic-triad geometry (internal Nuc–Base & Base–Acid) — CFG cuts
        _nb24, _ba24 = _num24("dist_nuc_base_internal"), _num24("dist_base_acid_internal")
        _nbm, _bam = CFG.TIER_NB_MAX, CFG.TIER_BA_MAX
        _triad_lbls_24 = ["Triad tight", "Triad moderate", "Triad loose"]
        d["triad_cat"] = np.where((_nb24 <= _nbm[CFG.TIER_TOP]) & (_ba24 <= _bam[CFG.TIER_TOP]), "Triad tight",
                          np.where((_nb24 <= _nbm[CFG.TIER_ORDER[2]]) & (_ba24 <= _bam[CFG.TIER_ORDER[2]]), "Triad moderate",
                                   "Triad loose"))
        # 5) SN2 attack angle — cuts = CFG.TIER_ANGLE_MIN
        _am = CFG.TIER_ANGLE_MIN
        _ang_ord_24 = [f"<{_am[CFG.TIER_ORDER[3]]:.0f}°",
                       f"{_am[CFG.TIER_ORDER[3]]:.0f}–{_am[CFG.TIER_ORDER[2]]:.0f}°",
                       f"{_am[CFG.TIER_ORDER[2]]:.0f}–{_am[CFG.TIER_ORDER[1]]:.0f}°",
                       f"{_am[CFG.TIER_ORDER[1]]:.0f}–{_am[CFG.TIER_TOP]:.0f}°",
                       f"≥{_am[CFG.TIER_TOP]:.0f}°"]
        d["ang_cat"] = (pd.cut(_num24("SN2_Attack_Angle", -1),
                               bins=[-np.inf, _am[CFG.TIER_ORDER[3]], _am[CFG.TIER_ORDER[2]], _am[CFG.TIER_ORDER[1]], _am[CFG.TIER_TOP], np.inf],
                               labels=_ang_ord_24)
                        .astype(str).fillna("Unknown")) if "SN2_Attack_Angle" in d.columns else "Unknown"
        # 6) Mechanistic score — cuts = CFG.TIER_MECH_MIN
        _mm = CFG.TIER_MECH_MIN
        _mech_lbls_24 = [f"<{_mm[CFG.TIER_ORDER[2]]:.1f}",
                         f"{_mm[CFG.TIER_ORDER[2]]:.1f}–{_mm[CFG.TIER_ORDER[1]]:.1f}",
                         f"{_mm[CFG.TIER_ORDER[1]]:.1f}–{_mm[CFG.TIER_TOP]:.1f}",
                         f"≥{_mm[CFG.TIER_TOP]:.1f}"]
        _mech_edges_24 = [-np.inf, _mm[CFG.TIER_ORDER[2]], _mm[CFG.TIER_ORDER[1]], _mm[CFG.TIER_TOP], np.inf]
        # Collapse bins whose CFG thresholds coincide (e.g. Tier_1A == Tier_1B mech floor),
        # keeping edges strictly increasing for pd.cut and dropping each vanished bin's label.
        _me24, _ml24 = [_mech_edges_24[0]], []
        for _i24, _lbl24 in enumerate(_mech_lbls_24):
            if _mech_edges_24[_i24 + 1] > _me24[-1]:
                _me24.append(_mech_edges_24[_i24 + 1])
                _ml24.append(_lbl24)
        _mech_edges_24, _mech_lbls_24 = _me24, _ml24
        d["mech_cat"] = (pd.cut(_num24("mechanistic_score", np.nan),
                                bins=_mech_edges_24,
                                labels=_mech_lbls_24)
                         .astype(str).fillna("Unknown")) if "mechanistic_score" in d.columns else "Unknown"

        # Order mirrors 02_Production tier cascade: Tier_1A → … → Tier_3 → Tier_4 → Tier_5_Decoy
        _tier_ord_24 = CFG.TIER_ORDER + ["Other"]
        d["tier_cat"] = d[CFG.COL_TIER].astype(str)
        d.loc[~d["tier_cat"].isin(_tier_ord_24[:-1]), "tier_cat"] = "Other"

        # ── Column definitions (best category first → drawn on top) ────────────
        _all_intermed_24 = [
            ("nuc_cat",   _nuc_lbls_24,                "Nuc Distance",         "(Å · ASP110)",          "#1B4D2E"),
            ("clamp_cat", ["Clamp intact", "Clamp broken"], "Carboxylate Clamp", "(Arg/Lys clamp)",   "#2C4A1E"),
            ("stab_cat",  ["Stabilised", "Unstabilised"],   "Halide Stabilisation", "(Trp / Tyr / polar)", "#1E4A3A"),
            ("triad_cat", _triad_lbls_24,              "Triad Geometry",       "(Nuc–Base · Base–Acid)", "#1E3A4A"),
            ("ang_cat",   list(reversed(_ang_ord_24)), "SN2 Angle",            "(degrees)",             "#1B3A5E"),
            ("mech_cat",  list(reversed(_mech_lbls_24)), "Mech Score",         "(score 0–1)",           "#3A1E4A"),
        ]
        _valid_im_24 = [(c, o, t, s, b) for c, o, t, s, b in _all_intermed_24
                        if (d[c] != "Unknown").any()]

        _col_seq_24  = ["all"] + [c for c, *_ in _valid_im_24] + ["tier_cat"]
        _ord_seq_24  = [None]  + [o for _, o, *_ in _valid_im_24] + [_tier_ord_24]
        _ttl_seq_24  = ["All Complexes"] + [t for _, _, t, *_ in _valid_im_24] + ["Final Tier"]
        _sub_seq_24  = ["(starting pool)"] + [s for _, _, _, s, *_ in _valid_im_24] + ["(degradation tier)"]
        _bg_seq_24   = ["#2E4053"] + [b for _, _, _, _, b in _valid_im_24] + ["#2E1B5E"]
        _n_cols_24   = len(_col_seq_24)

        # ── Colours — sourced from CFG § 8.10 (best=green … worst=red) ─────────
        _colmap_24 = {
            "nuc_cat":   dict(zip(_nuc_lbls_24, CFG.SANKEY_GRAD5)),
            "clamp_cat": dict(CFG.SANKEY_CLAMP_COLOUR),
            "stab_cat":  dict(CFG.SANKEY_STAB_COLOUR),
            "triad_cat": dict(CFG.SANKEY_TRIAD_COLOUR),
            # angle labels ascend (<… → ≥…); colour worst→best so the top ≥ band is green
            "ang_cat":   dict(zip(_ang_ord_24, list(reversed(CFG.SANKEY_GRAD5)))),
            "mech_cat":  dict(zip(_mech_lbls_24, CFG.SANKEY_MECH_GRAD)),
        }
        _tier_clr_24  = {t: TIER_PALETTE.get(t, "#999")
                         for t in [CFG.TIER_TOP,CFG.TIER_ORDER[1],CFG.TIER_ORDER[2],CFG.TIER_ORDER[3],CFG.TIER_ORDER[4],CFG.TIER_POOR,CFG.TIER_DECOY,"Other"]}
        _tier_alp_24  = {
            CFG.TIER_TOP: 0.72, CFG.TIER_ORDER[1]: 0.68, CFG.TIER_ORDER[2]: 0.62,
            CFG.TIER_ORDER[3]: 0.57, CFG.TIER_ORDER[4]: 0.48, CFG.TIER_POOR: 0.42, CFG.TIER_DECOY: 0.36, "Other": 0.28,
        }

        def _clr_24(cat_col, lbl):
            if cat_col == "tier_cat":
                return TIER_PALETTE.get(lbl, _tier_clr_24.get(lbl, "#CCC"))
            return _colmap_24.get(cat_col, {}).get(lbl, "#CCC")

        # ── Layout ────────────────────────────────────────────────────────────
        _BY0_24  = 0.07
        _BH_24   = 0.78
        _BGAP_24 = 0.012
        _bw0_24  = 0.025    # thin "All Complexes" bar
        _bw_24   = 0.093    # normal column width

        _n_norm_24 = _n_cols_24 - 1
        _gap_24    = (0.970 - _bw0_24 - _n_norm_24 * _bw_24) / max(_n_norm_24, 1)
        _xs_24     = [0.010]
        for _ci24 in range(1, _n_cols_24):
            _xs_24.append(_xs_24[-1] + (_bw0_24 if _ci24 == 1 else _bw_24) + _gap_24)

        _total_j_24 = max(len(d), 1)
        _rank_m_24  = {t: -i for i, t in enumerate(_tier_ord_24)}

        def _pos24(order, cnts):
            nz = {k: int(cnts.get(k, 0)) for k in order if int(cnts.get(k, 0)) > 0}
            if not nz: return {}
            total_h = _BH_24 - _BGAP_24 * max(len(nz) - 1, 0)
            total_c = max(sum(nz.values()), 1)
            """
            Stack TOP-DOWN so the first label in `order` sits at the top of the
            column (best category on top — Tier_1A, ≤3.0Å, ≥170°, …).
            """
            pos = {}
            y = _BY0_24 + _BH_24
            for lbl, cnt in nz.items():
                h = max(0.016, total_h * cnt / total_c)
                pos[lbl] = (y - h, y)
                y -= (h + _BGAP_24)
            return pos

        # Positions and counts for all non-'all' columns
        _cpos_24  = {}   # cat_col → {lbl: (y0, y1)}
        _ccnt_24  = {}   # cat_col → {lbl: int}
        for _cc, _co in zip(_col_seq_24[1:], _ord_seq_24[1:]):
            _vc = d[_cc].value_counts().to_dict()
            _ccnt_24[_cc] = _vc
            _cpos_24[_cc] = _pos24(_co, _vc)

        # ── Bezier ribbon ─────────────────────────────────────────────────────
        def _brib_24(x0, bw_s, x1, y0a, y1a, y0b, y1b, color, alpha, zorder):
            _ha = max(y1a - y0a, 0.0)
            _hb = max(y1b - y0b, 0.0)
            if _ha < 1e-5 and _hb < 1e-5: return
            _mv = 0.0015                       # smaller min height → less overflow past boxes
            _ha = max(_ha, _mv) if _ha > 0 else _ha
            _hb = max(_hb, _mv) if _hb > 0 else _hb
            """
            Tuck both ends a hair INTO the node boxes so the junction is flush —
            the boxes (higher zorder) cover the overlap, leaving no gap/overflow.
            """
            _pad24 = 0.002
            _xa = x0 + bw_s - _pad24
            _xe = x1 + _pad24
            cx = (_xa + _xe) / 2
            verts = [
                (_xa, y0a), (cx, y0a), (cx, y0b), (_xe, y0b),
                (_xe, y0b + _hb), (cx, y0b + _hb), (cx, y0a + _ha), (_xa, y0a + _ha),
                (_xa, y0a),
            ]
            codes = [_MplPath.MOVETO,
                     _MplPath.CURVE4, _MplPath.CURVE4, _MplPath.CURVE4,
                     _MplPath.LINETO,
                     _MplPath.CURVE4, _MplPath.CURVE4, _MplPath.CURVE4,
                     _MplPath.CLOSEPOLY]
            ax24.add_patch(_PathPatch(_MplPath(verts, codes),
                                      facecolor=color, alpha=alpha,
                                      edgecolor="none", zorder=zorder))
            _pct = max(_ha, _hb) / _BH_24 * 100
            if _ha > 0.033 and _hb > 0.033 and _pct >= 5.0:
                ax24.text((x0 + bw_s + x1) / 2,
                          ((y0a + _ha / 2) + (y0b + _hb / 2)) / 2,
                          f"{_pct:.0f}%",
                          ha="center", va="center", fontsize=14, fontweight="bold",
                          color="white", zorder=zorder + 1,
                          bbox=dict(boxstyle="round,pad=0.05", fc=color, ec="none", alpha=0.70))

        # ── Node box ──────────────────────────────────────────────────────────
        def _dbox_24(x, bw, y0, y1, label, count, box_color):
            _h = max(y1 - y0, 0.001)
            # Slightly blunt corners (pad=0 → no halo; small rounding softens edges)
            ax24.add_patch(_mp.FancyBboxPatch(
                (x, y0), bw, _h,
                boxstyle="round,pad=0,rounding_size=0.004", mutation_aspect=1.0,
                facecolor=box_color, edgecolor="white",
                linewidth=1.0, alpha=0.96, zorder=5))
            if _h < 0.006:
                return
            """
            Auto-contrast: white text on dark boxes, dark text on light boxes,
            chosen from the box's own luminance (fixes dark text on dark fills).
            """
            lbl_col = _auto_label_colour(box_color)
            """
            Show the actual COUNT (not %) — a rare class like Tier_1A (n=8) must
            read "8", never a rounded "0.0%".
            """
            _cnt_s = f"{int(count):,}"
            """
            Single line "label  count" (wrap only if long). Font scales with box
            height so the value FILLS the otherwise-empty top/bottom space; thin
            boxes keep a readable minimum and simply overlay the box ("on top").
            """
            _txt = _tw24.fill(f"{label}  {_cnt_s}", width=15)
            _fs  = float(np.clip(_h * 520, 9.0, 16.0))
            ax24.text(x + bw / 2, y0 + _h / 2, _txt,
                      ha="center", va="center", fontsize=_fs, fontweight="bold",
                      color=lbl_col, zorder=10, linespacing=1.0)

        # ── Figure ────────────────────────────────────────────────────────────
        fig24, ax24 = plt.subplots(figsize=(30, 16))
        ax24.set_xlim(0, 1); ax24.set_ylim(0, 1.05); ax24.axis("off")
        fig24.patch.set_facecolor("white")

        # ── Column 0: thin "All Complexes" bar (90°-rotated label) ────────────
        ax24.add_patch(_mp.FancyBboxPatch(
            (_xs_24[0], _BY0_24), _bw0_24, _BH_24,
            boxstyle="round,pad=0,rounding_size=0.004", mutation_aspect=1.0,
            facecolor="#2E4053", edgecolor="white", linewidth=1.0, alpha=0.96, zorder=5))
        ax24.text(_xs_24[0] + _bw0_24 / 2, _BY0_24 + _BH_24 / 2,
                  f"All Complexes  {_total_j_24:,}",
                  ha="center", va="center", fontsize=15, fontweight="bold",
                  color="white", zorder=6, rotation=90)

        # ── Node boxes for columns 1..n_cols-1 ───────────────────────────────
        for _ci24, _cc24 in enumerate(_col_seq_24[1:], start=1):
            for lbl, (y0, y1) in _cpos_24.get(_cc24, {}).items():
                _dbox_24(_xs_24[_ci24], _bw_24, y0, y1, lbl,
                         int(_ccnt_24.get(_cc24, {}).get(lbl, 0)),
                         _clr_24(_cc24, lbl))

        # ── Column header banners — extra gap so they clear the node content ──
        _hdr_gap_24 = 0.078   # wider gap → more separation between content and headers
        _hdr_y_24   = _BY0_24 + _BH_24 + _hdr_gap_24
        _hbox_24    = []      # (left, right) x-extent of each header banner
        for _ci24, (_htxt, _hsub, _hbg) in enumerate(zip(_ttl_seq_24, _sub_seq_24, _bg_seq_24)):
            if _ci24 == 0:
                """
                thin source bar → give its header a wider banner with horizontal,
                wrapped text so the full title shows (no rotation/clipping)
                """
                _cbw   = 0.066
                _hx    = max(_xs_24[0] + _bw0_24 / 2 - _cbw / 2, 0.0)
                _title = _tw24.fill(_htxt, width=8)
            else:
                _cbw   = _bw_24
                _hx    = _xs_24[_ci24]
                _title = _tw24.fill(_htxt, width=14)
            _hbox_24.append((_hx, _hx + _cbw))
            ax24.add_patch(_mp.FancyBboxPatch(
                (_hx, _hdr_y_24), _cbw, 0.060,
                boxstyle="round,pad=0,rounding_size=0.005", mutation_aspect=1.0,
                facecolor=_hbg, edgecolor="none", alpha=0.95, zorder=7))
            ax24.text(_hx + _cbw / 2, _hdr_y_24 + 0.038, _title,
                      ha="center", va="center", fontsize=14, fontweight="bold",
                      color="white", zorder=8, linespacing=1.0)
            ax24.text(_hx + _cbw / 2, _hdr_y_24 + 0.012, _hsub,
                      ha="center", va="center", fontsize=12, color="#CFE0EA", zorder=8)

        # ── Flow arrows between headers — anchored exactly at banner edges ────
        _ya24 = _hdr_y_24 + 0.030
        for _ci24 in range(_n_cols_24 - 1):
            """
            Tail flush on this banner's right edge, head flush on next banner's
            left edge — arrow touches both boxes (no gap, matches box shape).
            """
            _x_start = _hbox_24[_ci24][1]
            _x_end   = _hbox_24[_ci24 + 1][0]
            if _x_end <= _x_start:
                continue
            ax24.annotate("", xy=(_x_end, _ya24), xytext=(_x_start, _ya24),
                          arrowprops=dict(arrowstyle="->", color="#888", lw=1.4,
                                          mutation_scale=12, shrinkA=0, shrinkB=0),
                          zorder=9)

        # ── Ribbons col 0 → col 1 ─────────────────────────────────────────────
        _fc24    = _col_seq_24[1]
        _fo24    = _ord_seq_24[1]
        _flow01  = d.groupby([_fc24, "tier_cat"]).size().reset_index(name="count")
        _flow01["_nr"] = _flow01[_fc24].apply(lambda v: _fo24.index(v) if v in _fo24 else 999)
        _flow01["_tr"] = _flow01["tier_cat"].map(_rank_m_24).fillna(-99)
        _flow01  = _flow01.sort_values(["_nr", "_tr"])
        _used_all_24 = _BY0_24
        _dest01_24   = {k: v[0] for k, v in _cpos_24.get(_fc24, {}).items()}
        for _, r01 in _flow01.iterrows():
            if r01["count"] <= 0: continue
            nc, tc = r01[_fc24], r01["tier_cat"]
            if nc not in _cpos_24.get(_fc24, {}): continue
            _h_d = _cpos_24[_fc24][nc][1] - _cpos_24[_fc24][nc][0]
            _rhs = _BH_24 * r01["count"] / _total_j_24
            _rht = _h_d  * r01["count"] / max(int(_ccnt_24.get(_fc24, {}).get(nc, 1)), 1)
            y0a, y0b = _used_all_24, _dest01_24.get(nc, 0)
            _used_all_24 = y0a + _rhs
            _dest01_24[nc] = y0b + _rht
            _brib_24(_xs_24[0], _bw0_24, _xs_24[1],
                     y0a, y0a + _rhs, y0b, y0b + _rht,
                     TIER_PALETTE.get(tc, _tier_clr_24.get(tc, "#999")),
                     _tier_alp_24.get(tc, 0.40), zorder=2)

        # ── Ribbons col i → col i+1 (i = 1 .. n_cols-2) ─────────────────────
        for _ri in range(1, _n_cols_24 - 1):
            _sc24 = _col_seq_24[_ri]
            _tc24 = _col_seq_24[_ri + 1]
            _sp24 = _cpos_24.get(_sc24, {})
            _tp24 = _cpos_24.get(_tc24, {})
            _sn24 = _ccnt_24.get(_sc24, {})
            _tn24 = _ccnt_24.get(_tc24, {})
            _grp_cols_24 = [_sc24, _tc24] if _tc24 == "tier_cat" else [_sc24, _tc24, "tier_cat"]
            _frm  = (d.groupby(_grp_cols_24)
                     .size().reset_index(name="count"))
            if "tier_cat" not in _frm.columns:
                _frm["tier_cat"] = _frm[_tc24]
            _frm["_tr"] = _frm["tier_cat"].map(_rank_m_24).fillna(-99)
            _frm = _frm.sort_values("_tr")
            _so24 = {k: v[0] for k, v in _sp24.items()}
            _to24 = {k: v[0] for k, v in _tp24.items()}
            for _, rr in _frm.iterrows():
                if rr["count"] <= 0: continue
                s0, t0, tc = rr[_sc24], rr[_tc24], rr["tier_cat"]
                if s0 not in _sp24 or t0 not in _tp24: continue
                _hs = _sp24[s0][1] - _sp24[s0][0]
                _ht = _tp24[t0][1] - _tp24[t0][0]
                _rhs = _hs * rr["count"] / max(int(_sn24.get(s0, 1)), 1)
                _rht = _ht * rr["count"] / max(int(_tn24.get(t0, 1)), 1)
                y0a, y0b = _so24.get(s0, 0), _to24.get(t0, 0)
                _so24[s0] = y0a + _rhs
                _to24[t0] = y0b + _rht
                _brib_24(_xs_24[_ri], _bw_24, _xs_24[_ri + 1],
                         y0a, y0a + _rhs, y0b, y0b + _rht,
                         TIER_PALETTE.get(tc, _tier_clr_24.get(tc, "#999")),
                         _tier_alp_24.get(tc, 0.40), zorder=2 + _ri)

        """
        Legend removed — Final Tier column already carries the tier names/colours.
        Figure title removed per request.
        """
        plt.tight_layout(pad=0.3)
        plt.savefig(out_dir / "Figure_25_Sankey_Workflow.png",
                    dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight", facecolor="white")
        plt.close(fig24)

    except Exception as e:
        reporter.log(f"  ! Figure 25 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig


def _fig25_pfas_size(df: pd.DataFrame, out_dir: Path, reporter) -> None:
    """Figure 26a–25e — PFAS chain-length selectivity saved as five separate PNGs."""
    try:
        _req25 = ["total_fluorine_count", "SN2_Attack_Angle", "Boltz_Model_Confidence", "degrader_tier"]
        _miss25 = [c for c in _req25 if c not in df.columns]
        if _miss25:
            reporter.log(f"  ! Figure 26 skipped: missing columns {_miss25}")
            return

        _d25 = df[_req25].copy()
        for _c25 in _req25[:-1]:
            _d25[_c25] = pd.to_numeric(_d25[_c25], errors="coerce")
        _d25 = _d25.dropna(subset=_req25)
        _d25 = _d25[_d25["degrader_tier"].isin(TIER_ORDER_LOGIC)]
        if len(_d25) < 30:
            reporter.log("  ! Figure 26 skipped: insufficient data after filtering")
            return

        # PFAS chain-length bins (F-count as proxy: C4≈9F, C6≈13F, C8≈17F, C10≈21F) — from CFG (SSOT)
        _f25_bins   = list(CFG.VIS_PFAS_FCOUNT_BINS)
        # Labels paired with CFG.VIS_PFAS_FCOUNT_BINS (co-located in CFG so edits can't drift apart).
        _f25_labels       = list(CFG.VIS_PFAS_SIZE_LABELS)
        _f25_labels_clean = list(CFG.VIS_PFAS_SIZE_LABELS_SHORT)
        _d25["_size_bin"] = pd.cut(
            _d25["total_fluorine_count"], bins=_f25_bins,
            labels=_f25_labels, right=True
        )
        _d25 = _d25.dropna(subset=["_size_bin"])

        # Mechanistic outcome classification
        _ang25  = _d25["SN2_Attack_Angle"]
        _conf25 = _d25["Boltz_Model_Confidence"]
        _SA25, _IA25, _SC25 = CFG.SUBSTRATE_ANGLE_MIN, CFG.INHIBITOR_ANGLE_MAX, CFG.SUBSTRATE_CONF_MIN
        _d25["_outcome"] = np.where(
            (_ang25 >= _SA25) & (_conf25 >= _SC25), "Substrate",
            np.where(
                (_ang25 < _IA25) & (_conf25 >= _SC25), "Potential Inhibitor",
                np.where(
                    (_ang25 >= _SA25) & (_conf25 < _SC25), "Reactive (low conf)",
                    np.where((_ang25 >= _IA25) & (_ang25 < _SA25), "Borderline", "Non-reactive")
                )
            )
        )
        _oc25_order  = ["Substrate", "Borderline", "Reactive (low conf)",
                        "Non-reactive", "Potential Inhibitor"]
        _oc25_colors = dict(CFG.OUTCOME_COLOUR)   # sourced from CFG § 8.9

        # ── Panel A (25a): hexbin landscape + tier scatter + rolling median ──
        fig25a, axA = plt.subplots(figsize=(14, 8))
        axA.axhspan(CFG.SUBSTRATE_ANGLE_MIN, 185, color="#009E73", alpha=0.08, zorder=0)
        axA.axhspan(90, CFG.INHIBITOR_ANGLE_MAX, color="#D55E00", alpha=0.08, zorder=0)
        axA.axhline(CFG.SUBSTRATE_ANGLE_MIN, color="#009E73", lw=1.2, ls="--", alpha=0.7, zorder=2)
        axA.axhline(CFG.INHIBITOR_ANGLE_MAX, color="#D55E00", lw=1.2, ls="--", alpha=0.7, zorder=2)
        _hb25 = axA.hexbin(
            _d25["total_fluorine_count"], _d25["SN2_Attack_Angle"],
            gridsize=50, mincnt=1, cmap="YlOrRd", alpha=0.65, zorder=1, linewidths=0.2
        )
        _cb25 = fig25a.colorbar(_hb25, ax=axA, pad=0.01, aspect=30, shrink=0.85)
        _cb25.set_label("Count per hex", fontsize=9)
        _cb25.ax.tick_params(labelsize=8)
        _samp_A25 = _d25.sample(min(3000, len(_d25)), random_state=42)
        for _t25a in sorted(TIER_ORDER_LOGIC, key=lambda t: t == CFG.TIER_TOP):
            _is_pa25a = (_t25a == CFG.TIER_TOP)
            # Tier_1A is crucial and rare → plot ALL of its points (never sampled)
            _src25a = _d25 if _is_pa25a else _samp_A25
            _ts25 = _src25a[_src25a["degrader_tier"] == _t25a]
            if _ts25.empty:
                continue
            axA.scatter(
                _ts25["total_fluorine_count"], _ts25["SN2_Attack_Angle"],
                color=TIER_PALETTE.get(_t25a, "#999"),
                alpha=0.85 if _is_pa25a else 0.45,
                s=60 if _is_pa25a else 14,
                edgecolors="black" if _is_pa25a else "none",
                linewidths=0.8 if _is_pa25a else 0,
                zorder=5 if _is_pa25a else 3, label=_t25a
            )
        _roll25 = _d25[["total_fluorine_count", "SN2_Attack_Angle"]].sort_values("total_fluorine_count")
        if len(_roll25) >= 20:
            _win25 = max(20, len(_roll25) // 40)
            _roll25["_med"] = _roll25["SN2_Attack_Angle"].rolling(
                window=_win25, center=True, min_periods=5).median()
            axA.plot(_roll25["total_fluorine_count"], _roll25["_med"],
                     color="black", lw=2.2, zorder=6, label="Rolling median SN2")
        axA.text(0.99, 0.92, f"SUBSTRATE ZONE  (SN2 ≥ {CFG.SUBSTRATE_ANGLE_MIN:.0f}°)",
                 transform=axA.transAxes, ha="right", va="top", fontsize=9.5,
                 color="#007A52", fontweight="bold",
                 bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"))
        axA.text(0.99, 0.10, f"POTENTIAL INHIBITOR ZONE  (SN2 < {CFG.INHIBITOR_ANGLE_MAX:.0f}°)",
                 transform=axA.transAxes, ha="right", va="bottom", fontsize=9.5,
                 color="#C04000", fontweight="bold",
                 bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"))
        _bins25 = list(CFG.VIS_PFAS_FCOUNT_BINS)          # [0, 8, 13, 18, 24, inf]
        for _fb25 in _bins25[1:-1]:
            axA.axvline(_fb25, color="#444", lw=0.8, ls=":", alpha=0.5, zorder=2)
        _bin_toplbls25 = list(CFG.VIS_PFAS_SIZE_LABELS_SHORT)
        axA.set_ylim(85, 186)
        _fc_max25 = int(_d25["total_fluorine_count"].max()) + 2
        # Label midpoints derived from the CFG bin edges; the open-ended final bin uses the axis max.
        _edges25 = _bins25[:-1] + [max(_fc_max25, _bins25[-2] + 4)]
        _bin_xmids25 = [0.5 * (_edges25[i] + _edges25[i + 1]) for i in range(len(_edges25) - 1)]
        axA.set_xticks(range(0, _fc_max25 + 1, 2))
        for _bx25, _bl25 in zip(_bin_xmids25, _bin_toplbls25):
            axA.text(_bx25, 184, _bl25, ha="center", va="top",
                     fontsize=7.5, color="#444", style="italic")
        axA.set_xlabel("Total fluorine count  (proxy for carbon chain length)", fontsize=11)
        axA.set_ylabel("SN2 Attack Angle (°)", fontsize=11)
        axA.legend(loc="upper left", bbox_to_anchor=(0.0, 1.05),
                   ncol=len(TIER_ORDER_LOGIC) + 1, fontsize=7,
                   framealpha=0.92, fancybox=True, handlelength=0.8, borderpad=0.3)
        plt.tight_layout()
        _out25a = out_dir / "Figure_26a_PFAS_Size_Hexbin_Landscape.png"
        fig25a.savefig(_out25a, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close(fig25a)

        """
        Figure 26b: ONE figure, two stories per chain-length bin — a left
        stacked bar (mechanistic outcome) and a right stacked bar (degrader
        tier) placed side by side within each bin group.
        """
        _oc_ct25 = (
            _d25.groupby(["_size_bin", "_outcome"], observed=True)
            .size().unstack(fill_value=0)
        )
        _oc_ct25 = _oc_ct25.reindex(columns=[c for c in _oc25_order if c in _oc_ct25.columns])
        _oc_pct25 = _oc_ct25.div(_oc_ct25.sum(axis=1), axis=0) * 100

        _tier_ct25 = (
            _d25.groupby(["_size_bin", "degrader_tier"], observed=True)
            .size().unstack(fill_value=0)
        )
        _tier_ct25 = _tier_ct25.reindex(columns=[t for t in TIER_ORDER_LOGIC if t in _tier_ct25.columns])
        _tier_pct25 = _tier_ct25.div(_tier_ct25.sum(axis=1), axis=0) * 100

        _xlabels_25 = [_f25_labels_clean[_f25_labels.index(lb)] if lb in _f25_labels else str(lb)
                       for lb in _oc_pct25.index.tolist()]

        # Outcome palette + per-bin colours — sourced from CFG (§ 8.9, § 8.9b).
        _oc25_colors_b = dict(CFG.OUTCOME_COLOUR)
        _bin_ttlette_25 = list(CFG.PFAS_SIZE_BIN_COLOUR)
        _bin_cols_25 = [_bin_ttlette_25[i % len(_bin_ttlette_25)] for i in range(len(_oc_pct25))]

        fig25b, axB = plt.subplots(figsize=(15, 8))
        _xB, _wB, _offB = np.arange(len(_oc_pct25)), 0.27, 0.19   # bars 75% of prior width

        """
        Master container bar per x-category (like Fig 19a) — faint bin-coloured
        fill + matching coloured outline, holding both child stacked bars.
        """
        from matplotlib.colors import to_rgba as _to_rgba25
        for _xi, _bc in zip(_xB, _bin_cols_25):
            axB.bar(_xi, 104, 0.74, bottom=0, color=_to_rgba25(_bc, 0.10),
                    edgecolor=_bc, linewidth=1.4, zorder=2)

        def _stack_at_25(xpos, pct_df, ct_df, colour_fn, side):
            """
            Inside label when the segment is tall enough (≥4%); otherwise a small
            leader-arrow points to the thin segment with the value placed in the
            gap beside the bar (left bar → labels left, right bar → labels right).
            """
            _bot = np.zeros(len(pct_df))
            _dx  = -0.30 if side == "left" else 0.30
            _ha  = "right" if side == "left" else "left"
            _esign = -1 if side == "left" else 1
            for _col in pct_df.columns:
                _vals = pct_df[_col].values
                _ccol = colour_fn(_col)
                axB.bar(xpos, _vals, _wB, bottom=_bot, color=_ccol,
                        edgecolor="white", linewidth=0.6, zorder=3)
                for _xi, (_v, _b) in enumerate(zip(_vals, _bot)):
                    if _v <= 0:
                        continue
                    _ymid = _b + _v / 2
                    if _v >= 4:
                        axB.text(xpos[_xi], _ymid, f"{_v:.0f}%", ha="center",
                                 va="center", fontsize=8.5, color="white",
                                 fontweight="bold", zorder=4)
                    elif _v >= 0.8:
                        axB.annotate(f"{_v:.0f}%",
                                     xy=(xpos[_xi] + _esign * _wB / 2, _ymid),
                                     xytext=(xpos[_xi] + _dx, _ymid),
                                     ha=_ha, va="center", fontsize=8.5,
                                     color=_ccol, fontweight="bold", zorder=7,
                                     arrowprops=dict(arrowstyle="-", color=_ccol,
                                                     lw=0.5, shrinkA=1, shrinkB=1))
                _bot = _bot + _vals
            for _xi, _n in enumerate(ct_df.sum(axis=1).values):
                axB.text(xpos[_xi], 101.5, f"n={_n:,}", ha="center", va="bottom",
                         fontsize=6.3, color=_bin_cols_25[_xi], fontweight="bold")

        _stack_at_25(_xB - _offB, _oc_pct25, _oc_ct25, lambda c: _oc25_colors_b.get(c, "#999"), "left")
        _stack_at_25(_xB + _offB, _tier_pct25, _tier_ct25, lambda t: TIER_PALETTE.get(t, "#999"), "right")
        # Mini-tags under each paired bar so the two stories are unambiguous
        for _xi, _bc in zip(_xB, _bin_cols_25):
            axB.text(_xi - _offB, -1.5, "outcome", ha="center", va="top",
                     fontsize=8.5, color=_bc, style="italic", fontweight="bold")
            axB.text(_xi + _offB, -1.5, "tier", ha="center", va="top",
                     fontsize=8.5, color=_bc, style="italic", fontweight="bold")

        axB.set_ylim(0, 108)
        axB.set_xlim(-0.6, len(_oc_pct25) - 0.4)
        axB.set_xticks(_xB)
        axB.set_xticklabels(_xlabels_25, fontsize=9)
        for _tk, _bc in zip(axB.get_xticklabels(), _bin_cols_25):
            _tk.set_color(_bc)
            _tk.set_fontweight("bold")
        axB.tick_params(axis="x", pad=16)   # room for the outcome/tier mini-tags
        axB.set_ylabel("Percentage of complexes (%)", fontsize=10)
        axB.set_xlabel("PFAS chain-length bin", fontsize=10)
        axB.spines["top"].set_visible(False)
        axB.spines["right"].set_visible(False)
        axB.yaxis.grid(True, color="#ECECEC", linewidth=0.7, zorder=0)
        axB.set_axisbelow(True)

        from matplotlib.patches import Patch as _P25b
        _oc_h25 = [_P25b(facecolor=_oc25_colors_b.get(c, "#999"), edgecolor="white",
                         linewidth=0.5, label=c) for c in _oc_pct25.columns]
        _ti_h25 = [_P25b(facecolor=TIER_PALETTE.get(t, "#999"), edgecolor="white",
                         linewidth=0.5, label=t) for t in _tier_pct25.columns]
        """
        Two separate, clearly-titled legends so each bar maps to its OWN key —
        the left bar and the right bar use different category schemes (and some
        colours coincide), so a single merged legend was ambiguous.
        Both legends on ONE row: outcome (left bar) left-anchored, tier (right
        bar) right-anchored — saves vertical space.
        """
        _leg_oc25 = axB.legend(
            handles=_oc_h25, labels=[h.get_label() for h in _oc_h25],
            loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=len(_oc_h25),
            fontsize=8, title="Left bar — mechanistic outcome", title_fontsize=8,
            framealpha=0.92, fancybox=True, handlelength=0.9,
            handletextpad=0.35, columnspacing=0.7, borderpad=0.35)
        _leg_oc25.get_title().set_fontweight("bold")
        _leg_oc25._legend_box.align = "left"
        axB.add_artist(_leg_oc25)
        _leg_ti25 = axB.legend(
            handles=_ti_h25, labels=[h.get_label() for h in _ti_h25],
            loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=len(_ti_h25),
            fontsize=8, title="Right bar — degrader tier", title_fontsize=8,
            framealpha=0.92, fancybox=True, handlelength=0.9,
            handletextpad=0.35, columnspacing=0.7, borderpad=0.35)
        _leg_ti25.get_title().set_fontweight("bold")
        _leg_ti25._legend_box.align = "right"
        plt.tight_layout()
        _out25b = out_dir / "Figure_26b_PFAS_Size_Composition_Merged.png"
        fig25b.savefig(_out25b, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight",
                       bbox_extra_artists=(_leg_oc25, _leg_ti25))
        plt.close(fig25b)

        """
        Figure 26c: does FAcDs prefer smaller PFAS? The question is one of CATALYTIC
        competence, not model confidence, so the LEFT axis carries catalytic-competence
        readouts per CARBON-NUMBER group (x = C2…Cn, fluorine counts in parentheses):
          • soft_catalytic_score (box plots) — the continuous sigmoid composite of
            nucleophile reach, SN2 angle and triad relay; the distribution of mechanistic
            plausibility at each chain length.
          • degrader fraction (line) — the share of poses reaching a degrader tier
            (is_degrader), i.e. the share catalytically competent at that size.
        Boltz model confidence is plotted only as a CONTROL (faint grey median line): it
        measures structure-prediction self-certainty, not catalysis, and declines with
        ligand size as a prediction artefact (larger flexible chains are harder to place).
        The RIGHT axis carries median molecular weight as the size reference. Ligand
        nC/nF/MW are computed on the fly from the SMILES panel via RDKit
        (utils.compute_ligand_properties); nothing is hardcoded.
        """
        from matplotlib.patches import Patch as _P25
        _smi25 = out_dir.parents[1] / CFG.INPUT_SMILES
        if not _smi25.exists():
            _smi25 = Path.cwd() / CFG.INPUT_SMILES
        _LP25  = _utils_mod.compute_ligand_properties(_smi25)
        if not _LP25:
            reporter.log(f"  ! Figure 26c skipped: ligand SMILES not found ({_smi25})")
            raise RuntimeError("ligand properties unavailable")
        _lcol25 = next((c for c in [CFG.COL_LIG, "Ligand_Name", "ligand"]
                        if c in df.columns), None)
        if _lcol25 is None:
            reporter.log("  ! Figure 26c skipped: ligand column not found in data")
            raise RuntimeError("ligand column unavailable")
        """
        Catalytic-competence columns: the continuous soft_catalytic_score (0–1) and the
        per-pose degrader flag. Confidence is the labelled control. mechanistic_score is
        deliberately not used here — it saturates near 1.0 within the degrader tiers, so
        it has no spread to show a size trend (it is a gate, not a graded competence axis).
        """
        _soft_c25 = "soft_catalytic_score" if "soft_catalytic_score" in df.columns else None
        _conf_c25 = "Boltz_Model_Confidence" if "Boltz_Model_Confidence" in df.columns else None
        _deg_c25  = "is_degrader" if "is_degrader" in df.columns else None
        if _soft_c25 is None and _deg_c25 is None:
            reporter.log("  ! Figure 26c skipped: no catalytic-competence column (soft_catalytic_score / is_degrader)")
            raise RuntimeError("competence column unavailable")
        """
        _d25 was column-subset to numeric requirements (no ligand/score columns);
        pull them back from the source df by the surviving row index.
        """
        _cols25 = [_lcol25] + [c for c in (_soft_c25, _conf_c25, _deg_c25) if c]
        _d25c = df.loc[_d25.index, _cols25].copy()
        _d25c["_nC"] = _d25c[_lcol25].map(lambda l: _LP25.get(str(l), {}).get("nC"))
        _d25c["_mw"] = _d25c[_lcol25].map(lambda l: _LP25.get(str(l), {}).get("mw"))
        if _soft_c25:
            _d25c["_soft"] = pd.to_numeric(_d25c[_soft_c25], errors="coerce")
        if _conf_c25:
            _d25c["_conf"] = pd.to_numeric(_d25c[_conf_c25], errors="coerce")
        if _deg_c25:
            _d25c["_deg"] = (_d25c[_deg_c25].astype(str).str.strip().str.lower()
                             .isin(["true", "1", "1.0", "yes"]))
        _d25c = _d25c.dropna(subset=["_nC", "_mw"])
        if _d25c.empty:
            reporter.log("  ! Figure 26c skipped: no ligand property matches in data")
        else:
            _groups25 = sorted(_d25c["_nC"].unique())
            _present25 = set(_d25c[_lcol25].astype(str))

            def _xlabel25(nc):
                _fs = sorted({_LP25[l]["nF"] for l in _LP25
                              if _LP25[l]["nC"] == nc and l in _present25})
                return f"C{int(nc)}\n(" + ", ".join(f"{f}F" for f in _fs) + ")"

            _labels25 = [_xlabel25(g) for g in _groups25]
            _pos25    = np.arange(len(_groups25))
            _mw_by    = [_d25c[_d25c["_nC"] == g]["_mw"].dropna().values for g in _groups25]

            def _grp25(col):
                return [_d25c[_d25c["_nC"] == g][col].dropna().values for g in _groups25]

            fig25c, axL = plt.subplots(figsize=(15, 6))
            _SOFT_C, _SOFT_F = "#238B45", "#A1D99B"   # soft catalytic score (green)
            _DEG_C           = "#6A51A3"              # degrader fraction (purple)
            _CTRL_C          = "#9E9E9E"              # confidence control (grey)

            _handles25 = []
            # LEFT axis: soft_catalytic_score distribution per carbon group
            if "_soft" in _d25c.columns:
                _sb = _grp25("_soft")
                _bp = axL.boxplot(_sb, positions=_pos25, widths=0.42,
                                  patch_artist=True, showfliers=False, manage_ticks=False,
                                  medianprops=dict(color=_SOFT_C, linewidth=1.8))
                for _b in _bp["boxes"]:
                    _b.set(facecolor=_SOFT_F, alpha=0.85, edgecolor=_SOFT_C, linewidth=1.0)
                for _w in _bp["whiskers"] + _bp["caps"]:
                    _w.set(color=_SOFT_C, linewidth=1.0)
                _handles25.append(_P25(facecolor=_SOFT_F, edgecolor=_SOFT_C,
                                       label="Soft catalytic score (competence distribution)"))
                # Explicit legend entry for the in-box line so it is not confused with the
                # connected mean below — the box centre-line is the MEDIAN, not the mean.
                from matplotlib.lines import Line2D as _L2D25
                _handles25.append(_L2D25([0], [0], color=_SOFT_C, lw=1.8,
                                         label="Median (box centre line)"))

                # LEFT axis: connected MEAN of the competence boxes (size trend).
                # Dark-orange dashed diamonds — deliberately distinct from the purple
                # degrader-fraction line and the green box bodies.
                _smean25 = [float(np.mean(v)) if len(v) else np.nan for v in _sb]
                _mln25, = axL.plot(_pos25, _smean25, color="#D94801", lw=2.2, ls="--",
                                   marker="D", markersize=7, markeredgecolor="white",
                                   markeredgewidth=0.8, zorder=7,
                                   label="Mean competence (connected)")
                _handles25.append(_mln25)

            # LEFT axis: degrader fraction per carbon group (catalytic competence rate)
            if "_deg" in _d25c.columns:
                _dfrac = [float(_d25c[_d25c["_nC"] == g]["_deg"].mean())
                          if len(_d25c[_d25c["_nC"] == g]) else np.nan for g in _groups25]
                _dln, = axL.plot(_pos25, _dfrac, color=_DEG_C, lw=2.4, marker="o",
                                 markersize=6, zorder=6, label="Degrader fraction (is_degrader)")
                _handles25.append(_dln)

            # LEFT axis: Boltz confidence — CONTROL only (not catalytic)
            if "_conf" in _d25c.columns:
                _cmed = [float(np.median(v)) if len(v) else np.nan for v in _grp25("_conf")]
                _cln, = axL.plot(_pos25, _cmed, color=_CTRL_C, lw=1.6, ls=(0, (5, 3)),
                                 marker="D", markersize=4, alpha=0.85, zorder=4,
                                 label="Boltz confidence (control — prediction artefact, not catalytic)")
                _handles25.append(_cln)

            # Competence + degrader fraction carry the size story; molecular weight is not
            # plotted (it rises monotonically with carbon number, so it would be redundant
            # with the x-axis).

            # Per-group molecule count (unique ligands), inside each container bar near the top
            for _xi, g in zip(_pos25, _groups25):
                _nmol = int(_d25c[_d25c["_nC"] == g][_lcol25].nunique())
                axL.text(_xi, 1.075, f"{_nmol} mol", ha="center", va="top",
                         fontsize=8, fontweight="bold", color="#333", zorder=7)

            axL.set_xticks(_pos25)
            axL.set_xticklabels(_labels25, fontsize=8)
            axL.set_xlabel("PFAS carbon number  (fluorine counts present shown in parentheses)",
                           fontsize=10)
            axL.set_ylabel("Catalytic competence  (soft score · degrader fraction, 0–1)", fontsize=10)
            axL.spines["top"].set_visible(False)
            axL.set_ylim(0.0, 1.10)
            axL.set_xlim(-0.6, len(_groups25) - 0.4)

            """
            Master container bar per carbon group: very light tint + faint border
            — kept subtle so the box plots stay the visual focus.
            """
            from matplotlib.colors import to_rgba as _to_rgba25c
            _bin_pal25  = list(CFG.PFAS_SIZE_BIN_COLOUR)
            _grp_cols25 = [_bin_pal25[i % len(_bin_pal25)] for i in range(len(_groups25))]
            for _xi, _gc in zip(_pos25, _grp_cols25):
                axL.bar(_xi, 1.10, 0.74, bottom=0.0,
                        color=_to_rgba25c(_gc, 0.04), edgecolor=_to_rgba25c(_gc, 0.28),
                        linewidth=0.7, zorder=0)

            """
            Competence grid: green solid horizontal lines, faint grey vertical lines.
            """
            axL.yaxis.grid(True, color="#4CAF50", alpha=0.30, linewidth=0.7, zorder=0)
            axL.xaxis.grid(True, color="#EEEEEE", linewidth=0.5, zorder=0)
            axL.set_axisbelow(True)
            axL.legend(handles=_handles25, loc="lower left", bbox_to_anchor=(0.0, 1.01),
                       ncol=len(_handles25), fontsize=7.5, framealpha=0.95,
                       columnspacing=1.0, handletextpad=0.5, borderaxespad=0.0)
            plt.tight_layout()
            _out25c = out_dir / "Figure_26c_PFAS_Carbon_Confidence.png"
            fig25c.savefig(_out25c, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close(fig25c)

    except Exception as e:
        reporter.log(f"  ! Figure 26 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig


# =============================================================================
# SECTION 4D: DIAGNOSTIC & MULTI-MODEL TREND FIGURES  (folder 07)
# =============================================================================

def generate_additional_figures(df: pd.DataFrame, out_dir: Path,
                                reporter: ReportManager) -> None:
    """
    Pocket-fit and multi-model diagnostic figures written to
    <Run>/3_Validation_Figures/07_Diagnostic_and_MultiModel_Trends.

    These panels expose the steric and consensus signals that the main publication
    suite does not surface directly: whether each ligand physically fits its pocket
    (active-site volume vs ligand volume, occupancy, fit rate) and how the per-job
    Boltz-2 multi-model agreement (model_degrader_consensus) tracks both confidence
    and the final degrader tier. They are intentionally diagnostic — sourced from the
    pocket-fit columns added by 02_Production (active_site_volume, ligand_volume,
    pocket_occupancy, fit_ratio, ligand_fits) — and follow the same house style as the
    main figures: CFG palettes, median markers, 95% CI overlays and per-tier n counts.

    Each panel is guarded independently so a single missing column or empty subset
    skips only that figure, never the whole folder.
    """
    from scipy import stats as _sc_stats
    reporter.section("Step 5/5 — Diagnostic & Multi-Model Trend Figures  [writes folder 07]")
    out_dir.mkdir(parents=True, exist_ok=True)

    # CFG-sourced colours (single source of truth — config §8).
    _fit_col   = CFG.OUTCOME_COLOUR.get("Substrate", "#0E7C7B")
    _nofit_col = CFG.OUTCOME_COLOUR.get("Potential Inhibitor", "#6E2C00")

    existing_tiers = [t for t in TIER_ORDER_LOGIC if t in df.get("degrader_tier", pd.Series()).unique()]

    def _num(col: str) -> pd.Series:
        """Coerce a column to numeric, returning an all-NaN series if absent."""
        if col not in df.columns:
            return pd.Series(np.nan, index=df.index)
        return pd.to_numeric(df[col], errors="coerce")

    def _tier_ticklabels(ax, tiers, stats):
        """House-style coloured tier tick labels carrying μ ± CI, median and n."""
        ax.set_xticks(range(len(tiers)))
        labels = []
        for t in tiers:
            if t in stats:
                s = stats[t]
                lab = f"{t}\nμ={s['mean']:.3f}  ±{s['ci95']:.3f}"
                if abs(s["med"] - s["mean"]) > 0.005:
                    lab += f"\nmed={s['med']:.3f}"
                lab += f"\nn={s['n']:,}"
            else:
                lab = t
            labels.append(lab)
        ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=7.5)
        for tick, t in zip(ax.get_xticklabels(), tiers):
            tick.set_color(TIER_PALETTE.get(t, "black"))
            tick.set_fontweight("bold")

    def _per_tier_stats(values_by_tier):
        out = {}
        for t, v in values_by_tier.items():
            v = pd.Series(v).dropna()
            if len(v) >= 2:
                out[t] = {
                    "mean": float(v.mean()),
                    "ci95": float(_sc_stats.sem(v) * _sc_stats.t.ppf(0.975, len(v) - 1)),
                    "med":  float(v.median()),
                    "n":    int(len(v)),
                    "vals": v,
                }
        return out

    def _spear_box(ax, x, y, loc="lower right"):
        """Annotate an axis with the Spearman ρ, p-value and n of a relationship."""
        _xy = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"),
                            "y": pd.to_numeric(y, errors="coerce")}).dropna()
        if len(_xy) < 10:
            return
        _rho, _p = _sc_stats.spearmanr(_xy["x"], _xy["y"])
        if not np.isfinite(_rho):
            return
        _ptxt = "p < 0.001" if _p < 0.001 else f"p = {_p:.3f}"
        _txt = f"Spearman ρ = {_rho:+.2f}\n{_ptxt}   n = {len(_xy):,}"
        _xa, _ha = (0.97, "right") if "right" in loc else (0.03, "left")
        _ya, _va = (0.03, "bottom") if "lower" in loc else (0.97, "top")
        ax.text(_xa, _ya, _txt, transform=ax.transAxes, ha=_ha, va=_va, fontsize=8.5,
                fontweight="bold", color="#2C3E50",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#AAB7B8",
                          alpha=0.92, linewidth=0.8), zorder=20)

    def _spear_str(x, y):
        """Spearman ρ / p / n as a one-line string for embedding in a legend title; None if <10 pts."""
        _xy = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"),
                            "y": pd.to_numeric(y, errors="coerce")}).dropna()
        if len(_xy) < 10:
            return None
        _rho, _p = _sc_stats.spearmanr(_xy["x"], _xy["y"])
        if not np.isfinite(_rho):
            return None
        _ptxt = "p < 0.001" if _p < 0.001 else f"p = {_p:.3f}"
        return f"Spearman ρ = {_rho:+.2f}   {_ptxt}   n = {len(_xy):,}"

    # ── Figure 01: Pocket volume vs ligand volume (steric fit boundary) ─────────
    """
    Active-site cavity volume (x, Å³) against ligand molecular volume (y, Å³),
    one point per complex, coloured by degrader tier. The y = x diagonal is the
    geometric fit boundary: points above it describe ligands larger than the
    pocket that nominally encloses them (steric overfill), points below fit with
    headroom. Marginal median guide lines summarise each axis.
    """
    try:
        asv, lv = _num("active_site_volume"), _num("ligand_volume")
        sub = pd.DataFrame({"asv": asv, "lv": lv, "tier": df.get("degrader_tier")}).dropna(subset=["asv", "lv"])
        sub = sub[(sub["asv"] > 0) & (sub["lv"] > 0)]
        if len(sub) >= 5:
            fig, ax = plt.subplots(figsize=(8.4, 7.0))
            for t in existing_tiers:
                d = sub[sub["tier"] == t]
                if d.empty:
                    continue
                ax.scatter(d["asv"], d["lv"], s=14, alpha=0.30, linewidths=0,
                           color=TIER_PALETTE.get(t, "#999999"), label=t, zorder=3)
            # Highlight where the top tier (Tier_1A) lands in steric space.
            _t1a = sub[sub["tier"] == CFG.TIER_TOP]
            if len(_t1a) >= 3:
                _t1c = TIER_PALETTE.get(CFG.TIER_TOP, "#1B9E77")
                ax.scatter(_t1a["asv"], _t1a["lv"], s=48, facecolor=_t1c,
                           edgecolors="black", linewidths=0.8, alpha=0.95, zorder=6,
                           label=f"{CFG.TIER_TOP} (highlighted)")
                from matplotlib.patches import Rectangle as _Rect17
                _x05, _x95 = _t1a["asv"].quantile([0.05, 0.95])
                _y05, _y95 = _t1a["lv"].quantile([0.05, 0.95])
                ax.add_patch(_Rect17((_x05, _y05), _x95 - _x05, _y95 - _y05, fill=False,
                                     edgecolor=_t1c, linewidth=1.8, linestyle="--", zorder=7))
                _mx1, _my1 = float(_t1a["asv"].median()), float(_t1a["lv"].median())
                _mocc1 = (_my1 / _mx1) if _mx1 else float("nan")
                ax.annotate(f"{CFG.TIER_TOP}: median {_my1:.0f} Å³ ligand\nin {_mx1:.0f} Å³ pocket  "
                            f"(occupancy ≈ {_mocc1:.2f})",
                            xy=(_mx1, _my1), xycoords="data",
                            xytext=(0.50, 0.30), textcoords="axes fraction",
                            fontsize=8, fontweight="bold", color=_t1c, ha="left", va="top",
                            bbox=dict(boxstyle="round,pad=0.22", fc="none", ec=_t1c,
                                      alpha=0.95, linewidth=1.0),
                            arrowprops=dict(arrowstyle="->", color=_t1c, lw=1.4))
            _lim = float(np.ceil(max(sub["asv"].quantile(0.99), sub["lv"].quantile(0.99)) / 250.0) * 250.0)
            ax.plot([0, _lim], [0, _lim], color="#444444", linestyle="--", linewidth=1.3,
                    alpha=0.8, zorder=4, label="Fit boundary (ligand = pocket)")
            ax.axvline(sub["asv"].median(), color="#888888", linestyle=":", linewidth=1.0, alpha=0.6, zorder=2)
            ax.axhline(sub["lv"].median(), color="#888888", linestyle=":", linewidth=1.0, alpha=0.6, zorder=2)
            ax.set_xlim(0, _lim); ax.set_ylim(0, _lim)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("Active-site cavity volume  (Å³)", fontsize=11)
            ax.set_ylabel("Ligand molecular volume  (Å³)", fontsize=11)
            _ovf_n = int((sub["lv"] > sub["asv"]).sum())
            _ovf = _ovf_n / len(sub) * 100
            ax.text(0.02, 0.98, f"{_ovf:.1f}% of complexes overfill the pocket "
                    f"(n={_ovf_n:,} above boundary)",
                    transform=ax.transAxes, ha="left", va="top", fontsize=8.5, style="italic",
                    color="#6E2C00", bbox=dict(boxstyle="round,pad=0.18", fc="#FBEEE6",
                    ec=_nofit_col, alpha=0.9, linewidth=0.8))
            # Two-row legend above the axes (tiers + fit boundary); the Spearman ρ/p/n
            # rides in the legend title so it sits with the key, not as a floating box.
            _ss01 = _spear_str(sub["asv"], sub["lv"])
            ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=5, fontsize=7.5,
                      framealpha=0.92, columnspacing=0.9, handletextpad=0.4, borderaxespad=0.0,
                      title=_ss01, title_fontsize=8)
            plt.tight_layout()
            _o = out_dir / "01_Pocket_vs_Ligand_Volume.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 01 skipped: insufficient pocket/ligand volume data")
    except Exception as e:
        reporter.log(f"  ! Diag 01 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 02: Pocket occupancy distribution per PFAS carbon number ──────────
    """
    Pocket occupancy (ligand volume / cavity volume) grouped by PFAS carbon number
    (with the fluorine counts present per group shown in parentheses, matching the
    Fig 26c convention). Occupancy rises with chain length; the violin+box shows the
    spread, the red diamond marks the group median, and the occupancy = 1.0 line
    flags the point at which a ligand fully saturates the modelled cavity. Ligand
    nC/nF are computed on the fly from the SMILES panel via RDKit
    (utils.compute_ligand_properties); nothing is hardcoded.
    """
    try:
        # Ligand nC/nF map from the SMILES panel (same source as Fig 26c).
        _smi02 = out_dir.parents[1] / CFG.INPUT_SMILES
        if not _smi02.exists():
            _smi02 = Path.cwd() / CFG.INPUT_SMILES
        _LP02 = _utils_mod.compute_ligand_properties(_smi02) if _smi02.exists() else {}
        occ = _num("pocket_occupancy")
        sub = pd.DataFrame({"occ": occ, "lig": df.get("Ligand_Name")}).dropna(subset=["occ", "lig"])
        sub["_nC"] = sub["lig"].map(lambda l: _LP02.get(str(l), {}).get("nC"))
        sub = sub.dropna(subset=["_nC"])
        if not _LP02:
            reporter.log("  ! Diag 02 skipped: ligand SMILES panel not found for carbon grouping")
        elif len(sub) >= 5 and sub["_nC"].nunique() >= 2:
            sub["_nC"] = sub["_nC"].astype(int)
            groups = sorted(sub["_nC"].unique())
            _present02 = set(sub["lig"].astype(str))

            def _xlab02(nc):
                _fs = sorted({_LP02[l]["nF"] for l in _LP02
                              if _LP02[l].get("nC") == nc and l in _present02})
                return f"C{int(nc)}\n(" + ", ".join(f"{int(f)}F" for f in _fs) + ")"

            labels = [_xlab02(g) for g in groups]
            fig, ax = plt.subplots(figsize=(max(9.0, 0.95 * len(groups)), 6.4))
            _occ_cols = [plt.cm.viridis(_v) for _v in np.linspace(0.12, 0.9, len(groups))]
            sns.violinplot(data=sub, x="_nC", y="occ", order=groups, ax=ax, cut=0,
                           inner=None, density_norm="width", linewidth=0.6, palette=_occ_cols)
            sns.boxplot(data=sub, x="_nC", y="occ", order=groups, ax=ax, width=0.22,
                        showfliers=False, boxprops=dict(facecolor="white", alpha=0.9),
                        medianprops=dict(color="#C0392B", linewidth=1.6),
                        whiskerprops=dict(linewidth=0.8), capprops=dict(linewidth=0.8))
            _med_xs02, _med_ys02 = [], []
            for i, g in enumerate(groups):
                m = sub.loc[sub["_nC"] == g, "occ"].median()
                _med_xs02.append(i); _med_ys02.append(float(m))
            ax.plot(_med_xs02, _med_ys02, color="#C0392B", linewidth=1.8, zorder=5, alpha=0.85)
            ax.scatter(_med_xs02, _med_ys02, marker="D", s=42, color="#C0392B",
                       edgecolors="black", linewidths=0.8, zorder=6, label="Group median")
            ax.axhline(1.0, color=_nofit_col, linestyle="--", linewidth=1.2, alpha=0.8,
                       label="Full cavity saturation (occupancy = 1.0)")
            # A handful of poses report spurious occupancy (cavity-detection failures give
            # a near-zero cavity volume → occupancy ≫ 1). Clamp the view to the real
            # distribution so the chain-length trend is legible; note any hidden outliers.
            _yhi02 = float(min(max(1.05, sub["occ"].quantile(0.995) * 1.2), 3.0))
            _nhi02 = int((sub["occ"] > _yhi02).sum())
            ax.set_ylim(0, _yhi02)
            if _nhi02:
                ax.text(0.985, 0.97, f"{_nhi02:,} pose(s) with occupancy > {_yhi02:.1f} "
                        f"hidden (cavity mis-detection)",
                        transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                        style="italic", color="#566573")
            ax.set_xticks(range(len(groups))); ax.set_xticklabels(labels, fontsize=8)
            ax.set_xlabel("PFAS carbon number  (fluorine counts present shown in parentheses)", fontsize=10)
            ax.set_ylabel("Pocket occupancy  (ligand vol / cavity vol)", fontsize=11)
            ax.legend(loc="upper left", fontsize=8, framealpha=0.92, ncol=1)
            plt.tight_layout()
            _o = out_dir / "02_Pocket_Occupancy_by_Carbon_Number.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 02 skipped: insufficient occupancy/carbon data")
    except Exception as e:
        reporter.log(f"  ! Diag 02 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 03: Occupancy vs catalytic competence (binned median trend) ──────
    """
    Pocket occupancy (x) against catalytic competence score (y), coloured by tier.
    A binned-median trend line (10 equal-width occupancy bins) summarises whether
    competence falls as ligands crowd the cavity — the steric penalty signal that
    feeds the feasibility weighting.
    """
    try:
        occ, comp = _num("pocket_occupancy"), _num("competence_score")
        sub = pd.DataFrame({"occ": occ, "comp": comp, "tier": df.get("degrader_tier")}).dropna(subset=["occ", "comp"])
        if len(sub) >= 10:
            fig, ax = plt.subplots(figsize=(8.4, 6.4))
            for t in existing_tiers:
                d = sub[sub["tier"] == t]
                if d.empty:
                    continue
                ax.scatter(d["occ"], d["comp"], s=12, alpha=0.22, linewidths=0,
                           color=TIER_PALETTE.get(t, "#999999"), label=t, zorder=2)
            # Clamp the x-view to the real occupancy distribution (a few cavity-detection
            # failures push occupancy ≫ 1 and otherwise crush all data against x = 0).
            _xhi03 = float(min(max(1.0, sub["occ"].quantile(0.995) * 1.1), 3.0))
            _nhi03 = int((sub["occ"] > _xhi03).sum())
            _bins = np.linspace(0.0, _xhi03, 11)
            sub["_b"] = pd.cut(sub["occ"].clip(upper=_xhi03), _bins, include_lowest=True)
            med = sub.groupby("_b", observed=True).agg(x=("occ", "median"), y=("comp", "median"),
                                                       n=("comp", "size")).dropna()
            med = med[med["n"] >= 5]
            if len(med) >= 2:
                ax.plot(med["x"].clip(upper=_xhi03), med["y"], color="#C0392B", linewidth=2.4,
                        marker="D", markersize=7, markeredgecolor="black", markeredgewidth=0.8,
                        zorder=6, label="Binned median trend")
            ax.set_xlim(-_xhi03 * 0.02, _xhi03)
            if _nhi03:
                ax.text(0.985, 0.03, f"{_nhi03:,} pose(s) with occupancy > {_xhi03:.1f} hidden",
                        transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5,
                        style="italic", color="#566573")
            ax.set_xlabel("Pocket occupancy  (ligand vol / cavity vol)", fontsize=11)
            ax.set_ylabel("Catalytic competence score", fontsize=11)
            _spear_box(ax, sub["occ"], sub["comp"], "upper right")
            ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=12, fontsize=7,
                      framealpha=0.92, columnspacing=0.9, handletextpad=0.4, borderaxespad=0.0)
            plt.tight_layout()
            _o = out_dir / "03_Occupancy_vs_Competence.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 03 skipped: insufficient occupancy/competence data")
    except Exception as e:
        reporter.log(f"  ! Diag 03 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 04: Ligand fit rate per ligand ───────────────────────────────────
    """
    Fraction of complexes per ligand for which the steric fit test passed
    (ligand_fits == True). Bars are sorted by fit rate and labelled with the
    percentage and complex count; the colour ramps from fit (teal) to no-fit
    (maroon) using the CFG outcome palette.
    """
    try:
        if "ligand_fits" in df.columns and "Ligand_Name" in df.columns:
            fits = df["ligand_fits"]
            if fits.dtype == object:
                fits = fits.astype(str).str.strip().str.lower().map(
                    {"true": True, "false": False, "1": True, "0": False})
            fits = fits.astype("boolean")
            sub = pd.DataFrame({"fit": fits, "lig": df["Ligand_Name"]}).dropna(subset=["fit", "lig"])
            if len(sub) >= 5 and sub["lig"].nunique() >= 2:
                g = sub.groupby("lig")["fit"].agg(["mean", "size"]).sort_values("mean", ascending=True)
                fig, ax = plt.subplots(figsize=(9.0, max(5.0, 0.34 * len(g))))
                import matplotlib.colors as _mc
                _cmap = _mc.LinearSegmentedColormap.from_list("fit", [_nofit_col, "#E8E8E8", _fit_col])
                colours = [_cmap(v) for v in g["mean"]]
                ax.barh(range(len(g)), g["mean"] * 100, color=colours, edgecolor="black", linewidth=0.5)
                ax.set_yticks(range(len(g))); ax.set_yticklabels(g.index, fontsize=8)
                from matplotlib.patches import Patch as _Patch
                ax.legend(handles=[_Patch(facecolor=_fit_col, edgecolor="black", label="High fit rate"),
                                   _Patch(facecolor=_nofit_col, edgecolor="black", label="Low fit rate")],
                          loc="upper right", fontsize=8, framealpha=0.92, ncol=2)
                for i, (rate, n) in enumerate(zip(g["mean"], g["size"])):
                    ax.text(rate * 100 + 1, i, f"{rate*100:.0f}%  (n={int(n):,})",
                            va="center", ha="left", fontsize=7.5, color="#222222")
                ax.set_xlim(0, 108)
                ax.set_xlabel("Complexes passing steric fit test  (%)", fontsize=11)
                ax.set_ylabel("Ligand", fontsize=11)
                plt.tight_layout()
                _o = out_dir / "04_Ligand_Fit_Rate_by_Ligand.png"
                fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
                reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
            else:
                reporter.log("  ! Diag 04 skipped: insufficient ligand-fit data")
        else:
            reporter.log("  ! Diag 04 skipped: ligand_fits/Ligand_Name column absent")
    except Exception as e:
        reporter.log(f"  ! Diag 04 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 05: Multi-model degrader consensus by tier ───────────────────────
    """
    Per-tier distribution of the Boltz-2 multi-model degrader consensus (fraction
    of the 5 diffusion models that independently called the complex a degrader),
    with strip overlay, tier mean (diamond), 95% CI bar and median annotation —
    the multi-model agreement signal underlying tier assignment.
    """
    try:
        cons = _num("model_degrader_consensus")
        sub = pd.DataFrame({"cons": cons, "tier": df.get("degrader_tier")}).dropna(subset=["cons"])
        tiers = [t for t in existing_tiers if (sub["tier"] == t).sum() >= 2]
        if len(sub) >= 5 and tiers:
            stats = _per_tier_stats({t: sub.loc[sub["tier"] == t, "cons"] for t in tiers})
            fig, ax = plt.subplots(figsize=(max(8.0, 1.1 * len(tiers)), 6.4))
            sns.stripplot(data=sub[sub["tier"].isin(tiers)], x="tier", y="cons", order=tiers,
                          palette=TIER_PALETTE, alpha=0.15, size=3.0, jitter=0.28, ax=ax, zorder=1)
            for i, t in enumerate(tiers):
                if t not in stats:
                    continue
                s = stats[t]
                ax.plot([i, i], [s["mean"] - s["ci95"], s["mean"] + s["ci95"]],
                        color="black", linewidth=2.2, zorder=4)
                ax.scatter([i], [s["mean"]], color=TIER_PALETTE.get(t, "#999"), s=180,
                           zorder=5, edgecolors="black", linewidths=1.2)
            ax.set_ylim(-0.02, 1.06)
            ax.set_xlabel("Degrader tier", fontsize=11)
            ax.set_ylabel("Multi-model degrader consensus  (fraction of 5 models)", fontsize=10)
            _tier_ticklabels(ax, tiers, stats)
            ax.yaxis.grid(True, color="#DCDCDC", linewidth=0.6, alpha=0.7, zorder=0)
            ax.set_axisbelow(True)
            _lh = [
                Line2D([0], [0], marker="o", color="none", markerfacecolor="#888888",
                       markeredgecolor="none", markersize=5, alpha=0.35, label="Individual complex"),
                Line2D([0], [0], marker="D", color="none", markerfacecolor="#888888",
                       markeredgecolor="black", markersize=8, markeredgewidth=1.2, label="Tier mean (◆)"),
                Line2D([0], [0], color="black", linewidth=2.0, label="95% confidence interval"),
            ]
            _mdhc = _md_ready_stars_cat(ax, df, tiers, "model_degrader_consensus")
            if _mdhc is not None:
                _lh.append(_mdhc)
            # Legend above the axes (single row) so it never sits on the top data band.
            ax.legend(handles=_lh, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=4,
                      fontsize=8.5, framealpha=0.92, fancybox=True,
                      columnspacing=1.2, handletextpad=0.5, borderaxespad=0.0)
            plt.tight_layout()
            _o = out_dir / "05_MultiModel_Consensus_by_Tier.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 05 skipped: insufficient consensus data")
    except Exception as e:
        reporter.log(f"  ! Diag 05 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 06: Model confidence by multi-model consensus level ──────────────
    """
    Boltz-2 model confidence distribution at each multi-model degrader-consensus level
    (fraction of the 5 diffusion models that independently called the complex a
    degrader). Consensus is a discrete fraction, so it is the categorical axis and the
    confidence spread within each level is shown as a violin + box with the level median.
    A rising median confirms that cross-model agreement tracks AI confidence; the per-
    level n and Spearman ρ quantify the relationship. This is the robustness check that a
    raw confidence-vs-consensus scatter cannot show (it collapses onto horizontal bands).
    """
    try:
        conf = _num("Boltz_Model_Confidence")
        cons = _num("model_degrader_consensus")
        sub = pd.DataFrame({"conf": conf, "cons": cons}).dropna()
        _levels = sorted(sub["cons"].unique())
        if len(sub) >= 20 and len(_levels) >= 2:
            _pos = {c: i for i, c in enumerate(_levels)}
            sub["_x"] = sub["cons"].map(_pos)
            fig, ax = plt.subplots(figsize=(8.8, 6.4))
            _order06 = list(range(len(_levels)))
            _conf_cols = [plt.cm.YlGnBu(_v) for _v in np.linspace(0.25, 0.9, len(_order06))]
            sns.violinplot(data=sub, x="_x", y="conf", order=_order06, ax=ax, cut=0,
                           inner=None, density_norm="width", linewidth=0.5,
                           palette=_conf_cols, zorder=2)
            sns.boxplot(data=sub, x="_x", y="conf", order=_order06, ax=ax, width=0.18,
                        showfliers=False, boxprops=dict(facecolor="white", alpha=0.9),
                        medianprops=dict(color="#C0392B", linewidth=1.6),
                        whiskerprops=dict(linewidth=0.8), capprops=dict(linewidth=0.8),
                        zorder=3)
            _meds06 = [sub.loc[sub["cons"] == c, "conf"].median() for c in _levels]
            ax.plot(_order06, _meds06, color="#C0392B", linewidth=1.8, alpha=0.85, zorder=4,
                    marker="D", markersize=6, markeredgecolor="black", markeredgewidth=0.7,
                    label="Level median")
            ax.set_xticks(_order06)
            ax.set_xticklabels([f"{c:.2g}\n(n={int((sub['cons'] == c).sum()):,})" for c in _levels],
                               fontsize=8)
            ax.set_xlabel("Multi-model degrader consensus  (fraction of 5 models)", fontsize=11)
            ax.set_ylabel("Boltz-2 model confidence", fontsize=11)
            _ss06 = _spear_str(sub["cons"], sub["conf"])
            ax.legend(loc="lower right", fontsize=8.5, framealpha=0.92,
                      title=_ss06, title_fontsize=8.5)
            plt.tight_layout()
            _o = out_dir / "06_Confidence_vs_Consensus.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 06 skipped: insufficient confidence/consensus data")
    except Exception as e:
        reporter.log(f"  ! Diag 06 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 07: Quality & competence diagnostics (2-panel scatter) ───────────
    """
    Two tier-coloured scatter panels relating the AI, geometric and competence signals
    per complex:
      A. Boltz model confidence (x) vs catalytic competence (y) — whether confident
         predictions are also catalytically competent, or whether AI certainty and
         catalysis diverge.
      B. Mechanistic score (x, geometric/machinery quality) vs catalytic competence
         (y) — whether better active-site geometry translates into higher competence,
         the coupling that drives the within-tier ranking. Both axes are continuous
         0–1 scores, so the relationship reads as a clean cloud + binned-median trend.
    """
    fig = None
    try:
        _qc_panels = []
        if "Boltz_Model_Confidence" in df.columns and "competence_score" in df.columns:
            _qc_panels.append(("Boltz_Model_Confidence", "competence_score",
                               "Boltz model confidence", "Catalytic competence score", None))
        _func_c = "mechanistic_score" if "mechanistic_score" in df.columns else (
            "soft_catalytic_score" if "soft_catalytic_score" in df.columns else None)
        if _func_c and "competence_score" in df.columns:
            _func_lbl = ("Mechanistic score (geometry / machinery)"
                         if _func_c == "mechanistic_score"
                         else "Soft catalytic score (geometry / machinery)")
            _qc_panels.append((_func_c, "competence_score",
                               _func_lbl, "Catalytic competence score", None))
        if _qc_panels:
            fig, axes = plt.subplots(1, len(_qc_panels),
                                     figsize=(6.6 * len(_qc_panels), 5.6), squeeze=False)
            for _axi, (_xc, _yc, _xl, _yl, _mode) in zip(axes[0], _qc_panels):
                _xv = _num(_xc)
                _yv = _num(_yc)
                if _mode == "invert":
                    _yv = (1.0 - _yv).clip(0, 1)
                _pdat = pd.DataFrame({"x": _xv, "y": _yv, "tier": df.get("degrader_tier")}).dropna(subset=["x", "y"])
                for _t in existing_tiers:
                    _d = _pdat[_pdat["tier"] == _t]
                    if _d.empty:
                        continue
                    _axi.scatter(_d["x"], _d["y"], s=16, alpha=0.35, linewidths=0,
                                 color=TIER_PALETTE.get(_t, "#999999"), label=_t, zorder=3)
                # Binned-median trend + relationship statistic for analytical depth.
                if len(_pdat) >= 20:
                    _qb = np.linspace(_pdat["x"].min(), _pdat["x"].max(), 11)
                    _pdat["_b"] = pd.cut(_pdat["x"], _qb, include_lowest=True)
                    _qm = _pdat.groupby("_b", observed=True).agg(
                        x=("x", "median"), y=("y", "median"), n=("y", "size")).dropna()
                    _qm = _qm[_qm["n"] >= 5]
                    if len(_qm) >= 2:
                        _axi.plot(_qm["x"], _qm["y"], color="#C0392B", linewidth=2.2,
                                  marker="D", markersize=6, markeredgecolor="white",
                                  markeredgewidth=0.7, zorder=6, label="Binned median")
                # Panel-specific Spearman sits in the panel title (top, with the key row)
                # rather than as a floating stat box inside the data.
                _ss07 = _spear_str(_pdat["x"], _pdat["y"])
                if _ss07:
                    _axi.set_title(_ss07, fontsize=9, fontweight="bold", color="#2C3E50", pad=6)
                _axi.set_xlabel(_xl, fontsize=10)
                _axi.set_ylabel(_yl, fontsize=10)
            # Single-row legend spanning the top of the figure (above both panels).
            _h07, _l07 = axes[0][0].get_legend_handles_labels()
            fig.legend(_h07, _l07, loc="lower left", bbox_to_anchor=(0.02, 0.99),
                       ncol=12, fontsize=7.5, framealpha=0.92, columnspacing=1.0,
                       handletextpad=0.4)
            plt.tight_layout()
            _o = out_dir / "07_Quality_and_Competence_Diagnostics.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 07 skipped: quality/competence columns absent")
    except Exception as e:
        reporter.log(f"  ! Diag 07 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig
    finally:
        if fig is not None:
            plt.close(fig)

    # ── Figure 08: Size preference — competence distribution & viability vs ligand size ──
    """
    FAcD is a small-substrate (haloacetate) hydrolase. Single-panel size-preference
    summary in the Figure-26c visual grammar: per-size-bin distribution of the effective
    mechanistic score (green boxes), the connected bin means (dark-orange size trend), the
    Tier_1A viability floor, and — on the right axis — the catalytic hit-rate (fraction
    reaching Tier_2A+) and the mean catalytic-shell pocket_containment that drives the
    demotion. The reliable size axis is ligand_max_extent with catalytic-anchored
    pocket_containment, NOT the convex-hull pocket_occupancy proxy.
    """
    try:
        from matplotlib.patches import Patch as _P08
        _ext = _num("ligand_max_extent"); _eff = _num("mechanistic_score_effective")
        _con = _num("pocket_containment")
        _torder = ["Tier_5_Decoy", "Tier_4", "Tier_3", "Tier_2B", "Tier_2A", "Tier_1B", "Tier_1A"]
        _tv = df["degrader_tier"].astype(str).map({t: i for i, t in enumerate(_torder)})
        sub = pd.DataFrame({"ext": _ext, "eff": _eff, "con": _con, "tv": _tv}).dropna(subset=["ext", "eff", "tv"])
        if len(sub) >= 50:
            _hi = float(np.nanpercentile(sub["ext"], 99))
            _bins = np.linspace(float(sub["ext"].min()), _hi, 13)
            sub["_b"] = pd.cut(sub["ext"], _bins, include_lowest=True)
            _NFLOOR08 = 30
            _grpo = [g for g, d in sub.groupby("_b", observed=True) if len(d) >= _NFLOOR08]
            _pos = np.arange(len(_grpo))
            _cent = [iv.mid for iv in _grpo]
            _boxd = [sub.loc[sub["_b"].eq(g), "eff"].dropna().values for g in _grpo]
            _mean = [float(np.mean(v)) if len(v) else np.nan for v in _boxd]
            _i2a = _torder.index("Tier_2A")
            _hit = [float((sub.loc[sub["_b"].eq(g), "tv"] >= _i2a).mean()) for g in _grpo]
            _cmn = [float(sub.loc[sub["_b"].eq(g), "con"].mean()) for g in _grpo]

            fig, axL = plt.subplots(figsize=(15, 6))
            from matplotlib.colors import to_rgba as _rgba08
            # Horizontal tier-quality zones on the effective-mechanistic-score axis, coloured
            # by CFG.TIER_COLOUR and bounded by the CFG.TIER_MECH_MIN gates (single source of
            # truth). Only Tier_1A/1B (0.85) and Tier_2A (0.70) are mechanistic-score-gated;
            # Tier_2B and below are geometry/relay-gated, so the sub-0.70 band is a neutral
            # "below-2A" zone rather than a fabricated mech floor.
            _f1a = CFG.TIER_MECH_MIN["Tier_1A"]; _f2a = CFG.TIER_MECH_MIN["Tier_2A"]
            for _lo, _hi_z, _zc in [(_f1a, 1.06, CFG.TIER_COLOUR["Tier_1A"]),
                                    (_f2a, _f1a, CFG.TIER_COLOUR["Tier_2A"]),
                                    (0.0,  _f2a, CFG.TIER_COLOUR["Tier_2B"])]:
                axL.axhspan(_lo, _hi_z, color=_rgba08(_zc, 0.09), zorder=0)
            # Boxes: effective-mechanistic-score distribution per size bin (green).
            _SOFT_C, _SOFT_F = "#238B45", "#A1D99B"
            _bp = axL.boxplot(_boxd, positions=_pos, widths=0.5, patch_artist=True,
                              showfliers=False, manage_ticks=False,
                              medianprops=dict(color=_SOFT_C, linewidth=1.7))
            for _b in _bp["boxes"]:
                _b.set(facecolor=_SOFT_F, alpha=0.85, edgecolor=_SOFT_C, linewidth=1.0)
            for _w in _bp["whiskers"] + _bp["caps"]:
                _w.set(color=_SOFT_C, linewidth=1.0)
            # Connected box means (dark-orange dashed diamonds — the size trend).
            _mln08, = axL.plot(_pos, _mean, color="#D94801", lw=2.2, ls="--", marker="D",
                               markersize=7, markeredgecolor="white", markeredgewidth=0.8,
                               zorder=7, label="Mean effective mechanistic score")
            # Both real mechanistic-score tier floors (per CFG.TIER_MECH_MIN), labelled on the
            # right where the tier zones are otherwise empty.
            for _fy, _flbl, _fc in [(_f1a, f"Tier_1A / 1B floor {_f1a:g}", CFG.TIER_COLOUR["Tier_1A"]),
                                    (_f2a, f"Tier_2A floor {_f2a:g}", CFG.TIER_COLOUR["Tier_2A"])]:
                axL.axhline(_fy, ls="--", color=_fc, lw=1.3, zorder=3)
                axL.text(len(_pos) - 0.55, _fy + 0.012, _flbl, ha="right", va="bottom",
                         fontsize=8, color=_fc, fontweight="bold", zorder=4)
            axL.set_ylim(0.0, 1.06); axL.set_xlim(-0.6, len(_pos) - 0.4)
            axL.set_ylabel("Effective mechanistic score  (box = distribution)", fontsize=10)
            axL.set_xlabel("Ligand max extent (Å)  →  larger", fontsize=10)
            axL.set_xticks(_pos); axL.set_xticklabels([f"{c:.1f}" for c in _cent], rotation=90, fontsize=8)
            axL.spines["top"].set_visible(False)
            axL.yaxis.grid(True, color="#4CAF50", alpha=0.28, linewidth=0.7, zorder=0)
            axL.set_axisbelow(True)
            # Right axis: catalytic hit-rate + mean pocket containment (both 0–1).
            axR = axL.twinx()
            _hln08, = axR.plot(_pos, _hit, color="#2c7fb8", lw=2.0, marker="o", markersize=6,
                               zorder=6, label="Catalytic hit-rate (fraction ≥ Tier_2A)")
            _cln08, = axR.plot(_pos, _cmn, color="#1B9E9E", lw=1.8, ls=(0, (5, 3)), marker="s",
                               markersize=5, zorder=5, label="Mean pocket containment")
            axR.set_ylim(0.0, 1.06); axR.grid(False)
            axR.set_ylabel("Catalytic hit-rate  ·  pocket containment", fontsize=10)
            from matplotlib.lines import Line2D as _L2D08
            _mln08.set_label("Mean effective mech score")
            _hln08.set_label("Catalytic hit-rate (≥ Tier_2A)")
            _handles08 = [_P08(facecolor=_SOFT_F, edgecolor=_SOFT_C,
                               label="Effective mech score (distribution)"),
                          _L2D08([0], [0], color=_SOFT_C, lw=1.7, label="Median (box centre line)"),
                          _mln08, _hln08, _cln08]
            axL.legend(handles=_handles08, loc="lower left", bbox_to_anchor=(0.0, 1.01),
                       ncol=5, fontsize=7.5, framealpha=0.95, columnspacing=1.0,
                       handletextpad=0.4)
            plt.tight_layout()
            _o = out_dir / "08_Size_Preference_Containment.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 08 skipped: insufficient size/tier data")
    except Exception as e:
        reporter.log(f"  ! Diag 08 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 09: Reactive-centre engagement with the 8 catalytic residues vs size ──
    """
    Why FAcD favours small substrates is NOT that large PFAS cannot reach the catalytic
    machinery — they can (the reactive carbon still lands within van-der-Waals attack
    distance of the nucleophile). It is that the FULL catalytic constellation cannot fire
    at once: reactive-C→nucleophile contact AND a productive SN2 backside angle AND the
    effective-mechanistic viability floor must hold simultaneously, and that collapses with
    chain length. Left axis (Å): per-carbon-number distribution of the reactive-C →
    catalytic-nucleophile distance (green boxes + connected mean) inside the 2.5–3.5 Å
    productive attack-distance window, over the IQR spread of all eight catalytic-residue
    distances (999-sentinels dropped). Right axis (0–1): the strict 'properly positioned'
    fraction (collapses, tracks the catalytic hit-rate), the loose 'reaches nucleophile'
    fraction (stays high — necessary, not sufficient), mean effective mechanistic score,
    and the CFG.TIER_MECH_MIN floors. Carbon number is computed from the SMILES panel via
    utils.compute_ligand_properties; nothing is hardcoded.
    """
    try:
        from matplotlib.patches import Patch as _P09
        from matplotlib.lines import Line2D as _L2D09
        from matplotlib.ticker import MultipleLocator as _MLoc09
        _D8_09 = ["Dist_Nucleophile", "Dist_Base", "Dist_Acid", "Dist_Clamp1", "Dist_Clamp2",
                  "Dist_Stabiliser_H", "Dist_Stabiliser_W", "Dist_Stabiliser_Y"]
        _lcol09 = next((c for c in [CFG.COL_LIG, "Ligand_Name", "ligand"] if c in df.columns), None)
        _need09 = _D8_09 + ["SN2_Attack_Angle", "mechanistic_score_effective", "degrader_tier"]
        if _lcol09 is None or not all(c in df.columns for c in _need09):
            reporter.log("  ! Diag 09 skipped: reactive-geometry / ligand columns absent")
        else:
            _smi09 = out_dir.parents[1] / CFG.INPUT_SMILES
            if not _smi09.exists():
                _smi09 = Path.cwd() / CFG.INPUT_SMILES
            _LP09 = _utils_mod.compute_ligand_properties(_smi09)
            if not _LP09:
                reporter.log(f"  ! Diag 09 skipped: ligand SMILES not found ({_smi09})")
                raise RuntimeError("ligand properties unavailable")
            _PROD_LO, _PROD_HI, _READY, _SENT = 2.5, 3.5, 3.5, 20.0
            _ANG_MIN = 150.0
            _MECH_MIN = CFG.TIER_MECH_MIN["Tier_2A"]
            _d09 = df.copy()
            _d09["_nC"] = _d09[_lcol09].map(lambda l: _LP09.get(str(l), {}).get("nC"))
            _d8_09 = _d09[_D8_09].apply(pd.to_numeric, errors="coerce").mask(lambda s: s >= _SENT)
            _dnuc09 = _d8_09["Dist_Nucleophile"]
            _ang09 = _num("SN2_Attack_Angle")
            _eff09 = _num("mechanistic_score_effective")
            _torder09 = ["Tier_5_Decoy", "Tier_4", "Tier_3", "Tier_2B", "Tier_2A", "Tier_1B", "Tier_1A"]
            _tv09 = _d09["degrader_tier"].astype(str).map({t: i for i, t in enumerate(_torder09)})
            _dnraw09 = pd.to_numeric(_d09["Dist_Nucleophile"], errors="coerce")
            _ready09 = ((_dnraw09 <= _READY) & (_ang09 >= _ANG_MIN) & (_eff09 >= _MECH_MIN)).astype(float)
            _reach09 = (_dnraw09 <= _READY).astype(float)
            _sub09 = pd.DataFrame({"nC": _d09["_nC"], "dnuc": _dnuc09, "eff": _eff09,
                                   "tv": _tv09, "ready": _ready09, "reach": _reach09}).dropna(subset=["nC"])
            _grps09 = sorted(_sub09["nC"].unique())
            if len(_grps09) < 3:
                reporter.log("  ! Diag 09 skipped: insufficient carbon-number groups")
            else:
                _pos09 = np.arange(len(_grps09))
                _box09 = [_sub09.loc[_sub09.nC.eq(g), "dnuc"].dropna().values for g in _grps09]
                _shlo = [np.nanpercentile(_d8_09[_d09["_nC"].eq(g)].values, 25) for g in _grps09]
                _shhi = [np.nanpercentile(_d8_09[_d09["_nC"].eq(g)].values, 75) for g in _grps09]
                _rdy09 = [_sub09.loc[_sub09.nC.eq(g), "ready"].mean() for g in _grps09]
                _rch09 = [_sub09.loc[_sub09.nC.eq(g), "reach"].mean() for g in _grps09]
                _mm09 = [_sub09.loc[_sub09.nC.eq(g), "eff"].mean() for g in _grps09]
                _hit09 = [float((_sub09.loc[_sub09.nC.eq(g), "tv"] >= _torder09.index("Tier_2A")).mean())
                          for g in _grps09]

                fig, axL = plt.subplots(figsize=(15, 8.2))
                _BAND09 = "#B15FBF"     # 8-residue spread — distinct visible violet
                axL.fill_between(_pos09, _shlo, _shhi, color=_BAND09, alpha=0.22, zorder=0,
                                 label="8-residue distance spread (IQR)")
                axL.axhspan(_PROD_LO, _PROD_HI, color="#2ca02c", alpha=0.10, zorder=0)
                for _pe in (_PROD_LO, _PROD_HI):   # clear boundaries for the attack-distance window
                    axL.axhline(_pe, ls=(0, (5, 3)), color="#1a7a1a", lw=1.1, alpha=0.6, zorder=1)
                axL.text(len(_pos09) - 0.42, (_PROD_LO + _PROD_HI) / 2,
                         "productive attack distance (2.5–3.5 Å)", ha="center", va="center",
                         rotation=90, fontsize=7.5, color="#1a7a1a", fontweight="bold")
                _SOFT_C09, _SOFT_F09 = "#238B45", "#A1D99B"
                _bp09 = axL.boxplot(_box09, positions=_pos09, widths=0.42, patch_artist=True,
                                    showfliers=False, manage_ticks=False,
                                    medianprops=dict(color=_SOFT_C09, linewidth=1.7))
                for _b in _bp09["boxes"]:
                    _b.set(facecolor=_SOFT_F09, alpha=0.9, edgecolor=_SOFT_C09, linewidth=1.0)
                for _w in _bp09["whiskers"] + _bp09["caps"]:
                    _w.set(color=_SOFT_C09, linewidth=1.0)
                _bmean09 = [float(np.nanmean(b)) if len(b) else np.nan for b in _box09]
                _mbl09, = axL.plot(_pos09, _bmean09, color="#0B5345", lw=2.0, ls="--", marker="D",
                                   ms=6, mec="white", mew=0.7, zorder=6, label="Mean nuc. distance")
                axL.set_ylim(2.0, 5.0); axL.set_xlim(-0.6, len(_pos09) - 0.15)
                axL.yaxis.set_major_locator(_MLoc09(0.25))
                axL.set_ylabel("Reactive-carbon → catalytic-residue distance (Å)", fontsize=10.5)
                axL.set_xlabel("PFAS carbon number  →  larger ligand", fontsize=10)
                axL.set_xticks(_pos09); axL.set_xticklabels([f"C{int(g)}" for g in _grps09], fontsize=8)
                axL.spines["top"].set_visible(False)
                axL.yaxis.grid(True, color="#4C8BC9", alpha=0.28, lw=0.6, zorder=0)   # left (Å) grid
                axL.set_axisbelow(True)

                axR = axL.twinx()
                _rln09, = axR.plot(_pos09, _rdy09, color="#7B1FA2", lw=2.6, marker="o", ms=6,
                                   zorder=8, label="Properly positioned")
                _rcl09, = axR.plot(_pos09, _rch09, color="#9E9E9E", lw=1.5, ls=(0, (2, 2)),
                                   marker=".", ms=6, zorder=5, label="Reaches nucleophile")
                _mln09, = axR.plot(_pos09, _mm09, color="#D94801", lw=2.2, ls="--", marker="D",
                                   ms=6, mec="white", mew=0.7, zorder=6, label="Mean mech score")
                _hln09, = axR.plot(_pos09, _hit09, color="#2c7fb8", lw=2.0, ls=(0, (4, 2)),
                                   marker="s", ms=5, zorder=6, label="Hit-rate (≥Tier_2A)")
                for _fy, _fc, _lbl in [(CFG.TIER_MECH_MIN["Tier_1A"], CFG.TIER_COLOUR["Tier_1A"], "Tier 1A/1B"),
                                       (CFG.TIER_MECH_MIN["Tier_2A"], CFG.TIER_COLOUR["Tier_2A"], "Tier 2A")]:
                    axR.axhline(_fy, ls=":", color=_fc, lw=1.4, zorder=3)
                    axR.text(0.12, _fy + 0.004, f"{_lbl} floor ≥{_fy:g}", ha="left", va="bottom",
                             fontsize=7.5, color=_fc, fontweight="bold", zorder=9,
                             bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=_fc, lw=0.5, alpha=0.88))
                axR.set_ylim(0, 1.02); axR.set_ylabel("Catalytic metric (0–1)", fontsize=10)
                axR.yaxis.set_major_locator(_MLoc09(0.1))
                axR.yaxis.grid(True, color="#E08A3C", alpha=0.30, lw=0.6, ls=(0, (4, 3)), zorder=0)  # right grid
                _leg09 = [_P09(facecolor=_SOFT_F09, edgecolor=_SOFT_C09, label="Nuc. distance (box)"),
                          _L2D09([0], [0], color=_SOFT_C09, lw=1.7, label="Median"),
                          _mbl09,
                          _P09(facecolor=_BAND09, alpha=0.22, label="8-residue spread (IQR)"),
                          _rln09, _rcl09, _mln09, _hln09]
                axL.legend(handles=_leg09, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=8,
                           fontsize=6.3, framealpha=0.95, columnspacing=0.8, handletextpad=0.35)
                plt.tight_layout()
                _o = out_dir / "09_Reactive_Engagement.png"
                fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
                reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
    except Exception as e:
        reporter.log(f"  ! Diag 09 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig


def generate_ramachandran_figures(prod_dir: Path, out_dir: Path, reporter: ReportManager):
    """
    Backbone-geometry validation of the control predictions against the 3R3U crystal.
    Renders the crystal Ramachandran plot plus, for each Boltz-2 control (DeHa4 and 3R3U
    sequences), a standalone plot and a crystal-overlay comparison. All plots are written
    to <Run>/3_Validation_Figures/01_Ramachandran. The 3R3U crystal PDB is read from the
    input folder produced by 02_Production; nothing is downloaded here.
    """
    import gemmi
    reporter.section("Step 3/5 — Ramachandran Backbone-Geometry Validation  [writes folder 01]")
    # `out_dir` is already the 01_Ramachandran folder (created in main); use it directly.
    rama_dir = out_dir
    rama_dir.mkdir(parents=True, exist_ok=True)

    crystal = prod_dir / "1_Input_Data" / CFG.REFERENCE_PDB_FILE
    if not crystal.exists():
        reporter.log(f"  ! Skipped Ramachandran: 3R3U crystal not found at {crystal}")
        return

    dpi = int(getattr(CFG, "VIS_FIGURE_DPI", 300))

    def _load(path):
        st = gemmi.read_structure(str(path))
        st.setup_entities()
        return st

    def _control_cif(token):
        jobs_dir = prod_dir / "4_Prediction_Jobs"
        if not jobs_dir.exists():
            return None
        for job in sorted(jobs_dir.glob(f"*{token}*Fluoroacetate*")):
            bc = sorted(job.glob("Best_Complex/*.cif"))
            if bc:
                return bc[0]
            pred = sorted(job.glob("boltz_results_*/predictions/**/*model_0.cif"))
            if pred:
                return pred[0]
        return None

    rama_ref = _utils_mod.compute_ramachandran_angles(_load(crystal))
    _utils_mod.save_ramachandran_plot(rama_ref, "3R3U Crystal Structure",
                                      rama_dir / "Ramachandran_3R3U_Crystal.png", dpi=dpi)
    reporter.log("  ✔ Saved: Ramachandran_3R3U_Crystal.png")

    for label, token in [("DeHa4", "DeHa4_Control"), ("3R3U", "3R3U_Control")]:
        cif = _control_cif(token)
        if not cif:
            reporter.log(f"  ! Ramachandran: no {label} control CIF found")
            continue
        rama_con = _utils_mod.compute_ramachandran_angles(_load(cif))
        _utils_mod.save_ramachandran_plot(rama_con, f"{label} Control (Boltz-2)",
                                          rama_dir / f"Ramachandran_{label}_Control.png", dpi=dpi)
        _utils_mod.save_ramachandran_comparison(rama_ref, rama_con,
                                                "3R3U (Crystal)", f"{label} (Boltz-2)",
                                                rama_dir / f"Ramachandran_{label}_vs_Crystal.png", dpi=dpi)
        reporter.log(f"  ✔ Saved: Ramachandran_{label}_Control.png + comparison")


# =============================================================================
# SECTION 5: FIGURE DESCRIPTIONS & REPORTING
# =============================================================================

def write_figure_descriptions(out_dir: Path):
    """
    Writes a human-readable log file describing every figure produced by this pipeline.

    Entries follow the on-disk folder-by-folder layout:
      01_Ramachandran/                          backbone-geometry validation of controls
      02_Dataset_and_Alignment_Overview/        dataset + sequence-alignment overview
      03_AI_Confidence_Quality/                 Boltz-2 confidence metrics
      04_Catalytic_Geometry_and_Mechanism/      structural + mechanistic geometry
      05_Ligand_Interactions_and_Chemical_Space/ interaction profile + chemical space
      06_PFAS_Scope_and_Synthesis/              multi-metric synthesis + publication assembly
      07_Diagnostic_and_MultiModel_Trends/      pocket-fit + multi-model consensus
    """
    lines = [
        "=" * 80,
        "BOLTZ-2 VALIDATION FIGURES — DESCRIPTION LOG",
        "=" * 80,
        "",
        "This file describes every figure saved in this directory.",
        "Figures are ordered by scientific narrative (basic → advanced).",
        "Each entry lists the filename, chart type, axis meanings, and what to look for.",
        "",
        "=" * 80,
        "PART 1 — DATASET OVERVIEW",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 01 — Figure_01_Active_Site_Residue_Mapping_Coverage.png",
        "  Title   : Active-site residue mapping coverage across all FAcD variants",
        "  Type    : Vertical bar chart — one bar per canonical catalytic residue",
        "  X-axis  : The 8 canonical FAcD active-site residues, ordered by mechanistic role:",
        "            nucleophile → acid/base catalysis → carboxylate clamp → fluoride pocket",
        "  Y-axis  : Number of variants in which the residue was resolved during alignment",
        "  Colour  : Functional role group (CFG.ACTIVE_SITE_ROLE_GROUP_COLOUR); residues of the",
        "            same role share a colour; dashed line = total variants",
        "  Companion: 01_Active_Site_Residue_Mapping_Missed.csv — per residue, the",
        "            sequence IDs that failed to map it (one column per residue).",
        "  Look for: Which catalytic machinery is universally alignable versus divergent.",
        "            The fluoride-pocket residues and the second clamp arm typically map",
        "            least often — the family diverges most there.",
        "",
        "-" * 80,
        "Figure 02 — Figure_02_Tier_Distribution.png",
        "  Title   : Distribution of Catalytic Tiers",
        "  Type    : Vertical bar chart — count + percentage annotations inside bars",
        f"  X-axis  : Tier name (ordered best→worst: {CFG.TIER_TOP}, {CFG.TIER_ORDER[1]}, {CFG.TIER_ORDER[2]}, Tier_2B, {CFG.TIER_ORDER[4]}, {CFG.TIER_POOR})",
        "  Y-axis  : Count of protein–ligand complexes",
        "  Inset   : Pie chart showing Boltz-2 diffusion model selection frequency (model_0–model_4)",
        "  Look for: The fraction of candidates reaching each quality level.",
        "            Inset reveals whether one Boltz model dominates (~model_0 ≈ 70%).",
        "",
        "-" * 80,
        "Figure 03 — Figure_03_Alignment_Grades.png",
        "  Title   : Sequence Identity Distribution — Alignment Grade Counts + KDE by Tier  [2-panel]",
        "  Panel A : Horizontal bar chart of protein counts per identity grade band",
        "            (Grade A = 90–100% identity → Grade I = <20%); bars labelled with count and %",
        "  Panel B : Kernel density estimate of sequence identity per tier — filled curves overlaid",
        "            Vertical dashed lines mark Grade A (>=90%), D (>=60%), F (>=40%) thresholds",
        "  Look for: Most proteins cluster in Grade D-E (50-70% identity) — evolutionarily distant",
        "            yet sharing the active-site fold. Panel B reveals whether higher-tier proteins",
        "            cluster at higher identity (more conserved sequence) compared to lower tiers.",
        "",
        "-" * 80,
        "Figure 04 — Figure_04_Tier_Grade_Distribution.png",
        "  Title   : Tier x Alignment Grade — Stacked 100% Horizontal Bar Chart",
        "  Type    : 100% stacked horizontal bars (one row per tier)",
        "  Y-axis  : Degrader tiers (each row normalised to 100%)",
        "  X-axis  : Grade bands (A [>=90%] to I [<20%]); 10% intervals with pastel background stripes",
        "  Colour  : Grade-specific green-to-brown ramp (same palette as Figure 03)",
        "  Annotations: Count + % for segments >=2.5%",
        "  Look for: Do higher-tier proteins cluster at higher sequence identity grades?",
        "            If Tier_1/Tier_2 rows are dominated by Grade A-C, conservation tracks quality.",
        "",
        "=" * 80,
        "PART 2 — AI PREDICTION QUALITY",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 05a — Figure_05a_AI_Quality_Assessment.png",
        "  Title   : AI Quality Assessment — Confidence by Tier + pTM vs. ipTM  [2-panel]",
        "  Panel A : Box plot of Boltz Model Confidence Score per tier (Y-axis zoomed >=0.70)",
        "            Median annotated above each box; green dashed = high (>=0.90); yellow dotted = acceptable (>=0.80)",
        "  Panel B : Scatter of pTM (overall fold confidence) vs. ipTM (interface confidence) per tier",
        "            Points above diagonal = interface confidence exceeds fold confidence (preferred)",
        "            Points below diagonal = fold more confident than interface (decoy risk)",
        "            Diamond markers show per-tier median positions",
        "  Look for: Panel A: higher tiers should show higher median confidence.",
        "            Panel B: best-tier points should sit above the diagonal and cluster top-right.",
        "",
        "-" * 80,
        "Figure 05b — Figure_05b_TT_AI_Quality_Space.png",
        f"  Title   : {CFG.TIER_TOP} AI Quality Space (hexbin density + structure thumbnails)",
        "  Type    : Hexbin density scatter of pTM vs. ipTM; inset PyMOL structure thumbnails",
        "  Axes    : X = pTM (global fold confidence, 0-1); Y = ipTM (interface confidence, 0-1)",
        f"  Stars   : Mark {CFG.TIER_TOP} PFAS positions in the pTM/ipTM confidence space",
        f"  Insets  : PyMOL-rendered protein-ligand structures for each {CFG.TIER_TOP} representative",
        f"  Look for: {CFG.TIER_TOP} stars in the top-right quadrant (high pTM AND high ipTM).",
        "            Structure thumbnails show the active-site geometry at a glance.",
        "",
        "-" * 80,
        "Figure 06 — Figure_06_pTM_vs_ipTM_by_Tier.png",
        "  Title   : pTM vs. ipTM — 2-D Scatter by Tier",
        "  Type    : Scatter with diagonal reference line (pTM = ipTM) + coloured deviation zones",
        "  X-axis  : pTM — overall protein fold confidence (0-1); 1.0 = perfect fold prediction",
        "  Y-axis  : ipTM — interface confidence (0-1); 1.0 = perfect interface prediction",
        "  Diagonal: pTM = ipTM reference; +/-0.03 band; above diagonal = preferred (interface > fold)",
        "  Colour  : Catalytic tier; diamond markers show per-tier median positions",
        "  Zones   : Green band above diagonal = interface-dominant (preferred); red band below = fold-dominant",
        "  Look for: Best-tier points cluster above the diagonal and in the top-right corner.",
        f"            {CFG.TIER_DECOY} complexes often fall below the diagonal (fold confident, interface uncertain).",
        "",
        "=" * 80,
        "PART 3 — STRUCTURAL VALIDATION",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 07 — Figure_07_ActiveSite_RMSD_by_Tier.png",
        "  Title   : Active Site RMSD vs. Reference Control by Tier",
        "  Type    : Violin + strip plot with coloured quality zones",
        "  X-axis  : Catalytic tier (ordered best to worst)",
        "  Y-axis  : RMSD of 8 active-site residues vs. characterised DeHa4 control structure (Angstrom)",
        "            Lower = more structurally conserved active site; clipped at 99th percentile.",
        "  Zones   : Green <=1.0 Ang (Excellent), Yellow 1.0-2.0 Ang (Acceptable), Red >2.0 Ang (Diverged)",
        "  NOTE    : Complexes above the 99th percentile RMSD are excluded as failed alignments.",
        "            The % of data within each zone is annotated on the figure.",
        "  Look for: Whether higher-tier proteins show lower RMSD (tighter active-site conservation).",
        "            A flat RMSD profile across tiers confirms fold conservation regardless of tier.",
        "",
        "-" * 80,
        "Figure 08 — Figure_08_Feature_Correlations.png",
        "  Title   : Feature Correlation Matrix (Spearman rho) with Significance",
        "  Type    : Triangular heatmap (lower triangle only); pairwise-complete Spearman",
        "  Cells   : Spearman rank correlation coefficient (-1 to +1) + significance stars",
        "  Stars   : *** p<0.001   ** p<0.01   * p<0.05  (BH-FDR on unique pairs; no star = n.s.)",
        "  Order   : Features grouped by curated family (outcome, geometry, affinity, binding)",
        "  Look for: Strong red cells = features that rise and fall together (redundant or causal).",
        "            Strong blue = features that are inversely related.",
        "            AI Confidence + Binding Prob. cluster together — same underlying signal.",
        "            Pareto Rank inverts against physics/AI scores (lower rank = better = higher score).",
        "",
        "-" * 80,
        "Figure 09 — Figure_09_Tier_Quality_DotPlot.png",
        "  Title   : Tier Quality Summary — Normalised Multi-metric Cleveland Dot Plot",
        "  Type    : Horizontal Cleveland dot plot; one metric per row, tiers as coloured dots",
        "  Y-axis  : 6 normalised metrics (0=worst, 1=best within each metric):",
        "            AI Confidence, Mech. Score, Binding Prob., Seq. Identity, SN2 Angle, RMSD (inverted)",
        "  X-axis  : Normalised score (0-1 scale, within-metric min-max)",
        "  Dots    : One per tier, sized ~210 pt; raw median value annotated above each dot",
        "  Lines   : Grey range lines show min-max span across tiers for each metric",
        "  Look for: Metrics with the widest spread most discriminate between tiers.",
        "            Tiers that consistently appear at X=1.0 excel across all dimensions.",
        "",
        "=" * 80,
        "PART 4 — MECHANISTIC ANALYSIS",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 10 — Figure_10_Mech_State_CrossTab.png",
        "  Title   : Mechanistic State Cross-Tab (Halide Stabilisation x Carboxylate Clamp by Tier)",
        "  Type    : Heatmap (count + row %) — rows = tiers; columns = mechanistic state combinations",
        "  States  : Stabilised/Unstabilised (TRP/TYR aromatic shield) x Clamped/Unclamped (ASP carboxylate)",
        "  Look for: Virtually all complexes across all tiers show Stabilised+Clamped (98-100%).",
        "            This confirms the DEHA4 active-site architecture is robustly conserved.",
        "            Any tier deviating from 100% warrants structural inspection.",
        "",
        "-" * 80,
        "Figure 11 — Figure_11_Mechanistic_Score_by_Tier.png",
        "  Title   : Mechanistic Score by Tier",
        "  Type    : Mean +/- 95% CI dot plot with jittered individual data points behind",
        "  X-axis  : Catalytic tier (ordered best to worst)",
        "  Y-axis  : Mechanistic Score (0-1; CFG-weighted anchor set [nucleophile reach,",
        "            triad relay distances, carboxylate clamp, halide stabilisation] plus a",
        "            graded SN2 attack-angle term — config §5.1; saturates at 1.0 only for a",
        "            complete anchor set at an ideal 180 deg trajectory)",
        f"  Zones   : Green (>={CFG.MECH_FP_BAND_STRONG:.2f} = strong), Yellow ({CFG.MECH_FP_BAND_MODERATE:.2f}-{CFG.MECH_FP_BAND_STRONG:.2f} = moderate), Red (<{CFG.MECH_FP_BAND_MODERATE:.2f} = weak)",
        "  Annotations: Large dot = mean; vertical bar = 95% CI; % label = fraction in strong zone",
        "  Look for: Strong-zone complexes carry full machinery AND a productive angle; a high",
        "            mean with a low strong-zone fraction signals good anchors but poor trajectories.",
        f"            {CFG.TIER_POOR} tier has the widest spread.",
        "",
        "-" * 80,
        "Figure 12 — Figure_12_SN2_Angle_by_Tier.png",
        "  Title   : SN2 Attack Angle — Empirical Cumulative Distribution (ECDF) by Tier",
        "  Type    : ECDF step plot per tier",
        "  X-axis  : SN2 attack angle (degrees); 180 deg = ideal linear nucleophilic back-attack",
        f"            Vertical dashed lines = tier-specific angle thresholds ({', '.join(f'{v:.0f}' for v in sorted(CFG.TIER_ANGLE_MIN.values()))} deg)",
        "  Y-axis  : Cumulative fraction of complexes in each tier with angle <= X (read as %)",
        "  Dots    : Median angle per tier",
        "  How to read: Y% of complexes in this tier have SN2 angle <= X degrees",
        "  Look for: Tier_1/Tier_2 tiers: ECDF curves shifted right (most complexes at high angles).",
        f"            {CFG.TIER_ORDER[4]}/{CFG.TIER_POOR}: curves shifted left. Tier thresholds show what % meet each criterion.",
        "",
        "-" * 80,
        "Figure 13a — Figure_13a_Mechanism_Geometry_Scatter.png",
        "  Title   : SN2 Attack Angle vs. Nucleophile-Ligand Distance",
        "  Type    : Scatter plot with tier threshold lines and ideal-zone shading",
        "  X-axis  : ASP110 nucleophile to electrophilic carbon distance (Ang); shorter = reaction-ready",
        "            Clipped to <=5 Ang to remove outliers with no catalytic contact.",
        "  Y-axis  : SN2 attack angle (degrees); 180 = perfect linear back-attack geometry",
        "  Shape   : Circle = halide-stabilised (TRP/TYR aromatic shield present); X = not stabilised",
        "  Dashed lines: Vertical = distance thresholds per tier; Horizontal = angle thresholds.",
        "  Shading : Green zone = ideal geometry (short distance AND high angle).",
        "  Look for: Tier_1/Tier_2 tier points cluster top-left (short distance + high angle).",
        f"            {CFG.TIER_ORDER[4]}/{CFG.TIER_POOR} tier points scatter widely — geometry less constrained.",
        "",
        "-" * 80,
        "Figure 13b — Figure_13b_TT_Mechanistic_Quality_Space.png",
        f"  Title   : {CFG.TIER_TOP} Mechanistic Quality Space (hexbin density + structure thumbnails)",
        "  Type    : Hexbin density scatter; inset PyMOL protein-ligand structure thumbnails",
        "  Axes    : X = Dist_ASP110 (Ang); Y = SN2_Attack_Angle (degrees)",
        f"  Stars   : Highlight {CFG.TIER_TOP} PFAS complexes within the mechanistic space",
        f"  Insets  : PyMOL-rendered active-site views for each {CFG.TIER_TOP} representative protein",
        f"  Look for: {CFG.TIER_TOP} stars clustered in the ideal mechanistic zone (top-left corner).",
        "            Insets confirm the geometry seen in data is reflected in the 3-D structure.",
        "",
        "=" * 80,
        "PART 5 — LIGAND INTERACTIONS",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 14a — Figure_14a_Molecular_Interaction_Profile.png",
        "  Title   : Molecular Interaction Profile — Bond Types + Fluorine Engagement  [2-panel]",
        "  Panel A : Grouped bar chart — bond type on X-axis, tier as colour hue",
        "            Bond types: H-Bond, Salt Bridge, Halogen, F-Polar, F-Hydrophobic, Hydrophobic",
        "            Y-axis = mean count per complex; bars annotated with count values",
        "  Panel B : Violin plot of Fluorine Engagement Ratio (FER) per tier",
        "            FER = interacting F atoms / total F atoms in the PFAS ligand  (0-1 scale)",
        "            Reference lines at 50% and 100% engagement",
        "            NOTE: violin KDE may extend slightly above 1.0 (kernel smoothing artefact; raw data <=1.0)",
        "  Look for: Halogen contacts dominate (PFAS ligands carry many F atoms).",
        "            Higher tiers should show richer F-Polar / F-Hydrophobic engagement.",
        "            Panel B: Tier_1/Tier_2 tiers engage a higher fraction of ligand fluorine.",
        "",
        "-" * 80,
        "Figure 14b — Figure_14b_TT_Interaction_Quality_Space.png",
        f"  Title   : {CFG.TIER_TOP} Interaction Quality Space (hexbin density + structure thumbnails)",
        "  Type    : Hexbin density scatter of interaction density vs. total interactions; inset PyMOL structures",
        "  Axes    : X = Interaction Density Norm (interactions per ligand heavy atom); Y = Total Interactions (num_interactions, all contact types)",
        f"  Stars   : Mark {CFG.TIER_TOP} PFAS positions in the interaction-quality 2-D space",
        f"  Insets  : PyMOL-rendered protein-ligand structures for each {CFG.TIER_TOP} representative",
        f"  Look for: {CFG.TIER_TOP} stars clustered top-right (dense, numerous active-site contacts).",
        "",
        "-" * 80,
        "Figure 15 — Figure_15_Fluorine_Engagement_by_Tier.png",
        "  Title   : Fluorine Engagement Ratio by Tier (box + strip + median trend line)",
        "  Type    : Dual-axis (box + strip plot on left; median trend line on right)",
        "  Left axis : Box + strip of Fluorine Engagement Ratio (FER) per tier",
        f"              Quality zones: Green >={CFG.ENGAGEMENT_BAND_HIGH:.2f} (high), Yellow {CFG.ENGAGEMENT_BAND_MODERATE:.2f}-{CFG.ENGAGEMENT_BAND_HIGH:.2f} (moderate), Red <{CFG.ENGAGEMENT_BAND_MODERATE:.2f} (low)",
        "  Right axis: Median FER trend line per tier (blue diamonds with 95% CI band)",
        "  Annotations: med=X, n=Y badge above each box",
        "  Look for: Do higher tiers engage a greater fraction of ligand fluorine atoms?",
        f"            A rising trend line from {CFG.TIER_POOR} to Tier_1A confirms tier-quality tracks F-engagement.",
        "",
        "=" * 80,
        "PART 6 — BINDING ENERGETICS",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 16 — Figure_16_Catalytic_Quality_vs_Inhibition.png",
        "  Title   : Catalytic quality vs active-site contact density per tier  [dual-axis]",
        "  Panel A : Violin plot of soft_catalytic_score (continuous SN2-geometry composite) per tier",
        "            A genuinely varying catalytic-quality metric (Binding Probability is omitted — it",
        "            saturates near 0.99 across all tiers and carries no per-tier structure)",
        "            Median annotated inside each violin",
        "  Panel B : Connected dot-line of per-tier mean active-site contact density (interaction_density)",
        "            on the right axis; shaded band = 95% CI around the mean",
        f"  Look for: Panel A: catalytic quality rises with tier; {CFG.TIER_POOR} sits lowest/broadest.",
        "            Panel B: contact density reflects active-site engagement. Product inhibition is",
        "            assessed downstream by the Step-06 MM-GBSA stage, not here.",
        "",
        "-" * 80,
        "Figure 17 — Figure_17_ActiveSite_Contact_Density_by_Tier.png",
        "  Title   : Active-site contact density by tier",
        "  Type    : Box + strip plot",
        "  X-axis  : Catalytic tier (ordered best to worst)",
        "  Y-axis  : interaction_density (active-site contacts per complex; higher = richer engagement)",
        "  Annotations: Median badge per box",
        "  Stats   : Kruskal-Wallis H-test p-value annotated",
        "  Look for: how active-site contact richness varies across tiers. Product inhibition is",
        "            assessed downstream by the Step-06 MM-GBSA stage, not here.",
        "",
        "=" * 80,
        "PART 7 — CHEMICAL SPACE",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 18a — Figure_18a_Chemical_Space_Map.png",
        "  Title   : Chemical Space Map (UMAP Manifold)",
        "  Type    : 2D scatter (UMAP dimensionality-reduced embedding)",
        "  Axes    : UMAP dimensions 1 & 2 — arbitrary units, no physical meaning",
        "  Colour  : Catalytic tier (same palette as all other figures)",
        "  Markers : Star markers highlight Hidden Gems (physics-good / AI-missed complexes)",
        "  Look for: Clusters of the same colour = structurally/chemically similar complexes.",
        "            Tier islands separated in space = distinct mechanistic families.",
        "            Hidden Gems isolated from their tier cluster may have unusual chemistry.",
        "",
        "-" * 80,
        "Figure 18b — Figure_18b_TT_Chemical_Space_Landscape.png",
        f"  Title   : {CFG.TIER_TOP} Chemical Space Landscape (KDE density + structure thumbnails)",
        "  Type    : KDE density contour overlay on UMAP scatter; inset PyMOL structure thumbnails",
        "  Axes    : UMAP dimensions 1 & 2",
        f"  Stars   : Gold/green star markers indicate {CFG.TIER_TOP} PFAS positions",
        f"  Insets  : PyMOL-rendered protein-ligand structures for each {CFG.TIER_TOP} representative",
        f"  Look for: Density ridgelines isolating {CFG.TIER_TOP} from lower-tier complexes.",
        "            Structure thumbnails reveal active-site geometry at a glance.",
        f"            If {CFG.TIER_TOP} forms a tight cluster, they share structural/chemical features.",
        "",
        "=" * 80,
        "PART 8 — MULTI-METRIC SYNTHESIS",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 19a — Figure_19a_Radar_TopHits.png",
        f"  Title   : Candidate Radar — Top-5 Hits vs. {CFG.TIER_POOR}-Tier Baseline  [Radar]",
        "  Type    : Radar / spider chart with normalised metric spokes (0-1 each)",
        "  Spokes  : 5 normalised metrics: AI Confidence, ipTM, Interaction Density, mean pLDDT, Binding Prob.",
        f"  Lines   : Top-5 best-tier hits (coloured lines) vs. {CFG.TIER_POOR}-tier average baseline (red dashed)",
        "  Look for: Hits forming large polygons outperform on multiple axes simultaneously.",
        "            Spokes where hits touch 1.0 = this metric is at its best possible value.",
        "            Spokes where baseline and hits overlap = metric does not discriminate.",
        "",
        "-" * 80,
        "Figure 19b — Figure_19b_Radar_TierReps.png",
        "  Title   : Candidate Radar — One Representative per Tier  [Radar]",
        "  Type    : Radar / spider chart; one line per tier (best Boltz confidence per tier selected)",
        "  Spokes  : Same 5 normalised metrics as Figure 19a",
        "  Look for: Which tier consistently achieves the largest polygon (best on all metrics)?",
        "            Which metrics most differentiate tiers (largest gap between lines)?",
        "",
        "-" * 80,
        "Figure 20a — Figure_20a_Tier_Success_Rates.png",
        "  Title   : Tier success rates — substrate vs inhibitor geometry",
        "  Type    : Grouped horizontal bar chart (one tier per row)",
        "  X-axis  : Percentage of complexes (%)",
        f"  Bars    : Substrate (SN2≥{CFG.SUBSTRATE_ANGLE_MIN:.0f}° + Conf≥{CFG.SUBSTRATE_CONF_MIN:.2f}); SN2≥{CFG.SUBSTRATE_ANGLE_MIN:.0f}° (any conf); Potential inhibitor (SN2<{CFG.INHIBITOR_ANGLE_MAX:.0f}°)",
        "  Look for: Top tiers carry higher substrate-geometry pass rates and lower inhibitor fractions.",
        "",
        "Figure 20b — Figure_20b_Conf_SN2_Landscape.png",
        "  Title   : Tier medians in Confidence × SN2-angle space",
        "  Type    : Scatter of per-tier median diamonds (IQR error bars) over a hexbin density background",
        "  X-axis  : Boltz Model Confidence",
        "  Y-axis  : SN2 Attack Angle (°)",
        f"  Zones   : Substrate zone (SN2≥{CFG.SUBSTRATE_ANGLE_MIN:.0f}°, green); Inhibitor zone (SN2<{CFG.INHIBITOR_ANGLE_MAX:.0f}°, orange); Conf threshold {CFG.SUBSTRATE_CONF_MIN:.2f}",
        "  Diamonds: Per-tier median position, coloured by tier, annotated with median SN2 + confidence.",
        "  Look for: Top tiers sit inside the substrate zone, right of the 0.75 confidence threshold.",
        "",
        "-" * 80,
        "Figure 21 — Figure_21_Conflict_Composition.png",
        "  Title   : Conflict Category x Tier Composition — Stacked Bar",
        "  Type    : Side-by-side stacked bar (absolute counts left, % composition right)",
        f"  X-axis  : Conflict category (Consensus High, Hidden Gem, Ambiguous, Consensus Low, {CFG.TIER_DECOY})",
        "  Y-axis left : Absolute number of complexes in each category x tier combination",
        "  Y-axis right: Percentage composition of each conflict category by tier (sums to 100%)",
        "  Colour  : Degrader tier (same palette as all other figures)",
        "  Definitions:",
        "    Consensus High  = top 50% on both physics AND AI confidence",
        f"    {CFG.TIER_DECOY}           = high AI confidence but weak physics geometry (false positive risk)",
        "    Hidden Gem      = strong physics but low AI confidence (under-estimated by AI)",
        "    Ambiguous       = intermediate on both axes",
        "    Consensus Low   = bottom 50% on both axes",
        "  Look for: Consensus High dominated by Tier_1/Tier_2 tiers = metrics agree.",
        f"            {CFG.TIER_DECOY}s concentrated in {CFG.TIER_POOR}/{CFG.TIER_ORDER[4]} = AI over-confident on weaker candidates.",
        "            The % panel shows tier composition normalised — compare categories fairly.",
        "",
        "-" * 80,
        "Figure 22 — Figure_22_Hidden_Gems_DeepDive.png",
        f"  Title   : Hidden Gems — Physics-{CFG.TIER_ORDER[4]} / AI-Missed Conflicts",
        "  Type    : Horizontal dual-dot lollipop slope chart (only generated if Hidden Gems exist)",
        "  X-axis  : Percentile rank (0-100%; 0 = worst, 100 = best)",
        "  Y-axis  : Each hidden gem complex (sorted by AI rank, best at top)",
        "  Left dot (blue circle)    : Physics Pareto rank percentile",
        "  Right dot (orange diamond): AI Ensemble rank percentile",
        "  Line width : Proportional to |gap| between physics and AI ranks",
        "  Background zones: Tertiles (0-33 lower, 33-67 middle, 67-100 top) with colour shading",
        "  Annotations: Delta% label right of rightmost dot; tier-coloured y-axis labels",
        "  Look for: Gems that rank high on physics (left dot far right) but low on AI (right dot left).",
        "            Wide connecting lines = large AI/physics disagreement — worth manual rescue.",
        "",
        "-" * 80,
        "Figure 23 — Figure_23_Category_Overlap_Euler.png",
        "  Title   : Category Overlap — Euler / Venn Diagram",
        f"  Type    : Manual circle patches (3 main circles + {CFG.TIER_TOP} dashed overlay)",
        "  Sets    :",
        "    A (Physics-Strong)  : soft_catalytic_score / mechanistic_score / SN2 angle >= 67th percentile (top tercile catalytic geometry; confidence-independent)",
        "    B (AI-Strong)       : Boltz_Model_Confidence / ipTM / pTM >= 67th percentile (top 33% AI)",
        "    C (Tier-Strong)     : degrader_tier in {" + CFG.TIER_TOP + ", " + CFG.TIER_ORDER[1] + ", " + CFG.TIER_ORDER[2] + "}",
        "  Overlays:",
        f"    Dashed pink circle  : {CFG.TIER_TOP} — weighted to their position across A/B/C",
        "    Dotted grey circle  : All-three intersection highlight",
        "    Star markers        : 3R3U and DeHa4 reference controls (positioned in correct region)",
        "  Annotations: 7 region count boxes (count + row %); 4-column legend",
        "  Look for: How many top-tier candidates excel across all three metrics simultaneously?",
        "            The A+B+C intersection = Consensus Best on physics, AI, AND catalytic tier.",
        "            Reference controls (3R3U, DeHa4) should land in the A+B+C region.",
        "",
        "-" * 80,
        "Figure 24a — Figure_24a_Top25_Multitarget_Proteins.png",
        "  Type    : Horizontal stacked bar chart (one bar per protein, stacked by tier)",
        "  X-axis  : Number of unique PFAS ligands degraded (best tier per protein–ligand pair)",
        "  Y-axis  : Protein (top 25, ranked by quality-weighted degradation breadth)",
        f"  Colours : tier palette ({CFG.TIER_TOP} green → {CFG.TIER_DECOY} grey)",
        "  Look for: Wide bars = broad-spectrum PFAS degraders. Deep-green left stack = top-tier activity across many ligands.",
        "",
        "-" * 80,
        "Figure 24b — Figure_24b_TopTier_Protein_PFAS_Breakdown.png",
        f"  Type    : Horizontal stacked bar; one bar per protein that reaches {CFG.TIER_TOP} on ≥1 PFAS",
        "  X-axis  : Number of PFAS degraded (each block = one PFAS ligand)",
        "  Blocks  : Colour = PFAS ligand identity; printed label = that ligand's tier (T1A, T1B, T2A, …)",
        f"  Y-axis  : Top-tier proteins (★×N = number of {CFG.TIER_TOP} ligands)",
        "  Look for: Whether the elite proteins also degrade many OTHER PFAS, and at what tiers (substrate breadth of the best hits).",
        "",
        "-" * 80,
        "Figure 25 — Figure_25_Sankey_Workflow.png",
        "  Title   : Sankey Pathway — full FAcD tier-decision rules → Final Tier",
        "  Type    : Flat Sankey (ribbon polygons); node boxes labelled with COUNTS",
        "  Columns : Nuc Distance (CFG cuts) · Carboxylate Clamp (Arg/Lys) ·",
        "            Halide Stabilisation · Triad Geometry (Nuc–Base · Base–Acid) ·",
        "            SN2 Angle · Mechanistic Score · Final Tier",
        f"  Order   : best category on top in every column ({CFG.TIER_TOP}, ≤{CFG.TIER_NUC_DIST['Tier_1A']:.1f} Å, ≥{CFG.TIER_ANGLE_MIN['Tier_1A']:.0f}°, ≥{CFG.MECH_ELITE_HI:.2f})",
        "  Bins    : all cut points sourced from CFG (TIER_NUC_DIST / TIER_ANGLE_MIN /",
        "            TIER_NB_MAX / TIER_BA_MAX / TIER_MECH_MIN) — match 02_Production gates",
        "  Ribbons : Width proportional to complex count; coloured by destination tier",
        f"  Look for: {CFG.TIER_TOP} requires ALL rules to pass (tight nuc + clamp + stabilisation",
        f"            + tight triad + ≥{CFG.TIER_ANGLE_MIN['Tier_1A']:.0f}° + mech ≥{CFG.MECH_ELITE_HI:.2f}); the narrowing chain shows the attrition.",
        "",
        "-" * 80,
        "Figure 26a — Figure_26a_PFAS_Size_Hexbin_Landscape.png",
        "  Title   : PFAS Chain-Length vs SN2 Geometry — Substrate Preference and Inhibition Risk",
        "  Type    : Hexbin density (F-count × SN2 angle) with tier scatter overlay and rolling median",
        "  X-axis  : Total fluorine count (proxy for carbon chain length); ticks every 2 units",
        f"  Y-axis  : SN2 Attack Angle (°); green zone ≥{CFG.SUBSTRATE_ANGLE_MIN:.0f}° = substrate; red zone <{CFG.INHIBITOR_ANGLE_MAX:.0f}° = inhibitor",
        f"  Look for: Rolling median SN2 trend across chain lengths; tier scatter reveals {CFG.TIER_TOP}",
        "            distribution relative to substrate zone.",
        "",
        "-" * 80,
        "Figure 26b — Figure_26b_PFAS_Size_Composition_Merged.png",
        "  Title   : Chain-length Composition — Mechanistic Outcome + Degrader Tier (merged)",
        "  Type    : Single panel; per chain-length bin a LEFT stacked bar (outcome) and",
        "            a RIGHT stacked bar (degrader tier) side by side.",
        "  Left bar: Substrate (green), Borderline (amber), Reactive low conf (blue),",
        "            Non-reactive (grey), Potential Inhibitor (orange-red)",
        f"  Right bar: Degrader-tier composition ({CFG.TIER_TOP} → {CFG.TIER_DECOY} palette)",
        "  Look for: Substrate % peaks at short/medium chains; inhibition risk rises with chain length.",
        "",
        "-" * 80,
        "Figure 26c — Figure_26c_PFAS_Carbon_Confidence.png",
        "  Title   : Catalytic competence vs molecular weight by PFAS carbon number",
        "  Type    : Dual-axis per carbon-number group — LEFT (0–1): soft_catalytic_score box + degrader fraction line; RIGHT: median MW",
        "  X-axis  : Carbon number C2…Cn; fluorine counts present per group in parentheses (e.g. C2 (1F, 2F, 3F))",
        "  Control : Boltz confidence median plotted as a faint grey dashed line — it falls with size as a prediction artefact, NOT catalysis",
        "  Look for: catalytic competence (soft score / degrader fraction) is roughly flat across chain length, whereas confidence declines with size — i.e. any apparent 'small-molecule preference' from confidence is a prediction artefact, not a catalytic one.",
        "",
        "=" * 80,
        "BACKBONE GEOMETRY — RAMACHANDRAN (subfolder: Ramachandran/)",
        "=" * 80,
        "",
        "-" * 80,
        "Ramachandran_3R3U_Crystal.png — phi/psi of the 3R3U crystal reference.",
        "Ramachandran_<Control>_Control.png — phi/psi of each Boltz-2 control prediction",
        "  (DeHa4 and 3R3U sequences against fluoroacetate).",
        "Ramachandran_<Control>_vs_Crystal.png — control overlaid on the crystal.",
        "  Look for: control backbones occupying the same favoured basins as the crystal,",
        "            confirming the predictions reproduce native secondary structure.",
        "",
        "=" * 80,
        "DIAGNOSTIC & MULTI-MODEL TRENDS (subfolder: 07_Diagnostic_and_MultiModel_Trends/)",
        "=" * 80,
        "",
        "-" * 80,
        "01_Pocket_vs_Ligand_Volume.png",
        "  Title   : Pocket capacity vs ligand size — steric fit boundary",
        "  Type    : Scatter (one point per complex), coloured by tier, equal-aspect axes",
        "  X-axis  : Active-site cavity volume (Å³); Y-axis: ligand molecular volume (Å³)",
        "  Boundary: y = x dashed line — points above it are ligands larger than the pocket",
        "  Look for: ligands sit well below the boundary (occupancy << 1); only the longest",
        "            perfluoro chains approach it. Annotation reports the overfill count.",
        "",
        "-" * 80,
        "02_Pocket_Occupancy_by_Ligand.png",
        "  Title   : Pocket occupancy by PFAS carbon number",
        "  Type    : Violin + box per carbon group (C2…Cn, fluorine counts in parentheses);",
        "            red diamond = group median, joined by a median trend line",
        "  Y-axis  : Pocket occupancy = ligand volume / cavity volume; dashed line at 1.0",
        "  Look for: occupancy rises monotonically with carbon number (longest chains highest,",
        "            fluoroacetate lowest); tails crossing 1.0 flag steric saturation.",
        "",
        "-" * 80,
        "03_Occupancy_vs_Competence.png",
        "  Title   : Catalytic competence vs pocket occupancy",
        "  Type    : Scatter coloured by tier + binned-median trend (red diamonds)",
        "  Axes    : X = pocket occupancy; Y = catalytic competence score",
        "  Look for: whether competence falls as the cavity crowds — the steric signal that",
        "            informs the feasibility weighting.",
        "",
        "-" * 80,
        "04_Ligand_Fit_Rate_by_Ligand.png",
        "  Title   : Steric fit rate by ligand",
        "  Type    : Horizontal bars, sorted; colour ramps no-fit (maroon) → fit (teal)",
        "  X-axis  : % of complexes passing the steric fit test (ligand_fits == True)",
        "  Look for: small ligands pass universally; fit rate erodes for the largest PFAS.",
        "",
        "-" * 80,
        "05_MultiModel_Consensus_by_Tier.png",
        "  Title   : Multi-model degrader consensus by tier",
        "  Type    : Strip + tier mean (diamond) + 95% CI; coloured tier ticks with μ/median/n",
        "  Y-axis  : Fraction of the 5 Boltz-2 diffusion models that independently called a degrader",
        f"  Look for: consensus decreases monotonically from {CFG.TIER_TOP} to {CFG.TIER_DECOY}, "
        "confirming tier ordering reflects cross-model agreement.",
        "",
        "-" * 80,
        "06_Confidence_vs_Consensus.png",
        "  Title   : AI confidence vs cross-model consensus",
        "  Type    : Density hexbin (log count) + binned-median trend",
        "  Axes    : X = Boltz-2 model confidence; Y = multi-model degrader consensus",
        "  Look for: high confidence usually coincides with model agreement; confident calls that",
        "            still split the ensemble are robustness outliers.",
        "",
        "-" * 80,
        "07_Quality_and_Competence_Diagnostics.png",
        "  Title   : Quality & competence diagnostics (2-panel scatter, tier-coloured)",
        "  Panel A : Boltz model confidence (x) vs catalytic competence (y) — whether confident",
        "            predictions are also catalytically competent, or whether the two diverge",
        "  Panel B : Mechanistic score (x, functional quality) vs chemical penalty (y, 1 − feasibility",
        "            or active-site contact flag) — the functional-vs-flag trade-off per complex",
        "  Look for: top-tier complexes clustering at high confidence + high competence (A) and high",
        "            functional score + low penalty (B). Ported from the HADs extra-validation set.",
        "",
        "=" * 80,
        "ADDED CATALYTIC / BINDING FIGURES (main suite)",
        "=" * 80,
        "",
        "-" * 80,
        "04_Catalytic_Geometry_and_Mechanism/09_Mechanistic_Fingerprint.png",
        "  Title   : Mechanistic fingerprint — per-tier catalytic feature profile",
        "  Type    : Radar / spider chart; one feature spoke per axis,",
        "            one tier-coloured filled polygon across them",
        "  Axes    : Mechanistic features (config §5) — carboxylate head, halide",
        "            stabilisation, carboxylate clamp, SN2 alignment (angle/180°),",
        "            mechanistic score (summary)",
        "  Radius  : Tier-mean engagement 0–1 (boolean requirements → fraction engaged;",
        "            SN2 angle normalised by 180°; continuous scores → clipped mean)",
        "  Look for: which catalytic requirements each tier satisfies and how the profile",
        "            degrades down the tiers; the top tier rides near 1.0 on every axis.",
        "            Complements the scalar Figure 11 (mech-score mean ± CI).",
        "",
        "-" * 80,
        "05_Ligand_Interactions_and_Chemical_Space/08_Binding_Energetics.png",
        "  Title   : Binding energetics — binding probability by tier",
        "  Type    : Binding-probability violin per tier (red median bar)",
        "  Look for: how binding probability (sigmoid of SN2 geometry, computed in Step 02) varies",
        "            across tiers. Product inhibition is assessed downstream by the Step-06 MM-GBSA stage.",
        "",
        "=" * 80,
        "PART 9 — TWO-CRITERIA TIER LOGIC & DEAD-END FEASIBILITY",
        "=" * 80,
        "",
        "-" * 80,
        "04_Catalytic_Geometry_and_Mechanism/10_TwoCriteria_Tier_Logic.png",
        "  Title   : Two-criteria tier logic + dead-end feasibility (single 3-panel figure)",
        "  Panel a : Criterion A (active-site integrity, fraction of the 8 catalytic residues correctly",
        "            placed) vs Criterion B (catalytic constellation score) violins per A-bin, with the",
        "            median-B trend and the Tier_1A Criterion-B floor. Top row = share of complexes in",
        "            each A-bin clearing the B-floor.",
        "  Panel b : Criterion-B ECDF per tier with per-tier mean diamonds (arrowed values) — higher",
        "            tiers shifted right; the constellation cleanly separates the tiers.",
        "  Panel c : SN2 dead-end chemistry gate — scissile C-F BDE vs backside steric occlusion for the",
        "            ligands that actually carry the penalty (plus FA/DFA controls); the shaded quadrant",
        "            is the BDE x occlusion dead-end gate.",
        "  Look for: B clears the floor only once A is complete (a); tier order tracks Criterion B (b);",
        "            trifluoroacetate sits in the dead-end quadrant while FA/DFA stay feasible (c).",
        "",
        "=" * 80,
    ]
    # Repoint every internal "Figure_NN_*.png" mention at its real folder-relative
    # path so the description log matches the on-disk 7-folder layout. Entries stay
    # in thematic reading order; each now carries its folder path.
    _fig_paths = {
        "Figure_01_Active_Site_Residue_Mapping_Coverage.png": "02_Dataset_and_Alignment_Overview/01_Active_Site_Residue_Mapping_Coverage.png",
        "Figure_02_Tier_Distribution.png": "02_Dataset_and_Alignment_Overview/02_Tier_Distribution.png",
        "Figure_03_Alignment_Grades.png": "02_Dataset_and_Alignment_Overview/03_Alignment_Grades.png",
        "Figure_04_Tier_Grade_Distribution.png": "02_Dataset_and_Alignment_Overview/04_Tier_Grade_Distribution.png",
        "Figure_05a_AI_Quality_Assessment.png": "03_AI_Confidence_Quality/01_AI_Quality_Assessment.png",
        "Figure_05b_TT_AI_Quality_Space.png": "03_AI_Confidence_Quality/02_AI_Quality_Space.png",
        "Figure_06_pTM_vs_ipTM_by_Tier.png": "03_AI_Confidence_Quality/03_pTM_vs_ipTM_by_Tier.png",
        "Figure_07_ActiveSite_RMSD_by_Tier.png": "04_Catalytic_Geometry_and_Mechanism/01_ActiveSite_RMSD_by_Tier.png",
        "Figure_08_Feature_Correlations.png": "04_Catalytic_Geometry_and_Mechanism/07_Feature_Correlations.png",
        "Figure_09_Tier_Quality_DotPlot.png": "04_Catalytic_Geometry_and_Mechanism/08_Tier_Quality_DotPlot.png",
        "Figure_10_Mech_State_CrossTab.png": "04_Catalytic_Geometry_and_Mechanism/02_Mech_State_CrossTab.png",
        "Figure_11_Mechanistic_Score_by_Tier.png": "04_Catalytic_Geometry_and_Mechanism/03_Mechanistic_Score_by_Tier.png",
        "Figure_12_SN2_Angle_by_Tier.png": "04_Catalytic_Geometry_and_Mechanism/04_SN2_Angle_by_Tier.png",
        "Figure_13a_Mechanism_Geometry_Scatter.png": "04_Catalytic_Geometry_and_Mechanism/05_Mechanism_Geometry_Scatter.png",
        "Figure_13b_TT_Mechanistic_Quality_Space.png": "04_Catalytic_Geometry_and_Mechanism/06_Mechanistic_Quality_Space.png",
        "Figure_14a_Molecular_Interaction_Profile.png": "05_Ligand_Interactions_and_Chemical_Space/01_Molecular_Interaction_Profile.png",
        "Figure_14b_TT_Interaction_Quality_Space.png": "05_Ligand_Interactions_and_Chemical_Space/02_Interaction_Quality_Space.png",
        "Figure_15_Fluorine_Engagement_by_Tier.png": "05_Ligand_Interactions_and_Chemical_Space/03_Fluorine_Engagement_by_Tier.png",
        "Figure_16_Catalytic_Quality_vs_Inhibition.png": "05_Ligand_Interactions_and_Chemical_Space/04_Catalytic_Quality_vs_Inhibition.png",
        "Figure_17_ActiveSite_Contact_Density_by_Tier.png": "05_Ligand_Interactions_and_Chemical_Space/05_ActiveSite_Contact_Density_by_Tier.png",
        "Figure_18a_Chemical_Space_Map.png": "05_Ligand_Interactions_and_Chemical_Space/06_Chemical_Space_Map.png",
        "Figure_18b_TT_Chemical_Space_Landscape.png": "05_Ligand_Interactions_and_Chemical_Space/07_Chemical_Space_Landscape.png",
        "Figure_19a_Radar_TopHits.png": "06_PFAS_Scope_and_Synthesis/01_Radar_TopHits.png",
        "Figure_19b_Radar_TierReps.png": "06_PFAS_Scope_and_Synthesis/02_Radar_TierReps.png",
        "Figure_20a_Tier_Success_Rates.png": "06_PFAS_Scope_and_Synthesis/03_Tier_Success_Rates.png",
        "Figure_20b_Conf_SN2_Landscape.png": "06_PFAS_Scope_and_Synthesis/04_Conf_SN2_Landscape.png",
        "Figure_21_Conflict_Composition.png": "06_PFAS_Scope_and_Synthesis/05_Conflict_Composition.png",
        "Figure_22_Hidden_Gems_DeepDive.png": "06_PFAS_Scope_and_Synthesis/06_Hidden_Gems_DeepDive.png",
        "Figure_23_Category_Overlap_Euler.png": "06_PFAS_Scope_and_Synthesis/07_Category_Overlap_Euler.png",
        "Figure_24a_Top25_Multitarget_Proteins.png": "06_PFAS_Scope_and_Synthesis/08_Top25_Multitarget_Proteins.png",
        "Figure_24b_TopTier_Protein_PFAS_Breakdown.png": "06_PFAS_Scope_and_Synthesis/09_TopTier_Protein_PFAS_Breakdown.png",
        "Figure_25_Sankey_Workflow.png": "06_PFAS_Scope_and_Synthesis/10_Sankey_Workflow.png",
        "Figure_26a_PFAS_Size_Hexbin_Landscape.png": "06_PFAS_Scope_and_Synthesis/11_PFAS_Size_Hexbin_Landscape.png",
        "Figure_26b_PFAS_Size_Composition_Merged.png": "06_PFAS_Scope_and_Synthesis/12_PFAS_Size_Composition_Merged.png",
        "Figure_26c_PFAS_Carbon_Confidence.png": "06_PFAS_Scope_and_Synthesis/13_PFAS_Carbon_Confidence.png",
        "02_Pocket_Occupancy_by_Ligand.png": "02_Pocket_Occupancy_by_Carbon_Number.png",
    }
    _txt = "\n".join(lines)
    for _old_fp, _new_fp in _fig_paths.items():
        _txt = _txt.replace(_old_fp, _new_fp)
    desc_path = _aux_dir(out_dir) / "05_Figure_Descriptions.txt"
    desc_path.write_text(_txt, encoding="utf-8")
    return desc_path


# =============================================================================
# SECTION 6: MAIN EXECUTION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Boltz-2 Master Validation & Dendrogram Framework",
        usage="%(prog)s <run_folder>  (e.g. Boltz-2_Run_20260309T085406Z)"
    )
    parser.add_argument("run", help="Name of the Boltz-2 run folder (e.g. Boltz-2_Run_20260309T085406Z)")
    args = parser.parse_args()

    root_dir = Path.cwd()
    run_ttth = root_dir / args.run
    if not run_ttth.exists():
        print(f"Error: run folder not found: {run_ttth}")
        sys.exit(1)

    _utils_mod.print_script_banner(
        "03_Validation_Figures_FAcDs.py",
        "Multi-Objective Ranking  ·  Pareto Frontiers  ·  Publication Figures",
    )
    print(f"  Run Name : {run_ttth.name}", flush=True)

    prod_dir = run_ttth / "1_Boltz2_Production"
    out_dir = run_ttth / "3_Validation_Figures"

    # Create subdirectories (seven-folder structure, 01–07)
    rama_dir = out_dir / "01_Ramachandran"
    diag_dir = out_dir / "07_Diagnostic_and_MultiModel_Trends"

    out_dir.mkdir(parents=True, exist_ok=True)
    rama_dir.mkdir(parents=True, exist_ok=True)
    diag_dir.mkdir(parents=True, exist_ok=True)

    np.random.seed(42)

    reporter = _make_reporter(out_dir)

    try:
        # -------------------------------------------------------------------------------
        # Pipeline Execution Phase
        # -------------------------------------------------------------------------------
        df, features = load_and_prep_data(prod_dir, reporter)
        df = perform_advanced_ranking(df, features, out_dir, reporter)
        df = analyse_conflicts(df, out_dir, reporter)

        final_csv = _aux_dir(out_dir) / CFG.FILE_VALIDATED_MASTER
        df.to_csv(final_csv, index=False)
        reporter.log(f"Final Validated Dataset Saved: {final_csv.resolve()}")

        # Figures are generated folder-by-folder in narrative order (01 → 07).
        # 01_Ramachandran — control backbone-geometry validation.
        generate_ramachandran_figures(prod_dir, rama_dir, reporter)

        # 02–06 main validation suite, then 07 diagnostics; the two-criteria figure
        # lands in folder 04 (fig 10), all routed into their numbered folders from inside this call.
        generate_comprehensive_figures(df, features, out_dir, reporter)

        # Cleanup temporary structure thumbnails
        for p in out_dir.rglob("_tt_thumbnails"):
            if p.is_dir():
                shutil.rmtree(p)

        desc_path = write_figure_descriptions(out_dir)
        reporter.log(f"Figure Descriptions Log Saved: {desc_path.resolve()}")

    except Exception as e:
        print(f"\n[CRITICAL ERROR] Pipeline Failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    import time as _time
    _t0 = _time.perf_counter()
    main()
    _utils_mod.print_elapsed(_t0, "03_Validation_Figures_FAcDs.py")
