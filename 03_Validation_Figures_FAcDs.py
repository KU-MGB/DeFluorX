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
Date   : 20 July 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 03_Validation_Figures_FAcDs.py
  Role          : Scoring, ranking, and visual reporting of Boltz-2 predictions.
  Imports from  : 00_01_Project_Config_FAcDs.py  (CFG — tier colours, vis params)
                  00_02_Project_Utils_FAcDs.py   (ConsoleColours, setup_logging,
                                                 console_info, console_separator)
  Reads         : <Run>/1_Boltz2_Production/*_Ranked_*.csv  (falls back to *_Master_*.csv)
                  (Master CSV written by 02_Production_FAcDs.py; latest file selected)
  Writes        : <Run>/3_Validation_Figures/  (figures grouped folder-by-folder)
                    02_Ramachandran/                        Ramachandran_*.png
                    03_Dataset_and_Alignment_Overview/      01_*.png … 04_*.png
                    04_AI_Confidence_Quality/               01_*.png … 03_*.png
                    05_Catalytic_Geometry_and_Mechanism/    01_*.png … 10_*.png
                    06_Ligand_Interactions_and_Chemical_Space/ 01_*.png … 08_*.png
                    07_PFAS_Scope_and_Synthesis/            01_*.png … 13_*.png
                    08_Diagnostic_and_MultiModel_Trends/    01_*.png … 09_*.png  (pocket-fit + consensus + competence)
                  <Run>/3_Validation_Figures/01_Analysis_Data/03_Figure_Enriched_Dataset.csv
                  <Run>/3_Validation_Figures/01_Analysis_Data/04_ACTION_Rescue_Hidden_Gems.csv
                  <Run>/3_Validation_Figures/01_Analysis_Data/05_Figure_Descriptions.txt
                  <Run>/3_Validation_Figures/01_Analysis_Data/00_Validation_Figures.log
  Upstream      : 02_Production_FAcDs.py → writes the master ranked CSV (incl. the pocket-fit
                  columns active_site_volume, ligand_volume, pocket_occupancy, fit_ratio,
                  ligand_fits) consumed here
  Downstream    : 04_Dendrogram_FAcDs.py  → reads 03_Figure_Enriched_Dataset.csv (figure columns
                  only; the authoritative rank stays in 02's ranked CSV)
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
    [Data — written to 01_Analysis_Data/]
    • 03_Figure_Enriched_Dataset.csv        <-- ranked CSV + figure columns (PCA/UMAP/Pareto/conflict); NOT a rank source
    • 04_ACTION_Rescue_Hidden_Gems.csv      <-- MANUAL REVIEW LIST
    • 06_Statistical_Tests.csv              <-- per-test BH q-values (multiplicity paid once)
    • 07_Boltz2_MultiModel_QC_Variance.csv  <-- per-Boltz-model geometry variance (figs 13–15)
    • 05_Figure_Descriptions.txt            <-- Per-figure description log
    • 00_Validation_Figures.log             <-- Detailed Execution Log

    [Figures — grouped folder-by-folder; each folder numbered in narrative order]
    ── 02_Ramachandran/ ── backbone-geometry validation of the controls
    • Ramachandran_3R3U_Crystal.png · Ramachandran_{DeHa4,3R3U}_Control.png (+ _vs_Crystal)

    ── 03_Dataset_and_Alignment_Overview/ ── dataset + sequence-alignment overview
    • 01_Active_Site_Residue_Mapping_Coverage.png  <-- Per-residue alignment coverage (+ MISSED.csv)
    • 02_Tier_Distribution.png                     <-- Catalytic tier counts + model-selection pie
    • 03_Alignment_Grades.png                      <-- Sequence-identity grades + KDE by tier
    • 04_Tier_Grade_Distribution.png               <-- Tier × alignment-grade cross-tabulation

    ── 04_AI_Confidence_Quality/ ── Boltz-2 confidence metrics
    • 01_AI_Quality_Assessment.png                 <-- Boltz confidence boxes + pTM/ipTM panel
    • 02_AI_Quality_Space.png                      <-- {CFG.TIER_TOP} pTM/ipTM density + thumbnails
    • 03_pTM_vs_ipTM_by_Tier.png                   <-- pTM vs ipTM 2-D scatter by tier

    ── 05_Catalytic_Geometry_and_Mechanism/ ── structural + mechanistic geometry
    • 01_ActiveSite_RMSD_by_Tier.png               <-- Active-site RMSD vs reference control (median bar)
    • 02_Mech_State_CrossTab.png                   <-- Halide stabilisation × carboxylate clamp heatmap
    • 03_Mechanistic_Score_by_Tier.png             <-- Mechanistic score mean ± CI (median annotated)
    • 04_SN2_Angle_by_Tier.png                     <-- SN2 attack-angle ECDF by tier
    • 05_Mechanism_Geometry_Scatter.png            <-- SN2 angle vs nucleophile-distance scatter
    • 06_Mechanistic_Quality_Space.png             <-- {CFG.TIER_TOP} mechanistic density + thumbnails
    • 07_Feature_Correlations.png                  <-- Spearman ρ correlation heatmap of 17 features (incl. mechanistic/chemistry drivers), family-grouped
    • 08_Tier_Quality_DotPlot.png                  <-- Multi-metric tier-quality Cleveland dot plot
    • 09_Mechanistic_Fingerprint.png               <-- Per-tier catalytic feature profile (radar / spider; config §5 components)
    • 10_Criterion_A_Gates_B.png                   <-- active-site integrity gates the catalytic constellation
    • 11_Criterion_B_ECDF_by_Tier.png             <-- the constellation score separates the tiers (ECDF)
    • 12_SN2_DeadEnd_Gate.png                     <-- SN2 dead-end chemistry gate (C–F BDE × backside occlusion)

    ── 06_Ligand_Interactions_and_Chemical_Space/ ── interaction profile + chemical space
    • 01_Molecular_Interaction_Profile.png         <-- Bond-type profile + fluorine engagement
    • 02_Interaction_Quality_Space.png             <-- {CFG.TIER_TOP} interaction density + thumbnails
    • 03_Fluorine_Engagement_by_Tier.png           <-- Fluorine engagement ratio box + trend line
    • 04_Catalytic_Quality_vs_Inhibition.png       <-- Soft-catalytic-score violin + active-site contact density
    • 05_ActiveSite_Contact_Density_by_Tier.png    <-- Active-site contact density box + strip
    • 06_Binding_Energetics.png                    <-- Binding-probability violin by tier
    • 07_Chemical_Space_Map.png                    <-- UMAP chemical-space manifold, competence hexbin + {CFG.TIER_TOP} structure thumbnails
    • 08_Binding_Affinity_Metrics.png              <-- Binding-affinity distribution by tier (2-column legend + stats)

    ── 07_PFAS_Scope_and_Synthesis/ ── multi-metric synthesis + publication assembly
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

    ── 08_Diagnostic_and_MultiModel_Trends/ ── pocket-fit + multi-model consensus diagnostics
    • 01_Pocket_vs_Ligand_Volume.png               <-- Cavity vs ligand volume, y=x steric fit boundary
    • 02_Pocket_Occupancy_by_Carbon_Number.png     <-- Whole-cavity coverage (violin) + 8-residue active-site coverage (box) per PFAS carbon number, both means trended
    • 03_MultiModel_Consensus_by_Tier.png          <-- Multi-model degrader consensus, mean ± CI by tier
    • 04_Confidence_vs_Consensus.png               <-- Confidence vs cross-model consensus hexbin + trend
    • 05_Quality_and_Competence_Diagnostics.png    <-- Confidence×competence + mech-score×penalty scatter (tier)
    • 06_Size_Preference_Containment.png           <-- Effective-mech distribution + means + hit-rate + pocket containment vs ligand size
    • 07_Reactive_Engagement.png                   <-- Reactive-C→catalytic-residue distance + properly-positioned fraction (vs hit-rate) by carbon number
    • 08_Model_Agreement.png                       <-- Diffusion-sample consensus by tier: are the elite hits reproducible across samples

    • 01_Geometry_and_Uncertainty.png              <-- Nucleophile distance + SN2 angle with multi-model uncertainty
    • 02_Binding_Affinity_Metrics.png              <-- Binding-affinity distribution by tier (2-column legend + stats)
    • 03_Evolutionary_Phylogeny.png       <-- Evolutionary phylogeny, detailed variant
    • 04_Pillar_Divergence_by_Tier.png             <-- Divergence of the scoring pillars across tiers
    • 05_Mechanistic_Breakdown_by_Tier.png         <-- Mechanistic components broken down per tier
    • 06_Chain_Length_by_Tier.png                  <-- PFAS chain-length distribution by tier (stats top-left)
    • 07_Tier1A_Cross_Ligand_Heatmap.png           <-- Tier_1A proteins × ligands cross-tabulation heatmap

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
       00_01_Project_Config_FAcDs.py — see that module's Scientific References
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
import re
import json
import shutil
import time as _time
import concurrent.futures as cf
import argparse
import warnings
from pathlib import Path
import contextlib
from dataclasses import dataclass
from typing import Any, Iterable, List, Optional, Sequence, Tuple

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
import gemmi
import scipy.stats as _sc_stats
from scipy.stats import spearmanr, gaussian_kde, chi2_contingency, mannwhitneyu, t as _t_dist
from sklearn.preprocessing import MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA

# Bioinformatics Imports
from matplotlib.patches import ConnectionPatch
from matplotlib.lines import Line2D
import matplotlib.patheffects as pe
from PIL import Image, ImageDraw, ImageFont


# -------------------------------------------------------------------------------
# Step 1.3: Pipeline modules (00_01) via importlib
# -------------------------------------------------------------------------------
import importlib.util as _ilu

def _load_module(name: str, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod  = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_cfg_mod   = _load_module("ProjectConfig", Path(__file__).resolve().parent / "00_01_Project_Config_FAcDs.py")
_utils_mod = _load_module("ProjectUtils",  Path(__file__).resolve().parent / "00_02_Project_Utils_FAcDs.py")
CFG        = _cfg_mod.CFG()

ConsoleColours  = _utils_mod.ConsoleColours
latest_by_mtime = _utils_mod.latest_by_mtime   # newest ranked/master CSV by mtime (prefix-agnostic)
SEPARATOR_HEAVY = _utils_mod.SEPARATOR_HEAVY
SEPARATOR_LIGHT = _utils_mod.SEPARATOR_LIGHT
SEPARATOR_DASH  = _utils_mod.SEPARATOR_DASH
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
    "Error": CFG.VIS_ACCENT["error"],                # 03-specific analytics-failure indicator (not in CFG)
}

# Specific order for tiers to ensure logical plotting (Tier_1A to Tier_5_Decoy)
TIER_ORDER_LOGIC = list(CFG.TIER_ORDER)

# Consistent vertical separator between tier / class groups on every per-tier figure (one style
# everywhere). Colour is CFG-sourced (VIS_INK) so it stays in the single source of truth.
_TIER_SEP_KW = dict(color=CFG.VIS_INK["mid"], linestyle=":", alpha=0.5, linewidth=1.0, zorder=0)


def _tier_seps(ax, n=None):
    """Vertical dotted separators between the n categorical groups on a per-tier x-axis (lines at
    0.5, 1.5, … n-1.5). One consistent look so tier/class boundaries read the same in every figure.
    n is inferred from the axis tick count when omitted; lines outside the current xlim are skipped."""
    if n is None:
        n = len(ax.get_xticks())
    try:
        _n = int(n)
    except (TypeError, ValueError):
        return
    _x0, _x1 = ax.get_xlim()
    for _i in range(_n - 1):
        _x = _i + 0.5
        if _x0 <= _x <= _x1:
            ax.axvline(_x, **_TIER_SEP_KW)


CONFLICT_PALETTE = dict(CFG.CONFLICT_COLOUR)   # sourced from CFG § 8.7


# =============================================================================
# Inferential-statistics helpers.
# Each returns a short annotation string; _stat_box() draws it inside an axis.
# =============================================================================
def _stat_box(ax, text, loc="upper right", fontsize=8):
    """Draw a boxed statistics annotation inside an axis (no effect if text is empty)."""
    if not text:
        return
    xy = {"upper right": (0.985, 0.985, "right", "top"),
          "upper left": (0.015, 0.985, "left", "top"),
          "lower right": (0.985, 0.015, "right", "bottom"),
          "lower left": (0.015, 0.015, "left", "bottom"),
          "lower centre": (0.5, 0.015, "center", "bottom")}.get(loc, (0.985, 0.985, "right", "top"))
    ax.text(xy[0], xy[1], text, transform=ax.transAxes, ha=xy[2], va=xy[3], fontsize=fontsize,
            family="monospace", zorder=30,
            bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=CFG.VIS_INK["paler"], alpha=0.92))


"""
Every hypothesis test drawn on a panel is registered here, and the family is corrected ONCE, at the
end of the run, by Benjamini-Hochberg. Several tests are run across the same candidate set (one per
panel), so an uncorrected p-value overstates significance: with a handful of tests at alpha = 0.05 a
false positive is expected by chance alone. The panel keeps the RAW p (it is what that test
measured); the corrected q-value for the whole family goes to 01_Analysis_Data/06_Statistical_Tests.csv,
which is the value to quote.
"""
_PVALUES: list = []


def _register_p(test: str, panel: str, stat: float, n: int, p: float, **extra) -> None:
    """Register one test into the family. extra carries effect size, group ns, medians — whatever the
    test can honestly report. A p-value with no effect size is half a result: at n = 58,056 almost
    anything is 'significant', and the number that decides whether it MATTERS is the effect size."""
    _row = {"test": test, "panel": panel, "statistic": stat, "n": n, "p_raw": p}
    _row.update(extra)
    _PVALUES.append(_row)


def _statistical_battery(df: pd.DataFrame, reporter) -> None:
    """The full statistical battery, registered into the same BH-corrected family as the figures.

    The file was only ever recording the handful of tests a panel happened to run, which left the
    reader to take the paper's central claims on trust. Every claim the ranking rests on is now tested
    explicitly and corrected together:

      · Kruskal-Wallis across tiers for each metric  — does the metric separate the tiers at all?
      · Mann-Whitney, degraders vs non-degraders     — with rank-biserial r, the effect SIZE. At
                                                       n = 58,056 a p-value is nearly free; r is not.
      · Elite (Tier_1A) vs everything else           — the specific claim the headline makes.
      · Spearman rho between the ranking metrics     — are the pillars independent, or restating each
                                                       other? A screen built on three correlated
                                                       metrics has one metric and two echoes.

    Every test lands in 06_Statistical_Tests.csv with its q_BH, so the multiplicity is paid for once,
    across the whole family, rather than per figure.
    """
    from scipy.stats import kruskal as _kw, mannwhitneyu as _mw, spearmanr as _sr

    _metrics = [c for c in (
        "mechanistic_score_effective", CFG.COL_MECH_S, "competence_score",
        CFG.COL_SN2, "sn2_attack_angle_effective", "Dist_Nucleophile",
        "Binding_Affinity_Score", "Interaction_Density_Norm", CFG.COL_CONF,
        "active_site_plddt", "catalytic_constellation_score", CFG.COL_ID_PCT,
        "pocket_containment_cavity", "pocket_containment_site8", "total_fluorine_count",
    ) if c in df.columns]
    if not _metrics or CFG.COL_TIER not in df.columns:
        return

    _tier = df[CFG.COL_TIER].astype(str)
    _hq = set(getattr(CFG, "TIER_HIGH_QUALITY", []) or
              [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]])
    _is_deg = _tier.isin(_hq)
    _is_elite = _tier.eq(CFG.TIER_TOP)

    def _rb(_a, _b, _u):
        """Rank-biserial r from the Mann-Whitney U: the probability that a random member of A exceeds
        a random member of B, rescaled to [-1, 1]. Unlike p, it does not inflate with n."""
        _n = len(_a) * len(_b)
        return (2.0 * _u / _n - 1.0) if _n else float("nan")

    for _m in _metrics:
        _v = pd.to_numeric(df[_m], errors="coerce")

        # 1) Does it separate the tiers at all?
        _groups = [_v[_tier == _t].dropna().to_numpy()
                   for _t in CFG.TIER_ORDER if (_tier == _t).sum() >= 3]
        _groups = [g for g in _groups if len(g) >= 3]
        if len(_groups) >= 2:
            try:
                _h, _p = _kw(*_groups)
                _n_tot = int(sum(len(g) for g in _groups))
                # epsilon-squared: the share of rank variance the tier explains.
                _eps2 = (float(_h) - len(_groups) + 1) / (_n_tot - len(_groups)) if _n_tot > len(_groups) else float("nan")
                _register_p("Kruskal-Wallis across tiers", _m, float(_h), _n_tot, float(_p),
                            effect_size=round(_eps2, 4), effect_type="epsilon^2",
                            n_groups=len(_groups))
            except Exception as _e:                              # noqa: BLE001
                # A test that fails to run must SAY SO. Swallowed silently, it simply vanishes from the
                # results file, and a missing test looks exactly like a test that was never wanted.
                reporter.log(f"    ! Kruskal-Wallis skipped for {_m}: {type(_e).__name__}: {_e}")

        # 2) Degraders vs non-degraders, with the effect size.
        for _lbl, _mask in (("Degraders vs non-degraders", _is_deg),
                            (f"{CFG.TIER_TOP} vs rest", _is_elite)):
            _a = _v[_mask].dropna().to_numpy()
            _b = _v[~_mask].dropna().to_numpy()
            if len(_a) < 3 or len(_b) < 3:
                continue
            try:
                _u, _p = _mw(_a, _b, alternative="two-sided")
                _register_p(f"Mann-Whitney U — {_lbl}", _m, float(_u), len(_a) + len(_b), float(_p),
                            effect_size=round(_rb(_a, _b, float(_u)), 4),
                            effect_type="rank-biserial r", n_group_a=len(_a), n_group_b=len(_b),
                            median_a=round(float(np.median(_a)), 4),
                            median_b=round(float(np.median(_b)), 4))
            except Exception as _e:                              # noqa: BLE001
                reporter.log(f"    ! Mann-Whitney skipped for {_m} ({_lbl}): {type(_e).__name__}: {_e}")

    # 3) Are the ranking metrics independent, or restating one another?
    _corr = [c for c in ("mechanistic_score_effective", "competence_score",
                         "Binding_Affinity_Score", CFG.COL_CONF,
                         "Interaction_Density_Norm", CFG.COL_SN2) if c in df.columns]
    for _i in range(len(_corr)):
        for _j in range(_i + 1, len(_corr)):
            _x = pd.to_numeric(df[_corr[_i]], errors="coerce")
            _y = pd.to_numeric(df[_corr[_j]], errors="coerce")
            _ok = _x.notna() & _y.notna()
            if int(_ok.sum()) < 20:
                continue
            try:
                _rho, _p = _sr(_x[_ok], _y[_ok])
                _register_p("Spearman correlation", f"{_corr[_i]} vs {_corr[_j]}",
                            float(_rho), int(_ok.sum()), float(_p),
                            effect_size=round(float(_rho), 4), effect_type="Spearman rho")
            except Exception as _e:                              # noqa: BLE001
                reporter.log(f"    ! Spearman skipped for {_corr[_i]} vs {_corr[_j]}: "
                             f"{type(_e).__name__}: {_e}")

    reporter.log(f"  Statistical battery: {len(_PVALUES)} tests registered (BH-corrected together)")


def _write_statistical_tests(out_dir):
    """The whole test family with Benjamini-Hochberg q-values, ranked by q."""
    if not _PVALUES:
        return None
    from scipy.stats import false_discovery_control
    _df = pd.DataFrame(_PVALUES)
    # A metric can be registered by both the central battery and a figure helper (e.g. the confidence
    # Kruskal-Wallis). Collapse identical (test, panel) rows before BH so one test counts once — a
    # duplicate would inflate the family denominator and write a repeated row to the CSV.
    # keep="last": on a (test, panel) collision the battery registers AFTER the figure helpers, and its
    # row carries the effect size and the stricter min-group (≥3) — so the battery's row must win.
    _df = _df.drop_duplicates(subset=["test", "panel"], keep="last").reset_index(drop=True)
    # A non-finite p (a Kruskal-Wallis on all-identical groups, a zero-variance Spearman) would poison the
    # q-values for the WHOLE family — false_discovery_control propagates the NaN. Set those aside (recorded
    # with q_BH = NaN) so the finite tests are corrected among themselves.
    _finite = np.isfinite(pd.to_numeric(_df["p_raw"], errors="coerce").to_numpy(float))
    _df["q_BH"] = np.nan
    if _finite.any():
        _df.loc[_finite, "q_BH"] = false_discovery_control(
            _df.loc[_finite, "p_raw"].to_numpy(float), method="bh")
    _df["significant_q<0.05"] = _df["q_BH"] < 0.05
    _df = _df.sort_values("q_BH").reset_index(drop=True)
    _path = _aux_dir(out_dir) / "06_Statistical_Tests.csv"
    _df.to_csv(_path, index=False)
    return _path


def _kruskal_by_tier(df, value_col, tier_col=CFG.COL_TIER, tiers=None):
    """Kruskal–Wallis across tiers with epsilon-squared effect size. Returns annotation or ''."""
    from scipy.stats import kruskal as _kw
    if value_col not in df.columns or tier_col not in df.columns:
        return ""
    order = tiers or [t for t in TIER_ORDER_LOGIC if t in set(df[tier_col].dropna())]
    groups = [pd.to_numeric(df.loc[df[tier_col] == t, value_col], errors="coerce").dropna().values for t in order]
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 2:
        return ""
    try:
        H, p = _kw(*groups)
    except Exception:
        return ""
    N = sum(len(g) for g in groups); k = len(groups)
    eps2 = max(0.0, (H - k + 1) / (N - k)) if N > k else float("nan")
    _register_p("Kruskal-Wallis across tiers", value_col, float(H), int(N), float(p),
                effect_size=(round(float(eps2), 4) if np.isfinite(eps2) else None), effect_type="epsilon^2")
    p_str = "p < 1e-4" if p < 1e-4 else f"p = {p:.3g}"
    return f"Kruskal-Wallis across tiers: {p_str} | eps^2 = {eps2:.2f}"


def _wilcoxon_ptm_iptm(df):
    """Paired Wilcoxon signed-rank test of ipTM vs pTM across complexes."""
    from scipy.stats import wilcoxon as _wil
    if "ptm" not in df.columns or "iptm" not in df.columns:
        return ""
    sub = df.dropna(subset=["ptm", "iptm"])
    if len(sub) < 10:
        return ""
    try:
        _w, p = _wil(pd.to_numeric(sub["ptm"], errors="coerce"), pd.to_numeric(sub["iptm"], errors="coerce"))
    except Exception:
        return ""
    med_d = float((pd.to_numeric(sub["iptm"], errors="coerce") - pd.to_numeric(sub["ptm"], errors="coerce")).median())
    _register_p("Wilcoxon signed-rank (paired)", "ipTM vs pTM", float(_w), int(len(sub)), float(p))
    p_str = "p < 1e-4" if p < 1e-4 else f"p = {p:.3g}"
    return f"Wilcoxon ipTM vs pTM (paired): {p_str} | median d(ipTM-pTM) = {med_d:+.3f}"


def _umap_tier_separation(df, xcol="UMAP_X", ycol="UMAP_Y", tier_col=CFG.COL_TIER, n_perm=199, cap=2000):
    """Silhouette of tier labels in UMAP space with a label-permutation p-value.
    Subsampled to `cap` points so the O(N^2) silhouette stays tractable."""
    try:
        from sklearn.metrics import silhouette_score
    except Exception:
        return ""
    for c in (xcol, ycol, tier_col):
        if c not in df.columns:
            return ""
    sub = df.dropna(subset=[xcol, ycol, tier_col])
    if len(sub) < 30 or sub[tier_col].nunique() < 2:
        return ""
    rng = np.random.default_rng(0)
    if len(sub) > cap:
        sub = sub.iloc[rng.choice(len(sub), cap, replace=False)]
    X = sub[[xcol, ycol]].to_numpy(float)
    lab = sub[tier_col].astype("category").cat.codes.to_numpy()
    try:
        s0 = silhouette_score(X, lab)
        cnt = sum(1 for _ in range(n_perm) if silhouette_score(X, rng.permutation(lab)) >= s0)
    except Exception:
        return ""
    p = (cnt + 1) / (n_perm + 1)
    p_str = "p < 0.005" if p < 0.005 else f"p = {p:.3g}"
    _register_p("Silhouette label-permutation", "UMAP tier separation", float(s0), int(len(X)), float(p))
    return f"Tier separation in UMAP (silhouette perm-test): {p_str} | s = {s0:.3f}"

# MD-ready cohort — the MD-selected candidate complexes, overlaid as gold stars on tier plots
# (they may sit inside Tier_1A). The 3R3U × fluoroacetate positive control — the one control taken
# to MD — is drawn instead as a RED star with a GOLD border (an inversion of the MD-selected star)
# so the reference lands as a distinct landmark on any figure that already carries MD-selected stars.
# All three styles are single-sourced here + from CFG.
_MD_STAR_KW   = dict(marker="*", s=300, facecolor=CFG.VIS_ACCENT["star"],    edgecolor="black",                    linewidths=1.1, zorder=9)
_CTRL_STAR_KW = dict(marker="*", s=340, facecolor=CFG.VIS_ACCENT["control"], edgecolor=CFG.VIS_ACCENT["control_edge"], linewidths=1.4, zorder=10)

def _is_control_mask(df: pd.DataFrame) -> pd.Series:
    """Boolean mask of control rows (is_control flag from Step 02; safe when the column is absent)."""
    if "is_control" not in df.columns:
        return pd.Series(False, index=df.index)
    return df["is_control"].astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes"])

def _md_ready_df(df: pd.DataFrame) -> pd.DataFrame:
    """MD-selected CANDIDATE cohort (controls excluded) — the complexes taken to MD for the study.
    The 3R3U control is kept out here so it never appears as a generic MD-ready point; it is drawn
    only as its own distinct marker via _control_star_df on the star figures."""
    if "MD_Selected" not in df.columns:
        return df.iloc[0:0]
    _m = df["MD_Selected"].astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes"])
    return df[_m & ~_is_control_mask(df)]

def _control_star_df(df: pd.DataFrame) -> pd.DataFrame:
    """The single 3R3U × FA positive control taken to MD (is_control & Control_Ref == '3R3U' &
    MD_Selected). Empty when the columns are absent (older CSV) or no control was MD-selected."""
    if not {"MD_Selected", "Control_Ref"} <= set(df.columns):
        return df.iloc[0:0]
    _m = df["MD_Selected"].astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes"])
    _r3 = df["Control_Ref"].astype(str).str.strip() == "3R3U"
    return df[_m & _is_control_mask(df) & _r3]

def _ctrl_star(ax, xs, ys, size=None):
    """Overlay the 3R3U × FA positive control as a distinct red star / gold border at (xs, ys)
    (data coords). Used at every figure that draws MD-selected stars so the control never blends
    into the gold cohort. Caller passes the control's coordinates in that figure's own axes."""
    _kw = dict(_CTRL_STAR_KW)
    if size is not None:
        _kw["s"] = size
    ax.scatter(np.atleast_1d(xs), np.atleast_1d(ys), **_kw)


def _md_ready_stars_cat(ax, df: pd.DataFrame, tier_order, valcol: str):
    """Overlay the MD cohort on a categorical (per-tier) axis where x = tier index and y = valcol;
    same-tier stars are x-jittered so they never overlap. Candidates are gold stars; the 3R3U × FA
    control is a red star with a gold border. Returns a LIST of Line2D legend handles (0-2)."""
    _pos = {t: i for i, t in enumerate(tier_order)}

    def _overlay(_sub, _kw):
        _n = 0
        for _t, _grp in _sub.groupby(CFG.COL_TIER):
            if _t not in _pos:
                continue
            _vals = pd.to_numeric(_grp[valcol], errors="coerce").dropna().values
            _xs = np.linspace(-0.18, 0.18, len(_vals)) if len(_vals) > 1 else [0.0]
            for _dx, _v in zip(_xs, _vals):
                ax.scatter([_pos[_t] + _dx], [_v], **_kw)
                _n += 1
        return _n

    _handles = []
    _n_cand = _overlay(_md_ready_df(df), _MD_STAR_KW)
    if _n_cand:
        _handles.append(Line2D([0], [0], marker="*", ls="", markerfacecolor=CFG.VIS_ACCENT["star"],
                               markeredgecolor="black", markersize=12, label=f"MD-ready (n={_n_cand})"))
    _n_ctrl = _overlay(_control_star_df(df), _CTRL_STAR_KW)
    if _n_ctrl:
        _handles.append(Line2D([0], [0], marker="*", ls="", markerfacecolor=CFG.VIS_ACCENT["control"],
                               markeredgecolor=CFG.VIS_ACCENT["control_edge"], markersize=13,
                               label="3R3U × FA (positive control)"))
    return _handles


def _auto_label_colour(bg, threshold: float = 0.5) -> str:
    """The auto-contrast label colour, from utils — one luminance rule for every step's figures."""
    return _utils_mod.auto_label_colour(CFG, bg, threshold)

# Modern Plotting Theme (Publication Quality) — typography and canvas from CFG, via utils.
sns.set_theme(style="whitegrid", context="paper", font_scale=1.3)
_utils_mod.apply_figure_style(CFG)

SEPARATOR = "-" * 80


# =============================================================================
# SECTION 2: LOGGING & UTILITIES
# =============================================================================

logger = None



def console_info(msg: str) -> None:
    _console_info(msg, logger)



def _aux_dir(out_dir: Path) -> Path:
    """Subfolder for the non-figure data deliverables (CSV/TXT), keeping the
    3_Validation_Figures root clutter-free — only figure folders + this one."""
    d = out_dir / "01_Analysis_Data"
    d.mkdir(parents=True, exist_ok=True)
    return d


ReportManager = _utils_mod.ReportManager   # shared logger (00_02)


def _make_reporter(out_dir: Path):
    """Construct the shared ReportManager with this step's log path/header/logger."""
    return ReportManager(
        _aux_dir(out_dir) / "00_Validation_Figures.log",
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
    Load the ranked CSV from Step 02 (has Scientific_Rank + the tier/MD_Selected columns the
    figures key on). Ranked-only: no master-CSV fallback — a missing ranked CSV is a hard error,
    since every figure axis and the enriched dataset depend on Scientific_Rank.
    """
    # SSOT glob first, wildcard fallback; newest by mtime (a name sort can rank an older file last when the leading number differs).
    target = (latest_by_mtime(prod_dir.glob(CFG.GLOB_RANKED_CSV))
              or latest_by_mtime(prod_dir.glob("*_Ranked_*.csv")))
    if target is None: raise FileNotFoundError("No Ranked CSV found in Production folder.")
    reporter.log(f"Input Data Source: {target.resolve()}")
    df = pd.read_csv(target, low_memory=False)

    if "job_name" in df.columns:
        df = df.drop_duplicates(subset=["job_name"], keep="last")

    tier_map = CFG.TIER_RANK
    df["tier_numeric"] = df[CFG.COL_TIER].map(tier_map).fillna(0)

    col_map = {
        "binding_likelihood_computed": "Binding_Probability",
        "Binding_Probability_Score": "Binding_Probability",
        "custom_affinity_score": "Chemical_Affinity_Score",
        "confidence_score": CFG.COL_CONF,
        CFG.COL_IDENS: "Interaction_Density_Norm",
        "sn2_attack_angle": CFG.COL_SN2,
        "Active_Site_RMSD_to_Control": "Active_Site_RMSD",
    }
    df.rename(columns=col_map, inplace=True)

    # tier_numeric is the target label — it is deliberately NOT an input feature,
    # else PCA/UMAP would align their principal axes with the tier and manufacture
    # a circular, guaranteed tier separation (target leakage). It is used only to
    # orient the PC1 sign after projection and to colour points post hoc.
    candidate_features = [
        CFG.COL_CONF, "iptm", "mean_plddt",
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


def write_residue_mapping_missed_csv(out_dir: Path, reporter) -> None:
    """Step 1 (01_Analysis_Data): per active-site residue, list the query proteins whose alignment
    failed to map it → 01_Active_Site_Residue_Mapping_Missed.csv. Written here so the whole
    01_Analysis_Data set is finalised before the figure steps; the Step-3 coverage figure recomputes
    its own bar counts from the same alignment table rather than this CSV."""
    import json as _json
    _aln_csv = (out_dir.parent / "1_Boltz2_Production" / "3_Sequence_Reference_Data"
                / "Active_Site_Alignments" / CFG.FILE_ALIGNMENT_STATS)
    if not _aln_csv.exists():
        return
    _aa3 = {"ASP": "Asp", "ARG": "Arg", "HIS": "His", "TRP": "Trp", "TYR": "Tyr"}
    _ref = CFG.REF_ACTIVE_SITE_MAP
    _keys = [k for k in CFG.ACTIVE_SITE_ROLE_ORDER if k in _ref]
    _adf = pd.read_csv(_aln_csv)
    if "protein" in _adf.columns:   # controls (…_Control) are excluded — coverage reflects the query set
        _adf = _adf[~_adf["protein"].astype(str).str.endswith("_Control")].drop_duplicates("protein")
    _missed = {k: [] for k in _keys}
    for _prot, _s in zip(_adf.get("protein", range(len(_adf))), _adf["active_site_mapping"].fillna("")):
        try:
            _m = _json.loads(_s) if _s else None
        except Exception:
            _m = None
        if not isinstance(_m, dict):
            # An absent or unparseable mapping produced NO residues — this is the coverage-evidence
            # CSV, so count the protein as missed for every role. Skipping the row would let an
            # alignment that mapped nothing read as fully mapped (missing is not zero-missed).
            for k in _keys:
                _missed[k].append(str(_prot))
            continue
        for k in _keys:
            if _m.get(k) is None:
                _missed[k].append(str(_prot))
    _col = {k: f"{_aa3.get(_ref[k]['res'], _ref[k]['res'].title())}{_ref[k]['id']}_{k}" for k in _keys}
    _maxlen = max((len(v) for v in _missed.values()), default=0)
    pd.DataFrame({_col[k]: _missed[k] + [""] * (_maxlen - len(_missed[k])) for k in _keys}).to_csv(
        _aux_dir(out_dir) / "01_Active_Site_Residue_Mapping_Missed.csv", index=False)
    reporter.log("  ✔ Saved: 01_Analysis_Data/01_Active_Site_Residue_Mapping_Missed.csv")


def perform_advanced_ranking(df: pd.DataFrame, features: list[str], out_dir: Path, reporter: ReportManager):
    """Multi-objective ranking of complexes.

    Scales the feature matrix (RobustScaler — it centres a zero-IQR column instead of
    dividing by zero), derives a PCA PC1 score (sign-aligned to tier_numeric), a CFG-weighted
    composite score, Pareto fronts (confidence vs binding probability) and a UMAP
    embedding. Adds the Score/Ensemble/Pareto/UMAP columns to df in place and
    returns it.
    """
    reporter.section("Step 1/8 — 01_Analysis_Data")
    reporter.log("  Multi-Objective Ranking (PCA & Pareto)")

    x = df[features].dropna()
    '''
    RobustScaler already guards its own zero-IQR case: sklearn's _handle_zeros_in_scale sets a zero
    scale to 1.0 per column, so a constant / heavily-skewed feature is merely centred, never divided by
    zero — no inf/NaN reaches PCA or UMAP. A whole-matrix StandardScaler would re-scale EVERY column the
    moment one degenerate column appeared, throwing away the outlier robustness the ranking depends on, so
    RobustScaler is used unconditionally.
    '''
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
        check_col = next((c for c in (CFG.COL_MECH_S, "competence_score") if c in x.columns), features[0])
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
        df.loc[x.index, "Score_PCA_Data"] = CFG.SCORE_PCA_NEUTRAL

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
    if CFG.COL_CONF in df.columns and "Binding_Probability" in df.columns:
        _mech_obj = next((c for c in ("mechanistic_score_effective", CFG.COL_MECH_S) if c in df.columns), None)
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
            _pr = calculate_pareto_fronts(df.loc[_viable], [CFG.COL_CONF, "Binding_Probability"], [True, True])
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
            reducer = umap.UMAP(n_neighbors=min(15, len(x)-1), min_dist=0.1, random_state=int(CFG.ANALYSIS_SEED))
            umap_map = reducer.fit_transform(x_scaled)
            df.loc[x.index, "UMAP_X"] = umap_map[:, 0]
            df.loc[x.index, "UMAP_Y"] = umap_map[:, 1]
        except Exception as e:
            reporter.log(f"  ! UMAP Calculation Skipped: {e}")
    else:
        reporter.log("  ! Skipping UMAP: Insufficient data points (< 5).")

    if len(x) >= 2:
        pd.DataFrame(pca.components_.T, columns=["PC1", "PC2"], index=features).to_csv(_aux_dir(out_dir) / "02_PCA_Loadings.csv", index=True)
        reporter.log("  ✔ Saved: 01_Analysis_Data/02_PCA_Loadings.csv")

    return df

def analyse_conflicts(df: pd.DataFrame, out_dir: Path, reporter: ReportManager):
    """Flag conflict/opportunity cases and write the validated master CSV.

    Classifies each complex (e.g. hidden gems = high mechanistic merit but
    low confidence; decoys = high confidence but poor mechanism), logs the
    per-class counts, and saves the final validated master table to out_dir.
    """
    reporter.log("")
    reporter.log("  Conflict & Opportunity Analysis")

    """
    The classification is four mutually exclusive tests on two columns, so it is expressed as vector
    masks rather than a Python callback invoked once per row: np.select evaluates the same conditions
    in the same order (first match wins, as in the original if-chain) across the whole frame at once.
    On 58,000 rows the row-wise apply was the single slowest statement in this step.
    """
    _high_quality = [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]
    _hi, _lo = CFG.CONFLICT_CONF_HIGH, CFG.CONFLICT_CONF_LOW

    _tier = df[CFG.COL_TIER] if CFG.COL_TIER in df.columns else pd.Series(CFG.TIER_DECOY, index=df.index)
    _conf = (pd.to_numeric(df[CFG.COL_CONF], errors="coerce").fillna(0.0)
             if CFG.COL_CONF in df.columns else pd.Series(0.0, index=df.index))

    _is_hq = _tier.isin(_high_quality)
    _is_lowtier = _tier.isin([CFG.TIER_POOR, CFG.TIER_DECOY])

    df["Conflict_Category"] = np.select(
        [
            _is_hq & (_conf >= _hi),
            _tier.isin([CFG.TIER_POOR, CFG.TIER_DECOY, "Error"]) & (_conf < _lo),
            _is_hq & (_conf < _hi),
            _is_lowtier & (_conf >= _hi),
        ],
        ["Consensus High", "Consensus Low", "Hidden Gem", "Decoy"],
        default="Ambiguous",
    )

    gems = df[df["Conflict_Category"] == "Hidden Gem"].sort_values("Pareto_Rank")
    gems.to_csv(_aux_dir(out_dir) / "04_ACTION_Rescue_Hidden_Gems.csv", index=False)
    reporter.log("  ✔ Saved: 01_Analysis_Data/04_ACTION_Rescue_Hidden_Gems.csv")

    reporter.log(f"Hidden Gems (Rescue Target) : {len(gems)}")
    reporter.log(f"Decoys (Potential Artifacts): {len(df[df['Conflict_Category'] == 'Decoy'])}")

    return df

# =============================================================================
# SECTION 4B: Tier_1A Landscape Companion Figures (the *b variants: 05b, 11b, 13b,
#             14b, 17b, 19b, 20b, 24b, 26b)
# =============================================================================

"""
Companion figures overlaying Tier_1A structural thumbnails onto the same
chemical spaces as their parent figures, showing where the Tier_1A hits
sit relative to the full 58 k-complex dataset.
"""

# ── Design constants ──────────────────────────────────────────────────────────
_TT_ENTRY_COLS = [CFG.VIS_ACCENT_DEEP["brick"], CFG.VIS_ACCENT_DEEP["jade"], CFG.VIS_ACCENT_DEEP["emerald_alt"], CFG.VIS_RAMP["purple"][2], CFG.VIS_ACCENT_DEEP["tangerine"]]
_TT_STAR_FILL  = CFG.VIS_ACCENT_DEEP["emerald"]
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
    cifs = sorted(bc.glob("*.cif")) if bc.exists() else sorted(jd.rglob("*.cif"))
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
        sub = df[df[CFG.COL_TIER] == tier]
        if sub.empty: continue
        ax.scatter(sub["_X"], sub["_Y"],
                   c=TIER_PALETTE.get(tier, CFG.VIS_INK["paler"]),
                   s=_TT_TIER_S.get(tier, 4), alpha=_TT_TIER_A.get(tier, 0.2),
                   linewidths=0, rasterized=True, zorder=2)


def _tt_draw_stars(ax, pa):
    for local_i, (_, row) in enumerate(pa.iterrows()):
        # The 3R3U × FA positive control (present in _pa via MD_Selected) is drawn as a distinct
        # red star with a gold border + gold halo, so it never blends into the per-entry palette.
        if str(row.get("is_control", "")).strip().lower() in ("true", "1", "1.0", "yes"):
            ax.scatter(row["_X"], row["_Y"], c=CFG.VIS_ACCENT["control"], s=int(_TT_STAR_S * 1.15),
                       marker="*", edgecolors=CFG.VIS_ACCENT["control_edge"], linewidths=2.2, zorder=12)
            continue
        ec   = _TT_ENTRY_COLS[local_i % len(_TT_ENTRY_COLS)]
        fill = CFG.VIS_ACCENT_DEEP["sun"] if local_i == 4 else _TT_STAR_FILL
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
        # The 3R3U × FA positive control gets a gold border (CFG.VIS_ACCENT["control_edge"]) —
        # deliberately outside the per-entry palette so its thumbnail reads instantly as the control.
        if str(row.get("is_control", "")).strip().lower() in ("true", "1", "1.0", "yes"):
            colour = CFG.VIS_ACCENT["control_edge"]
        else:
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
        prot = str(row.get(CFG.COL_PROT, "")).replace("_Control","").split("_")[0][:11]
        lig  = str(row.get(CFG.COL_LIG, ""))
        lig  = lig[:13] + "…" if len(lig) > 13 else lig
        sn2  = row.get(CFG.COL_SN2,        float("nan"))
        conf = row.get(CFG.COL_CONF,   float("nan"))
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
    n_tt = len(df[df[CFG.COL_TIER] == CFG.TIER_TOP])
    # Symbol-type header entries so the reader understands both glyphs
    h = [
        Line2D([0],[0], marker="o", color="w",
               markerfacecolor=CFG.VIS_INK["ghost"], markeredgecolor="none",
               markersize=7, alpha=0.55, label="Individual complex  (scatter dot)"),
        Line2D([0],[0], marker="*", color="w",
               markerfacecolor=_TT_STAR_FILL, markeredgecolor=_TT_ENTRY_COLS[0],
               markersize=14, label=f"{CFG.TIER_TOP} ★ highlighted  (n={n_tt})"),
        Line2D([0],[0], marker="*", color="w",
               markerfacecolor=CFG.VIS_ACCENT_DEEP["sun"], markeredgecolor=_TT_ENTRY_COLS[4],
               markersize=11, label=f"{CFG.TIER_TOP} ★ gold (TFA duplicate)"),
    ]
    # 3R3U × FA positive control — distinct red star / gold border (only when present in the frame)
    if not _control_star_df(df).empty:
        h.append(Line2D([0],[0], marker="*", color="w",
                        markerfacecolor=CFG.VIS_ACCENT["control"], markeredgecolor=CFG.VIS_ACCENT["control_edge"],
                        markersize=15, markeredgewidth=1.8, label="3R3U × FA  (positive control)"))
    for t in [t for t in CFG.TIER_ORDER if t != CFG.TIER_TOP]:
        sub = df[df[CFG.COL_TIER] == t]
        if sub.empty:
            continue
        h.append(Line2D([0],[0], marker="o", color="w",
                        markerfacecolor=TIER_PALETTE.get(t, CFG.VIS_INK["paler"]),
                        markersize=7, alpha=0.85, label=f"{t}  (n={len(sub):,})"))
    return h


def _tt_add_legend(fig, handles):
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.045), ncol=4,   # small gap below the x-axis label — not over it, not far
               title=f"Degrader tier  (★ = {CFG.TIER_TOP} highlighted; ● = scatter background)",
               )


def _tt_new_fig():
    fig, ax = plt.subplots(figsize=(18, 12))
    fig.subplots_adjust(**_TT_MARGINS)
    ax.set_facecolor(CFG.VIS_INK["panel"])
    return fig, ax


def _tt_style(ax, xlabel, ylabel, title=None):
    # `title` accepted for call-site compatibility but intentionally not rendered
    # (figures carry no titles).
    ax.set_xlabel(xlabel,  labelpad=6)
    ax.set_ylabel(ylabel,  labelpad=6)
    ax.tick_params(labelsize=9)
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.18, color=CFG.VIS_INK["pale"], linewidth=0.5)


# ── Figure 13b — SN2 Angle × Confidence + PA Landscape ───────────────────────
def _fig_13b_tt_mechanistic(df, pa, imgs, out_dir: Path, reporter):
    try:
        dv = df.dropna(subset=[CFG.COL_SN2,CFG.COL_CONF]).copy()
        dv["_X"] = pd.to_numeric(dv[CFG.COL_SN2],      errors="coerce")
        dv["_Y"] = pd.to_numeric(dv[CFG.COL_CONF], errors="coerce")
        dv = dv.dropna(subset=["_X","_Y"])
        if dv.empty:
            reporter.log(f"  ! Skipped: {_fig_path('13b')} — no rows with both axes present")
            return
        pax = pa.copy()
        pax["_X"] = pd.to_numeric(pa[CFG.COL_SN2],      errors="coerce")
        pax["_Y"] = pd.to_numeric(pa[CFG.COL_CONF], errors="coerce")
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
            ax.text(angle, yhi - 0.001, f"≥{angle:g}° {tlbl.replace('Tier_', '')}", rotation=90,
                    rotation_mode="anchor",
                    fontsize=7.0, color=tcol, fontweight="bold", va="center", ha="right",
                    zorder=7,
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", alpha=0.85, ec=tcol,
                              linewidth=0.5))
        ax.axhline(dv["_Y"].median(), color=CFG.VIS_INK["muted"], ls=":", lw=0.9, alpha=0.55, zorder=3)
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
        reporter.log(f"  ! Skipped: {_fig_path('13b')} — {e}")
        plt.close("all")   # release the figure left open by the failed savefig


# ── Figure 05b — Confidence × ipTM + PA Landscape ────────────────────────────
def _fig_05b_tt_ai_quality(df, pa, imgs, out_dir: Path, reporter):
    try:
        dv = df.dropna(subset=[CFG.COL_CONF,"iptm"]).copy()
        if dv.empty:
            reporter.log(f"  ! Skipped: {_fig_path('05b')} — no rows with both confidence and ipTM")
            return
        dv["_X"] = dv[CFG.COL_CONF]; dv["_Y"] = dv["iptm"]
        pax = pa.copy()
        pax["_X"] = pa[CFG.COL_CONF]; pax["_Y"] = pa["iptm"]
        xlo,xhi = dv["_X"].quantile(0.01)-0.01, dv["_X"].max()+0.005
        ylo,yhi = dv["_Y"].quantile(0.01)-0.01, dv["_Y"].max()+0.005
        fig, ax = _tt_new_fig()
        hb = ax.hexbin(dv["_X"], dv["_Y"], gridsize=60, cmap="YlOrBr",
                       mincnt=1, linewidths=0.15, alpha=0.85, zorder=1,
                       extent=[xlo,xhi,ylo,yhi])
        cb = fig.colorbar(hb, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("Complex count per bin", fontsize=9.5)
        for v, fn in [(dv["_X"].median(), ax.axvline),(dv["_Y"].median(), ax.axhline)]:
            fn(v, color=CFG.VIS_INK["muted"], ls="--", lw=0.9, alpha=0.55, zorder=3)
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
        reporter.log(f"  ! Skipped: {_fig_path('05b')} — {e}")
        plt.close("all")   # release the figure left open by the failed savefig


# ── Figure 14b — Interaction Density × Count + PA Landscape ──────────────────
def _fig_14b_tt_interactions(df, pa, imgs, out_dir: Path, reporter):
    try:
        dv = df.dropna(subset=["Interaction_Density_Norm","num_interactions"]).copy()
        dv["_X"] = pd.to_numeric(dv["Interaction_Density_Norm"], errors="coerce")
        dv["_Y"] = pd.to_numeric(dv["num_interactions"],          errors="coerce")
        dv = dv.dropna(subset=["_X","_Y"])
        if dv.empty:
            reporter.log(f"  ! Skipped: {_fig_path('14b')} — no rows with both interaction metrics")
            return
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
        ax.axvline(xmed, color=CFG.VIS_INK["muted"], ls="--", lw=0.9, alpha=0.55, zorder=3)
        ax.axhline(ymed, color=CFG.VIS_INK["muted"], ls="--", lw=0.9, alpha=0.55, zorder=3)
        ax.text(xmed+0.05, 0.5, f"median {xmed:.2f}", fontsize=_JF_FA, color=CFG.VIS_INK["muted"])
        ax.text(0.1, ymed+0.3, f"median {ymed:.0f}",  fontsize=_JF_FA, color=CFG.VIS_INK["muted"])
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
        reporter.log(f"  ! Skipped: {_fig_path('14b')} — {e}")
        plt.close("all")   # release the figure left open by the failed savefig


_FIG_MAPPING = {
        # ── 03_Dataset_and_Alignment_Overview ──
        "Figure_01_Active_Site_Residue_Mapping_Coverage.png": "03_Dataset_and_Alignment_Overview/01_Active_Site_Residue_Mapping_Coverage.png",
        "Figure_02_Tier_Distribution.png": "03_Dataset_and_Alignment_Overview/02_Tier_Distribution.png",
        "Figure_03_Alignment_Grades.png": "03_Dataset_and_Alignment_Overview/03_Alignment_Grades.png",
        "Figure_04_Tier_Grade_Distribution.png": "03_Dataset_and_Alignment_Overview/04_Tier_Grade_Distribution.png",
        # ── 04_AI_Confidence_Quality ──
        "Figure_05a_AI_Quality_Assessment.png": "04_AI_Confidence_Quality/01_AI_Quality_Assessment.png",
        "Figure_05b_TT_AI_Quality_Space.png": "04_AI_Confidence_Quality/02_AI_Quality_Space.png",
        "Figure_06_pTM_vs_ipTM_by_Tier.png": "04_AI_Confidence_Quality/03_pTM_vs_ipTM_by_Tier.png",
        # ── 05_Catalytic_Geometry_and_Mechanism ──
        "Figure_07_ActiveSite_RMSD_by_Tier.png": "05_Catalytic_Geometry_and_Mechanism/01_ActiveSite_RMSD_by_Tier.png",
        "Figure_11b_Mechanistic_Fingerprint.png": "05_Catalytic_Geometry_and_Mechanism/06_Mechanistic_Fingerprint.png",
        "Figure_08_Feature_Correlations.png": "05_Catalytic_Geometry_and_Mechanism/02_Feature_Correlations.png",
        "Figure_09_Tier_Quality_DotPlot.png": "05_Catalytic_Geometry_and_Mechanism/03_Tier_Quality_DotPlot.png",
        "Figure_10_Mech_State_CrossTab.png": "05_Catalytic_Geometry_and_Mechanism/04_Mech_State_CrossTab.png",
        "Figure_11_Mechanistic_Score_by_Tier.png": "05_Catalytic_Geometry_and_Mechanism/05_Mechanistic_Score_by_Tier.png",
        "Figure_12_SN2_Angle_by_Tier.png": "05_Catalytic_Geometry_and_Mechanism/07_SN2_Angle_by_Tier.png",
        "Figure_13a_Mechanism_Geometry_Scatter.png": "05_Catalytic_Geometry_and_Mechanism/08_Mechanism_Geometry_Scatter.png",
        "Figure_13b_TT_Mechanistic_Quality_Space.png": "05_Catalytic_Geometry_and_Mechanism/09_Mechanistic_Quality_Space.png",
        # ── 06_Ligand_Interactions_and_Chemical_Space ──
        "Figure_14a_Molecular_Interaction_Profile.png": "06_Ligand_Interactions_and_Chemical_Space/01_Molecular_Interaction_Profile.png",
        "Figure_17b_Binding_Energetics.png": "06_Ligand_Interactions_and_Chemical_Space/06_Binding_Energetics.png",
        "Figure_14b_TT_Interaction_Quality_Space.png": "06_Ligand_Interactions_and_Chemical_Space/02_Interaction_Quality_Space.png",
        "Figure_15_Fluorine_Engagement_by_Tier.png": "06_Ligand_Interactions_and_Chemical_Space/03_Fluorine_Engagement_by_Tier.png",
        "Figure_16_Catalytic_Quality_vs_Inhibition.png": "06_Ligand_Interactions_and_Chemical_Space/04_Catalytic_Quality_vs_Inhibition.png",
        "Figure_17_ActiveSite_Contact_Density_by_Tier.png": "06_Ligand_Interactions_and_Chemical_Space/05_ActiveSite_Contact_Density_by_Tier.png",
        "Figure_18a_Chemical_Space_Map.png": "06_Ligand_Interactions_and_Chemical_Space/07_Chemical_Space_Map.png",
        # ── 07_PFAS_Scope_and_Synthesis ──
        "Figure_19a_Radar_TopHits.png": "07_PFAS_Scope_and_Synthesis/01_Radar_TopHits.png",
        "Figure_19b_Radar_TierReps.png": "07_PFAS_Scope_and_Synthesis/02_Radar_TierReps.png",
        "Figure_20a_Tier_Success_Rates.png": "07_PFAS_Scope_and_Synthesis/03_Tier_Success_Rates.png",
        "Figure_20b_Conf_SN2_Landscape.png": "07_PFAS_Scope_and_Synthesis/04_Conf_SN2_Landscape.png",
        "Figure_21_Conflict_Composition.png": "07_PFAS_Scope_and_Synthesis/05_Conflict_Composition.png",
        "Figure_22_Hidden_Gems_DeepDive.png": "07_PFAS_Scope_and_Synthesis/06_Hidden_Gems_DeepDive.png",
        "Figure_23_Category_Overlap_Euler.png": "07_PFAS_Scope_and_Synthesis/07_Category_Overlap_Euler.png",
        "Figure_24a_Top25_Multitarget_Proteins.png": "07_PFAS_Scope_and_Synthesis/08_Top25_Multitarget_Proteins.png",
        "Figure_24b_TopTier_Protein_PFAS_Breakdown.png": "07_PFAS_Scope_and_Synthesis/09_TopTier_Protein_PFAS_Breakdown.png",
        "Figure_25_Sankey_Workflow.png": "07_PFAS_Scope_and_Synthesis/10_Sankey_Workflow.png",
        "Figure_26a_PFAS_Size_Hexbin_Landscape.png": "07_PFAS_Scope_and_Synthesis/11_PFAS_Size_Hexbin_Landscape.png",
        "Figure_26b_PFAS_Size_Composition_Merged.png": "07_PFAS_Scope_and_Synthesis/12_PFAS_Size_Composition_Merged.png",
        "Figure_26c_PFAS_Carbon_Confidence.png": "07_PFAS_Scope_and_Synthesis/13_PFAS_Carbon_Confidence.png",
        # ── Two-criteria tier logic + dead-end feasibility (single 3-panel figure,
        #    filed in the catalytic folder) ──
        "Figure_30a_Criterion_A_Gates_B.png": "05_Catalytic_Geometry_and_Mechanism/10_Criterion_A_Gates_B.png",
        "Figure_30b_Criterion_B_ECDF_by_Tier.png": "05_Catalytic_Geometry_and_Mechanism/11_Criterion_B_ECDF_by_Tier.png",
        "Figure_30c_SN2_DeadEnd_Gate.png": "05_Catalytic_Geometry_and_Mechanism/12_SN2_DeadEnd_Gate.png",
}


# Folded per-folder analysis panels are routed by _panel (not _FIG_MAPPING); their skip logs map here.
_PANEL_FIG_PATHS = {
    "01C": "05_Catalytic_Geometry_and_Mechanism/13_Geometry_and_Uncertainty.png",
    "05b": "05_Catalytic_Geometry_and_Mechanism/14_Mechanistic_Breakdown_by_Tier.png",
    "02A": "06_Ligand_Interactions_and_Chemical_Space/08_Binding_Affinity_Metrics.png",
    "04A": "03_Dataset_and_Alignment_Overview/05_Evolutionary_Phylogeny.png",
    "05c": "07_PFAS_Scope_and_Synthesis/14_Chain_Length_by_Tier.png",
    "06A": "07_PFAS_Scope_and_Synthesis/15_Tier1A_Cross_Ligand_Heatmap.png",
    "05a": "08_Diagnostic_and_MultiModel_Trends/09_Pillar_Divergence_by_Tier.png",
}


def _fig_path(fig_key: str) -> str:
    """Map a figure key (e.g. '22', '13b', '01C') to its folder/NN_name.png for skip logs — the folded
    panels via _PANEL_FIG_PATHS, the main suite via _FIG_MAPPING. Falls back to 'Figure <key>'."""
    if fig_key in _PANEL_FIG_PATHS:
        return _PANEL_FIG_PATHS[fig_key]
    m = next((v for k, v in _FIG_MAPPING.items() if k.startswith(f"Figure_{fig_key}_")), None)
    return m if m else f"Figure {fig_key}"


@contextlib.contextmanager
def _redirect_savefig(out_dir: Path, reporter):
    """Route each Figure_NN_*.png save to its numbered folder (per _FIG_MAPPING)
    for the duration of the block, restoring the original savefig on exit — even
    on exception, so a failure mid-suite can never leak the patched savefig into
    later routines or the rest of the process. Target folders are created on entry.
    """
    import matplotlib
    _orig_plt = plt.savefig
    _orig_fig = matplotlib.figure.Figure.savefig
    def _redirect(fname):
        name = fname.name if isinstance(fname, Path) else Path(fname).name
        if name in _FIG_MAPPING:
            _target = out_dir / _FIG_MAPPING[name]
            _target.parent.mkdir(parents=True, exist_ok=True)   # folder created lazily, at its turn
            reporter.log(f'  ✔ Saved: {_FIG_MAPPING[name]}')
            return _target
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



# ---------------------------------------------------------------------------
# Shared helper hoisted from the figure suite (auto-contrast label colour).
# ---------------------------------------------------------------------------
def _text_color(hex_bg: str, threshold: float = 0.5) -> str:
    """Auto-contrast label colour for hex_bg (delegates to _auto_label_colour)."""
    return _auto_label_colour(hex_bg, threshold)


# =============================================================================
# ADDITIONAL JOURNAL FIGURE TYPES (folded panels, routed through _panel)
# Seven figure types 03 does not otherwise produce — volcano, Manhattan, swimmer,
# feature-importance + tier streams, circos, metric network, treemap.
# Each is a folded panel: fn(df, folder_dir, reporter) that saves folder_dir/NN_name.png,
# so _panel logs it '  ✔ Saved: folder/NN_name.png' in the owning step's block, exactly like
# every other figure. Styling stays inside a private rc-context so 03's global figure style is
# untouched; names are prefixed _jf_/_JF_ to avoid any collision with the rest of the module.
# =============================================================================
from matplotlib.patches import Rectangle as _JF_Rect, PathPatch as _JF_PathPatch, Patch as _JF_Patch
from matplotlib.lines import Line2D as _JF_Line
from matplotlib.path import Path as _JF_MPath

# All colour, font and DPI values below are SOURCED FROM CFG (single source of truth) — the same
# palette, fonts and quality every other 03 figure uses. No titles (none is set).
_JF_TORDER = ["Tier_1A", "Tier_1B", "Tier_2A", "Tier_2B", "Tier_3", "Tier_4", "Tier_5_Decoy"]
_JF_TCOL = dict(CFG.TIER_COLOUR)                                    # tier palette from CFG
_JF_OKABE = [CFG.VIS_ACCENT[k] for k in ("blue", "vermillion", "green", "amber",
                                         "magenta", "sky", "yellow")] + [CFG.VIS_INK["faint"]]
_JF_SIG = {"up": CFG.VIS_ACCENT["bad"], "down": CFG.VIS_ACCENT["blue"], "ns": CFG.VIS_INK["palest"]}
# rc-context that re-asserts CFG's own figure style (fonts/family/sizes/grid) for these panels, so a
# stray global change cannot alter them; identical settings to _utils_mod.apply_figure_style(CFG).
_JF_RC = {"figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
          "font.family": "sans-serif", "font.sans-serif": list(CFG.VIS_FONT_FAMILY),
          "axes.labelsize": CFG.VIS_FONT_AXIS_LABEL, "xtick.labelsize": CFG.VIS_FONT_TICK,
          "ytick.labelsize": CFG.VIS_FONT_TICK, "legend.fontsize": CFG.VIS_FONT_LEGEND,
          "axes.edgecolor": "black", "axes.linewidth": 0.5, "axes.grid": True, "axes.axisbelow": True,
          "grid.color": CFG.VIS_INK["grid"], "grid.linewidth": 0.3, "legend.frameon": False}
_JF_FA = float(CFG.VIS_FONT_ANNOT)          # in-figure annotations (from CFG)
_JF_FD = float(CFG.VIS_FONT_ANNOT) - 1.0    # dense categorical tick labels (many ligands/metrics)
_JF_LIGSHORT = {"Fluoroacetate": "FA", "Difluoroacetate": "DFA"}
# Canonical short PFAS display label — the ONE shortener every figure routes through so ligand
# names read uniformly: strip the "NN_" ordering prefix and collapse the long acid names to
# FA / DFA / TFA. Display-only; the underlying data columns keep their full values.
_LIG_SHORT_MAP = {"fluoroacetate": "FA", "monofluoroacetate": "FA", "mono-fluoroacetate": "FA", "mfa": "FA",
                  "difluoroacetate": "DFA", "dfa": "DFA",
                  "trifluoroacetate": "TFA", "tfa": "TFA"}


def _lig_short(name) -> str:
    if name is None or (isinstance(name, float) and pd.isna(name)):
        return ""
    s = re.sub(r"^\d+_", "", str(name).strip())
    return _LIG_SHORT_MAP.get(s.lower().replace("_ref", "").strip(), s)


def _lig_short_series(s):
    return s.astype(str).map(_lig_short)


_JF_MLAB = {"iptm": "iPTM", "ptm": "pTM", "mean_plddt": "pLDDT", "Boltz_Model_Confidence": "Confidence",
            "Binding_Probability_Score": "Binding P", "Binding_Probability": "Binding P",
            "custom_affinity_score": "Affinity", "Chemical_Affinity_Score": "Affinity",
            "interaction_density": "Int. density", "Interaction_Density_Norm": "Int. density",
            "competence_score": "Competence", "feasibility_factor": "Feasibility",
            "mechanistic_score_effective": "Mech. score", "soft_catalytic_score": "Soft cat.",
            "catalytic_constellation_score": "Cat. constel.", "halide_stabilisation_score": "Halide stab.",
            "carboxylate_clamp_integrity": "Cbx clamp", "pocket_occupancy": "Pocket occ.",
            "num_interactions": "N interact.", "count_hydrogen_bond": "H-bonds",
            "count_salt_bridge": "Salt bridges", "SN2_Attack_Angle": "S$_N$2 angle",
            "scissile_cf_bde": "C–F BDE", "sn2_backside_occlusion": "Backside occ."}
_JF_METRIC_CANDS = ["iptm", "ptm", "mean_plddt", "Boltz_Model_Confidence", "Binding_Probability_Score",
                    "Binding_Probability", "custom_affinity_score", "Chemical_Affinity_Score",
                    "interaction_density", "Interaction_Density_Norm", "competence_score",
                    "feasibility_factor", "mechanistic_score_effective", "soft_catalytic_score",
                    "catalytic_constellation_score", "halide_stabilisation_score",
                    "carboxylate_clamp_integrity", "pocket_occupancy", "num_interactions",
                    "count_hydrogen_bond", "count_salt_bridge", "SN2_Attack_Angle", "scissile_cf_bde",
                    "sn2_backside_occlusion"]


def _jf_ml(c):
    return _JF_MLAB.get(c, str(c).replace("_", " "))


def _jf_ligshort(df):
    return _lig_short_series(df["Ligand_Name"])


def _jf_tiers(df, ts):
    return [t for t in _JF_TORDER if (ts == t).any()]


def _jf_metrics(df):
    seen = set(); out = []
    for c in _JF_METRIC_CANDS:
        lab = _jf_ml(c)
        if c in df.columns and lab not in seen and pd.to_numeric(df[c], errors="coerce").notna().sum() > 100:
            out.append(c); seen.add(lab)
    return out


def _jf_metric_matrix(df):
    cols = _jf_metrics(df)
    return df[cols].apply(pd.to_numeric, errors="coerce")


def _jf_box(ax):
    ax.grid(which="minor", visible=False)
    for s in ax.spines.values():
        s.set_visible(True); s.set_color("black"); s.set_linewidth(0.5)


def _jf_save(fig, folder_dir, name):
    Path(folder_dir).mkdir(parents=True, exist_ok=True)
    fig.savefig(Path(folder_dir) / name, dpi=int(CFG.VIS_FIGURE_DPI), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _jf_fit_fs(ax, cx, cy, dxd, dyd, lines, base, fmin):
    """Largest font (≤ base pt) at which `lines` fit inside the box (data coords cx±dxd/2, cy±dyd/2);
    returns None when that would fall below fmin — the caller then skips the label rather than render
    it unreadably small. Sizing uses the box's on-screen extent so it adapts to each rectangle."""
    p0 = ax.transData.transform((cx - dxd / 2.0, cy - dyd / 2.0))
    p1 = ax.transData.transform((cx + dxd / 2.0, cy + dyd / 2.0))
    bw = abs(p1[0] - p0[0]); bh = abs(p1[1] - p0[1])        # box size, display px
    if bw <= 1 or bh <= 1:
        return None
    ppp = ax.figure.dpi / 72.0                              # px per point
    maxchars = max((len(str(s)) for s in lines), default=1) or 1
    fs = min(base, (bw * 0.86) / (maxchars * 0.58 * ppp), (bh * 0.80) / (len(lines) * 1.30 * ppp))
    return fs if fs >= fmin else None


def _jf_volcano(df, folder_dir, reporter):
    """Multi-group volcano — per-metric effect size (degrader − non-degrader) vs Mann–Whitney significance."""
    from scipy import stats as _sps
    with plt.rc_context(_JF_RC):
        fig = plt.figure(figsize=(9, 5.5)); ax = fig.subplots()
        deg = df["is_degrader"] == True; rows = []
        for c in _jf_metrics(df):
            a = pd.to_numeric(df.loc[deg, c], errors="coerce").dropna()
            b = pd.to_numeric(df.loc[~deg, c], errors="coerce").dropna()
            if len(a) < 30 or len(b) < 30:
                continue
            _u, p = _sps.mannwhitneyu(a, b, alternative="two-sided")
            pooled = np.sqrt(((a.std()**2)+(b.std()**2))/2) or 1e-9
            _d = (a.mean()-b.mean())/pooled
            # Register into the shared statistical family so these tests are BH-corrected together with
            # the battery and land in 06_Statistical_Tests.csv linked to this figure. Distinct from the
            # battery's Mann-Whitney: that groups by tier membership, this by the is_degrader flag.
            _register_p(f"Mann-Whitney U — is_degrader — {_jf_ml(c)}",
                        "10_Volcano_Metric_Significance", float(_u), int(len(a) + len(b)), float(p),
                        effect_size=round(float(_d), 4), effect_type="Cohen d",
                        n_group_a=int(len(a)), n_group_b=int(len(b)),
                        median_a=round(float(np.median(a)), 4), median_b=round(float(np.median(b)), 4))
            rows.append((c, _d, p))
        v = pd.DataFrame(rows, columns=["metric", "effect", "p"]); v["nlp"] = -np.log10(v["p"].clip(lower=1e-300))
        v["col"] = np.where((v.p < 0.05) & (v.effect > 0.1), _JF_SIG["up"],
                   np.where((v.p < 0.05) & (v.effect < -0.1), _JF_SIG["down"], _JF_SIG["ns"]))
        ax.scatter(v.effect, v.nlp, c=v.col, s=48, edgecolor="black", lw=0.5, rasterized=True)
        ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color=CFG.VIS_INK["faint"]); ax.axvline(0, ls="-", lw=0.5, color=CFG.VIS_INK["faint"])
        # Labels in a clean column in the empty right margin, ordered by effect so the thin leaders
        # fan out without crossing the cloud.
        _emax = float(v.effect.max()); _emin = float(v.effect.min()); _espan = (_emax - _emin) or 1.0
        top = v.sort_values("nlp", ascending=False).head(8).sort_values("effect", ascending=False)
        _xr = _emax + 0.14 * _espan
        _yr = np.linspace(v.nlp.max() * 0.95, v.nlp.max() * 0.28, len(top))
        for (_, r), _yy in zip(top.iterrows(), _yr):
            ax.annotate(_jf_ml(r.metric), (r.effect, r.nlp), (_xr, _yy), fontsize=_JF_FA,
                        ha="left", va="center", color=r.col if r.col != _JF_SIG["ns"] else CFG.VIS_INK["near_black"],
                        arrowprops=dict(arrowstyle="-", lw=0.5, color=CFG.VIS_INK["palest"]))
        ax.set_xlim(right=_xr + 0.55 * _espan)
        ax.set_xlabel("Effect size  (Cohen's d)"); ax.set_ylabel(r"$-\log_{10}$ p  (Mann–Whitney)"); _jf_box(ax)
        _jf_save(fig, folder_dir, "10_Volcano_Metric_Significance.png")


def _jf_manhattan(df, folder_dir, reporter):
    """Manhattan / TWAS — each candidate a locus grouped by ligand; y = rank-based −log10 significance."""
    with plt.rc_context(_JF_RC):
        fig = plt.figure(figsize=(8.5, 4.5)); ax = fig.subplots()
        d = df.copy(); d["_lg"] = _jf_ligshort(d)
        d["nlp"] = -np.log10((d["Scientific_Rank"].rank(pct=True)).clip(lower=1e-6))
        d = d.sort_values(["_lg", "Scientific_Rank"]); ligs = list(d["_lg"].unique()); xpos = 0; ticks = []
        for i, lg in enumerate(ligs):
            s = d[d._lg == lg]; x = np.arange(xpos, xpos + len(s))
            ax.scatter(x, s["nlp"], s=3, color=_JF_OKABE[i % len(_JF_OKABE)], alpha=0.7, edgecolor="none", rasterized=True)
            ticks.append(xpos + len(s)/2); xpos += len(s) + 40
        ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color=_JF_SIG["up"])
        ax.set_xticks(ticks); ax.set_xticklabels(ligs, rotation=90, fontsize=_JF_FD)
        ax.set_ylabel(r"$-\log_{10}$ p  (rank-based)"); ax.set_xlabel("Ligand locus"); ax.margins(x=0.01); _jf_box(ax)
        _jf_save(fig, folder_dir, "11_Manhattan_Candidate_Significance.png")


def _jf_importance(df, folder_dir, reporter):
    """Feature importance (PC1 loading, or metric SD fallback) + tier-composition streamgraph across ligands."""
    with plt.rc_context(_JF_RC):
        fig = plt.figure(figsize=(11, 5)); axL, axR = fig.subplots(1, 2, gridspec_kw={"width_ratios": [1, 1.5], "wspace": 0.22})
        pca = None
        _pcap = Path(folder_dir).parent / "01_Analysis_Data" / "02_PCA_Loadings.csv"
        if _pcap.exists():
            try:
                pca = pd.read_csv(_pcap, index_col=0)
            except Exception:
                pca = None
        if pca is not None and "PC1" in pca.columns:
            imp = pca["PC1"].abs().sort_values()
            axL.barh(range(len(imp)), imp.values, edgecolor="black", lw=0.4,
                     color=plt.cm.viridis(np.linspace(0.12, 0.9, len(imp))))   # distinct colour per feature
            axL.set_yticks(range(len(imp))); axL.set_yticklabels([_jf_ml(c) for c in imp.index], fontsize=_JF_FA)
            axL.set_xlabel("|PC1 loading|")
        else:
            var = _jf_metric_matrix(df).std().sort_values()
            axL.barh(range(len(var)), var.values, edgecolor="black", lw=0.4,
                     color=plt.cm.viridis(np.linspace(0.12, 0.9, len(var))))   # distinct colour per feature
            axL.set_yticks(range(len(var))); axL.set_yticklabels([_jf_ml(c) for c in var.index], fontsize=_JF_FA)
            axL.set_xlabel("Metric SD (importance)")
        _jf_box(axL)
        d = df.copy(); d["_lg"] = _jf_ligshort(d)
        order = d.groupby("_lg", observed=True)["iptm"].mean().sort_values(ascending=False).index.tolist()
        ct = (d.groupby(["_lg", "degrader_tier"], observed=True).size().unstack(fill_value=0).reindex(order).fillna(0))
        ct = ct[[t for t in _JF_TORDER if t in ct.columns]]; x = np.arange(len(ct)); base = -ct.sum(axis=1).values/2
        for t in ct.columns:
            axR.fill_between(x, base, base + ct[t].values, color=_JF_TCOL[t],
                             label=t.replace("Tier_", "").replace("_Decoy", "-D"), lw=0); base = base + ct[t].values
        axR.set_xticks(x); axR.set_xticklabels(ct.index, rotation=90, fontsize=_JF_FA)
        axR.set_ylabel("Candidates (stream)"); axR.set_xlabel("Ligand")
        axR.legend(ncol=len(ct.columns), fontsize=_JF_FA - 0.5, loc="lower center",
                   bbox_to_anchor=(0.5, 1.0), handlelength=0.9, columnspacing=1.0,
                   borderpad=0.3); axR.grid(False)   # framealpha/handletextpad inherit the CFG legend SSOT
        for s in axR.spines.values():
            s.set_visible(False)
        _jf_save(fig, folder_dir, "12_Feature_Importance_Tier_Streams.png")


def _jf_swimmer(df, folder_dir, reporter):
    """Swimmer — top 3 candidates per tier (Tier_1A → Tier_5_Decoy); bar = mechanistic score, colour = tier."""
    with plt.rc_context(_JF_RC):
        fig = plt.figure(figsize=(8, 6.5)); ax = fig.subplots()
        d = df.copy(); d["_lg"] = _jf_ligshort(d); tiers = _jf_tiers(d, d["degrader_tier"])
        picks = [d[d.degrader_tier == t].sort_values("Scientific_Rank").head(3) for t in tiers]
        s = pd.concat(picks).copy(); s["len"] = pd.to_numeric(s["mechanistic_score_effective"], errors="coerce").fillna(0)
        s = s.iloc[::-1]; y = np.arange(len(s))
        ax.barh(y, s["len"], color=[_JF_TCOL.get(str(t), CFG.VIS_INK["faint"]) for t in s["degrader_tier"]],
                edgecolor="black", lw=0.4, height=0.7)
        md = (s["MD_Selected"] == True).values if "MD_Selected" in s else np.zeros(len(s), bool)
        ax.scatter(s["len"].values[md] + 0.012, y[md], marker=">", s=34, color="black", zorder=3)
        deg = (s["is_degrader"] == True).values
        ax.scatter(np.full(deg.sum(), -0.012), y[deg], marker="o", s=18, color=_JF_SIG["up"], zorder=3)
        ax.set_yticks(y); ax.set_yticklabels([f"{str(t).replace('Tier_','').replace('_Decoy','-D')} · {lg}"
                                              for t, lg in zip(s["degrader_tier"], s["_lg"])], fontsize=_JF_FA)
        for _tl, _t in zip(ax.get_yticklabels(), s["degrader_tier"]):
            _tl.set_color(_JF_TCOL.get(str(_t), CFG.VIS_INK["near_black"]))   # tick label matches its bar's tier colour
        ax.set_xlabel("Mechanistic score (effective)"); ax.set_ylabel("Tier · ligand  (top 3 per tier)")
        leg = ax.legend(handles=[_JF_Line([0], [0], marker=">", ls="", mfc="black", mec="black", label="MD-selected"),
                                 _JF_Line([0], [0], marker="o", ls="", mfc=_JF_SIG["up"], mec=_JF_SIG["up"], label="degrader")],
                        loc="lower right", fontsize=_JF_FA, frameon=True, edgecolor="none")  # framealpha ← CFG SSOT
        leg.get_frame().set_facecolor("white"); ax.margins(y=0.01); _jf_box(ax)
        _jf_save(fig, folder_dir, "16_Swimmer_Top_Per_Tier.png")


def _jf_circos(df, folder_dir, reporter):
    """Circos / chord — ligand → tier assignment; ribbon width ∝ candidate count."""
    def _bez(ax, a0, a1, color, lw):
        p0 = np.array([np.cos(a0), np.sin(a0)]); p1 = np.array([np.cos(a1), np.sin(a1)])
        ax.add_patch(_JF_PathPatch(_JF_MPath([p0, p0*0.15, p1*0.15, p1],
                     [_JF_MPath.MOVETO, _JF_MPath.CURVE4, _JF_MPath.CURVE4, _JF_MPath.CURVE4]),
                     fc="none", ec=color, lw=lw, alpha=0.5))
    from matplotlib.patches import Wedge as _Wedge
    _LIG_ARC = CFG.VIS_INK["faint"]          # ligands share ONE neutral arc colour; tiers keep their own
    with plt.rc_context(_JF_RC):
        fig = plt.figure(figsize=(8, 8)); ax = fig.subplots(); ax.set_aspect("equal"); ax.axis("off")
        d = df.copy(); d["_lg"] = _jf_ligshort(d)
        ct = d.groupby(["_lg", "degrader_tier"], observed=True).size().unstack(fill_value=0)
        ct = ct[[t for t in _JF_TORDER if t in ct.columns]]; ligs = ct.index.tolist(); tiers = ct.columns.tolist()
        nodes = ligs + tiers; N = len(nodes); ang = {n: 2*np.pi*i/N for i, n in enumerate(nodes)}
        _half = np.degrees(np.pi / N) * 0.90    # half angular slot, small gap between wedges
        for n, a in ang.items():
            is_t = n in _JF_TCOL; col = _JF_TCOL.get(n, _LIG_ARC)
            ax.plot([0.97*np.cos(a), 1.03*np.cos(a)], [0.97*np.sin(a), 1.03*np.sin(a)],
                    color=col, lw=8 if is_t else 3.5, solid_capstyle="butt")
            de = np.degrees(a) % 360; rot = np.degrees(a); rot = rot-180 if 90 < de < 270 else rot
            # Every node carries an outer arc so the ring reads as an elegant band: each TIER in its
            # own colour, all LIGANDS in one shared neutral colour. Labels sit just beyond the arc.
            ax.add_patch(_Wedge((0, 0), 1.19, de - _half, de + _half, width=0.055,
                                facecolor=col if is_t else _LIG_ARC, edgecolor="none",
                                alpha=0.95 if is_t else 0.85, zorder=2))
            ax.text(1.235*np.cos(a), 1.235*np.sin(a), str(n).replace("Tier_", "T").replace("_Decoy", "5-D"),
                    ha="left" if de <= 90 or de >= 270 else "right", va="center", fontsize=_JF_FA,
                    rotation=rot, rotation_mode="anchor", fontweight="bold" if is_t else "normal",
                    color=col if is_t else CFG.VIS_INK["near_black"])   # tier labels match their arc colour
        mx = ct.values.max()
        for lg in ligs:
            for t in tiers:
                w = ct.loc[lg, t]
                if w > 0:
                    _bez(ax, ang[lg], ang[t], _JF_TCOL[t], 0.3 + 3.0*w/mx)
        ax.set_xlim(-1.5, 1.5); ax.set_ylim(-1.5, 1.5)
        _jf_save(fig, folder_dir, "17_Circos_Ligand_Tier_Assignment.png")


def _jf_network(df, folder_dir, reporter):
    """Metric co-variation — a circular chord diagram. Metrics sit on a ring; a chord joins any
    pair with |Pearson r| ≥ 0.5, coloured by the SIGN of r (red = positive, blue = negative) and
    widened by |r|. Node size ∝ how many strong correlations the metric has; labels sit outside the
    ring so nothing overlaps."""
    from matplotlib.colors import Normalize as _Norm
    THR = 0.5
    with plt.rc_context(_JF_RC):
        fig = plt.figure(figsize=(9, 9)); ax = fig.subplots(); ax.set_aspect("equal"); ax.axis("off")
        mm = _jf_metric_matrix(df); m = mm.sample(min(5000, len(mm)), random_state=3)
        corr = m.corr(); cols = list(corr.columns); N = len(cols)
        ang = {c: 2 * np.pi * i / N for i, c in enumerate(cols)}
        P = {c: np.array([np.cos(a), np.sin(a)]) for c, a in ang.items()}
        _deg = {c: 0 for c in cols}
        _edges = []
        for i in range(N):
            for j in range(i + 1, N):
                r = corr.iloc[i, j]
                if np.isfinite(r) and abs(r) >= THR:
                    _edges.append((cols[i], cols[j], float(r))); _deg[cols[i]] += 1; _deg[cols[j]] += 1
        _norm = _Norm(-1, 1); _cm = plt.get_cmap("RdBu_r")
        # weakest chords first so the strong ones read on top
        for u, v, r in sorted(_edges, key=lambda e: abs(e[2])):
            p0, p1 = P[u], P[v]
            ax.add_patch(_JF_PathPatch(_JF_MPath([p0, p0 * 0.28, p1 * 0.28, p1],
                         [_JF_MPath.MOVETO, _JF_MPath.CURVE4, _JF_MPath.CURVE4, _JF_MPath.CURVE4]),
                         fc="none", ec=_cm(_norm(r)), lw=0.7 + 3.4 * (abs(r) - THR) / (1 - THR),
                         alpha=0.72, zorder=1))
        for c in cols:
            a = ang[c]; p = P[c]
            ax.scatter(*p, s=70 + 60 * _deg[c], color=CFG.VIS_ACCENT["blue"],
                       edgecolor="white", lw=1.1, zorder=3)
            de = np.degrees(a) % 360; rot = np.degrees(a); rot = rot - 180 if 90 < de < 270 else rot
            ax.text(1.10 * np.cos(a), 1.10 * np.sin(a), _jf_ml(c),
                    ha="left" if de <= 90 or de >= 270 else "right", va="center",
                    fontsize=_JF_FA, rotation=rot, rotation_mode="anchor")
        from matplotlib.lines import Line2D as _L15
        ax.legend(handles=[_L15([0], [0], color=_cm(_norm(0.85)), lw=3.2, label="positive correlation"),
                           _L15([0], [0], color=_cm(_norm(-0.85)), lw=3.2, label="negative correlation")],
                  loc="lower center", bbox_to_anchor=(0.5, -0.02), ncol=2, frameon=False,
                  title=f"chord = |r| ≥ {THR:g}  ·  width ∝ |r|  ·  node size ∝ # strong links")
        ax.set_xlim(-1.45, 1.45); ax.set_ylim(-1.45, 1.45)
        _jf_save(fig, folder_dir, "15_Metric_CoVariation_Network.png")


def _jf_treemap(df, folder_dir, reporter):
    """Treemap — candidates by tier × ligand; rectangle area ∝ count, colour = tier (squarified)."""
    def _sq(sizes, x, y, dx, dy):
        sizes = list(sizes); total = float(sum(sizes))
        if total <= 0:
            return []
        sizes = [s*dx*dy/total for s in sizes]
        def worst(row, L): s = sum(row); return max((L**2)*max(row)/(s**2), (s**2)/((L**2)*min(row)))
        def lay(row, x, y, dx, dy):
            cov = sum(row); r = []
            if dx >= dy:
                w = cov/dy; yy = y
                for v in row: h = v/w; r.append({"x": x, "y": yy, "dx": w, "dy": h}); yy += h
                return r, x+w, y, dx-w, dy
            h = cov/dx; xx = x
            for v in row: w = v/h; r.append({"x": xx, "y": y, "dx": w, "dy": h}); xx += w
            return r, x, y+h, dx, dy-h
        out = []; row = []; rx, ry, rdx, rdy = x, y, dx, dy
        for s in sizes:
            L = min(rdx, rdy)
            if not row or worst(row+[s], L) <= worst(row, L): row.append(s)
            else: rr, rx, ry, rdx, rdy = lay(row, rx, ry, rdx, rdy); out += rr; row = [s]
        if row: rr, *_ = lay(row, rx, ry, rdx, rdy); out += rr
        return out
    with plt.rc_context(_JF_RC):
        fig = plt.figure(figsize=(9, 6.5)); ax = fig.subplots(); ax.axis("off")
        d = df.copy(); d["_lg"] = _jf_ligshort(d)
        cnt = d.groupby(["degrader_tier", "_lg"], observed=True).size().reset_index(name="n")
        cnt = cnt[cnt.n > 0].sort_values("n", ascending=False)
        ax.set_xlim(0, 100); ax.set_ylim(0, 100)   # set before labels so the fit calc sees real box sizes
        for r, (_, row) in zip(_sq(cnt["n"].values, 0, 0, 100, 100), cnt.iterrows()):
            c = _JF_TCOL.get(str(row["degrader_tier"]), CFG.VIS_INK["faint"])
            ax.add_patch(_JF_Rect((r["x"], r["y"]), r["dx"], r["dy"], facecolor=c, edgecolor="white", lw=0.8))
            _lines = [str(row["_lg"]), str(int(row["n"]))]
            # base = font-size cap (big boxes never exceed it); 3.0 = floor (smaller → label skipped)
            _fs = _jf_fit_fs(ax, r["x"]+r["dx"]/2, r["y"]+r["dy"]/2, r["dx"], r["dy"], _lines, _JF_FA + 1.0, 3.0)
            if _fs is not None:   # font sized to the box; skipped when it would be unreadably small
                ax.text(r["x"]+r["dx"]/2, r["y"]+r["dy"]/2, "\n".join(_lines), ha="center", va="center",
                        fontsize=_fs, color=_text_color(c))
        ax.legend(handles=[_JF_Patch(color=_JF_TCOL[t], label=t.replace("Tier_", "T").replace("_Decoy", "5-D"))
                           for t in _JF_TORDER], ncol=7, fontsize=_JF_FA, loc="upper center",
                  bbox_to_anchor=(0.5, -0.01), handlelength=1)
        _jf_save(fig, folder_dir, "06_Treemap_Tier_Ligand_Composition.png")


# =============================================================================
# STEP 3/8 — 03_Dataset_and_Alignment_Overview
# =============================================================================
def _fig_folder03_dataset(df, features, out_dir, reporter, existing_tiers):
    """Folder 03_Dataset_and_Alignment_Overview — coverage, tiers, alignment grades."""
    reporter.section("Step 3/8 — 03_Dataset_and_Alignment_Overview · coverage, tiers, alignment grades")
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
        # The companion 01_Active_Site_Residue_Mapping_Missed.csv is written up front in Step 1
        # (write_residue_mapping_missed_csv); this figure only draws the coverage bars.

        fig, ax = plt.subplots(figsize=(11, 6))
        # Colour each bar by its functional role group (CFG single source) so residues
        # of the same mechanistic role share a colour.
        _groups01 = [CFG.ACTIVE_SITE_ROLE_GROUP.get(k, "") for k in _keys]
        _colours01 = [CFG.ACTIVE_SITE_ROLE_GROUP_COLOUR.get(g, CFG.VIS_INK["faint"]) for g in _groups01]
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
        # Each tick label takes its own bar's role colour: the residue, its bar and the legend
        # entry for its mechanistic role then carry one identity the eye can follow.
        for _tl01, _c01 in zip(ax.get_xticklabels(), _colours01):
            _tl01.set_color(_c01)
        ax.set_ylabel("Variants with residue mapped")
        ax.set_ylim(0, _n01 * 1.08 if _n01 else 1)
        # Legend: one entry per role group (in mechanistic order) + the all-variants line,
        # laid out as a single row in the top-left.
        from matplotlib.patches import Patch as _Patch01
        _seen01 = []
        for g in _groups01:
            if g not in _seen01:
                _seen01.append(g)
        _handles01 = [_Patch01(facecolor=CFG.ACTIVE_SITE_ROLE_GROUP_COLOUR.get(g, CFG.VIS_INK["faint"]),
                               edgecolor="white", label=g) for g in _seen01]
        _handles01.append(Line2D([0], [0], color="grey", linestyle="--", label=f"All variants (n={_n01:,})"))
        ax.legend(handles=_handles01, loc="upper left", ncol=len(_handles01),
                   frameon=True)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_01_Active_Site_Residue_Mapping_Coverage.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()
    else:
        reporter.log(f"  ! Skipped Figure 01: Alignment_Stats.csv not found at {_aln_csv}")

    # --- Figure 02: Tier Distribution + Model Selection (pie inset) ---
    if CFG.COL_TIER in df.columns:
        total_complexes = len(df)
        model_col = "best_model_name" if "best_model_name" in df.columns else None
        counts6 = [len(df[df[CFG.COL_TIER] == t]) for t in existing_tiers]

        # Adaptive figure height: scale with the tallest bar so the pie doesn't float
        _max_cnt6 = max(counts6) if counts6 else 1
        _fig_h6 = max(5.5, min(10.0, 4.5 + _max_cnt6 / 12000))

        fig, ax = plt.subplots(figsize=(11, _fig_h6))
        bar_x6 = np.arange(len(existing_tiers))

        ax.bar(bar_x6, counts6,
               color=[TIER_PALETTE.get(t, CFG.VIS_INK["faint"]) for t in existing_tiers],
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
        ax.set_xlabel("Degrader Tier", )
        ax.set_ylabel("Number of complexes", )
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
            """
            The pie says which Boltz diffusion model produced each pose. A gold star marks the wedges
            the MD-SELECTED complexes were actually drawn from, so the reader can see at a glance
            whether the cohort committed to simulation came from one model or several. A cohort that
            all came from a single diffusion sample is a different claim from one spread across
            models: the first could be a single-model artefact, the second is model-independent.
            """
            _md6 = _md_ready_df(df)
            _md_models6 = set()
            if not _md6.empty and model_col in _md6.columns:
                _md_models6 = {str(_v) for _v in _md6[model_col].dropna().unique()}
            _ctrl6 = _control_star_df(df)
            _ctrl_models6 = ({str(_v) for _v in _ctrl6[model_col].dropna().unique()}
                             if (not _ctrl6.empty and model_col in _ctrl6.columns) else set())

            for _w6, _lab6, _col6, _at6, _key6 in zip(wedges6, _model_labels6, _pie_colours,
                                                      autotexts, model_counts.index):
                _tc6 = _text_color(_col6)
                _at6.set_fontsize(7.5)
                _at6.set_color(_tc6)
                _ang6 = np.deg2rad((_w6.theta1 + _w6.theta2) / 2.0)
                ax_pie.text(0.82 * np.cos(_ang6), 0.82 * np.sin(_ang6), _lab6,
                            ha="center", va="center", fontsize=_JF_FA,
                            fontweight="bold", color=_tc6)
                if str(_key6) in _md_models6:
                    """
                    ONE STAR PER MD-SELECTED COMPLEX. The count is the number of stars — no text
                    label, because the mark already carries the whole message: how many of the
                    simulated complexes this diffusion model produced. Wedge borders are left
                    untouched; outlining only the starred wedges reads as a rendering artefact,
                    since a wedge edge is shared with its neighbour and the highlight bleeds.
                    """
                    _n_md6 = int((_md6[model_col].astype(str) == str(_key6)).sum())
                    _rad6 = 0.42
                    _sx6, _sy6 = _rad6 * np.cos(_ang6), _rad6 * np.sin(_ang6)
                    # Stars laid out along the wedge's tangent so they never sit on the % text.
                    _tx6, _ty6 = -np.sin(_ang6), np.cos(_ang6)
                    _gap6 = 0.115
                    _off6 = (np.arange(_n_md6) - (_n_md6 - 1) / 2.0) * _gap6
                    ax_pie.scatter(_sx6 + _off6 * _tx6, _sy6 + _off6 * _ty6,
                                   marker="*", s=150, facecolor=CFG.VIS_ACCENT["star"],
                                   edgecolor="black", linewidths=0.9, zorder=6)
                # 3R3U × FA control: red/gold star at the wedge of the model it came from
                # (drawn even when that model produced no candidate, so the control is never hidden).
                # Placed at r=0.26 — the empty inner zone, clear of the % autotext (r≈0.6) and the
                # candidate stars (r=0.42) so it never sits on the percentage label.
                if str(_key6) in _ctrl_models6:
                    _crad6 = 0.26
                    ax_pie.scatter([_crad6 * np.cos(_ang6)], [_crad6 * np.sin(_ang6)],
                                   marker="*", s=175, facecolor=CFG.VIS_ACCENT["control"],
                                   edgecolor=CFG.VIS_ACCENT["control_edge"], linewidths=1.2, zorder=7)
            ax_pie.set_frame_on(False)

        # Gold star on the tier the MD-ready cohort sits in — the bar the pipeline acts on.
        _md6_bar = _md_ready_df(df)
        if not _md6_bar.empty and CFG.COL_TIER in _md6_bar.columns:
            _md6_tiers = _md6_bar[CFG.COL_TIER].astype(str).value_counts()
            for _t6, _n6 in _md6_tiers.items():
                if _t6 not in existing_tiers:
                    continue
                _idx6 = existing_tiers.index(_t6)
                _xi6, _yi6 = bar_x6[_idx6], float(counts6[_idx6])
                ax.scatter([_xi6], [_yi6 + _y_ceil6 * 0.115], **_MD_STAR_KW)
                _cbar6 = _control_star_df(df)
                if not _cbar6.empty and _t6 in _cbar6[CFG.COL_TIER].astype(str).values:
                    _cbw6 = (bar_x6[1] - bar_x6[0]) if len(bar_x6) > 1 else 0.8
                    _ctrl_star(ax, _xi6 + _cbw6 * 0.28, _yi6 + _y_ceil6 * 0.115, size=330)
                ax.text(_xi6, _yi6 + _y_ceil6 * 0.155, f"MD-selected\n(n={_n6})",
                        ha="center", va="bottom", fontsize=7.0, fontweight="bold",
                        color=CFG.VIS_ACCENT["star_edge"], zorder=9)

        plt.tight_layout()
        plt.savefig(out_dir / "Figure_02_Tier_Distribution.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # --- Figure 03: Sequence Identity — Histogram (grade bands) overlaid with tier KDE lines ---
    """
    Single chart: histogram bars show the overall grade distribution (how many proteins
    per identity band); KDE lines per tier overlay on the same x-axis, revealing whether
    different catalytic tiers cluster at different identity levels.
    """
    if CFG.COL_ID_PCT in df.columns:
        # Use existing single-letter Alignment_Grade from CSV; only compute if missing
        if CFG.COL_ALN_G not in df.columns or df[CFG.COL_ALN_G].isnull().all():
            def _simple_grade(v):
                try:
                    val = float(v)
                except (ValueError, TypeError):
                    return "I"
                return _utils_mod.get_alignment_grade(val, CFG)
            df[CFG.COL_ALN_G] = df[CFG.COL_ID_PCT].apply(_simple_grade)
        # Ensure single-letter format (convert "Grade A (...)" → 'A' if needed)
        _ag = df[CFG.COL_ALN_G].astype(str)
        if _ag.str.startswith("Grade ").any():
            df[CFG.COL_ALN_G] = _ag.str.extract(r"Grade\s+([A-I])", expand=False).fillna("I")
        id_plot = df.dropna(subset=[CFG.COL_ID_PCT]).copy()
        id_plot[CFG.COL_ID_PCT] = pd.to_numeric(id_plot[CFG.COL_ID_PCT], errors="coerce").clip(0, 100)
        id_plot = id_plot.dropna(subset=[CFG.COL_ID_PCT])

        fig, ax = plt.subplots(figsize=(13, 7))
        ax2 = ax.twinx()   # secondary y-axis for KDE density

        # Grade band background shading — sourced from CFG § 8.8
        grade_bands = list(CFG.GRADE_BANDS)
        # Per-protein grade counts (one row per unique protein)
        _prot_dedup02 = id_plot.drop_duplicates(subset=[CFG.COL_PROT]) if CFG.COL_PROT in id_plot.columns else id_plot
        _n_total_prot02 = len(_prot_dedup02)
        for lo, hi, col, lbl in grade_bands:
            ax.axvspan(lo, hi, alpha=0.20, color=col, zorder=0)
            ax.text((lo + hi) / 2, 0.972, f"Grade {lbl}", ha="center", va="top",
                    fontsize=_JF_FA, color=col, fontweight="bold",
                    transform=ax.get_xaxis_transform(),
                    bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none",
                              alpha=0.70), zorder=11)
            _cnt_g02 = int((_prot_dedup02[CFG.COL_ALN_G] == lbl).sum())
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
            sub_t = id_plot[id_plot[CFG.COL_TIER] == tier]
            if sub_t.empty:
                continue
            tier_hist_data.append(sub_t[CFG.COL_ID_PCT].values)
            tier_hist_labels.append(f"{tier}  (n={len(sub_t):,})")
            tier_hist_colors.append(TIER_PALETTE.get(tier, CFG.VIS_INK["faint"]))
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
                            if len(id_plot[id_plot[CFG.COL_TIER] == t]) >= 10]
        _ridge_step = 0.008   # vertical offset per tier index
        _tt_kde_med_xy = None   # captured KDE-median dot position of the top tier (on ax2)
        for _tier_idx, tier in enumerate(_kde_tiers_valid):
            sub_t = id_plot[id_plot[CFG.COL_TIER] == tier]
            kde_x = np.linspace(0, 100, 400)
            try:
                kde_y = _gkde(sub_t[CFG.COL_ID_PCT])(kde_x)
            except (np.linalg.LinAlgError, ValueError):
                # Zero-variance tier (all sequences share one identity) → singular
                # covariance; skip this ridge rather than crash the whole figure.
                continue
            _tier_offset = _tier_idx * _ridge_step
            kde_y_shifted = kde_y + _tier_offset
            _tcol = TIER_PALETTE.get(tier, CFG.VIS_INK["faint"])
            ax2.plot(kde_x, kde_y_shifted, color=_tcol,
                     linewidth=2.5, alpha=0.75, zorder=5 + _tier_idx,
                     linestyle="-")
            # Fill under each shifted KDE curve for ridge visual separation
            ax2.fill_between(kde_x, _tier_offset, kde_y_shifted,
                             color=_tcol, alpha=0.15, zorder=4 + _tier_idx)
            # Median tick mark on the KDE line
            med_x = float(sub_t[CFG.COL_ID_PCT].median())
            med_y = float(np.atleast_1d(_gkde(sub_t[CFG.COL_ID_PCT])(np.array([med_x])))[0])
            ax2.scatter([med_x], [med_y + _tier_offset], color=_tcol,
                        s=55, zorder=6 + _tier_idx, edgecolors="black", linewidths=0.8)
            if tier == CFG.TIER_TOP:
                _tt_kde_med_xy = (med_x, med_y + _tier_offset)

        # Anchor KDE axis at zero; show full density range with 10% headroom
        _kde_ymax = ax2.get_ylim()[1]
        ax2.set_ylim(0, _kde_ymax * 1.10)

        # Grade threshold vertical lines
        for pct_v, _lbl_v, col_v in [(90, "Grade A boundary", CFG.VIS_ACCENT["good"]),
                                      (60, "Grade D boundary", CFG.VIS_ACCENT["amber"]),
                                      (40, "Grade F boundary", CFG.VIS_ACCENT["vermillion"])]:
            ax.axvline(x=pct_v, color=col_v, linestyle="--", alpha=0.65,
                       linewidth=1.3, zorder=4)

        ax.set_xlim(0, 100)
        ax.set_xticks(range(0, 101, 10))
        ax.tick_params(axis="x", labelsize=9)
        # Left axis (count): steel-blue tick labels matching the grey count gridlines
        ax.tick_params(axis="y", labelsize=9, labelcolor=CFG.VIS_ACCENT["axis_left"])
        ax.yaxis.label.set_color(CFG.VIS_INK["near_black"])
        # Right axis (KDE density): amber tick labels matching the amber density gridlines
        ax2.tick_params(axis="y", labelsize=9, labelcolor=CFG.VIS_ACCENT["axis_right"])
        ax2.yaxis.label.set_color(CFG.VIS_INK["near_black"])
        ax.set_xlabel("Sequence Identity to Reference DeHa4 Control  (%)", )
        ax.set_ylabel("Number of Protein–Ligand Complexes  (stacked by tier)", )
        ax2.set_ylabel("KDE Density  (probability density)",  rotation=270, labelpad=14)

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
        ax.yaxis.grid(True, color=CFG.VIS_ACCENT["axis_left"], linewidth=0.55, linestyle="-", alpha=0.25, zorder=0)
        ax.xaxis.grid(False)
        ax2.yaxis.grid(True, color=CFG.VIS_ACCENT["axis_right"], linewidth=0.55, linestyle="--", alpha=0.20, zorder=0)

        _ymax02 = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 200
        # Convergence identity is COMPUTED from the data (median identity of the degrader tiers —
        # where the tier distributions actually overlap), not a hard-coded value, and the arrow is
        # anchored to that x so it points at the real convergence, not an arbitrary 60%.
        _deg_t02 = [t for t in CFG.TIER_ORDER if t not in (CFG.TIER_DECOY, CFG.TIER_POOR, "Tier_4")]
        _conv_src02 = (id_plot[id_plot[CFG.COL_TIER].astype(str).isin(_deg_t02)]
                       if CFG.COL_TIER in id_plot.columns else id_plot)
        _conv02 = float(pd.to_numeric(_conv_src02[CFG.COL_ID_PCT], errors="coerce").median())
        if _conv02 == _conv02:   # guard NaN
            ax.annotate(f"Tiers converge\nnear {_conv02:.0f}% identity",
                        xy=(_conv02, _ymax02 * 0.55),
                        xytext=(min(_conv02 + 24, 80), _ymax02 * 0.82),
                        fontsize=7.0, color=CFG.VIS_INK["soft"], style="italic", ha="left",
                        arrowprops=dict(arrowstyle="->", color=CFG.VIS_INK["ghost"], lw=0.9))
        # Tier_1A visibility annotation — arrow points to the Tier_1A KDE-curve median
        # (a visible dot), not the near-invisible count bar. Drawn on the KDE axis and
        # lifted above the ridge lines with an opaque box so it never hides behind them.
        _tt_sub02 = id_plot[id_plot[CFG.COL_TIER] == CFG.TIER_TOP]
        if len(_tt_sub02) > 0 and _tt_kde_med_xy is not None:
            _tt_med_x = float(_tt_sub02[CFG.COL_ID_PCT].median())
            ax2.annotate(
                f"{CFG.TIER_TOP}  (n={len(_tt_sub02):,})\nmedian identity {_tt_med_x:.0f}%",
                xy=_tt_kde_med_xy, xycoords="data",
                xytext=(0.80, 0.40), textcoords=ax2.transAxes,
                fontsize=_JF_FA, ha="left", va="center",
                color=TIER_PALETTE.get(CFG.TIER_TOP, CFG.VIS_BAND["high"]),
                fontweight="bold",
                arrowprops=dict(arrowstyle="->", lw=1.4,
                                color=TIER_PALETTE.get(CFG.TIER_TOP, CFG.VIS_BAND["high"])),
                bbox=dict(boxstyle="round,pad=0.25", fc="white",
                          ec=TIER_PALETTE.get(CFG.TIER_TOP, CFG.VIS_BAND["high"]),
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
            _tcol  = TIER_PALETTE.get(_tname, CFG.VIS_INK["faint"])
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
                _tcol  = TIER_PALETTE.get(_tname, CFG.VIS_INK["faint"])
                _kdot  = _L7([0],[0], marker="o", color=_tcol, linewidth=1.5,
                             markerfacecolor=_tcol, markeredgecolor="black",
                             markeredgewidth=0.5, markersize=5)
                _pair_handles.append((_bh, _kdot))
                _pair_labels.append(_tname)
            """
            The MD-selected cohort is drawn BEFORE the legend is built, so its star can join the
            existing tier legend as one more entry. A second legend box would say nothing the first
            one could not, and two boxes stacked in the same corner cost more space than the marks
            they explain.
            """
            _md7 = _md_ready_df(df)
            if not _md7.empty and CFG.COL_ID_PCT in _md7.columns:
                _mdx7 = pd.to_numeric(_md7[CFG.COL_ID_PCT], errors="coerce").dropna()
                if len(_mdx7):
                    _y7 = ax.get_ylim()[1] * 0.035
                    ax.scatter(_mdx7.to_numpy(), np.full(len(_mdx7), _y7), **_MD_STAR_KW)
                    for _xv7 in _mdx7:
                        ax.annotate("", xy=(float(_xv7), 0), xytext=(float(_xv7), _y7),
                                    arrowprops=dict(arrowstyle="-", color=CFG.VIS_ACCENT["star_edge"],
                                                    lw=0.9, alpha=0.85), zorder=8)
                    _pair_handles.append(_L7([0], [0], marker="*", linestyle="none",
                                             markersize=13, markerfacecolor=CFG.VIS_ACCENT["star"],
                                             markeredgecolor="black", markeredgewidth=1.0))
                    _pair_labels.append(f"MD-selected (n={len(_mdx7)})")
                    _ctrl7 = _control_star_df(df)
                    if not _ctrl7.empty and CFG.COL_ID_PCT in _ctrl7.columns:
                        _cx7 = pd.to_numeric(_ctrl7[CFG.COL_ID_PCT], errors="coerce").dropna()
                        if len(_cx7):
                            _ctrl_star(ax, _cx7.to_numpy(), np.full(len(_cx7), _y7))
                            _pair_handles.append(_L7([0], [0], marker="*", linestyle="none",
                                                     markersize=14, markerfacecolor=CFG.VIS_ACCENT["control"],
                                                     markeredgecolor=CFG.VIS_ACCENT["control_edge"], markeredgewidth=1.5))
                            _pair_labels.append("3R3U × FA (control)")

            _leg7 = ax.legend(_pair_handles, _pair_labels,
                              loc="upper left", ncol=3,
                               fancybox=True,
                              handler_map={tuple: _HT02(ndivide=None, pad=0.2)},
                              title="Tier (bars = counts, KDE lines ● median)",
                              
                              bbox_to_anchor=(0.01, 0.93))
        except Exception:
            _leg7 = ax.legend(_bar_hdls_02, _bar_lbls_02,
                              loc="upper left", ncol=3,
                               fancybox=True,
                              title="Tier (bars = counts, KDE lines ● median)",
                              
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
    if CFG.COL_ID_PCT in df.columns and CFG.COL_TIER in df.columns:
        try:
            f11_df = df.dropna(subset=[CFG.COL_ID_PCT]).copy()
            f11_df[CFG.COL_ID_PCT] = pd.to_numeric(f11_df[CFG.COL_ID_PCT], errors="coerce").clip(0, 100)
            f11_df = f11_df.dropna(subset=[CFG.COL_ID_PCT])
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
            f11_df["Grade"] = pd.cut(f11_df[CFG.COL_ID_PCT], bins=grade_cuts,
                                     labels=grade_lbls, right=False)
            ct11 = pd.crosstab(f11_df[CFG.COL_TIER], f11_df["Grade"])
            ct11 = ct11.reindex(index=[t for t in existing_tiers if t in ct11.index],
                                columns=grade_lbls, fill_value=0)
            ct11_pct = ct11.div(ct11.sum(axis=1), axis=0).fillna(0) * 100

            from scipy.stats import chi2_contingency as _chi2_cont
            try:
                _chi2_f11, _p_f11, _dof_f11, _ = _chi2_cont(ct11.values)
                _chi2_str = (f"χ²({_dof_f11}) = {_chi2_f11:.1f},  "
                             f"p {'< 0.001' if _p_f11 < 0.001 else f'= {_p_f11:.3f}'}")
                _register_p("Chi-square — grade × tier", "11_Grade_by_Tier",
                            float(_chi2_f11), int(ct11.values.sum()), float(_p_f11), dof=int(_dof_f11))
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
                        color=grade_colors.get(grade, CFG.VIS_INK["faint"]),
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
                    txt_col = "white" if grade_colors.get(grade, CFG.VIS_INK["faint"])[1:] < "AAAAAA" else CFG.VIS_INK["near_black"]
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
                        ha="center", va="bottom", fontsize=_JF_FA, color=col,
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

            # X-axis ticks every 10% — identical to Fig02
            ax.set_xlim(0, 100)
            ax.set_xticks(range(0, 101, 10))
            ax.set_xticklabels([f"{v}%" for v in range(0, 101, 10)],
                               fontsize=8, rotation=45, ha="right")
            ax.set_xlabel("Proportion of complexes in each sequence identity grade  (%)",
                          )
            ax.set_ylabel("Degrader Tier", labelpad=-8)   # sit closer to the tier tick labels
            ax.xaxis.grid(True, color=CFG.VIS_INK["tick"], linewidth=0.65, zorder=0)
            ax.set_axisbelow(True)

            """
            Gold stars mark the MD-selected cohort on its own tier row, positioned at the cumulative
            proportion of the identity grade each complex falls in. The row shows how that tier is
            spread across sequence-identity grades; the star shows where inside that spread the
            complexes actually taken to MD sit — which is the part of the distribution the downstream
            claims rest on.
            """
            _md11 = _md_ready_df(f11_df)
            if not _md11.empty and "Grade" in _md11.columns:
                """
                ONE STAR PER MD-SELECTED COMPLEX, placed in the GAP ABOVE the bar rather than on it.
                The segments already carry their percentage and count; a star dropped on the row lands
                on those numbers and obscures the very values it is meant to annotate. Sitting in the
                gap, the stars point at the segment they belong to without covering anything, and
                their number carries the count — no text label needed.
                """
                _tier_row11 = {t: bar_y[i] for i, t in enumerate(tiers_f11)}
                _grades11 = list(ct11.columns)
                _pct11 = ct11.div(ct11.sum(axis=1), axis=0) * 100.0
                _star_seen11 = 0
                for _t11, _grp11 in _md11.groupby(_md11[CFG.COL_TIER].astype(str)):
                    if _t11 not in _tier_row11:
                        continue
                    for _g11, _sub11 in _grp11.groupby(_grp11["Grade"].astype(str)):
                        if _g11 not in _grades11:
                            continue
                        _gi11 = _grades11.index(_g11)
                        # Mid-point of this grade's segment within the tier's stacked row.
                        _cum11 = float(_pct11.loc[_t11].iloc[:_gi11].sum())
                        _seg11 = float(_pct11.loc[_t11].iloc[_gi11])
                        _cx11 = _cum11 + _seg11 / 2.0
                        _n11 = len(_sub11)
                        _dx11 = (np.arange(_n11) - (_n11 - 1) / 2.0) * 2.2   # % units, side by side
                        ax.scatter(_cx11 + _dx11,
                                   np.full(_n11, _tier_row11[_t11] + 0.40),
                                   **({**_MD_STAR_KW, "s": 240}))
                        _star_seen11 += _n11
                _ctrl11 = _control_star_df(f11_df)
                if not _ctrl11.empty and "Grade" in _ctrl11.columns:
                    for _ct11, _cgrp11 in _ctrl11.groupby(_ctrl11[CFG.COL_TIER].astype(str)):
                        if _ct11 not in _tier_row11:
                            continue
                        for _cg11, _csub11 in _cgrp11.groupby(_cgrp11["Grade"].astype(str)):
                            if _cg11 not in _grades11:
                                continue
                            _cgi11 = _grades11.index(_cg11)
                            _ccx11 = (float(_pct11.loc[_ct11].iloc[:_cgi11].sum())
                                      + float(_pct11.loc[_ct11].iloc[_cgi11]) / 2.0)
                            _ctrl_star(ax, [_ccx11], [_tier_row11[_ct11] + 0.62], size=260)
                            _star_seen11 += 1
                if _star_seen11:
                    from matplotlib.lines import Line2D as _L2D11
                    # The band below the bottom row is already empty — the legend goes there as it
                    # is. Extending the axis to make room would open white space the figure does not
                    # need, which costs more than the legend it houses.
                    ax.legend([_L2D11([], [], marker="*", linestyle="none", markersize=13,
                                      markerfacecolor=CFG.VIS_ACCENT["star"], markeredgecolor="black")],
                              [f"MD-selected (n={_star_seen11})  ·  one star per complex"],
                              loc="lower left", bbox_to_anchor=(0.005, 0.005),
                                fancybox=True).set_zorder(11)

            if _chi2_str:
                ax.text(0.99, 0.01, _chi2_str, transform=ax.transAxes,
                        ha="right", va="bottom", fontsize=8.5, style="italic",
                        color=CFG.VIS_INK["dark"],
                        bbox=dict(boxstyle="round,pad=0.25", fc="white",
                                  ec=CFG.VIS_INK["palest"], alpha=0.88, linewidth=0.6))
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_04_Tier_Grade_Distribution.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('04')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig



# =============================================================================
# STEP 4/8 — 04_AI_Confidence_Quality
# =============================================================================
def _fig_folder04_ai_confidence(df, features, out_dir, reporter, existing_tiers, _pa, _imgs, _tt_has_imgs):
    """Folder 04_AI_Confidence_Quality — Boltz-2 confidence + pTM/ipTM."""
    reporter.section("Step 4/8 — 04_AI_Confidence_Quality · Boltz-2 confidence + pTM/ipTM")
    # --- Figure 05: AI Quality — Boltz Confidence (boxes) + pTM & ipTM (lines) overlaid ---
    """
    Single chart: box plots show per-tier distribution of Boltz Model Confidence (primary
    y-axis). Overlaid connected dot-lines show per-tier median pTM and ipTM on the same
    scale (both 0–1), making agreement and divergence between the three AI metrics visible.
    """
    _has_conf  = CFG.COL_CONF in df.columns and CFG.COL_TIER in df.columns
    _has_ptm   = "ptm" in df.columns and "iptm" in df.columns and CFG.COL_TIER in df.columns
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
            if _has_conf and not df.empty:
                sns.boxplot(data=df, x=CFG.COL_TIER, y=CFG.COL_CONF,
                            order=existing_tiers, palette=TIER_PALETTE, linewidth=1.2,
                            flierprops=dict(marker="o", markersize=2, alpha=0.25),
                            width=0.55, zorder=2, ax=ax)
                # Median text above each box
                tier_order_f10 = list(existing_tiers)
                medians_f10 = []
                for i, tier in enumerate(tier_order_f10):
                    med = df.loc[df[CFG.COL_TIER] == tier, CFG.COL_CONF].median()
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
                                fontsize=_JF_FA, color=CFG.VIS_ACCENT["warn"], style="italic",
                                arrowprops=dict(arrowstyle="->", color=CFG.VIS_ACCENT["warn"], lw=0.9),
                                bbox=dict(boxstyle="round,pad=0.2", fc=CFG.VIS_ACCENT_DEEP["rose_edge"],
                                          ec=CFG.VIS_ACCENT["warn"], alpha=0.9, lw=0.8))

            # pTM — connected median line with CI band
            if _has_ptm:
                ptm_xs, ptm_meds, ptm_lo, ptm_hi = [], [], [], []
                iptm_xs, iptm_meds, iptm_lo, iptm_hi = [], [], [], []
                for tier in existing_tiers:
                    sub = df.loc[df[CFG.COL_TIER] == tier].dropna(subset=["ptm", "iptm"])
                    if len(sub) < 3:
                        continue
                    x_pos = tier_x[tier]
                    for metric, xs_l, meds_l, lo_l, hi_l, col in [
                        ("ptm", ptm_xs, ptm_meds, ptm_lo, ptm_hi, CFG.VIS_ACCENT["blue"]),
                        ("iptm", iptm_xs, iptm_meds, iptm_lo, iptm_hi, CFG.VIS_ACCENT["magenta"]),
                    ]:
                        vals = sub[metric].dropna()
                        med_v = float(vals.median())
                        lo_ci, hi_ci = _median_ci95(vals)
                        xs_l.append(x_pos)
                        meds_l.append(med_v)
                        lo_l.append(lo_ci)
                        hi_l.append(hi_ci)

                for xs_l, meds_l, lo_l, hi_l, col, lbl, mk in [
                    (ptm_xs, ptm_meds, ptm_lo, ptm_hi, CFG.VIS_ACCENT["blue"], "pTM  (fold confidence)", "s"),
                    (iptm_xs, iptm_meds, iptm_lo, iptm_hi, CFG.VIS_ACCENT["magenta"], "ipTM  (interface confidence)", "^"),
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
                all_vals_f10.extend(df[CFG.COL_CONF].dropna().tolist())
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
                _bax_kw = dict(transform=ax.transAxes, color=CFG.VIS_INK["dark"], clip_on=False,
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
                        fontsize=_JF_FA, color=CFG.VIS_INK["muted"], ha="left", va="bottom",
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
            # The LABEL below already reads CFG.CONF_BAND_HIGH; the shading must start at the same
            # value or the panel states a threshold it does not draw.
            _hc_vis_lo10 = max(float(CFG.CONF_BAND_HIGH), y_lo_f10)
            if _hc_vis_lo10 < _y_top10 - 0.005:
                ax.text(0.01, 0.98, f"High confidence (≥{CFG.CONF_BAND_HIGH:.2f})", color=CFG.VIS_BAND["high"],
                        fontsize=8, ha="left", va="top", fontweight="bold",
                        style="italic", transform=ax.transAxes, zorder=6,
                        bbox=dict(boxstyle="round,pad=0.15", fc=CFG.VIS_TINT["green"], ec=CFG.CONF_BAND_COLOURS["high"],
                                  alpha=0.85, linewidth=0.6))
            # Acceptable zone: 0.80–0.90
            _acc_vis_lo10 = max(0.80, y_lo_f10)
            _acc_vis_hi10 = min(0.90, _y_top10)
            if _acc_vis_hi10 > _acc_vis_lo10 + 0.005:
                _acc_mid10 = (_acc_vis_lo10 + _acc_vis_hi10) / 2
                ax.text(0.01, _acc_mid10, f"Acceptable ({CFG.CONF_BAND_ACCEPTABLE:.2f}–{CFG.CONF_BAND_HIGH:.2f})", color=CFG.VIS_BAND["moderate"],
                        fontsize=8, ha="left", va="center", fontweight="bold",
                        style="italic", transform=_yt10, zorder=6,
                        bbox=dict(boxstyle="round,pad=0.15", fc=CFG.VIS_TINT["cream"], ec=CFG.CONF_BAND_COLOURS["acceptable"],
                                  alpha=0.85, linewidth=0.6))
            # Below threshold zone: y_lo–0.80
            _bel_vis_hi10 = min(0.80, _y_top10)
            if _bel_vis_hi10 > y_lo_f10 + 0.005:
                _bel_mid10 = (y_lo_f10 + _bel_vis_hi10) / 2
                ax.text(0.01, _bel_mid10, f"Below threshold (<{CFG.CONF_BAND_ACCEPTABLE:.2f})", color=CFG.VIS_ACCENT_DEEP["orange"],
                        fontsize=8, ha="left", va="center", fontweight="bold",
                        style="italic", transform=_yt10, zorder=6,
                        bbox=dict(boxstyle="round,pad=0.15", fc=CFG.VIS_TINT["red"], ec=CFG.CONF_BAND_COLOURS["below"],
                                  alpha=0.85, linewidth=0.6))

            ax.set_xticks(range(len(existing_tiers)))
            ax.set_xticklabels(existing_tiers, rotation=40, ha="right", fontsize=9)
            for _tick10x, _tier10x in zip(ax.get_xticklabels(), existing_tiers):
                _tick10x.set_color(TIER_PALETTE.get(_tier10x, "black"))
            ax.tick_params(axis="y", labelsize=9, labelcolor=CFG.VIS_ACCENT["axis_left"])
            ax.yaxis.label.set_color(CFG.VIS_INK["near_black"])
            ax.yaxis.grid(True, color=CFG.VIS_ACCENT["axis_left"], linewidth=0.5, linestyle="-", alpha=0.20, zorder=0)
            ax.set_axisbelow(True)
            ax.set_xlabel("Degrader Tier", )
            ax.set_ylabel("AI Quality Score  (Boltz Confidence / pTM / ipTM;  0–1 scale)", )

            # Legend: single horizontal row, top-right, above data
            from matplotlib.patches import Patch as _Patch10
            conf_patch = _Patch10(facecolor=CFG.VIS_INK["pale"], edgecolor="black", linewidth=0.8,
                                  label="Boltz Confidence  (box)")
            handles_f10, labels_f10 = ax.get_legend_handles_labels()
            _all_h10 = [conf_patch] + handles_f10
            _all_l10 = ["Boltz Confidence  (box)"] + labels_f10
            _leg10 = ax.legend(_all_h10, _all_l10,
                               loc="lower right", bbox_to_anchor=(1.0, 1.01),   # lifted just above the panel
                                fancybox=True,
                               ncol=len(_all_h10))   # all in one row
            _leg10.set_zorder(20)

            plt.tight_layout()
            _stat_box(ax, _kruskal_by_tier(df, CFG.COL_CONF), "lower centre")
            _tier_seps(plt.gca())   # consistent tier separators
            plt.savefig(out_dir / "Figure_05a_AI_Quality_Assessment.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('05')} — {e}")
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
    if "ptm" in df.columns and "iptm" in df.columns and CFG.COL_TIER in df.columns:
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
                                 facecolor=CFG.VIS_BAND["high_fill"], alpha=0.45, zorder=0, linewidth=0))
            # Zone 3: below diagonal (pTM > ipTM)
            _verts_dn = ([(x, x - _band) for x in _diag_pts] +
                         [(_hi18, _lo18), (_lo18, _lo18)])
            ax.add_patch(_Poly18(_verts_dn, closed=True,
                                 facecolor=CFG.VIS_BAND["low_fill"], alpha=0.45, zorder=0, linewidth=0))
            # Zone 2 (diagonal band) is implicitly the white strip between the two coloured zones

            # Identity diagonal line
            ax.plot(_diag_pts, _diag_pts, color=CFG.VIS_INK["muted"], linestyle="--", linewidth=1.3,
                    alpha=0.60, zorder=2)
            # Band boundary dashes
            ax.plot(_diag_pts, _diag_pts + _band, color=CFG.VIS_INK["muted"], linestyle=":",
                    linewidth=0.7, alpha=0.40, zorder=2)
            ax.plot(_diag_pts, _diag_pts - _band, color=CFG.VIS_INK["muted"], linestyle=":",
                    linewidth=0.7, alpha=0.40, zorder=2)

            """
            Zone text — placed BELOW the upper-left legend to avoid overlap.
            Legend with ncol=2 is ~2 rows tall; empirically offset to 0.60 axes-y.
            """
            ax.text(0.04, 0.60, "ipTM > pTM\nInterface well-modelled",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=8.5, color=CFG.VIS_BAND["high"], fontweight="bold", style="italic",
                    bbox=dict(boxstyle="round,pad=0.22", fc=CFG.VIS_BAND["high_fill"], ec=CFG.VIS_BAND["high"],
                              alpha=0.80, linewidth=0.7))
            ax.text(0.97, 0.97, "pTM = ipTM\n(±0.03 band)",
                    transform=ax.transAxes, ha="right", va="top",
                    fontsize=_JF_FA, color=CFG.VIS_INK["muted"], style="italic")
            ax.text(0.96, 0.10, "pTM > ipTM\nFold > interface\n(potential decoy)",
                    transform=ax.transAxes, ha="right", va="bottom",
                    fontsize=8.5, color=CFG.VIS_ACCENT_DEEP["red"], fontweight="bold", style="italic",
                    bbox=dict(boxstyle="round,pad=0.22", fc=CFG.VIS_BAND["low_fill"], ec=CFG.VIS_ACCENT_DEEP["red"],
                              alpha=0.80, linewidth=0.7))

            """
            Rendering strategy:
            Tier_3 tier (n~34k) → hexbin background density (YlOrBr)
            All other tiers   → scatter with tier colour, drawn on top in order
            This keeps tier colours distinguishable while showing the density of the bulk.
            """
            for tier in sorted(existing_tiers, key=lambda t: t == CFG.TIER_TOP):
                sub18 = f18_df[f18_df[CFG.COL_TIER] == tier]
                if sub18.empty:
                    continue
                ptm_vals  = sub18["ptm"].values
                iptm_vals = sub18["iptm"].values
                col18     = TIER_PALETTE.get(tier, CFG.VIS_INK["faint"])
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
                               markerfacecolor=CFG.VIS_INK["ghost"], markeredgecolor="none",
                               markersize=7, alpha=0.55,
                               label="Individual complex  (hexbin/scatter)"))
            _leg18.append(_L18([0], [0], marker="D", color="w",
                               markerfacecolor=CFG.VIS_INK["ghost"], markeredgecolor="black",
                               markersize=9, markeredgewidth=1.0,
                               label="Tier median  (◆ diamond)"))
            for tier in existing_tiers:
                sub18 = f18_df[f18_df[CFG.COL_TIER] == tier]
                if len(sub18) < 2:
                    continue
                mx18, my18 = float(sub18["ptm"].median()), float(sub18["iptm"].median())
                col18 = TIER_PALETTE.get(tier, CFG.VIS_INK["faint"])
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
            existing_tiers18 = [t for t in existing_tiers if t in f18_df[CFG.COL_TIER].values]
            _kde_x18 = np.linspace(_lo18, _hi18, 300)
            for tier in existing_tiers18:
                sub18_m = f18_df[f18_df[CFG.COL_TIER] == tier]
                _min_kde = 2 if tier == CFG.TIER_TOP else 10   # gaussian_kde needs ≥ 2 points
                if len(sub18_m) < _min_kde:
                    continue
                _col18m = TIER_PALETTE.get(tier, CFG.VIS_INK["faint"])
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
            _tt_col_f05 = TIER_PALETTE.get(CFG.TIER_TOP, CFG.VIS_ACCENT["green"])
            _tt_sub_f05 = f18_df[f18_df[CFG.COL_TIER] == CFG.TIER_TOP]
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
            ax_top.set_ylabel("pTM density",  labelpad=2)
            ax_right.set_xlabel("ipTM density",  labelpad=2)
            ax_top.set_axisbelow(True)
            ax_right.set_axisbelow(True)
            ax_top.tick_params(axis="y", labelsize=7.0)
            ax_right.tick_params(axis="x", labelsize=7.0)
            for _sp in list(ax_top.spines.values()) + list(ax_right.spines.values()):
                _sp.set_visible(False)

            ax.set_xlim(_lo18, _hi18)
            ax.set_ylim(_lo18, _hi18)
            ax.set_xlabel("pTM — Global Fold Confidence  (0–1)", )
            ax.set_ylabel("ipTM — Interface Confidence  (0–1)", )
            ax.tick_params(axis="both", labelsize=9)
            ax.grid(True, color=CFG.VIS_INK["grid"], linewidth=0.5, alpha=0.6, zorder=1)
            # Legend: placed in the top KDE marginal (left side) to keep scatter uncluttered
            ax_top.legend(handles=_leg18, loc="upper left", ncol=3,
                           fancybox=True,
                          title="Tier  (◆ = median pTM, ipTM)", 
                          bbox_to_anchor=(0.0, 1.0))
            plt.tight_layout()
            _stat_box(ax, _wilcoxon_ptm_iptm(df), "lower right")
            plt.savefig(out_dir / "Figure_06_pTM_vs_ipTM_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('06')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig



# =============================================================================
# STEP 5/8 — 05_Catalytic_Geometry_and_Mechanism
# =============================================================================
def _fig_folder05_catalytic(df, features, out_dir, reporter, existing_tiers, _pa, _imgs, _tt_has_imgs):
    """Folder 05_Catalytic_Geometry_and_Mechanism — structural + mechanistic geometry."""
    reporter.section("Step 5/8 — 05_Catalytic_Geometry_and_Mechanism · structural + mechanistic geometry")
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
        reporter.log(f"  ! Skipped: {_fig_path('07')} — no RMSD column found. "
                     f"All column names tried: {_rmsd17_aliases}. "
                     f"Partial matches in CSV (containing 'rmsd'/'active_site'): {_rmsd_candidates}. "
                     f"Full column list: {list(df.columns)}")
    if _rmsd17_col and CFG.COL_TIER in df.columns:
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
            valid_tiers_f17 = [t for t in existing_tiers if t in f17_plot[CFG.COL_TIER].values]
            if f17_df.empty or not valid_tiers_f17:
                # No numeric Active_Site_RMSD rows → the % zone maths would divide by zero and seaborn
                # would be asked to plot an empty frame. Skip cleanly via the block's own handler below.
                raise RuntimeError("no rows with a numeric Active_Site_RMSD")
            _sn_max17 = 0.0
            for _t17s in valid_tiers_f17:
                _sv17 = f17_plot.loc[f17_plot[CFG.COL_TIER] == _t17s, "Active_Site_RMSD"]
                if 0 < len(_sv17) <= 30:
                    _sn_max17 = max(_sn_max17, float(_sv17.max()))
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
            sns.violinplot(data=f17_plot, x=CFG.COL_TIER, y="Active_Site_RMSD",
                           order=valid_tiers_f17, palette=TIER_PALETTE,
                           inner="quartile", cut=0, linewidth=1.1, ax=ax, zorder=3)
            # Connected median trend across tiers (markers joined by a line; value
            # labels boxed just below each node) — mirrors the Fig 16 penalty line.
            _med_xs17, _med_ys17, _med_raw17 = [], [], []
            for _im17, _tier17m in enumerate(valid_tiers_f17):
                _med17m = f17_df.loc[f17_df[CFG.COL_TIER] == _tier17m,
                                     "Active_Site_RMSD"].median()
                if np.isnan(_med17m):
                    continue
                _med_xs17.append(_im17)
                _med_ys17.append(min(float(_med17m), _violin_cap))
                _med_raw17.append(float(_med17m))
            if _med_xs17:
                ax.plot(_med_xs17, _med_ys17, color=CFG.VIS_ACCENT["alert"], linewidth=2.2,
                        marker="D", markersize=9, markeredgecolor="black",
                        markeredgewidth=0.9, zorder=8, label="Median RMSD (per tier)")
                for _xi17, _yi17, _rv17 in zip(_med_xs17, _med_ys17, _med_raw17):
                    ax.text(_xi17, _yi17 - 0.13, f"{_rv17:.2f}", ha="center", va="top",
                            fontsize=_JF_FA, color=CFG.VIS_ACCENT["alert"], fontweight="bold", zorder=9,
                            bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                      ec=CFG.VIS_ACCENT["alert"], alpha=0.9, linewidth=0.6))
            # Very sparse strip for large-n tiers only
            _large_n_tiers17 = [t for t in valid_tiers_f17
                                 if len(f17_plot[f17_plot[CFG.COL_TIER] == t]) > 30]
            if _large_n_tiers17:
                _f17_large = f17_plot[f17_plot[CFG.COL_TIER].isin(_large_n_tiers17)]
                sns.stripplot(data=_f17_large, x=CFG.COL_TIER, y="Active_Site_RMSD",
                              order=valid_tiers_f17, color="black", alpha=0.04,
                              size=1.5, jitter=True, ax=ax, zorder=2)
            """
            For tiers with very few data points (≤30), draw individual dots prominently
            so every value is visible — avoids them being lost inside the violin body
            """
            for _i17, _tier17 in enumerate(valid_tiers_f17):
                _sub17 = f17_plot.loc[f17_plot[CFG.COL_TIER] == _tier17, "Active_Site_RMSD"].dropna()
                if len(_sub17) <= 30:
                    _col17 = TIER_PALETTE.get(_tier17, CFG.VIS_INK["faint"])
                    _jit17 = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(_sub17))
                    ax.scatter(_i17 + _jit17, _sub17.values,
                               color=_col17, edgecolors="black", s=55,
                               linewidths=0.8, alpha=0.88, zorder=7)
            # Outlier scatter: actual RMSD values above 3 Å (unclipped)
            for _ti, tier in enumerate(valid_tiers_f17):
                sub_f6 = f17_df[f17_df[CFG.COL_TIER] == tier]
                _out = sub_f6[sub_f6["Active_Site_RMSD"] > _violin_cap]["Active_Site_RMSD"].values
                if len(_out) == 0:
                    continue
                _jit_f6 = np.random.default_rng(17 + _ti).uniform(-0.22, 0.22, size=len(_out))
                ax.scatter(_ti + _jit_f6, _out,
                           color=TIER_PALETTE.get(tier, CFG.VIS_INK["faint"]), s=38, alpha=0.65,
                           marker="o", edgecolors="white", linewidths=0.4, zorder=4)
                if len(_out) >= 1:
                    _lbl_col = TIER_PALETTE.get(tier, CFG.VIS_INK["faint"])
                    ax.text(_ti, y_ceil + 0.04, f"+{len(_out)}",
                            ha="center", va="bottom", fontsize=9.5,
                            color=_lbl_col, fontweight="bold",
                            zorder=30, clip_on=False,
                            bbox=dict(boxstyle="round,pad=0.18", fc="white",
                                      ec=_lbl_col, alpha=1.0, linewidth=1.0))

            # Capped axis, fixed ticks, and dashed break marker at 3 Å
            ax.axhline(3.0, color=CFG.VIS_INK["ghost"], linestyle=":", linewidth=0.9, alpha=0.7)
            ax.set_ylim(0, y_ceil)
            _y_ticks_f6 = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
            _y_extra_f6 = [t for t in np.arange(3.5, y_ceil, 0.5) if t <= y_ceil - 0.1]
            _y_ticks_f6 += _y_extra_f6
            ax.set_yticks(_y_ticks_f6)
            ax.set_yticklabels([f"{t:.1f}" for t in _y_ticks_f6])
            ax.tick_params(axis="y", labelsize=7)
            ax.set_xlabel("Degrader Tier", )
            ax.set_ylabel("Active Site RMSD (Å)  — lower = better", )
            ax.text(0.99, 0.99, " ",
                    ha="right", va="top", fontsize=7, color=CFG.VIS_INK["ghost"], style="italic",
                    transform=ax.transAxes)
            # Build two-line x-tick labels: tier name + n / median stats (using unclipped values)
            _xtick_labels17 = []
            for _tier17 in valid_tiers_f17:
                _sub17 = f17_df.loc[f17_df[CFG.COL_TIER] == _tier17,
                                    "Active_Site_RMSD"].dropna()
                if len(_sub17) > 0:
                    _med17 = float(_sub17.median())
                    _xtick_labels17.append(
                        f"{_tier17}\nn={len(_sub17):,}  med={_med17:.2f}Å"
                    )
                else:
                    _xtick_labels17.append(_tier17)
            ax.set_xticks(range(len(valid_tiers_f17)))
            ax.set_xticklabels(_xtick_labels17, rotation=40, ha="right", fontsize=_JF_FA)
            for _tick17, _tier17 in zip(ax.get_xticklabels(), valid_tiers_f17):
                _tick17.set_color(TIER_PALETTE.get(_tier17, "black"))

            # Per-zone count labels (90°-rotated) inside each RMSD zone band
            _re, _ra, _rd = CFG.RMSD_BAND_EXCELLENT, CFG.RMSD_BAND_ACCEPTABLE, CFG.RMSD_BAND_DIVERGED
            _cbc = CFG.CONF_BAND_COLOURS
            _zone_count_defs = [
                (0.0, _re,    _cbc["high"],       f"<{_re:.0f} Å"),
                (_re, _ra,    _cbc["acceptable"], f"{_re:.0f}–{_ra:.0f} Å"),
                (_ra, _rd,    _cbc["below"],      f"{_ra:.0f}–{_rd:.0f} Å"),
                (_rd, np.inf, CFG.VIS_ACCENT_DEEP["red_dark"],          f">{_rd:.0f} Å"),
            ]
            for _zlo, _zhi, _zcol, _zlbl in _zone_count_defs:
                _zmask = (f17_df["Active_Site_RMSD"] >= _zlo) & (f17_df["Active_Site_RMSD"] < _zhi)
                _zn = int(_zmask.sum())
                _zy = min(_zhi, y_ceil) * 0.98
                ax.text(len(valid_tiers_f17) - 0.3, _zy,
                        f"n={_zn:,}  {_zlbl}",
                        ha="center", va="top", fontsize=7.0, color=_zcol, fontweight="bold",
                        rotation=90, zorder=25, clip_on=False,
                        bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=_zcol,
                                  alpha=1.0, linewidth=0.5))

            # Zone band legend (top-left) for the shaded quality bands
            from matplotlib.patches import Patch as _P6z
            _zone_legend = [
                _P6z(facecolor=CFG.VIS_ACCENT_DEEP["mint_fill"], alpha=0.55, edgecolor="none",
                     label=f"Excellent  (< {CFG.RMSD_BAND_EXCELLENT:g} Å)  {pct_green:.0f}%"),
                _P6z(facecolor=CFG.VIS_ACCENT_DEEP["cream_fill"], alpha=0.55, edgecolor="none",
                     label=f"Acceptable  ({CFG.RMSD_BAND_EXCELLENT:g}–{CFG.RMSD_BAND_ACCEPTABLE:g} Å)  {pct_yellow:.0f}%"),
                _P6z(facecolor=CFG.VIS_ACCENT_DEEP["rose_fill"], alpha=0.55, edgecolor="none",
                     label=f"Diverged  (> {CFG.RMSD_BAND_ACCEPTABLE:g} Å)  {pct_red:.0f}%"),
            ]
            _tier_handles, _tier_labels = ax.get_legend_handles_labels()
            _combined_handles = _zone_legend + _tier_handles
            _combined_labels  = [h.get_label() for h in _zone_legend] + _tier_labels
            _mdh17 = _md_ready_stars_cat(ax, f17_df, valid_tiers_f17, "Active_Site_RMSD")
            for _h17 in _mdh17:
                _combined_handles = _combined_handles + [_h17]
                _combined_labels  = _combined_labels + [_h17.get_label()]
            ax.legend(_combined_handles, _combined_labels,
                      loc="upper left")
            plt.tight_layout()
            _tier_seps(plt.gca())   # consistent vertical tier separators
            plt.savefig(out_dir / "Figure_07_ActiveSite_RMSD_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('07')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 08: Feature Correlation Matrix (Spearman ρ) ---
    _corr_col_labels = {
        CFG.COL_CONF:      "Conf",
        "iptm":                        "ipTM",
        "Binding_Probability":         "BindP",
        "Interaction_Density_Norm":    "IntDen",
        "Chemical_Affinity_Score":     "ChemAff",
        "tier_numeric":                "Tier",
        "Pareto_Rank":                 "Pareto",
        CFG.COL_SN2:            "SN2°",
        CFG.COL_ID_PCT:                "SeqID%",
        "SN2_Trajectory_Deviation_A":  "TrajDev",
        "soft_catalytic_score":        "SoftCat",
        # Mechanistic / chemistry drivers — real matrix cells: their actual pairwise Spearman with
        # every other feature.
        "scissile_cf_bde":             "C-F BDE",
        "sn2_backside_occlusion":      "Backside",
        "halide_stabilisation_score":  "Halide",
        "competence_score":            "Compet.",
        "carboxylate_clamp_integrity": "Cbx clamp",
        "mechanistic_score_effective": "Mech.eff",
    }
    # Map tier_numeric alias — try every plausible name written by different pipeline versions
    for _alias in ["tier_val", "tier_value", "Tier_Score", "tier_score",
                   "Tier_Numeric", "degrader_tier_numeric", "Degrader_Tier_Numeric"]:
        if _alias in df.columns and "tier_numeric" not in df.columns:
            df["tier_numeric"] = pd.to_numeric(df[_alias], errors="coerce")
    # Last-resort: derive numeric rank from the categorical degrader_tier string
    if "tier_numeric" not in df.columns and CFG.COL_TIER in df.columns:
        _tier_rank_map = {t: i for i, t in enumerate(TIER_ORDER_LOGIC)}
        df["tier_numeric"] = df[CFG.COL_TIER].map(_tier_rank_map)
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

        """
        This is a SECOND, INDEPENDENT Benjamini-Hochberg family — the feature-correlation matrix — and it
        is not the family written to 06_Statistical_Tests.csv (which corrects the panel tests and the
        statistical battery together). Two families are legitimate: these are different questions asked of
        different hypotheses, and pooling them would inflate m and cost power on both.

        But a Spearman that appears in both is corrected against a different m in each, so its q here and
        its q in the CSV WILL NOT MATCH. That is not an error; it is a fact a reader must be told, and it
        is stated on the panel rather than left to be discovered.
        """
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
        # Benjamini-Hochberg from the library rather than a hand-rolled step-up: same correction,
        # one fewer place for an off-by-one in the rank index to hide.
        from scipy.stats import false_discovery_control
        _adj7 = np.clip(false_discovery_control(np.asarray(_pv_arr, dtype=float), method="bh"), 0, 1)
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
            ("Chemistry / mechanism", ["C-F BDE", "Backside", "Halide", "Compet.", "Cbx clamp", "Mech.eff"]),
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
        _, ax_f3 = plt.subplots(figsize=(14, 12.5))
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
                   fontsize=8, color=CFG.VIS_INK["muted"], style="italic",
                   bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=CFG.VIS_INK["palest"],
                             alpha=0.88, linewidth=0.7))
        # Labels are short (≤7 chars) — no truncation needed
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
                    _txt_art.set_color("white" if abs(_rho_val) > 0.5 else CFG.VIS_INK["ink_pure"])
                except ValueError:
                    pass

        # ── Self-explanatory overlays ─────────────────────────────────────────
        from matplotlib.patches import Rectangle as _Rect7
        _N7 = len(_display_names)
        _fam_lookup_f7 = {l: g for g, ls in _fam_groups_f7 for l in ls}
        _fams_seq_f7   = [_fam_lookup_f7.get(l, "Other") for l in _display_names]

        # Colour the feature tick labels (and the family brackets below) by family group, so the
        # four families read at a glance rather than as one undifferentiated black axis.
        _FAM_COL7 = {"Outcome": CFG.VIS_ACCENT["green"], "Catalytic geometry": CFG.VIS_ACCENT["blue"],
                     "Chemistry / mechanism": CFG.VIS_ACCENT["vermillion"],
                     "Affinity / seq / rank": CFG.VIS_ACCENT["amber"],
                     "Binding & confidence": CFG.VIS_ACCENT["magenta"]}
        for _tl7, _fam7 in zip(ax_f3.get_xticklabels(), _fams_seq_f7):
            _tl7.set_color(_FAM_COL7.get(_fam7, CFG.VIS_INK["near_black"]))
        for _tl7, _fam7 in zip(ax_f3.get_yticklabels(), _fams_seq_f7):
            _tl7.set_color(_FAM_COL7.get(_fam7, CFG.VIS_INK["near_black"]))

        """
        1) Family separator lines at each group boundary — clipped to the lower
        (data) triangle so they never run through the empty upper-right area
        that holds the key + take-away note.
        """
        for _i7 in range(1, _N7):
            if _fams_seq_f7[_i7] != _fams_seq_f7[_i7 - 1]:
                ax_f3.plot([0, _i7], [_i7, _i7], color=CFG.VIS_INK["outline"], linewidth=1.7, zorder=5)
                ax_f3.plot([_i7, _i7], [_i7, _N7], color=CFG.VIS_INK["outline"], linewidth=1.7, zorder=5)

        # 2) Family bracket labels under the matrix (one per group).
        for _gname, _glbls in _fam_groups_f7:
            _gi = [i for i, l in enumerate(_display_names) if l in _glbls]
            if not _gi:
                continue
            ax_f3.annotate(_gname, xy=((min(_gi) + max(_gi) + 1) / 2.0, _N7 + 1.15),
                           xycoords="data", ha="center", va="top", fontsize=8.5,
                           fontweight="bold", color=_FAM_COL7.get(_gname, CFG.VIS_INK["dark"]),
                           annotation_clip=False)

        """
        3) Spotlight the Tier column — every feature's correlation with the
        degradation outcome is the question the figure answers.
        """
        if "Tier" in _display_names:
            _ti7 = _display_names.index("Tier")
            ax_f3.add_patch(_Rect7((_ti7, _ti7), 1, _N7 - _ti7, fill=False,
                            edgecolor=CFG.VIS_INK["near_black"], linewidth=2.4, zorder=6))


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
        CFG.COL_CONF: "AI Confidence", CFG.COL_MECH_S: "Mech. Score",
        "Binding_Probability": "Binding Prob.", CFG.COL_ID_PCT: "Seq. Identity (%)",
        CFG.COL_SN2: "SN2 Angle (°)", "Active_Site_RMSD": "RMSD (Å, inv.)",
        "SN2_Trajectory_Deviation_A": "SN2 Traj. Dev. (Å)",
        "soft_catalytic_score": "Soft Catalytic Score",
    }
    _f08_invert = {"Active_Site_RMSD", "SN2_Trajectory_Deviation_A"}
    _f08_cols = [c for c in _f08_metric_map if c in df.columns and CFG.COL_TIER in df.columns]
    if len(_f08_cols) >= 2:
        try:
            _f08_rows = []
            for tier in existing_tiers:
                sub = df[df[CFG.COL_TIER] == tier]
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
                    _means = np.mean(_samples, axis=1)          # the panel plots tier MEANS
                    return float(np.percentile(_means, 2.5)), float(np.percentile(_means, 97.5))

                # Precompute per-tier bootstrap CIs for each metric (raw scale)
                _f08_ci = {}  # (tier, metric) -> (lo_norm, hi_norm)
                for col in _f08_cols:
                    disp = _f08_metric_map[col]
                    mn_raw = float(_f08_df[disp].min())
                    mx_raw = float(_f08_df[disp].max())
                    _inv = col in _f08_invert
                    for tier in tiers_f08:
                        sub_ci = df[df[CFG.COL_TIER] == tier]
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
                    ax08.xaxis.grid(True, color=CFG.VIS_INK["grid"], linewidth=0.7, linestyle="-", zorder=0)
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
                            _dot_data08.append((nv, rv, TIER_PALETTE.get(tier, CFG.VIS_INK["faint"]), tier))
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
                                ax08.hlines(_base_y, 0.0, 1.0, colors=CFG.VIS_ACCENT_DEEP["gold"],
                                            linewidth=1.0, linestyle="--", alpha=0.55, zorder=1)
                            else:
                                ax08.hlines(_base_y, float(np.nanmin(norm_vals)),
                                            float(np.nanmax(norm_vals)),
                                            colors=CFG.VIS_INK["palest"], linewidth=2.5, zorder=1)

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
                                          ha="center", va="bottom", fontsize=7.0,
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
                    )

                from matplotlib.lines import Line2D as _L08
                _leg08 = [_L08([0],[0], marker="o", color="w",
                               markerfacecolor=TIER_PALETTE.get(t,CFG.VIS_INK["faint"]),
                               markeredgecolor="black", markeredgewidth=0.8,
                               markersize=10, label=t) for t in tiers_f08]
                _leg08.append(_L08([0],[0], color="none", linewidth=0,
                                   label="★ RMSD metrics: inverted — lower RMSD = better score"))
                _leg08.append(_L08([0],[0], color=CFG.VIS_INK["faint"], linewidth=5.5, alpha=0.3,
                                   label="95% bootstrap CI on mean (1 000 resamples)"))
                if _any_tied08:
                    _leg08.append(_L08([0],[0], color=CFG.VIS_ACCENT_DEEP["gold"], linewidth=1.5,
                                       linestyle="--",
                                       label="⚡ Amber dashed = metric tied across all tiers"))
                    for _tm_name, _tm_val in _tied_metrics08:
                        _leg08.append(_L08([0],[0], color="none", linewidth=0,
                                           label=f"   {_tm_name}: all tiers = {_tm_val:.2f}"))
                axes_f08[0].legend(handles=_leg08, loc="upper left",
                                   bbox_to_anchor=(0.0, 1.0),
                                   ncol=max(1, (len(_leg08) + 1) // 2),
                                     fancybox=True)
                plt.tight_layout()
                plt.savefig(out_dir / "Figure_09_Tier_Quality_DotPlot.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
                plt.close()
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('09')} — {e}")
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
        reporter.log(f"  ! Skipped: {_fig_path('10')} — columns not found: {_mech_missing}. "
                     f"Partial matches in CSV: {_mech_candidates}. "
                     f"Full column list: {list(df.columns)}")
    if "Halide_Stabilisation" in df.columns and "Carboxylate_Clamp" in df.columns and CFG.COL_TIER in df.columns:
        try:
            mech_df = df.copy()
            mech_df["HS"] = mech_df["Halide_Stabilisation"].map(
                lambda x: "Stabilised" if str(x).strip().lower() in ("true", "1", "yes") else "Unstabilised")
            mech_df["CC"] = mech_df["Carboxylate_Clamp"].map(
                lambda x: "Clamped" if str(x).strip().lower() in ("true", "1", "yes") else "Unclamped")
            mech_df["Mech_State"] = mech_df["HS"] + "\n" + mech_df["CC"]

            state_order = ["Stabilised\nClamped", "Stabilised\nUnclamped",
                           "Unstabilised\nClamped", "Unstabilised\nUnclamped"]
            cross = pd.crosstab(mech_df[CFG.COL_TIER], mech_df["Mech_State"])
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
                _register_p("Chi-square — interaction-density band × tier", "19_Interaction_Density_by_Tier",
                            float(_chi2_09), int(_cross_raw.values.sum()), float(_p_09), dof=int(_dof_09))
            except Exception:
                _chi2_stat = ""

            _, ax19 = plt.subplots(figsize=(11, 7))
            sns.heatmap(cross_pct, annot=annot_labels, fmt="", cmap="YlOrRd", linewidths=0.5,
                        annot_kws={"size": 10, "weight": "bold"},
                        cbar_kws={"label": "% of tier"}, ax=ax19)
            if _chi2_stat:
                ax19.text(0.99, 0.02, _chi2_stat, transform=ax19.transAxes,
                          ha="right", va="bottom", fontsize=8.5, style="italic",
                          color=CFG.VIS_INK["dark"],
                          bbox=dict(boxstyle="round,pad=0.25", fc="white",
                                    ec=CFG.VIS_INK["palest"], alpha=0.95, linewidth=0.6))
            ax19.set_xlabel("Mechanistic State", )
            ax19.set_ylabel("Degrader Tier", )
            ax19.set_xticklabels(ax19.get_xticklabels(), rotation=15, ha="right", fontsize=9)
            # Colour y-tick labels by tier (same palette as all other figures)
            for tick in ax19.get_yticklabels():
                tier_name = tick.get_text()
                tick.set_color(TIER_PALETTE.get(tier_name, "black"))
                tick.set_rotation(0)
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_10_Mech_State_CrossTab.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('10')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 11: Mechanistic Score — Mean±CI dot plot ---
    if CFG.COL_MECH_S in df.columns and CFG.COL_TIER in df.columns:
        from scipy import stats as _scipy_stats
        f13_data = df.dropna(subset=[CFG.COL_MECH_S])
        fig, ax = plt.subplots(figsize=(13, 8))
        _mfs, _mfm = CFG.MECH_FP_BAND_STRONG, CFG.MECH_FP_BAND_MODERATE   # CFG single source
        _cbc = CFG.CONF_BAND_COLOURS
        ax.axhspan(_mfs, 1.01, alpha=0.10, color=_cbc["high"], zorder=0, label="_nolegend_")
        ax.axhspan(_mfm, _mfs, alpha=0.09, color=_cbc["acceptable"], zorder=0, label="_nolegend_")
        ax.axhspan(0.00, _mfm, alpha=0.08, color=_cbc["below"], zorder=0, label="_nolegend_")
        ax.text(0.99, 0.942, f"Strong zone  (≥{_mfs:.2f})", color=CFG.VIS_BAND["high"],
                fontsize=8, ha="right", va="center", fontweight="bold",
                style="italic", transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.15", fc=CFG.VIS_TINT["green"], ec=CFG.CONF_BAND_COLOURS["high"],
                          alpha=0.85, linewidth=0.6))
        ax.text(0.99, 0.798, f"Moderate zone  ({_mfm:.2f}–{_mfs:.2f})", color=CFG.VIS_BAND["moderate"],
                fontsize=8, ha="right", va="center", fontweight="bold",
                style="italic", transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.15", fc=CFG.VIS_TINT["cream"], ec=CFG.CONF_BAND_COLOURS["acceptable"],
                          alpha=0.85, linewidth=0.6))
        ax.text(0.99, 0.435, f"Weak zone  (<{_mfm:.2f})", color=CFG.VIS_ACCENT_DEEP["orange"],
                fontsize=8, ha="right", va="center", fontweight="bold",
                style="italic", transform=ax.transAxes,
                bbox=dict(boxstyle="round,pad=0.15", fc=CFG.VIS_TINT["red"], ec=CFG.CONF_BAND_COLOURS["below"],
                          alpha=0.85, linewidth=0.6))
        ax.axhline(y=_mfs, color=_cbc["high"], linestyle="--", linewidth=1.2, alpha=0.7)
        ax.axhline(y=_mfm, color=_cbc["below"], linestyle=":", linewidth=1.0, alpha=0.7)
        sns.stripplot(data=f13_data, x=CFG.COL_TIER, y=CFG.COL_MECH_S,
                      order=existing_tiers, palette=TIER_PALETTE, alpha=0.15, size=3.0,
                      jitter=0.28, ax=ax, zorder=1)
        _f10_stats = {}
        for _t10 in existing_tiers:
            _v10 = f13_data.loc[f13_data[CFG.COL_TIER] == _t10,
                                 CFG.COL_MECH_S].dropna()
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
            ax.scatter([i], [mean_v], color=TIER_PALETTE.get(tier, CFG.VIS_INK["faint"]), s=180,
                       zorder=5, edgecolors="black", linewidths=1.2)
            pct_strong   = (vals >= _mfs).sum() / len(vals) * 100
            pct_moderate = ((vals >= _mfm) & (vals < _mfs)).sum() / len(vals) * 100
            pct_weak     = (vals < _mfm).sum() / len(vals) * 100
            if pct_strong > 0:
                ax.text(i, (_mfs + 1.0) / 2, f"{pct_strong:.0f}%", ha="center", va="center",
                        fontsize=7, color=CFG.VIS_BAND["high"], fontweight="bold", zorder=10,
                        bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["high"],
                                  alpha=0.75, linewidth=0.5))
            if pct_moderate > 0:
                ax.text(i, (_mfm + _mfs) / 2, f"{pct_moderate:.0f}%", ha="center", va="center",
                        fontsize=7, color=CFG.VIS_BAND["moderate"], fontweight="bold", zorder=10,
                        bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["acceptable"],
                                  alpha=0.75, linewidth=0.5))
            if pct_weak > 0:
                ax.text(i, _mfm / 2, f"{pct_weak:.0f}%", ha="center", va="center",
                        fontsize=7, color=CFG.VIS_ACCENT_DEEP["orange"], fontweight="bold", zorder=10,
                        bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["below"],
                                  alpha=0.75, linewidth=0.5))
        ax.set_ylim(-0.02, 1.06)
        ax.set_yticks(sorted({0.0, 0.2, 0.4, 0.6, 0.8, 1.0, round(_mfm, 2), round(_mfs, 2)}))
        ax.set_xlabel("Degrader Tier", )
        ax.set_ylabel("Mechanistic Score  (anchor set + graded SN2 angle; config §5.1)", )
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
        ax.set_xticklabels(_xtlbl10, rotation=40, ha="right", fontsize=_JF_FA)
        for _tick13x, _tier13x in zip(ax.get_xticklabels(), existing_tiers):
            _tick13x.set_color(TIER_PALETTE.get(_tier13x, "black"))
        ax.yaxis.grid(True, color=CFG.VIS_INK["tick"], linewidth=0.6, alpha=0.70, zorder=0)
        ax.set_axisbelow(True)
        _means10 = [_f10_stats[t]["mean"] for t in existing_tiers if t in _f10_stats]
        if _means10 and (max(_means10) - min(_means10)) < 0.05:
            _mean10_all = float(np.mean(_means10))
            ax.axhline(y=_mean10_all, color=CFG.VIS_ACCENT_DEEP["gold"], linewidth=1.5,
                       linestyle="--", alpha=0.65, zorder=2)
            ax.text(0.99, 0.52,
                    f"⚡ All tier means within {max(_means10)-min(_means10):.3f} "
                    f"— minimal inter-tier separation\n"
                    f"(mech-score cluster near {_mean10_all:.3f}; see also Fig 08 tied-metric row)",
                    transform=ax.transAxes, ha="right", va="center", fontsize=8,
                    color=CFG.VIS_ACCENT_DEEP["gold_deep"], fontweight="bold", style="italic",
                    bbox=dict(boxstyle="round,pad=0.15", fc=CFG.VIS_TINT["amber"], ec=CFG.VIS_ACCENT_DEEP["gold"],
                              alpha=0.90, linewidth=0.8), zorder=8)
        from matplotlib.lines import Line2D as _L10
        _lh10 = [
            _L10([0],[0], marker="o", color="none", markerfacecolor=CFG.VIS_INK["ghost"],
                 markeredgecolor="none", markersize=5, alpha=0.35,
                 label="Individual complex"),
            _L10([0],[0], marker="D", color="none", markerfacecolor=CFG.VIS_INK["ghost"],
                 markeredgecolor="black", markersize=8, markeredgewidth=1.2,
                 label="Tier mean (◆)"),
            _L10([0],[0], color="black", linewidth=2.0, label="95% confidence interval"),
        ]
        _mdh10 = _md_ready_stars_cat(ax, df, existing_tiers, CFG.COL_MECH_S)
        _lh10.extend(_mdh10)
        ax.legend(handles=_lh10, loc="lower left",
                   fancybox=True,
                  ncol=4).set_zorder(20)
        plt.tight_layout()
        _tier_seps(plt.gca())   # consistent vertical tier separators
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
        ("Mechanistic\nscore",         CFG.COL_MECH_S,        "score"),
    ]
    _fp_avail = [(lbl, col, kind) for (lbl, col, kind) in _fp_feats if col in df.columns]
    if _fp_avail and CFG.COL_TIER in df.columns:
        fig = None
        try:
            _fp_tiers = [t for t in existing_tiers if (df[CFG.COL_TIER] == t).any()]
            _fp_rows, _fp_n = [], []
            for _ft in _fp_tiers:
                _fp_sub = df[df[CFG.COL_TIER] == _ft]
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
                ax.plot(_fp_ring_ang, [_rl] * 300, color=CFG.VIS_BAND["low"], linewidth=0.55,
                        alpha=0.45, zorder=1, linestyle="-", solid_capstyle="round")
            for _ti, _t in enumerate(_fp_tiers):
                _yv = _fpm.loc[_t].values.tolist()
                _vals = _yv + [_yv[0]]
                _c = TIER_PALETTE.get(_t, CFG.VIS_INK["faint"])
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
                               color=CFG.VIS_BAND["low"], size=7.5, fontweight="bold")
            # Radial tick labels parked in an empty wedge between spokes (no data line there).
            ax.set_rlabel_position(np.degrees(np.pi / float(_N_fp)))
            ax.set_ylim(0, 1.10)
            ax.yaxis.grid(False)
            ax.xaxis.grid(True, color=CFG.VIS_INK["palest"], linewidth=0.6, alpha=0.7)
            _leg11b = ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.18),
                                ncol=min(3, len(_fp_tiers)),
                                 fancybox=True)
            _leg11b.set_zorder(20)
            for _lt, _t in zip(_leg11b.get_texts(), _fp_tiers):
                _lt.set_color(TIER_PALETTE.get(_t, "black"))
            plt.tight_layout(rect=[0, 0.15, 1, 1])
            plt.savefig(out_dir / "Figure_11b_Mechanistic_Fingerprint.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('11b')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig
        finally:
            if fig is not None:
                plt.close(fig)

    # --- Figure 12: SN2 Attack Angle — ECDF by Tier ---
    if CFG.COL_SN2 in df.columns and CFG.COL_TIER in df.columns:
        plot_df = df[df[CFG.COL_SN2] > 0].copy()
        valid_tiers_f14 = [t for t in existing_tiers if t in plot_df[CFG.COL_TIER].values]
        fig, ax = plt.subplots(figsize=(12, 7))
        # Shade tier quality zones + count labels inside each zone band.
        # Zone lower bounds derive from CFG.TIER_ANGLE_MIN (single source); the upper
        # bound of each zone is the next-higher tier's minimum (top zone capped at 180°).
        _zt11   = [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]
        _zcol11 = [CFG.VIS_ACCENT["green"], CFG.VIS_ACCENT["sky"], CFG.VIS_ACCENT["blue"], CFG.VIS_ACCENT["magenta"]]
        _zlo11  = [CFG.TIER_ANGLE_MIN[t] for t in _zt11]
        _zhi11  = [180.0] + _zlo11[:-1]
        _zone_defs11 = [(_zlo11[i], _zhi11[i], _zcol11[i]) for i in range(len(_zt11))]
        for lo, hi, col in _zone_defs11:
            ax.axvspan(lo, hi, alpha=0.07, color=col)
            ax.text((lo + hi) / 2, 1.025, f"≥{lo}°", ha="center", va="bottom",
                    fontsize=_JF_FA, color=col, fontweight="bold", transform=ax.get_xaxis_transform())
            _zn11 = int(((plot_df[CFG.COL_SN2] >= lo) & (plot_df[CFG.COL_SN2] < hi)).sum())
            ax.text((lo + hi) / 2, 0.10, f"n={_zn11:,}",
                    ha="center", va="center", fontsize=7.0, color=col, fontweight="bold",
                    rotation=90, zorder=10,
                    transform=ax.get_xaxis_transform(),
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.68))
        # ECDF per tier — median value labels are collected here and placed after the
        # loop in alternating vertical bands so clustered medians never overlap.
        _med_pts = []
        for tier in valid_tiers_f14:
            tier_angles = np.sort(plot_df.loc[plot_df[CFG.COL_TIER] == tier, CFG.COL_SN2].values)
            ecdf_y = np.arange(1, len(tier_angles) + 1) / len(tier_angles)
            col = TIER_PALETTE.get(tier, CFG.VIS_INK["faint"])
            ax.step(tier_angles, ecdf_y, where="post", color=col, linewidth=2.5,
                    label=f"{tier}  (n={len(tier_angles):,})")
            _n_ecdf = len(tier_angles)
            _eps_dkw = np.sqrt(np.log(2 / 0.05) / (2 * _n_ecdf))
            _ecdf_y_lo = np.clip(ecdf_y - _eps_dkw, 0, 1)
            _ecdf_y_hi = np.clip(ecdf_y + _eps_dkw, 0, 1)
            ax.fill_between(tier_angles, _ecdf_y_lo, _ecdf_y_hi,
                            color=TIER_PALETTE.get(tier, CFG.VIS_INK["faint"]), alpha=0.15, zorder=1)
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
        # MD-ready overlay — the top-performing MD-selected complexes (incl. those inside
        # Tier_1A) marked as gold stars on their own tier's ECDF curve.
        _mdsel14 = plot_df[plot_df.get("MD_Selected", pd.Series(False, index=plot_df.index))
                           .astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes"])]
        for _, _mr in _mdsel14.iterrows():
            _ma = float(_mr[CFG.COL_SN2]); _mt = _mr[CFG.COL_TIER]
            _ta = np.sort(plot_df.loc[plot_df[CFG.COL_TIER] == _mt, CFG.COL_SN2].values)
            if len(_ta) == 0:
                continue
            _myv = float(np.interp(_ma, _ta, np.arange(1, len(_ta) + 1) / len(_ta)))
            ax.scatter([_ma], [_myv], marker="*", s=300, facecolor=CFG.VIS_ACCENT["star"],
                       edgecolor="black", linewidths=1.1, zorder=9)
        for _, _cr14 in _control_star_df(df).iterrows():
            _ca14 = float(_cr14[CFG.COL_SN2]); _ct14 = _cr14[CFG.COL_TIER]
            _cta14 = np.sort(plot_df.loc[plot_df[CFG.COL_TIER] == _ct14, CFG.COL_SN2].values)
            if len(_cta14) == 0:
                continue
            _cyv14 = float(np.interp(_ca14, _cta14, np.arange(1, len(_cta14) + 1) / len(_cta14)))
            _ctrl_star(ax, [_ca14], [_cyv14])
        if len(_mdsel14):
            ax.scatter([], [], marker="*", s=180, facecolor=CFG.VIS_ACCENT["star"], edgecolor="black",
                       linewidths=1.0, label=f"MD-ready hits (n={len(_mdsel14)})")
        if not _control_star_df(df).empty:
            ax.scatter([], [], marker="*", s=200, facecolor=CFG.VIS_ACCENT["control"],
                       edgecolor=CFG.VIS_ACCENT["control_edge"], linewidths=1.4, label="3R3U × FA (control)")

        # Vertical threshold lines
        # Threshold lines sourced from CFG so they always match the tier gates
        # (e.g. Tier_1A = 170°, not a stale hardcoded literal).
        _ecdf_tiers = [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]
        _ecdf_angles = [CFG.TIER_ANGLE_MIN[t] for t in _ecdf_tiers]
        for angle, tier_key in zip(_ecdf_angles, _ecdf_tiers):
            ax.axvline(x=angle, color=TIER_PALETTE.get(tier_key, CFG.VIS_INK["faint"]), linestyle="--", alpha=0.65, linewidth=1.0)
        ax.set_xlim(30, 185)
        ax.set_ylim(0, 1.09)
        ax.set_xticks([30, 45, 90, 130] + sorted(set(_ecdf_angles)) + [180])
        ax.set_yticks([i/10 for i in range(0, 11)])
        ax.set_yticklabels([f"{i*10}%" for i in range(0, 11)])
        ax.set_xlabel("SN2 Attack Angle (°)  — 180° = ideal linear nucleophilic back-attack", )
        ax.set_ylabel("Cumulative Fraction of Complexes in Tier", )
        # upper left — ECDF lines fan right so top-left is always clear
        from matplotlib.lines import Line2D as _L11
        _hdr11 = _L11([0],[0], color="none", linewidth=0, label="Tier (●=median)")
        _h11, _l11 = ax.get_legend_handles_labels()
        _leg14 = ax.legend(handles=[_hdr11] + _h11, labels=["Tier (●=median)"] + _l11,
                           loc="upper left", bbox_to_anchor=(0.01, 0.99),
                             fancybox=True,
                           ncol=max(2, (len(_h11) + 1) // 2))
        _leg14.set_zorder(20)
        # Thin horizontal grid lines at every 10% ECDF level for easy reading
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color=CFG.VIS_INK["tick"], linewidth=0.6, linestyle="-", alpha=0.85, zorder=0)
        ax.xaxis.grid(True, color=CFG.VIS_INK["grid"], linewidth=0.5, linestyle=":", alpha=0.7, zorder=0)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_12_SN2_Angle_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # --- Figure 13a: SN2 Mechanism Geometry Scatter ---
    if "Dist_Nucleophile" in df.columns and CFG.COL_SN2 in df.columns:
        _n_raw9 = len(df)
        plot_df12 = df[
            (df["Dist_Nucleophile"] > 0) &
            (df["Dist_Nucleophile"] < 5.0) &   # plot display window (Å), not a scientific gate
            (df[CFG.COL_SN2] >= 0)
        ].copy()
        n_excluded12 = _n_raw9 - len(plot_df12)
        plot_df12["Stabilised"] = plot_df12.get("halide_stabilisation_score",
                                                 pd.Series(0.0, index=plot_df12.index)) > 0.5
        plot_df12["Status"] = plot_df12["Stabilised"].map({True: "Aromatic Shield", False: "Unstabilised"})

        fig, ax = plt.subplots(figsize=(12, 8))
        markers12 = {"Aromatic Shield": "o", "Unstabilised": "X"}
        plot_df12_sorted = plot_df12.copy()
        plot_df12_sorted["_tier_ord"] = plot_df12_sorted[CFG.COL_TIER].map(
            {t: i for i, t in enumerate(reversed(existing_tiers))})
        plot_df12_sorted = plot_df12_sorted.sort_values("_tier_ord")
        sns.scatterplot(data=plot_df12_sorted, x="Dist_Nucleophile", y=CFG.COL_SN2,
                        hue=CFG.COL_TIER, hue_order=existing_tiers, style="Status",
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
                _td12 = plot_df12[plot_df12[CFG.COL_TIER] == tier]
                if len(_td12) < 30:
                    continue
                _col12 = TIER_PALETTE.get(tier, CFG.VIS_INK["faint"])
                try:
                    _kde12 = _gkde12a(np.vstack([
                        _td12["Dist_Nucleophile"].values,
                        _td12[CFG.COL_SN2].values
                    ]))
                    _ZZ12 = _kde12(_grid12).reshape(_XX12.shape)
                    _zmax12 = _ZZ12.max()
                    if _zmax12 <= 0:              # all-zero KDE → skip (dividing gives NaN into contour)
                        continue
                    _ZZ12 = _ZZ12 / _zmax12
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
        ax.axvspan(0, _nuc_cfg12[CFG.TIER_TOP], alpha=0.06, color=CFG.VIS_ACCENT["green"], label="_nolegend_")
        ax.axhspan(_ang_cfg12[CFG.TIER_TOP], 180, alpha=0.06, color=CFG.VIS_ACCENT["green"], label="_nolegend_")

        """
        The distance cut-offs are written UP their own vertical lines, not across the top. Tiers that
        share a cut-off (Tier_1B and Tier_2A are both ≤ 3.2 Å) draw one line and carry one merged
        label, so no label is ever displaced from the line it names.
        """
        _dist_groups12 = {}
        for tier_key in [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]:
            _dist_groups12.setdefault(_nuc_cfg12[tier_key], []).append(tier_key)

        for dist, _tiers12 in _dist_groups12.items():
            col = TIER_PALETTE.get(_tiers12[0], CFG.VIS_INK["faint"])
            ax.axvline(x=dist, color=col, linestyle="--", alpha=0.65, linewidth=1.2)
            _lbl12 = "/".join(t.replace("Tier_", "") for t in _tiers12)
            ax.text(dist, 182, f"{_lbl12} ≤{dist}Å", color=col,
                    fontsize=7.0, rotation=90, rotation_mode="anchor",
                    ha="left", va="center", alpha=0.95, zorder=7,
                    bbox=dict(boxstyle="round,pad=0.18", fc="white", alpha=0.85, ec=col,
                              linewidth=0.5))

        for tier_key in [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]:
            angle = _ang_cfg12[tier_key]
            col = TIER_PALETTE.get(tier_key, CFG.VIS_INK["faint"])
            ax.axhline(y=angle, color=col, linestyle=":", alpha=0.65, linewidth=1.2, zorder=1)
            ax.text(0.012, angle + 0.8, f"≥{angle:g}°  {tier_key}", color=col,
                    fontsize=7.0, ha="left", va="bottom", alpha=0.95, zorder=7,
                    transform=ax.get_yaxis_transform(),
                    bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none",
                              alpha=0.70, linewidth=0))

        # Start the x-axis at the data, not 0 — only ~0.3% of complexes sit below 2 Å, so a 0-start
        # opens a wide empty band. Use the 0.5th percentile (robust to the sparse sub-2 Å tail),
        # floored to 0.1 Å; that lands the axis at ~2 Å where the population actually begins.
        _xmin12 = float(pd.to_numeric(plot_df12["Dist_Nucleophile"], errors="coerce").quantile(0.005))
        ax.set_xlim(max(0.0, np.floor(_xmin12 * 10) / 10), 5.0)
        # Headroom above 180° carries the rotated cut-off labels. It is kept to the minimum they
        # need: any more and the panel opens a band of empty white above the data.
        ax.set_ylim(0, 200)
        ax.set_yticks(range(0, 181, 20))
        # The distance is measured to EACH variant's own mapped catalytic aspartate (Mapped_Nucleophile:
        # Asp109 / Asp110 / Asp112 / Asp104 …), never a fixed residue number — a FAcD variant's nucleophile
        # can sit a little before or after the reference Asp110. So the axis names the ROLE, not one number.
        ax.set_xlabel('Mapped catalytic Asp (nucleophile) → Carbon Distance (Å)  — shorter = closer to reaction geometry',
                      )
        ax.set_ylabel("SN2 Attack Angle (°)  — 180° = perfect linear back-attack", )
        handles12, labels12 = ax.get_legend_handles_labels()
        clean_labels12, clean_handles12 = [], []
        for h12, l12 in zip(handles12, labels12):
            if l12 in (CFG.COL_TIER, "Status"):
                continue
            clean_handles12.append(h12)
            clean_labels12.append(l12)
        # Legend inside the lower-left (the sparse substrate zone), with the population count as its
        # TITLE so the "n = … shown / excluded" line and the tier key read as one merged block.
        _leg12 = ax.legend(clean_handles12, clean_labels12,
                           loc="lower left", bbox_to_anchor=(0.008, 0.008),
                           fancybox=True, ncol=max(2, (len(clean_handles12) + 1) // 2),
                           title=f"n = {len(plot_df12):,} complexes shown  ·  {n_excluded12:,} excluded (dist = 0 or ≥ 5 Å)")
        _leg12.get_title().set_fontsize(_JF_FA)
        _leg12.get_title().set_color(CFG.VIS_INK["mid"])
        ax.set_axisbelow(True)
        ax.grid(color=CFG.VIS_INK["wash"], linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_13a_Mechanism_Geometry_Scatter.png",
                    dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    # Figure 13b: PA companion (mechanistic scatter overlay)
    if _tt_has_imgs:
        _fig_13b_tt_mechanistic(df, _pa, _imgs, out_dir, reporter)


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
                    CFG.COL_TIER, "scissile_cf_bde", "sn2_backside_occlusion",
                    "sn2_dead_end", CFG.COL_LIG]
        if not all(_c in df.columns for _c in _need_tc):
            reporter.log("  ! Two-criteria figure skipped: required columns absent")
        else:
            _BF = CFG.TIER_CONSTELLATION_MIN["Tier_1A"]
            _BDE_MAX = CFG.SCISSILE_CF_BDE_MAX
            _OCC_MAX = CFG.SN2_BACKSIDE_OCCL_MAX
            _tord = [t for t in TIER_ORDER_LOGIC if t in df[CFG.COL_TIER].unique()]
            """
            Three INDEPENDENT figures, not one three-panel strip. Each answers its own question, and
            each carries its own axis labels and legend at a readable size — sharing a strip lets the
            widest panel dictate the size of the others.
            """
            _figa = plt.figure(figsize=(13, 7))

            # ── Figure 10 : Criterion A (integrity) gates Criterion B (constellation) ──
            _axa = _figa.add_subplot(111)
            _A = pd.to_numeric(df["active_site_residues_correct"], errors="coerce") / 8.0
            _B = pd.to_numeric(df["catalytic_constellation_score"], errors="coerce")
            _sa = pd.DataFrame({"A": _A, "B": _B}).dropna()
            _abins = sorted(_sa["A"].round(3).unique())
            _adata = [_sa.loc[_sa["A"].round(3).eq(b), "B"].values for b in _abins]
            _ans = [len(v) for v in _adata]
            _parts = _axa.violinplot(_adata, positions=range(len(_abins)), widths=0.85, showextrema=False)
            for _pc, _col in zip(_parts["bodies"], plt.cm.RdYlGn(np.linspace(0.1, 0.9, len(_abins)))):
                _pc.set(facecolor=_col, alpha=0.7, edgecolor=CFG.VIS_INK["muted"], linewidth=0.5)
            _amed = [float(np.median(v)) if len(v) else np.nan for v in _adata]
            _axa.plot(range(len(_abins)), _amed, "-D", color=CFG.VIS_ACCENT["bad"], lw=2, ms=6, zorder=6, label="Median B per bin")
            _axa.axhline(_BF, ls="--", color=CFG.VIS_ACCENT["bad"], lw=1.6, label=f"Criterion-B floor for Tier_1A ({_BF:g})")
            for _i, _b in enumerate(_abins):
                _pct = 100.0 * (_sa.loc[_sa["A"].round(3).eq(_b), "B"] >= _BF).mean()
                _axa.text(_i, 1.03, f"{_pct:.0f}%", ha="center", va="bottom", fontsize=8.5, fontweight="bold", color=CFG.VIS_ACCENT["axis_left"])
            _axa.text(0.015, 0.965, f"Top % = share with B ≥ {_BF:g}", transform=_axa.transAxes,
                      fontsize=8, color=CFG.VIS_ACCENT["axis_left"], va="top", ha="left", fontweight="bold")
            _axa.set_xticks(range(len(_abins)))
            _axa.set_xticklabels([f"{b:.3g}\n(n={n:,})" for b, n in zip(_abins, _ans)], fontsize=8)
            for _tla, _cla in zip(_axa.get_xticklabels(), plt.cm.RdYlGn(np.linspace(0.1, 0.9, len(_abins)))):
                _tla.set_color(_cla)   # each bin's tick label matches its violin colour
            _axa.set_ylim(0, 1.18)
            _axa.set_xlabel("Active-site integrity  (fraction of the 8 catalytic residues placed)", )
            _axa.set_ylabel(f"Catalytic constellation score  (0–1 vs {CFG.REFERENCE_PDB_ID} crystal)", )
            _axa.legend(loc="upper right",  ncol=2)
            _axa.grid(True, color=CFG.VIS_INK["palest"], linewidth=0.6, alpha=0.75, zorder=0); _axa.set_axisbelow(True)

            # ── panel b : Criterion-B ECDF by tier + arrowed per-tier mean values ──
            _figb = plt.figure(figsize=(12, 7))
            _axb = _figb.add_subplot(111)
            # Shaded per-tier constellation zones (mirrors the SN2-angle ECDF, Fig 12): each zone's
            # lower bound is that tier's TIER_CONSTELLATION_MIN gate, its upper bound the next-higher
            # tier's gate (top capped at 1.0). Each band carries its ≥threshold label and the count of
            # complexes whose score falls inside it — so the reader sees the gate AND the population.
            _zt_cb  = [CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2], CFG.TIER_ORDER[3]]
            _zlo_cb = [CFG.TIER_CONSTELLATION_MIN[t] for t in _zt_cb]
            _zhi_cb = [1.0] + _zlo_cb[:-1]
            _score_all_cb = pd.to_numeric(df["catalytic_constellation_score"], errors="coerce")
            for _lo, _hi, _tk in zip(_zlo_cb, _zhi_cb, _zt_cb):
                _zc = TIER_PALETTE.get(_tk, CFG.VIS_INK["faint"])
                _axb.axvspan(_lo, _hi, alpha=0.07, color=_zc, zorder=0)
                _axb.text((_lo + _hi) / 2, 1.02, f"≥{_lo:.2f}", ha="center", va="bottom",
                          fontsize=_JF_FA, color=_zc, fontweight="bold", transform=_axb.get_xaxis_transform())
                _ncb = int(((_score_all_cb >= _lo) & (_score_all_cb < _hi)).sum())
                _axb.text((_lo + _hi) / 2, 0.10, f"n={_ncb:,}", ha="center", va="center",
                          fontsize=7.0, color=_zc, fontweight="bold", rotation=90, zorder=10,
                          transform=_axb.get_xaxis_transform(),
                          bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.68))
            # ECDF per tier + DKW 95% band + median dot (median labels placed after the loop in
            # alternating vertical bands so clustered medians never overlap).
            _med_cb = []
            for _t in _tord:
                _v = pd.to_numeric(df.loc[df[CFG.COL_TIER].eq(_t), "catalytic_constellation_score"],
                                   errors="coerce").dropna().sort_values().values
                if len(_v) < 5:
                    continue
                _yv = np.arange(1, len(_v) + 1) / len(_v)
                _cc = TIER_PALETTE.get(_t, CFG.VIS_INK["faint"])
                _axb.step(_v, _yv, where="post", color=_cc, lw=2.2, label=f"{_t}  (n={len(_v):,})", zorder=3)
                _eps_cb = np.sqrt(np.log(2 / 0.05) / (2 * len(_v)))
                _axb.fill_between(_v, np.clip(_yv - _eps_cb, 0, 1), np.clip(_yv + _eps_cb, 0, 1),
                                  color=_cc, alpha=0.15, zorder=1)
                _mdv = float(np.median(_v)); _mdy = float(np.interp(_mdv, _v, _yv))
                _axb.scatter([_mdv], [_mdy], color=_cc, s=70, zorder=6, edgecolors="black", linewidths=0.8)
                _med_cb.append((_mdv, _mdy, _cc))
            _band_cb = [40, -40, 24, -24]
            for _k, (_mx, _my, _mc) in enumerate(sorted(_med_cb, key=lambda z: z[0])):
                _oy = _band_cb[_k % len(_band_cb)]
                _axb.annotate(f"{_mx:.2f}", (_mx, _my), xytext=(0, _oy), textcoords="offset points",
                              fontsize=8, color=_mc, fontweight="bold", ha="center",
                              va="bottom" if _oy > 0 else "top", zorder=7,
                              bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=_mc, alpha=0.85, linewidth=0.6),
                              arrowprops=dict(arrowstyle="->", color=_mc, lw=0.7, shrinkA=0, shrinkB=3))
            # MD-ready overlay: top-performing MD-selected complexes on their tier's B-ECDF curve.
            _mdb = _md_ready_df(df)
            for _, _mr in _mdb.iterrows():
                _bx = pd.to_numeric(pd.Series([_mr.get("catalytic_constellation_score")]), errors="coerce").iloc[0]
                _bt = _mr.get(CFG.COL_TIER)
                if pd.isna(_bx):
                    continue
                _bv = pd.to_numeric(df.loc[df[CFG.COL_TIER].eq(_bt), "catalytic_constellation_score"],
                                    errors="coerce").dropna().sort_values().values
                if len(_bv) == 0:
                    continue
                _axb.scatter([_bx], [float(np.interp(_bx, _bv, np.linspace(0, 1, len(_bv))))], **_MD_STAR_KW)
            for _, _cr_b in _control_star_df(df).iterrows():
                _cbx = pd.to_numeric(pd.Series([_cr_b.get("catalytic_constellation_score")]), errors="coerce").iloc[0]
                _cbt = _cr_b.get(CFG.COL_TIER)
                if pd.isna(_cbx):
                    continue
                _cbv = pd.to_numeric(df.loc[df[CFG.COL_TIER].eq(_cbt), "catalytic_constellation_score"],
                                     errors="coerce").dropna().sort_values().values
                if len(_cbv) == 0:
                    continue
                _ctrl_star(_axb, [_cbx], [float(np.interp(_cbx, _cbv, np.linspace(0, 1, len(_cbv))))])
            if len(_mdb):
                _axb.scatter([], [], marker="*", s=140, facecolor=CFG.VIS_ACCENT["star"], edgecolor="black",
                             linewidths=0.9, label=f"MD-ready (n={len(_mdb)})")
            if not _control_star_df(df).empty:
                _axb.scatter([], [], marker="*", s=160, facecolor=CFG.VIS_ACCENT["control"],
                             edgecolor=CFG.VIS_ACCENT["control_edge"], linewidths=1.4, label="3R3U × FA (control)")
            # Per-tier constellation threshold lines (colours match the tiers/zones); the Tier_1A line
            # is the Criterion-B floor. A dotted line marks the crystal-grade constellation that
            # compensates a near-elite mech for Tier_1A (MECH_ELITE_CONSTELLATION).
            for _t in _zt_cb[1:]:   # Tier_1A's line is the Criterion-B floor, drawn once below
                _axb.axvline(CFG.TIER_CONSTELLATION_MIN[_t],
                             color=TIER_PALETTE.get(_t, CFG.VIS_INK["faint"]), ls="--", alpha=0.65, lw=1.0, zorder=2)
            _axb.axvline(_BF, ls="--", color=CFG.VIS_ACCENT["bad"], lw=1.4, alpha=0.7,
                         label=f"Criterion-B floor for {CFG.TIER_TOP} ({_BF:g})")
            _axb.axvline(CFG.MECH_ELITE_CONSTELLATION, ls=":", color=CFG.VIS_INK["dark"], lw=1.3, alpha=0.8,
                         zorder=2, label=f"crystal-grade ({CFG.MECH_ELITE_CONSTELLATION:g})")
            _axb.set_xlabel("Catalytic constellation score  (reactive-geometry match, 0–1)", )
            _axb.set_ylabel("Cumulative fraction", )
            """
            The legend is stacked in one narrow column at the top left — the one corner no ECDF
            crosses, since every tier is still flat at zero below a constellation score of ~0.17.
            Strung out along a single row (one column per tier) the entries would shrink below
            legibility and still collide. The criterion floor carries its own label rather than
            being an unexplained red line.
            """
            _axb.legend(loc="upper left", ncol=1)
            # y every 0.1; x-ticks carry the base 0.2 grid PLUS every tier gate + the crystal-grade
            # mark, so each coloured threshold line is read off a labelled tick (0.05 minor between).
            from matplotlib.ticker import MultipleLocator as _MLoc_cb
            _axb.set_xlim(0.0, 1.0)
            _axb.set_ylim(0.0, 1.06)
            _axb.set_xticks(sorted(set([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
                                       + [round(v, 2) for v in _zlo_cb]
                                       + [round(CFG.MECH_ELITE_CONSTELLATION, 2)])))
            _axb.yaxis.set_major_locator(_MLoc_cb(0.1))
            _axb.xaxis.set_minor_locator(_MLoc_cb(0.05))
            _axb.tick_params(axis="both", which="major", labelsize=8)
            _axb.grid(True, color=CFG.VIS_INK["palest"], linewidth=0.6, alpha=0.75, zorder=0); _axb.set_axisbelow(True)

            """
            ── Figure 12 : the SN2 dead-end chemistry gate ──────────────────────────────────────
            A scatter of BDE against occlusion cannot work here, and no amount of label-nudging
            fixes it: the ligands SHARE coordinates. Five perfluoroalkyl sulfonates have the same
            C–F bond-dissociation energy and the same backside occlusion, so they land on one dot,
            and their names have to be stacked or fanned around it. The figure was spending all its
            ink on the collision and none on the chemistry.

            One row per ligand removes the collision by construction. Each ligand is a row; the two
            criteria are read side by side against their own gate; and the verdict — dead end, or
            feasible α-attack — is the row's colour. Nothing can overlap, and the question the
            figure exists to answer (WHICH criterion killed this ligand) is now answerable at a
            glance: a bar crossing into the shaded zone is the one that failed.
            """
            _DEAD, _FEAS = CFG.VIS_ACCENT_DEEP["red_bright"], CFG.VIS_ACCENT["blue"]
            _dd = df.dropna(subset=["scissile_cf_bde", "sn2_backside_occlusion"])
            _gg = _dd.groupby(CFG.COL_LIG).agg(
                bde=("scissile_cf_bde", "median"), occ=("sn2_backside_occlusion", "median"),
                dead=("sn2_dead_end", lambda s: pd.to_numeric(s, errors="coerce").fillna(0).max())).reset_index()
            _gg["pen"] = (_gg["dead"] > 0) | (_gg["bde"] > _BDE_MAX) | (_gg["occ"] > _OCC_MAX)
            _keep = _gg[_gg["pen"] | _gg[CFG.COL_LIG].isin(["26_Fluoroacetate", "27_Difluoroacetate"])].copy()

            _figc = plt.figure(figsize=(12, 7))
            _gsc = _figc.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.06)
            _axc = _figc.add_subplot(_gsc[0, 0])          # criterion 1 — bond strength
            _axc2 = _figc.add_subplot(_gsc[0, 1], sharey=_axc)   # criterion 2 — backside access

            if len(_keep) >= 2:
                _keep["lig"] = _lig_short_series(_keep[CFG.COL_LIG])
                _keep["isdead"] = _keep["dead"] > 0
                # order by the verdict, then by bond strength: the feasible controls sit together at
                # the foot of the chart, the dead ends above them.
                _keep = _keep.sort_values(["isdead", "bde", "occ"], ascending=[True, True, True])
                _y = np.arange(len(_keep))
                _cols = [_DEAD if _d else _FEAS for _d in _keep["isdead"]]

                # left — scissile C–F BDE
                _axc.axvspan(_BDE_MAX, _keep.bde.max() + 2.5, color=_DEAD, alpha=0.07, zorder=0)
                _axc.hlines(_y, _keep.bde.min() - 2.0, _keep["bde"], color=_cols, lw=1.4, alpha=0.55, zorder=2)
                _axc.scatter(_keep["bde"], _y, s=110, color=_cols, edgecolor="white", lw=0.8, zorder=3)
                _axc.axvline(_BDE_MAX, ls="--", color=_DEAD, lw=1.4)
                for _yy, _v in zip(_y, _keep["bde"]):
                    _axc.annotate(f"{_v:.1f}", (_v, _yy), xytext=(7, 0), textcoords="offset points",
                                  va="center", fontsize=8, color=CFG.VIS_INK["dark"])
                _axc.set_xlim(_keep.bde.min() - 2.0, _keep.bde.max() + 2.5)
                _axc.set_xlabel(f"Scissile C–F BDE  (kcal/mol)   ·   gate > {_BDE_MAX:g}", )
                _axc.set_yticks(_y)
                _axc.set_yticklabels(_keep["lig"], fontsize=9.5, fontweight="bold")
                for _tl, _cc in zip(_axc.get_yticklabels(), _cols):
                    _tl.set_color(_cc)

                # right — backside steric occlusion
                _axc2.axvspan(_OCC_MAX, max(_keep.occ.max() + 0.35, _OCC_MAX + 0.35),
                              color=_DEAD, alpha=0.07, zorder=0)
                _axc2.hlines(_y, 0, _keep["occ"], color=_cols, lw=1.4, alpha=0.55, zorder=2)
                _axc2.scatter(_keep["occ"], _y, s=110, color=_cols, edgecolor="white", lw=0.8, zorder=3)
                _axc2.axvline(_OCC_MAX, ls="--", color=_DEAD, lw=1.4)
                for _yy, _v in zip(_y, _keep["occ"]):
                    _axc2.annotate(f"{_v:.2f}", (_v, _yy), xytext=(7, 0), textcoords="offset points",
                                   va="center", fontsize=8, color=CFG.VIS_INK["dark"])
                _axc2.set_xlim(-0.05, max(_keep.occ.max() + 0.5, _OCC_MAX + 0.5))
                _axc2.set_xlabel(f"Backside steric occlusion  (Σ vdW, Å)   ·   gate > {_OCC_MAX:g}", )
                _axc2.tick_params(labelleft=False)
                _axc.set_ylim(-0.8, len(_keep) - 0.2)

            _axc.legend(handles=[_L2Dtc([], [], marker="o", ls="", color=_DEAD, label="SN2 dead-end (flagged)"),
                                 _L2Dtc([], [], marker="o", ls="", color=_FEAS, label="feasible α-attack (control)"),
                                 _L2Dtc([], [], ls="--", color=_DEAD, label="gate (shaded = fails it)")],
                        loc="lower right")
            for _a in (_axc, _axc2):
                _a.grid(True, axis="x", color=CFG.VIS_INK["palest"], linewidth=0.6, alpha=0.75, zorder=0)
                _a.set_axisbelow(True)

            _tier_seps(_axa, len(_abins))   # consistent separators (replaces ad-hoc grid)
            _figa.savefig(out_dir / "Figure_30a_Criterion_A_Gates_B.png",
                          dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            _figb.savefig(out_dir / "Figure_30b_Criterion_B_ECDF_by_Tier.png",
                          dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            _figc.savefig(out_dir / "Figure_30c_SN2_DeadEnd_Gate.png",
                          dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close(_figa); plt.close(_figb); plt.close(_figc)
    except Exception as e:
        reporter.log(f"  ! Two-criteria figure skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig



# =============================================================================
# STEP 6/8 — 06_Ligand_Interactions_and_Chemical_Space
# =============================================================================
def _fig_folder06_ligand(df, features, out_dir, reporter, existing_tiers, _pa, _imgs, _tt_has_imgs, int_label_map, present_int_cols):
    """Folder 06_Ligand_Interactions_and_Chemical_Space — interaction profile + chemical space."""
    reporter.section("Step 6/8 — 06_Ligand_Interactions_and_Chemical_Space · interaction profile + chemical space")
    # --- Figure 14: Molecular Interaction Profile — Stacked bars (bond types) + F-engagement line ---
    """
    Single chart, dual y-axis: stacked bars per tier show the total interaction count and
    its bond-type composition (left y-axis); overlaid connected dot-line shows the fluorine
    engagement ratio per tier on the right y-axis, revealing whether richer interaction
    profiles correlate with higher fluorine utilisation.
    """
    _has_int   = bool(present_int_cols) and CFG.COL_TIER in df.columns
    _has_f16   = ("interacting_fluorine_count" in df.columns and
                  "total_fluorine_count" in df.columns and CFG.COL_TIER in df.columns)
    if _has_int or _has_f16:
        try:
            # Pre-compute tier count for dynamic figure height
            n_tiers_f13 = len([t for t in existing_tiers if t in df[CFG.COL_TIER].values]) if CFG.COL_TIER in df.columns else 6
            fig, (ax, ax_hm) = plt.subplots(1, 2, figsize=(20, max(5.5, n_tiers_f13 * 0.95 + 2.5)),
                                             gridspec_kw={"width_ratios": [3, 2], "wspace": 0.35})
            ax_r15 = ax.twinx()   # right y-axis for fluorine engagement line

            # Left axis: stacked bar chart — one bar per tier, stacked by bond type
            if _has_int:
                tier_int = df.groupby(CFG.COL_TIER)[[present_int_cols[k] for k in present_int_cols]].mean()
                tier_int = tier_int.reindex([t for t in existing_tiers if t in tier_int.index])
                tier_int.columns = [int_label_map.get(k, k) for k in present_int_cols]
                tier_int = tier_int.fillna(0)
                # Normalise to 100 % so each bar shows interaction-type composition
                _tier_int_abs     = tier_int.copy()
                _row_totals_abs   = _tier_int_abs.sum(axis=1)
                tier_int          = _tier_int_abs.div(_row_totals_abs, axis=0).fillna(0).mul(100)
                # Stacked bars: each bar = one tier; stack = bond type proportion
                bond_palette = dict(CFG.BOND_TYPE_COLOUR)
                bar_bottom = np.zeros(len(tier_int))
                bar_x = np.arange(len(tier_int))
                tier_labels = list(tier_int.index)
                for bond_type in list(int_label_map.values()):
                    if bond_type not in tier_int.columns:
                        continue
                    vals_pct = tier_int[bond_type].values           # percent for bar height
                    vals_abs = _tier_int_abs[bond_type].values      # absolute mean for annotation
                    ax.bar(bar_x, vals_pct, bottom=bar_bottom,
                           color=bond_palette.get(bond_type, CFG.VIS_INK["faint"]),
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
                ax.tick_params(axis="y", labelsize=9, labelcolor=CFG.VIS_ACCENT["axis_left"])
                ax.yaxis.label.set_color(CFG.VIS_INK["near_black"])
                ax.set_xlabel("Degrader Tier", )
                ax.set_ylabel("Interaction type proportion  (%)  —  annotation = mean count", )

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
                    if tier not in f16_df[CFG.COL_TIER].values:
                        continue
                    vals = f16_df.loc[f16_df[CFG.COL_TIER] == tier, "FER"].dropna()
                    if len(vals) < 2:
                        continue
                    med_v = float(vals.median())
                    lo_ci, hi_ci = _median_ci95(vals)
                    fer_xs.append(i)
                    fer_meds.append(med_v)
                    fer_lo_ci.append(max(0, lo_ci))
                    fer_hi_ci.append(min(1, hi_ci))
                if fer_xs:
                    ax_r15.plot(fer_xs, fer_meds, color=CFG.VIS_ACCENT["amber"], linewidth=2.6,
                                marker="o", markersize=9, markeredgecolor="black",
                                markeredgewidth=0.9, zorder=8,
                                label="Fluorine Engagement Ratio  (median ± 95% CI)")
                    ax_r15.fill_between(fer_xs, fer_lo_ci, fer_hi_ci,
                                        color=CFG.VIS_ACCENT["amber"], alpha=0.22, zorder=7)
                    # Filled zone backgrounds on right axis
                    ax_r15.axhspan(0.75, 1.20, alpha=0.15, color=CFG.VIS_ACCENT["green"], zorder=0)
                    ax_r15.axhspan(0.50, 0.75, alpha=0.13, color=CFG.VIS_ACCENT["blue"], zorder=0)
                    ax_r15.axhspan(0.00, 0.50, alpha=0.12, color=CFG.VIS_ACCENT["vermillion"], zorder=0)
                    ax_r15.axhline(y=0.5, color=CFG.VIS_ACCENT["blue"], linestyle="--",
                                   alpha=0.70, linewidth=1.3, zorder=6)
                    ax_r15.axhline(y=1.0, color=CFG.VIS_ACCENT["green"], linestyle=":",
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
                                        ha="center", va="bottom", fontsize=_JF_FA,
                                        color=CFG.VIS_ACCENT["axis_right"], fontweight="bold", zorder=10,
                                        arrowprops=dict(arrowstyle="->", color=CFG.VIS_ACCENT["amber"],
                                                        lw=0.5, mutation_scale=6,
                                                        shrinkA=0, shrinkB=2),
                                        bbox=dict(boxstyle="round,pad=0.09", fc="white",
                                                  ec=CFG.CONF_BAND_COLOURS["acceptable"], alpha=0.94, linewidth=0.5))
                ax_r15.set_ylim(0, 1.30)
                ax_r15.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
                ax_r15.set_yticklabels(["0%", "25%", "50%", "75%", "100%"], fontsize=9)
                ax_r15.set_ylabel("Fluorine Engagement Ratio  (interacting F / total F per ligand)",
                                   rotation=270, labelpad=14, color=CFG.VIS_ACCENT["axis_right"])
                ax_r15.tick_params(axis="y", labelcolor=CFG.VIS_ACCENT["axis_right"], labelsize=9)
                """
                Rotated zone labels — left side of figure, written bottom-to-top
                ylim=(0, 1.30): Low midpoint=0.25→0.192; Mid midpoint=0.625→0.481; High midpoint=0.875→0.673
                """
                ax_r15.text(0.99, 0.192, "Low engagement (<50%)",
                            ha="center", va="center", rotation=90,
                            fontsize=_JF_FA, color=CFG.VIS_ACCENT["vermillion"], style="italic", fontweight="bold",
                            transform=ax_r15.transAxes,
                            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["below"],
                                      alpha=0.80, linewidth=0.5))
                ax_r15.text(0.99, 0.481, "50% engagement",
                            ha="center", va="center", rotation=90,
                            fontsize=_JF_FA, color=CFG.VIS_ACCENT["blue"], style="italic", fontweight="bold",
                            transform=ax_r15.transAxes,
                            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.VIS_ACCENT["blue"],
                                      alpha=0.80, linewidth=0.5))
                ax_r15.text(0.99, 0.673, "Full engagement",
                            ha="center", va="center", rotation=90,
                            fontsize=_JF_FA, color=CFG.VIS_ACCENT["green"], style="italic", fontweight="bold",
                            transform=ax_r15.transAxes,
                            bbox=dict(boxstyle="round,pad=0.10", fc="white", ec=CFG.CONF_BAND_COLOURS["high"],
                                      alpha=0.80, linewidth=0.5))

            """
            Dual-colour grid system: left = steel-blue (matches left tick labels),
            right = amber (matches right tick labels). Both BEHIND bars via set_axisbelow.
            """
            ax.set_axisbelow(True)
            ax.yaxis.grid(True, color=CFG.VIS_ACCENT["axis_left"], linewidth=0.55, linestyle="-",
                          alpha=0.20, zorder=0)
            ax.xaxis.grid(False)
            ax_r15.yaxis.grid(True, color=CFG.VIS_ACCENT["axis_right"], linewidth=0.55, linestyle="--",
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
                                 fancybox=True)
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
                ax_hm.set_xlabel("Degrader Tier  —  Interaction Type Proportion (%)", )
                ax_hm.set_ylabel("Interaction Type", )
                ax_hm.tick_params(axis="x", labelrotation=45, labelsize=8.5)
                ax_hm.tick_params(axis="y", labelrotation=0,  labelsize=8)
            else:
                ax_hm.set_visible(False)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_14a_Molecular_Interaction_Profile.png",
                        dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('14')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # Figure 14b: PA companion (interaction profile overlay)
    if _tt_has_imgs:
        _fig_14b_tt_interactions(df, _pa, _imgs, out_dir, reporter)

    # --- Figure 15: Fluorine Engagement Ratio by Tier (box + strip + median trend line) ---
    if "interacting_fluorine_count" in df.columns and "total_fluorine_count" in df.columns and CFG.COL_TIER in df.columns:
        try:
            f16_df = df.copy()
            f16_df["FER"] = np.where(
                f16_df["total_fluorine_count"] > 0,
                f16_df["interacting_fluorine_count"] / f16_df["total_fluorine_count"],
                np.nan)
            f16_df = f16_df.dropna(subset=["FER"])
            valid_t16 = [t for t in existing_tiers if t in f16_df[CFG.COL_TIER].values]
            fig, ax = plt.subplots(figsize=(12, 7))

            ax.axhspan(0.75, 1.01, alpha=0.07, color=CFG.VIS_ACCENT["green"], zorder=0)
            ax.axhspan(0.50, 0.75, alpha=0.06, color=CFG.VIS_ACCENT["amber"], zorder=0)
            ax.axhspan(0.00, 0.50, alpha=0.06, color=CFG.VIS_ACCENT["vermillion"], zorder=0)
            ax.axhline(y=0.75, color=CFG.VIS_ACCENT["green"], linestyle="--", alpha=0.6, linewidth=1.1)
            ax.axhline(y=0.50, color=CFG.VIS_ACCENT["amber"], linestyle=":", alpha=0.6, linewidth=1.1)
            ax.set_axisbelow(True)
            # KDE violin + IQR box overlay (no boxplot — keeps figure clean)
            _skipped14 = []   # tiers with too few points for a KDE violin (kept visible via a note)
            for _ti14, tier14 in enumerate(valid_t16):
                _vals14 = f16_df.loc[f16_df[CFG.COL_TIER] == tier14, "FER"].dropna().values
                if len(_vals14) < 4:
                    # Too few complexes for a stable KDE — draw the raw points so the
                    # tier is not silently omitted, and record it for the caption.
                    if len(_vals14):
                        ax.scatter([_ti14] * len(_vals14), _vals14, s=26, zorder=6,
                                   color=TIER_PALETTE.get(tier14, CFG.VIS_INK["faint"]),
                                   edgecolors="black", linewidths=0.5)
                    _skipped14.append((tier14, len(_vals14)))
                    continue
                from scipy.stats import gaussian_kde as _gkde14
                try:
                    _kde14  = _gkde14(_vals14)
                    _y14    = np.linspace(np.min(_vals14), np.max(_vals14), 200)
                    _w14    = _kde14(_y14)
                    _wmax14 = _w14.max()
                    if _wmax14 <= 0:             # all-zero KDE → skip (division would NaN the violin)
                        continue
                    _w14   /= _wmax14
                    _w14   *= 0.35   # half-width of violin
                    ax.fill_betweenx(_y14, _ti14 - _w14, _ti14 + _w14,
                                     alpha=0.55, color=TIER_PALETTE.get(tier14, CFG.VIS_INK["faint"]))
                except (np.linalg.LinAlgError, ValueError):
                    # Degenerate (collinear) FER values give a singular KDE covariance;
                    # render the raw points as a jittered strip in place of the violin.
                    _jit14 = np.random.default_rng(0).uniform(-0.12, 0.12, size=len(_vals14))
                    ax.scatter(_ti14 + _jit14, _vals14, s=12, alpha=0.45,
                               color=TIER_PALETTE.get(tier14, CFG.VIS_INK["faint"]),
                               edgecolors="none", zorder=2)
                # IQR box on top
                _q25, _q75 = np.percentile(_vals14, [25, 75])
                _med14 = np.median(_vals14)
                ax.vlines(_ti14, _q25, _q75, color=CFG.VIS_INK["dark"], linewidth=2.5, zorder=5)
                ax.scatter([_ti14], [_med14], color="white", s=30, zorder=6,
                           edgecolors=CFG.VIS_INK["dark"], linewidths=1.2)

            # Per-tier median statistics badge
            _med16_xs, _med16_ys = [], []
            for i, tier in enumerate(valid_t16):
                sub16 = f16_df.loc[f16_df[CFG.COL_TIER] == tier, "FER"].dropna()
                if len(sub16) == 0:
                    continue
                med = float(sub16.median())
                n16 = len(sub16)
                _med16_xs.append(i)
                _med16_ys.append(med)
            # FER connecting line across tier medians
            if len(_med16_xs) >= 2:
                ax.plot(_med16_xs, _med16_ys, color=CFG.VIS_ACCENT["amber"], linewidth=2.2,
                        linestyle="--", alpha=0.85, zorder=4,
                        marker="D", markersize=6, markeredgecolor=CFG.VIS_INK["dark"],
                        markeredgewidth=0.7, label="FER trend (median)")
            for i, tier in enumerate(valid_t16):
                sub16 = f16_df.loc[f16_df[CFG.COL_TIER] == tier, "FER"].dropna()
                if len(sub16) == 0:
                    continue
                med = float(sub16.median())
                n16 = len(sub16)
                ax.text(i, med + 0.022, f"med={med:.2f}\nn={n16:,}",
                        ha="center", va="bottom", fontsize=7.0, fontweight="bold",
                        zorder=7, color=CFG.VIS_INK["near_black"],
                        bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=CFG.VIS_INK["palest"],
                                  linewidth=0.6, alpha=0.85))

            ax.set_ylim(0, 1.08)
            ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
            ax.set_yticklabels(["0%", "25%", "50%", "75%", "100%"])
            # Color y-tick labels to match zone backgrounds
            _ytick_colors16 = {
                "0%":   CFG.VIS_RAMP["orange"][3],
                "25%":  CFG.VIS_RAMP["orange"][3],
                "50%":  CFG.VIS_BAND["moderate"],
                "75%":  CFG.VIS_BAND["high"],
                "100%": CFG.VIS_BAND["high"],
            }
            for tick16 in ax.get_yticklabels():
                tick16.set_color(_ytick_colors16.get(tick16.get_text(), CFG.VIS_INK["near_black"]))
            """
            Start-y (axes fraction) = bottom border of each label's OWN zone, so
            every label sits on its category's lower boundary line (va='bottom',
            text rises into its zone). Low's zone bottom coincides with the axes
            floor, which is expected — High/Moderate now match that anchoring.
            """
            for (_zy14_start, _ztxt14, _zcol14, _zec14) in [
                (0.73, f"High engagement (≥{CFG.ENGAGEMENT_BAND_HIGH*100:.0f}%)", CFG.VIS_BAND["high"], CFG.VIS_ACCENT["green"]),
                (0.50, f"Moderate ({CFG.ENGAGEMENT_BAND_MODERATE*100:.0f}–{CFG.ENGAGEMENT_BAND_HIGH*100:.0f}%)", CFG.VIS_BAND["moderate"], CFG.VIS_ACCENT["amber"]),
                (0.10, f"Low engagement (<{CFG.ENGAGEMENT_BAND_MODERATE*100:.0f}%)", CFG.VIS_RAMP["orange"][3], CFG.VIS_ACCENT["vermillion"]),
            ]:
                ax.text(0.985, _zy14_start, _ztxt14,
                        ha="center", va="bottom", clip_on=True,
                        fontsize=7.0, color=_zcol14, style="italic",
                        fontweight="bold", rotation=90, transform=ax.transAxes,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=_zec14,
                                  alpha=0.80, linewidth=0.5))
            ax.set_xlabel("Degrader Tier", labelpad=-18)   # lift into the whitespace above the angled tier ticks
            ax.set_ylabel("Fluorine Engagement Ratio  (interacting F / total F)", )
            ax.set_xticks(range(len(valid_t16)))
            ax.set_xticklabels(valid_t16, rotation=35, ha="right", fontsize=9)
            for tick16x, tier16x in zip(ax.get_xticklabels(), valid_t16):
                tick16x.set_color(TIER_PALETTE.get(tier16x, "black"))

            # Legend — violin + IQR box description
            from matplotlib.patches import Patch as _P16
            from matplotlib.lines import Line2D as _L16
            _leg16_h = [
                _P16(facecolor=CFG.VIS_INK["ghost"], alpha=0.55, label="KDE violin (per tier)"),
                _L16([0], [0], color=CFG.VIS_INK["dark"], linewidth=2.5, label="IQR (25–75%)"),
                _L16([0], [0], color="white", marker="o", markersize=6,
                     markeredgecolor=CFG.VIS_INK["dark"], markeredgewidth=1.2,
                     linestyle="None", label="Median"),
            ]
            _mdh16 = _md_ready_stars_cat(ax, f16_df, valid_t16, "FER")
            _leg16_h.extend(_mdh16)
            _leg16 = ax.legend(handles=_leg16_h,
                               loc="lower left", bbox_to_anchor=(0.01, 0.01),
                                 fancybox=True, ncol=4)
            _leg16.set_zorder(20)

            if _skipped14:
                _sk_txt = ", ".join(f"{t} (n={n})" for t, n in _skipped14)
                ax.text(0.99, 0.015, f"n<4, shown as points (no KDE): {_sk_txt}",
                        transform=ax.transAxes, ha="right", va="bottom",
                        fontsize=7.0, style="italic", color=CFG.VIS_INK["muted"], zorder=20)

            plt.tight_layout()
            _tier_seps(plt.gca())   # consistent vertical tier separators
            plt.savefig(out_dir / "Figure_15_Fluorine_Engagement_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('15')} — {e}")
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
    _qual_col  = next((c for c in ["soft_catalytic_score", CFG.COL_SN2, CFG.COL_MECH_S] if c in df.columns), None)
    _qual_lbl  = {"soft_catalytic_score": "Soft catalytic score  (SN2-geometry composite)",
                  CFG.COL_SN2: "SN2 attack angle (°)",
                  CFG.COL_MECH_S: "Mechanistic score"}.get(_qual_col, str(_qual_col))
    _qual_fmt  = "{:.1f}" if _qual_col == CFG.COL_SN2 else "{:.3f}"
    _has_bind  = _qual_col is not None and CFG.COL_TIER in df.columns
    _has_dens  = "Interaction_Density_Norm" in df.columns and CFG.COL_TIER in df.columns
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
                sns.violinplot(data=f12_data, x=CFG.COL_TIER, y=_qual_col, order=existing_tiers,
                               palette=TIER_PALETTE, inner="box", linewidth=1.2, cut=0,
                               alpha=0.75, ax=ax)
                for i, tier in enumerate(existing_tiers):
                    med = f12_data.loc[f12_data[CFG.COL_TIER] == tier, _qual_col].median()
                    if not np.isnan(med):
                        # Red median line across the violin width
                        ax.hlines(med, i - 0.22, i + 0.22, colors=CFG.VIS_ACCENT["alert"],
                                  linewidth=2.0, zorder=7, linestyle="-")
                        ax.text(i + 0.25, med, "  " + _qual_fmt.format(med), va="center", ha="left",
                                fontsize=7, color=CFG.VIS_ACCENT["alert"], fontweight="bold", zorder=8)
                # Top extends to 1.0 (0–1 composite) so the upper violin tails are
                # not clipped; floor keeps the zoomed lower bound.
                _qtop12 = min(1.0, max(_qhi12 + _qpad12, float(_qv12.max()) + _qpad12))
                if _qtop12 >= 0.9:
                    _qtop12 = 1.0
                ax.set_ylim(max(0.0, _qlo12 - _qpad12), _qtop12)
                ax.set_ylabel(_qual_lbl, )
            else:
                ax.set_visible(False)

            # Right axis: per-tier mean active-site contact density connected dot-line
            if _has_dens:
                f20_df = df.dropna(subset=["Interaction_Density_Norm"]).copy()
                valid_tiers_f20 = [t for t in existing_tiers if t in f20_df[CFG.COL_TIER].values]
                pip_xs, pip_meds, pip_lo_ci, pip_hi_ci = [], [], [], []
                for tier in valid_tiers_f20:
                    x_pos = existing_tiers.index(tier) if tier in existing_tiers else 0
                    vals  = f20_df.loc[f20_df[CFG.COL_TIER] == tier, "Interaction_Density_Norm"].dropna()
                    if len(vals) < 2:
                        continue
                    mean_v = float(vals.mean())
                    ci95  = float(_sc_stats12.sem(vals) * _sc_stats12.t.ppf(0.975, len(vals) - 1))
                    pip_xs.append(x_pos)
                    pip_meds.append(mean_v)
                    pip_lo_ci.append(mean_v - ci95)
                    pip_hi_ci.append(mean_v + ci95)
                if pip_xs:
                    ax_r.plot(pip_xs, pip_meds, color=CFG.VIS_ACCENT["vermillion"], linewidth=2.5,
                              marker="D", markersize=9, markeredgecolor="black",
                              markeredgewidth=0.9, zorder=7, label="Active-site contact density\n(mean ± 95 % CI)")
                    ax_r.fill_between(pip_xs, pip_lo_ci, pip_hi_ci, color=CFG.VIS_ACCENT["vermillion"],
                                      alpha=0.20, zorder=6)
                    # Mean annotation
                    for xi, yi in zip(pip_xs, pip_meds):
                        ax_r.text(xi, yi + (max(pip_meds) - min(pip_meds)) * 0.04,
                                  f"{yi:.2f}", ha="center", va="bottom", fontsize=8,
                                  color=CFG.VIS_ACCENT["vermillion"], fontweight="bold", zorder=8)
                ax_r.set_ylabel("Active-site contact density  (interactions per complex)",
                                 rotation=270, labelpad=14, color=CFG.VIS_ACCENT["vermillion"])
                ax_r.tick_params(axis="y", labelcolor=CFG.VIS_ACCENT["vermillion"], labelsize=9)

            ax.set_xticks(range(len(existing_tiers)))
            ax.set_xticklabels(existing_tiers, rotation=40, ha="right", fontsize=9)
            for _tick12x, _tier12x in zip(ax.get_xticklabels(), existing_tiers):
                _tick12x.set_color(TIER_PALETTE.get(_tier12x, "black"))
            # Left axis (binding prob.): steel-blue labels matching blue gridlines
            ax.tick_params(axis="y", labelsize=9, labelcolor=CFG.VIS_ACCENT["axis_left"])
            ax.yaxis.label.set_color(CFG.VIS_INK["near_black"])
            # Right axis (contact density): amber labels matching amber gridlines
            ax_r.tick_params(axis="y", labelsize=9, labelcolor=CFG.VIS_ACCENT["axis_right"])
            ax_r.yaxis.label.set_color(CFG.VIS_INK["near_black"])
            ax.set_xlabel("Degrader Tier", labelpad=-18)   # lift into the whitespace above the angled tier ticks

            # Dual-colour grids: left-axis = steel-blue; right-axis = amber
            ax.set_axisbelow(True)
            ax.yaxis.grid(True, color=CFG.VIS_ACCENT["axis_left"], linewidth=0.55, linestyle="-", alpha=0.20, zorder=0)
            ax.xaxis.grid(False)
            ax_r.yaxis.grid(True, color=CFG.VIS_ACCENT["axis_right"], linewidth=0.55, linestyle="--", alpha=0.18, zorder=0)
            ax_r.set_axisbelow(True)

            # Legend — bottom-left (violins tend to be taller on right)
            if _has_dens and pip_xs:
                from matplotlib.lines import Line2D as _Line12
                _pip_handle = _Line12([0], [0], color=CFG.VIS_ACCENT["vermillion"], linewidth=2.5,
                                      marker="D", markersize=8,
                                      label="Active-site contact density  (mean ± 95% CI)")
                _leg12 = ax_r.legend(handles=[_pip_handle], loc="lower left",
                                       fancybox=True)
                _leg12.set_zorder(20)

            plt.tight_layout()
            _tier_seps(plt.gca())   # consistent tier separators
            plt.savefig(out_dir / "Figure_16_Catalytic_Quality_vs_Inhibition.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('16')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 17: Active-site contact density by tier ---
    if "Interaction_Density_Norm" in df.columns and CFG.COL_TIER in df.columns:
        try:
            f20_df = df.dropna(subset=["Interaction_Density_Norm"]).copy()
            valid_t20 = [t for t in existing_tiers if t in f20_df[CFG.COL_TIER].values]
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
                ax.axhspan(_z_hi20, _d_max, alpha=0.07, color=CFG.VIS_ACCENT["green"], zorder=0)
                ax.axhspan(_z_lo20, _z_hi20, alpha=0.06, color=CFG.VIS_ACCENT["amber"], zorder=0)
                ax.axhspan(0.0, _z_lo20, alpha=0.06, color=CFG.VIS_ACCENT["vermillion"], zorder=0)
                ax.axhline(_z_hi20, color=CFG.VIS_ACCENT["green"], linestyle="--", alpha=0.55, linewidth=1.0)
                ax.axhline(_z_lo20, color=CFG.VIS_ACCENT["vermillion"], linestyle=":", alpha=0.55, linewidth=1.0)
                sns.boxplot(data=f20_df, x=CFG.COL_TIER, y="Interaction_Density_Norm",
                            order=valid_t20, palette=TIER_PALETTE, linewidth=1.2,
                            showfliers=False, ax=ax)
                sns.stripplot(data=f20_df, x=CFG.COL_TIER, y="Interaction_Density_Norm",
                              order=valid_t20, color="black", alpha=0.10, size=2.0, jitter=True, ax=ax)
                # Per-tier mean trend line (diamond markers) over the median boxes.
                _mean20_xs, _mean20_ys = [], []
                for i, tier in enumerate(valid_t20):
                    sub20 = f20_df.loc[f20_df[CFG.COL_TIER] == tier, "Interaction_Density_Norm"].dropna()
                    if len(sub20) == 0:
                        continue
                    _mean20_xs.append(i)
                    _mean20_ys.append(float(sub20.mean()))
                if len(_mean20_xs) >= 2:
                    ax.plot(_mean20_xs, _mean20_ys, color=CFG.VIS_ACCENT_DEEP["purple_deep"], linewidth=2.0,
                            linestyle="--", alpha=0.9, zorder=6, marker="D", markersize=6,
                            markeredgecolor="white", markeredgewidth=0.7,
                            label="Mean trend")
                for i, tier in enumerate(valid_t20):
                    sub20 = f20_df.loc[f20_df[CFG.COL_TIER] == tier, "Interaction_Density_Norm"].dropna()
                    if len(sub20) == 0:
                        continue
                    med = float(sub20.median())
                    _bg_col20 = TIER_PALETTE.get(tier, CFG.VIS_INK["mid"])
                    # Adaptive text: white on dark backgrounds, near-black on light ones
                    _txt_col20 = _text_color(_bg_col20)
                    ax.text(i, med, f"{med:.1f}", ha="center", va="center",
                            fontsize=8.5, fontweight="bold", zorder=7,
                            color=_txt_col20,
                            bbox=dict(boxstyle="round,pad=0.15", fc=_bg_col20,
                                      ec="white", linewidth=0.6, alpha=0.88))
                ax.set_xlabel("Degrader Tier", labelpad=-18)   # lift into the whitespace above the angled tier ticks
                ax.set_ylabel("Active-site contact density  (interactions per complex)", )
                ax.set_xticks(range(len(valid_t20)))
                ax.set_xticklabels(valid_t20, rotation=35, ha="right", fontsize=9)
                for _tick20x, _tier20x in zip(ax.get_xticklabels(), valid_t20):
                    _tick20x.set_color(TIER_PALETTE.get(_tier20x, "black"))
                ax.set_axisbelow(True)
                ax.yaxis.grid(True, color=CFG.VIS_ACCENT["axis_left"], linewidth=0.55, linestyle="-",
                              alpha=0.22, zorder=0)
                ax.xaxis.grid(False)
                ax.tick_params(axis="y", labelcolor=CFG.VIS_ACCENT["axis_left"], labelsize=9)
                ax.yaxis.label.set_color(CFG.VIS_INK["near_black"])
                ax.set_ylim(bottom=0, top=_d_max)
                # Zone labels parked at the right margin.
                for _zy20, _zt20, _zc20 in [((_z_hi20 + _d_max) / 2, "High", CFG.VIS_BAND["high"]),
                                            ((_z_lo20 + _z_hi20) / 2, "Moderate", CFG.VIS_BAND["moderate"]),
                                            (_z_lo20 / 2, "Low", CFG.VIS_RAMP["orange"][3])]:
                    ax.text(0.995, _zy20, _zt20, transform=ax.get_yaxis_transform(),
                            ha="right", va="center", fontsize=_JF_FA, color=_zc20,
                            fontweight="bold", alpha=0.85, zorder=8)
                # Stats notes go bottom-left: the top-left corner carries the highest boxes (Tier_1A
                # and Tier_1B), so text placed there sits on the data it is describing.
                if _n_hidden20 > 0:
                    ax.text(0.01, 0.055,
                            f"{_n_hidden20:,} complexes > {_d_max:.1f} hidden (axis clamped at p99)",
                            transform=ax.transAxes, ha="left", va="bottom",
                            fontsize=7, color=CFG.VIS_INK["grey"], style="italic")
                ax.legend(loc="upper right")

                # Kruskal-Wallis significance note
                try:
                    from scipy.stats import kruskal as _kw20
                    _kw_groups = [f20_df.loc[f20_df[CFG.COL_TIER] == t,
                                             "Interaction_Density_Norm"].dropna().values
                                  for t in valid_t20 if len(f20_df[f20_df[CFG.COL_TIER]==t]) >= 3]
                    if len(_kw_groups) >= 2:
                        _kw_h, _kw_p = _kw20(*_kw_groups)
                        # Same family as the other tests: it must be corrected WITH them, or the
                        # FDR control is applied to an incomplete set.
                        _register_p("Kruskal-Wallis across tiers", "Interaction_Density_Norm",
                                    float(_kw_h), int(sum(len(g) for g in _kw_groups)), float(_kw_p))
                        """
                        The annotation is marked UNCORRECTED. The Benjamini-Hochberg q depends on the
                        whole test family, which is complete only after every figure has run, so no q
                        exists at draw time — and printing a bare p next to a significance claim reads
                        as though the multiplicity had been accounted for. The corrected q_BH for this
                        exact test is in 06_Statistical_Tests.csv, which is the reportable value.
                        """
                        _pstr = ("p<0.001" if _kw_p < 0.001 else f"p={_kw_p:.3f}")
                        ax.text(0.01, 0.015,
                                f"Kruskal–Wallis {_pstr} uncorrected  (across tiers; "
                                f"BH q in 06_Statistical_Tests.csv)",
                                transform=ax.transAxes, ha="left", va="bottom",
                                fontsize=8, color=CFG.VIS_INK["dark"], style="italic")
                except Exception:
                    pass

                plt.tight_layout()
                _tier_seps(plt.gca())   # consistent vertical tier separators
                plt.savefig(out_dir / "Figure_17_ActiveSite_Contact_Density_by_Tier.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
                plt.close()
            else:
                reporter.log(f"  ! Skipped: {_fig_path('17')} — no valid tier groups for Interaction_Density_Norm.")
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('17')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig
    else:
        reporter.log(f"  ! Skipped: {_fig_path('17')} — column 'Interaction_Density_Norm' or 'degrader_tier' absent.")

    # --- Figure 17b: Binding Energetics — Binding Probability (violin) per tier ---
    """
    Per-tier read-out of binding quality. The Binding_Probability distribution is shown
    as a violin (zoomed to where values cluster) with a red median bar per tier. Product
    inhibition is not plotted here; it is assessed downstream by the Step-06 MM-GBSA stage.
    """
    _be_bind = "Binding_Probability" if "Binding_Probability" in df.columns else None
    _be_has_bind = _be_bind is not None and CFG.COL_TIER in df.columns
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
            sns.violinplot(data=_bedata, x=CFG.COL_TIER, y=_be_bind, order=existing_tiers,
                           palette=TIER_PALETTE, inner="box", linewidth=1.2, cut=0,
                           alpha=0.75, ax=ax)
            _be_mxs, _be_mys = [], []
            for _i, _t in enumerate(existing_tiers):
                _m = _bedata.loc[_bedata[CFG.COL_TIER] == _t, _be_bind].median()
                if not np.isnan(_m):
                    _be_mxs.append(_i); _be_mys.append(float(_m))
            # Single connected median trend (diamonds joined) — no per-violin bars.
            if len(_be_mxs) >= 2:
                ax.plot(_be_mxs, _be_mys, color=CFG.VIS_ACCENT["alert"], linewidth=2.2, linestyle="--",
                        alpha=0.9, marker="D", markersize=7, markeredgecolor="black",
                        markeredgewidth=0.7, zorder=8, label="Median trend")
                # Boxed value labels just below each node (Fig 03 style).
                _be_span = (_be_hi - _be_lo) if _be_hi > _be_lo else 1.0
                for _mx, _my in zip(_be_mxs, _be_mys):
                    ax.text(_mx, _my - _be_span * 0.03, f"{_my:.3f}", ha="center", va="top",
                            fontsize=7.0, fontweight="bold", color=CFG.VIS_ACCENT["alert"], zorder=9,
                            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec=CFG.VIS_INK["palest"],
                                      linewidth=0.5, alpha=0.85))
            if _be_hi > _be_lo:
                ax.set_ylim(_be_lo, _be_hi)
            ax.set_ylabel("Binding probability  (sigmoid of interaction density, cross-PAE, confidence)",  color=CFG.VIS_ACCENT["axis_left"])
            ax.tick_params(axis="y", labelcolor=CFG.VIS_ACCENT["axis_left"], labelsize=9)
            ax.set_xticks(range(len(existing_tiers)))
            ax.set_xticklabels(existing_tiers, rotation=40, ha="right", fontsize=9)
            for _txk, _tt in zip(ax.get_xticklabels(), existing_tiers):
                _txk.set_color(TIER_PALETTE.get(_tt, "black"))
            ax.set_xlabel("Degrader tier", labelpad=-18)   # lift into the whitespace above the angled tier ticks
            ax.set_axisbelow(True)
            ax.yaxis.grid(True, color=CFG.VIS_ACCENT["axis_left"], linewidth=0.55, alpha=0.20, zorder=0)
            ax.xaxis.grid(False)
            from matplotlib.patches import Patch as _Patch17b
            from matplotlib.lines import Line2D as _Line17b
            _lh17b = [
                _Patch17b(facecolor=CFG.VIS_RAMP["blue"][1], edgecolor=CFG.VIS_ACCENT["axis_left"],
                          label="Binding probability (per-tier violin)"),
                _Line17b([0], [0], color=CFG.VIS_ACCENT["alert"], linewidth=2.2,
                         label="Binding probability median"),
            ]
            # Single row, anchored top-left above the axes (keeps it off the violins).
            ax.legend(handles=_lh17b, loc="lower left", bbox_to_anchor=(0.0, 1.01),
                      ncol=len(_lh17b),   fancybox=True).set_zorder(20)
            plt.tight_layout()
            _stat_box(ax, _kruskal_by_tier(_bedata, _be_bind, tiers=existing_tiers), "lower left")
            _tier_seps(plt.gca())   # consistent tier separators
            plt.savefig(out_dir / "Figure_17b_Binding_Energetics.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('17b')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig
        finally:
            if fig is not None:
                plt.close(fig)

    # --- Figure 18a: Chemical Space UMAP Manifold ---
    if "UMAP_X" in df.columns:
        # Degradability landscape: a hexbin over the UMAP chemical space coloured by the MEAN competence
        # per bin, so chemical regions enriched for good degraders light up; the MD-ready / elite hits are
        # overlaid as stars, each with its PyMOL active-site thumbnail in the right margin. The main axes
        # sit in the middle band and the right margin holds the thumbnails — no tight_layout, it would
        # fight the margins.
        fig, ax = plt.subplots(figsize=(11, 7.5))
        fig.subplots_adjust(**_TT_MARGINS)
        _u = df.dropna(subset=["UMAP_X", "UMAP_Y"]).copy()
        _ux, _uy = _u["UMAP_X"].values, _u["UMAP_Y"].values
        _cmetric = next((c for c in ("competence_score", CFG.COL_MECH_S,
                                     CFG.COL_CONF) if c in _u.columns), None)
        if _cmetric:
            _cvals = pd.to_numeric(_u[_cmetric], errors="coerce").fillna(0.0).values
            hb = ax.hexbin(_ux, _uy, C=_cvals, reduce_C_function=np.mean,
                           gridsize=45, cmap="viridis", mincnt=1, linewidths=0.15)
            cb = fig.colorbar(hb, ax=ax, fraction=0.025, pad=0.01)
            cb.set_label(f"Mean {_cmetric.replace('_', ' ')} per bin  (brighter = more degradable)",
                         fontsize=9)
        else:
            ax.hexbin(_ux, _uy, gridsize=45, cmap="Greys", mincnt=1)

        _md_col = getattr(CFG, "MD_SELECTED_COL", "MD_Selected")
        if _md_col in _u.columns and _u[_md_col].astype(str).str.lower().isin(["true", "1", "1.0"]).any():
            _elite = _u[_u[_md_col].astype(str).str.lower().isin(["true", "1", "1.0"])]
            _elite_lbl = f"MD-ready hits (n={len(_elite)})"
        else:
            _elite = _u[_u[CFG.COL_TIER] == CFG.TIER_TOP]
            _elite_lbl = f"{CFG.TIER_TOP} hits (n={len(_elite)})"
        # Build the thumbnail cohort ONCE (same order _tt_draw_thumbnails uses) and map each entry to its
        # box/arrow colour, so a star can be tinted to match its OWN thumbnail rather than every star being
        # an indistinguishable gold. Matched by (protein, ligand) identity.
        _pax = None
        _thumb_colour = {}
        if _tt_has_imgs and _imgs and _pa is not None and not _pa.empty:
            _pax = _pa.dropna(subset=["UMAP_X", "UMAP_Y"]).copy()
            if not _pax.empty:
                _pax["_X"] = _pax["UMAP_X"]; _pax["_Y"] = _pax["UMAP_Y"]
                _n_thumb = min(len(_pax), len(_imgs), CFG.VIS_MAX_THUMBNAILS)
                for _ti, (_, _tr) in enumerate(_pax.head(_n_thumb).iterrows()):
                    _tk = (str(_tr.get(CFG.COL_PROT, "")), str(_tr.get(CFG.COL_LIG, "")))
                    _thumb_colour[_tk] = _TT_ENTRY_COLS[_ti % len(_TT_ENTRY_COLS)]

        if not _elite.empty:
            # Each star takes its thumbnail's colour (fill); an elite hit with no thumbnail stays gold.
            # No ligand labels on the map — the thumbnail box already names protein + ligand.
            for _, _er in _elite.iterrows():
                # The 3R3U × FA positive control is drawn as the distinct red star / gold border;
                # candidate hits keep their per-thumbnail tint (so each star matches its own box).
                if str(_er.get("is_control", "")).strip().lower() in ("true", "1", "1.0", "yes"):
                    _ctrl_star(ax, _er["UMAP_X"], _er["UMAP_Y"], size=210)
                    continue
                _ek = (str(_er.get(CFG.COL_PROT, "")), str(_er.get(CFG.COL_LIG, "")))
                _sc = _thumb_colour.get(_ek, CFG.VIS_ACCENT["star"])
                ax.scatter(_er["UMAP_X"], _er["UMAP_Y"], marker="*", s=185,
                           color=_sc, edgecolors="black", linewidths=0.9, zorder=8)
            ax.scatter([], [], marker="*", s=185, color=CFG.VIS_ACCENT["star"],
                       edgecolors="black", linewidths=0.9, label=_elite_lbl)   # single neutral legend key

        ax.set_xlabel("UMAP Dimension 1  (distances reflect chemical similarity)", )
        ax.set_ylabel("UMAP Dimension 2", )
        ax.grid(False)
        ax.legend(loc="lower left")
        _stat_box(ax, _umap_tier_separation(_u), "upper left")
        # PyMOL active-site thumbnails for the Tier_1A representatives: each thumbnail's arrow and box
        # share the colour of the star it points to.
        if _pax is not None and not _pax.empty:
            _tt_draw_thumbnails(fig, ax, _pax, _imgs)
        plt.savefig(out_dir / "Figure_18a_Chemical_Space_Map.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()



# =============================================================================
# STEP 7/8 — 07_PFAS_Scope_and_Synthesis
# =============================================================================
def _fig_folder07_pfas(df, features, out_dir, reporter, existing_tiers):
    """Folder 07_PFAS_Scope_and_Synthesis — synthesis + publication assembly."""
    reporter.section("Step 7/8 — 07_PFAS_Scope_and_Synthesis · synthesis + publication assembly")
    # --- Figure 19: Candidate Radar: 18a (top-5 hits) + 18b (one per tier) ---
    _radar_labels = {
        CFG.COL_CONF:   "AI Conf.",
        "iptm":                     "ipTM",
        "Interaction_Density_Norm": "Int.Den",
        "mean_plddt":               "pLDDT",
        "Binding_Probability":      "Bind.Prob",
        CFG.COL_SN2:         "SN2(°)",
        CFG.COL_MECH_S: "Mech",
        "Active_Site_RMSD":         "RMSD(Å)",
        CFG.COL_ID_PCT:             "Seq.ID%",
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
        colours_r = [CFG.VIS_BAND["high"], CFG.VIS_RADAR_SERIES[1], CFG.VIS_ACCENT["amber"], CFG.VIS_ACCENT["magenta"], CFG.VIS_ACCENT["blue"],
                     CFG.VIS_ACCENT["sky"], CFG.VIS_ACCENT["yellow"], CFG.VIS_ACCENT["green"], CFG.VIS_ACCENT["vermillion"], CFG.VIS_ACCENT["magenta"]]
        lname_col = next((c for c in [CFG.COL_LIG, "ligand"] if c in rows_df.columns), None)

        _ring_levels = [0.2, 0.4, 0.6, 0.8, 1.0]
        _ring_angles = np.linspace(0, 2 * np.pi, 300)
        for _rl in _ring_levels:
            ax_r.plot(_ring_angles, [_rl] * 300,
                      color=CFG.VIS_BAND["low"], linewidth=0.55, alpha=0.45, zorder=1,
                      linestyle="-", solid_capstyle="round")

        import re as _re18
        for i in range(len(rows_df)):
            vals  = _scaler_vals.iloc[i].values.flatten().tolist() + [_scaler_vals.iloc[i].values[0]]
            _raw  = str(rows_df.iloc[i][lname_col]) if lname_col else f"Hit {i+1}"
            # Canonical short label (FA/DFA/TFA, "NN_" prefix stripped) — matches every other figure.
            _clean = _lig_short(_raw)
            tier  = str(rows_df.iloc[i].get(CFG.COL_TIER, ""))
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
            ax_r.plot(angles, bvals, linewidth=2, linestyle="--", color=CFG.VIS_INK["faint"],
                      label=f"{CFG.TIER_DECOY} avg.  (reference)", zorder=3)
            ax_r.fill(angles, bvals, color=CFG.VIS_INK["faint"], alpha=0.05, zorder=2)
            ax_r.scatter(angles[:-1], bvals[:-1],
                         color=CFG.VIS_INK["faint"], s=45, zorder=5, edgecolors="white", linewidths=0.8)

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
                             color=CFG.VIS_BAND["low"], size=7.5, fontweight="bold")
        ax_r.set_rlabel_position(45)
        ax_r.set_ylim(0, 1.10)
        ax_r.yaxis.grid(False)
        ax_r.xaxis.grid(True, color=CFG.VIS_INK["palest"], linewidth=0.6, alpha=0.7)

        n_lines = len(rows_df) + (0 if baseline_df.empty else 1)
        ax_r.legend(loc="lower center", bbox_to_anchor=(0.5, -0.18),
                    ncol=min(3, n_lines), fontsize=_radar_fs,
                     fancybox=True)
        plt.tight_layout(rect=[0, 0.15, 1, 1])
        plt.savefig(fig_path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close()

    try:
        metrics_r = [m for m in _radar_labels if m in df.columns]
        if metrics_r:
            worst_tier = existing_tiers[-1] if existing_tiers else None
            worst_avg  = (df[df[CFG.COL_TIER] == worst_tier].mean(numeric_only=True).to_frame().T
                          if worst_tier else pd.DataFrame())

            """
            The radar profiles the MD-READY COHORT, not the whole elite tier.

            A radar is a per-complex figure: every hit is one polygon, and the reader compares them
            spoke by spoke. That only works for a handful of complexes. Tier_1A now holds ~100
            members, so drawing the tier would either overplot into an unreadable web or — worse —
            silently truncate to an arbitrary top-N and present it as though it were the tier.

            The MD-selected cohort is the set this figure is actually about: the few complexes the
            pipeline commits to simulate, one per ligand. Those are the polygons worth comparing.
            Falls back to the top Tier_1A ranks only when no MD selection exists yet (e.g. a run
            stopped before Step 02 wrote the column).
            """
            _pa18a = _md_ready_df(df).copy()
            _radar_src = "MD-selected cohort"
            if _pa18a.empty:
                _pa18a = df[df[CFG.COL_TIER] == CFG.TIER_TOP].copy()
                _radar_src = f"top {CFG.TIER_TOP} ranks (no MD selection in this run)"
            _sort_col18a = "Scientific_Rank" if "Scientific_Rank" in _pa18a.columns \
                           else CFG.COL_CONF
            _asc18a = _sort_col18a == "Scientific_Rank"
            _pa18a  = _pa18a.sort_values(_sort_col18a, ascending=_asc18a)
            if _pa18a.empty:
                _pa18a = df.sort_values(CFG.COL_CONF, ascending=False).head(8)
            # Cap plotted hits so the radar stays readable (top-N by ranking key).
            _pa18a = _pa18a.head(CFG.VIS_RADAR_MAX_HITS)
            _draw_radar(_pa18a, worst_avg, metrics_r,
                        out_dir / "Figure_19a_Radar_TopHits.png",
                        f"{_radar_src} vs. worst-tier baseline",
                        label_mode="rank")

            # --- 18b: One representative per tier (best Boltz confidence per tier) ---
            reps_b = []
            for t in existing_tiers:
                sub_t = df[df[CFG.COL_TIER] == t]
                if not sub_t.empty:
                    reps_b.append(sub_t.sort_values(CFG.COL_CONF, ascending=False).iloc[[0]])
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
    For each tier: the % meeting SN2 ≥ 165° (geometry) against the % meeting
    Boltz_Model_Confidence ≥ 0.75 (AI), which places each tier in the biologically
    meaningful substrate-vs-inhibitor space. Ranking on Binding_Probability cannot
    do this — it is ~constant at 0.99 across the whole dataset and separates nothing.
    """
    if CFG.COL_SN2 in df.columns and CFG.COL_CONF in df.columns:
        try:

            _d19 = df[[CFG.COL_SN2, CFG.COL_CONF, CFG.COL_TIER]].copy()
            _d19[CFG.COL_SN2]      = pd.to_numeric(_d19[CFG.COL_SN2],      errors="coerce")
            _d19[CFG.COL_CONF] = pd.to_numeric(_d19[CFG.COL_CONF], errors="coerce")
            _d19 = _d19.dropna()
            _tiers19 = [t for t in TIER_ORDER_LOGIC if t in _d19[CFG.COL_TIER].unique()]

            # Per-tier metrics
            _rows19 = []
            for _t19 in _tiers19:
                _sub = _d19[_d19[CFG.COL_TIER] == _t19]
                _n   = len(_sub)
                _pct_sn2  = float((_sub[CFG.COL_SN2] >= CFG.SUBSTRATE_ANGLE_MIN).mean() * 100)
                _pct_conf = float((_sub[CFG.COL_CONF] >= CFG.SUBSTRATE_CONF_MIN).mean() * 100)
                _pct_both = float(((_sub[CFG.COL_SN2] >= CFG.SUBSTRATE_ANGLE_MIN) & (_sub[CFG.COL_CONF] >= CFG.SUBSTRATE_CONF_MIN)).mean() * 100)
                _pct_inh  = float((_sub[CFG.COL_SN2] < CFG.INHIBITOR_ANGLE_MAX).mean() * 100)
                _med_sn2  = float(_sub[CFG.COL_SN2].median())
                _med_conf = float(_sub[CFG.COL_CONF].median())
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
            _master_ec19 = [TIER_PALETTE.get(t, CFG.VIS_INK["silver"]) for t in _df19["tier"].tolist()]
            """
            Faint tint of each tier's own colour as the track fill (alpha baked in
            so the coloured outline stays vivid) — cohesive, not flat grey.
            """
            _master_fc19 = [_to_rgba19(c, 0.13) for c in _master_ec19]
            ax19L.barh(_y19, 105, height=0.82, color=_master_fc19,
                       edgecolor=_master_ec19, linewidth=1.4, zorder=2)
            ax19L.barh(_y19 + _bar_h,  _df19["pct_substrate"].values,  height=_bar_h,
                       color=CFG.VIS_ACCENT["green"], alpha=0.85, label=f"Substrate (SN2≥{CFG.SUBSTRATE_ANGLE_MIN:.0f}° + Conf≥{CFG.SUBSTRATE_CONF_MIN:.2f})", zorder=3)
            ax19L.barh(_y19,            _df19["pct_sn2_ok"].values,     height=_bar_h,
                       color=CFG.VIS_ACCENT["sky"], alpha=0.75, label=f"SN2≥{CFG.SUBSTRATE_ANGLE_MIN:.0f}° (any conf)", zorder=3)
            ax19L.barh(_y19 - _bar_h,  _df19["pct_inhibitor"].values,  height=_bar_h,
                       color=CFG.VIS_ACCENT["vermillion"], alpha=0.80, label=f"Potential inhibitor (SN2<{CFG.INHIBITOR_ANGLE_MAX:.0f}°)", zorder=3)

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
                                   va="center", ha="center", fontsize=7.0,
                                   color="white", fontweight="bold", zorder=5)
                    else:
                        ax19L.text(_pct_val + 1.0, _yi + _yo, f"{_pct_val:.0f}%",
                                   va="center", ha="left", fontsize=7.0,
                                   color=CFG.VIS_INK["dark"], zorder=5)

            ax19L.set_yticks(_y19)
            ax19L.set_yticklabels(
                [f'{r["tier"]}  (n={r["n"]:,})' for _, r in _df19.iterrows()],
                fontsize=7
            )
            for tick, tier in zip(ax19L.get_yticklabels(), _df19["tier"].tolist()):
                tick.set_color(TIER_PALETTE.get(tier, "black"))
            ax19L.set_xlim(0, 105)
            ax19L.axvline(25, color=CFG.VIS_INK["faint"], lw=0.7, ls=":", zorder=2.4)
            ax19L.axvline(50, color=CFG.VIS_INK["faint"], lw=0.7, ls=":", zorder=2.4)
            ax19L.axvline(75, color=CFG.VIS_INK["faint"], lw=0.7, ls=":", zorder=2.4)
            ax19L.set_xlabel("Percentage of complexes (%)", )
            ax19L.tick_params(axis="x", labelsize=8)
            # Legend: single row, 70% font (7.5→5.25), inside bottom-right
            ax19L.legend(loc="lower right", ncol=3,
                          fancybox=True, fontsize=CFG.VIS_FONT_LEGEND - 1.5)   # shrink to fit within the panel width
            ax19L.invert_yaxis()

            # ── Right panel: scatter of tier medians in Conf × SN2 space ─────
            # Background hexbin of all complexes (shows density landscape)
            ax19R.hexbin(_d19[CFG.COL_CONF], _d19[CFG.COL_SN2],
                                  gridsize=40, cmap="Greys", mincnt=1, alpha=0.4, linewidths=0.2, zorder=1)

            # Substrate / inhibitor zone shading — bounds from CFG
            ax19R.axhspan(CFG.SUBSTRATE_ANGLE_MIN, 185, color=CFG.VIS_ACCENT["green"], alpha=0.10, zorder=0)
            ax19R.axhspan(0, CFG.INHIBITOR_ANGLE_MAX, color=CFG.VIS_ACCENT["vermillion"], alpha=0.08, zorder=0)
            ax19R.axhline(CFG.SUBSTRATE_ANGLE_MIN, color=CFG.VIS_ACCENT["green"], lw=1.2, ls="--", alpha=0.7, zorder=2)
            ax19R.axhline(CFG.INHIBITOR_ANGLE_MAX, color=CFG.VIS_ACCENT["vermillion"], lw=1.2, ls="--", alpha=0.7, zorder=2)
            ax19R.axvline(CFG.SUBSTRATE_CONF_MIN, color=CFG.VIS_ACCENT["blue"], lw=1.2, ls="--", alpha=0.7, zorder=2)

            # Tier median diamonds with IQR error bars — collect positions first
            _tier_pts19 = []   # (tier, cx, cy, colour)
            for _t19 in _tiers19:
                _sub = _d19[_d19[CFG.COL_TIER] == _t19]
                if _sub.empty:
                    continue
                _cx = float(_sub[CFG.COL_CONF].median())
                _cy = float(_sub[CFG.COL_SN2].median())
                _ex_lo = _cx - float(_sub[CFG.COL_CONF].quantile(0.25))
                _ex_hi = float(_sub[CFG.COL_CONF].quantile(0.75)) - _cx
                _ey_lo = _cy - float(_sub[CFG.COL_SN2].quantile(0.25))
                _ey_hi = float(_sub[CFG.COL_SN2].quantile(0.75)) - _cy
                _col19 = TIER_PALETTE.get(_t19, CFG.VIS_INK["faint"])
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
                               fontsize=7.0, color=_col19, fontweight="bold",
                               arrowprops=dict(arrowstyle="-|>", color=_col19, lw=1.3,
                                               mutation_scale=15, shrinkA=3, shrinkB=5,
                                               connectionstyle="arc3,rad=0.12"),
                               bbox=dict(boxstyle="round,pad=0.18", fc="white",
                                         ec=_col19, alpha=0.92, linewidth=0.8),
                               zorder=8)

            ax19R.set_xlabel("Boltz Model Confidence", )
            ax19R.set_ylabel("SN2 Attack Angle (°)", )
            ax19R.tick_params(axis="both", labelsize=8)

            # Zone / threshold labels — all inside the axes, single line each
            ax19R.text(0.87, 181, "Substrate zone", color=CFG.VIS_BAND["high"], fontsize=7,
                       style="italic", ha="center", va="center")
            ax19R.text(0.87, 8, "Inhibitor zone", color=CFG.VIS_ACCENT_DEEP["orange_hot"], fontsize=7,
                       style="italic", ha="center", va="center")
            ax19R.text(0.745, 95, f"Conf threshold = {CFG.SUBSTRATE_CONF_MIN:.2f}", color=CFG.VIS_ACCENT["blue"], fontsize=7.0,
                       ha="center", va="center", style="italic", rotation=90)

            # Legend: 3 columns, shrunk icons, fixed outside above the top-left corner
            ax19R.legend(loc="lower left", bbox_to_anchor=(0.0, 1.02),
                         ncol=2,   fancybox=True,
                         markerscale=0.45)

            fig19a.savefig(out_dir / "Figure_20a_Tier_Success_Rates.png",
                           dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            fig19b.savefig(out_dir / "Figure_20b_Conf_SN2_Landscape.png",
                           dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close(fig19a)
            plt.close(fig19b)
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('20')} — {e}")
            plt.close("all")   # release the figure left open by the failed savefig

    # --- Figure 21: Conflict Composition — 100% Stacked Horizontal Bar per Tier ---
    """
    One bar per tier (best → worst on y-axis).  Left portion = favourable categories
    (Consensus High, Hidden Gem); right portion = concerning categories.
    A thin vertical divider at the "favourable total" position separates the two sides.
    Segment labels show % only for segments ≥ 5%.  Y-axis labels include tier n-count.
    """
    if "Conflict_Category" in df.columns and CFG.COL_TIER in df.columns:
        cat_colours_f2 = {
            "Consensus High": CFG.VIS_ACCENT["green"], "Hidden Gem": CFG.VIS_ACCENT["magenta"],
            CFG.TIER_DECOY: CFG.VIS_ACCENT["vermillion"], "Consensus Low": CFG.VIS_INK["grey"], "Ambiguous": CFG.VIS_ACCENT["amber"]}
        # Canonical plotting order: favourable → concerning, left → right
        cat_order_f2 = ["Consensus High", "Hidden Gem", "Ambiguous", "Consensus Low", "Decoy"]
        present_cats_f2 = [c for c in cat_order_f2 if c in df["Conflict_Category"].unique()]

        ct_f2 = pd.crosstab(df[CFG.COL_TIER], df["Conflict_Category"])
        tier_rows_f2 = [t for t in existing_tiers if t in ct_f2.index]
        ct_f2 = ct_f2.reindex(index=tier_rows_f2, columns=present_cats_f2, fill_value=0)
        tier_n_f2   = ct_f2.sum(axis=1)
        # Grand total = all rows with a valid degrader_tier (not filtered by Conflict_Category)
        _total_rows_f2 = len(df)
        grand_total    = int(df[CFG.COL_TIER].notna().sum())
        pct_f2         = ct_f2.div(tier_n_f2, axis=0).fillna(0) * 100
        # Detect control cases — generic mask first, then specific 3R3U / DeHa4
        _ctrl_mask = pd.Series(False, index=df.index)
        for _ctrl_col in ["is_control", "control", "Control", "is_Control"]:
            if _ctrl_col in df.columns:
                _ctrl_mask = df[_ctrl_col].astype(bool)
                break
        if not _ctrl_mask.any() and "job_name" in df.columns:
            _ctrl_mask = df["job_name"].str.lower().str.contains("control|ctrl|deha4_ref|3r3u", na=False)
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

        # Y positions (top tier at top)
        tier_y_f2   = np.arange(len(tier_rows_f2) - 1, -1, -1)
        bar_h       = 0.60

        # Height: fit bars tightly — no whitespace padding below lowest tier
        _fh_f2 = max(4.0, len(tier_rows_f2) * 0.65 + 1.8)
        fig_f2, ax_f2 = plt.subplots(figsize=(13, _fh_f2))
        ax_f2.set_facecolor(CFG.VIS_INK["panel"])

        left_f2 = np.zeros(len(tier_rows_f2))
        for cat in present_cats_f2:
            vals   = pct_f2[cat].values
            counts = ct_f2[cat].values
            ax_f2.barh(tier_y_f2, vals, left=left_f2, height=bar_h,
                       color=cat_colours_f2.get(cat, CFG.VIS_INK["faint"]),
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
        ax_f2.axvline(50, color=CFG.VIS_INK["muted"], linewidth=1.1, linestyle="--", alpha=0.5, zorder=1)
        ax_f2.text(50.5, len(tier_rows_f2) - 0.5, "50%",
                   va="top", ha="left", fontsize=8, color=CFG.VIS_INK["muted"])

        # Y-axis labels: tier name + n-count in matching colour
        ax_f2.set_yticks(tier_y_f2)
        ax_f2.set_yticklabels(
            [f"{t}  (n={tier_n_f2.loc[t]:,})" for t in tier_rows_f2],
            fontsize=10)
        for tick, tier in zip(ax_f2.get_yticklabels(), tier_rows_f2):
            tick.set_color(TIER_PALETTE.get(tier, "black"))

        ax_f2.set_xlim(0, 100)
        # Clamp y-range to exactly the number of tiers — prevents blank canvas below
        ax_f2.set_ylim(-0.55, len(tier_rows_f2) - 0.45)
        ax_f2.set_xlabel("Percentage of complexes within tier  (%)", )
        ax_f2.set_xticks([0, 25, 50, 75, 100])
        ax_f2.set_xticklabels(["0%", "25%", "50%", "75%", "100%"])
        ax_f2.set_axisbelow(True)
        ax_f2.xaxis.grid(True, alpha=0.30, linestyle=":", color=CFG.VIS_INK["ghost"], zorder=0)

        # Top-right badge: total jobs | tiered count — lifted above the axes frame
        ax_f2.text(0.99, 1.06,
                   f"Total jobs: {_total_rows_f2:,}   |   Tiered: {grand_total:,}",
                   transform=ax_f2.transAxes, ha="right", va="bottom",
                   fontsize=8.5, color=CFG.VIS_INK["dark"], fontweight="bold",
                   bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=CFG.VIS_INK["palest"],
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
                       color=CFG.VIS_ACCENT["green"], fontweight="bold",
                       transform=ax_f2.get_xaxis_transform(),
                       bbox=dict(boxstyle="round,pad=0.12", fc="white", ec=CFG.CONF_BAND_COLOURS["high"],
                                 alpha=0.80, linewidth=0.5))
            ax_f2.text(mid_con, 1.015, "Concerning ◀",
                       ha="center", va="bottom", fontsize=9.5,
                       color=CFG.VIS_ACCENT["vermillion"], fontweight="bold",
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
                    col_arr = cat_colours_f2.get(cat, CFG.VIS_INK["faint"])
                    ax_f2.annotate(f"n={int(ni):,}",
                                   xy=(mid_x, yi),
                                   xytext=(mid_x, yi + 0.38),
                                   ha="center", va="bottom", fontsize=7.0,
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
            if CFG.COL_LIG in df.columns:
                _ln = str(df.at[_idx, CFG.COL_LIG])
            if (not _ln or _ln.lower() == "nan") and "job_name" in df.columns:
                _jn = str(df.at[_idx, "job_name"])
                _ln = _jn.split("Control_", 1)[1] if "Control_" in _jn else _jn
            # Canonical short label (FA/DFA/TFA, "NN_" prefix stripped) — matches every other figure.
            return _lig_short(_ln)

        _ctrl_specs_f20 = [(_ctrl_3R3U_mask, "3R3U", CFG.VIS_ACCENT_DEEP["red_deep"]),
                           (_ctrl_deha4_mask, "DeHa4", CFG.VIS_ACCENT_DEEP["blue"])]
        _ctrl_by_tier_f20 = {}   # tier -> [(label, colour), ...]
        for _cmask, _cprot, _ccol in _ctrl_specs_f20:
            if not (_cmask.any() and "job_name" in df.columns):
                continue
            for _cidx in df.index[_cmask]:
                _ctier = df.at[_cidx, CFG.COL_TIER]
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
                               ha="right", va="center", fontsize=7.0,
                               color=_ccol, fontweight="bold",
                               arrowprops=dict(arrowstyle="->", color=_ccol, lw=0.9,
                                               shrinkA=1, shrinkB=2),
                               bbox=dict(boxstyle="round,pad=0.12", fc="white",
                                         ec=_ccol, alpha=0.90, linewidth=0.5),
                               zorder=7)

        from matplotlib.patches import Patch as _P02
        _leg02 = [_P02(facecolor=cat_colours_f2.get(c, CFG.VIS_INK["faint"]), edgecolor="white",
                       linewidth=0.5, label=c) for c in present_cats_f2]
        # Single-row legend below chart — ncol = number of categories so all fit on one line
        _leg02_obj = ax_f2.legend(handles=_leg02,
                                  loc="upper center", bbox_to_anchor=(0.5, -0.10),
                                  ncol=len(present_cats_f2),
                                    fancybox=True,
                                  title="Conflict Category", )
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
    _phys_col21_abs = next((c for c in [CFG.COL_SN2, "Scientific_Rank", "Pareto_Rank"] if c in df.columns), None)
    _ai_col21_abs   = next((c for c in [CFG.COL_CONF, "Ensemble_Data_Rank"] if c in df.columns), None)
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

                lname_col21 = next((c for c in ["complex_id", CFG.COL_LIG, "ligand",
                                                CFG.COL_PROT, "protein"]
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
                ax21.set_facecolor(CFG.VIS_INK["panel"])
                ax21.set_axisbelow(True)
                ax21.xaxis.grid(True, color=CFG.VIS_INK["tick"], linewidth=0.7, zorder=0)

                # Background zone bands — tertile reference lines
                ax21.axvspan(0,  33, alpha=0.04, color=CFG.VIS_ACCENT["vermillion"], zorder=0)  # lower third
                ax21.axvspan(33, 67, alpha=0.03, color=CFG.VIS_ACCENT["amber"], zorder=0)  # middle
                ax21.axvspan(67, 100, alpha=0.04, color=CFG.VIS_ACCENT["green"], zorder=0)  # top third
                for _xv21 in [33, 67]:
                    ax21.axvline(_xv21, color=CFG.VIS_INK["palest"], linewidth=0.8,
                                 linestyle="--", zorder=1)

                for _gi21, (_, row21) in enumerate(gems_s21.iterrows()):
                    _phys21 = float(row21["_phys_pct21"])
                    _ai21   = float(row21["_ai_pct21"])
                    _gap21  = float(row21["_gap21"])
                    _y21    = _gi21

                    # Connecting line — width proportional to |gap|
                    _lw21 = np.clip(0.8 + abs(_gap21) / 25, 0.8, 3.5)
                    _lo_x21, _hi_x21 = sorted([_phys21, _ai21])
                    ax21.hlines(_y21, _lo_x21, _hi_x21,
                                colors=CFG.VIS_INK["palest"], linewidth=_lw21,
                                alpha=0.75, zorder=2)

                    # Physics dot (blue circle, red outline)
                    ax21.scatter(_phys21, _y21, color=CFG.VIS_ACCENT["blue"], s=130, zorder=4,
                                 marker="o", edgecolors="red", linewidths=1.2)
                    ax21.text(_phys21, _y21 + 0.34, f"{_phys21:.0f}%",
                              va="bottom", ha="center", fontsize=7.0, color=CFG.VIS_ACCENT["blue"],
                              fontweight="bold", zorder=5)
                    # AI dot (golden diamond, no border)
                    ax21.scatter(_ai21, _y21, color=CFG.VIS_ACCENT["amber"], s=130, zorder=4,
                                 marker="D", edgecolors="none")
                    ax21.text(_ai21, _y21 - 0.34, f"{_ai21:.0f}%",
                              va="top", ha="center", fontsize=7.0, color=CFG.VIS_ACCENT["amber"],
                              fontweight="bold", zorder=5)

                    # Value labels just outside the right dot
                    _rightmost21 = max(_phys21, _ai21)
                    ax21.text(_rightmost21 + 1.2, _y21,
                              f"Δ={abs(_gap21):.1f}%",
                              va="center", ha="left", fontsize=_JF_FA,
                              color=CFG.VIS_INK["muted"], fontweight="bold")

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
                    _t21 = (row21.get(CFG.COL_TIER, "")
                            if CFG.COL_TIER in gems.columns else "")
                    _tick21.set_color(TIER_PALETTE.get(_t21, CFG.VIS_INK["dark"]))

                ax21.set_xlim(-1, _x_hi21 + 2)
                _x21_phys_lbl = ("SN2 Attack Angle" if _phys_col21_abs == CFG.COL_SN2
                                 else "Physics Rank")
                _x21_ai_lbl   = ("Boltz Confidence" if _ai_col21_abs == CFG.COL_CONF
                                 else "AI Rank")
                ax21.set_xlabel(
                    f"Normalised Score  ({_x21_phys_lbl} / {_x21_ai_lbl})  "
                    "(0 = lowest · 100 = highest)",
                    )
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
                    (16.5,  "Lower third", CFG.VIS_ACCENT["vermillion"], 33),
                    (50.0,  "Middle third", CFG.VIS_BAND["moderate"], 67),
                    (83.5,  "Top third",   CFG.VIS_ACCENT["green"], 100),
                ]
                for _zx21, _ztxt21, _zcol21, _zmax21 in _zone_lbls21:
                    if _zx21 <= _x_hi21 + 1:
                        ax21.text(_zx21, 0.97, _ztxt21, ha="center", va="top",
                                  fontsize=_JF_FA, color=_zcol21, style="italic",
                                  transform=ax21.get_xaxis_transform())

                # Legend
                _lh21 = [
                    _L21([0],[0], marker="o", color="w", markerfacecolor=CFG.VIS_ACCENT["blue"],
                         markeredgecolor="red", markeredgewidth=1.2,
                         markersize=9, label=f"Physics score ({_x21_phys_lbl})"),
                    _L21([0],[0], marker="D", color="w", markerfacecolor=CFG.VIS_ACCENT["amber"],
                         markeredgecolor="none",
                         markersize=9, label=f"AI score ({_x21_ai_lbl})"),
                ]
                # Add tier entries if multiple tiers
                _tier_in_gems = [t for t in existing_tiers
                                 if CFG.COL_TIER in gems.columns and
                                 t in gems[CFG.COL_TIER].values]
                _lh21 += [_L21([0],[0], color=TIER_PALETTE.get(t,CFG.VIS_INK["faint"]), linewidth=3,
                                label=f'{t}  (n={len(gems[gems[CFG.COL_TIER]==t]):,})')
                          for t in _tier_in_gems]
                ax21.legend(handles=_lh21, loc="lower right",
                             fancybox=True, ncol=2,
                            title="Marker type  |  Tier (line colour)",
                            )

                # Badge: summary stats
                n_gems21 = len(gems)
                med_gap21 = float(gems["_gap21"].median())
                ax21.text(0.99, 0.99,
                          f"Hidden Gems: {n_gems21:,}\nMedian |Δ|: {abs(med_gap21):.1f}%",
                          transform=ax21.transAxes, va="top", ha="right",
                          fontsize=8.5, color=CFG.VIS_INK["dark"],
                          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec=CFG.VIS_INK["palest"],
                                    alpha=0.88, linewidth=0.5))

                plt.tight_layout()
                df.drop(columns=["_phys_pct21", "_ai_pct21"], errors="ignore", inplace=True)
                plt.savefig(out_dir / "Figure_22_Hidden_Gems_DeepDive.png",
                            dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
                plt.close(fig21)
            else:
                reporter.log(f"  ! Skipped: {_fig_path('22')} — No Hidden Gems found in dataset.")
        except Exception as e:
            reporter.log(f"  ! Skipped: {_fig_path('22')} — {e}")
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
        _phys_col22 = next((c for c in ["soft_catalytic_score", CFG.COL_MECH_S, CFG.COL_SN2] if c in df.columns), None)
        _ai_col22   = next((c for c in [CFG.COL_CONF, "iptm", "ptm"] if c in df.columns), None)
        _tier_strong22 = {CFG.TIER_TOP, CFG.TIER_ORDER[1], CFG.TIER_ORDER[2]}

        if _phys_col22 and _ai_col22 and CFG.COL_TIER in df.columns:
            _n22 = len(df)
            _phys_vals22   = pd.to_numeric(df[_phys_col22], errors="coerce")
            _phys_thresh22 = float(_phys_vals22.quantile(0.67))   # top tercile of catalytic geometry (higher = better)
            _set_A = set(df.index[_phys_vals22 >= _phys_thresh22])
            _ai_thresh22   = float(pd.to_numeric(df[_ai_col22], errors="coerce").quantile(0.67))
            _set_B = set(df.index[pd.to_numeric(df[_ai_col22], errors="coerce") >= _ai_thresh22])
            _set_C = set(df.index[df[CFG.COL_TIER].isin(_tier_strong22)])

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
                    (_cx_A, _cy_A, CFG.VIS_ACCENT["blue"]),
                    (_cx_B, _cy_B, CFG.VIS_ACCENT["vermillion"]),
                    (_cx_C, _cy_C, CFG.VIS_ACCENT["green"])]:
                ax22.add_patch(_Circ22((_cx22, _cy22), _r22,
                                       facecolor=_col22, alpha=0.20,
                                       edgecolor=_col22, linewidth=2.2, zorder=2))

            _abc_cx, _abc_cy, _r_abc = 5.0, 5.10, 0.75
            ax22.add_patch(_Circ22((_abc_cx, _abc_cy), _r_abc,
                                   facecolor=CFG.VIS_INK["pale"], alpha=0.15,
                                   edgecolor=CFG.VIS_INK["soft"], linewidth=2.0,
                                   linestyle=":", zorder=4))

            _set_PA22 = set(df.index[df[CFG.COL_TIER] == CFG.TIER_TOP])
            _n_PA22   = len(_set_PA22)
            _tt_cx22, _tt_cy22, _r_PA22, _tt_col22 = _cx_C, _cy_C, 0.50, CFG.VIS_ACCENT["magenta"]
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
                              bbox=dict(boxstyle="round,pad=0.20", fc=CFG.VIS_RAMP["purple"][0],
                                        ec=_tt_col22, alpha=0.93, linewidth=1.3),
                              arrowprops=dict(arrowstyle="->", color=_tt_col22,
                                              lw=1.1, shrinkA=0, shrinkB=3))

            # Primary set labels — uniform size, bold, distinct filled-and-outlined box
            _main_lbl_kw22 = dict(ha="center", va="center", fontsize=11,
                                  fontweight="bold", linespacing=1.25, zorder=6)
            ax22.text(2.15, 7.25, "Physics-Strong", color=CFG.VIS_ACCENT["blue"], **_main_lbl_kw22,
                      bbox=dict(boxstyle="round,pad=0.40", fc=CFG.VIS_RAMP["blue"][0], ec=CFG.VIS_ACCENT["blue"],
                                alpha=0.95, linewidth=1.8))
            ax22.text(7.45, 7.25, "AI-Strong", color=CFG.VIS_ACCENT["vermillion"], **_main_lbl_kw22,
                      bbox=dict(boxstyle="round,pad=0.40", fc=CFG.VIS_RAMP["orange"][0], ec=CFG.CONF_BAND_COLOURS["below"],
                                alpha=0.95, linewidth=1.8))
            ax22.text(5.00, 1.55, "Tier-Strong", color=CFG.VIS_ACCENT["green"], **_main_lbl_kw22,
                      bbox=dict(boxstyle="round,pad=0.40", fc=CFG.VIS_BAND["high_fill"], ec=CFG.CONF_BAND_COLOURS["high"],
                                alpha=0.95, linewidth=1.8))

            _region_data22 = [
                (1.85, 5.80, _only_A,    CFG.VIS_ACCENT["blue"], "Physics only"),
                (8.15, 5.80, _only_B,    CFG.VIS_ACCENT["vermillion"], "AI only"),
                (5.00, 2.20, _only_C,    CFG.VIS_ACCENT["green"], "Tier only"),
                (5.00, 6.90, _AB_not_C,  CFG.VIS_ACCENT_DEEP["blue_alt"], "Physics ∩ AI"),
                (2.80, 4.40, _AC_not_B,  CFG.VIS_ACCENT_DEEP["green"], "Physics ∩ Tier"),
                (7.20, 4.40, _BC_not_A,  CFG.VIS_ACCENT_DEEP["orange_alt"], "AI ∩ Tier"),
                (5.00, 5.30, _ABC,       CFG.VIS_INK["outline"],    "All three"),
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
                          bbox=dict(boxstyle="square,pad=0.32", fc=CFG.VIS_INK["canvas"], ec=_rc22,
                                    alpha=0.92, linewidth=0.9), zorder=5)

            _ctrl_items22 = [
                (_ctrl_3R3U_idx22,  "3R3U (control)",  CFG.VIS_ACCENT_DEEP["red_deep"]),
                (_ctrl_deha4_idx22, "DeHa4 (control)", CFG.VIS_ACCENT_DEEP["blue"]),
            ]
            _ctrl_rpos22 = {
                (True,  True,  True):  (5.00, 5.10),
                (True,  True,  False): (5.00, 6.90),
                (True,  False, True):  (2.80, 4.40),
                (False, True,  True):  (7.20, 4.40),
                (True,  False, False): (1.85, 5.80),
                (False, True,  False): (8.15, 5.80),
                (False, False, True):  (5.00, 2.20),
                # A control in NONE of the sets sits outside every circle, with margin from the
                # corner so the placement reads as deliberate.
                (False, False, False): (1.70, 1.45),
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
                _Patch22(facecolor=CFG.VIS_ACCENT["blue"], alpha=0.45, edgecolor=CFG.VIS_ACCENT["blue"],
                         linewidth=1.2, label=f"Physics-Strong  (top-tercile {_phys_col22})"),
                _Patch22(facecolor=CFG.VIS_ACCENT["vermillion"], alpha=0.45, edgecolor=CFG.VIS_ACCENT["vermillion"],
                         linewidth=1.2, label="AI-Strong  (top-tercile confidence)"),
                _Patch22(facecolor=CFG.VIS_ACCENT["green"], alpha=0.45, edgecolor=CFG.VIS_ACCENT["green"],
                         linewidth=1.2, label=f"Tier-Strong  ({CFG.TIER_TOP}, {CFG.TIER_ORDER[1]}, {CFG.TIER_ORDER[2]})"),
            ]
            if _n_PA22 > 0:
                _leg22.append(_L22([0], [0], color=CFG.VIS_ACCENT["magenta"], linewidth=2.2,
                                   linestyle="--",
                                   label=f"{CFG.TIER_TOP} overlay  (n={_n_PA22:,})"))
            _leg22.append(_L22([0], [0], color=CFG.VIS_INK["soft"], linewidth=1.8,
                               linestyle=":", marker="o", markersize=12,
                               markerfacecolor="none", markeredgecolor=CFG.VIS_INK["soft"],
                               markeredgewidth=1.5,
                               label=f"All-three intersection  ({_ABC:,})"))
            _leg22.extend(_ctrl_legend_h22)
            # Legend sits close under the diagram
            ax22.legend(handles=_leg22,
                        loc="upper center", bbox_to_anchor=(0.5, 0.02),
                        ncol=3,   fancybox=True)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_23_Category_Overlap_Euler.png",
                        dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close()
        else:
            reporter.log(f"  ! Skipped: {_fig_path('23')} — requires a catalytic-geometry column (soft_catalytic_score/mechanistic_score/SN2_Attack_Angle), a confidence column, and degrader_tier")
    except Exception as e:
        reporter.log(f"  ! Skipped: {_fig_path('23')} — {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    _fig23_multitarget(df, out_dir, reporter)
    _fig24_sankey(df, out_dir, reporter)
    _fig25_pfas_size(df, out_dir, reporter)


# ── Multi-model per-Boltz-model geometry variance (folder 08, figs 13–15) ──────────────────────
# Ported from the standalone TEST/generate_all_figures.py into 03's conventions (CFG palette / DPI /
# fonts, apply_figure_style, consistent tier separators). Shows how the 5 Boltz diffusion models
# spread the SN2 distance/angle within each tier, and the same for the high-confidence subsets.
def _mm_load_variance_df(fig_root: Path, df: pd.DataFrame):
    """Per-model variance CSV (one row per model per complex) joined to the pipeline tier from df.
    Returns the merged, tier-ordered frame, or None when the variance CSV is absent."""
    vpath = _aux_dir(fig_root) / "07_Boltz2_MultiModel_QC_Variance.csv"
    if not vpath.exists():
        return None
    try:
        var = pd.read_csv(vpath)
    except Exception:
        return None
    var = var[pd.to_numeric(var.get("sn2_distance_A"), errors="coerce") < 100]      # drop ~999 Å sentinel
    _jc = _xn__col(df, "job_name") or "job_name"
    if _jc not in df.columns or CFG.COL_TIER not in df.columns or "complex_id" not in var.columns:
        return None
    tdf = (df[[_jc, CFG.COL_TIER]].rename(columns={_jc: "complex_id"}).drop_duplicates("complex_id"))
    m = var.merge(tdf, on="complex_id", how="left")
    m = m[m[CFG.COL_TIER].isin(TIER_ORDER_LOGIC)].copy()
    if m.empty:
        return None
    m[CFG.COL_TIER] = pd.Categorical(m[CFG.COL_TIER], categories=TIER_ORDER_LOGIC, ordered=True)
    return m


def _mm_nac_lines(ax_d, ax_a):
    """Strict-NAC reference lines shared by the multi-model panels (CFG-sourced)."""
    _dl = float(CFG.NAC_DIST_STRICT)
    _al = float(CFG.NAC_ANGLE_STRICT)
    ax_d.axhline(_dl, color=CFG.VIS_ACCENT["blue"], ls="--", lw=1.2, label=f"strict NAC {_dl:g} Å")
    ax_a.axhline(_al, color=CFG.VIS_ACCENT["magenta"], ls="--", lw=1.2, label=f"strict NAC {_al:g}°")


def _mm_variance_by_model(df, out_dir, reporter):
    """Fig 13 — per-Boltz-model SN2 distance/angle spread by tier (are the 5 models consistent?)."""
    mv = _mm_load_variance_df(out_dir.parent, df)
    if mv is None:
        reporter.log("  ! Skipped: 08_Diagnostic_and_MultiModel_Trends/13_MultiModel_Geometry_Variance_by_Model.png — variance CSV unavailable.")
        return
    _tiers = [t for t in TIER_ORDER_LOGIC if t in set(mv[CFG.COL_TIER].dropna())]
    _models = sorted(str(m) for m in mv["model_name"].dropna().unique())
    _acc = CFG.VIS_ACCENT
    _mpal = [_acc[k] for k in ("green", "sky", "blue", "magenta", "amber", "vermillion") if k in _acc]
    _mpal = (_mpal * (len(_models) // max(1, len(_mpal)) + 1))[:max(1, len(_models))]
    fig, (ax_d, ax_a) = plt.subplots(2, 1, figsize=(13, 8.5), sharex=True)
    sns.boxplot(data=mv, x=CFG.COL_TIER, y="sn2_distance_A", hue="model_name", order=_tiers,
                hue_order=_models, palette=_mpal, ax=ax_d, fliersize=0, linewidth=0.8)
    ax_d.set_ylabel("Nucleophile distance (Å)"); ax_d.set_xlabel("")
    sns.boxplot(data=mv, x=CFG.COL_TIER, y="sn2_angle_deg", hue="model_name", order=_tiers,
                hue_order=_models, palette=_mpal, ax=ax_a, fliersize=0, linewidth=0.8)
    ax_a.set_ylabel("SN2 attack angle (°)"); ax_a.set_xlabel("Degrader tier")
    _mm_nac_lines(ax_d, ax_a)
    for _ax in (ax_d, ax_a):
        _tier_seps(_ax, len(_tiers)); _ax.spines["top"].set_visible(False); _ax.spines["right"].set_visible(False); _ax.grid(axis="y", alpha=0.3)
        if _ax.get_legend():
            _ax.get_legend().remove()
    ax_a.tick_params(axis="x", rotation=30)
    # Model legend INSIDE the top panel (its upper band is empty — data sits at 2–4.5 Å); no title.
    # The strict-NAC reference lines are drawn by _mm_nac_lines but sit outside seaborn's hue legend,
    # so add them explicitly (CFG values, matching the dashed-line colours) — otherwise the lines show
    # with no key.
    _dl = float(CFG.NAC_DIST_STRICT); _al = float(CFG.NAC_ANGLE_STRICT)
    _h = [Line2D([0], [0], marker="s", ls="", color=_mpal[i % len(_mpal)]) for i in range(len(_models))]
    _h += [Line2D([0], [0], color=CFG.VIS_ACCENT["blue"], ls="--", lw=1.2),
           Line2D([0], [0], color=CFG.VIS_ACCENT["magenta"], ls="--", lw=1.2)]
    _labs = list(_models) + [f"strict NAC {_dl:g} Å", f"strict NAC {_al:g}°"]
    ax_d.legend(_h, _labs, loc="upper left", ncol=len(_labs),
                fontsize=CFG.VIS_FONT_LEGEND, columnspacing=1.0)   # framealpha/handletextpad ← CFG SSOT
    fig.tight_layout()
    fig.savefig(out_dir / "13_MultiModel_Geometry_Variance_by_Model.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def _mm_variance_cut(df, out_dir, reporter, *, filt_col, thr, num, fname):
    """Shared box+strip of per-model geometry for a high-quality subset (confidence / ipTM ≥ thr)."""
    mv = _mm_load_variance_df(out_dir.parent, df)
    if mv is None:
        reporter.log(f"  ! Skipped: 08_Diagnostic_and_MultiModel_Trends/{num}_{fname}.png — variance CSV unavailable.")
        return
    if filt_col in mv.columns:
        mv = mv[pd.to_numeric(mv[filt_col], errors="coerce") >= thr]
    if mv.empty:
        reporter.log(f"  ! Skipped: {num}_{fname}.png — no rows with {filt_col} ≥ {thr}.")
        return
    _tiers = [t for t in TIER_ORDER_LOGIC if t in set(mv[CFG.COL_TIER].dropna())]
    _tpal = [TIER_PALETTE.get(t, CFG.VIS_INK["faint"]) for t in _tiers]   # each tier its own CFG colour
    fig, (ax_d, ax_a) = plt.subplots(2, 1, figsize=(11, 9.5), sharex=True)
    for _ax, _y, _ylab in ((ax_d, "sn2_distance_A", "Nucleophile distance (Å)"),
                           (ax_a, "sn2_angle_deg", "SN2 attack angle (°)")):
        sns.boxplot(data=mv, x=CFG.COL_TIER, y=_y, order=_tiers, hue=CFG.COL_TIER, palette=_tpal,
                    dodge=False, legend=False, fliersize=0, ax=_ax, linewidth=0.8)
        sns.stripplot(data=mv, x=CFG.COL_TIER, y=_y, order=_tiers, color=CFG.VIS_INK["dark"],
                      alpha=0.12, jitter=0.28, size=1.2, ax=_ax)
        _ax.set_ylabel(_ylab)
        _tier_seps(_ax, len(_tiers)); _ax.spines["top"].set_visible(False); _ax.spines["right"].set_visible(False); _ax.grid(axis="y", alpha=0.3)
    ax_d.set_xlabel(""); ax_a.set_xlabel("Degrader tier"); ax_a.tick_params(axis="x", rotation=30)
    for _lab, _t in zip(ax_a.get_xticklabels(), _tiers):   # tick label inherits its tier's box colour
        _lab.set_color(TIER_PALETTE.get(_t, CFG.VIS_INK["dark"]))
    _mm_nac_lines(ax_d, ax_a)
    ax_d.legend(loc="upper right", fontsize=CFG.VIS_FONT_LEGEND)
    ax_a.legend(loc="lower right", fontsize=CFG.VIS_FONT_LEGEND)
    fig.tight_layout()
    fig.savefig(out_dir / f"{num}_{fname}.png", dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)


def _mm_variance_conf90(df, out_dir, reporter):
    """Fig 14 — multi-model geometry variance, confidence ≥ 0.90 subset."""
    _mm_variance_cut(df, out_dir, reporter, filt_col="confidence_score", thr=0.90,
                     num="14", fname="MultiModel_Geometry_Variance_conf90")


def _mm_variance_iptm90(df, out_dir, reporter):
    """Fig 15 — multi-model geometry variance, ipTM ≥ 0.90 subset."""
    _mm_variance_cut(df, out_dir, reporter, filt_col="iptm", thr=0.90,
                     num="15", fname="MultiModel_Geometry_Variance_iptm90")


def _generate_comprehensive_figures_impl(df: pd.DataFrame, features: list[str], out_dir: Path, reporter: ReportManager):
    """Orchestrate the publication figure suite (folders 03-07): build the shared prologue
    (tier order, companion thumbnails, interaction-column maps) once, then render each folder.
    Runs inside the _redirect_savefig context opened by generate_comprehensive_figures."""
    global _xn__PROD_DIR
    _xn__PROD_DIR = out_dir.parent / "1_Boltz2_Production"   # per-folder analysis panels read the model CIFs from here
    _utils_mod.apply_figure_style(CFG)
    existing_tiers = [t for t in TIER_ORDER_LOGIC if t in df[CFG.COL_TIER].unique()]
    # ── Companion thumbnail setup (shared by Figs 04b, 12b, 13b, 17b) ─────────
    # Source the PyMOL structure thumbnails from the MD-ready cohort (MD_Selected,
    # SECTION 18) so only the ~9-10 complexes actually taken to MD are rendered —
    # not all Tier_1A. Falls back to the top tier if the flag is absent (older CSV).
    _md_col = getattr(CFG, "MD_SELECTED_COL", "MD_Selected")
    if _md_col in df.columns:
        _sel = df[_md_col].astype(str).str.lower().isin(["true", "1", "1.0"])
        _pa = df[_sel].reset_index(drop=True)
    elif CFG.COL_TIER in df.columns:
        _pa = df[df[CFG.COL_TIER] == CFG.TIER_TOP].reset_index(drop=True)
    else:
        _pa = pd.DataFrame()
    # Rank-sort so the capped thumbnail set (CFG.VIS_MAX_THUMBNAILS) takes the best
    # entries; Scientific_Rank ascending = best first, else Boltz confidence.
    if not _pa.empty:
        if "Scientific_Rank" in _pa.columns:
            _pa = _pa.sort_values("Scientific_Rank", ascending=True).reset_index(drop=True)
        elif CFG.COL_CONF in _pa.columns:
            _pa = _pa.sort_values(CFG.COL_CONF, ascending=False).reset_index(drop=True)
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
    def _panel(fn, folder, name):
        """Render one folded per-folder analysis panel into its folder and buffer the save; a failure is
        logged and the rest still run. The save is flushed (in numeric order) with the folder's block."""
        try:
            (out_dir / folder).mkdir(parents=True, exist_ok=True)   # folder created at its turn
            fn(df, out_dir / folder, reporter)
            if (out_dir / folder / f"{name}.png").exists():
                reporter.log(f"  ✔ Saved: {folder}/{name}.png")
        except Exception as _e:                                  # noqa: BLE001
            reporter.log(f"  ! Skipped: {folder}/{name}.png — {type(_e).__name__}: {_e}")
            plt.close("all")

    _fig_folder03_dataset(df, features, out_dir, reporter, existing_tiers)
    _panel(_xo__fig_04A_evolutionary_phylogeny, "03_Dataset_and_Alignment_Overview", "05_Evolutionary_Phylogeny")
    _panel(_jf_treemap, "03_Dataset_and_Alignment_Overview", "06_Treemap_Tier_Ligand_Composition")

    _fig_folder04_ai_confidence(df, features, out_dir, reporter, existing_tiers, _pa, _imgs, _tt_has_imgs)

    _fig_folder05_catalytic(df, features, out_dir, reporter, existing_tiers, _pa, _imgs, _tt_has_imgs)
    _panel(_xn__fig_01C_geometry_and_uncertainty, "05_Catalytic_Geometry_and_Mechanism", "13_Geometry_and_Uncertainty")
    _panel(_xo__fig_05b_mechanistic_size_modified, "05_Catalytic_Geometry_and_Mechanism", "14_Mechanistic_Breakdown_by_Tier")
    _panel(_jf_network, "05_Catalytic_Geometry_and_Mechanism", "15_Metric_CoVariation_Network")

    _fig_folder06_ligand(df, features, out_dir, reporter, existing_tiers, _pa, _imgs, _tt_has_imgs, int_label_map, present_int_cols)
    _panel(_xo__fig_02A_binding_affinity_metrics, "06_Ligand_Interactions_and_Chemical_Space", "08_Binding_Affinity_Metrics")

    _fig_folder07_pfas(df, features, out_dir, reporter, existing_tiers)
    _panel(_xo__fig_05c_size_by_tier_modified, "07_PFAS_Scope_and_Synthesis", "14_Chain_Length_by_Tier")
    _panel(_xn_figure_06a, "07_PFAS_Scope_and_Synthesis", "15_Tier1A_Cross_Ligand_Heatmap")
    _panel(_jf_swimmer, "07_PFAS_Scope_and_Synthesis", "16_Swimmer_Top_Per_Tier")
    _panel(_jf_circos, "07_PFAS_Scope_and_Synthesis", "17_Circos_Ligand_Tier_Assignment")

    # Step 8/8 — 08_Diagnostic_and_MultiModel_Trends (its own step, last folder on disk)
    _diag_dir = out_dir / "08_Diagnostic_and_MultiModel_Trends"
    _diag_dir.mkdir(parents=True, exist_ok=True)
    generate_additional_figures(df, _diag_dir, reporter)
    _panel(_xn__fig_05a_pillar_divergence_modified, "08_Diagnostic_and_MultiModel_Trends", "09_Pillar_Divergence_by_Tier")
    _panel(_jf_volcano, "08_Diagnostic_and_MultiModel_Trends", "10_Volcano_Metric_Significance")
    _panel(_jf_manhattan, "08_Diagnostic_and_MultiModel_Trends", "11_Manhattan_Candidate_Significance")
    _panel(_jf_importance, "08_Diagnostic_and_MultiModel_Trends", "12_Feature_Importance_Tier_Streams")
    _panel(_mm_variance_by_model, "08_Diagnostic_and_MultiModel_Trends", "13_MultiModel_Geometry_Variance_by_Model")
    _panel(_mm_variance_conf90, "08_Diagnostic_and_MultiModel_Trends", "14_MultiModel_Geometry_Variance_conf90")
    _panel(_mm_variance_iptm90, "08_Diagnostic_and_MultiModel_Trends", "15_MultiModel_Geometry_Variance_iptm90")





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
        tier_colors = {t: TIER_PALETTE.get(t, CFG.VIS_INK["faint"]) for t in tier_order}
        d["tier_rank"] = d["tier"].map(CFG.TIER_RANK).fillna(0)
        d["prot"] = d[CFG.COL_PROT] if CFG.COL_PROT in d.columns else d.get(CFG.COL_PROT, d.get("protein", ""))
        d["lig"] = d[CFG.COL_LIG] if CFG.COL_LIG in d.columns else d.get(CFG.COL_LIG, d.get("ligand", ""))
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
            )
        ax.set_ylabel("Protein  (top 25, ranked by quality-weighted degradation breadth)",
                      )
        ax.tick_params(axis="y", labelsize=9)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_visible(False)
        ax.xaxis.grid(True, color=CFG.VIS_INK["hairline"], linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        for i, prot in enumerate(bar_data.index):
            total  = int(bar_data.loc[prot].sum())
            pa_cnt = int(bar_data.loc[prot, CFG.TIER_TOP]) \
                     if CFG.TIER_TOP in bar_data.columns else 0
            _lbl = f"{total}" if pa_cnt == 0 else f"{total}  ★×{pa_cnt}"
            ax.text(total + 0.15, i, _lbl, va="center", fontsize=8.5,
                    color=tier_colors[CFG.TIER_TOP] if pa_cnt > 0 else CFG.VIS_INK["soft"],
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
                      ncol=len(legend_h23), frameon=False,
                      bbox_to_anchor=(0.5, 0.045),   # lift up toward the panel, clear of the x-axis label
                      title=_title23, )
        plt.tight_layout(rect=[0, 0.07, 1, 1])
        plt.savefig(out_dir / "Figure_24a_Top25_Multitarget_Proteins.png",
                    dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close(fig23a)

        _fig23b_toptier_breakdown(best_pairs, tier_order, tier_colors, out_dir, reporter)

    except Exception as e:
        reporter.log(f"  ! Skipped: {_fig_path('24')} — {e}")
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
            reporter.log(f"  ! Skipped: {_fig_path('24b')} — no proteins reached {top_tier}")
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
                col = lig_colour.get(r["lig"], CFG.VIS_INK["faint"])
                ax.barh(yi, 1.0, left=left, height=_BH, color=col, edgecolor="white",
                        linewidth=0.8, zorder=3)
                ax.text(left + 0.5, yi, _tshort(r["tier"]), ha="center", va="center",
                        fontsize=7.0, fontweight="bold", color=_label_col(col), zorder=5)
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
                _tc = tier_colors.get(_tier_seq[_s], CFG.VIS_INK["mid"])
                ax.add_patch(_mp.Rectangle((_s, yi - _CAT_H / 2), (_e - _s + 1), _CAT_H,
                             fill=False, edgecolor=_tc, linewidth=2.3, zorder=4))
                _s = _e + 1

            ntop = int(order.loc[prot, "_ntop"])
            ax.text(left + 0.25, yi, f"{int(left)}  ★×{ntop}", va="center",
                    fontsize=8, fontweight="bold",
                    color=tier_colors.get(top_tier, CFG.VIS_INK["soft"]))
        ax.set_xlim(-0.35, float(order["_ntot"].max()) + 3.0)

        ax.set_yticks(range(len(prots)))
        ax.set_yticklabels(prots, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Number of PFAS degraded  (each block = one PFAS; block label = its tier)",
                      )
        ax.set_ylabel(f"Top-tier proteins  (reach {_tshort(top_tier)} on ≥1 PFAS;  ★×N = {_tshort(top_tier)} count)",
                      )
        for _s in ("top", "right", "left"):
            ax.spines[_s].set_visible(False)
        ax.xaxis.grid(True, color=CFG.VIS_INK["hairline"], linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)

        shown = [l for l in lig_order
                 if l in set(deg[deg["prot"].isin(prots)]["lig"])]
        leg = [_mp.Patch(color=lig_colour[l], label=l) for l in shown]
        ax.legend(handles=leg, loc="center left", bbox_to_anchor=(1.01, 0.5),
                   frameon=False, title="PFAS ligand", )
        plt.tight_layout()
        out = out_dir / "Figure_24b_TopTier_Protein_PFAS_Breakdown.png"
        plt.savefig(out, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
        plt.close(fig)
    except Exception as e:
        reporter.log(f"  ! Skipped: {_fig_path('24b')} — {e}")
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
        # Duplicate column labels make every subsequent assignment raise "cannot insert X, already
        # exists" — pandas cannot address a name that resolves to more than one column. Collapse them
        # before anything is written.
        d = d.loc[:, ~d.columns.duplicated()].copy()
        for _stale24 in ["g_alpha", "g_nuc", "g_triad", "g_ang", "g_mech",
                         "tier_cat", "md_cat"]:
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

        # The score the tier gate actually keys on: raw geometry MINUS the graded chemistry and
        # containment penalties. Binning the raw mechanistic_score would make the column a
        # non-sequitur — a complex could sit in the top band and still land in a low tier, with the
        # penalty that demoted it invisible to the reader.
        _mech_col_24 = ("mechanistic_score_effective"
                        if "mechanistic_score_effective" in d.columns else CFG.COL_MECH_S)

        """
        ── DECISION FUNNEL ────────────────────────────────────────────────────────────────────────
        The diagram follows the pipeline's actual decision, so the stream NARROWS: every complex
        enters at the left, each gate eliminates those that fail it, and what survives to the right
        is the MD-selected cohort. The eliminated complexes leave the diagram at the gate that
        rejected them — they are not re-binned into the next column, which would draw a
        classification rather than a screen and read as though nothing was ever rejected. The one
        number the figure exists to explain is how the screen gets from 58,056 complexes to a handful
        of MD systems, and that number is only visible if the stream is allowed to narrow.

        A complex eliminated at any gate joins the REJECTED band and stays there: it is never
        re-admitted, so the band only grows and the surviving stream only shrinks. The stage at which
        the band widens is the stage that did the work. Gates are applied in the ladder's own order
        and every threshold comes from CFG.
        """
        _sc24 = _num24("scissile_is_alpha", 0.0).fillna(0.0)
        _dn24 = _num24("Dist_Nucleophile", np.inf)
        _nbv24, _bav24 = _num24("dist_nuc_base_internal", np.inf), _num24("dist_base_acid_internal", np.inf)
        _angv24 = _num24("sn2_attack_angle_effective",
                         np.nan) if "sn2_attack_angle_effective" in d.columns else _num24(CFG.COL_SN2, np.nan)
        _mv24 = _num24(_mech_col_24, np.nan)

        _LOOSE = CFG.TIER_ORDER[3]          # Tier_2B — the loosest degrader rung
        _gates_24 = [
            ("g_alpha", _sc24 >= 1.0,
             "α-Carbon Attack", "(reactive centre)", CFG.VIS_ACCENT_DEEP["green_alt"], "No α-attack"),
            ("g_nuc", _dn24 <= CFG.TIER_NUC_DIST[_LOOSE],
             "Nucleophile Reach", f"(≤ {CFG.TIER_NUC_DIST[_LOOSE]:g} Å)", CFG.VIS_ACCENT_DEEP["moss"], "Out of reach"),
            ("g_triad", (_nbv24 <= CFG.TIER_NB_MAX[_LOOSE]) & (_bav24 <= CFG.TIER_BA_MAX[_LOOSE]),
             "Catalytic Relay", "(Nuc–Base · Base–Acid)", CFG.VIS_ACCENT_DEEP["pine"], "Relay broken"),
            ("g_ang", _angv24 >= CFG.TIER_ANGLE_MIN[_LOOSE],
             "SN2 Attack Angle", f"(effective ≥ {CFG.TIER_ANGLE_MIN[_LOOSE]:g}°)", CFG.VIS_RAMP["blue"][4], "Wrong trajectory"),
            ("g_mech", _mv24 >= CFG.TIER_MECH_MIN[CFG.TIER_ORDER[2]],
             "Mechanistic Score", f"(≥ {CFG.TIER_MECH_MIN[CFG.TIER_ORDER[2]]:.2f})", CFG.VIS_RAMP["purple"][4], "Machinery incomplete"),
        ]

        """
        A rejected complex EXITS at the gate that killed it. It is drawn once, as a dead-end node in
        that gate's column, and then it is gone — no node and no ribbon in any later column, because
        the ribbon loop skips a category that has no node. Carrying a rejected band all the way to
        the right would say the opposite of what the pipeline does: it would show 58,056 complexes
        arriving at the end, which is precisely the impression this figure has to destroy.

        Each gate therefore has exactly two nodes — the survivors, who continue, and the rejects,
        who stop — and the surviving stream is visibly thinner at every step.
        """
        """
        EVERY complex is tiered — the ladder assigns a tier to all 58,056, including the ones that
        fail a catalytic gate (they land in Tier_3/4/5_Decoy). So the gates do not delete anyone: a
        complex that fails one drops into a FADED 'already rejected' lane, carries on through the
        remaining gate columns in that lane, and still arrives at its real tier. The dead-end node at
        each gate names WHICH gate rejected it and how many.

        The tier column therefore holds all 58,056, distributed across the tiers — which is what the
        pipeline actually produces. Only the MD-ready cohort continues past it; everything else stops
        at its tier, so the final column is the handful of complexes taken to simulation.

        Colour carries the survival story: the passing stream is drawn dark and saturated, the
        rejected lanes light and washed out, so the eye follows the candidates that are still alive.
        """
        _OUT24 = "__out__"                   # no node exists for this label → the complex stops here
        _REJLANE24 = "Rejected earlier"
        _alive24 = pd.Series(True, index=d.index)
        _stage_cols_24 = []
        for _key24, _pass24, _ttl24, _sub24, _bg24, _faillbl24 in _gates_24:
            _pass_now = _alive24 & _pass24.fillna(False)
            _fail_now = _alive24 & ~_pass_now          # failed HERE, not earlier
            _n_fail24 = int(_fail_now.sum())
            # The node drawer already prints the count; adding it here prints it twice.
            _rej_lbl24 = f"✗ {_faillbl24}"
            d[_key24] = np.where(_pass_now, "Pass",
                                 np.where(_fail_now, _rej_lbl24, _REJLANE24))
            _alive24 = _pass_now
            _ord24 = ["Pass"] + ([_rej_lbl24] if _n_fail24 else []) \
                + ([_REJLANE24] if (d[_key24] == _REJLANE24).any() else [])
            _stage_cols_24.append((_key24, _ord24, _ttl24, _sub24, _bg24, _rej_lbl24))

        # Every complex carries its real tier — the ladder assigns one to all of them.
        d["tier_cat"] = d[CFG.COL_TIER].astype(str)
        _tier_ord_24 = [t for t in CFG.TIER_ORDER if (d["tier_cat"] == t).any()]

        # Terminal column: only the MD-ready cohort continues. Everything else stops at its tier.
        _mdmask24 = _md_ready_df(d).index
        _MDSEL24 = "MD-selected"
        d["md_cat"] = np.where(d.index.isin(_mdmask24), _MDSEL24, _OUT24)
        _md_ord_24 = [_MDSEL24] if (d["md_cat"] == _MDSEL24).any() else []

        _col_seq_24  = ["all"] + [c for c, *_ in _stage_cols_24] + ["tier_cat", "md_cat"]
        _ord_seq_24  = [None]  + [o for _, o, *_ in _stage_cols_24] + [_tier_ord_24, _md_ord_24]
        _ttl_seq_24  = ["All Complexes"] + [t for _, _, t, *_ in _stage_cols_24] + ["Degrader Tier", "MD Cohort"]
        _sub_seq_24  = ["(starting pool)"] + [s for _, _, _, s, *_ in _stage_cols_24] \
            + ["(all complexes tiered)", "(taken to MD)"]
        _bg_seq_24   = [CFG.VIS_INK["ink_navy"]] + [b for _, _, _, _, b, _ in _stage_cols_24] + [CFG.VIS_ACCENT_DEEP["purple_dark"], CFG.VIS_ACCENT_DEEP["gold_dark"]]
        _n_cols_24   = len(_col_seq_24)

        # ── Colours ────────────────────────────────────────────────────────────
        # One survivor colour per gate (best-of-CFG's green end) and one rejected colour throughout,
        # so the eye follows a single narrowing stream against a single growing dead-end band.
        _PASS_CLR_24 = CFG.SANKEY_GRAD5[0]     # survivors: saturated
        _REJ_CLR_24  = CFG.SANKEY_GRAD5[-1]    # rejected AT this gate: red dead-end
        _LANE_CLR_24 = CFG.VIS_INK["steel"]               # already-rejected lane: washed out, recedes
        _colmap_24 = {}
        for _key24, _ords24, _t24, _s24, _b24, _rej24 in _stage_cols_24:
            _colmap_24[_key24] = {"Pass": _PASS_CLR_24, _rej24: _REJ_CLR_24,
                                  _REJLANE24: _LANE_CLR_24}
        _colmap_24["md_cat"] = {_MDSEL24: CFG.VIS_ACCENT_DEEP["gold_hot"]}
        _tier_clr_24  = {t: TIER_PALETTE.get(t, CFG.VIS_INK["faint"])
                         for t in [CFG.TIER_TOP,CFG.TIER_ORDER[1],CFG.TIER_ORDER[2],CFG.TIER_ORDER[3],CFG.TIER_ORDER[4],CFG.TIER_POOR,CFG.TIER_DECOY,"Other"]}
        _tier_alp_24  = {
            CFG.TIER_TOP: 0.72, CFG.TIER_ORDER[1]: 0.68, CFG.TIER_ORDER[2]: 0.62,
            CFG.TIER_ORDER[3]: 0.57, CFG.TIER_ORDER[4]: 0.48, CFG.TIER_POOR: 0.42, CFG.TIER_DECOY: 0.36, "Other": 0.28,
        }

        def _clr_24(cat_col, lbl):
            if str(lbl).startswith("✗"):          # a dead-end node: the complexes rejected at this gate
                return _REJ_CLR_24
            if cat_col == "tier_cat":
                return TIER_PALETTE.get(lbl, _tier_clr_24.get(lbl, CFG.VIS_INK["palest"]))
            return _colmap_24.get(cat_col, {}).get(lbl, CFG.VIS_INK["palest"])

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
            """
            Node heights are scaled against the WHOLE corpus, not against the column's own sum.
            Every intermediate column holds all 58,056 complexes, so the two agree there — but the
            final column holds only the MD-selected handful, and normalising it to itself would blow
            4 complexes up to the full height of the figure, drawn exactly as large as the 22,184 of
            Tier_3 beside it. The point of the last column is how FEW survive; it must be a sliver.
            """
            total_c = max(_total_j_24, 1)
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
            facecolor=CFG.VIS_INK["ink_navy"], edgecolor="white", linewidth=1.0, alpha=0.96, zorder=5))
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
                      ha="center", va="center", fontsize=12, color=CFG.VIS_TINT["blue"], zorder=8)

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
                          arrowprops=dict(arrowstyle="->", color=CFG.VIS_INK["ghost"], lw=1.4,
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
                     TIER_PALETTE.get(tc, _tier_clr_24.get(tc, CFG.VIS_INK["faint"])),
                     _tier_alp_24.get(tc, 0.40), zorder=2)

        # ── Ribbons col i → col i+1 (i = 1 .. n_cols-2) ─────────────────────
        for _ri in range(1, _n_cols_24 - 1):
            _sc24 = _col_seq_24[_ri]
            _tc24 = _col_seq_24[_ri + 1]
            _sp24 = _cpos_24.get(_sc24, {})
            _tp24 = _cpos_24.get(_tc24, {})
            _sn24 = _ccnt_24.get(_sc24, {})
            _tn24 = _ccnt_24.get(_tc24, {})
            # De-duplicated: on the last hop the source column IS tier_cat, so appending it again for
            # the ribbon ordering key names it twice and reset_index cannot insert it a second time.
            _grp_cols_24 = list(dict.fromkeys([_sc24, _tc24, "tier_cat"]))
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
                """
                The ribbon's weight follows the complex's fate. A flow that is still ALIVE — leaving a
                Pass node and not entering a dead end — is drawn in its tier colour at full strength.
                A flow that has been rejected, whether it is falling into a dead end now or drifting
                along in the already-rejected lane, is drawn pale and translucent so it recedes. The
                reader then follows the surviving candidates without having to read a single label.
                """
                _dead_src = (s0 == _REJLANE24) or str(s0).startswith("✗")
                _dead_tgt = (t0 == _REJLANE24) or str(t0).startswith("✗")
                if _dead_src or _dead_tgt:
                    _rc24, _ra24 = _LANE_CLR_24, 0.30
                else:
                    _rc24 = TIER_PALETTE.get(tc, _tier_clr_24.get(tc, CFG.VIS_INK["faint"]))
                    _ra24 = min(0.92, _tier_alp_24.get(tc, 0.40) + 0.28)
                _brib_24(_xs_24[_ri], _bw_24, _xs_24[_ri + 1],
                         y0a, y0a + _rhs, y0b, y0b + _rht,
                         _rc24, _ra24, zorder=2 + _ri)

        """
        No legend or title — the Final Tier column already carries the tier names and colours.
        """
        plt.tight_layout(pad=0.3)
        plt.savefig(out_dir / "Figure_25_Sankey_Workflow.png",
                    dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight", facecolor="white")
        plt.close(fig24)

    except Exception as e:
        reporter.log(f"  ! Skipped: {_fig_path('25')} — {e}")
        plt.close("all")   # release the figure left open by the failed savefig


def _fig25_pfas_size(df: pd.DataFrame, out_dir: Path, reporter) -> None:
    """Figure 26a–25e — PFAS chain-length selectivity saved as five separate PNGs."""
    try:
        _req25 = ["total_fluorine_count", CFG.COL_SN2, CFG.COL_CONF, CFG.COL_TIER]
        _miss25 = [c for c in _req25 if c not in df.columns]
        if _miss25:
            reporter.log(f"  ! Skipped: {_fig_path('26')} — missing columns {_miss25}")
            return

        # MD_Selected is carried through: subsetting to _req25 alone drops it, and the MD stars then
        # silently never draw — the figure looks finished and is simply missing its point.
        _keep25 = _req25 + [c for c in ("MD_Selected", "Control_Ref", "is_control") if c in df.columns]
        _d25 = df[_keep25].copy()
        for _c25 in _req25[:-1]:
            _d25[_c25] = pd.to_numeric(_d25[_c25], errors="coerce")
        _d25 = _d25.dropna(subset=_req25)
        _d25 = _d25[_d25[CFG.COL_TIER].isin(TIER_ORDER_LOGIC)]
        if len(_d25) < 30:
            reporter.log(f"  ! Skipped: {_fig_path('26')} — insufficient data after filtering")
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
        _ang25  = _d25[CFG.COL_SN2]
        _conf25 = _d25[CFG.COL_CONF]
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

        # ── Panel A (25a): hexbin landscape + tier scatter + rolling median ──
        fig25a, axA = plt.subplots(figsize=(14, 8))
        axA.axhspan(CFG.SUBSTRATE_ANGLE_MIN, 185, color=CFG.VIS_ACCENT["green"], alpha=0.08, zorder=0)
        axA.axhspan(90, CFG.INHIBITOR_ANGLE_MAX, color=CFG.VIS_ACCENT["vermillion"], alpha=0.08, zorder=0)
        axA.axhline(CFG.SUBSTRATE_ANGLE_MIN, color=CFG.VIS_ACCENT["green"], lw=1.2, ls="--", alpha=0.7, zorder=2)
        axA.axhline(CFG.INHIBITOR_ANGLE_MAX, color=CFG.VIS_ACCENT["vermillion"], lw=1.2, ls="--", alpha=0.7, zorder=2)
        _hb25 = axA.hexbin(
            _d25["total_fluorine_count"], _d25[CFG.COL_SN2],
            gridsize=50, mincnt=1, cmap="YlOrRd", alpha=0.65, zorder=1, linewidths=0.2
        )
        _cb25 = fig25a.colorbar(_hb25, ax=axA, pad=0.01, aspect=30, shrink=0.85)
        _cb25.set_label("Count per hex", fontsize=9)
        _cb25.ax.tick_params(labelsize=8)
        _samp_A25 = _d25.sample(min(3000, len(_d25)), random_state=int(CFG.ANALYSIS_SEED))
        for _t25a in sorted(TIER_ORDER_LOGIC, key=lambda t: t == CFG.TIER_TOP):
            _is_pa25a = (_t25a == CFG.TIER_TOP)
            # Tier_1A is crucial and rare → plot ALL of its points (never sampled)
            _src25a = _d25 if _is_pa25a else _samp_A25
            _ts25 = _src25a[_src25a[CFG.COL_TIER] == _t25a]
            if _ts25.empty:
                continue
            axA.scatter(
                _ts25["total_fluorine_count"], _ts25[CFG.COL_SN2],
                color=TIER_PALETTE.get(_t25a, CFG.VIS_INK["faint"]),
                alpha=0.85 if _is_pa25a else 0.45,
                s=60 if _is_pa25a else 14,
                edgecolors="black" if _is_pa25a else "none",
                linewidths=0.8 if _is_pa25a else 0,
                zorder=5 if _is_pa25a else 3, label=_t25a
            )
        _roll25 = _d25[["total_fluorine_count", CFG.COL_SN2]].sort_values("total_fluorine_count")
        if len(_roll25) >= 20:
            _win25 = max(20, len(_roll25) // 40)
            _roll25["_med"] = _roll25[CFG.COL_SN2].rolling(
                window=_win25, center=True, min_periods=5).median()
            axA.plot(_roll25["total_fluorine_count"], _roll25["_med"],
                     color="black", lw=2.2, zorder=6, label="Rolling median SN2")

        # The MD-selected cohort, starred on the size/angle landscape: this figure argues that
        # attack geometry degrades with chain length, and these are the complexes that argument is
        # being acted on.
        _md25 = _md_ready_df(_d25)
        if not _md25.empty:
            axA.scatter(_md25["total_fluorine_count"], _md25[CFG.COL_SN2],
                        label=f"MD-selected (n={len(_md25)})",
                        **({**_MD_STAR_KW, "s": 340}))
        _ctrl25 = _control_star_df(_d25)
        if not _ctrl25.empty and "total_fluorine_count" in _ctrl25.columns:
            _ctrl_star(axA, pd.to_numeric(_ctrl25["total_fluorine_count"], errors="coerce").to_numpy(),
                       pd.to_numeric(_ctrl25[CFG.COL_SN2], errors="coerce").to_numpy(), size=360)
            axA.scatter([], [], marker="*", s=200, facecolor=CFG.VIS_ACCENT["control"],
                        edgecolor=CFG.VIS_ACCENT["control_edge"], linewidths=1.4, label="3R3U × FA (control)")
        axA.text(0.99, 0.92, f"SUBSTRATE ZONE  (SN2 ≥ {CFG.SUBSTRATE_ANGLE_MIN:.0f}°)",
                 transform=axA.transAxes, ha="right", va="top", fontsize=9.5,
                 color=CFG.VIS_BAND["high"], fontweight="bold",
                 bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"))
        axA.text(0.99, 0.10, f"POTENTIAL INHIBITOR ZONE  (SN2 < {CFG.INHIBITOR_ANGLE_MAX:.0f}°)",
                 transform=axA.transAxes, ha="right", va="bottom", fontsize=9.5,
                 color=CFG.VIS_ACCENT_DEEP["orange_hot"], fontweight="bold",
                 bbox=dict(facecolor="white", alpha=0.75, edgecolor="none"))
        _bins25 = list(CFG.VIS_PFAS_FCOUNT_BINS)          # [0, 8, 13, 18, 24, inf]
        for _fb25 in _bins25[1:-1]:
            axA.axvline(_fb25, color=CFG.VIS_INK["soft"], lw=0.8, ls=":", alpha=0.5, zorder=2)
        _bin_toplbls25 = list(CFG.VIS_PFAS_SIZE_LABELS_SHORT)
        axA.set_ylim(85, 186)
        _fc_max25 = int(_d25["total_fluorine_count"].max()) + 2
        # Label midpoints derived from the CFG bin edges; the open-ended final bin uses the axis max.
        _edges25 = _bins25[:-1] + [max(_fc_max25, _bins25[-2] + 4)]
        _bin_xmids25 = [0.5 * (_edges25[i] + _edges25[i + 1]) for i in range(len(_edges25) - 1)]
        axA.set_xticks(range(0, _fc_max25 + 1, 2))
        for _bx25, _bl25 in zip(_bin_xmids25, _bin_toplbls25):
            axA.text(_bx25, 184, _bl25, ha="center", va="top",
                     fontsize=_JF_FA, color=CFG.VIS_INK["soft"], style="italic")
        axA.set_xlabel("Total fluorine count  (proxy for carbon chain length)", )
        axA.set_ylabel("SN2 Attack Angle (°)", )
        # One row: the entry count is the tiers plus the rolling median plus the MD stars, so the
        # column count has to include all three or the legend wraps. markerscale shrinks the MD star
        # to the legend's own scale — the plotted star is deliberately large to be findable in a
        # 55,000-point field, which is the wrong size for a legend swatch.
        _nleg25 = len(axA.get_legend_handles_labels()[0])
        axA.legend(loc="upper left", bbox_to_anchor=(0.0, 1.05),
                   ncol=max(_nleg25, 1),  markerscale=0.55,
                    fancybox=True)
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
            _d25.groupby(["_size_bin", CFG.COL_TIER], observed=True)
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
                         fontsize=7.0, color=_bin_cols_25[_xi], fontweight="bold")

        _stack_at_25(_xB - _offB, _oc_pct25, _oc_ct25, lambda c: _oc25_colors_b.get(c, CFG.VIS_INK["faint"]), "left")
        _stack_at_25(_xB + _offB, _tier_pct25, _tier_ct25, lambda t: TIER_PALETTE.get(t, CFG.VIS_INK["faint"]), "right")
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
        axB.tick_params(axis="x", pad=16)   # room for the outcome/tier mini-tags
        axB.set_ylabel("Percentage of complexes (%)", )
        axB.set_xlabel("PFAS chain-length bin", )
        axB.spines["top"].set_visible(False)
        axB.spines["right"].set_visible(False)
        axB.yaxis.grid(True, color=CFG.VIS_INK["mist"], linewidth=0.7, zorder=0)
        axB.set_axisbelow(True)

        from matplotlib.patches import Patch as _P25b
        _oc_h25 = [_P25b(facecolor=_oc25_colors_b.get(c, CFG.VIS_INK["faint"]), edgecolor="white",
                         linewidth=0.5, label=c) for c in _oc_pct25.columns]
        _ti_h25 = [_P25b(facecolor=TIER_PALETTE.get(t, CFG.VIS_INK["faint"]), edgecolor="white",
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
             title="Left bar — mechanistic outcome", 
             fancybox=True)
        _leg_oc25.get_title().set_fontweight("bold")
        _leg_oc25._legend_box.align = "left"
        axB.add_artist(_leg_oc25)
        _leg_ti25 = axB.legend(
            handles=_ti_h25, labels=[h.get_label() for h in _ti_h25],
            loc="lower right", bbox_to_anchor=(1.0, 1.01), ncol=len(_ti_h25),
             title="Right bar — degrader tier", 
             fancybox=True)
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
            reporter.log(f"  ! Skipped: {_fig_path('26c')} — ligand SMILES not found ({_smi25})")
            raise RuntimeError("ligand properties unavailable")
        _lcol25 = next((c for c in [CFG.COL_LIG, CFG.COL_LIG, "ligand"]
                        if c in df.columns), None)
        if _lcol25 is None:
            reporter.log(f"  ! Skipped: {_fig_path('26c')} — ligand column not found in data")
            raise RuntimeError("ligand column unavailable")
        """
        Catalytic-competence columns: the continuous soft_catalytic_score (0–1) and the
        per-pose degrader flag. Confidence is the labelled control. mechanistic_score is
        deliberately not used here — it saturates near 1.0 within the degrader tiers, so
        it has no spread to show a size trend (it is a gate, not a graded competence axis).
        """
        _soft_c25 = "soft_catalytic_score" if "soft_catalytic_score" in df.columns else None
        _conf_c25 = CFG.COL_CONF if CFG.COL_CONF in df.columns else None
        _deg_c25  = "is_degrader" if "is_degrader" in df.columns else None
        if _soft_c25 is None and _deg_c25 is None:
            reporter.log(f"  ! Skipped: {_fig_path('26c')} — no catalytic-competence column (soft_catalytic_score / is_degrader)")
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
            reporter.log(f"  ! Skipped: {_fig_path('26c')} — no ligand property matches in data")
        else:
            _groups25 = sorted(_d25c["_nC"].unique())
            _present25 = set(_d25c[_lcol25].astype(str))

            def _xlabel25(nc):
                _fs = sorted({_LP25[l]["nF"] for l in _LP25
                              if _LP25[l]["nC"] == nc and l in _present25})
                return f"C{int(nc)}\n(" + ", ".join(f"{f}F" for f in _fs) + ")"

            _labels25 = [_xlabel25(g) for g in _groups25]
            _pos25    = np.arange(len(_groups25))

            def _grp25(col):
                return [_d25c[_d25c["_nC"] == g][col].dropna().values for g in _groups25]

            fig25c, axL = plt.subplots(figsize=(15, 6))
            _SOFT_C, _SOFT_F = CFG.VIS_RAMP["green"][2], CFG.VIS_RAMP["green"][1]   # soft catalytic score (green)
            _DEG_C           = CFG.VIS_RAMP["purple"][3]              # degrader fraction (purple)
            _CTRL_C          = CFG.VIS_INK["slate"]              # confidence control (grey)

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
                _mln25, = axL.plot(_pos25, _smean25, color=CFG.VIS_RAMP["orange"][2], lw=2.2, ls="--",
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
                         fontsize=8, fontweight="bold", color=CFG.VIS_INK["dark"], zorder=7)

            axL.set_xticks(_pos25)
            axL.set_xticklabels(_labels25, fontsize=8)
            axL.set_xlabel("PFAS carbon number  (fluorine counts present shown in parentheses)",
                           )
            axL.set_ylabel("Catalytic competence  (soft score · degrader fraction, 0–1)", )
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
            axL.yaxis.grid(True, color=CFG.VIS_ACCENT_DEEP["leaf"], alpha=0.30, linewidth=0.7, zorder=0)
            axL.xaxis.grid(True, color=CFG.VIS_INK["smoke"], linewidth=0.5, zorder=0)
            axL.set_axisbelow(True)
            axL.legend(handles=_handles25, loc="lower left", bbox_to_anchor=(0.0, 1.01),
                       ncol=len(_handles25))
            plt.tight_layout()
            _out25c = out_dir / "Figure_26c_PFAS_Carbon_Confidence.png"
            fig25c.savefig(_out25c, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            plt.close(fig25c)

    except Exception as e:
        reporter.log(f"  ! Skipped: {_fig_path('26')} — {e}")
        plt.close("all")   # release the figure left open by the failed savefig


# =============================================================================
# SECTION 4D: DIAGNOSTIC & MULTI-MODEL TREND FIGURES  (folder 07)
# =============================================================================

def generate_additional_figures(df: pd.DataFrame, out_dir: Path,
                                reporter: ReportManager) -> None:
    """
    Pocket-fit and multi-model diagnostic figures written to
    <Run>/3_Validation_Figures/08_Diagnostic_and_MultiModel_Trends.

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
    reporter.section("Step 8/8 — 08_Diagnostic_and_MultiModel_Trends · diagnostic + multi-model trends")
    out_dir.mkdir(parents=True, exist_ok=True)

    # CFG-sourced colours (single source of truth — config §8).
    _fit_col   = CFG.OUTCOME_COLOUR.get("Substrate", CFG.VIS_ACCENT_DEEP["teal"])
    _nofit_col = CFG.OUTCOME_COLOUR.get("Potential Inhibitor", CFG.VIS_ACCENT_DEEP["orange_deepest"])

    existing_tiers = [t for t in TIER_ORDER_LOGIC if t in df.get(CFG.COL_TIER, pd.Series()).unique()]

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
        ax.set_xticklabels(labels, rotation=40, ha="right", fontsize=_JF_FA)
        for tick, t in zip(ax.get_xticklabels(), tiers):
            tick.set_color(TIER_PALETTE.get(t, "black"))

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


    def _spear_str(x, y, panel: str = "extended figure"):
        """Spearman ρ / p / n as a one-line string for embedding in a legend title; None if <10 pts.

        The test is registered into the shared Benjamini-Hochberg family, and the p printed on the
        panel is marked uncorrected: the family is complete only after every figure has run, so no q
        exists at draw time. The corrected q_BH for this exact test is in 06_Statistical_Tests.csv,
        which is the reportable value.
        """
        _xy = pd.DataFrame({"x": pd.to_numeric(x, errors="coerce"),
                            "y": pd.to_numeric(y, errors="coerce")}).dropna()
        if len(_xy) < 10:
            return None
        _rho, _p = _sc_stats.spearmanr(_xy["x"], _xy["y"])
        if not np.isfinite(_rho):
            return None
        _register_p("Spearman correlation", panel, float(_rho), int(len(_xy)), float(_p), rho=float(_rho))
        _ptxt = "p < 0.001" if _p < 0.001 else f"p = {_p:.3f}"
        return f"Spearman ρ = {_rho:+.2f}   {_ptxt} uncorrected   n = {len(_xy):,}"

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
        sub = pd.DataFrame({"asv": asv, "lv": lv, "tier": df.get(CFG.COL_TIER)}).dropna(subset=["asv", "lv"])
        sub = sub[(sub["asv"] > 0) & (sub["lv"] > 0)]
        if len(sub) >= 5:
            fig, ax = plt.subplots(figsize=(8.4, 7.0))
            for t in existing_tiers:
                d = sub[sub["tier"] == t]
                if d.empty:
                    continue
                ax.scatter(d["asv"], d["lv"], s=14, alpha=0.30, linewidths=0,
                           color=TIER_PALETTE.get(t, CFG.VIS_INK["faint"]), label=t, zorder=3)
            # Highlight where the top tier (Tier_1A) lands in steric space.
            _t1a = sub[sub["tier"] == CFG.TIER_TOP]
            if len(_t1a) >= 3:
                _t1c = TIER_PALETTE.get(CFG.TIER_TOP, CFG.VIS_ACCENT_DEEP["sea"])
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
            """
            Each axis is scaled to ITS OWN data. Forcing a shared limit and an equal aspect — which a
            y = x boundary line invites — sized the ligand axis by the POCKET volumes: cavities reach
            ~1,700 Å³ while almost every ligand sits below 500 Å³, so two thirds of the panel was
            empty and the data was crushed into a band along the bottom.

            The y = x boundary is a locus, not a 45° line: it marks where the ligand exactly fills the
            cavity, and it remains exactly that under independent scales — it simply becomes steep,
            which is the honest picture when ligands are an order of magnitude smaller than the
            pockets holding them. That is the finding, not a plotting artefact to be hidden by
            padding the axis until the line looks diagonal.
            """
            _xlim = float(np.ceil(sub["asv"].quantile(0.995) / 250.0) * 250.0)
            _ylim = float(np.ceil(sub["lv"].quantile(0.995) * 1.08 / 50.0) * 50.0)
            _lim = max(_xlim, _ylim)
            ax.plot([0, _lim], [0, _lim], color=CFG.VIS_INK["soft"], linestyle="--", linewidth=1.3,
                    alpha=0.8, zorder=4, label="Fit boundary (ligand = pocket)")
            ax.axvline(sub["asv"].median(), color=CFG.VIS_INK["ghost"], linestyle=":", linewidth=1.0, alpha=0.6, zorder=2)
            ax.axhline(sub["lv"].median(), color=CFG.VIS_INK["ghost"], linestyle=":", linewidth=1.0, alpha=0.6, zorder=2)
            ax.set_xlim(0, _xlim); ax.set_ylim(0, _ylim)
            ax.set_xlabel("Active-site cavity volume  (Å³)", )
            ax.set_ylabel("Ligand molecular volume  (Å³)", )
            _ovf_n = int((sub["lv"] > sub["asv"]).sum())
            _ovf = _ovf_n / len(sub) * 100
            ax.text(0.02, 0.98, f"{_ovf:.1f}% of complexes overfill the pocket "
                    f"(n={_ovf_n:,} above boundary)",
                    transform=ax.transAxes, ha="left", va="top", fontsize=8.5, style="italic",
                    color=CFG.VIS_ACCENT_DEEP["orange_deepest"], bbox=dict(boxstyle="round,pad=0.18", fc=CFG.VIS_ACCENT_DEEP["peach_fill"],
                    ec=_nofit_col, alpha=0.9, linewidth=0.8))
            # Two-row legend above the axes (tiers + fit boundary); the Spearman ρ/p/n
            # rides in the legend title so it sits with the key, not as a floating box.
            _ss01 = _spear_str(sub["asv"], sub["lv"], "XN-01 pocket volume vs ligand volume")
            ax.legend(loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=5,

                      title=_ss01, )
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
        """
        TWO coverage measures per carbon number, on one panel:

          · pocket_containment_cavity — the fraction of the ligand the PROTEIN CAVITY encloses,
            by ray-cast buriedness against every protein heavy atom. Drawn as the VIOLIN: it is the
            full distribution, and its shape is the story (a long low tail = the chain is spilling).
          · pocket_containment_site8 — the fraction of the ligand inside the contact shell of the
            EIGHT mapped catalytic residues (CFG §2.5, positions taken from the ranked sheet's
            Mapped_* columns, never hardcoded). Drawn as the BOX inside the violin: it is the
            reactive-shell engagement, and its median/IQR is what the tier gate cares about.

        Both means are traced as trend lines, so the divergence is readable at a glance: the cavity
        can still hold a long chain while the catalytic shell has already lost it. That gap IS the
        selectivity of the enzyme, and no single containment number can show it.

        Falls back to the pocket_occupancy view when the two containment columns are absent from the
        ranked sheet, so the panel still draws on a sheet that carries only the occupancy measure.
        """
        _has_cav02 = "pocket_containment_cavity" in df.columns
        _has_s802 = "pocket_containment_site8" in df.columns
        _two02 = _has_cav02 and _has_s802

        if _two02:
            sub = pd.DataFrame({"cav": _num("pocket_containment_cavity"),
                                "s8": _num("pocket_containment_site8"),
                                "lig": df.get(CFG.COL_LIG)}).dropna(subset=["cav", "s8", "lig"])
        else:
            sub = pd.DataFrame({"occ": _num("pocket_occupancy"),
                                "lig": df.get(CFG.COL_LIG)}).dropna(subset=["occ", "lig"])
        sub["_nC"] = sub["lig"].map(lambda l: _LP02.get(str(l), {}).get("nC"))
        sub = sub.dropna(subset=["_nC"])
        if not _LP02:
            reporter.log("  ! Diag 02 skipped: ligand SMILES panel not found for carbon grouping")
        elif _two02 and len(sub) >= 5 and sub["_nC"].nunique() >= 2:
            sub["_nC"] = sub["_nC"].astype(int)
            groups = sorted(sub["_nC"].unique())
            _present02 = set(sub["lig"].astype(str))

            def _xlab02b(nc):
                _fs = sorted({_LP02[l]["nF"] for l in _LP02
                              if _LP02[l].get("nC") == nc and l in _present02})
                return f"C{int(nc)}\n(" + ", ".join(f"{int(f)}F" for f in _fs) + ")"

            labels = [_xlab02b(g) for g in groups]
            fig, ax = plt.subplots(figsize=(max(9.5, 1.0 * len(groups)), 6.6))
            _CAV02, _S802 = CFG.VIS_ACCENT_DEEP["blue_light"], CFG.VIS_ACCENT_DEEP["tangerine"]

            _cav_by = [sub.loc[sub["_nC"] == g, "cav"].to_numpy() for g in groups]
            _s8_by = [sub.loc[sub["_nC"] == g, "s8"].to_numpy() for g in groups]
            _pos02 = np.arange(len(groups))

            _vp02 = ax.violinplot(_cav_by, positions=_pos02, widths=0.86,
                                  showmeans=False, showextrema=False)
            for _b in _vp02["bodies"]:
                _b.set_facecolor(_CAV02); _b.set_alpha(0.35)
                _b.set_edgecolor(_CAV02); _b.set_linewidth(0.9)

            _bp02 = ax.boxplot(_s8_by, positions=_pos02, widths=0.24, showfliers=False,
                               patch_artist=True,
                               medianprops=dict(color=CFG.VIS_RAMP["orange"][4], linewidth=1.5))
            for _b in _bp02["boxes"]:
                _b.set(facecolor=_S802, alpha=0.55, edgecolor=CFG.VIS_RAMP["orange"][4], linewidth=0.9)
            for _w in _bp02["whiskers"] + _bp02["caps"]:
                _w.set(color=CFG.VIS_RAMP["orange"][4], linewidth=0.9)

            _mcav = [float(np.mean(v)) if len(v) else np.nan for v in _cav_by]
            _ms8 = [float(np.mean(v)) if len(v) else np.nan for v in _s8_by]
            ax.plot(_pos02, _mcav, color=_CAV02, lw=2.2, marker="o", ms=6,
                    mec="white", mew=0.8, zorder=7, label="Mean — whole-cavity coverage")
            ax.plot(_pos02, _ms8, color=CFG.VIS_ACCENT_DEEP["orange_mid"], lw=2.2, ls="--", marker="D", ms=6,
                    mec="white", mew=0.8, zorder=7,
                    label="Mean — 8-residue active-site coverage")

            from matplotlib.patches import Patch as _P02
            ax.legend(handles=[
                _P02(facecolor=_CAV02, alpha=0.35, edgecolor=_CAV02,
                     label="Whole-cavity coverage (violin — full distribution)"),
                _P02(facecolor=_S802, alpha=0.55, edgecolor=CFG.VIS_RAMP["orange"][4],
                     label="8-residue active-site coverage (box — median · IQR)"),
                *ax.get_legend_handles_labels()[0],
            ], loc="lower left",   ncol=2)

            ax.set_ylim(0, 1.12)
            ax.set_yticks(np.arange(0, 1.01, 0.2))
            ax.axhline(1.0, color=CFG.VIS_INK["shadow"], ls=":", lw=1.1, alpha=0.8)
            ax.text(len(groups) - 0.45, 1.012, "ligand fully contained", ha="right", va="bottom",
                    fontsize=_JF_FA, color=CFG.VIS_INK["shadow"], style="italic")
            ax.set_xticks(_pos02); ax.set_xticklabels(labels, fontsize=8)
            ax.set_xlabel("PFAS carbon number  (fluorine counts present shown in parentheses)",
                          )
            ax.set_ylabel("Fraction of the ligand contained  (0–1)", )
            _tier_seps(ax, len(_pos02))   # consistent per-carbon separators
            ax.yaxis.grid(True, ls=":", alpha=0.35)
            ax.set_axisbelow(True)
            plt.tight_layout()
            _o = out_dir / "02_Pocket_Occupancy_by_Carbon_Number.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
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
                        medianprops=dict(color=CFG.VIS_ACCENT["bad"], linewidth=1.6),
                        whiskerprops=dict(linewidth=0.8), capprops=dict(linewidth=0.8))
            _med_xs02, _med_ys02 = [], []
            for i, g in enumerate(groups):
                m = sub.loc[sub["_nC"] == g, "occ"].median()
                _med_xs02.append(i); _med_ys02.append(float(m))
            ax.plot(_med_xs02, _med_ys02, color=CFG.VIS_ACCENT["bad"], linewidth=1.8, zorder=5, alpha=0.85)
            ax.scatter(_med_xs02, _med_ys02, marker="D", s=42, color=CFG.VIS_ACCENT["bad"],
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
                        transform=ax.transAxes, ha="right", va="top", fontsize=_JF_FA,
                        style="italic", color=CFG.VIS_INK["charcoal"])
            ax.set_xticks(range(len(groups))); ax.set_xticklabels(labels, fontsize=8)
            ax.set_xlabel("PFAS carbon number  (fluorine counts present shown in parentheses)", )
            ax.set_ylabel("Pocket occupancy  (ligand vol / cavity vol)", )
            _tier_seps(ax, len(groups))   # consistent per-carbon separators
            ax.legend(loc="upper left",   ncol=1)
            plt.tight_layout()
            _o = out_dir / "02_Pocket_Occupancy_by_Carbon_Number.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 02 skipped: insufficient occupancy/carbon data")
    except Exception as e:
        reporter.log(f"  ! Diag 02 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 03: Multi-model degrader consensus by tier ───────────────────────
    """
    Per-tier distribution of the Boltz-2 multi-model degrader consensus (fraction
    of the 5 diffusion models that independently called the complex a degrader),
    with strip overlay, tier mean (diamond), 95% CI bar and median annotation —
    the multi-model agreement signal underlying tier assignment.
    """
    try:
        cons = _num("model_degrader_consensus")
        sub = pd.DataFrame({"cons": cons, "tier": df.get(CFG.COL_TIER)}).dropna(subset=["cons"])
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
                ax.scatter([i], [s["mean"]], color=TIER_PALETTE.get(t, CFG.VIS_INK["faint"]), s=180,
                           zorder=5, edgecolors="black", linewidths=1.2)
            ax.set_ylim(-0.02, 1.06)
            ax.set_xlabel("Degrader tier", )
            ax.set_ylabel("Multi-model degrader consensus  (fraction of 5 models)", )
            _tier_ticklabels(ax, tiers, stats)
            ax.yaxis.grid(True, color=CFG.VIS_INK["tick"], linewidth=0.6, alpha=0.7, zorder=0)
            ax.set_axisbelow(True)
            _lh = [
                Line2D([0], [0], marker="o", color="none", markerfacecolor=CFG.VIS_INK["ghost"],
                       markeredgecolor="none", markersize=5, alpha=0.35, label="Individual complex"),
                Line2D([0], [0], marker="D", color="none", markerfacecolor=CFG.VIS_INK["ghost"],
                       markeredgecolor="black", markersize=8, markeredgewidth=1.2, label="Tier mean (◆)"),
                Line2D([0], [0], color="black", linewidth=2.0, label="95% confidence interval"),
            ]
            _mdhc = _md_ready_stars_cat(ax, df, tiers, "model_degrader_consensus")
            _lh.extend(_mdhc)
            # Legend above the axes (single row) so it never sits on the top data band.
            ax.legend(handles=_lh, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=4,
                        fancybox=True)
            plt.tight_layout()
            _tier_seps(plt.gca())   # consistent vertical tier separators
            _o = out_dir / "03_MultiModel_Consensus_by_Tier.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 03 skipped: insufficient consensus data")
    except Exception as e:
        reporter.log(f"  ! Diag 03 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 04: Model confidence by multi-model consensus level ──────────────
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
        conf = _num(CFG.COL_CONF)
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
                        medianprops=dict(color=CFG.VIS_ACCENT["bad"], linewidth=1.6),
                        whiskerprops=dict(linewidth=0.8), capprops=dict(linewidth=0.8),
                        zorder=3)
            _meds06 = [sub.loc[sub["cons"] == c, "conf"].median() for c in _levels]
            ax.plot(_order06, _meds06, color=CFG.VIS_ACCENT["bad"], linewidth=1.8, alpha=0.85, zorder=4,
                    marker="D", markersize=6, markeredgecolor="black", markeredgewidth=0.7,
                    label="Level median")
            ax.set_xticks(_order06)
            ax.set_xticklabels([f"{c:.2g}\n(n={int((sub['cons'] == c).sum()):,})" for c in _levels],
                               fontsize=8)
            for _tl06, _c06 in zip(ax.get_xticklabels(), _conf_cols):
                _tl06.set_color(_c06)   # each level's tick label matches its violin colour
            ax.set_xlabel("Multi-model degrader consensus  (fraction of 5 models)", )
            ax.set_ylabel("Boltz-2 model confidence", )
            _ss06 = _spear_str(sub["cons"], sub["conf"], "XN-06 pose consensus vs confidence")
            ax.legend(loc="lower right",
                      title=_ss06, )
            plt.tight_layout()
            _tier_seps(plt.gca())   # consistent vertical tier separators
            _o = out_dir / "04_Confidence_vs_Consensus.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 04 skipped: insufficient confidence/consensus data")
    except Exception as e:
        reporter.log(f"  ! Diag 04 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 05: Quality & competence diagnostics (2-panel scatter) ───────────
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
        if CFG.COL_CONF in df.columns and "competence_score" in df.columns:
            _qc_panels.append((CFG.COL_CONF, "competence_score",
                               "Boltz model confidence", "Catalytic competence score", None))
        _func_c = CFG.COL_MECH_S if CFG.COL_MECH_S in df.columns else (
            "soft_catalytic_score" if "soft_catalytic_score" in df.columns else None)
        if _func_c and "competence_score" in df.columns:
            _func_lbl = ("Mechanistic score (geometry / machinery)"
                         if _func_c == CFG.COL_MECH_S
                         else "Soft catalytic score (geometry / machinery)")
            _qc_panels.append((_func_c, "competence_score",
                               _func_lbl, "Catalytic competence score", None))
        if _qc_panels:
            """
            Both panels plot the SAME quantity on y — the catalytic competence score, on the same
            0–1 scale — so the axis is drawn once and shared. Repeating an identical scale beside
            itself spends horizontal space on nothing and invites the reader to check whether the two
            axes really are the same. Sharing it also locks the two clouds to a common vertical
            reference, which is the only way the panels can honestly be compared by eye.
            """
            fig, axes = plt.subplots(1, len(_qc_panels), sharey=True,
                                     figsize=(6.1 * len(_qc_panels), 5.6), squeeze=False)
            for _pi, (_axi, (_xc, _yc, _xl, _yl, _mode)) in enumerate(zip(axes[0], _qc_panels)):
                if _pi > 0:
                    _yl = None                 # the shared axis is labelled once, on the left
                _xv = _num(_xc)
                _yv = _num(_yc)
                if _mode == "invert":
                    _yv = (1.0 - _yv).clip(0, 1)
                _pdat = pd.DataFrame({"x": _xv, "y": _yv, "tier": df.get(CFG.COL_TIER)}).dropna(subset=["x", "y"])
                for _t in existing_tiers:
                    _d = _pdat[_pdat["tier"] == _t]
                    if _d.empty:
                        continue
                    _axi.scatter(_d["x"], _d["y"], s=16, alpha=0.35, linewidths=0,
                                 color=TIER_PALETTE.get(_t, CFG.VIS_INK["faint"]), label=_t, zorder=3)
                # Binned-median trend + relationship statistic for analytical depth.
                if len(_pdat) >= 20:
                    _qb = np.linspace(_pdat["x"].min(), _pdat["x"].max(), 11)
                    _pdat["_b"] = pd.cut(_pdat["x"], _qb, include_lowest=True)
                    _qm = _pdat.groupby("_b", observed=True).agg(
                        x=("x", "median"), y=("y", "median"), n=("y", "size")).dropna()
                    _qm = _qm[_qm["n"] >= 5]
                    if len(_qm) >= 2:
                        _axi.plot(_qm["x"], _qm["y"], color=CFG.VIS_ACCENT["bad"], linewidth=2.2,
                                  marker="D", markersize=6, markeredgecolor="white",
                                  markeredgewidth=0.7, zorder=6, label="Binned median")
                # Panel-specific Spearman sits INSIDE the panel (top-left whitespace), small and
                # unbold, so the two panels' statistics never collide above the axes.
                _ss07 = _spear_str(_pdat["x"], _pdat["y"], "XN-07 inter-model geometry spread")
                if _ss07:
                    _axi.text(0.035, 0.975, _ss07.replace("   ", "\n").replace("  ", " "),
                              transform=_axi.transAxes, ha="left", va="top",
                              fontsize=CFG.VIS_FONT_ANNOT - 0.5, fontweight="normal",
                              color=CFG.VIS_INK["ink_deep"], linespacing=1.3, zorder=8)
                _axi.set_xlabel(_xl, )
                if _yl:
                    _axi.set_ylabel(_yl, )
            # Single-row legend spanning the top of the figure (above both panels).
            _h07, _l07 = axes[0][0].get_legend_handles_labels()
            fig.legend(_h07, _l07, loc="lower left", bbox_to_anchor=(0.02, 0.99),
                       ncol=12)
            plt.tight_layout()
            for _qax in axes[0]: _tier_seps(_qax)   # consistent tier separators
            _o = out_dir / "05_Quality_and_Competence_Diagnostics.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 05 skipped: quality/competence columns absent")
    except Exception as e:
        reporter.log(f"  ! Diag 05 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig
    finally:
        if fig is not None:
            plt.close(fig)

    # ── Figure 06: Size preference — competence distribution & viability vs ligand size ──
    """
    FAcD is a small-substrate (haloacetate) hydrolase. Single-panel size-preference
    summary in the Figure-26c visual grammar: per-size-bin distribution of the effective
    mechanistic score (green boxes), the connected bin means (dark-orange size trend), the
    Tier_1A viability floor, and — on the right axis — the catalytic hit-rate (fraction
    reaching Tier_2A+) and the mean pocket_containment_cavity that drives the demotion. The
    reliable size axis is ligand_max_extent against the protein-aware cavity containment, NOT
    the convex-hull pocket_occupancy proxy.
    """
    try:
        from matplotlib.patches import Patch as _P08
        _ext = _num("ligand_max_extent"); _eff = _num("mechanistic_score_effective")
        _con = _num("pocket_containment_cavity")
        _torder = ["Tier_5_Decoy", "Tier_4", "Tier_3", "Tier_2B", "Tier_2A", "Tier_1B", "Tier_1A"]
        _tv = df[CFG.COL_TIER].astype(str).map({t: i for i, t in enumerate(_torder)})
        sub = pd.DataFrame({"ext": _ext, "eff": _eff, "con": _con, "tv": _tv}).dropna(subset=["ext", "eff", "tv"])
        if len(sub) >= CFG.VIS_DIAG_MIN_N_FOR_BINNING:
            _hi = float(np.nanpercentile(sub["ext"], 99))
            _bins = np.linspace(float(sub["ext"].min()), _hi, 13)
            sub["_b"] = pd.cut(sub["ext"], _bins, include_lowest=True)
            _NFLOOR08 = CFG.VIS_DIAG_MIN_N_PER_BIN
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
            _SOFT_C, _SOFT_F = CFG.VIS_RAMP["green"][2], CFG.VIS_RAMP["green"][1]
            _bp = axL.boxplot(_boxd, positions=_pos, widths=0.5, patch_artist=True,
                              showfliers=False, manage_ticks=False,
                              medianprops=dict(color=_SOFT_C, linewidth=1.7))
            for _b in _bp["boxes"]:
                _b.set(facecolor=_SOFT_F, alpha=0.85, edgecolor=_SOFT_C, linewidth=1.0)
            for _w in _bp["whiskers"] + _bp["caps"]:
                _w.set(color=_SOFT_C, linewidth=1.0)
            # Connected box means (dark-orange dashed diamonds — the size trend).
            _mln08, = axL.plot(_pos, _mean, color=CFG.VIS_RAMP["orange"][2], lw=2.2, ls="--", marker="D",
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
            axL.set_ylabel("Effective mechanistic score  (box = distribution)", )
            axL.set_xlabel("Ligand max extent (Å)  →  larger", )
            axL.set_xticks(_pos); axL.set_xticklabels([f"{c:.1f}" for c in _cent], rotation=90, fontsize=8)
            axL.spines["top"].set_visible(False)
            axL.yaxis.grid(True, color=CFG.VIS_ACCENT_DEEP["leaf"], alpha=0.28, linewidth=0.7, zorder=0)
            axL.set_axisbelow(True)
            # Right axis: catalytic hit-rate + mean pocket containment (both 0–1).
            axR = axL.twinx()
            _hln08, = axR.plot(_pos, _hit, color=CFG.VIS_ACCENT_DEEP["blue_mid"], lw=2.0, marker="o", markersize=6,
                               zorder=6, label="Catalytic hit-rate (fraction ≥ Tier_2A)")
            _cln08, = axR.plot(_pos, _cmn, color=CFG.VIS_ACCENT_DEEP["cyan"], lw=1.8, ls=(0, (5, 3)), marker="s",
                               markersize=5, zorder=5, label="Mean pocket containment")
            axR.set_ylim(0.0, 1.06); axR.grid(False)
            axR.set_ylabel("Catalytic hit-rate  ·  pocket containment", )
            from matplotlib.lines import Line2D as _L2D08
            _mln08.set_label("Mean effective mech score")
            _hln08.set_label("Catalytic hit-rate (≥ Tier_2A)")
            _handles08 = [_P08(facecolor=_SOFT_F, edgecolor=_SOFT_C,
                               label="Effective mech score (distribution)"),
                          _L2D08([0], [0], color=_SOFT_C, lw=1.7, label="Median (box centre line)"),
                          _mln08, _hln08, _cln08]
            axL.legend(handles=_handles08, loc="lower left", bbox_to_anchor=(0.0, 1.01),
                       ncol=5)
            plt.tight_layout()
            _o = out_dir / "06_Size_Preference_Containment.png"
            fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
            reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
        else:
            reporter.log("  ! Diag 06 skipped: insufficient size/tier data")
    except Exception as e:
        reporter.log(f"  ! Diag 06 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    # ── Figure 07: Reactive-centre engagement with the 8 catalytic residues vs size ──
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
        _lcol09 = next((c for c in [CFG.COL_LIG, CFG.COL_LIG, "ligand"] if c in df.columns), None)
        _need09 = _D8_09 + [CFG.COL_SN2, "mechanistic_score_effective", CFG.COL_TIER]
        if _lcol09 is None or not all(c in df.columns for c in _need09):
            reporter.log("  ! Diag 07 skipped: reactive-geometry / ligand columns absent")
        else:
            _smi09 = out_dir.parents[1] / CFG.INPUT_SMILES
            if not _smi09.exists():
                _smi09 = Path.cwd() / CFG.INPUT_SMILES
            _LP09 = _utils_mod.compute_ligand_properties(_smi09)
            if not _LP09:
                reporter.log(f"  ! Diag 07 skipped: ligand SMILES not found ({_smi09})")
                raise RuntimeError("ligand properties unavailable")
            _PROD_LO, _PROD_HI = CFG.VIS_DIAG_PROD_LO_A, CFG.VIS_DIAG_PROD_HI_A
            _READY, _SENT = CFG.VIS_DIAG_READY_DIST_A, CFG.VIS_DIAG_DIST_SENTINEL_A
            _ANG_MIN = CFG.VIS_DIAG_ANGLE_MIN_DEG
            _MECH_MIN = CFG.TIER_MECH_MIN["Tier_2A"]
            _d09 = df.copy()
            _d09["_nC"] = _d09[_lcol09].map(lambda l: _LP09.get(str(l), {}).get("nC"))
            _d8_09 = _d09[_D8_09].apply(pd.to_numeric, errors="coerce").mask(lambda s: s >= _SENT)
            _dnuc09 = _d8_09["Dist_Nucleophile"]
            _ang09 = _num(CFG.COL_SN2)
            _eff09 = _num("mechanistic_score_effective")
            _torder09 = ["Tier_5_Decoy", "Tier_4", "Tier_3", "Tier_2B", "Tier_2A", "Tier_1B", "Tier_1A"]
            _tv09 = _d09[CFG.COL_TIER].astype(str).map({t: i for i, t in enumerate(_torder09)})
            _dnraw09 = pd.to_numeric(_d09["Dist_Nucleophile"], errors="coerce")
            _ready09 = ((_dnraw09 <= _READY) & (_ang09 >= _ANG_MIN) & (_eff09 >= _MECH_MIN)).astype(float)
            _reach09 = (_dnraw09 <= _READY).astype(float)
            _sub09 = pd.DataFrame({"nC": _d09["_nC"], "dnuc": _dnuc09, "eff": _eff09,
                                   "tv": _tv09, "ready": _ready09, "reach": _reach09}).dropna(subset=["nC"])
            _grps09 = sorted(_sub09["nC"].unique())
            if len(_grps09) < 3:
                reporter.log("  ! Diag 07 skipped: insufficient carbon-number groups")
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
                _BAND09 = CFG.VIS_RAMP["purple"][1]     # 8-residue spread — distinct visible violet
                axL.fill_between(_pos09, _shlo, _shhi, color=_BAND09, alpha=0.22, zorder=0,
                                 label="8-residue distance spread (IQR)")
                axL.axhspan(_PROD_LO, _PROD_HI, color=CFG.VIS_ACCENT_DEEP["clover"], alpha=0.10, zorder=0)
                for _pe in (_PROD_LO, _PROD_HI):   # clear boundaries for the attack-distance window
                    axL.axhline(_pe, ls=(0, (5, 3)), color=CFG.VIS_BAND["high"], lw=1.1, alpha=0.6, zorder=1)
                axL.text(len(_pos09) - 0.42, (_PROD_LO + _PROD_HI) / 2,
                         "productive attack distance (2.5–3.5 Å)", ha="center", va="center",
                         rotation=90, fontsize=_JF_FA, color=CFG.VIS_BAND["high"], fontweight="bold")
                _SOFT_C09, _SOFT_F09 = CFG.VIS_RAMP["green"][2], CFG.VIS_RAMP["green"][1]
                _bp09 = axL.boxplot(_box09, positions=_pos09, widths=0.42, patch_artist=True,
                                    showfliers=False, manage_ticks=False,
                                    medianprops=dict(color=_SOFT_C09, linewidth=1.7))
                for _b in _bp09["boxes"]:
                    _b.set(facecolor=_SOFT_F09, alpha=0.9, edgecolor=_SOFT_C09, linewidth=1.0)
                for _w in _bp09["whiskers"] + _bp09["caps"]:
                    _w.set(color=_SOFT_C09, linewidth=1.0)
                _bmean09 = [float(np.nanmean(b)) if len(b) else np.nan for b in _box09]
                _mbl09, = axL.plot(_pos09, _bmean09, color=CFG.VIS_RAMP["green"][4], lw=2.0, ls="--", marker="D",
                                   ms=6, mec="white", mew=0.7, zorder=6, label="Mean nuc. distance")
                axL.set_ylim(2.0, 5.0); axL.set_xlim(-0.6, len(_pos09) - 0.15)
                axL.yaxis.set_major_locator(_MLoc09(0.25))
                axL.set_ylabel("Reactive-carbon → catalytic-residue distance (Å)", )
                axL.set_xlabel("PFAS carbon number  →  larger ligand", )
                axL.set_xticks(_pos09); axL.set_xticklabels([f"C{int(g)}" for g in _grps09], fontsize=8)
                axL.spines["top"].set_visible(False)
                axL.yaxis.grid(True, color=CFG.VIS_RAMP["blue"][2], alpha=0.28, lw=0.6, zorder=0)   # left (Å) grid
                axL.set_axisbelow(True)

                axR = axL.twinx()
                _rln09, = axR.plot(_pos09, _rdy09, color=CFG.VIS_ACCENT_DEEP["purple"], lw=2.6, marker="o", ms=6,
                                   zorder=8, label="Properly positioned")
                _rcl09, = axR.plot(_pos09, _rch09, color=CFG.VIS_INK["slate"], lw=1.5, ls=(0, (2, 2)),
                                   marker=".", ms=6, zorder=5, label="Reaches nucleophile")
                _mln09, = axR.plot(_pos09, _mm09, color=CFG.VIS_RAMP["orange"][2], lw=2.2, ls="--", marker="D",
                                   ms=6, mec="white", mew=0.7, zorder=6, label="Mean mech score")
                _hln09, = axR.plot(_pos09, _hit09, color=CFG.VIS_ACCENT_DEEP["blue_mid"], lw=2.0, ls=(0, (4, 2)),
                                   marker="s", ms=5, zorder=6, label="Hit-rate (≥Tier_2A)")
                for _fy, _fc, _lbl in [(CFG.TIER_MECH_MIN["Tier_1A"], CFG.TIER_COLOUR["Tier_1A"], "Tier 1A/1B"),
                                       (CFG.TIER_MECH_MIN["Tier_2A"], CFG.TIER_COLOUR["Tier_2A"], "Tier 2A")]:
                    axR.axhline(_fy, ls=":", color=_fc, lw=1.4, zorder=3)
                    axR.text(0.12, _fy + 0.004, f"{_lbl} floor ≥{_fy:g}", ha="left", va="bottom",
                             fontsize=_JF_FA, color=_fc, fontweight="bold", zorder=9,
                             bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=_fc, lw=0.5, alpha=0.88))
                """
                MD-selected complexes, starred on the left (distance) axis at their own carbon
                number. This panel's claim is that reactive engagement collapses as the ligand grows;
                the starred points are the cohort that claim is being acted on, so they are shown as
                individuals rather than dissolved into the box for their carbon bin.
                """
                _md09 = _md_ready_df(_d09)
                if not _md09.empty and "_nC" in _md09.columns:
                    _grp_pos09 = {int(g): p for g, p in zip(_grps09, _pos09)}
                    _mdc09 = pd.to_numeric(_md09["_nC"], errors="coerce")
                    _mdd09 = pd.to_numeric(_md09["Dist_Nucleophile"], errors="coerce")
                    _mx09, _my09 = [], []
                    for _cn09, _dv09 in zip(_mdc09, _mdd09):
                        if pd.isna(_cn09) or pd.isna(_dv09):
                            continue
                        if int(_cn09) in _grp_pos09:
                            _mx09.append(_grp_pos09[int(_cn09)])
                            _my09.append(float(_dv09))
                    if _mx09:
                        axL.scatter(_mx09, _my09, label=f"MD-selected (n={len(_mx09)})",
                                    **({**_MD_STAR_KW, "s": 300}))
                    _ctrl09 = _control_star_df(_d09)
                    if not _ctrl09.empty and "_nC" in _ctrl09.columns:
                        _cx09, _cy09 = [], []
                        for _ccn09, _cdv09 in zip(pd.to_numeric(_ctrl09["_nC"], errors="coerce"),
                                                  pd.to_numeric(_ctrl09["Dist_Nucleophile"], errors="coerce")):
                            if pd.isna(_ccn09) or pd.isna(_cdv09) or int(_ccn09) not in _grp_pos09:
                                continue
                            _cx09.append(_grp_pos09[int(_ccn09)]); _cy09.append(float(_cdv09))
                        if _cx09:
                            _ctrl_star(axL, _cx09, _cy09, size=320)
                            axL.scatter([], [], marker="*", s=200, facecolor=CFG.VIS_ACCENT["control"],
                                        edgecolor=CFG.VIS_ACCENT["control_edge"], linewidths=1.4, label="3R3U × FA (control)")

                axR.set_ylim(0, 1.02); axR.set_ylabel("Catalytic metric (0–1)", )
                axR.yaxis.set_major_locator(_MLoc09(0.1))
                axR.yaxis.grid(True, color=CFG.VIS_RAMP["orange"][1], alpha=0.30, lw=0.6, ls=(0, (4, 3)), zorder=0)  # right grid
                _leg09 = [_P09(facecolor=_SOFT_F09, edgecolor=_SOFT_C09, label="Nuc. distance (box)"),
                          _L2D09([0], [0], color=_SOFT_C09, lw=1.7, label="Median"),
                          _mbl09,
                          _P09(facecolor=_BAND09, alpha=0.22, label="8-residue spread (IQR)"),
                          _rln09, _rcl09, _mln09, _hln09]
                axL.legend(handles=_leg09, loc="lower left", bbox_to_anchor=(0.0, 1.01), ncol=8)
                plt.tight_layout()
                _o = out_dir / "07_Reactive_Engagement.png"
                fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
                reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
    except Exception as e:
        reporter.log(f"  ! Diag 07 skipped: {e}")
        plt.close("all")   # release the figure left open by the failed savefig

    _diag10_model_agreement(df, out_dir, reporter)


def _diag10_model_agreement(df: pd.DataFrame, out_dir: Path, reporter) -> None:
    """
    Diagnostic 10 — does the diffusion ensemble AGREE about the elite hits?

    Every complex is predicted several times (CFG.BOLTZ_DIFFUSION_SAMPLES independent diffusion
    samples). model_degrader_consensus is the fraction of those samples that independently reach a
    degrader tier. It answers the question the tier alone cannot: is a Tier_1A hit a reproducible
    property of the protein-ligand pair, or one lucky sample out of five?

    This matters more than it looks. The tier is assigned to the single representative pose, so a
    complex whose consensus is 0.2 got its elite label from ONE sample while four others disagreed —
    a coin-flip dressed as a result. If the elite tiers are not consensus-backed, the headline claim
    is a sampling artefact, and the reader is entitled to see that in a figure rather than infer it.

    Drawn as the consensus distribution per tier (violin + quartile box) with the per-tier mean
    annotated, so a tier whose hits are unanimous is visibly separated from one whose hits are not.
    """
    try:
        _c10 = "model_degrader_consensus"
        if _c10 not in df.columns or CFG.COL_TIER not in df.columns:
            reporter.log(f"  ! Diag 08 skipped: needs {_c10} + {CFG.COL_TIER}")
            return
        d10 = df[[_c10, CFG.COL_TIER]].copy()
        d10[_c10] = pd.to_numeric(d10[_c10], errors="coerce")
        d10 = d10.dropna(subset=[_c10])
        _t10 = [t for t in CFG.TIER_ORDER if (d10[CFG.COL_TIER] == t).sum() >= 3]
        if not _t10:
            reporter.log("  ! Diag 08 skipped: no tier with ≥3 complexes")
            return

        fig, ax = plt.subplots(figsize=(11, 6))
        _data10 = [d10.loc[d10[CFG.COL_TIER] == t, _c10].to_numpy() for t in _t10]
        _pos10 = np.arange(len(_t10))
        _vp10 = ax.violinplot(_data10, positions=_pos10, widths=0.8,
                              showmeans=False, showextrema=False)
        for _b10, _t in zip(_vp10["bodies"], _t10):
            _b10.set_facecolor(TIER_PALETTE.get(_t, CFG.VIS_INK["faint"]))
            _b10.set_alpha(0.55)
            _b10.set_edgecolor(CFG.VIS_INK["dark"])
            _b10.set_linewidth(0.7)
        _bp10 = ax.boxplot(_data10, positions=_pos10, widths=0.16, showfliers=False,
                           patch_artist=True, medianprops=dict(color="black", linewidth=1.4))
        for _b in _bp10["boxes"]:
            _b.set(facecolor="white", alpha=0.9, edgecolor=CFG.VIS_INK["dark"], linewidth=0.8)

        for _i10, (_t, _v10) in enumerate(zip(_t10, _data10)):
            _m10 = float(np.mean(_v10))
            ax.scatter([_i10], [_m10], marker="D", s=42, color=CFG.VIS_ACCENT["star_edge"],
                       edgecolor="white", linewidths=0.7, zorder=6)
            ax.text(_i10, 1.04, f"{_m10:.2f}\nn={len(_v10):,}", ha="center", va="bottom",
                    fontsize=_JF_FA, fontweight="bold", color=TIER_PALETTE.get(_t, CFG.VIS_INK["dark"]))

        ax.axhline(0.5, ls=":", lw=1.2, color=CFG.VIS_INK["grey"], zorder=2)
        ax.text(len(_t10) - 0.45, 0.51, "majority of samples agree", ha="right", va="bottom",
                fontsize=_JF_FA, color=CFG.VIS_INK["grey"], style="italic")
        ax.set_xticks(_pos10)
        ax.set_xticklabels(_t10, rotation=30, ha="right", fontsize=9)
        for _tk, _t in zip(ax.get_xticklabels(), _t10):
            _tk.set_color(TIER_PALETTE.get(_t, "black"))
        ax.set_ylim(0, 1.18)
        ax.set_yticks(np.arange(0, 1.01, 0.2))
        ax.set_ylabel("Diffusion-sample consensus\n(fraction of samples independently reaching a degrader tier)",
                      )
        ax.set_xlabel("Degrader Tier", )
        ax.yaxis.grid(True, ls=":", alpha=0.35)
        ax.set_axisbelow(True)
        _md10 = _md_ready_df(df)
        if not _md10.empty and _c10 in _md10.columns:
            for _t, _g in _md10.groupby(_md10[CFG.COL_TIER].astype(str)):
                if _t not in _t10:
                    continue
                _x10 = _t10.index(_t)
                _y10 = pd.to_numeric(_g[_c10], errors="coerce").dropna()
                if len(_y10):
                    ax.scatter(np.full(len(_y10), _x10), _y10.to_numpy(),
                               **({**_MD_STAR_KW, "s": 200}))
        _ctrl10 = _control_star_df(df)
        if not _ctrl10.empty and _c10 in _ctrl10.columns:
            for _t, _g in _ctrl10.groupby(_ctrl10[CFG.COL_TIER].astype(str)):
                if _t not in _t10:
                    continue
                _cy10 = pd.to_numeric(_g[_c10], errors="coerce").dropna()
                if len(_cy10):
                    _ctrl_star(ax, np.full(len(_cy10), _t10.index(_t)), _cy10.to_numpy(), size=220)
        plt.tight_layout()
        _tier_seps(plt.gca())   # consistent vertical tier separators
        _o = out_dir / "08_Model_Agreement.png"
        fig.savefig(_o, dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight"); plt.close(fig)
        reporter.log(f"  ✔ Saved: {_o.parent.name}/{_o.name}")
    except Exception as e:                                   # noqa: BLE001
        reporter.log(f"  ! Diag 08 skipped: {e}")
        plt.close("all")


# ===============================================================================
# SECTION 4E: PER-FOLDER ANALYSIS PANELS  (folded into their folder steps by the orchestrator)
# ===============================================================================
"""
Eight panels merged in from the two more_Plots prototypes, each taken from whichever prototype drew
it best (both were rendered and reviewed side by side before selection).

The prototypes define helpers of the SAME NAME with DIFFERENT bodies (_save, _legend_with_stats), and
each defines a phylogeny figure — of which BOTH are wanted. Everything is therefore namespaced by
origin (_xn_ = newer prototype, _xo_ = older) rather than folded into one set of helpers: a shared
namespace would silently hand one version's figure the other version's helper, and the result would
look plausible and be wrong.

The prototypes' raw filenames (Figure_1_…, Figure_05C_…, mixed numbering and letter suffixes) are
replaced by the pipeline's convention — sequential 01-08 with descriptive names.
"""

_xn_STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

_xn_WATER_NAMES = {"HOH", "WAT", "SOL", "TIP3", "TIP3P", "SPC", "SPCE", "OPC"}

_xn_NUCLEOPHILE_SEARCH_RADIUS_A = 6.0

_xn_CF_BOND_CUTOFF_A = 1.7

_xn_IDEAL_SN2_ANGLE_DEG = 180.0

"""
FAcD attacks with an ASPARTATE and nothing else. The carboxylate's two oxygens are the only atoms that
can be the nucleophile; WHICH aspartate it is varies by variant, and 02 records that per complex as
Mapped_Nucleophile in the ranked CSV.

There is no table of candidate residue types here on purpose. A generic ASP/GLU/SER/THR/CYS list — the
kind a general-purpose structure tool ships with — let a SERINE win the geometry contest on 4 of 30
Tier_4 complexes (SER147, SER149, SER163). A serine hydroxyl is a different enzyme's mechanism. The
angle reported was geometrically real and mechanistically meaningless, which is the worst kind of wrong:
it looks like data.
"""
_xn_NUCLEOPHILE_ATTACK_ATOMS = {"ASP": ("OD1", "OD2"), "ASH": ("OD1", "OD2")}

@dataclass
class _xn_AtomRecord:
    """An atom plus the minimum metadata needed to locate it: residue, sequence id, chain."""
    atom: gemmi.Atom
    residue_name: str
    residue_seqid: int
    chain_name: str

@dataclass
class _xn_GeometryResult:
    """The best SN2 geometry found in one model."""
    sn2_distance_A: float
    sn2_angle_deg: float
    sn2_angle_deviation_deg: float
    nucleophile_resname: str
    nucleophile_resseq: int
    nucleophile_chain: str
    nucleophile_atom: str
    reactive_carbon_atom: str
    leaving_fluorine_atom: str
    cf_distance_A: float
    candidate_count: int
    status: str

def _xn_calculate_angle(p1: gemmi.Position, p2: gemmi.Position, p3: gemmi.Position) -> float:
    """The p1-p2-p3 angle in degrees.

    For the SN2: p1 = the attacking Oδ of the catalytic aspartate, p2 = the ligand's α-carbon,
    p3 = the leaving fluorine. 180° is a perfect backside attack.
    """
    v1 = np.array([p1.x - p2.x, p1.y - p2.y, p1.z - p2.z], dtype=float)
    v2 = np.array([p3.x - p2.x, p3.y - p2.y, p3.z - p2.z], dtype=float)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    if norm_v1 == 0 or norm_v2 == 0:
        return 0.0
    cosine_angle = np.dot(v1, v2) / (norm_v1 * norm_v2)
    angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    return float(np.degrees(angle))

def _xn_is_ligand_residue(residue: gemmi.Residue) -> bool:
    """Anything that is not a standard amino acid or water is treated as the ligand."""
    name = residue.name.strip().upper()
    return name not in _xn_STANDARD_AA and name not in _xn_WATER_NAMES

def _xn_iter_atoms_from_cif(cif_path: Path) -> Tuple[List[_xn_AtomRecord], List[_xn_AtomRecord]]:
    """
    Split a Boltz-2 CIF into its ligand atoms and its protein atoms.

    A parse failure is propagated to the caller, which skips the
    singolo modello senza interrompere l'intero batch.
    """
    doc = gemmi.cif.read_file(str(cif_path))
    structure = gemmi.make_structure_from_block(doc.sole_block())
    ligand_atoms: List[_xn_AtomRecord] = []
    protein_atoms: List[_xn_AtomRecord] = []

    for model in structure:
        for chain in model:
            for residue in chain:
                bucket = ligand_atoms if _xn_is_ligand_residue(residue) else protein_atoms
                for atom in residue:
                    # Hydrogens play no part in this geometry; skipping them keeps the parse cheap.
                    if atom.element.name == "H":
                        continue
                    bucket.append(_xn_AtomRecord(
                        atom=atom,
                        residue_name=residue.name.strip().upper(),
                        residue_seqid=int(residue.seqid.num) if residue.seqid.num is not None else -1,
                        chain_name=chain.name,
                    ))
    return ligand_atoms, protein_atoms

def _xn_find_reactive_cf_pairs(lig_atoms: Sequence[_xn_AtomRecord]) -> List[Tuple[_xn_AtomRecord, _xn_AtomRecord]]:
    """Every bonded C–F pair in the ligand (closer than the C–F bond cutoff).

    The scissile bond is one of these: the α-carbon holds the fluorine that leaves.
    """
    carbons = [a for a in lig_atoms if a.atom.element.name == "C"]
    fluorines = [a for a in lig_atoms if a.atom.element.name == "F"]
    pairs: List[Tuple[_xn_AtomRecord, _xn_AtomRecord]] = []
    for carbon in carbons:
        for fluorine in fluorines:
            if carbon.atom.pos.dist(fluorine.atom.pos) < _xn_CF_BOND_CUTOFF_A:
                pairs.append((carbon, fluorine))
    return pairs

def _xn_iter_nucleophile_atoms(protein_atoms: Sequence[_xn_AtomRecord]) -> Iterable[_xn_AtomRecord]:
    """The Oδ atoms of the catalytic aspartate. `protein_atoms` is already confined to that residue."""
    for atom_record in protein_atoms:
        allowed_atoms = _xn_NUCLEOPHILE_ATTACK_ATOMS.get(atom_record.residue_name)
        if not allowed_atoms:
            continue
        if atom_record.atom.name.strip() in allowed_atoms:
            yield atom_record

def _xn_geometry_rank(distance_A: float, angle_deg: float) -> Tuple[float, float]:
    """Rank the two carboxylate oxygens: short distance AND an angle near 180°.

    The first term combines the distance with the normalised angular deviation; the second keeps the
    distance as a stable tie-break when two candidates score alike.
    """
    angle_deviation = abs(_xn_IDEAL_SN2_ANGLE_DEG - angle_deg)
    combined = distance_A + (angle_deviation / 45.0)
    return combined, distance_A

def _xn_compute_sn2_geometry(cif_path: Path, nuc_resnum: int | None = None) -> _xn_GeometryResult:
    """The best SN2 geometry in one model, measured AT THE CATALYTIC NUCLEOPHILE.

    THE NUCLEOPHILE IS ONE RESIDUE. FAcD attacks with an ASPARTATE — the residue 02 mapped for this
    protein by global alignment (Mapped_Nucleophile: Asp110, Asp112, Asp109 …). Scanning the whole
    protein for whichever ASP/GLU/SER/THR/CYS scores best on distance and angle lets a residue that does
    no catalysis win the contest: measured on the Tier_4 cohort it picked a SERINE in 4 of 30 complexes
    (SER147, SER149, SER163, at 3-4 A). A serine hydroxyl is a different enzyme's mechanism. The angle
    it reports is geometrically real and mechanistically meaningless.

    So the residue is passed in and the search is confined to its two carboxylate oxygens. Without it the
    function REFUSES: there is no fallback scan, because the fallback scan is what produced the serine.
    """
    if nuc_resnum is None:
        # No mapped aspartate, no measurement. Falling back to a whole-protein scan is precisely what
        # put a serine in the results; an absent number is honest, a confident wrong one is not.
        raise ValueError("no mapped catalytic nucleophile for this complex")
    lig_atoms, protein_atoms = _xn_iter_atoms_from_cif(cif_path)
    protein_atoms = [a for a in protein_atoms if a.residue_seqid == int(nuc_resnum)]
    cf_pairs = _xn_find_reactive_cf_pairs(lig_atoms)
    if not lig_atoms:
        return _xn_empty_geometry("no_ligand_atoms")
    if not cf_pairs:
        return _xn_empty_geometry("no_reactive_cf_pair")

    best: Optional[Tuple[Tuple[float, float], _xn_GeometryResult]] = None
    candidate_count = 0

    for nuc in _xn_iter_nucleophile_atoms(protein_atoms):
        for carbon, fluorine in cf_pairs:
            distance_A = nuc.atom.pos.dist(carbon.atom.pos)
            if distance_A > _xn_NUCLEOPHILE_SEARCH_RADIUS_A:
                continue
            angle_deg = _xn_calculate_angle(nuc.atom.pos, carbon.atom.pos, fluorine.atom.pos)
            angle_deviation = abs(_xn_IDEAL_SN2_ANGLE_DEG - angle_deg)
            cf_distance_A = carbon.atom.pos.dist(fluorine.atom.pos)
            candidate_count += 1
            result = _xn_GeometryResult(
                sn2_distance_A=distance_A,
                sn2_angle_deg=angle_deg,
                sn2_angle_deviation_deg=angle_deviation,
                nucleophile_resname=nuc.residue_name,
                nucleophile_resseq=nuc.residue_seqid,
                nucleophile_chain=nuc.chain_name,
                nucleophile_atom=nuc.atom.name.strip(),
                reactive_carbon_atom=carbon.atom.name.strip(),
                leaving_fluorine_atom=fluorine.atom.name.strip(),
                cf_distance_A=cf_distance_A,
                candidate_count=candidate_count,
                status="ok",
            )
            ranked = _xn_geometry_rank(distance_A, angle_deg)
            if best is None or ranked < best[0]:
                best = (ranked, result)

    if best is None:
        return _xn_empty_geometry("no_nucleophile_within_6A")

    final = best[1]
    return _xn_GeometryResult(
        **{**final.__dict__, "candidate_count": candidate_count}
    )

def _xn_empty_geometry(status: str) -> _xn_GeometryResult:
    """The sentinel result for a model with no measurable SN2 geometry — every field NaN, never 0.0."""
    return _xn_GeometryResult(
        sn2_distance_A=math.nan,
        sn2_angle_deg=math.nan,
        sn2_angle_deviation_deg=math.nan,
        nucleophile_resname="",
        nucleophile_resseq=-1,
        nucleophile_chain="",
        nucleophile_atom="",
        reactive_carbon_atom="",
        leaving_fluorine_atom="",
        cf_distance_A=math.nan,
        candidate_count=0,
        status=status,
    )

_xn_ELITE_TIERS = ["Tier_1A", "Tier_1B"]

# Panel width. The seven tier names are printed in full and upright; at a narrower width
# 'Tier_4' and 'Tier_5_Decoy' overlap and the axis becomes unreadable at the one place the
# figure is decided. The width is set so the longest name fits without rotating it.
_xn_COL_DOUBLE_IN = 14.5

# Set from --no-variance in main(). Declared here so the geometry figure can read it even when a panel
# is invoked directly rather than through a full run.
_ALLOW_VARIANCE_COMPUTE = True

_xn_STRIP_MAX_PER_GROUP = 250

_xn__RNG = np.random.default_rng(int(CFG.ANALYSIS_SEED))   # CFG has no RANDOM_SEED — that branch was always 42



def _xn__nuc_distance(df: pd.DataFrame) -> pd.Series:
    """Nucleophile–substrate distance, resolving the several historical names.

    The Step-02 writer renames ``dist_Nuc`` to ``Dist_Nucleophile_ASP10`` via
    ``CFG.COLUMN_RENAMING_MAP``; older/raw exports expose the live geometry as
    ``best_nucleophile_distance``. All three are treated as equivalent.
    """
    col = _xn__col(df, 'Dist_Nucleophile_ASP10', 'Dist_Nucleophile', 'best_nucleophile_distance', 'dist_Nuc')
    return _xn__num(df, col)

_xn_PILLAR_ALIASES = {
    'Model_Quality_Score': ['Model_Quality_Score', 'Model_Quality_ScoreNormalised', CFG.COL_CONF],
    'Binding_Affinity_Score': ['Binding_Affinity_Score', 'Binding_Affinity_ScoreNormalised', 'Chemical_Affinity_Score', 'Binding_Probability', 'custom_affinity_score', 'Binding_Probability_Score'],
    'Catalytic_Competence_Score': ['Catalytic_Competence_Score', 'Catalytic_Competence_ScoreNormalised', 'competence_score', 'soft_catalytic_score'],
    'Evolutionary_Fingerprint_Score': ['Evolutionary_Fingerprint_Score', 'Evolutionary_Fingerprint_ScoreNormalised', CFG.COL_LIKE_S, CFG.COL_ID_PCT],
    'Final_Unified_Score': ['Final_Unified_Score', 'Ranking_Score_Calc'],
}

def _xn__col(df: pd.DataFrame, *candidates: str) -> str | None:
    """Return the first candidate column that exists in ``df`` (else ``None``)."""
    for c in candidates:
        if c in df.columns:
            return c
    return None

def _xn__pillar_col(df: pd.DataFrame, pillar: str) -> str | None:
    """Resolve a headline-pillar name to whichever concrete column is present."""
    return _xn__col(df, *_xn_PILLAR_ALIASES.get(pillar, [pillar]))

def _xn__num(df: pd.DataFrame, col: str | None) -> pd.Series:
    """Coerce a column to numeric, returning an empty series when it is absent."""
    if col is None or col not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[col], errors='coerce')


def _xn__fmt_p(p: float) -> str:
    """Compact, publication-style p-value formatting."""
    if p is None or not np.isfinite(p):
        return 'p = n/a'
    if p < 0.0001:
        return 'p < 1e-4'
    return f'p = {p:.3g}'



def _ext_match_03_style(fig) -> None:
    """Bring an extended-analysis panel onto 03's typography before it is written.

    The merged prototypes never set label or tick fonts — they inherited whatever rcParams were in
    force, which in their own sandbox was seaborn's scaled theme and here is 03's. The result is a
    folder whose axis labels are a different size and colour from every other figure in the set, and
    a reader flicking between folders sees two typefaces and wonders which is authoritative.

    The values are the ones 03's own figures pass explicitly: 11 pt axis labels, 9 pt ticks, black.
    Applied at save time so it reaches every axis of every panel, including the ones the prototypes
    build through seaborn and never touch again.
    """
    _LBL, _TICK, _COL = 11.0, 9.0, CFG.VIS_INK["black"]
    _LEG = float(getattr(CFG, "VIS_FONT_LEGEND", 8.5))
    for _ax in fig.get_axes():
        # A twin (secondary) axis may be deliberately colour-coded so its quantity cannot be misread as
        # the other axis's: it carries a `_twin_axis_colour` tag and keeps that colour on its y-label and
        # y-ticks. The x-axis stays black everywhere (tier ticks are recoloured separately below).
        _ycol = getattr(_ax, "_twin_axis_colour", None) or _COL
        _ax.xaxis.label.set_fontsize(_LBL); _ax.xaxis.label.set_color(_COL); _ax.xaxis.label.set_fontweight("normal")
        _ax.yaxis.label.set_fontsize(_LBL); _ax.yaxis.label.set_color(_ycol); _ax.yaxis.label.set_fontweight("normal")
        """
        A long y-label set at the same size as a short one runs the height of the panel and forces the
        reader to track it vertically. Labels past the threshold step down a point size; short ones are
        untouched, so the axis-label size stays uniform across the set wherever it can be.
        """
        _yl = _ax.yaxis.label.get_text()
        if len(_yl) > 34:
            _ax.yaxis.label.set_fontsize(_LBL - 1.5)
        _ax.tick_params(axis="x", labelsize=_TICK, labelcolor=_COL)
        _ax.tick_params(axis="y", labelsize=_TICK, labelcolor=_ycol)
        for _t in _ax.get_yticklabels():
            _t.set_fontsize(_TICK)
            _t.set_color(_ycol)

        """
        A tier tick takes its TIER'S colour, as it does in every other figure in the set. The tick
        and the box above it then say the same thing twice, in the same language: the reader can go
        from a colour in the plot to a name on the axis without consulting a key. Ticks that are not
        tier names (a carbon number, a ligand) stay black — colouring them would imply an encoding
        that does not exist.
        """
        for _t in _ax.get_xticklabels():
            _t.set_fontsize(_TICK)
            _raw = _t.get_text().split("\n")[0].strip()
            if _raw in TIER_PALETTE:
                _t.set_color(TIER_PALETTE[_raw])
            else:
                _t.set_color(_COL)

        # Legend text on the pipeline's legend size, so a legend here and a legend in 03 are the
        # same object at the same scale.
        _lg = _ax.get_legend()
        if _lg is not None:
            for _txt in _lg.get_texts():
                _txt.set_fontsize(_LEG)
            if _lg.get_title() is not None:
                _lg.get_title().set_fontsize(_LEG)

    """
    When every panel of a multi-panel figure carries the SAME x-label, the label is a property of the
    figure, not of each panel. Printing it under all of them repeats one fact n times and eats the
    width the tier names need. It is lifted to a single figure-level label instead.

    Only when they are identical: panels that measure different quantities keep their own labels.
    """
    _axes = [_a for _a in fig.get_axes() if _a.get_visible() and _a.get_xlabel()]
    if len(_axes) > 1:
        _labels = {_a.get_xlabel() for _a in _axes}
        if len(_labels) == 1:
            for _a in _axes:
                _a.set_xlabel("")
            fig.supxlabel(_labels.pop(), fontsize=_LBL, color=_COL)


def _xn__save(fig, out_dir: Path, name: str, reporter) -> None:
    """Persist a figure as a PNG at the pipeline's publication resolution, then free its memory.

    The DPI comes from CFG.VIS_FIGURE_DPI — the same value every other step renders at, so a
    sandbox figure and a Step-03 figure are the same physical object at the same scale.
    """
    _ext_match_03_style(fig)
    stem = Path(name).stem
    png = out_dir / f'{stem}.png'
    fig.savefig(png, dpi=CFG.VIS_FIGURE_DPI)
    plt.close(fig)

def _xn__panel(ax, letter: str) -> None:
    """No-op: panel letters are not drawn.

    No other figure in this set carries an (A)/(B) tag, and a lettering convention that appears in one
    figure and nowhere else is noise — the reader assumes it means something.
    """
    return

def _xn__fmt_n(n: int) -> str:
    """A sample size that fits the ~0.5 in a tier occupies on a two-panel figure.

    Seven tiers across a double-column axis leave no room for '22,184' under each tick; the
    thousands are what a reader takes from it, not the units.
    """
    n = int(n)
    return f'{n:,}' if n < 1000 else f'{n / 1000:.1f}k'

def _xn__stat_header(ax, text: str) -> None:
    """The panel's test result, on its own line ABOVE the axes.

    Inside the axes it has nowhere to go that is not on top of either the data or the tail of a
    violin, and a boxed annotation nudged out of the way by constrained_layout collides with the
    neighbouring panel instead. Above the frame it always has room and never covers a mark.
    """
    if not text:
        return
    ax.text(0.015, 0.015, text, transform=ax.transAxes, ha='left', va='bottom',
            fontsize=CFG.VIS_FONT_ANNOT - 0.5, color=CFG.VIS_INK["soft"],
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=CFG.VIS_LEGEND_FRAME_ALPHA,
                      edgecolor=CFG.VIS_INK["palest"], linewidth=0.6), zorder=9)

def _xn__thin(sub: pd.DataFrame, group_col: str, cap: int = _xn_STRIP_MAX_PER_GROUP) -> pd.DataFrame:
    """A per-group random subsample for strip overlays (see _xn_STRIP_MAX_PER_GROUP)."""
    parts = []
    for _, g in sub.groupby(group_col, observed=True):
        if len(g) > cap:
            g = g.iloc[_xn__RNG.choice(len(g), cap, replace=False)]
        parts.append(g)
    return pd.concat(parts, ignore_index=True) if parts else sub

def _xn__kruskal(sub: pd.DataFrame, group_col: str, val_col: str, order, reg_name: str = None) -> str:
    """Kruskal-Wallis across the ordered groups, formatted for an on-panel annotation.

    reg_name is the name the test is REGISTERED under (defaults to val_col). Callers that pass a
    renamed working column (e.g. every pillar copied to 'y') must pass the real column here, or all
    their tests collapse onto one registry key ('Kruskal-Wallis — y across tier') and only the last
    survives the BH-family dedup."""
    _reg = reg_name or val_col
    groups = [sub.loc[sub[group_col] == t, val_col].dropna().values for t in order]
    groups = [g for g in groups if len(g) >= 2]
    if len(groups) < 2:
        return ''
    H, p = _sc_stats.kruskal(*groups)
    # Epsilon-squared: the share of rank variance the grouping explains. H alone grows with n
    # and says nothing about how large the tier separation actually is.
    n = sum(len(g) for g in groups)
    eps2 = (H - len(groups) + 1) / (n - len(groups)) if n > len(groups) else np.nan
    # Register so this on-panel p pays the same BH multiplicity toll as the battery's tests and appears
    # in 06_Statistical_Tests.csv — otherwise a reader cannot discover it was run. reg_name keeps the
    # call sites (geometry / binding / the four pillars) distinct in the family.
    _register_p(f"Kruskal-Wallis — {_reg} across {group_col}", "03_Extended_Panels",
                float(H), int(n), float(p),
                eps_sq=(round(float(eps2), 4) if np.isfinite(eps2) else None))
    return f'Kruskal–Wallis  H = {H:,.0f}   {_xn__fmt_p(p)}   ε² = {eps2:.2f}'

# The extended panels are filed in the thematic subfolders alongside the figures they belong with,
# numbered to continue each destination folder's own sequence.
_xn_FIG_NAMES = {
    'geometry':   '13_Geometry_and_Uncertainty',
    'binding':    '08_Binding_Affinity_Metrics',
    'phylogeny':  '05_Evolutionary_Phylogeny',
    'pillars':    '09_Pillar_Divergence_by_Tier',
    'size_mech':  '14_Mechanistic_Breakdown_by_Tier',
    'size_tier':  '14_Chain_Length_by_Tier',
    'validation': '15_Tier1A_Cross_Ligand_Heatmap',
}

def _xn__gate_lines_distance(ax) -> None:
    """The tier distance gates, drawn where the tiers are actually decided (CFG.TIER_NUC_DIST).

    A box plot shows where the data sit; it does not show the cut the classifier applied. With
    the gate drawn, a reader can see for themselves why a complex landed in the tier it did.
    """
    gates = getattr(CFG, 'TIER_NUC_DIST', {}) or {}
    # The 1A and 2A gates are only 0.2 A apart: centred labels would print on top of each other,
    # so consecutive labels sit alternately below and above their own rule.
    _gate_keys = []
    for i, tier in enumerate((CFG.TIER_TOP, 'Tier_2A', 'Tier_2B')):
        v = gates.get(tier)
        if v is None or not np.isfinite(float(v)) or float(v) > 12.0:
            continue
        ax.axhline(float(v), ls=':', lw=1.0, color=TIER_PALETTE.get(tier, CFG.VIS_INK["grey"]),
                   alpha=0.9, zorder=4)
        """
        The rule carries a LETTER, not a sentence. 'Tier_1A ≤3.0 Å' written along its line is wider
        than the gap between the lines it has to fit into — Tier_1A (3.0 Å) and Tier_2A (3.2 Å) are
        0.2 Å apart — so whatever it clears, it lands on something else. A single letter fits between
        them, and the key at the foot of the panel says what each letter means.
        """
        _tag = chr(ord('A') + i)
        _gate_keys.append(f'{_tag} = {tier} ≤ {float(v):.1f} Å')
        ax.text(0.995, float(v), _tag, transform=ax.get_yaxis_transform(),
                ha='right', va='center', fontsize=CFG.VIS_FONT_ANNOT,
                fontweight='bold', color=TIER_PALETTE.get(tier, CFG.VIS_INK["grey"]),
                bbox=dict(boxstyle='circle,pad=0.16', facecolor='white', alpha=0.92,
                          edgecolor=TIER_PALETTE.get(tier, CFG.VIS_INK["grey"]), linewidth=0.8),
                zorder=7)

    # The key is NOT drawn here. Two boxes — one naming the tier colours, one naming the gate tags —
    # say two halves of the same thing and cost the panel twice the space. The strings are handed to
    # the caller, which folds them into the single legend.
    ax._gate_keys = _gate_keys

def _xn__gate_lines_angle(ax) -> None:
    """The tier SN2-angle gates (CFG.TIER_ANGLE_MIN); see _xn__gate_lines_distance."""
    gates = getattr(CFG, 'TIER_ANGLE_MIN', {}) or {}
    _gate_keys = []
    for i, tier in enumerate((CFG.TIER_TOP, 'Tier_2A', 'Tier_2B')):
        v = gates.get(tier)
        if v is None or not np.isfinite(float(v)):
            continue
        ax.axhline(float(v), ls=':', lw=1.0, color=TIER_PALETTE.get(tier, CFG.VIS_INK["grey"]),
                   alpha=0.9, zorder=4)
        _tag = chr(ord('A') + i)
        _gate_keys.append(f'{_tag} = {tier} ≥ {float(v):.0f}°')
        ax.text(0.995, float(v), _tag, transform=ax.get_yaxis_transform(),
                ha='right', va='center', fontsize=CFG.VIS_FONT_ANNOT,
                fontweight='bold', color=TIER_PALETTE.get(tier, CFG.VIS_INK["grey"]),
                bbox=dict(boxstyle='circle,pad=0.16', facecolor='white', alpha=0.92,
                          edgecolor=TIER_PALETTE.get(tier, CFG.VIS_INK["grey"]), linewidth=0.8),
                zorder=7)

    # The key is NOT drawn here. Two boxes — one naming the tier colours, one naming the gate tags —
    # say two halves of the same thing and cost the panel twice the space. The strings are handed to
    # the caller, which folds them into the single legend.
    ax._gate_keys = _gate_keys

def _xn__tidy_tier_ticks(ax, counts=None) -> None:
    """Tier labels upright and shortened, optionally carrying their own n.

    Rotated labels cost a reader a head-tilt per panel, and the 'Tier_' stem repeats on every
    tick without adding anything the axis title lacks. Where the top of the axis is already
    spoken for by a legend, the sample size rides on the tick label rather than fighting it.
    """
    labels = []
    for t in ax.get_xticklabels():
        raw = t.get_text()
        # The full tier name, as every other figure in the set prints it. Abbreviating to '1A' saves
        # a few pixels and costs the reader a translation on every glance between figures.
        lbl = raw
        if counts is not None:
            lbl += f'\n{_xn__fmt_n(counts.get(raw, 0))}'
        labels.append(lbl)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, rotation=0, fontsize=CFG.VIS_FONT_TICK)
    if counts is not None:
        ax.set_xlabel(f'{ax.get_xlabel()}   (n below each tier)')

def _xn__mean_trend(ax, sub, tiers, val_col, *, group_col='tier', label='Mean (trend)'):
    """Connect the per-tier MEANS with a line.

    Seven boxes side by side ask the reader to compare seven distributions by eye and infer whether the
    quantity rises or falls with tier. The trend line states that in one stroke — and it must be the MEAN,
    not the median already drawn inside each box, or the line would only restate a mark that is there.
    """
    _m = [pd.to_numeric(sub.loc[sub[group_col] == _t, val_col], errors='coerce').mean() for _t in tiers]
    ax.plot(range(len(tiers)), _m, color=CFG.VIS_ACCENT["bad"], lw=2.0, marker='D', ms=5.5,
            mec='white', mew=0.8, zorder=9, label=label)
    ax.legend(loc='best')   # font + frame inherited from apply_figure_style (SSOT)
    return ax


def _xn__tier_boxstrip(ax, sub, tiers, val_col, ylabel, *, group_col='tier'):
    """The tier-versus-value panel used by several figures: box + capped strip + n + test.

    Draw order matters. The strip goes down FIRST at low zorder and the box on top of it with a
    semi-transparent face, so the median and the quartiles stay readable through the points
    instead of being buried under them.
    """
    thin = _xn__thin(sub, group_col)
    sns.stripplot(data=thin, x=group_col, y=val_col, order=tiers, color=CFG.VIS_INK["dark"],
                  size=1.8, alpha=0.28, jitter=0.28, ax=ax, zorder=1, legend=False)
    sns.boxplot(data=sub, x=group_col, y=val_col, order=tiers, hue=group_col,
                palette=TIER_PALETTE, legend=False, fliersize=0, ax=ax,
                width=0.62, linewidth=1.0, zorder=3,
                boxprops=dict(alpha=0.85), medianprops=dict(color=CFG.VIS_INK["near_black"], linewidth=1.6))
    _xn__mean_trend(ax, sub, tiers, val_col, group_col=group_col)
    ax.set_xlabel('Degrader tier')
    ax.set_ylabel(ylabel)
    ax.set_axisbelow(True)
    _xn__tidy_tier_ticks(ax, counts=sub[group_col].value_counts())
    _xn__stat_header(ax, _xn__kruskal(sub, group_col, val_col, tiers))
    return ax

def _xn__variance_rows_for_job(args: tuple) -> tuple:
    """One complex: its 5 model CIFs parsed to rows. Runs in a worker process.

    Module-level and self-contained so it pickles for the process pool. Returns (rows, n_geom_fail)
    rather than logging: a worker writing to the run log would interleave its output with the others.
    """
    job_dir, _nuc = args
    _mre = re.compile(r'_(model_\d+)\.json$')
    _rows, _fail = [], 0
    for cj in sorted(job_dir.glob('boltz_results_*/predictions/*/confidence_*_model_*.json')):
        _mm = _mre.search(cj.name)
        if not _mm:
            continue
        cif_path = cj.with_name(cj.name.replace('confidence_', '', 1).replace('.json', '.cif'))
        if not cif_path.exists():
            continue
        try:
            with open(cj) as fh:
                cd = json.load(fh)
        except (OSError, ValueError):
            cd = {}
        try:
            geo = _xn_compute_sn2_geometry(cif_path, _nuc)
            sn2_dist, sn2_ang = float(geo.sn2_distance_A), float(geo.sn2_angle_deg)
        except Exception:                                        # noqa: BLE001
            # A model with no nucleophile in range has no SN2 geometry; that is a real outcome, not an
            # error. It is COUNTED and reported by the caller rather than silently becoming a NaN.
            sn2_dist, sn2_ang, _fail = float('nan'), float('nan'), _fail + 1
        _rows.append({'complex_id': job_dir.name, 'model_name': _mm.group(1),
                      'ptm': cd.get('ptm', np.nan), 'iptm': cd.get('iptm', np.nan),
                      'ligand_iptm': cd.get('ligand_iptm', np.nan),
                      'confidence_score': cd.get('confidence_score', np.nan),
                      'sn2_distance_A': sn2_dist, 'sn2_angle_deg': sn2_ang})
    return _rows, _fail


def _xn__ensure_multimodel_variance_csv(prod_dir: Path, out_dir: Path, reporter) -> Path:
    """Return the per-model variance CSV, computing it if absent.

    It READS the model CIFs from the production folder (02's output) and WRITES the CSV into 03's own
    01_Analysis_Data. A step writes its artefacts into its own folder: 1_Boltz2_Production is 02's
    record of what it produced, and 03 dropping a derived table into it would leave an analysis
    product filed as though the production run had made it.

    The geometry is recomputed by parsing every model .cif under 4_Prediction_Jobs. The SN2 geometry
    (nucleophile → reactive carbon → leaving fluorine) is delegated to the internal gemmi engine
    (_xn_compute_sn2_geometry) so the recomputed values are numerically identical to the pipeline's
    own analysis rather than a divergent re-implementation. Per-model ptm/iptm/ligand_iptm/
    confidence_score are read from the sibling confidence_*.json.
    """
    jobs_dir = prod_dir / '4_Prediction_Jobs'
    target = _aux_dir(out_dir) / '07_Boltz2_MultiModel_QC_Variance.csv'
    if not jobs_dir.exists():
        reporter.log(f'  ! Multi-model variance: {jobs_dir} not found; cannot compute.')
        return target
    job_folders = sorted((p for p in jobs_dir.iterdir() if p.is_dir()))
    n_jobs = len(job_folders)
    _all_ids = {p.name for p in job_folders}

    """
    An existing CSV is REUSED, but only after it is checked against the job folders it claims to
    describe. `exists()` alone is not a validation: a file truncated by a killed run, or one written
    before more predictions were added, would be trusted for ever and the uncertainty panels drawn
    from a corpus that is not the corpus.

    Three things are checked — the schema, the coverage, and whether any geometry actually resolved
    (an all-NaN table is the failure mode this build had, and it looks complete from the outside).
    Whatever it already covers is KEPT: only the complexes missing from it are parsed, so a build
    interrupted at 40,000 of 58,056 resumes from 40,000 rather than starting again.
    """
    _done_ids: set = set()
    _prior = None
    if target.exists():
        try:
            _prior = pd.read_csv(target, low_memory=False)
        except Exception as _e:                                  # noqa: BLE001
            reporter.log(f'  ! Existing variance CSV unreadable ({type(_e).__name__}); rebuilding.')
            _prior = None

    if _prior is not None:
        _needed = {'complex_id', 'model_name', 'sn2_distance_A', 'sn2_angle_deg'}
        if not _needed.issubset(_prior.columns):
            reporter.log('  ! Existing variance CSV has the wrong schema; rebuilding.')
            _prior = None       # rejected: it must not be carried into the rebuilt table
        elif int(_prior['sn2_distance_A'].notna().sum()) == 0:
            reporter.log('  ! Existing variance CSV has NO resolved geometry (every row NaN); rebuilding.')
            _prior = None       # rejected: concatenating it back would restore the rows just rejected
        else:
            # A complex counts as done only when all 5 of its models are present; a half-written
            # complex from a killed run is re-parsed rather than half-trusted.
            _per = _prior.groupby('complex_id')['model_name'].nunique()
            _done_ids = set(_per[_per >= 5].index) & _all_ids
            _prior = _prior[_prior['complex_id'].isin(_done_ids)]
            if _done_ids >= _all_ids:
                _geo = _prior['sn2_distance_A'].notna().mean()
                reporter.log(f'  ✔ Found: 01_Analysis_Data/{target.name}  '
                             f'(complete for all {n_jobs:,} complexes · geometry {_geo:.1%}) — reusing, nothing recomputed.')
                return target
            reporter.log(f'  ⧗ Variance CSV covers {len(_done_ids):,}/{n_jobs:,} complexes — '
                         f'resuming; only the {n_jobs - len(_done_ids):,} missing will be parsed.')

    job_folders = [p for p in job_folders if p.name not in _done_ids]
    n_jobs = len(job_folders)

    """
    Each complex is measured at ITS OWN catalytic aspartate — the residue 02 mapped for that protein by
    global alignment. Without it the geometry engine scans the whole protein and can settle on a residue
    that does no catalysis (on Tier_4 it chose a SERINE in 4 of 30), reporting an attack angle that is
    geometrically real and mechanistically meaningless.
    """
    _nuc_by_job: dict = {}
    _prot_by_job: dict = {}
    _rank_csv = (latest_by_mtime(prod_dir.glob(CFG.GLOB_RANKED_CSV))
                 or latest_by_mtime(prod_dir.glob("*_Ranked_*.csv")))
    if _rank_csv is not None:
        try:
            _rk = pd.read_csv(_rank_csv, low_memory=False,
                              usecols=["job_name", "Mapped_Nucleophile", "Protein_Name"])
            for _jn, _mn, _pn in zip(_rk["job_name"], _rk["Mapped_Nucleophile"], _rk["Protein_Name"]):
                _prot_by_job[str(_jn)] = str(_pn)
                _m = re.search(r"(\d+)\s*$", str(_mn))
                if _m:
                    _nuc_by_job[str(_jn)] = int(_m.group(1))
        except Exception as _e:                                  # noqa: BLE001
            # There is NO whole-protein fallback — _xn_compute_sn2_geometry refuses without the
            # mapped aspartate, because the fallback scan is what put a serine in the results.
            # Complexes with no mapped nucleophile are SKIPPED, and the message must say so.
            reporter.log(f'  ! Could not read Mapped_Nucleophile ({type(_e).__name__}); '
                         f'those complexes will be SKIPPED — geometry is measured only at the '
                         f'mapped catalytic aspartate, never by a whole-protein scan.')
    _missing_jobs = [p.name for p in job_folders if p.name not in _nuc_by_job]
    _missing = len(_missing_jobs)
    if _missing:
        # The gap is per-protein, not per-complex: when the alignment cannot place a protein's catalytic
        # aspartate, ALL of that protein's ligand complexes lack it. State both counts so the number here
        # reconciles with the per-variant coverage figure (81 complexes = 3 protein variants × 27 ligands),
        # rather than reading as a contradiction of it.
        _missing_prots = {_prot_by_job.get(j) for j in _missing_jobs if _prot_by_job.get(j)}
        _npr = len(_missing_prots) or "?"
        reporter.log(f'  ! {_missing:,} complex(es) have no mapped nucleophile — these are {_npr} protein '
                     f'variant(s) (each × its ligands) whose catalytic aspartate the alignment could not '
                     f'place; see 01_Active_Site_Residue_Mapping_Coverage. SKIPPED — geometry is measured '
                     f'only at the mapped catalytic aspartate, never by a whole-protein scan.')
    _tasks = [(p, _nuc_by_job.get(p.name)) for p in job_folders]

    """
    The work is one independent gemmi parse per model CIF — five per complex, tens of thousands of
    complexes — with no shared state, so it parallelises cleanly across cores. Two cores are left free
    so the machine stays usable and a GPU job's feeder process is never starved.
    """
    _n_proc = max(1, (_os.cpu_count() or 4) - 2)
    reporter.log(f'  ⧗ Building 07_Boltz2_MultiModel_QC_Variance.csv — parsing {n_jobs:,} complexes '
                 f'× 5 model CIFs across {_n_proc} cores…')

    """
    A three-hour job that prints one line every few thousand complexes is indistinguishable from a
    hung one. The console gets a live counter rewritten in place (\\r) carrying rate and ETA, so the
    run can be watched; the run LOG gets a plain milestone line every 5,000 complexes, because a log
    file full of carriage returns is unreadable. The counter is only drawn to a terminal — piped to a
    file it would be noise.
    """
    _t0 = _time.time()
    _tty = sys.stdout.isatty()

    def _tick(done: int) -> None:
        _el = max(1e-6, _time.time() - _t0)
        _rate = done / _el
        _eta = (n_jobs - done) / _rate if _rate > 0 else 0.0
        if _tty:
            _bar_n = 28
            _fill = int(_bar_n * done / max(1, n_jobs))
            sys.stdout.write(
                f'\r    · variance: [{"█" * _fill}{"·" * (_bar_n - _fill)}] '
                f'{done:,}/{n_jobs:,} ({done / max(1, n_jobs):5.1%})  '
                f'{_rate:6.1f} cx/s  elapsed {_el / 60:5.1f}m  ETA {_eta / 60:5.1f}m   ')
            sys.stdout.flush()
        if done % 5000 == 0 or done == n_jobs:
            _milestone = (f'    · variance progress: {done:,}/{n_jobs:,} complexes '
                          f'({done / max(1, n_jobs):.0%})  ·  elapsed {_el / 60:.1f} min  ·  ETA {_eta / 60:.1f} min')
            _file_only = getattr(reporter, 'log_file_only', None)
            if _tty and _file_only is not None:
                # The console already carries this in the bar; printing it again would break the bar.
                _file_only(_milestone)
            else:
                reporter.log(_milestone)

    target.parent.mkdir(parents=True, exist_ok=True)

    def _flush(rows_so_far: list) -> None:
        """Write everything resolved so far — the prior rows plus the new ones — atomically.

        A build this long must be able to die without losing the hours it already paid for. The write
        goes to a temporary file and is renamed over the target, because rename is atomic within a
        filesystem: a run killed mid-write leaves either the old complete file or the new one, never a
        truncated hybrid that the next run would read as truth and resume from.
        """
        _out = pd.DataFrame(rows_so_far)
        if _prior is not None and len(_prior):
            _out = pd.concat([_prior, _out], ignore_index=True)
        _tmp = target.with_suffix('.csv.tmp')
        _out.to_csv(_tmp, index=False)
        _tmp.replace(target)

    rows: list = []
    _n_geom_fail = 0
    with cf.ProcessPoolExecutor(max_workers=_n_proc) as _ex:
        _done = 0
        for _job_rows, _fails in _ex.map(_xn__variance_rows_for_job, _tasks, chunksize=16):
            rows.extend(_job_rows)
            _n_geom_fail += _fails
            _done += 1
            if _done % 50 == 0 or _done == n_jobs:
                _tick(_done)
            if _done % 5000 == 0:
                _flush(rows)

    # Close the bar's line, so the summary below is not written over the last frame of it.
    if _tty:
        sys.stdout.write('\n')
        sys.stdout.flush()

    _flush(rows)
    var_df = pd.read_csv(target, low_memory=False)

    """
    The geometry failures are REPORTED, not swallowed. A silent `except: nan` here is what let the
    whole table be written with an empty geometry column: every row present, every distance and angle
    NaN, and nothing in the log to say so. If the geometry is missing for most models the CSV is not
    usable and the run must say that out loud.
    """
    _n_rows = len(var_df)
    _n_geo = int(var_df['sn2_distance_A'].notna().sum()) if _n_rows else 0
    reporter.log(f"  ✔ Saved: 01_Analysis_Data/{target.name}  "
                 f"({var_df['complex_id'].nunique():,} complexes · {_n_rows:,} model rows · "
                 f"geometry resolved {_n_geo:,} = {_n_geo / max(1, _n_rows):.1%})")
    if _n_geo == 0:
        reporter.log('  ! Multi-model variance: NO geometry resolved on any model — the uncertainty '
                     'panels would be empty. Treat this CSV as unusable.')
    elif _n_geom_fail:
        reporter.log(f'    · {_n_geom_fail:,} of {_n_rows:,} model poses ({_n_geom_fail / max(1, _n_rows):.1%}) '
                     f'had no resolvable SN2 geometry — the mapped aspartate was absent or beyond attack '
                     f'range in that model (this includes the no-nucleophile complexes above, × 5 models).')
    return target

def _xn__fig_01C_geometry_and_uncertainty(df, out_dir, reporter):
    """Geometry and inter-model uncertainty by degrader tier (2x2 grid).

    Identical to Figure_01c_Geometry_and_Uncertainty, except the 180 degree
    "ideal in-line" reference line (and its legend) is removed from the
    SN2-angle panel (top-right).

    Uses the variance CSV's own ``degrader_tier`` column (identical across a
    complex's 5 models — it is the pipeline's authoritative classification)
    rather than re-deriving a tier from raw distance/angle alone. A purely
    geometric re-derivation cannot represent Tier_5_Decoy here: the 02b
    geometry engine only searches for a nucleophile within a 6 A radius
    (``geometry_status == "no_nucleophile_within_6A"`` otherwise), so
    ``sn2_distance_A`` never exceeds ~6 A in practice — well inside the
    Tier_4 cutoff (8 A). Combined with taking the best of 5 models per
    complex, that made Tier_5_Decoy essentially unreachable even though it
    is a genuine, populated tier in the source data.
    """
    """
    The production directory sits beside 3_Validation_Figures under the RUN root — not beside the
    figure subfolder that out_dir points at. It is taken from the module global the dispatcher sets;
    the fallback walks up from out_dir rather than assuming a fixed depth, because deriving it as
    out_dir.parent silently yields 3_Validation_Figures/1_Boltz2_Production — a path that cannot
    exist, and the variance CSV is then reported 'unavailable' while the CIFs sit there all along.
    """
    _prod = globals().get("_xn__PROD_DIR") or globals().get("_PROD_DIR")
    if _prod is None or not Path(_prod).exists():
        _prod = next((_c for _c in (out_dir.parent.parent / '1_Boltz2_Production',
                                    out_dir.parent / '1_Boltz2_Production')
                      if _c.exists()), out_dir.parent / '1_Boltz2_Production')
    """
    The inter-model uncertainty panels need the per-model variance CSV. Building it from scratch
    means parsing 5 model CIFs for each of ~58,000 complexes — hours of gemmi. That is a
    deliberate, opt-in job (--variance), not something a figure refresh should trigger silently:
    without the flag, an absent CSV falls through to the two absolute-geometry panels, which
    carry the same tier claim and are drawn from the ranked CSV in seconds.
    """
    # The variance CSV is built up front in Step 1 (01_Analysis_Data); here we only READ the cache.
    _fig_root = out_dir if out_dir.name.startswith('3_') else out_dir.parent
    var_path = next((p for p in (_aux_dir(_fig_root) / '07_Boltz2_MultiModel_QC_Variance.csv',
                                 _prod / '4_Prediction_Jobs' / '07_Boltz2_MultiModel_QC_Variance.csv')
                     if p.exists()), None)
    if var_path is None and not globals().get('_ALLOW_VARIANCE_COMPUTE', False):
        reporter.log('  · no per-model variance CSV; drawing the two geometry panels. '
                     'Pass --variance to build it from the CIFs (slow) and get the uncertainty panels too.')
    if var_path is None or not Path(var_path).exists():
        # FAcDs has no per-model variance CSV (needs the 02b reanalysis engine), so the two
        # inter-model uncertainty panels cannot be drawn. Plot the two absolute-geometry panels
        # (nucleophile distance, SN2 attack angle) by tier from the ranked CSV instead.
        reporter.log('  ! Figure 01C: variance CSV unavailable — plotting geometry-only (nucleophile distance + SN2 angle) from ranked CSV.')
        _dist = _xn__nuc_distance(df).where(lambda s: s < 20.0)  # drop ~999/1000 Å "no nucleophile" sentinel
        _ang = _xn__num(df, _xn__col(df, CFG.COL_SN2))
        _tcol = _xn__col(df, CFG.COL_TIER)
        if _tcol is None or _dist.empty or _ang.empty:
            reporter.log(f"  ! Skipped: {_fig_path('01C')} — geometry columns unavailable.")
            return
        gdf = pd.DataFrame({'tier': df[_tcol].values, 'dist': _dist.values, 'ang': _ang.values}).dropna(subset=['dist', 'ang'])
        if gdf.empty:
            reporter.log(f"  ! Skipped: {_fig_path('01C')} — no geometry rows.")
            return
        _tiers = [t for t in TIER_ORDER_LOGIC if t in set(gdf['tier'])]
        fig, (axd, axa) = plt.subplots(1, 2, figsize=(_xn_COL_DOUBLE_IN, 0.45 * _xn_COL_DOUBLE_IN))
        _xn__tier_boxstrip(axd, gdf, _tiers, 'dist', 'Nucleophile distance  (Å)')
        _xn__tier_boxstrip(axa, gdf, _tiers, 'ang', 'SN2 attack angle  (°)')
        _xn__gate_lines_distance(axd)
        _xn__gate_lines_angle(axa)
        # The nucleophile distance is long-tailed: Tier_5_Decoy reaches ~17 Å while every
        # catalytically meaningful difference sits between 2 and 6 Å. On a linear axis spanning
        # the tail, the tiers the paper is about collapse into one flat line.
        axd.set_ylim(0, float(np.nanpercentile(gdf['dist'], 99)) * 1.05)
        axa.set_ylim(0, 185)
        axa.set_yticks([0, 45, 90, 135, 180])
        _xn__panel(axd, 'a')
        _xn__panel(axa, 'b')
        """
        ONE tier legend for BOTH panels, inside panel A. The two panels share the same colour
        encoding, so repeating the key beside each of them would state the same thing twice and eat
        the space the data needs. Placed inside the axes rather than beside the figure, it also
        cannot be cropped away from the panels it explains.
        """
        """
        ONE legend, top-left of panel A: the tier colours AND the gate tags.

        The tags mean the same tier in both panels (A = Tier_1A, B = Tier_2A, C = Tier_2B) and differ
        only in the threshold each panel enforces, so a single entry can carry both — the distance
        gate and the angle gate side by side. Two separate boxes stated the tier twice and spent the
        panel's space saying it.
        """
        from matplotlib.patches import Patch as _P1C
        from matplotlib.lines import Line2D as _L1C
        _h = [_P1C(facecolor=TIER_PALETTE.get(_t, CFG.VIS_INK["faint"]), alpha=0.85,
                   edgecolor=CFG.VIS_INK["outline"], linewidth=0.7, label=_t)
              for _t in _tiers]
        _l = list(_tiers)
        _dk = {k.split(' = ')[0]: k.split(' = ')[1] for k in getattr(axd, '_gate_keys', [])}
        _ak = {k.split(' = ')[0]: k.split(' = ')[1] for k in getattr(axa, '_gate_keys', [])}
        for _tag in ('A', 'B', 'C'):
            if _tag not in _dk:
                continue
            _tier_name = _dk[_tag].split(' ≤')[0]
            _dist_gate = _dk[_tag].split(' ≤')[-1]
            _ang_gate = _ak.get(_tag, '').split(' ≥')[-1] if _tag in _ak else ''
            _h.append(_L1C([], [], ls=':', lw=1.4,
                           color=TIER_PALETTE.get(_tier_name, CFG.VIS_INK["grey"])))
            _l.append(f'({_tag}) gate: ≤{_dist_gate}  ·  ≥{_ang_gate}')
        axd.legend(_h, _l, loc='upper left', bbox_to_anchor=(0.012, 0.985), ncol=2,
                     fancybox=True,

                   title='Degrader tier   ·   gate tags (A/B/C on the rules)',
                   )
        _xn__save(fig, out_dir, _xn_FIG_NAMES['geometry'], reporter)
        return
    try:
        vdf = pd.read_csv(var_path)
    except FileNotFoundError:
        reporter.log(f"  ! Skipped: {_fig_path('01C')} — could not open {var_path}.")
        return
    id_col = next((c for c in ('complex_id', 'job_name') if c in vdf.columns), None)
    if id_col is None or 'sn2_angle_deg' not in vdf.columns or 'sn2_distance_A' not in vdf.columns:
        reporter.log(f"  ! Skipped: {_fig_path('01C')} — variance CSV lacks id / 'sn2_distance_A' / 'sn2_angle_deg'.")
        return
    """
    The tier is a property of the COMPLEX, not of a diffusion sample, so it belongs to the ranked sheet
    and not to a per-model table. It is joined in here on the complex id rather than duplicated into
    every one of the 290,000 model rows — and joining also means an existing variance CSV stays valid
    when the tiering changes, instead of 165 minutes of parsing being thrown away because a column it
    never needed to carry is missing.
    """
    if CFG.COL_TIER not in vdf.columns:
        _rk = (latest_by_mtime(Path(_prod).glob(CFG.GLOB_RANKED_CSV))
               or latest_by_mtime(Path(_prod).glob("*_Ranked_*.csv")))
        if _rk is not None:
            try:
                _tier_map = pd.read_csv(_rk, low_memory=False,
                                        usecols=["job_name", CFG.COL_TIER])
                vdf = vdf.merge(_tier_map, how="left", left_on=id_col, right_on="job_name")
                _n_tier = int(vdf[CFG.COL_TIER].notna().sum())
                reporter.log(f'  · tier joined from the ranked sheet for {_n_tier:,} of '
                             f'{len(vdf):,} model rows.')
            except Exception as _e:                              # noqa: BLE001
                reporter.log(f"  ! Figure 01C: could not join the tier ({type(_e).__name__}: {_e}).")
    if CFG.COL_TIER not in vdf.columns or vdf[CFG.COL_TIER].notna().sum() == 0:
        reporter.log(f"  ! Skipped: {_fig_path('01C')} — no degrader_tier available (not in the variance "
                     "CSV, and the ranked sheet could not supply it).")
        return
    vdf['sn2_distance_A'] = pd.to_numeric(vdf['sn2_distance_A'], errors='coerce')
    vdf['sn2_angle_deg'] = pd.to_numeric(vdf['sn2_angle_deg'], errors='coerce')
    per_complex = []
    for cid, g in vdf.groupby(id_col):
        g_valid = g.dropna(subset=['sn2_distance_A', 'sn2_angle_deg'])
        if g_valid.empty:
            continue
        best_row = g_valid.loc[g_valid['sn2_distance_A'].idxmin()]
        per_complex.append({'complex_id': cid, 'best_geo_tier': g[CFG.COL_TIER].iloc[0], 'abs_distance': float(best_row['sn2_distance_A']), 'abs_angle': float(best_row['sn2_angle_deg']), 'distance_std': float(g['sn2_distance_A'].std(ddof=1)), 'angle_std': float(g['sn2_angle_deg'].std(ddof=1))})
    cdf = pd.DataFrame(per_complex)
    if cdf.empty:
        reporter.log(f"  ! Skipped: {_fig_path('01C')} — no complexes after geometric tiering.")
        return
    # Top row plots the SELECTED (representative) model the pipeline finally reports, taken from the
    # ranked sheet (df) per complex — not a per-model distance minimum. Falls back to the closest-of-5
    # only when the ranked geometry columns are unavailable.
    _jc = _xn__col(df, "job_name") or "job_name"
    _dc = _xn__col(df, "Dist_Nucleophile") or "Dist_Nucleophile"
    _ac = _xn__col(df, CFG.COL_SN2) or CFG.COL_SN2
    _sel_ok = all(c in df.columns for c in (_jc, _dc, _ac))
    if _sel_ok:
        _seldf = (df[[_jc, _dc, _ac]]
                  .rename(columns={_jc: "complex_id", _dc: "sel_distance", _ac: "sel_angle"})
                  .drop_duplicates("complex_id"))
        _seldf["sel_distance"] = pd.to_numeric(_seldf["sel_distance"], errors="coerce")
        _seldf["sel_angle"] = pd.to_numeric(_seldf["sel_angle"], errors="coerce")
        cdf = cdf.merge(_seldf, on="complex_id", how="left")
    _use_sel = _sel_ok and "sel_distance" in cdf.columns and cdf["sel_distance"].notna().any()
    _top_d, _top_a = ("sel_distance", "sel_angle") if _use_sel else ("abs_distance", "abs_angle")
    _sfx = "Selected Model" if _use_sel else "closest of 5 models (ranked geometry unavailable)"
    tiers = [t for t in TIER_ORDER_LOGIC if t in set(cdf['best_geo_tier'])]
    # Tight inter-row spacing: the top row carries no x tick labels, so the default gap leaves a wide
    # empty band between the rows. A small hspace pulls the two rows together, matching the other figures.
    fig, axes = plt.subplots(2, 2, figsize=(_xn_COL_DOUBLE_IN, 0.82 * _xn_COL_DOUBLE_IN),
                             gridspec_kw={'hspace': 0.06})
    (ax_tl, ax_tr), (ax_bl, ax_br) = axes
    _xn__tier_boxstrip(ax_tl, cdf, tiers, _top_d, f'Nucleophile distance  (Å)\n{_sfx}', group_col='best_geo_tier')
    _xn__gate_lines_distance(ax_tl)
    ax_tl.set_ylim(1.5, max(5.5, float(np.nanmax(cdf[_top_d])) * 1.02))   # 1.5 Å floor (no data below) → ~5.5 Å so the full scatter shows
    _xn__tier_boxstrip(ax_tr, cdf, tiers, _top_a, f'SN2 attack angle  (°)\n{_sfx}', group_col='best_geo_tier')
    _xn__gate_lines_angle(ax_tr)
    ax_tr.set_ylim(0, 185)
    ax_tr.set_yticks([0, 45, 90, 135, 180])
    # The top and bottom rows share the same tier x-axis; the bottom row already carries the tier names
    # and their n, so the top row's tick labels only repeat them — suppress them (labels shown once, below).
    for _axt in (ax_tl, ax_tr):
        _axt.tick_params(labelbottom=False)
        _axt.set_xlabel('')
    for _ax, _y, _lab in [(ax_bl, 'distance_std', 'Distance s.d. across 5 models  (Å)'),
                          (ax_br, 'angle_std', 'Angle s.d. across 5 models  (°)')]:
        sns.violinplot(data=cdf, x='best_geo_tier', y=_y, order=tiers, hue='best_geo_tier',
                       palette=TIER_PALETTE, legend=False, cut=0, inner='quartile',
                       linewidth=0.9, ax=_ax)
        for _c in _ax.collections:
            _c.set_alpha(0.85)
        # The violin's inner quartiles show the SPREAD of the spread; the trend line answers the question
        # the panel is actually asked — does inter-model disagreement rise or fall with tier?
        _xn__mean_trend(_ax, cdf, tiers, _y, group_col='best_geo_tier')
        _ax.set_xlabel('Degrader tier')
        _ax.set_ylabel(_lab)
        _ax.set_axisbelow(True)
        _xn__tidy_tier_ticks(_ax, counts=cdf['best_geo_tier'].value_counts())
        _xn__stat_header(_ax, _xn__kruskal(cdf, 'best_geo_tier', _y, tiers))
    # The tier names and their n already sit on the bottom tick labels; the figure-level
    # 'Degrader tier (n below each tier)' caption only restates that, so clear every panel's x-label
    # (no shared x-label left ⇒ _ext_match_03_style lifts no supxlabel).
    for _ax in (ax_bl, ax_br):
        _ax.set_xlabel('')
    for _ax in (ax_tl, ax_tr, ax_bl, ax_br):
        _tier_seps(_ax, len(tiers))
    for _ax, _l in [(ax_tl, 'a'), (ax_tr, 'b'), (ax_bl, 'c'), (ax_br, 'd')]:
        _xn__panel(_ax, _l)
    _xn__save(fig, out_dir, _xn_FIG_NAMES['geometry'], reporter)


def _xn__fig_05a_pillar_divergence_modified(df, out_dir, reporter):
    # Tier-resolved pillar divergence: one violin panel per scoring pillar, each showing the
    # pillar-score distribution across degrader tiers (green→red gradient). Replaces the flat
    # size-only line plot so each tier's contribution is visible.
    pillars = [('Model_Quality_Score', 'Model Quality'),
               ('Binding_Affinity_Score', 'Binding Affinity'),
               ('Catalytic_Competence_Score', 'Catalytic Competence'),
               ('Evolutionary_Fingerprint_Score', 'Evolutionary Fingerprint')]
    tcol = _xn__col(df, CFG.COL_TIER)
    if tcol is None:
        return
    tiers = [t for t in TIER_ORDER_LOGIC if t in set(df[tcol].dropna())]
    if not tiers:
        return
    fig, axes = plt.subplots(2, 2, figsize=(_xn_COL_DOUBLE_IN, 0.85 * _xn_COL_DOUBLE_IN),
                             gridspec_kw={'hspace': 0.06})
    counts = df[tcol].value_counts()
    for ax, letter, (pil, lab) in zip(axes.ravel(), ['a', 'b', 'c', 'd'], pillars):
        pc = _xn__pillar_col(df, pil)
        if pc is None or pc not in df.columns:
            ax.set_visible(False)
            continue
        sub = pd.DataFrame({'tier': df[tcol].values, 'y': pd.to_numeric(df[pc], errors='coerce').values}).dropna()
        sub = sub[sub['tier'].isin(tiers)]
        if sub.empty:
            ax.set_visible(False)
            continue
        sns.violinplot(data=sub, x='tier', y='y', order=tiers, hue='tier', palette=TIER_PALETTE,
                       legend=False, cut=0, inner='box', linewidth=0.9, ax=ax, zorder=2)
        for _c in ax.collections:
            _c.set_alpha(0.75)
        # The median trend across tiers. Four violins side by side show four distributions; the
        # question the panel is asked is whether the pillar rises or falls with tier, and only a
        # connected median answers that without the reader eyeballing four white bars.
        med = sub.groupby('tier', observed=True)['y'].median().reindex(tiers)
        ax.plot(range(len(tiers)), med.values, '-', color=CFG.VIS_INK["ink_pure"], lw=1.2, alpha=0.65,
                marker='o', markersize=3.2, markerfacecolor='white', markeredgewidth=0.9,
                zorder=6, label='Median trend')
        ax.set_xlabel('Degrader tier')
        ax.set_ylabel(lab)
        ax.set_axisbelow(True)
        _xn__tidy_tier_ticks(ax, counts=counts)
        _xn__stat_header(ax, _xn__kruskal(sub, 'tier', 'y', tiers, reg_name=pc))
        _xn__panel(ax, letter)
    # Both rows share the tier x-axis: the top row's tick labels only repeat the bottom row's, and the
    # tier names + n already sit on the bottom ticks — so suppress the top row's labels and clear every
    # x-label (no shared x-label left ⇒ no redundant 'Degrader tier (n below each tier)' supxlabel).
    for _ax in axes[0]:
        _ax.tick_params(labelbottom=False)
    for _ax in axes.ravel():
        _ax.set_xlabel('')
        _tier_seps(_ax)   # consistent tier separators
    _xn__save(fig, out_dir, _xn_FIG_NAMES['pillars'], reporter)

def _xn_figure_06a(df: pd.DataFrame, out_dir: Path, reporter) -> None:
    """Create Tier_1A enzyme × ligand heatmap for lab validation targets."""
    import matplotlib.colors as mcolors
    import re
    enzyme_col = CFG.COL_PROT if CFG.COL_PROT in df.columns else next((c for c in df.columns if c.lower() == 'job_name'), None)
    ligand_col = CFG.COL_LIG if CFG.COL_LIG in df.columns else next((c for c in df.columns if c.lower() in ('ligand', 'ligand_name')), None)
    tier_col = CFG.COL_TIER if CFG.COL_TIER in df.columns else None
    if not enzyme_col or not ligand_col or (not tier_col):
        if reporter:
            reporter.log(f"  ! Skipped: {_fig_path('06A')} — required columns missing")
        return
    df = df[df[ligand_col] != 'Fluoroacetate_Ref']
    df_hm = df[[enzyme_col, ligand_col, tier_col]].copy()
    df_hm[enzyme_col] = df_hm[enzyme_col].astype(str).str.strip()
    df_hm[ligand_col] = df_hm[ligand_col].astype(str).str.strip()
    df_hm[tier_col] = df_hm[tier_col].astype(str).fillna(CFG.TIER_DECOY)
    tier1a_enzymes = df_hm.loc[df_hm[tier_col] == CFG.TIER_TOP, enzyme_col].unique()
    if len(tier1a_enzymes) == 0:
        if reporter:
            reporter.log(f"  ! Skipped: {_fig_path('06A')} — no Tier_1A enzymes found")
        return
    tier_rank_map = CFG.TIER_RANK
    df_hm['tier_rank'] = df_hm[tier_col].map(tier_rank_map).fillna(max(tier_rank_map.values()) + 1).astype(int)
    df_best = df_hm.sort_values('tier_rank').drop_duplicates(subset=[enzyme_col, ligand_col], keep='first')

    def get_versatility_score(enz):
        sub = df_best[df_best[enzyme_col] == enz]
        score = 0
        for lig in ['MFA', 'DFA', 'TFA', 'Fluoroacetate']:
            tier_val = sub.loc[sub[ligand_col].str.lower() == lig.lower(), tier_col]
            if not tier_val.empty and tier_val.values[0] == CFG.TIER_TOP:
                score += 10
            elif not tier_val.empty and tier_val.values[0] == CFG.TIER_ORDER[1]:
                score += 5
        total_1a = (sub[tier_col] == CFG.TIER_TOP).sum()
        return (score, total_1a)
    enzymes = sorted(tier1a_enzymes, key=get_versatility_score, reverse=True)

    def _extract_number(l):
        m = re.search('\\d+', l)
        return int(m.group(0)) if m else float('inf')
    ligands = sorted(df_best[ligand_col].unique(), key=_extract_number)
    if not ligands:
        if reporter:
            reporter.log(f"  ! Skipped: {_fig_path('06A')} — no ligands found")
        return

    def _canonical_ligand_name(raw_name: object) -> str:
        # Delegate to the module-level shortener so this heatmap's labels match every other figure
        # (FA / DFA / TFA, with the "NN_" ordering prefix stripped).
        return _lig_short(raw_name)
    df_best['ligand_display'] = df_best[ligand_col].map(_canonical_ligand_name)
    ligands_display_ordered = []
    for l in ligands:
        dl = _canonical_ligand_name(l)
        if dl not in ligands_display_ordered:
            ligands_display_ordered.append(dl)
    pivot = pd.pivot_table(df_best, values='tier_rank', index=enzyme_col, columns='ligand_display', aggfunc='min', fill_value=np.nan).reindex(index=enzymes, columns=ligands_display_ordered)
    inv_rank_map = {v: '5' if k == 'Tier_5_Decoy' else k.replace('Tier_', '') for k, v in tier_rank_map.items()}
    pivot_labels = pivot.map(lambda x: inv_rank_map.get(int(x), '') if pd.notna(x) else '')
    fig_h = max(3.0, len(enzymes) * 0.25 + 1.8)
    fig_w = max(_xn_COL_DOUBLE_IN, len(ligands_display_ordered) * 0.42 + 2.4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    unique_ranks = sorted([r for r in tier_rank_map.values()])
    color_list = []
    for rank in unique_ranks:
        tier_name = [k for k, v in tier_rank_map.items() if v == rank][0]
        if tier_name == CFG.TIER_TOP:
            color_list.append(CFG.VIS_ACCENT["green"])
        else:
            color_list.append(TIER_PALETTE.get(tier_name, CFG.VIS_INK["palest"]))
    cmap = mcolors.ListedColormap(color_list)
    bounds = np.arange(min(unique_ranks) - 0.5, max(unique_ranks) + 1.5, 1)
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    # Cell labels are written in whichever of black/white survives the fill they sit on. A fixed
    # black label is legible on the pale mid-tiers and disappears into the dark Tier_1A green.
    _label_colours = pivot.map(
        lambda x: (_utils_mod.auto_label_colour(CFG, cmap(norm(int(x))))
                   if (pd.notna(x) and _utils_mod is not None) else CFG.VIS_INK["near_black"]))
    sns.heatmap(pivot, ax=ax, cmap=cmap, norm=norm, cbar=False, linewidths=0.5, linecolor='white',
                annot=pivot_labels, fmt='s', annot_kws={'fontsize': 7.5, 'weight': 'bold'},
                square=False, mask=pivot.isna())
    for _txt in ax.texts:
        _r, _c = int(_txt.get_position()[1] - 0.5), int(_txt.get_position()[0] - 0.5)
        try:
            _txt.set_color(_label_colours.iloc[_r, _c])
        except Exception:
            pass
    ax.set_facecolor(CFG.VIS_INK["grid"])
    ax.set_ylabel('Enzyme')
    ax.set_xlabel('PFAS ligand')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0, fontsize=CFG.VIS_FONT_TICK - 1.0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=CFG.VIS_FONT_TICK - 1.0)
    ax.tick_params(length=0)
    handles = []
    for rank, color in zip(unique_ranks, color_list):
        tier_name = [k for k, v in tier_rank_map.items() if v == rank][0]
        handles.append(plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=color,
                                  markersize=9, label=tier_name.replace('Tier_', '')))
    _grey = plt.Line2D([0], [0], marker='s', color='w', markerfacecolor=CFG.VIS_INK["grid"],
                       markeredgecolor=CFG.VIS_INK["paler"], markersize=9, label='not modelled')
    # ABOVE the heatmap. Below it, the legend has to share a strip with the rotated ligand labels and
    # the x-axis title, and every placement that clears one collides with the other. The top edge is
    # empty, so the legend sits there — attached to the panel, over nothing.
    ax.legend(handles=handles + [_grey], title='Degradation tier  (cell label = tier)',
              loc='lower center', bbox_to_anchor=(0.5, 1.015), ncol=len(handles) + 1,
              fontsize=CFG.VIS_FONT_LEGEND,  frameon=False)
    _xn__save(fig, out_dir, _xn_FIG_NAMES['validation'], reporter)



# The ELITE catalytic tiers (Tier_1A/1B), split OUT from the rest. This panel family contrasts elite vs
# rest — a different, deliberately narrower question than the battery's degrader/non-degrader split
# (CFG.TIER_HIGH_QUALITY, 4 tiers) — and is labelled "Elite vs rest" so the two never collide under one
# name. Keeping it at the elite pair preserves the affinity≠catalysis reading: the elite catalytic tiers
# (the FA/DFA-like near-attack poses) are LOW binding affinity, the high-affinity binders are not elite.
_xo_ELITE_TIERS = ["Tier_1A", "Tier_1B"]


def _xo__minmax(s):
    s = pd.to_numeric(s, errors='coerce')
    lo, hi = (np.nanmin(s), np.nanmax(s))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return pd.Series(np.nan, index=s.index)
    return (s - lo) / (hi - lo)

def _xo__nuc_distance(df: pd.DataFrame) -> pd.Series:
    """Nucleophile–substrate distance, resolving the several historical names.

    The Step-02 writer renames ``dist_Nuc`` to ``Dist_Nucleophile_ASP10`` via
    ``CFG.COLUMN_RENAMING_MAP``; older/raw exports expose the live geometry as
    ``best_nucleophile_distance``. All three are treated as equivalent.
    """
    col = _xo__col(df, 'Dist_Nucleophile_ASP10', 'Dist_Nucleophile', 'best_nucleophile_distance', 'dist_Nuc')
    return _xo__num(df, col)

def _xo_get_col(df, name, fallbacks=[]):
    if name in df.columns:
        return name
    for f in fallbacks:
        if f in df.columns:
            return f
    return None

_xo_PILLAR_ALIASES = {
    'Model_Quality_Score': ['Model_Quality_Score', 'Model_Quality_ScoreNormalised', CFG.COL_CONF],
    'Binding_Affinity_Score': ['Binding_Affinity_Score', 'Binding_Affinity_ScoreNormalised', 'Chemical_Affinity_Score', 'Binding_Probability', 'custom_affinity_score', 'Binding_Probability_Score'],
    'Catalytic_Competence_Score': ['Catalytic_Competence_Score', 'Catalytic_Competence_ScoreNormalised', 'competence_score', 'soft_catalytic_score'],
    'Evolutionary_Fingerprint_Score': ['Evolutionary_Fingerprint_Score', 'Evolutionary_Fingerprint_ScoreNormalised', CFG.COL_LIKE_S, CFG.COL_ID_PCT],
    'Final_Unified_Score': ['Final_Unified_Score', 'Ranking_Score_Calc'],
}

def _xo__col(df: pd.DataFrame, *candidates: str) -> str | None:
    """Return the first candidate column that exists in ``df`` (else ``None``)."""
    for c in candidates:
        if c in df.columns:
            return c
    return None

def _xo__pillar_col(df: pd.DataFrame, pillar: str) -> str | None:
    """Resolve a headline-pillar name to whichever concrete column is present."""
    return _xo__col(df, *_xo_PILLAR_ALIASES.get(pillar, [pillar]))

def _xo__num(df: pd.DataFrame, col: str | None) -> pd.Series:
    """Coerce a column to numeric, returning an empty series when it is absent."""
    if col is None or col not in df.columns:
        return pd.Series(dtype=float)
    return pd.to_numeric(df[col], errors='coerce')

def _xo__tiers_present(df: pd.DataFrame, tier_col: str=CFG.COL_TIER) -> list[str]:
    """Return the canonical tier ordering restricted to tiers actually present."""
    if tier_col not in df.columns:
        return []
    have = set(df[tier_col].dropna().unique())
    return [t for t in TIER_ORDER_LOGIC if t in have]

def _xo__fmt_p(p: float) -> str:
    """Compact, publication-style p-value formatting."""
    if p is None or not np.isfinite(p):
        return 'p = n/a'
    if p < 0.0001:
        return 'p < 1e-4'
    return f'p = {p:.3g}'

def _xo__annotate(ax, text: str, loc: str='upper right') -> None:
    """Place a boxed statistics annotation on an axis."""
    xy = {'upper right': (0.98, 0.97, 'right', 'top'), 'upper left': (0.02, 0.97, 'left', 'top'), 'lower right': (0.98, 0.03, 'right', 'bottom'), 'lower left': (0.02, 0.03, 'left', 'bottom')}.get(loc, (0.98, 0.97, 'right', 'top'))
    ax.text(xy[0], xy[1], text, transform=ax.transAxes, ha=xy[2], va=xy[3], fontsize=CFG.VIS_FONT_ANNOT,
            bbox=dict(boxstyle='round,pad=0.4', fc='white', ec=CFG.VIS_INK["paler"], alpha=CFG.VIS_LEGEND_FRAME_ALPHA))

def _xo__legend_with_stats(ax, handles, labels, stat_lines, loc, ncol=1):
    """One combined box: the legend entries, then the statistics lines as blank-handle rows,
    so the legend and the stats annotation read as a single unit rather than two boxes.

    ncol lays the entries out in columns. Stacked in a single column they grow into a tall strip
    down the side of the panel and start covering the data they describe; two or three columns give
    the same information in a fraction of the height. Font, size, box transparency and spacing all
    inherit apply_figure_style's rcParams, so this legend matches every other legend in the set."""
    from matplotlib.lines import Line2D as _L2D
    _blank = lambda: _L2D([], [], linestyle='', marker='', color='none')
    h = list(handles) + [_blank()] + [_blank() for _ in stat_lines]
    l = list(labels) + [''] + list(stat_lines)
    # A tuple loc is an axes-fraction anchor; matplotlib takes it via bbox_to_anchor, not loc.
    _kw = dict(ncol=max(1, int(ncol)))
    if isinstance(loc, (tuple, list)):
        ax.legend(h, l, loc='upper left', bbox_to_anchor=tuple(loc), **_kw)
    else:
        ax.legend(h, l, loc=loc, **_kw)

def _xo__save(fig, out_dir: Path, name: str, reporter) -> None:
    """Persist a figure at publication resolution and free its memory."""
    _ext_match_03_style(fig)
    path = out_dir / name
    fig.savefig(path, dpi=CFG.VIS_FIGURE_DPI, bbox_inches='tight')
    plt.close(fig)

def _xo__fig_02A_binding_affinity_metrics(df, out_dir, reporter, controls=None):
    """Binding_Affinity_Score violin + Affinity/Pocket-ratio/Density mean ± 95% CI lines.

    Mirrors Figure_01A's layout (green Degrader / orange Non-Degrader zones,
    Mann-Whitney stats box). The violin is Binding_Affinity_Score; the
    overlay lines are Affinity (custom_affinity_score), Pocket ratio
    (Pocket_Tightness_Score), and Interaction density, each min-max
    normalised to 0-1 and connected across tiers by a coloured line. Each
    line also gets a dashed horizontal reference at the DehH2+FA control's
    value (same colour). Rank-biserial r is signed so r > 0 whenever
    Degraders (Tier_1A/1B) exceed Non-Degraders on that metric.
    """
    controls = controls or {}
    ba_col = _xo__pillar_col(df, 'Binding_Affinity_Score')
    if ba_col is None:
        reporter.log(f"  ! Skipped: {_fig_path('02A')} — Binding_Affinity_Score column missing.")
        return
    aff_col = _xo__col(df, 'custom_affinity_score', 'Chemical_Affinity_Score', 'custom_affinity_calc')
    pocket_col = _xo__col(df, 'Pocket_Tightness_Score', 'pocket_enclosure_ratio')
    dens_col = _xo__col(df, CFG.COL_IDENS, 'Interaction_Density_Norm', 'interaction_density_calc')
    tiers = _xo__tiers_present(df)
    if not tiers:
        reporter.log(f"  ! Skipped: {_fig_path('02A')} — no tiers present.")
        return
    xpos = {t: i for i, t in enumerate(tiers)}
    fig, ax = plt.subplots(figsize=(10, 6))
    deg_idx = [xpos[t] for t in tiers if t in _xo_ELITE_TIERS]
    non_idx = [xpos[t] for t in tiers if t not in _xo_ELITE_TIERS]
    sns.violinplot(data=df, x=CFG.COL_TIER, y=ba_col, order=tiers, hue=CFG.COL_TIER, palette=TIER_PALETTE, legend=False, cut=0, inner='box', ax=ax, zorder=2)

    """
    The per-tier MEAN affinity, traced across the tiers. The violins carry the distributions;
    the line carries the point — affinity does NOT order the tiers, and a reader should be able
    to see that without integrating seven shapes by eye.
    """
    _means2 = [float(pd.to_numeric(df.loc[df[CFG.COL_TIER] == _t, ba_col],
                                   errors='coerce').mean()) for _t in tiers]
    ax.plot(range(len(tiers)), _means2, color=CFG.VIS_ACCENT["bad"], lw=2.2, marker='D', ms=6,
            mec='white', mew=0.8, zorder=8, label='Mean affinity')

    for _coll in ax.collections:
        _coll.set_alpha(0.6)

    def _mean_ci(values):
        v = np.asarray(values, float)
        v = v[np.isfinite(v)]
        if len(v) < 2:
            return (np.nan, np.nan) if len(v) == 0 else (float(v[0]), 0.0)
        mean = float(np.mean(v))
        sem = float(np.std(v, ddof=1) / np.sqrt(len(v)))
        ci = float(sem * _t_dist.ppf(0.975, len(v) - 1))
        return (mean, ci)
    trend_specs = [(aff_col, CFG.VIS_ACCENT["blue"], 'Affinity (95% CI)', 'o'), (pocket_col, CFG.VIS_ACCENT["amber"], 'Pocket ratio (95% CI)', 's'), (dens_col, CFG.VIS_TREND_SERIES[2], 'Int. density (95% CI)', '^')]
    df_norm = df.copy()
    for tcol, colour, tlabel, mk in trend_specs:
        if tcol is None:
            continue
        df_norm['_nrm'] = _xo__minmax(df_norm[tcol])
        xs, means, cis = ([], [], [])
        for t in tiers:
            vals = df_norm.loc[df_norm[CFG.COL_TIER] == t, '_nrm'].dropna().values
            mval, cval = _mean_ci(vals)
            if np.isfinite(mval):
                xs.append(xpos[t])
                means.append(mval)
                cis.append(cval)
        if xs:
            ax.errorbar(xs, means, yerr=cis, color=colour, marker=mk, markersize=6, lw=2.0, capsize=3, markeredgecolor='black', markeredgewidth=0.6, label=tlabel, zorder=6)
        ctrl_row = controls.get('DehH2+FA')
        if ctrl_row is not None:
            raw = pd.to_numeric(df[tcol], errors='coerce')
            lo, hi = (np.nanmin(raw), np.nanmax(raw))
            cval_raw = pd.to_numeric(pd.Series([ctrl_row.get(tcol, np.nan)]), errors='coerce').iloc[0]
            if np.isfinite(lo) and np.isfinite(hi) and (hi > lo) and pd.notna(cval_raw):
                cval_norm = float(np.clip((cval_raw - lo) / (hi - lo), 0, 1))
                ax.axhline(cval_norm, ls='--', lw=1.4, color=colour, alpha=0.85, zorder=1.5, label=f"{tlabel.split('  (')[0]} — DehH2+FA control = {cval_norm:.2f}")

    def _mw_signed(mcol, label):
        if mcol is None or mcol not in df.columns:
            return None
        v = pd.to_numeric(df[mcol], errors='coerce')
        is_deg = df[CFG.COL_TIER].isin(_xo_ELITE_TIERS)
        a = v[is_deg].dropna().values
        b = v[~is_deg].dropna().values
        if len(a) < 2 or len(b) < 2:
            return None
        U, p = mannwhitneyu(a, b, alternative='two-sided')
        r = 2.0 * U / (len(a) * len(b)) - 1.0
        # The panel prints the raw p; register it so its q_BH is paid in the same family.
        _register_p(f"Mann-Whitney U — {label} — Elite (Tier_1A/1B) vs rest",
                    "08_Binding_Affinity_Metrics", float(U), len(a) + len(b), float(p),
                    effect_size_r=round(float(r), 4))
        return (float(p), float(r))
    stat_lines = ['Elite (Tier_1A/1B) vs rest (Mann–Whitney U; q_BH in CSV; r>0 = Elite higher)']
    metric_map = [('BA_Score', ba_col), ('Affinity', aff_col), ('Pocket', pocket_col), ('IntDens', dens_col)]
    for label, mcol in metric_map:
        res = _mw_signed(mcol, label)
        if res is None:
            stat_lines.append(f'{label:<9}: n/a')
        else:
            p, r = res
            stat_lines.append(f'{label:<9}: {_xo__fmt_p(p)} | r = {r:+.2f}')
    y_max = np.nanmax(pd.to_numeric(df[ba_col], errors='coerce').values) if ba_col in df.columns else 1.0
    y_top = max(1.0, float(y_max)) * 1.08
    ax.set_xlabel('Catalytic degrader tier')
    ax.set_ylabel('Binding_Affinity_Score / normalised components (0–1)')
    _tier_seps(ax)   # consistent tier separators
    ax.set_ylim(0, y_top)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch as _PatchA
    handles, labels = ax.get_legend_handles_labels()
    handles = [_PatchA(facecolor=CFG.VIS_INK["paler"], alpha=0.6, label='BA_Score (violin)')] + handles
    labels = ['BA_Score (violin)'] + labels
    # Anchored inside the top-left corner: the default 'upper left' placement drifts out to the
    # frame and reads as a separate object floating beside the panel. Semi-transparent box so the
    # violins behind it stay visible.
    _xo__legend_with_stats(ax, handles, labels, stat_lines, (0.012, 0.985), ncol=2)
    _xo__save(fig, out_dir, '08_Binding_Affinity_Metrics.png', reporter)

def _xo__fig_04A_evolutionary_phylogeny(df, out_dir, reporter):
    idc = _xo__col(df, CFG.COL_ID_PCT, 'Identity_to_Control')
    evo = _xo__pillar_col(df, 'Evolutionary_Fingerprint_Score')
    mech = _xo__col(df, 'Mechanistic_Fingerprint_Score', 'Catalytic_Fingerprint_Score')
    rmsd = _xo__col(df, 'Active_Site_RMSD_to_Control')
    tiers = _xo__tiers_present(df)
    if not tiers or (idc is None and evo is None):
        reporter.log(f"  ! Skipped: {_fig_path('04A')} — necessary columns or tiers missing.")
        return
    xpos = {t: i for i, t in enumerate(tiers)}
    deg_idx = [xpos[t] for t in tiers if t in _xo_ELITE_TIERS]
    non_idx = [xpos[t] for t in tiers if t not in _xo_ELITE_TIERS]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))

    def _mw_signed_p1(mcol, label):
        if mcol is None or mcol not in df.columns:
            return None
        v = pd.to_numeric(df[mcol], errors='coerce')
        is_deg = df[CFG.COL_TIER].isin(_xo_ELITE_TIERS)
        a = v[is_deg].dropna().values
        b = v[~is_deg].dropna().values
        if len(a) < 2 or len(b) < 2:
            return None
        U, p = mannwhitneyu(a, b, alternative='two-sided')
        r = 2.0 * U / (len(a) * len(b)) - 1.0
        # The panel prints the raw p; register it so its q_BH is paid in the same family.
        _register_p(f"Mann-Whitney U — {label} — Elite (Tier_1A/1B) vs rest",
                    "05_Evolutionary_Phylogeny", float(U), len(a) + len(b), float(p),
                    effect_size_r=round(float(r), 4))
        return (float(p), float(r))
    if idc is not None:
        sns.violinplot(data=df, x=CFG.COL_TIER, y=idc, order=tiers, hue=CFG.COL_TIER, palette=TIER_PALETTE, legend=False, cut=0, inner='quartile', ax=ax1, zorder=2)
        for _coll in ax1.collections:
            _coll.set_alpha(0.6)
        # Seven violins ask the reader to infer the trend by eye; the mean line states it.
        _xn__mean_trend(ax1, df, tiers, idc, group_col=CFG.COL_TIER, label='Mean identity (trend)')
        ax1.set_xlabel('Catalytic degrader tier')
        ax1.set_ylabel('Sequence identity to control (%)')
        _tier_seps(ax1)   # consistent tier separators
        # Both panels share the same tier x-axis; the top panel's tick labels only repeat the bottom's,
        # so they are suppressed and the coloured tier names are shown once, under the bottom panel.
        ax1.tick_params(labelbottom=False)
        res = _mw_signed_p1(idc, 'Sequence Identity')
        stat_text = 'Elite (Tier_1A/1B) vs rest (Mann–Whitney U;  uncorrected p, q_BH in 06_Statistical_Tests.csv)\n'
        if res:
            p, r = res
            stat_text += f'Sequence Identity: {_xo__fmt_p(p)} | r = {r:+.2f}'
        else:
            stat_text += 'Sequence Identity: n/a'
        _xo__annotate(ax1, stat_text, loc='lower left')
    if evo is not None:
        sns.violinplot(data=df, x=CFG.COL_TIER, y=evo, order=tiers, hue=CFG.COL_TIER, palette=TIER_PALETTE, legend=False, cut=0, inner='box', ax=ax2, zorder=2)
        for _coll in ax2.collections:
            _coll.set_alpha(0.6)

        def _mean_ci(values):
            v = np.asarray(values, float)
            v = v[np.isfinite(v)]
            if len(v) < 2:
                return (np.nan, np.nan) if len(v) == 0 else (float(v[0]), 0.0)
            mean = float(np.mean(v))
            sem = float(np.std(v, ddof=1) / np.sqrt(len(v)))
            ci = float(sem * _t_dist.ppf(0.975, len(v) - 1))
            return (mean, ci)
        trend_specs = [(mech, CFG.VIS_ACCENT["blue"], 'Mechanistic fingerprint  (norm., mean ± 95% CI)', 'o'), (rmsd, CFG.VIS_ACCENT["amber"], 'Active site RMSD  (norm., mean ± 95% CI)', 's')]
        """
        THE NORMALISED TRACES GET THEIR OWN AXIS.

        The violins are the evolutionary fingerprint score, which runs 0-100. The overlaid traces are
        min-max normalised to 0-1. Drawn on the same axis, a 0-1 series against a 0-100 scale is a flat
        line pinned to the floor: it was plotted, it was legible as a colour, and it carried no
        information at all. They now share the x-axis and nothing else.

        The RMSD column also carries a 999.0 sentinel for 'undefined' (max = 999.00 in the ranked sheet).
        Min-max normalising through that sentinel crushes every real value into the bottom of the range,
        so it is dropped before the scaling rather than clipped after it — a clip to 5 A still lets a
        sentinel-derived 5.0 masquerade as a genuinely bad fold.
        """
        # Twin-axis colour coding: the two y-axes carry different quantities (evolutionary score 0–100 on
        # the left, normalised traces 0–1 on the right), so each takes its own CFG colour on its label,
        # ticks, spine and grid — a reader can then never trace a value onto the wrong scale. The right
        # grid is dashed so, where the two colour grids overlap, they stay tellable apart.
        _c_left, _c_right = CFG.VIS_ACCENT["axis_left"], CFG.VIS_ACCENT["axis_right"]
        ax2r = ax2.twinx()
        ax2r.set_ylim(0, 1.05)
        ax2r._twin_axis_colour = _c_right          # honoured by _ext_match_03_style (keeps this colour)
        ax2r.set_ylabel('normalised (0–1)   ·   RMSD inverted, 1.0 = best', color=_c_right)
        ax2r.tick_params(axis='y', labelsize=CFG.VIS_FONT_TICK, colors=_c_right)
        ax2r.spines['right'].set_color(_c_right)
        ax2r.set_axisbelow(True)
        ax2r.grid(True, axis='y', color=_c_right, alpha=CFG.VIS_GRID_ALPHA,
                  linewidth=CFG.VIS_GRID_LINEWIDTH, linestyle='--')

        df_norm = df.copy()
        for tcol, colour, tlabel, mk in trend_specs:
            if tcol is None:
                continue
            if 'RMSD' in tcol:
                rmsd_vals = pd.to_numeric(df_norm[tcol], errors='coerce')
                rmsd_vals = rmsd_vals.where(rmsd_vals < float(CFG.SENTINEL_UNDEFINED))  # drop 999 = undefined
                rmsd_norm = _xo__minmax(rmsd_vals.clip(upper=5.0))
                df_norm['_nrm'] = 1.0 - rmsd_norm
                tlabel = tlabel.replace('norm.', 'norm. inverted, 1.0=best')
            else:
                df_norm['_nrm'] = _xo__minmax(df_norm[tcol])
            xs, means, cis = ([], [], [])
            for t in tiers:
                vals = df_norm.loc[df_norm[CFG.COL_TIER] == t, '_nrm'].dropna().values
                mval, cval = _mean_ci(vals)
                if np.isfinite(mval):
                    xs.append(xpos[t])
                    means.append(mval)
                    cis.append(cval)
            if xs:
                ax2r.errorbar(xs, means, yerr=cis, color=colour, marker=mk, markersize=6, lw=2.0,
                              capsize=3, markeredgecolor='black', markeredgewidth=0.6, label=tlabel,
                              zorder=6)
        stat_lines = ['Elite (Tier_1A/1B) vs rest  (Mann–Whitney U;  uncorrected p, q_BH in 06_Statistical_Tests.csv;  r > 0 = Elite higher)']
        metric_map = [('Evo_Score', evo), ('Mech_Fpt', mech), ('Active_RMDA', rmsd)]
        for label, mcol in metric_map:
            res = _mw_signed_p1(mcol, label)
            if res is None:
                stat_lines.append(f'{label:<12}: n/a')
            else:
                p, r = res
                stat_lines.append(f'{label:<12}: {_xo__fmt_p(p)} | r = {r:+.2f}')
        y_max = np.nanmax(pd.to_numeric(df[evo], errors='coerce').values) if evo in df.columns else 1.0
        y_top = max(1.0, float(y_max)) * 1.08
        _xn__mean_trend(ax2, df, tiers, evo, group_col=CFG.COL_TIER, label='Mean evo score (trend)')
        ax2.set_xlabel('Catalytic degrader tier')
        ax2.set_ylabel('Evolutionary fingerprint score')
        _tier_seps(ax2)   # consistent tier separators
        ax2.set_ylim(0, y_top)
        # Left axis (evolutionary score) takes the paired left colour on its label, ticks, spine and grid,
        # so it reads as the counterpart of the amber right axis. Solid grid vs the right's dashed grid.
        ax2._twin_axis_colour = _c_left            # honoured by _ext_match_03_style (keeps this colour)
        ax2.spines['left'].set_color(_c_left)
        ax2.tick_params(axis='y', colors=_c_left)
        ax2.set_axisbelow(True)
        ax2.grid(True, axis='y', color=_c_left, alpha=CFG.VIS_GRID_ALPHA,
                 linewidth=CFG.VIS_GRID_LINEWIDTH, linestyle='-')
        from matplotlib.lines import Line2D
        from matplotlib.patches import Patch as _PatchA
        handles, labels = ax2.get_legend_handles_labels()
        handles = [_PatchA(facecolor=CFG.VIS_INK["paler"], alpha=0.6, label='Evolutionary_Fingerprint_Score (violin)')] + handles
        labels = ['Evolutionary_Fingerprint_Score (violin)'] + labels
        _xo__legend_with_stats(ax2, handles, labels, stat_lines, 'lower left')
    plt.tight_layout()
    _xo__save(fig, out_dir, '05_Evolutionary_Phylogeny.png', reporter)

def _xo__fig_05b_mechanistic_size_modified(df, out_dir, reporter):
    fcol = _xo__col(df, 'total_fluorine_count') if '_xo__col' in globals() else _xo_get_col(df, 'total_fluorine_count')
    if fcol is None:
        return
    m_all = pd.Series(True, index=df.index)
    inter = _xo__col(df, 'interacting_fluorine_count') if '_xo__col' in globals() else _xo_get_col(df, 'interacting_fluorine_count')
    tot = fcol
    if inter and tot:
        df['FER_computed'] = pd.to_numeric(df[inter], errors='coerce') / pd.to_numeric(df[tot], errors='coerce')
    else:
        df['FER_computed'] = np.nan
    rows = [(_xo__col(df, CFG.COL_SN2) if '_xo__col' in globals() else _xo_get_col(df, CFG.COL_SN2), 'SN2 attack angle (°)', m_all), ('nuc_dist', 'Nucleophile distance (Å)', m_all), ('FER_computed', 'Fluorine Engagement Ratio', m_all)]
    fig, axes = plt.subplots(3, 1, figsize=(8.5, 12), sharex=True)
    x_all = pd.to_numeric(df[fcol], errors='coerce')
    # Drop the "no nucleophile found" sentinel (~999/1000 Å) so the panel shows real
    # catalytic distances (2–8 Å) rather than being crushed to the axis floor.
    nuc = _xo__nuc_distance(df).where(lambda s: s < 20.0)
    unique_x = sorted(x_all.dropna().unique())
    for i, (col, ylab, mask) in enumerate(rows):
        ax = axes[i]
        if col == 'nuc_dist':
            y = nuc
        elif col is None or col not in df.columns:
            ax.set_visible(False)
            continue
        else:
            y = pd.to_numeric(df[col], errors='coerce')
        xx = x_all[mask]
        yy = y[mask]
        ok = xx.notna() & yy.notna()
        xx, yy = (xx[ok], yy[ok])
        if len(xx) < 5:
            ax.set_visible(False)
            continue
        temp_df = pd.DataFrame({'x': xx, 'y': yy})
        sns.boxplot(data=temp_df, x='x', y='y', ax=ax, color=CFG.VIS_ACCENT["blue"], fliersize=1)
        grp = temp_df.groupby('x')['y'].mean()
        valid_x = sorted(temp_df['x'].unique())
        x_indices = [unique_x.index(x_val) for x_val in valid_x]
        ax.plot(x_indices, grp.values, '-', color=CFG.VIS_ACCENT["vermillion"], lw=2, marker='o')
        ax.set_ylabel(ylab)
        """
        The predictor is a COUNT, so it is tested as one.

        Splitting fluorine count at F≤3 vs F>3 and running Mann-Whitney throws away everything the
        count knows: a difluoro and a perfluorodecyl land in the same bin, the monotonic trend the panel
        is drawn to show is discarded, and the power lost to the binarisation is paid for nothing. The
        cut point was also arbitrary — no threshold in CFG corresponds to it.

        Spearman's rho on the continuous count answers the question the figure asks — does the metric
        move monotonically with chain length — and reports the DIRECTION and STRENGTH of that move, not
        merely whether two arbitrary halves differ. It is registered like every other test, so it joins
        the Benjamini-Hochberg family instead of being a p-value printed on a panel and corrected
        nowhere.
        """
        if col in ['nuc_dist', rows[0][0]]:
            _fin = xx.notna() & yy.notna()
            if int(_fin.sum()) >= 3 and xx[_fin].nunique() > 1:
                _rho, _p = _sc_stats.spearmanr(xx[_fin], yy[_fin])
                if np.isfinite(_rho) and np.isfinite(_p):
                    _register_p("Spearman correlation", f"{ylab} vs total fluorine count",
                                float(_rho), int(_fin.sum()), float(_p),
                                effect_size=round(float(_rho), 4), effect_type="spearman rho")
                    p_str = 'p < 0.001' if _p < 0.001 else f'p = {_p:.3f}'
                    # Marked uncorrected: the BH family closes only after every figure has run, so no
                    # q exists at draw time. The corrected q_BH is in 06_Statistical_Tests.csv.
                    ax.text(0.95, 0.95, f'Spearman ρ = {_rho:+.2f}  ·  {p_str} uncorrected  ·  n = {int(_fin.sum()):,}',
                            transform=ax.transAxes, ha='right', va='top',
                            fontsize=CFG.VIS_FONT_ANNOT,
                            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    axes[-1].set_xticks(np.arange(len(unique_x)))
    axes[-1].set_xticklabels([str(int(val)) for val in unique_x])
    axes[-1].set_xlabel('Total fluorine count')
    for _ax in axes:
        _tier_seps(_ax)   # consistent per-class separators
    fig.tight_layout()
    _xo__save(fig, out_dir, '14_Mechanistic_Breakdown_by_Tier.png', reporter)

def _xo__fig_05c_size_by_tier_modified(df, out_dir, reporter):
    fcol = _xo__col(df, 'total_fluorine_count') if '_xo__col' in globals() else _xo_get_col(df, 'total_fluorine_count')
    if fcol is None:
        return
    tiers = _xo__tiers_present(df) if '_xo__tiers_present' in globals() else sorted(df[CFG.COL_TIER].dropna().unique())
    fig, ax = plt.subplots(figsize=(10, 6))
    sns.violinplot(data=df, x=CFG.COL_TIER, y=fcol, order=tiers, hue=CFG.COL_TIER, palette=TIER_PALETTE if 'TIER_PALETTE' in globals() else None, legend=False, cut=0, inner='box', ax=ax)
    degrader_tiers = [t for t in tiers if t in _xo_ELITE_TIERS] if '_xo_ELITE_TIERS' in globals() else [t for t in tiers if t in ['Tier_1A', 'Tier_1B']]
    non_degrader_tiers = [t for t in tiers if t not in degrader_tiers]
    xpos = {t: i for i, t in enumerate(tiers)}
    deg_idx = [xpos[t] for t in degrader_tiers]
    non_idx = [xpos[t] for t in non_degrader_tiers]
    """
    The per-tier MEAN, traced across the tiers. A violin shows each tier's shape but leaves the
    reader to compare seven of them by eye; the trend line states the claim the figure exists to
    make — chain length rises monotonically as the tier falls — in one stroke.
    """
    _means7 = [float(df.loc[df[CFG.COL_TIER] == _t, fcol].mean()) for _t in tiers]
    _xs7 = list(range(len(tiers)))
    ax.plot(_xs7, _means7, color=CFG.VIS_ACCENT["bad"], lw=2.2, marker='D', ms=6,
            mec='white', mew=0.8, zorder=8, label='Mean (trend)')
    ax.legend(loc='upper right')

    ax.set_xlabel('Catalytic degrader tier')
    ax.set_ylabel('Total fluorine count')
    _tier_seps(ax)   # consistent tier separators
    df_deg = df[df[CFG.COL_TIER].isin(degrader_tiers)][fcol].dropna()
    df_non = df[df[CFG.COL_TIER].isin(non_degrader_tiers)][fcol].dropna()
    if len(df_deg) > 0 and len(df_non) > 0:
        stat, p = mannwhitneyu(df_deg, df_non, alternative='two-sided')
        n1, n2 = (len(df_deg), len(df_non))
        r = 1 - 2 * stat / (n1 * n2)
        _register_p("Mann-Whitney U", "Total fluorine count: degrader vs non-degrader tiers",
                    float(stat), int(n1 + n2), float(p),
                    effect_size_r=float(r), n_degrader=int(n1), n_non_degrader=int(n2),
                    median_degrader=float(df_deg.median()), median_non_degrader=float(df_non.median()))
        p_str = f'p < 0.001' if p < 0.001 else f'p = {p:.3f}'
        # Marked uncorrected: the BH family closes only after every figure has run, so no q exists at
        # draw time. The corrected q_BH is in 06_Statistical_Tests.csv. One line, top-left: three
        # stacked lines in the top-right corner sit over the widest violins.
        ax.text(0.015, 0.97, f'Mann-Whitney U · {p_str} uncorrected · effect size r = {r:.2f}',
                transform=ax.transAxes, ha='left', va='top', fontsize=8.5,
                bbox=dict(facecolor='white', alpha=0.85, edgecolor=CFG.VIS_INK["palest"], linewidth=0.6,
                          boxstyle='round,pad=0.3'))
    _xo__save(fig, out_dir, '14_Chain_Length_by_Tier.png', reporter)


"""
Each extended panel is filed with the figures it belongs with, not in a holding folder of its own:
the phylogeny sits with the cohort overview, the two size panels with the PFAS scope, and so on. The
number continues the destination folder's existing sequence, so nothing already numbered moves.

The numbering is APPENDED rather than gap-filled. 07_PFAS_Scope_and_Synthesis is missing a 06 — a
retired figure — and reusing that slot would silently point an old citation at a new figure.
"""


def generate_ramachandran_figures(prod_dir: Path, out_dir: Path, reporter: ReportManager):
    """
    Backbone-geometry validation of the control predictions against the 3R3U crystal.
    Renders the crystal Ramachandran plot plus, for each Boltz-2 control (DeHa4 and 3R3U
    sequences), a standalone plot and a crystal-overlay comparison. All plots are written
    to <Run>/3_Validation_Figures/02_Ramachandran. The 3R3U crystal PDB is read from the
    input folder produced by 02_Production; nothing is downloaded here.
    """
    import gemmi
    reporter.section("Step 2/8 — 02_Ramachandran · Backbone-Geometry Validation")
    # `out_dir` is already the 02_Ramachandran folder (created in main); use it directly.
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
                                      rama_dir / "01_Ramachandran_3R3U_Crystal.png", dpi=dpi)
    reporter.log("  ✔ Saved: 02_Ramachandran/01_Ramachandran_3R3U_Crystal.png")

    _rn = 1   # running plot number within 02_Ramachandran (crystal was 01)
    for label, token in [("DeHa4", "DeHa4_Control"), ("3R3U", "3R3U_Control")]:
        cif = _control_cif(token)
        if not cif:
            reporter.log(f"  ! Ramachandran: no {label} control CIF found")
            continue
        rama_con = _utils_mod.compute_ramachandran_angles(_load(cif))
        _rn += 1; _c_con = f"{_rn:02d}_Ramachandran_{label}_Control.png"
        _utils_mod.save_ramachandran_plot(rama_con, f"{label} Control (Boltz-2)",
                                          rama_dir / _c_con, dpi=dpi)
        _rn += 1; _c_cmp = f"{_rn:02d}_Ramachandran_{label}_vs_Crystal.png"
        _utils_mod.save_ramachandran_comparison(rama_ref, rama_con,
                                                "3R3U (Crystal)", f"{label} (Boltz-2)",
                                                rama_dir / _c_cmp, dpi=dpi)
        reporter.log(f"  ✔ Saved: 02_Ramachandran/{_c_con}")
        reporter.log(f"  ✔ Saved: 02_Ramachandran/{_c_cmp}")


# =============================================================================
# SECTION 5: FIGURE DESCRIPTIONS & REPORTING
# =============================================================================

"""
Why a defined figure can legitimately draw nothing. Keyed by filename; used by the legend writer so a
skipped figure is reported with its cause instead of a bare absence.

A Hidden Gem is a complex in a high-quality tier whose Boltz-2 confidence is nonetheless below
CFG.CONFLICT_CONF_HIGH — physics says yes, the AI hedges. When the panel is empty it means no such
disagreement exists: on this corpus every surviving elite complex also carries high model confidence,
so there is nothing to rescue. An empty rescue list is the strong outcome, not a missing figure.
"""
_ABSENT_FIGURE_REASON: dict = {
    "06_Hidden_Gems_DeepDive.png":
        "No Hidden Gems: no complex sits in a high-quality tier while falling below the "
        "confidence threshold. Physics and model confidence agree across every elite hit.",
}


def write_figure_descriptions(out_dir: Path):
    """
    Writes a human-readable log file describing every figure produced by this pipeline.

    Entries follow the on-disk folder-by-folder layout:
      02_Ramachandran/                          backbone-geometry validation of controls
      03_Dataset_and_Alignment_Overview/        dataset + sequence-alignment overview
      04_AI_Confidence_Quality/                 Boltz-2 confidence metrics
      05_Catalytic_Geometry_and_Mechanism/      structural + mechanistic geometry
      06_Ligand_Interactions_and_Chemical_Space/ interaction profile + chemical space
      07_PFAS_Scope_and_Synthesis/              multi-metric synthesis + publication assembly
      08_Diagnostic_and_MultiModel_Trends/      pocket-fit + multi-model consensus
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
        "  Features: 17 metrics grouped by curated family (outcome, catalytic geometry, chemistry /",
        "            mechanism, affinity / seq / rank, binding & confidence). The chemistry / mechanism",
        "            drivers (C-F BDE, Backside, Halide, Compet., Cbx clamp, Mech.eff) are real matrix",
        "            cells — their actual pairwise Spearman with every other feature. Each family's",
        "            tick labels + bracket are colour-coded.",
        "  Look for: Strong red cells = features that rise and fall together (redundant or causal).",
        "            Strong blue = features that are inversely related.",
        "            AI Confidence + Binding Prob. cluster together — same underlying signal.",
        "            Pareto Rank inverts against physics/AI scores (lower rank = better = higher score).",
        "            Competence + Mech.eff track the Tier outcome; C-F BDE opposes it.",
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
        "  X-axis  : mapped catalytic Asp (nucleophile, per-variant) to electrophilic carbon distance (Ang); shorter = reaction-ready",
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
        "  Axes    : X = SN2_Attack_Angle (degrees, higher = near-ideal); Y = Boltz_Model_Confidence (AI structural quality)",
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
        "  Title   : Chemical Space Map (UMAP manifold, competence hexbin)",
        f"  Type    : Hexbin over the UMAP embedding coloured by mean competence per bin; {CFG.TIER_TOP} ★ overlaid, with inset PyMOL active-site thumbnails",
        "  Axes    : UMAP dimensions 1 & 2 — distances reflect chemical similarity",
        "  Colour  : Mean competence score per bin (brighter = more degradable)",
        f"  Markers : ★ = MD-ready / {CFG.TIER_TOP} hits, each ★ tinted to match its own thumbnail's box + arrow; the thumbnail names the protein + ligand and shows its active-site geometry",
        f"  Look for: Bright chemical regions enriched for degraders, with the {CFG.TIER_TOP} stars clustering there.",
        "            Thumbnails reveal the active-site geometry behind each top hit at a glance.",
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
        "DIAGNOSTIC & MULTI-MODEL TRENDS (subfolder: 08_Diagnostic_and_MultiModel_Trends/)",
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
        "02_Pocket_Occupancy_by_Carbon_Number.png",
        "  Title   : Ligand containment by PFAS carbon number",
        "  Type    : Per carbon group (C2…Cn, fluorine counts in parentheses) — whole-cavity",
        "            coverage as a violin (full distribution) with the 8-residue active-site",
        "            coverage as a box (median · IQR) inside it; both group means are trended",
        "  Y-axis  : Fraction of the ligand contained (0–1); dashed line at 1.0 = fully contained",
        "  Look for: both measures fall as the chain lengthens, and the active-site box drops",
        "            faster than the cavity violin — the cavity still holds a long chain after the",
        "            catalytic shell has lost it. That divergence is the enzyme's size selectivity.",
        "",
        "-" * 80,
        "  Title   : Catalytic competence vs pocket occupancy",
        "  Type    : Scatter coloured by tier + binned-median trend (red diamonds)",
        "  Axes    : X = pocket occupancy; Y = catalytic competence score",
        "  Look for: whether competence falls as the cavity crowds — the steric signal that",
        "            informs the feasibility weighting.",
        "",
        "-" * 80,
        "  Title   : Steric fit rate by ligand",
        "  Type    : Horizontal bars, sorted; colour ramps no-fit (maroon) → fit (teal)",
        "  X-axis  : % of complexes passing the steric fit test (ligand_fits == True)",
        "  Look for: small ligands pass universally; fit rate erodes for the largest PFAS.",
        "",
        "-" * 80,
        "03_MultiModel_Consensus_by_Tier.png",
        "  Title   : Multi-model degrader consensus by tier",
        "  Type    : Strip + tier mean (diamond) + 95% CI; coloured tier ticks with μ/median/n",
        "  Y-axis  : Fraction of the 5 Boltz-2 diffusion models that independently called a degrader",
        f"  Look for: consensus decreases monotonically from {CFG.TIER_TOP} to {CFG.TIER_DECOY}, "
        "confirming tier ordering reflects cross-model agreement.",
        "",
        "-" * 80,
        "04_Confidence_vs_Consensus.png",
        "  Title   : AI confidence vs cross-model consensus",
        "  Type    : Density hexbin (log count) + binned-median trend",
        "  Axes    : X = Boltz-2 model confidence; Y = multi-model degrader consensus",
        "  Look for: high confidence usually coincides with model agreement; confident calls that",
        "            still split the ensemble are robustness outliers.",
        "",
        "-" * 80,
        "05_Quality_and_Competence_Diagnostics.png",
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
        "06_Size_Preference_Containment.png",
        "  Title   : Size preference — effective mechanistic score vs ligand size",
        "  Type    : Per-size-bin distribution + bin means + catalytic hit-rate + mean containment",
        "  Look for: competence falls with ligand size; containment falls with it too. FAcD is a",
        "            small-substrate hydrolase, and this is where that shows.",
        "",
        "-" * 80,
        "07_Reactive_Engagement.png",
        "  Title   : Reactive-carbon engagement with the catalytic residues by carbon number",
        "  Type    : Nucleophile-distance box + 8-residue spread + properly-positioned fraction",
        "  Look for: the properly-positioned fraction collapses beyond ~C6 while the ligand is",
        "            still nominally bound — reach is not engagement.",
        "",
        "-" * 80,
        "08_Model_Agreement.png",
        "  Title   : Diffusion-sample consensus by tier — are the elite hits reproducible?",
        "  Type    : Consensus distribution per tier (violin + box), per-tier mean, MD-selected starred",
        "  X-axis  : Degrader tier; Y-axis: fraction of diffusion samples independently reaching a",
        "            degrader tier (model_degrader_consensus)",
        "  Look for: whether Tier_1A is consensus-backed. A tier whose hits sit BELOW 0.5 earned",
        "            its label from a minority of samples — a best-of-N over the diffusion ensemble,",
        "            not a reproducible property of the complex.",
        "",
        "=" * 80,
        "PER-FOLDER ANALYSIS PANELS (rendered within each folder's step)",
        "=" * 80,
        "",
        "-" * 80,
        "05_Catalytic_Geometry_and_Mechanism/13_Geometry_and_Uncertainty.png",
        "  Title   : Reaction geometry with multi-model uncertainty",
        "  Look for: nucleophile distance and SN2 angle per tier, with the spread across the",
        "            diffusion samples — a tight tier is a reproducible one.",
        "",
        "-" * 80,
        "06_Ligand_Interactions_and_Chemical_Space/08_Binding_Affinity_Metrics.png",
        "  Title   : Binding affinity by tier",
        "  Look for: affinity does NOT order the tiers — a high-affinity binder that presents the",
        "            wrong face to the mapped catalytic aspartate is not a degrader. This figure is the evidence.",
        "",
        "-" * 80,
        "03_Dataset_and_Alignment_Overview/05_Evolutionary_Phylogeny.png",
        "  Title   : Evolutionary phylogeny of the cohort",
        "  Look for: whether the degrader tiers cluster phylogenetically or are scattered across",
        "            the tree. Scattered = catalytic competence is not a clade property.",
        "",
        "-" * 80,
        "08_Diagnostic_and_MultiModel_Trends/09_Pillar_Divergence_by_Tier.png",
        "  Title   : Divergence of the scoring pillars across tiers",
        "  Look for: which pillar actually separates the tiers, and which merely follows.",
        "",
        "-" * 80,
        "05_Catalytic_Geometry_and_Mechanism/14_Mechanistic_Breakdown_by_Tier.png",
        "  Title   : Mechanistic components per tier",
        "  Look for: the component that collapses first as the tier falls.",
        "",
        "-" * 80,
        "07_PFAS_Scope_and_Synthesis/14_Chain_Length_by_Tier.png",
        "  Title   : PFAS chain length by tier",
        "  Look for: the elite tiers are short-chain. FAcD is a small-substrate hydrolase.",
        "",
        "-" * 80,
        "07_PFAS_Scope_and_Synthesis/15_Tier1A_Cross_Ligand_Heatmap.png",
        "  Title   : Tier_1A proteins x ligands",
        "  Look for: whether an elite protein is elite for ONE ligand or several — a protein that",
        "            is Tier_1A across ligands is a genuinely promiscuous defluorinase.",
        "",
        "-" * 80,
        "05_Catalytic_Geometry_and_Mechanism/06_Mechanistic_Fingerprint.png",
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
        "06_Ligand_Interactions_and_Chemical_Space/06_Binding_Energetics.png",
        "  Title   : Binding energetics — binding probability by tier",
        "  Type    : Binding-probability violin per tier (red median bar)",
        "  Look for: how binding probability (sigmoid of interaction density, cross-PAE and",
        "            across tiers. Product inhibition is assessed downstream by the Step-06 MM-GBSA stage.",
        "",
        "=" * 80,
        "PART 9 — TWO-CRITERIA TIER LOGIC & DEAD-END FEASIBILITY",
        "=" * 80,
        "",
        "-" * 80,
        "05_Catalytic_Geometry_and_Mechanism/10_Criterion_A_Gates_B.png",
        "05_Catalytic_Geometry_and_Mechanism/11_Criterion_B_ECDF_by_Tier.png",
        "05_Catalytic_Geometry_and_Mechanism/12_SN2_DeadEnd_Gate.png",
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
    # in thematic reading order; each carries its folder path.
    _fig_paths = {
        "Figure_01_Active_Site_Residue_Mapping_Coverage.png": "03_Dataset_and_Alignment_Overview/01_Active_Site_Residue_Mapping_Coverage.png",
        "Figure_02_Tier_Distribution.png": "03_Dataset_and_Alignment_Overview/02_Tier_Distribution.png",
        "Figure_03_Alignment_Grades.png": "03_Dataset_and_Alignment_Overview/03_Alignment_Grades.png",
        "Figure_04_Tier_Grade_Distribution.png": "03_Dataset_and_Alignment_Overview/04_Tier_Grade_Distribution.png",
        "Figure_05a_AI_Quality_Assessment.png": "04_AI_Confidence_Quality/01_AI_Quality_Assessment.png",
        "Figure_05b_TT_AI_Quality_Space.png": "04_AI_Confidence_Quality/02_AI_Quality_Space.png",
        "Figure_06_pTM_vs_ipTM_by_Tier.png": "04_AI_Confidence_Quality/03_pTM_vs_ipTM_by_Tier.png",
        "Figure_07_ActiveSite_RMSD_by_Tier.png": "05_Catalytic_Geometry_and_Mechanism/01_ActiveSite_RMSD_by_Tier.png",
        "Figure_08_Feature_Correlations.png": "05_Catalytic_Geometry_and_Mechanism/02_Feature_Correlations.png",
        "Figure_09_Tier_Quality_DotPlot.png": "05_Catalytic_Geometry_and_Mechanism/03_Tier_Quality_DotPlot.png",
        "Figure_10_Mech_State_CrossTab.png": "05_Catalytic_Geometry_and_Mechanism/04_Mech_State_CrossTab.png",
        "Figure_11_Mechanistic_Score_by_Tier.png": "05_Catalytic_Geometry_and_Mechanism/05_Mechanistic_Score_by_Tier.png",
        "Figure_12_SN2_Angle_by_Tier.png": "05_Catalytic_Geometry_and_Mechanism/07_SN2_Angle_by_Tier.png",
        "Figure_13a_Mechanism_Geometry_Scatter.png": "05_Catalytic_Geometry_and_Mechanism/08_Mechanism_Geometry_Scatter.png",
        "Figure_13b_TT_Mechanistic_Quality_Space.png": "05_Catalytic_Geometry_and_Mechanism/09_Mechanistic_Quality_Space.png",
        "Figure_14a_Molecular_Interaction_Profile.png": "06_Ligand_Interactions_and_Chemical_Space/01_Molecular_Interaction_Profile.png",
        "Figure_14b_TT_Interaction_Quality_Space.png": "06_Ligand_Interactions_and_Chemical_Space/02_Interaction_Quality_Space.png",
        "Figure_15_Fluorine_Engagement_by_Tier.png": "06_Ligand_Interactions_and_Chemical_Space/03_Fluorine_Engagement_by_Tier.png",
        "Figure_16_Catalytic_Quality_vs_Inhibition.png": "06_Ligand_Interactions_and_Chemical_Space/04_Catalytic_Quality_vs_Inhibition.png",
        "Figure_17_ActiveSite_Contact_Density_by_Tier.png": "06_Ligand_Interactions_and_Chemical_Space/05_ActiveSite_Contact_Density_by_Tier.png",
        "Figure_18a_Chemical_Space_Map.png": "06_Ligand_Interactions_and_Chemical_Space/07_Chemical_Space_Map.png",
        "Figure_19a_Radar_TopHits.png": "07_PFAS_Scope_and_Synthesis/01_Radar_TopHits.png",
        "Figure_19b_Radar_TierReps.png": "07_PFAS_Scope_and_Synthesis/02_Radar_TierReps.png",
        "Figure_20a_Tier_Success_Rates.png": "07_PFAS_Scope_and_Synthesis/03_Tier_Success_Rates.png",
        "Figure_20b_Conf_SN2_Landscape.png": "07_PFAS_Scope_and_Synthesis/04_Conf_SN2_Landscape.png",
        "Figure_21_Conflict_Composition.png": "07_PFAS_Scope_and_Synthesis/05_Conflict_Composition.png",
        "Figure_22_Hidden_Gems_DeepDive.png": "07_PFAS_Scope_and_Synthesis/06_Hidden_Gems_DeepDive.png",
        "Figure_23_Category_Overlap_Euler.png": "07_PFAS_Scope_and_Synthesis/07_Category_Overlap_Euler.png",
        "Figure_24a_Top25_Multitarget_Proteins.png": "07_PFAS_Scope_and_Synthesis/08_Top25_Multitarget_Proteins.png",
        "Figure_24b_TopTier_Protein_PFAS_Breakdown.png": "07_PFAS_Scope_and_Synthesis/09_TopTier_Protein_PFAS_Breakdown.png",
        "Figure_25_Sankey_Workflow.png": "07_PFAS_Scope_and_Synthesis/10_Sankey_Workflow.png",
        "Figure_26a_PFAS_Size_Hexbin_Landscape.png": "07_PFAS_Scope_and_Synthesis/11_PFAS_Size_Hexbin_Landscape.png",
        "Figure_26b_PFAS_Size_Composition_Merged.png": "07_PFAS_Scope_and_Synthesis/12_PFAS_Size_Composition_Merged.png",
        "Figure_26c_PFAS_Carbon_Confidence.png": "07_PFAS_Scope_and_Synthesis/13_PFAS_Carbon_Confidence.png",
        "02_Pocket_Occupancy_by_Ligand.png": "02_Pocket_Occupancy_by_Carbon_Number.png",
    }
    # Additional journal figure types (folded panels). Colours, fonts and DPI are the CFG values every
    # other figure uses; no titles. Full folder paths are given so each entry sits with its figure.
    lines += [
        "",
        "=" * 80,
        "PART 8 — ADDITIONAL FIGURE TYPES",
        "=" * 80,
        "",
        "-" * 80,
        "Figure — 03_Dataset_and_Alignment_Overview/06_Treemap_Tier_Ligand_Composition.png",
        "  Type    : Treemap (squarified) — one rectangle per tier × ligand cell",
        "  Area    : Number of protein–ligand candidates in that tier for that ligand",
        "  Colour  : Degrader tier (CFG.TIER_COLOUR)",
        "  Look for: Which ligands dominate each tier by volume; the largest blocks are the",
        "            most-populated tier/ligand combinations.",
        "",
        "-" * 80,
        "Figure — 05_Catalytic_Geometry_and_Mechanism/15_Metric_CoVariation_Network.png",
        "  Type    : Force-directed network of the scoring metrics",
        "  Nodes   : Metrics; edges drawn where |Pearson r| ≥ 0.5; node size ∝ degree",
        "  Colour  : Connected component (CFG Okabe-Ito accents)",
        "  Look for: Clusters of metrics that co-vary — redundant vs. independent signal.",
        "",
        "-" * 80,
        "Figure — 07_PFAS_Scope_and_Synthesis/16_Swimmer_Top_Per_Tier.png",
        "  Type    : Swimmer (horizontal bars) — top 3 candidates per tier (Tier_1A → Tier_5_Decoy)",
        "  X-axis  : Effective mechanistic score;  bar colour = tier (CFG.TIER_COLOUR)",
        "  Markers : ▶ MD-selected, ● degrader (CFG.VIS_ACCENT)",
        "  Look for: The score gradient down the tiers and which representatives were MD-selected.",
        "",
        "-" * 80,
        "Figure — 07_PFAS_Scope_and_Synthesis/17_Circos_Ligand_Tier_Assignment.png",
        "  Type    : Circos / chord diagram",
        "  Arcs    : Ligand → tier assignment; ribbon width ∝ candidate count",
        "  Colour  : Tier (CFG.TIER_COLOUR); grey ticks = ligands",
        "  Look for: Which tiers each ligand feeds, and the dominant ligand→tier flows.",
        "",
        "-" * 80,
        "Figure — 08_Diagnostic_and_MultiModel_Trends/10_Volcano_Metric_Significance.png",
        "  Type    : Volcano plot — one point per scoring metric",
        "  X-axis  : Effect size (Cohen's d, degrader − non-degrader)",
        "  Y-axis  : −log10 p (Mann–Whitney, degrader vs non-degrader)",
        "  Colour  : Higher in degraders / lower / n.s. (CFG.VIS_ACCENT + CFG.VIS_INK)",
        "  Look for: Which metrics most separate degraders from non-degraders.",
        "",
        "-" * 80,
        "Figure — 08_Diagnostic_and_MultiModel_Trends/11_Manhattan_Candidate_Significance.png",
        "  Type    : Manhattan plot — one point per candidate, grouped by ligand 'locus'",
        "  X-axis  : Candidates ordered by ligand;  Y-axis : rank-based −log10 significance",
        "  Colour  : Ligand locus (CFG Okabe-Ito accents); dashed line = p 0.05",
        "  Look for: Ligands with many high-significance candidates rising above the line.",
        "",
        "-" * 80,
        "Figure — 08_Diagnostic_and_MultiModel_Trends/12_Feature_Importance_Tier_Streams.png",
        "  Type    : Feature-importance bars (left) + tier-composition streamgraph (right)",
        "  Left    : |PC1 loading| per metric (importance); Right : tier stream across ligands",
        "  Colour  : Degrader tier (CFG.TIER_COLOUR)",
        "  Look for: Which metrics drive PC1, and how tier composition varies across ligands.",
        "",
        "-" * 80,
        "Figure — 08_Diagnostic_and_MultiModel_Trends/13_MultiModel_Geometry_Variance_by_Model.png",
        "  Type    : Per-Boltz-model SN2 distance (top) + angle (bottom) box-by-tier",
        "  X-axis  : Degrader tier;  Hue : the 5 Boltz diffusion models;  dashed = strict-NAC references",
        "  Look for: Whether the 5 models agree on the geometry within a tier (tight boxes = consistent).",
        "",
        "-" * 80,
        "Figure — 08_Diagnostic_and_MultiModel_Trends/14_MultiModel_Geometry_Variance_conf90.png",
        "  Type    : Per-tier SN2 distance/angle box + strip, restricted to confidence ≥ 0.90",
        "  Look for: Does the geometry–tier trend hold when only high-confidence models are kept?",
        "",
        "-" * 80,
        "Figure — 08_Diagnostic_and_MultiModel_Trends/15_MultiModel_Geometry_Variance_iptm90.png",
        "  Type    : As figure 14 but restricted to ipTM ≥ 0.90 (high interface confidence)",
        "  Look for: Whether the high-interface-confidence subset preserves the per-tier geometry trend.",
        "",
    ]
    _txt = "\n".join(lines)
    for _old_fp, _new_fp in _fig_paths.items():
        _txt = _txt.replace(_old_fp, _new_fp)

    """
    Folder 01's figures are named, not numbered, so the numbered blocks above never reach them and the
    Ramachandran panels — the only structural validation against a crystal in the whole pipeline — went
    undescribed. They are appended here rather than renumbered: the filenames are what the backbone
    validation refers to.
    """
    _txt += "\n" + "\n".join([
        "=" * 80,
        "BACKBONE GEOMETRY VALIDATION (subfolder: 02_Ramachandran/)",
        "=" * 80,
        "",
        "-" * 80,
        "02_Ramachandran/Ramachandran_3R3U_Crystal.png",
        "  Title   : Ramachandran plot — 3R3U crystal structure",
        "  Type    : phi/psi scatter over the favoured/allowed contour map",
        "  Look for: the reference. Every Boltz-2 prediction is judged against this distribution.",
        "",
        "-" * 80,
        "02_Ramachandran/Ramachandran_DeHa4_Control.png",
        "02_Ramachandran/Ramachandran_3R3U_Control.png",
        "  Title   : Ramachandran plot — Boltz-2 control prediction (DeHa4 / 3R3U sequence)",
        "  Look for: outliers in the disallowed regions. A prediction with good ipTM but a strained",
        "            backbone has bought its confidence with geometry the protein cannot hold.",
        "",
        "-" * 80,
        "02_Ramachandran/Ramachandran_DeHa4_vs_Crystal.png",
        "02_Ramachandran/Ramachandran_3R3U_vs_Crystal.png",
        "  Title   : Prediction overlaid on the crystal",
        "  Look for: whether the predicted backbone occupies the SAME basins as the crystal, not",
        "            merely allowed ones. This is the check that Boltz-2 folded the enzyme rather",
        "            than something plausible-looking of the same sequence.",
        "",
    ])

    """
    The legend must describe the figures the run actually produced. A figure can be skipped for a
    legitimate, data-driven reason — the Hidden Gems panel draws nothing when no complex is
    physics-strong yet AI-doubted — and a legend that still promises it sends the reader to a file
    that is not there. Blocks whose figure is absent are moved to a closing section that says so, and
    why, rather than being silently deleted: 'this figure has no data' is itself a finding.
    """
    _blocks = _txt.split("-" * 80)
    _kept, _absent = [], []
    for _b in _blocks:
        _pngs = re.findall(r"([0-9A-Za-z_]+/[0-9A-Za-z_]+\.png)", _b)
        if _pngs and not any((out_dir / _p).exists() for _p in _pngs):
            _absent.append((_pngs[0], _b))
        else:
            _kept.append(_b)
    _txt = ("-" * 80).join(_kept)
    if _absent:
        _txt += "\n" + "\n".join([
            "=" * 80,
            "NOT PRODUCED IN THIS RUN",
            "=" * 80,
            "",
            "These figures are defined by the pipeline but drew no data on this dataset. That is a",
            "result, not an omission — the reason is given per figure.",
            "",
        ])
        for _p, _b in _absent:
            _why = _ABSENT_FIGURE_REASON.get(Path(_p).name, "no rows matched this figure's criteria.")
            _txt += f"\n  · {_p}\n      {_why}\n"

    desc_path = _aux_dir(out_dir) / "05_Figure_Descriptions.txt"
    desc_path.write_text(_txt, encoding="utf-8")
    _stats_path = _write_statistical_tests(out_dir)
    if _stats_path is not None:
        # console_info (not bare print) so this line carries the same 2-space indent as every other
        # '✔ Saved:' line the reporter emits — otherwise it sits two spaces to the left of the rest.
        console_info(f"  ✔ Saved: 01_Analysis_Data/{_stats_path.name}")
    return desc_path


# =============================================================================
# SECTION 6: MAIN EXECUTION
# =============================================================================

def main():
    # Collapse stacked separator rules: a caller prints a rule, a helper prints its own,
    # and the log grows triple bars with nothing between them.
    _utils_mod.install_console_rule_filter()
    parser = argparse.ArgumentParser(
        description="Boltz-2 Master Validation & Dendrogram Framework",
        usage="%(prog)s <run_folder>  (e.g. Boltz-2_Run_20260309T085406Z)"
    )
    parser.add_argument("run", help="Name of the Boltz-2 run folder (e.g. Boltz-2_Run_20260309T085406Z)")
    parser.add_argument("--no-variance", action="store_true",
                        help="Skip building the per-model variance CSV. By default it is built once "
                             "(re-parsing the 5 model CIFs per complex) and cached in "
                             "01_Analysis_Data; later runs reuse it. Skipping it costs the geometry "
                             "figure its two inter-model uncertainty panels.")
    args = parser.parse_args()

    """
    The variance CSV is BUILT BY DEFAULT, once. It is the only artefact behind the uncertainty panels,
    and a figure that quietly omits half of itself on every run is worse than a run that pays the cost
    once: the build re-parses five CIFs per complex, but the result is cached in 01_Analysis_Data, so
    only the first run after a fresh production pays for it.
    """
    global _ALLOW_VARIANCE_COMPUTE
    _ALLOW_VARIANCE_COMPUTE = not bool(args.no_variance)

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

    # Each numbered figure folder is created lazily — only when its step first writes into it (see
    # _redirect_savefig, _panel, and the Step-2 Ramachandran / Step-8 diagnostic calls). Folders then
    # appear on disk in step order, never pre-created empty ahead of their turn.
    out_dir.mkdir(parents=True, exist_ok=True)
    rama_dir = out_dir / "02_Ramachandran"

    np.random.seed(int(CFG.ANALYSIS_SEED))

    reporter = _make_reporter(out_dir)

    try:
        # -------------------------------------------------------------------------------
        # Pipeline Execution Phase
        # -------------------------------------------------------------------------------
        df, features = load_and_prep_data(prod_dir, reporter)
        df = perform_advanced_ranking(df, features, out_dir, reporter)
        df = analyse_conflicts(df, out_dir, reporter)

        enriched_csv = _aux_dir(out_dir) / CFG.FILE_VALIDATED_MASTER
        # Atomic: this is the ~370 MB deliverable 04 reads with no row-count check, so a killed run must
        # not leave a truncated CSV consumed as truth. tmp + replace (as 05/06/07 do via this helper).
        _utils_mod.atomic_write_csv(df, enriched_csv, index=False)
        reporter.log(f"  ✔ Saved: 01_Analysis_Data/{CFG.FILE_VALIDATED_MASTER}")

        # Finalise 01_Analysis_Data before the figure steps: build the multi-model variance CSV here
        # (parse the 5 model CIFs per complex) so the whole 01_Analysis_Data set is complete before
        # Step 2. The geometry panel in Step 5 then only READS this cache — it never rebuilds.
        if _ALLOW_VARIANCE_COMPUTE:
            reporter.log("")
            reporter.log("  Multi-Model Variance QC")
            try:
                _xn__ensure_multimodel_variance_csv(prod_dir, out_dir, reporter)
            except Exception as _e:                               # noqa: BLE001
                reporter.log(f"  ! Multi-model variance skipped: {type(_e).__name__}: {_e}")

        reporter.log("")
        reporter.log("  Active-Site Residue Mapping QC")
        try:
            write_residue_mapping_missed_csv(out_dir, reporter)
        except Exception as _e:                                   # noqa: BLE001
            reporter.log(f"  ! Residue-mapping QC skipped: {type(_e).__name__}: {_e}")

        # Figures are generated folder-by-folder in narrative order (02 → 08).
        # 02_Ramachandran — control backbone-geometry validation.
        rama_dir.mkdir(parents=True, exist_ok=True)   # 02_Ramachandran created at its turn (Step 2)
        generate_ramachandran_figures(prod_dir, rama_dir, reporter)

        # 02–06 main validation suite, then 07 diagnostics; the two-criteria figure
        # lands in folder 04 (fig 10), all routed into their numbered folders from inside this call.
        generate_comprehensive_figures(df, features, out_dir, reporter)

        # Cleanup temporary structure thumbnails
        for p in sorted(out_dir.rglob("_tt_thumbnails")):
            if p.is_dir():
                shutil.rmtree(p)

        # The battery runs after every figure has registered its own tests, so the whole family —
        # the figures' tests and the battery's — is Benjamini-Hochberg corrected together, once.
        # It must run BEFORE write_figure_descriptions, which is what writes the stats CSV out.
        try:
            _statistical_battery(df, reporter)
        except Exception as _e:                               # noqa: BLE001
            reporter.log(f"  ! Statistical battery skipped: {type(_e).__name__}: {_e}")

        desc_path = write_figure_descriptions(out_dir)
        reporter.log(f"  ✔ Saved: 01_Analysis_Data/{desc_path.name}")

    except Exception as e:
        print(f"\n[CRITICAL ERROR] Pipeline Failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    _t0 = _time.perf_counter()
    main()
    _utils_mod.print_elapsed(_t0, "03_Validation_Figures_FAcDs.py")
