#!/usr/bin/env python3

"""
===============================================================================
FAcDs Pipeline  |  Step 03  |  Validation, Ranking & Visualisation
===============================================================================
The definitive "Judge": merges physics-based structural validation with
confidence metrics from Boltz-2, performs Pareto optimisation, rescues
hidden-gem candidates, and generates a 24-figure publication figure suite.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 June 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 03_Validation_Figures_FAcDs.py
  Role          : Scoring, ranking, and visual reporting of Boltz-2 predictions.
  Imports from  : 00_02_Project_Config_FAcDs.py  (CFG — tier colours, vis params)
                  00_02_Project_Utils_FAcDs.py   (ConsoleColours, setup_logging,
                                                 console_info, console_separator)
  Reads         : <Run>/1_Boltz2_Production/*_Ranked_*.csv  (falls back to *_Master_*.csv)
                  (Master CSV written by 02_Production_FAcDs.py; latest file selected)
  Writes        : <Run>/3_Validation_Figures/Figure_01_*.png … Figure_24_*.png
                  <Run>/3_Validation_Figures/03_Final_Validated_Master.csv
                  <Run>/3_Validation_Figures/04_ACTION_Rescue_Hidden_Gems.csv
                  <Run>/3_Validation_Figures/01_Analysis_Log.txt
  Upstream      : 02_Production_FAcDs.py → writes the master ranked CSV consumed here
  Downstream    : 04_Phylogeny_FAcDs.py  → reads 03_Final_Validated_Master.csv
───────────────────────────────────────────────────────────────────────────────

# ── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Memory overhead: Generating 24+ high-resolution figures simultaneously can
     spike RAM usage; recommended 32GB+ for large (>10k row) datasets.
  2. CSV Schema Sensitivity: Relies on the exact column naming convention from
     Step 02; custom CSV modifications will break the scoring logic.
  3. Pareto Convergence: Non-dominated sorting complexity scales O(N log N);
     runs with extremely high objective-conflict may take several minutes.
  4. Visualization Headless: Matplotlib MUST use the 'Agg' backend; the script
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

    Phylogenetic analysis is handled by the downstream 04_Phylogeny_FAcDs.py script.

-------------------------------------------------------------------------------
Outputs (Saved in <Run_Folder>/3_Validation_Figures/):
    [Data]
    • 03_Final_Validated_Master.csv  <-- THE FINAL DATASET
    • 04_ACTION_Rescue_Hidden_Gems.csv      <-- MANUAL REVIEW LIST
    • 01_Analysis_Log.txt                   <-- Detailed Execution Log

    [Figures — ordered by scientific narrative: basic → advanced]
    ── Part 1: Dataset Overview ──
    • Figure_01_Tier_Distribution.png                 <-- Catalytic tier counts + model selection pie
    • Figure_02_Alignment_Grades.png                  <-- Sequence identity grades + KDE by tier
    • Figure_03_Tier_Grade_Distribution.png           <-- Tier × alignment grade cross-tabulation
    ── Part 2: AI Prediction Quality ──
    • Figure_04a_AI_Quality_Assessment.png            <-- Boltz confidence boxes + pTM/ipTM panel
    • Figure_04b_PA_AI_Quality_Space.png              <-- Perfect_A pTM/ipTM density + thumbnails
    • Figure_05_pTM_vs_ipTM_by_Tier.png              <-- pTM vs ipTM 2-D scatter by tier
    ── Part 3: Structural Validation ──
    • Figure_06_ActiveSite_RMSD_by_Tier.png          <-- Active-site RMSD vs reference control
    • Figure_07_Feature_Correlations.png              <-- Spearman ρ feature correlation heatmap
    • Figure_08_Tier_Quality_DotPlot.png              <-- Multi-metric tier quality Cleveland dot plot
    ── Part 4: Mechanistic Analysis ──
    • Figure_09_Mech_State_CrossTab.png               <-- Halide stabilisation × carboxylate clamp heatmap
    • Figure_10_Mechanistic_Fingerprint_by_Tier.png   <-- Catalytic fingerprint score mean ± CI
    • Figure_11_SN2_Angle_by_Tier.png                 <-- SN2 attack angle ECDF by tier
    • Figure_12a_Mechanism_Geometry_Scatter.png       <-- SN2 angle vs nucleophile distance scatter
    • Figure_12b_PA_Mechanistic_Quality_Space.png     <-- Perfect_A mechanistic density + thumbnails
    ── Part 5: Ligand Interactions ──
    • Figure_13a_Molecular_Interaction_Profile.png    <-- Bond type profile + fluorine engagement
    • Figure_13b_PA_Interaction_Quality_Space.png     <-- Perfect_A interaction density + thumbnails
    • Figure_14_Fluorine_Engagement_by_Tier.png       <-- Fluorine engagement ratio box + trend line
    ── Part 6: Binding Energetics ──
    • Figure_15_Binding_Energetics.png                <-- Binding probability violin + product inhibition
    • Figure_16_Product_Inhibition_by_Tier.png        <-- Product inhibition penalty box + strip
    ── Part 7: Chemical Space ──
    • Figure_17a_Chemical_Space_Map.png               <-- UMAP chemical space manifold by tier
    • Figure_17b_PA_Chemical_Space_Landscape.png      <-- Perfect_A KDE density + structure thumbnails
    ── Part 8: Multi-metric Synthesis ──
    • Figure_18a_Fingerprint_TopHits.png              <-- Radar: top-5 hits vs worst-tier baseline
    • Figure_18b_Fingerprint_TierReps.png             <-- Radar: one representative per tier
    • Figure_19a_Tier_Success_Rates.png               <-- Tier % success: substrate vs inhibitor geometry
    • Figure_19b_Conf_SN2_Landscape.png               <-- Per-tier medians in Confidence × SN2-angle space
    • Figure_20_Conflict_Composition.png              <-- Conflict category × tier stacked bar
    • Figure_21_Hidden_Gems_DeepDive.png              <-- Physics-good / AI-missed slope chart
    • Figure_22_Category_Overlap_Euler.png            <-- Multi-set overlap Euler diagram
    ── Part 9: Publication Assembly ──
    • Figure_23_Top25_Multitarget_Proteins.png        <-- Top-25 multitarget proteins stacked bar
    • Figure_23b_Ligand_Network.png                   <-- PFAS radial ligand network (tier-coloured)
    • Figure_24_Sankey_Workflow.png                   <-- Sankey: nuc-dist → SN2 angle → final tier
    • Figure_25a_PFAS_Size_Hexbin_Landscape.png        <-- Chain-length hexbin + tier scatter + rolling median
    • Figure_25b_PFAS_Size_Composition_Merged.png      <-- Per bin: outcome (left) + degrader-tier (right) stacked bars
    • Figure_25c_PFAS_Size_Confidence_Boxplot.png       <-- Boltz confidence boxplot by chain-length bin


===============================================================================
"""

# ===============================================================================
# SECTION 1: SYSTEM CONFIGURATION & IMPORTS
# ===============================================================================

# -------------------------------------------------------------------------------
# Step 1.1: Standard Library Imports
# -------------------------------------------------------------------------------
import sys
import os
import shutil
import argparse
import warnings
from pathlib import Path
from datetime import datetime

# -------------------------------------------------------------------------------
# Step 1.2: Scientific Stack Imports
# -------------------------------------------------------------------------------
import numpy as np
import pandas as pd
import seaborn as sns
import umap
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, gaussian_kde
from sklearn.preprocessing import MinMaxScaler, RobustScaler
from sklearn.decomposition import PCA

# Bioinformatics Imports
from matplotlib.patches import ConnectionPatch
from matplotlib.lines import Line2D
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
_utils_mod = _load_module("ProjectUtils",  Path(__file__).resolve().parent / "00_02_Project_Utils_FAcDs.py")
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

# Tier palette — sourced from CFG.TIER_COLOUR (Okabe-Ito colourblind-safe).
# Local overrides: "Decoy" uses a softer grey for unlabelled entries;
# "Error" is a 03-specific indicator for analytics failures (not in CFG).
TIER_PALETTE = {
    **CFG.TIER_COLOUR,
    "Decoy":  "#999999",
    "Error": "#FF6B6B",
}

# Specific order for tiers to ensure logical plotting (Best to Worst)
TIER_ORDER_LOGIC = [
    "Perfect_A", "Perfect_B",
    "Best_A", "Best_B",
    "Good", "Poor", "Decoy"
]

CONFLICT_PALETTE = dict(CFG.CONFLICT_COLOUR)   # sourced from CFG § 9.7

# Modern Plotting Theme (Publication Quality)
sns.set_theme(style="whitegrid", context="paper", font_scale=1.3)
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], # Nature preference
    "axes.grid": True,
    "grid.color": "#F0F0F0",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "savefig.dpi": 300,
    "savefig.bbox": "tight"
})

SEPARATOR = "-" * 80


# ===============================================================================
# SECTION 2: LOGGING & UTILITIES
# ===============================================================================

logger = None

def console_title(msg: str) -> None:
    _console_title(msg, logger)

def console_info(msg: str) -> None:
    _console_info(msg, logger)

def console_separator() -> None:
    _console_sep(logger, heavy=True)

class ReportManager:
    """Manages writing simultaneous logs to console and file."""
    def __init__(self, out_dir: Path):
        self.path = out_dir / "01_Analysis_Log.txt"
        with open(self.path, "w") as f:
            f.write(f"BOLTZ-2 VALIDATION & PHYLOGENY REPORT \n")
            f.write(f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * 80 + "\n\n")

    def log(self, text: str):
        console_info(text)
        with open(self.path, "a") as f: f.write(f"[LOG] {text}\n")

    def section(self, title: str):
        print(f"\n{ConsoleColours.BOLD}{title}{ConsoleColours.ENDC}", flush=True)
        print(SEPARATOR_LIGHT, flush=True)
        with open(self.path, "a") as f: f.write(f"\n--- {title} ---\n")

def calculate_alignment_grade(identity):
    """Maps identity percentage to Grades A-I (matching Production logic)."""
    try:
        val = float(identity)
    except (ValueError, TypeError):
        return "Grade I (< 20%)"

    if val >= 90.0: return "Grade A (90-100%)"
    if val >= 80.0: return "Grade B (80-90%)"
    if val >= 70.0: return "Grade C (70-80%)"
    if val >= 60.0: return "Grade D (60-70%)"
    if val >= 50.0: return "Grade E (50-60%)"
    if val >= 40.0: return "Grade F (40-50%)"
    if val >= 30.0: return "Grade G (30-40%)"
    if val >= 20.0: return "Grade H (20-30%)"
    return "Grade I (< 20%)"


# ===============================================================================
# SECTION 3: ANALYSIS ALGORITHMS (PARETO OPTIMISATION)
# ===============================================================================

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
            
            min_y_seen = float('inf')
            current_front_indices = []
            
            for i in range(len(sorted_data)):
                y_val = sorted_data[i, 1]
                if y_val < min_y_seen:
                    min_y_seen = y_val
                    current_front_indices.append(sorted_ids[i])
                    
            pareto_ranks[current_front_indices] = current_rank
            mask_keep = np.isin(population_ids, current_front_indices, invert=True)
            population_ids = population_ids[mask_keep]

        else:
            # Vectorised dominance check — no Python loop.
            # boe[i,j,m]: subset[j,m] <= subset[i,m]  (j better-or-equal to i on obj m)
            # b[i,j,m]:   subset[j,m] <  subset[i,m]  (j strictly better on obj m)
            # j dominates i iff all-objectives boe AND any-objective b.
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

# ===============================================================================
# SECTION 4: DATA PROCESSING & FIGURES (PIPELINE)
# ===============================================================================

def load_and_prep_data(prod_dir: Path, reporter: ReportManager) -> tuple[pd.DataFrame, list[str]]:
    # Prefer ranked CSV (has Scientific_Rank + Ranking_Score_Calc used by downstream steps).
    # Fall back to master CSV if ranked not yet generated.
    candidates = sorted(list(prod_dir.glob("*_Ranked_*.csv")))
    if not candidates: raise FileNotFoundError("No Ranked CSV found in Production folder.")
        
    target = candidates[-1]
    reporter.log(f"Input Data Source: {target.resolve()}")
    df = pd.read_csv(target, low_memory=False)
    
    if 'job_name' in df.columns: 
        df = df.drop_duplicates(subset=['job_name'], keep='last')
    
    tier_map = CFG.TIER_RANK
    df['tier_numeric'] = df['degrader_tier'].map(tier_map).fillna(0)
    
    col_map = {
        "binding_likelihood_computed": "Binding_Probability",
        "Binding_Probability_Score": "Binding_Probability",
        "custom_affinity_score": "Chemical_Affinity_Score",
        "confidence_score": "Boltz_Model_Confidence",
        "interaction_density": "Interaction_Density_Norm",
        "sn2_attack_angle": "SN2_Attack_Angle",
        "dist_Nuc": "Dist_Nucleophile_ASP110",
        "dist_Stab_W": "Dist_Stabiliser_TRP156",
        "dist_Stab_Y": "Dist_Stabiliser_TYR219",
    }
    df.rename(columns=col_map, inplace=True)

    candidate_features = [
        "tier_numeric", "Boltz_Model_Confidence", "iptm", "mean_plddt",
        "Interaction_Density_Norm", "Chemical_Affinity_Score", "Binding_Probability",
        "count_hydrogen_bond", "count_salt_bridge"
    ]
    
    valid_features = []
    for f in candidate_features:
        if f in df.columns:
            df[f] = pd.to_numeric(df[f], errors='coerce').fillna(0)
            if df[f].std() > 0.001: valid_features.append(f)
            
    reporter.log(f"Loaded {len(df)} complexes. Valid Features for Analysis: {len(valid_features)}")
    return df, valid_features

def perform_advanced_ranking(df: pd.DataFrame, features: list[str], out_dir: Path, reporter: ReportManager):
    reporter.section("Step 1: Multi-Objective Ranking (PCA & Pareto)")
    
    x = df[features].dropna()
    x_scaled = RobustScaler().fit_transform(x)
    
    # Default in case PCA is skipped (len(x) < 2)
    score_pca = np.full(len(x), 50.0)

    # PCA
    if len(x) >= 2:
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(x_scaled)
        
        check_col = "tier_numeric" if "tier_numeric" in features else features[0]
        corr, _ = spearmanr(pcs[:, 0], df.loc[x.index, check_col])

        score_pca = MinMaxScaler(feature_range=(0, 100)).fit_transform(pcs[:, 0].reshape(-1, 1)).flatten()
        if corr < 0:
            score_pca = 100.0 - score_pca
            reporter.log(
                f"  -> PC1 direction: Spearman ρ(PC1, {check_col}) = {corr:.4f} — NEGATIVE."
                f"  PC1 scores mathematically inverted so high PC1 = high degradability tier."
                f"  [Methods note: PC1 sign chosen by Spearman correlation with tier_numeric; "
                f"  inversion applied when ρ < 0. Report this in supplementary methods.]"
            )
        else:
            reporter.log(
                f"  -> PC1 direction: Spearman ρ(PC1, {check_col}) = {corr:.4f} — POSITIVE."
                f"  PC1 preserved (high PC1 = high degradability tier, no inversion needed)."
            )
            
        df.loc[x.index, 'Score_PCA_Data'] = score_pca
    else:
        reporter.log("  ! Skipping PCA: Insufficient data points (< 2).")
        df.loc[x.index, 'Score_PCA_Data'] = 50.0 

    # Weighted Score
    score_weighted = np.zeros(len(x))
    total_w = 0.0
    for col, w in WEIGHTS.items():
        if col in df.columns:
            vals = df.loc[x.index, col].values.reshape(-1, 1)
            norm_col = MinMaxScaler().fit_transform(vals).flatten()
            score_weighted += norm_col * w
            total_w += w
    if total_w > 0: score_weighted = (score_weighted / total_w) * 100
    df.loc[x.index, 'Score_Weighted'] = score_weighted
    
    # Pareto
    if "Boltz_Model_Confidence" in df.columns and "Binding_Probability" in df.columns:
        reporter.log("  -> Calculating Pareto Frontiers (Confidence vs Binding Probability)...")
        df['Pareto_Rank'] = calculate_pareto_fronts(df, ["Boltz_Model_Confidence", "Binding_Probability"], [True, True])
    else:
        df['Pareto_Rank'] = 0 
        
    df.loc[x.index, 'Ensemble_Score'] = (0.5 * score_pca) + (0.5 * score_weighted)
    df['Ensemble_Score'] = df['Ensemble_Score'].fillna(0)
    df['Ensemble_Data_Rank'] = df['Ensemble_Score'].rank(ascending=False, method='min')
    
    # UMAP
    if len(x) >= 5:
        try:
            reporter.log("  -> Computing UMAP Manifold for Chemical Space Map...")
            reducer = umap.UMAP(n_neighbors=min(15, len(x)-1), min_dist=0.1, random_state=42)
            umap_map = reducer.fit_transform(x_scaled)
            df.loc[x.index, 'UMAP_X'] = umap_map[:, 0]
            df.loc[x.index, 'UMAP_Y'] = umap_map[:, 1]
        except Exception as e:
            reporter.log(f"  ! UMAP Calculation Skipped: {e}")
    else:
        reporter.log("  ! Skipping UMAP: Insufficient data points (< 5).")

    if len(x) >= 2:
        pd.DataFrame(pca.components_.T, columns=['PC1', 'PC2'], index=features).to_csv(out_dir / "02_PCA_Loadings.csv", index=True)
    
    return df

def analyse_conflicts(df: pd.DataFrame, out_dir: Path, reporter: ReportManager):
    reporter.section("Step 2: Conflict & Opportunity Analysis")
    
    def classify(row):
        tier = row.get('degrader_tier', 'Decoy')
        conf = row.get('Boltz_Model_Confidence', 0.0)

        _high_quality = [CFG.T_PA, CFG.T_PB, CFG.T_BA, CFG.T_BB]
        _hi = CFG.CONFLICT_CONF_HIGH
        _lo = CFG.CONFLICT_CONF_LOW
        if tier in _high_quality and conf >= _hi: return "Consensus High"
        if tier in [CFG.T_PR, CFG.T_DY, 'Error'] and conf < _lo: return "Consensus Low"
        if tier in _high_quality and conf < _hi: return "Hidden Gem"
        if tier in [CFG.T_PR, CFG.T_DY] and conf >= _hi: return "Decoy"
        return "Ambiguous"

    df['Conflict_Category'] = df.apply(classify, axis=1)
    
    gems = df[df['Conflict_Category'] == "Hidden Gem"].sort_values("Pareto_Rank")
    gems.to_csv(out_dir / "04_ACTION_Rescue_Hidden_Gems.csv", index=False)
    
    reporter.log(f"Hidden Gems (Rescue Target) : {len(gems)}")
    reporter.log(f"Decoys (Potential Artifacts): {len(df[df['Conflict_Category'] == 'Decoy'])}")
    
    return df

# ===============================================================================
# SECTION 4B: Perfect_A Landscape Companion Figures (05b, 09b, 10b, 15b)
# Companion figures overlaying Perfect_A structural thumbnails onto the same
# chemical spaces as their parent figures, showing where the Perfect_A hits
# sit relative to the full 58 k-complex dataset.
# ===============================================================================

# ── Design constants ──────────────────────────────────────────────────────────
_PA_ENTRY_COLS = ["#E74C3C", "#17A589", "#27AE60", "#8E44AD", "#E67E22"]
_PA_STAR_FILL  = "#2ECC71"
_PA_STAR_S     = CFG.VIS_PA_STAR_SIZE
_PA_TIER_S     = CFG.VIS_TIER_SIZES
_PA_TIER_A     = CFG.VIS_TIER_ALPHAS
_PA_STRUCT_COLS= ["tv_red","teal","forest","purple","tv_orange"]
_PA_SHRINK_B   = CFG.VIS_PA_SHRINK_BORDER
_PA_IMG_PX     = CFG.VIS_PA_IMG_PX
_PA_TEXT_PX    = CFG.VIS_PA_TEXT_PX
_PA_BORDER_PX  = CFG.VIS_PA_BORDER_PX
_PA_PAD        = CFG.VIS_PA_PAD
_PA_FONT_SIZE  = CFG.VIS_PA_FONT_SIZE
_PA_FONT_PATH  = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
_PA_AX_W       = CFG.VIS_PA_AX_WIDTH
_PA_LEFT_X0    = 0.058
_PA_RIGHT_X0   = 0.808
_PA_INSET_Y    = [0.60, 0.22, 0.69, 0.43, 0.17]
_PA_MARGINS    = dict(left=0.22, right=0.74, top=0.84, bottom=0.16)


def _pa_find_cif(pred_jobs: Path, job_name: str):
    jd = pred_jobs / job_name
    if not jd.exists():
        return None
    bc = jd / "Best_Complex"
    cifs = list(bc.glob("*.cif")) if bc.exists() else list(jd.rglob("*.cif"))
    return cifs[0] if cifs else None


def _pa_render_one(cif_path: Path, out_png: Path, prot_col: str, width=800, height=800) -> bool:
    try:
        import pymol
        pymol.finish_launching(["pymol", "-cq"])
        from pymol import cmd
        cmd.reinitialize()
        cmd.feedback('disable', 'all', 'everything')
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
        cmd.zoom("pocket", buffer=5)
        cmd.bg_color("black")
        cmd.viewport(width, height)
        cmd.set("ray_shadows",           "off")
        cmd.set("ambient",               0.6)
        cmd.set("ray_opaque_background", 1)
        cmd.png(str(out_png), width=width, height=height, ray=1, quiet=1)
        return out_png.exists()
    except Exception as e:
        return False


def _pa_render_all(pa: pd.DataFrame, pred_jobs: Path, thumb_dir: Path, reporter) -> list:
    thumb_dir.mkdir(exist_ok=True)
    imgs = []
    n_rendered = n_cached = n_failed = 0
    for local_i, (_, row) in enumerate(pa.iterrows()):
        png = thumb_dir / f"pa_{local_i}.png"
        if not png.exists():
            cif = _pa_find_cif(pred_jobs, str(row.get("job_name", "")))
            if cif:
                col = _PA_STRUCT_COLS[local_i % len(_PA_STRUCT_COLS)]
                ok = _pa_render_one(cif, png, col)
                if ok:
                    n_rendered += 1
                else:
                    n_failed += 1
            else:
                n_failed += 1
        else:
            n_cached += 1
        imgs.append(np.array(Image.open(png).convert("RGBA")) if png.exists() else None)
    return imgs


def _pa_make_composite(img_arr, label_str: str, colour_hex: str) -> np.ndarray:
    PAD = _PA_PAD
    br, bg, bb = (int(colour_hex[i:i+2], 16) for i in (1, 3, 5))
    W, H_s = _PA_IMG_PX, _PA_TEXT_PX
    src = Image.fromarray(img_arr).convert("RGBA")
    inner = _PA_IMG_PX - 2 * PAD
    src_small = src.resize((inner, inner), Image.LANCZOS)
    canvas = Image.new("RGBA", (W, W + H_s), (15, 15, 25, 255))
    canvas.paste(src_small, (PAD, PAD))
    draw = ImageDraw.Draw(canvas)
    draw.line([(  _PA_BORDER_PX, W), (W - _PA_BORDER_PX, W)],
              fill=(br, bg, bb, 200), width=4)
    usable_w = W - 2 * _PA_BORDER_PX - 8
    lines = label_str.split("\n")
    font = ImageFont.load_default()
    for fs in range(_PA_FONT_SIZE, 14, -2):
        try:
            f = ImageFont.truetype(_PA_FONT_PATH, fs)
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
    draw.rectangle([_PA_BORDER_PX // 2, _PA_BORDER_PX // 2,
                    W - _PA_BORDER_PX // 2 - 1, W + H_s - _PA_BORDER_PX // 2 - 1],
                   outline=(br, bg, bb, 255), width=_PA_BORDER_PX)
    return np.array(canvas)


def _pa_draw_scatter(ax, df):
    for tier in [t for t in reversed(CFG.TIER_ORDER) if t != CFG.T_PA]:
        sub = df[df["degrader_tier"] == tier]
        if sub.empty: continue
        ax.scatter(sub["_X"], sub["_Y"],
                   c=TIER_PALETTE.get(tier, "#BBB"),
                   s=_PA_TIER_S.get(tier, 4), alpha=_PA_TIER_A.get(tier, 0.2),
                   linewidths=0, rasterized=True, zorder=2)


def _pa_draw_stars(ax, pa):
    for local_i, (_, row) in enumerate(pa.iterrows()):
        ec   = _PA_ENTRY_COLS[local_i % len(_PA_ENTRY_COLS)]
        fill = "#F1C40F" if local_i == 4 else _PA_STAR_FILL
        sz   = int(_PA_STAR_S * 0.82) if local_i == 4 else _PA_STAR_S
        ax.scatter(row["_X"], row["_Y"], c=fill, s=sz, marker="*",
                   edgecolors=ec, linewidths=1.2, zorder=10)


def _pa_draw_thumbnails(fig, ax, pa, imgs):
    fw, fh = fig.get_size_inches()
    # Show ALL Perfect_A entries — distribute evenly between left and right columns
    n_max = min(len(pa), len(imgs))
    _ax_w = 0.085 if n_max > 5 else _PA_AX_W
    ax_h = _ax_w * (fw / fh) * (_PA_IMG_PX + _PA_TEXT_PX) / _PA_IMG_PX
    # Compute dynamic Y positions: split into left (even indices) and right (odd indices)
    # Each column gets ceil(n_max/2) entries; space them evenly within [0.10, 0.88]
    n_left  = (n_max + 1) // 2
    n_right = n_max - n_left
    _y_margin, _y_top = 0.10, 0.88
    def _ycols(n):
        if n == 0:
            return []
        step = (_y_top - _y_margin) / n
        return [_y_margin + i * step for i in range(n)]
    _y_left  = _ycols(n_left)
    _y_right = _ycols(n_right)
    for local_i, (_, row) in enumerate(pa.head(n_max).iterrows()):
        colour = _PA_ENTRY_COLS[local_i % len(_PA_ENTRY_COLS)]
        side   = "left" if local_i % 2 == 0 else "right"
        x0     = _PA_LEFT_X0 if side == "left" else _PA_RIGHT_X0
        col_idx = local_i // 2
        y0      = (_y_left[col_idx] if side == "left" else _y_right[col_idx])
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
        comp  = _pa_make_composite(img, label, colour)
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
            mutation_scale=22, shrinkA=4, shrinkB=_PA_SHRINK_B, zorder=11,
        ))


def _pa_legend_handles(df):
    n_pa = len(df[df["degrader_tier"] == "Perfect_A"])
    # Symbol-type header entries so the reader understands both glyphs
    h = [
        Line2D([0],[0], marker="o", color="w",
               markerfacecolor="#888888", markeredgecolor="none",
               markersize=7, alpha=0.55, label="Individual complex  (scatter dot)"),
        Line2D([0],[0], marker="*", color="w",
               markerfacecolor=_PA_STAR_FILL, markeredgecolor=_PA_ENTRY_COLS[0],
               markersize=14, label=f"Perfect_A ★ highlighted  (n={n_pa})"),
        Line2D([0],[0], marker="*", color="w",
               markerfacecolor="#F1C40F", markeredgecolor=_PA_ENTRY_COLS[4],
               markersize=11, label="Perfect_A ★ gold (TFA duplicate)"),
    ]
    for t in [t for t in CFG.TIER_ORDER if t != CFG.T_PA]:
        sub = df[df["degrader_tier"] == t]
        if sub.empty:
            continue
        h.append(Line2D([0],[0], marker="o", color="w",
                        markerfacecolor=TIER_PALETTE.get(t, "#BBB"),
                        markersize=7, alpha=0.85, label=f"{t}  (n={len(sub):,})"))
    return h


def _pa_add_legend(fig, handles):
    fig.legend(handles=handles, loc="lower center",
               bbox_to_anchor=(0.5, 0.02), ncol=4,
               fontsize=8.5, framealpha=0.95,
               title="Degrader tier  (★ = Perfect_A highlighted; ● = scatter background)",
               title_fontsize=9.0,
               borderpad=0.8, labelspacing=0.6, handletextpad=0.5)


def _pa_new_fig():
    fig, ax = plt.subplots(figsize=(18, 12))
    fig.subplots_adjust(**_PA_MARGINS)
    ax.set_facecolor("#FAFAFA")
    return fig, ax


def _pa_style(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel, fontsize=11, labelpad=6)
    ax.set_ylabel(ylabel, fontsize=11, labelpad=6)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=10)
    ax.tick_params(labelsize=9)
    ax.set_axisbelow(True)
    ax.grid(True, alpha=0.18, color="#AAAAAA", linewidth=0.5)


# ── Figure 17b — UMAP Chemical Space + PA Landscape ──────────────────────────
def _fig_17b_pa_landscape(df, pa, imgs, out_dir: Path, reporter):
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
        F  = np.clip(-0.596*np.log(np.clip(Zi/Zi.max(),1e-6,None)), 0, 15)
        fig, ax = _pa_new_fig()
        cf = ax.contourf(Xi,Yi,F, levels=45, cmap="Blues_r", vmin=0,vmax=15, alpha=0.90)
        ax.contour(Xi,Yi,F, levels=18, colors="white", linewidths=0.30, alpha=0.38)
        cb = fig.colorbar(cf, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("KDE Population Density  (darker = more complexes)", fontsize=9.5)
        cb.set_ticks([0,3,6,9,12,15])
        _pa_draw_scatter(ax, dv); _pa_draw_stars(ax, pax)
        _pa_draw_thumbnails(fig, ax, pax, imgs)
        _pa_style(ax, "UMAP 1", "UMAP 2", "")
        _pa_add_legend(fig, _pa_legend_handles(dv))
        ax.set_xlim(x0, x1); ax.set_ylim(y0, y1)
        out = out_dir / "Figure_17b_PA_Chemical_Space_Landscape.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        reporter.log(f"  ✔ Saved: {out.resolve()}")
    except Exception as e:
        reporter.log(f"  ! Figure 17b skipped: {e}")


# ── Figure 12b — SN2 Angle × Confidence + PA Landscape ───────────────────────
def _fig_12b_pa_mechanistic(df, pa, imgs, out_dir: Path, reporter):
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
        fig, ax = _pa_new_fig()
        hb = ax.hexbin(dv["_X"], dv["_Y"], gridsize=60, cmap="YlOrBr",
                       mincnt=1, linewidths=0.15, alpha=0.85, zorder=1,
                       extent=[dv["_X"].min(), dv["_X"].max(), ylo, yhi])
        cb = fig.colorbar(hb, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("Complex count per bin", fontsize=9.5)
        for angle, tlbl, tcol in [(175, "Perfect_A", _PA_STAR_FILL),
                                   (165, "Perfect_B", TIER_PALETTE["Perfect_B"]),
                                   (155, "Best_A",    TIER_PALETTE["Best_A"])]:
            ax.axvline(angle, color=tcol, ls="--", lw=1.6, alpha=0.85, zorder=4)
            ax.text(angle+0.3, yhi-0.002, f"≥{angle}°\n{tlbl}",
                    fontsize=5.5, color=tcol, fontweight="bold", va="top")
        ax.axhline(dv["_Y"].median(), color="#555", ls=":", lw=0.9, alpha=0.55, zorder=3)
        _pa_draw_scatter(ax, dv); _pa_draw_stars(ax, pax)
        _pa_draw_thumbnails(fig, ax, pax, imgs)
        _pa_style(ax,
                  "SN2 Attack Angle (°) — higher = near-ideal nucleophilic trajectory",
                  "Boltz Model Confidence (AI structural quality, higher = better)",
                  "")
        _pa_add_legend(fig, _pa_legend_handles(dv))
        ax.set_ylim(ylo, yhi)
        out = out_dir / "Figure_12b_PA_Mechanistic_Quality_Space.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        reporter.log(f"  ✔ Saved: {out.resolve()}")
    except Exception as e:
        reporter.log(f"  ! Figure 12b skipped: {e}")


# ── Figure 04b — Confidence × ipTM + PA Landscape ────────────────────────────
def _fig_04b_pa_ai_quality(df, pa, imgs, out_dir: Path, reporter):
    try:
        dv = df.dropna(subset=["Boltz_Model_Confidence","iptm"]).copy()
        dv["_X"] = dv["Boltz_Model_Confidence"]; dv["_Y"] = dv["iptm"]
        pax = pa.copy()
        pax["_X"] = pa["Boltz_Model_Confidence"]; pax["_Y"] = pa["iptm"]
        xlo,xhi = dv["_X"].quantile(0.01)-0.01, dv["_X"].max()+0.005
        ylo,yhi = dv["_Y"].quantile(0.01)-0.01, dv["_Y"].max()+0.005
        fig, ax = _pa_new_fig()
        hb = ax.hexbin(dv["_X"], dv["_Y"], gridsize=60, cmap="YlOrBr",
                       mincnt=1, linewidths=0.15, alpha=0.85, zorder=1,
                       extent=[xlo,xhi,ylo,yhi])
        cb = fig.colorbar(hb, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("Complex count per bin", fontsize=9.5)
        for v, fn in [(dv["_X"].median(), ax.axvline),(dv["_Y"].median(), ax.axhline)]:
            fn(v, color="#555", ls="--", lw=0.9, alpha=0.55, zorder=3)
        _pa_draw_scatter(ax, dv); _pa_draw_stars(ax, pax)
        _pa_draw_thumbnails(fig, ax, pax, imgs)
        _pa_style(ax,
                  "Boltz Model Confidence (higher = better folding quality)",
                  "ipTM — Interface Predicted TM-score",
                  "")
        _pa_add_legend(fig, _pa_legend_handles(dv))
        ax.set_xlim(xlo, xhi); ax.set_ylim(ylo, yhi)
        out = out_dir / "Figure_04b_PA_AI_Quality_Space.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        reporter.log(f"  ✔ Saved: {out.resolve()}")
    except Exception as e:
        reporter.log(f"  ! Figure 04b skipped: {e}")


# ── Figure 13b — Interaction Density × Count + PA Landscape ──────────────────
def _fig_13b_pa_interactions(df, pa, imgs, out_dir: Path, reporter):
    try:
        dv = df.dropna(subset=["Interaction_Density_Norm","num_interactions"]).copy()
        dv["_X"] = pd.to_numeric(dv["Interaction_Density_Norm"], errors="coerce")
        dv["_Y"] = pd.to_numeric(dv["num_interactions"],          errors="coerce")
        dv = dv.dropna(subset=["_X","_Y"])
        pax = pa.copy()
        pax["_X"] = pd.to_numeric(pa["Interaction_Density_Norm"], errors="coerce")
        pax["_Y"] = pd.to_numeric(pa["num_interactions"],          errors="coerce")
        xhi = max(dv["_X"].quantile(0.99)+0.5, pax["_X"].max()+0.5)
        yhi = max(dv["_Y"].quantile(0.99)+2,   pax["_Y"].max()+2)
        dv  = dv[(dv["_X"] <= xhi) & (dv["_Y"] <= yhi)]
        fig, ax = _pa_new_fig()
        hb = ax.hexbin(dv["_X"], dv["_Y"], gridsize=60, cmap="YlOrBr",
                       mincnt=1, linewidths=0.15, alpha=0.85, zorder=1)
        cb = fig.colorbar(hb, ax=ax, fraction=0.025, pad=0.01)
        cb.set_label("Complex count per bin", fontsize=9.5)
        xmed, ymed = dv["_X"].median(), dv["_Y"].median()
        ax.axvline(xmed, color="#555", ls="--", lw=0.9, alpha=0.55, zorder=3)
        ax.axhline(ymed, color="#555", ls="--", lw=0.9, alpha=0.55, zorder=3)
        ax.text(xmed+0.05, 0.5, f"median {xmed:.2f}", fontsize=7.5, color="#555")
        ax.text(0.1, ymed+0.3, f"median {ymed:.0f}",  fontsize=7.5, color="#555")
        _pa_draw_scatter(ax, dv); _pa_draw_stars(ax, pax)
        _pa_draw_thumbnails(fig, ax, pax, imgs)
        _pa_style(ax,
                  "Interaction Density Norm (interactions per ligand heavy atom)",
                  "Total Interactions (num_interactions, all contact types)",
                  "")
        _pa_add_legend(fig, _pa_legend_handles(dv))
        ax.set_xlim(0, xhi); ax.set_ylim(0, yhi)
        out = out_dir / "Figure_13b_PA_Interaction_Quality_Space.png"
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        reporter.log(f"  ✔ Saved: {out.resolve()}")
    except Exception as e:
        reporter.log(f"  ! Figure 13b skipped: {e}")


def generate_comprehensive_figures(df: pd.DataFrame, features: list[str], out_dir: Path, reporter: ReportManager):
    reporter.section("Step 3: Generating Scientific Figures (Publication Quality)")
    existing_tiers = [t for t in TIER_ORDER_LOGIC if t in df['degrader_tier'].unique()]

    # ── Global publication-quality rcParams ───────────────────────────────────
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'axes.titlesize': 13, 'axes.titleweight': 'bold',
        'axes.labelsize': 11,
        'xtick.labelsize': 9, 'ytick.labelsize': 9,
        'legend.fontsize': 8.5, 'legend.framealpha': 0.92,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True,
        'grid.color': '#EBEBEB', 'grid.linewidth': 0.6,
    })

    def _text_color(hex_bg: str, threshold: float = 0.45) -> str:
        """Return 'white' for dark backgrounds, near-black for light backgrounds."""
        try:
            r = int(hex_bg[1:3], 16) / 255
            g = int(hex_bg[3:5], 16) / 255
            b = int(hex_bg[5:7], 16) / 255
            luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
            return 'white' if luminance < threshold else '#1a1a1a'
        except Exception:
            return 'white'
    
    # --- Figure 01: Tier Distribution + Model Selection (pie inset) ---
    if 'degrader_tier' in df.columns:
        total_complexes = len(df)
        model_col = 'best_model_name' if 'best_model_name' in df.columns else None
        counts6 = [len(df[df['degrader_tier'] == t]) for t in existing_tiers]

        # Adaptive figure height: scale with the tallest bar so the pie doesn't float
        _max_cnt6 = max(counts6) if counts6 else 1
        _fig_h6 = max(5.5, min(10.0, 4.5 + _max_cnt6 / 12000))

        fig, ax = plt.subplots(figsize=(11, _fig_h6))
        bar_x6 = np.arange(len(existing_tiers))

        bars6 = ax.bar(bar_x6, counts6,
                       color=[TIER_PALETTE.get(t, '#999') for t in existing_tiers],
                       edgecolor='white', linewidth=0.8, width=0.65, zorder=2)
        for xi, cnt in zip(bar_x6, counts6):
            pct = cnt / total_complexes * 100 if total_complexes else 0
            _pct_str6 = (f'{pct:.3f}' if pct < 0.1 else
                         f'{pct:.2f}' if pct < 1.0 else
                         f'{pct:.1f}')
            ax.text(xi, cnt + _max_cnt6 * 0.004,
                    f'{cnt:,}\n({_pct_str6}%)',
                    ha='center', va='bottom', fontsize=9, fontweight='bold', zorder=5)

        # Y-axis ceiling: round up to the nearest 5000, then add 30% headroom for labels
        import math as _math6
        _y_ceil6 = _math6.ceil(_max_cnt6 / 5000) * 5000
        ax.set_ylim(0, _y_ceil6 * 1.28)
        ax.set_xticks(bar_x6)
        ax.set_xticklabels(existing_tiers, rotation=35, ha='right', fontsize=10)
        for tick, tier in zip(ax.get_xticklabels(), existing_tiers):
            tick.set_color(TIER_PALETTE.get(tier, 'black'))
            tick.set_fontweight('bold')
        ax.set_xlabel("Degrader Tier", fontsize=11)
        ax.set_ylabel("Number of complexes", fontsize=11)
        ax.yaxis.grid(True, linestyle=':', alpha=0.3, zorder=0)
        ax.set_axisbelow(True)

        # Pie inset: lower-left corner, bigger, with legend immediately to its right
        if model_col:
            model_counts = df[model_col].value_counts().sort_index()
            _pie_colours = ['#8ECFC9', '#A5C8E1', '#FABEBE', '#FFF4B8', '#C9E8C4',
                            '#F8C8A8', '#D5C8E8', '#B8D4E8'][:len(model_counts)]
            # Position: x=0.01, y=0.28 (lower — using the empty space above the bars)
            # Size: 0.34 wide × 0.50 tall — larger than before
            ax_pie = ax.inset_axes([0.01, 0.28, 0.34, 0.50])
            wedges, texts, autotexts = ax_pie.pie(
                model_counts.values,
                colors=_pie_colours,
                autopct='%1.0f%%',
                startangle=90,
                pctdistance=0.68,
                labeldistance=None,
                wedgeprops=dict(edgecolor='white', linewidth=1.6))
            for at in autotexts:
                at.set_fontsize(7.5)
                at.set_color('#111')
                at.set_fontweight('bold')
            ax_pie.set_frame_on(False)
            # Explicit model legend — placed immediately to the right of the pie inset
            _model_labels6 = []
            for raw_key in model_counts.index:
                raw_name = str(raw_key)
                if raw_name.lower().startswith('model_'):
                    raw_name = 'Model ' + raw_name.split('_', 1)[1]
                elif raw_name.lower().startswith('model'):
                    raw_name = 'Model ' + raw_name[5:].lstrip('_')
                _model_labels6.append(raw_name)
            from matplotlib.patches import Patch as _P06
            _pie_leg = [_P06(facecolor=_pie_colours[i], edgecolor='white',
                             label=_model_labels6[i]) for i in range(len(model_counts))]
            # Legend overlaps slightly with pie right edge to sit as close as possible.
            # Pie inset right edge = 0.01+0.34 = 0.35; anchor at 0.32 pushes legend
            # flush against the pie with no visible gap.
            ax.legend(handles=_pie_leg, loc='center left',
                      bbox_to_anchor=(0.32, 0.53),
                      fontsize=7.5, framealpha=0.90, fancybox=True,
                      title='Model selection', title_fontsize=7.5,
                      ncol=1, borderpad=0.25, labelspacing=0.20,
                      handlelength=0.8, handletextpad=0.35)

        plt.tight_layout()
        plt.savefig(out_dir / "Figure_01_Tier_Distribution.png", dpi=300, bbox_inches='tight')
        plt.close()
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_01_Tier_Distribution.png').resolve()}")

    # --- Figure 02: Sequence Identity — Histogram (grade bands) overlaid with tier KDE lines ---
    # Single chart: histogram bars show the overall grade distribution (how many proteins
    # per identity band); KDE lines per tier overlay on the same x-axis, revealing whether
    # different catalytic tiers cluster at different identity levels.
    if 'identity_pct' in df.columns:
        # Use existing single-letter Alignment_Grade from CSV; only compute if missing
        if 'Alignment_Grade' not in df.columns or df['Alignment_Grade'].isnull().all():
            _grade_map = {(90,100):'A',(80,90):'B',(70,80):'C',(60,70):'D',(50,60):'E',
                          (40,50):'F',(30,40):'G',(20,30):'H',(0,20):'I'}
            def _simple_grade(v):
                try:
                    v = float(v)
                except: return 'I'
                for (lo, hi), g in _grade_map.items():
                    if v >= lo: return g
                return 'I'
            df['Alignment_Grade'] = df['identity_pct'].apply(_simple_grade)
        # Ensure single-letter format (convert "Grade A (...)" → 'A' if needed)
        _ag = df['Alignment_Grade'].astype(str)
        if _ag.str.startswith('Grade ').any():
            df['Alignment_Grade'] = _ag.str.extract(r'Grade\s+([A-I])', expand=False).fillna('I')
        id_plot = df.dropna(subset=['identity_pct']).copy()
        id_plot['identity_pct'] = pd.to_numeric(id_plot['identity_pct'], errors='coerce').clip(0, 100)
        id_plot = id_plot.dropna(subset=['identity_pct'])

        fig, ax = plt.subplots(figsize=(13, 7))
        ax2 = ax.twinx()   # secondary y-axis for KDE density

        # Grade band background shading — sourced from CFG § 9.8
        grade_bands = list(CFG.GRADE_BANDS)
        # Per-protein grade counts (one row per unique protein)
        _prot_dedup02 = id_plot.drop_duplicates(subset=['Protein_Name']) if 'Protein_Name' in id_plot.columns else id_plot
        _n_total_prot02 = len(_prot_dedup02)
        for lo, hi, col, lbl in grade_bands:
            ax.axvspan(lo, hi, alpha=0.20, color=col, zorder=0)
            ax.text((lo + hi) / 2, 0.972, f'Grade {lbl}', ha='center', va='top',
                    fontsize=7.5, color=col, fontweight='bold',
                    transform=ax.get_xaxis_transform(),
                    bbox=dict(boxstyle='round,pad=0.08', fc='white', ec='none',
                              alpha=0.70), zorder=11)
            _cnt_g02 = int((_prot_dedup02['Alignment_Grade'] == lbl).sum())
            ax.text((lo + hi) / 2, 0.952, f'{_cnt_g02}/{_n_total_prot02}', ha='center', va='top',
                    fontsize=5.0, color=col,
                    transform=ax.get_xaxis_transform(), zorder=11)
        ax.set_axisbelow(True)

        # Primary: stacked histogram of proteins by tier (all tiers stacked per bin)
        bin_edges = np.arange(0, 101, 5)  # 5-pct-wide bins
        tier_hist_data = []
        tier_hist_labels = []
        tier_hist_colors = []
        for tier in existing_tiers:
            sub_t = id_plot[id_plot['degrader_tier'] == tier]
            if sub_t.empty:
                continue
            tier_hist_data.append(sub_t['identity_pct'].values)
            tier_hist_labels.append(f'{tier}  (n={len(sub_t):,})')
            tier_hist_colors.append(TIER_PALETTE.get(tier, '#999'))
        if tier_hist_data:
            ax.hist(tier_hist_data, bins=bin_edges, stacked=True,
                    color=tier_hist_colors, label=tier_hist_labels,
                    edgecolor='white', linewidth=0.4, alpha=0.75, zorder=2)

        # Secondary: KDE density line per tier — ridge (joyplot) approach with vertical offsets
        # Each tier is staggered upward so curves do not overlap (ridge/joyplot style).
        from scipy.stats import gaussian_kde as _gkde
        _kde_tiers_valid = [t for t in existing_tiers
                            if len(id_plot[id_plot['degrader_tier'] == t]) >= 10]
        _ridge_step = 0.008   # vertical offset per tier index
        for _tier_idx, tier in enumerate(_kde_tiers_valid):
            sub_t = id_plot[id_plot['degrader_tier'] == tier]
            kde_x = np.linspace(0, 100, 400)
            kde_y = _gkde(sub_t['identity_pct'])(kde_x)
            _tier_offset = _tier_idx * _ridge_step
            kde_y_shifted = kde_y + _tier_offset
            _tcol = TIER_PALETTE.get(tier, '#999')
            ax2.plot(kde_x, kde_y_shifted, color=_tcol,
                     linewidth=2.5, alpha=0.75, zorder=5 + _tier_idx,
                     linestyle='-')
            # Fill under each shifted KDE curve for ridge visual separation
            ax2.fill_between(kde_x, _tier_offset, kde_y_shifted,
                             color=_tcol, alpha=0.15, zorder=4 + _tier_idx)
            # Median tick mark on the KDE line
            med_x = float(sub_t['identity_pct'].median())
            med_y = float(np.atleast_1d(_gkde(sub_t['identity_pct'])(np.array([med_x])))[0])
            ax2.scatter([med_x], [med_y + _tier_offset], color=_tcol,
                        s=55, zorder=6 + _tier_idx, edgecolors='black', linewidths=0.8)

        # Anchor KDE axis at zero; show full density range with 10% headroom
        _kde_ymax = ax2.get_ylim()[1]
        ax2.set_ylim(0, _kde_ymax * 1.10)

        # Grade threshold vertical lines
        for pct_v, lbl_v, col_v in [(90, 'Grade A boundary', '#1B7837'),
                                      (60, 'Grade D boundary', '#E69F00'),
                                      (40, 'Grade F boundary', '#D55E00')]:
            ax.axvline(x=pct_v, color=col_v, linestyle='--', alpha=0.65,
                       linewidth=1.3, zorder=4)

        ax.set_xlim(0, 100)
        ax.set_xticks(range(0, 101, 10))
        ax.tick_params(axis='x', labelsize=9)
        # Left axis (count): steel-blue tick labels matching the grey count gridlines
        ax.tick_params(axis='y', labelsize=9, labelcolor='#2C6FAC')
        ax.yaxis.label.set_color('#2C6FAC')
        # Right axis (KDE density): amber tick labels matching the amber density gridlines
        ax2.tick_params(axis='y', labelsize=9, labelcolor='#A06000')
        ax2.yaxis.label.set_color('#A06000')
        ax.set_xlabel("Sequence Identity to Reference DeHa4 Control  (%)", fontsize=11)
        ax.set_ylabel("Number of Protein–Ligand Complexes  (stacked by tier)", fontsize=11)
        ax2.set_ylabel("KDE Density  (probability density)", fontsize=9, rotation=270, labelpad=14)

        # Dual-colour grid system:
        #   Left-axis (count) → steel-blue lines to match left tick labels
        #   Right-axis (KDE density) → amber lines to match right tick labels
        # Grid strictly behind everything: draw grids at zorder=0, bars at zorder=2+
        # ax.set_axisbelow only affects the primary axis; force zorder on both axes
        ax.set_axisbelow(True)
        ax2.set_axisbelow(True)
        # Raise ax2 (KDE lines) above ax (bars) so KDE lines render on top
        ax2.set_zorder(ax.get_zorder() + 1)
        ax2.patch.set_visible(False)   # keep ax background visible
        ax.yaxis.grid(True, color='#2C6FAC', linewidth=0.55, linestyle='-', alpha=0.25, zorder=0)
        ax.xaxis.grid(False)
        ax2.yaxis.grid(True, color='#A06000', linewidth=0.55, linestyle='--', alpha=0.20, zorder=0)

        _ymax02 = ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 200
        ax.annotate("Tiers converge\nnear 60% identity",
                    xy=(60, _ymax02 * 0.45),
                    xytext=(70, _ymax02 * 0.60),
                    fontsize=5.0, color='#444', style='italic',
                    arrowprops=dict(arrowstyle='->', color='#888', lw=0.9))
        # Perfect_A visibility annotation — arrow points to peak bar of the bottom stack layer
        _pa_sub02 = id_plot[id_plot['degrader_tier'] == 'Perfect_A']
        if len(_pa_sub02) > 0:
            _pa_hist02, _ = np.histogram(_pa_sub02['identity_pct'], bins=bin_edges)
            _pa_peak_idx = int(np.argmax(_pa_hist02))
            _pa_peak_x   = float(bin_edges[_pa_peak_idx]) + 2.5
            _pa_peak_cnt = int(_pa_hist02[_pa_peak_idx])
            _pa_med_x    = float(_pa_sub02['identity_pct'].median())
            ax.annotate(
                f'Perfect_A  (n={len(_pa_sub02):,})\nmedian identity {_pa_med_x:.0f}%',
                xy=(_pa_peak_x, max(_pa_peak_cnt, 1)),
                xytext=(min(_pa_peak_x + 12, 82), _ymax02 * 0.32),
                fontsize=5.0,
                color=TIER_PALETTE.get('Perfect_A', '#2E7D52'),
                fontweight='bold',
                arrowprops=dict(arrowstyle='->', lw=1.3,
                                color=TIER_PALETTE.get('Perfect_A', '#2E7D52')),
                bbox=dict(boxstyle='round,pad=0.22', fc='white',
                          ec=TIER_PALETTE.get('Perfect_A', '#2E7D52'),
                          alpha=0.90, linewidth=1.0),
                zorder=15
            )
            pass  # Fig02 PA diag suppressed
        # Single unified legend: compound handles (bar patch + KDE line with median dot)
        from matplotlib.lines import Line2D as _L7
        import matplotlib.patches as _mp02
        _bar_hdls_02, _bar_lbls_02 = ax.get_legend_handles_labels()
        _unified_handles = []
        _unified_labels  = []
        for _bh, _bl in zip(_bar_hdls_02, _bar_lbls_02):
            # Extract tier name from label like "Perfect_A  (n=8)"
            _tname = _bl.split('  ')[0] if '  ' in _bl else _bl
            _tcol  = TIER_PALETTE.get(_tname, '#999')
            # Compound: bar patch (facecolor) + dot marker
            _unified_handles.append(_bh)
            _kdot = _L7([0],[0], marker='o', color=_tcol, linewidth=1.5,
                        markerfacecolor=_tcol, markeredgecolor='black',
                        markeredgewidth=0.5, markersize=5)
            _unified_handles.append(_kdot)
            _unified_labels.append(_bl)
            _unified_labels.append('')  # blank label for the dot entry pairs
        # Merge pairs: interleave bar and dot into the same rows (ncol=2 pairs)
        # Simpler: build single-entry handles with dual marker using HandlerTuple
        try:
            from matplotlib.legend_handler import HandlerTuple as _HT02
            _pair_handles = []
            _pair_labels  = []
            for _bh, _bl in zip(_bar_hdls_02, _bar_lbls_02):
                _tname = _bl.split('  ')[0] if '  ' in _bl else _bl
                _tcol  = TIER_PALETTE.get(_tname, '#999')
                _kdot  = _L7([0],[0], marker='o', color=_tcol, linewidth=1.5,
                             markerfacecolor=_tcol, markeredgecolor='black',
                             markeredgewidth=0.5, markersize=5)
                _pair_handles.append((_bh, _kdot))
                _pair_labels.append(_tname)
            _leg7 = ax.legend(_pair_handles, _pair_labels,
                              loc='upper left', ncol=3, fontsize=8.5,
                              framealpha=0.92, fancybox=True,
                              handler_map={tuple: _HT02(ndivide=None, pad=0.2)},
                              title='Tier (bars = counts, KDE lines ● median)',
                              title_fontsize=8.0,
                              bbox_to_anchor=(0.01, 0.93))
        except Exception:
            _leg7 = ax.legend(_bar_hdls_02, _bar_lbls_02,
                              loc='upper left', ncol=3, fontsize=8.5,
                              framealpha=0.92, fancybox=True,
                              title='Tier (bars = counts, KDE lines ● median)',
                              title_fontsize=8.0,
                              bbox_to_anchor=(0.01, 0.93))
        _leg7.set_zorder(20)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_02_Alignment_Grades.png", dpi=300, bbox_inches='tight')
        plt.close()
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_02_Alignment_Grades.png').resolve()}")

    # --- Figure 03: Tier × Alignment Grade — Stacked 100% Horizontal Bar Chart ---
    # Each tier is a horizontal bar whose width = 100%.  The bar is divided into grade
    # bands (I worst → A best) using the same green-to-brown colour ramp as Fig 07,
    # making cross-figure reading intuitive.  Absolute counts are annotated inside wide
    # segments; narrow segments get small external labels.
    if 'identity_pct' in df.columns and 'degrader_tier' in df.columns:
        try:
            f11_df = df.dropna(subset=['identity_pct']).copy()
            f11_df['identity_pct'] = pd.to_numeric(f11_df['identity_pct'], errors='coerce').clip(0, 100)
            f11_df = f11_df.dropna(subset=['identity_pct'])
            grade_cuts = [0, 20, 30, 40, 50, 60, 70, 80, 90, 100]
            grade_lbls = ['I (<20%)', 'H (20–30%)', 'G (30–40%)', 'F (40–50%)',
                          'E (50–60%)', 'D (60–70%)', 'C (70–80%)', 'B (80–90%)', 'A (≥90%)']
            # Grade colours — sourced from CFG § 9.8
            grade_colors = dict(CFG.GRADE_COLOUR_FULL)
            f11_df['Grade'] = pd.cut(f11_df['identity_pct'], bins=grade_cuts,
                                     labels=grade_lbls, right=False)
            ct11 = pd.crosstab(f11_df['degrader_tier'], f11_df['Grade'])
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
                        color=grade_colors.get(grade, '#999'),
                        edgecolor='white', linewidth=0.6, zorder=2)
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
                        txt = f'{vi:.0f}%\n(n={int(ni):,})'
                        fs = 6.5
                    else:
                        txt = f'{vi:.0f}%'
                        fs = 6
                    txt_col = 'white' if grade_colors.get(grade, '#999')[1:] < 'AAAAAA' else '#111'
                    ax.text(mid_x, yi, txt, ha='center', va='center',
                            fontsize=fs, color=txt_col, fontweight='bold', zorder=3,
                            linespacing=1.1)
                left11 = left11 + vals

            # ── Grade bands identical to Fig02 — fixed identity-% positions ──
            # grade_bands = [(lo, hi, col, lbl), ...] with same positions as Fig02.
            # Both figures share a 0–100 x-axis so the visual layout matches exactly.
            _g11_lbl_to_full = {g.split(' ')[0]: g for g in grade_lbls}
            _g11_total       = ct11.sum(axis=0)
            _n_total11       = int(ct11.values.sum())
            for lo, hi, col, lbl in grade_bands:
                # Background shading — same alpha/style as Fig02
                ax.axvspan(lo, hi, alpha=0.20, color=col, zorder=0)
                # Dashed vertical divider at left boundary (skip x=0)
                if lo > 0:
                    ax.axvline(lo, color=col, lw=1.0, ls='--', alpha=0.55, zorder=3)
                # Grade header label at midpoint — same style as Fig02
                _mid = (lo + hi) / 2
                ax.text(_mid, 1.015, f'Grade {lbl}',
                        ha='center', va='bottom', fontsize=7.5, color=col,
                        fontweight='bold', transform=ax.get_xaxis_transform(),
                        clip_on=False,
                        bbox=dict(boxstyle='round,pad=0.08', fc='white',
                                  ec='none', alpha=0.70), zorder=11)
                _full_lbl = _g11_lbl_to_full.get(lbl, '')
                _gn = int(_g11_total.get(_full_lbl, 0))
                ax.text(_mid, 1.001, f'{_gn:,}/{_n_total11:,}',
                        ha='center', va='bottom', fontsize=5.0, color=col,
                        transform=ax.get_xaxis_transform(), clip_on=False, zorder=11)

            # Tier row labels
            _tier_n11 = ct11.sum(axis=1)
            ax.set_yticks(bar_y)
            ax.set_yticklabels([f'{t}  (n={_tier_n11.loc[t]:,})' for t in tiers_f11], fontsize=9.5)
            for tick, tier in zip(ax.get_yticklabels(), tiers_f11):
                tick.set_color(TIER_PALETTE.get(tier, 'black'))
                tick.set_fontweight('bold')

            # X-axis ticks every 10% — identical to Fig02
            ax.set_xlim(0, 100)
            ax.set_xticks(range(0, 101, 10))
            ax.set_xticklabels([f'{v}%' for v in range(0, 101, 10)],
                               fontsize=8, rotation=45, ha='right')
            ax.set_xlabel("Proportion of complexes in each sequence identity grade  (%)",
                          fontsize=11)
            ax.set_ylabel("Degrader Tier", fontsize=11)
            ax.xaxis.grid(True, color='#DCDCDC', linewidth=0.65, zorder=0)
            ax.set_axisbelow(True)

            if _chi2_str:
                ax.text(0.99, 0.01, _chi2_str, transform=ax.transAxes,
                        ha='right', va='bottom', fontsize=8.5, style='italic',
                        color='#333',
                        bbox=dict(boxstyle='round,pad=0.25', fc='white',
                                  ec='#ccc', alpha=0.88, linewidth=0.6))
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_03_Tier_Grade_Distribution.png", dpi=300, bbox_inches='tight')
            plt.close()
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_03_Tier_Grade_Distribution.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 03 skipped: {e}")

    # ── Perfect_A companion setup (shared by Figs 04b, 12b, 13b, 17b) ──────────
    _pa   = (df[df["degrader_tier"] == "Perfect_A"].reset_index(drop=True)
             if "degrader_tier" in df.columns else pd.DataFrame())
    _imgs = []
    try:
        if not _pa.empty:
            _pred_jobs = out_dir.parent / "1_Boltz2_Production" / "4_Prediction_Jobs"
            _thumb_dir = out_dir / "_pa_thumbnails"
            _thumb_dir.mkdir(parents=True, exist_ok=True)
            _imgs = _pa_render_all(_pa, _pred_jobs, _thumb_dir, reporter)
    except Exception:
        pass
    _pa_has_imgs = not _pa.empty and sum(i is not None for i in _imgs) > 0

    # --- Figure 04: AI Quality — Boltz Confidence (boxes) + pTM & ipTM (lines) overlaid ---
    # Single chart: box plots show per-tier distribution of Boltz Model Confidence (primary
    # y-axis). Overlaid connected dot-lines show per-tier median pTM and ipTM on the same
    # scale (both 0–1), making agreement and divergence between the three AI metrics visible.
    _has_conf  = 'Boltz_Model_Confidence' in df.columns and 'degrader_tier' in df.columns
    _has_ptm   = 'ptm' in df.columns and 'iptm' in df.columns and 'degrader_tier' in df.columns
    if _has_conf or _has_ptm:
        try:
            from scipy import stats as _sc_stats10
            fig, ax = plt.subplots(figsize=(13, 7))

            # Quality zone bands — distinct fills and matching boundary lines
            ax.axhspan(0.90, 1.01, alpha=0.10, color='#009E73', zorder=0)
            ax.axhspan(0.80, 0.90, alpha=0.09, color='#E69F00', zorder=0)
            ax.axhspan(0.00, 0.80, alpha=0.07, color='#D55E00', zorder=0)
            ax.axhline(y=0.90, color='#009E73', linestyle='--', linewidth=1.3, alpha=0.80, zorder=1)
            ax.axhline(y=0.80, color='#E69F00', linestyle=':', linewidth=1.3, alpha=0.80, zorder=1)
            # Zone labels — placed in DATA y-coordinates (get_yaxis_transform) AFTER ylim is set,
            # so they land at the correct zone midpoints regardless of the dynamic y-axis range.
            # (Labels are defined here as closures and drawn below after ax.set_ylim is called.)

            tier_x = {tier: i for i, tier in enumerate(existing_tiers)}

            # Boltz Confidence — box plots (primary, fills the background)
            if _has_conf:
                sns.boxplot(data=df, x='degrader_tier', y='Boltz_Model_Confidence',
                            order=existing_tiers, palette=TIER_PALETTE, linewidth=1.2,
                            flierprops=dict(marker='o', markersize=2, alpha=0.25),
                            width=0.55, zorder=2, ax=ax)
                # Median text above each box
                tier_order_f10 = list(existing_tiers)
                medians_f10 = []
                for i, tier in enumerate(tier_order_f10):
                    med = df.loc[df['degrader_tier'] == tier, 'Boltz_Model_Confidence'].median()
                    medians_f10.append(float(med) if not np.isnan(med) else float('nan'))
                    if not np.isnan(med):
                        ax.text(i, med + 0.003, f'{med:.3f}', ha='center', va='bottom',
                                fontsize=8, color='black', fontweight='bold', zorder=7)
                # Decoy overconfidence annotation
                _decoy_vals = [v for t, v in zip(tier_order_f10, medians_f10) if t == 'Decoy']
                if _decoy_vals and _decoy_vals[0] > 0.95:
                    _di = tier_order_f10.index('Decoy')
                    ax.annotate("Decoy\noverconfidence",
                                xy=(_di, _decoy_vals[0]),
                                xytext=(_di + 0.4, _decoy_vals[0] - 0.02),
                                fontsize=7.5, color='#B22222', style='italic',
                                arrowprops=dict(arrowstyle='->', color='#B22222', lw=0.9),
                                bbox=dict(boxstyle='round,pad=0.2', fc='#FFF5F5',
                                          ec='#B22222', alpha=0.9, lw=0.8))

            # pTM — connected median line with CI band
            if _has_ptm:
                ptm_xs, ptm_meds, ptm_lo, ptm_hi = [], [], [], []
                iptm_xs, iptm_meds, iptm_lo, iptm_hi = [], [], [], []
                for tier in existing_tiers:
                    sub = df.loc[df['degrader_tier'] == tier].dropna(subset=['ptm', 'iptm'])
                    if len(sub) < 3:
                        continue
                    x_pos = tier_x[tier]
                    for metric, xs_l, meds_l, lo_l, hi_l, col in [
                        ('ptm', ptm_xs, ptm_meds, ptm_lo, ptm_hi, '#0072B2'),
                        ('iptm', iptm_xs, iptm_meds, iptm_lo, iptm_hi, '#CC79A7'),
                    ]:
                        vals = sub[metric].dropna()
                        med_v = float(vals.median())
                        ci95  = float(_sc_stats10.sem(vals) * _sc_stats10.t.ppf(0.975, len(vals) - 1))
                        xs_l.append(x_pos)
                        meds_l.append(med_v)
                        lo_l.append(med_v - ci95)
                        hi_l.append(med_v + ci95)

                for xs_l, meds_l, lo_l, hi_l, col, lbl, mk in [
                    (ptm_xs, ptm_meds, ptm_lo, ptm_hi, '#0072B2', 'pTM  (fold confidence)', 's'),
                    (iptm_xs, iptm_meds, iptm_lo, iptm_hi, '#CC79A7', 'ipTM  (interface confidence)', '^'),
                ]:
                    if not xs_l:
                        continue
                    ax.plot(xs_l, meds_l, color=col, linewidth=2.0, zorder=5,
                            marker=mk, markersize=8, markeredgecolor='black',
                            markeredgewidth=0.8, label=lbl)
                    ax.fill_between(xs_l, lo_l, hi_l, color=col, alpha=0.18, zorder=4)

            # Y-axis zoom: start just below data minimum
            all_vals_f10 = []
            if _has_conf:
                all_vals_f10.extend(df['Boltz_Model_Confidence'].dropna().tolist())
            if _has_ptm:
                all_vals_f10.extend(df['ptm'].dropna().tolist())
                all_vals_f10.extend(df['iptm'].dropna().tolist())
            if all_vals_f10:
                y_lo_f10 = max(0.0, float(np.percentile(all_vals_f10, 1)) - 0.02)
                ax.set_ylim(y_lo_f10, 1.01)
            else:
                y_lo_f10 = 0.0
            # Broken-axis indicator at y_lo_f10 (axis is truncated — not starting at 0)
            if y_lo_f10 > 0.01:
                _bax_kw = dict(transform=ax.transAxes, color='#333', clip_on=False,
                               linewidth=2.2)
                _d = 0.022   # diagonal tick half-length in axes coords
                _bax_y = 0.0  # bottom of axes
                # Left y-axis break marks
                ax.plot((-_d, +_d), (_bax_y - _d, _bax_y + _d), **_bax_kw)
                ax.plot((-_d * 2.2, +_d * 0.2), (_bax_y - _d, _bax_y + _d), **_bax_kw)
                # Right side break marks
                ax.plot((1 - _d, 1 + _d), (_bax_y - _d, _bax_y + _d), **_bax_kw)
                ax.plot((1 - _d * 2.2, 1 + _d * 0.2), (_bax_y - _d, _bax_y + _d), **_bax_kw)
                ax.text(0.01, 0.015, "↕ axis break (starts at 0.80)",
                        transform=ax.transAxes,
                        fontsize=7.5, color='#555', ha='left', va='bottom',
                        style='italic')
            # Zone labels in DATA coordinates — each label at the midpoint of its
            # VISIBLE zone slice (zones clipped to the actual ylim after set_ylim).
            # get_yaxis_transform(): x=axes fraction [0,1], y=data coordinates
            _yt10 = ax.get_yaxis_transform()
            _y_top10 = 1.01
            # High confidence zone: 0.90–1.01 — fixed at top-right via transAxes so it is
            # always visible regardless of ylim and never lands on Perfect_A box bodies.
            _hc_vis_lo10 = max(0.90, y_lo_f10)
            if _hc_vis_lo10 < _y_top10 - 0.005:
                ax.text(0.01, 0.98, 'High confidence (≥0.90)', color='#007A50',
                        fontsize=8, ha='left', va='top', fontweight='bold',
                        style='italic', transform=ax.transAxes, zorder=6,
                        bbox=dict(boxstyle='round,pad=0.15', fc='#D6F5EB', ec='#009E73',
                                  alpha=0.85, linewidth=0.6))
            # Acceptable zone: 0.80–0.90
            _acc_vis_lo10 = max(0.80, y_lo_f10)
            _acc_vis_hi10 = min(0.90, _y_top10)
            if _acc_vis_hi10 > _acc_vis_lo10 + 0.005:
                _acc_mid10 = (_acc_vis_lo10 + _acc_vis_hi10) / 2
                ax.text(0.01, _acc_mid10, 'Acceptable (0.80–0.90)', color='#8A6000',
                        fontsize=8, ha='left', va='center', fontweight='bold',
                        style='italic', transform=_yt10, zorder=6,
                        bbox=dict(boxstyle='round,pad=0.15', fc='#FFF3CC', ec='#E69F00',
                                  alpha=0.85, linewidth=0.6))
            # Below threshold zone: y_lo–0.80
            _bel_vis_hi10 = min(0.80, _y_top10)
            if _bel_vis_hi10 > y_lo_f10 + 0.005:
                _bel_mid10 = (y_lo_f10 + _bel_vis_hi10) / 2
                ax.text(0.01, _bel_mid10, 'Below threshold (<0.80)', color='#A03000',
                        fontsize=8, ha='left', va='center', fontweight='bold',
                        style='italic', transform=_yt10, zorder=6,
                        bbox=dict(boxstyle='round,pad=0.15', fc='#FDECEA', ec='#D55E00',
                                  alpha=0.85, linewidth=0.6))

            ax.set_xticks(range(len(existing_tiers)))
            ax.set_xticklabels(existing_tiers, rotation=40, ha='right', fontsize=9)
            for _tick10x, _tier10x in zip(ax.get_xticklabels(), existing_tiers):
                _tick10x.set_color(TIER_PALETTE.get(_tier10x, 'black'))
                _tick10x.set_fontweight('bold')
            ax.tick_params(axis='y', labelsize=9, labelcolor='#2C6FAC')
            ax.yaxis.label.set_color('#2C6FAC')
            ax.yaxis.grid(True, color='#2C6FAC', linewidth=0.5, linestyle='-', alpha=0.20, zorder=0)
            ax.set_axisbelow(True)
            ax.set_xlabel("Degrader Tier", fontsize=11)
            ax.set_ylabel("AI Quality Score  (Boltz Confidence / pTM / ipTM;  0–1 scale)", fontsize=11)

            # Legend: single horizontal row, top-right, above data
            from matplotlib.patches import Patch as _Patch10
            conf_patch = _Patch10(facecolor='#AAAAAA', edgecolor='black', linewidth=0.8,
                                  label='Boltz Confidence  (box)')
            handles_f10, labels_f10 = ax.get_legend_handles_labels()
            _all_h10 = [conf_patch] + handles_f10
            _all_l10 = ['Boltz Confidence  (box)'] + labels_f10
            _leg10 = ax.legend(_all_h10, _all_l10,
                               loc='upper right', fontsize=8,
                               framealpha=0.92, fancybox=True,
                               ncol=len(_all_h10))   # all in one row
            _leg10.set_zorder(20)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_04a_AI_Quality_Assessment.png", dpi=300, bbox_inches='tight')
            plt.close()
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_04a_AI_Quality_Assessment.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 04 skipped: {e}")

    # Figure 04b: PA companion (Confidence × ipTM landscape)
    if _pa_has_imgs:
        _fig_04b_pa_ai_quality(df, _pa, _imgs, out_dir, reporter)

    # --- Figure 05: pTM vs ipTM — 2-D Scatter (REDESIGN) ---
    # Each complex is a point: X = pTM (fold confidence), Y = ipTM (interface confidence).
    # Colour = degrader tier. The identity diagonal (pTM = ipTM) divides the space:
    #   ABOVE diagonal → interface more confident than fold = preferred state
    #   BELOW diagonal → fold more confident than interface = potential decoy risk
    # Per-tier median diamonds summarise cluster positions.
    # This is fundamentally more informative than paired box plots.
    if 'ptm' in df.columns and 'iptm' in df.columns and 'degrader_tier' in df.columns:
        try:
            f18_df = df.dropna(subset=['ptm', 'iptm']).copy()
            f18_df['ptm']  = pd.to_numeric(f18_df['ptm'],  errors='coerce')
            f18_df['iptm'] = pd.to_numeric(f18_df['iptm'], errors='coerce')
            f18_df = f18_df.dropna(subset=['ptm', 'iptm'])

            all_vals18 = np.concatenate([f18_df['ptm'].values, f18_df['iptm'].values])
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
            # Zone 1 (above diagonal): ipTM > pTM — preferred state, interface well modelled
            # Zone 2 (diagonal band ±0.03): pTM ≈ ipTM — balanced confidence
            # Zone 3 (below diagonal): pTM > ipTM — fold stronger than interface (decoy risk)
            _diag_pts = np.linspace(_lo18, _hi18, 300)
            from matplotlib.patches import Polygon as _Poly18
            _band = 0.03
            # Zone 1: above diagonal (ipTM > pTM)
            _verts_up = ([(x, x + _band) for x in _diag_pts] +
                         [(_hi18, _hi18), (_lo18, _hi18)])
            ax.add_patch(_Poly18(_verts_up, closed=True,
                                 facecolor='#D4EFDF', alpha=0.45, zorder=0, linewidth=0))
            # Zone 3: below diagonal (pTM > ipTM)
            _verts_dn = ([(x, x - _band) for x in _diag_pts] +
                         [(_hi18, _lo18), (_lo18, _lo18)])
            ax.add_patch(_Poly18(_verts_dn, closed=True,
                                 facecolor='#FADBD8', alpha=0.45, zorder=0, linewidth=0))
            # Zone 2 (diagonal band) is implicitly the white strip between the two coloured zones

            # Identity diagonal line
            ax.plot(_diag_pts, _diag_pts, color='#555', linestyle='--', linewidth=1.3,
                    alpha=0.60, zorder=2)
            # Band boundary dashes
            ax.plot(_diag_pts, _diag_pts + _band, color='#555', linestyle=':',
                    linewidth=0.7, alpha=0.40, zorder=2)
            ax.plot(_diag_pts, _diag_pts - _band, color='#555', linestyle=':',
                    linewidth=0.7, alpha=0.40, zorder=2)

            # Zone text — placed BELOW the upper-left legend to avoid overlap.
            # Legend with ncol=2 is ~2 rows tall; empirically offset to 0.60 axes-y.
            ax.text(0.04, 0.60, 'ipTM > pTM\nInterface well-modelled',
                    transform=ax.transAxes, ha='left', va='top',
                    fontsize=8.5, color='#1A7A4A', fontweight='bold', style='italic',
                    bbox=dict(boxstyle='round,pad=0.22', fc='#D4EFDF', ec='#1A7A4A',
                              alpha=0.80, linewidth=0.7))
            ax.text(0.97, 0.97, 'pTM = ipTM\n(±0.03 band)',
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=7.5, color='#555', style='italic')
            ax.text(0.96, 0.10, 'pTM > ipTM\nFold > interface\n(potential decoy)',
                    transform=ax.transAxes, ha='right', va='bottom',
                    fontsize=8.5, color='#A93226', fontweight='bold', style='italic',
                    bbox=dict(boxstyle='round,pad=0.22', fc='#FADBD8', ec='#A93226',
                              alpha=0.80, linewidth=0.7))

            # Rendering strategy:
            #   Good tier (n~34k) → hexbin background density (YlOrBr)
            #   All other tiers   → scatter with tier colour, drawn on top in order
            # This keeps tier colours distinguishable while showing the density of the bulk.
            for tier in sorted(existing_tiers, key=lambda t: t == 'Perfect_A'):
                sub18 = f18_df[f18_df['degrader_tier'] == tier]
                if sub18.empty:
                    continue
                ptm_vals  = sub18['ptm'].values
                iptm_vals = sub18['iptm'].values
                col18     = TIER_PALETTE.get(tier, '#999')
                if tier == 'Good':
                    ax.hexbin(ptm_vals, iptm_vals, gridsize=45, cmap='YlOrBr',
                              mincnt=1, alpha=0.65, zorder=1)
                else:
                    _is_pa   = tier == 'Perfect_A'
                    _is_top  = tier in ('Perfect_A', 'Perfect_B', 'Best_A')
                    ax.scatter(ptm_vals, iptm_vals,
                               c=col18,
                               alpha=0.92 if _is_pa else (0.55 if _is_top else 0.20),
                               s=110 if _is_pa else (40 if _is_top else 10),
                               edgecolors='black' if _is_pa else 'none',
                               linewidths=1.2 if _is_pa else 0,
                               marker='*' if _is_pa else 'o',
                               zorder=10 if _is_pa else (5 if _is_top else 3),
                               label=None)

            # Per-tier median diamonds — NO inline text annotation;
            # coordinates are embedded in the legend label instead
            from matplotlib.lines import Line2D as _L18
            _leg18 = []
            # Header entries: show both symbol types so the reader knows what each glyph means
            _leg18.append(_L18([0], [0], marker='o', color='w',
                               markerfacecolor='#888888', markeredgecolor='none',
                               markersize=7, alpha=0.55,
                               label='Individual complex  (hexbin/scatter)'))
            _leg18.append(_L18([0], [0], marker='D', color='w',
                               markerfacecolor='#888888', markeredgecolor='black',
                               markersize=9, markeredgewidth=1.0,
                               label='Tier median  (◆ diamond)'))
            for tier in existing_tiers:
                sub18 = f18_df[f18_df['degrader_tier'] == tier]
                if len(sub18) < 2:
                    continue
                mx18, my18 = float(sub18['ptm'].median()), float(sub18['iptm'].median())
                col18 = TIER_PALETTE.get(tier, '#999')
                ax.scatter([mx18], [my18], marker='D', s=160,
                           color=col18, edgecolors='black', linewidths=1.3, zorder=8)
                # Perfect_A callout arrow on main scatter
                if tier == 'Perfect_A':
                    ax.annotate(
                        f'Perfect_A\n(n={len(sub18)}, ★)',
                        xy=(mx18, my18),
                        xytext=(mx18 - 0.035, my18 + 0.015),
                        fontsize=8, color=col18, fontweight='bold', zorder=12,
                        arrowprops=dict(arrowstyle='->', color=col18, lw=1.2,
                                        shrinkA=6, shrinkB=4),
                        bbox=dict(boxstyle='round,pad=0.25', fc='white',
                                  ec=col18, alpha=0.92, linewidth=1.2))
                # Median coords live in the legend label — no floating text cluttering the plot
                _leg18.append(_L18([0], [0], marker='D', color='w',
                                   markerfacecolor=col18, markeredgecolor='black',
                                   markersize=9, markeredgewidth=1.0,
                                   label=f'{tier}  ◆({mx18:.3f}, {my18:.3f})  n={len(sub18):,}'))

            # Marginal distributions — KDE lines (smoother than histograms)
            from scipy.stats import gaussian_kde as _gkde18
            existing_tiers18 = [t for t in existing_tiers if t in f18_df['degrader_tier'].values]
            _kde_x18 = np.linspace(_lo18, _hi18, 300)
            for tier in existing_tiers18:
                sub18_m = f18_df[f18_df['degrader_tier'] == tier]
                _min_kde = 1 if tier == 'Perfect_A' else 10
                if len(sub18_m) < _min_kde:
                    continue
                _col18m = TIER_PALETTE.get(tier, '#999')
                _is_pa18 = tier == 'Perfect_A'
                _lw18    = 2.8 if _is_pa18 else (2.0 if tier == 'Perfect_B' else 1.3)
                _alp18   = 0.95 if _is_pa18 else 0.80
                _zo18    = 10  if _is_pa18 else 5
                try:
                    _kde_ptm  = _gkde18(sub18_m['ptm'].values)(_kde_x18)
                    _kde_iptm = _gkde18(sub18_m['iptm'].values)(_kde_x18)
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
            # Annotate Perfect_A KDE peak on both marginal axes
            _pa_col_f05 = TIER_PALETTE.get('Perfect_A', '#009E73')
            _pa_sub_f05 = f18_df[f18_df['degrader_tier'] == 'Perfect_A']
            if len(_pa_sub_f05) >= 2:
                try:
                    _kde_pa_ptm  = _gkde18(_pa_sub_f05['ptm'].values)(_kde_x18)
                    _kde_pa_iptm = _gkde18(_pa_sub_f05['iptm'].values)(_kde_x18)
                    _peak_ptm_x  = float(_kde_x18[np.argmax(_kde_pa_ptm)])
                    _peak_ptm_y  = float(np.max(_kde_pa_ptm))
                    _peak_iptm_x = float(_kde_x18[np.argmax(_kde_pa_iptm)])
                    _peak_iptm_y = float(np.max(_kde_pa_iptm))
                    # Top KDE: annotation placed BELOW the peak (away from legend box)
                    ax_top.annotate(
                        'Perfect_A ▲',
                        xy=(_peak_ptm_x, _peak_ptm_y),
                        xytext=(_peak_ptm_x + 0.008, _peak_ptm_y * 0.40),
                        fontsize=7, color=_pa_col_f05, fontweight='bold',
                        ha='left', zorder=15,
                        arrowprops=dict(arrowstyle='->', color=_pa_col_f05, lw=1.0,
                                        shrinkA=2, shrinkB=2),
                        bbox=dict(boxstyle='round,pad=0.14', fc='white',
                                  ec=_pa_col_f05, alpha=0.95, linewidth=0.8))
                    # Right KDE: annotation to the left of peak where density is low
                    ax_right.annotate(
                        'Perfect_A ▲',
                        xy=(_peak_iptm_y, _peak_iptm_x),
                        xytext=(_peak_iptm_y * 0.45, _peak_iptm_x + 0.008),
                        fontsize=7, color=_pa_col_f05, fontweight='bold',
                        ha='right', zorder=15,
                        arrowprops=dict(arrowstyle='->', color=_pa_col_f05, lw=1.0,
                                        shrinkA=2, shrinkB=2),
                        bbox=dict(boxstyle='round,pad=0.14', fc='white',
                                  ec=_pa_col_f05, alpha=0.95, linewidth=0.8))
                except Exception:
                    pass
            plt.setp(ax_top.get_xticklabels(), visible=False)
            plt.setp(ax_right.get_yticklabels(), visible=False)
            # Label each marginal clearly so readers know which metric is shown
            ax_top.set_ylabel("Density", fontsize=7, labelpad=2)
            ax_top.set_title("pTM distribution per tier", fontsize=7,
                             color='#444', style='italic', pad=2)
            ax_right.set_xlabel("Density", fontsize=7, labelpad=2)
            ax_right.set_title("ipTM\ndistribution", fontsize=7,
                               color='#444', style='italic', pad=2, loc='left')
            ax_top.set_axisbelow(True)
            ax_right.set_axisbelow(True)
            ax_top.tick_params(axis='y', labelsize=6)
            ax_right.tick_params(axis='x', labelsize=6)
            for _sp in list(ax_top.spines.values()) + list(ax_right.spines.values()):
                _sp.set_visible(False)

            ax.set_xlim(_lo18, _hi18)
            ax.set_ylim(_lo18, _hi18)
            ax.set_xlabel("pTM — Global Fold Confidence  (0–1)", fontsize=11)
            ax.set_ylabel("ipTM — Interface Confidence  (0–1)", fontsize=11)
            ax.tick_params(axis='both', labelsize=9)
            ax.grid(True, color='#EBEBEB', linewidth=0.5, alpha=0.6, zorder=1)
            # Legend: placed in the top KDE marginal (left side) to keep scatter uncluttered
            ax_top.legend(handles=_leg18, loc='upper left', ncol=3, fontsize=7.0,
                          framealpha=0.92, fancybox=True,
                          title='Tier  (◆ = median pTM, ipTM)', title_fontsize=7.0,
                          bbox_to_anchor=(0.0, 1.0))
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_05_pTM_vs_ipTM_by_Tier.png", dpi=300, bbox_inches='tight')
            plt.close()
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_05_pTM_vs_ipTM_by_Tier.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 05 skipped: {e}")

    # --- Figure 06: Active Site RMSD vs. Reference Control by Tier ---
    # Column alias detection: try multiple plausible names before skipping
    _rmsd17_aliases = [
        'Active_Site_RMSD', 'active_site_rmsd', 'ActiveSite_RMSD', 'RMSD_ActiveSite',
        'rmsd_active_site', 'binding_site_rmsd', 'Binding_Site_RMSD',
        'RMSD', 'rmsd', 'structural_rmsd', 'Structural_RMSD',
        'site_rmsd', 'Site_RMSD', 'AS_RMSD', 'as_rmsd',
        'reference_rmsd', 'Reference_RMSD', 'RMSD_to_reference',
    ]
    # Broad partial-match fallback: any column with 'rmsd' in the name
    _rmsd17_col = next((c for c in _rmsd17_aliases if c in df.columns), None)
    if _rmsd17_col is None:
        _rmsd17_col = next((c for c in df.columns if 'rmsd' in c.lower()), None)
    if _rmsd17_col is None:
        _rmsd_candidates = [c for c in df.columns if 'rmsd' in c.lower() or 'active_site' in c.lower()]
        reporter.log(f"  ! Figure 06 skipped: no RMSD column found. "
                     f"All column names tried: {_rmsd17_aliases}. "
                     f"Partial matches in CSV (containing 'rmsd'/'active_site'): {_rmsd_candidates}. "
                     f"Full column list: {list(df.columns)}")
    if _rmsd17_col and 'degrader_tier' in df.columns:
        try:
            if _rmsd17_col != 'Active_Site_RMSD':
                df = df.rename(columns={_rmsd17_col: 'Active_Site_RMSD'})
            f17_df = df.dropna(subset=['Active_Site_RMSD']).copy()
            f17_df['Active_Site_RMSD'] = pd.to_numeric(f17_df['Active_Site_RMSD'], errors='coerce')
            f17_df = f17_df.dropna(subset=['Active_Site_RMSD'])
            # Show ALL data — no percentile clip so small-n tiers (e.g. Perfect_A n=5)
            # retain every value.  y_ceil caps the display at 4 Å, but is extended if any
            # small-n tier has a value beyond that so its dots always land inside the axes.
            f17_plot = f17_df.copy()
            n_outliers = 0
            valid_tiers_f17 = [t for t in existing_tiers if t in f17_plot['degrader_tier'].values]
            _sn_max17 = 0.0
            for _t17s in valid_tiers_f17:
                _sv17 = f17_plot.loc[f17_plot['degrader_tier'] == _t17s, 'Active_Site_RMSD']
                if 0 < len(_sv17) <= 30:
                    _sn_max17 = max(_sn_max17, float(_sv17.max()))
            _raw_max = float(f17_df['Active_Site_RMSD'].max())
            _violin_cap = 3.0   # violin body capped here
            _raw_pct98  = float(f17_df['Active_Site_RMSD'].quantile(0.98)) if len(f17_df) > 0 else 3.0
            y_ceil      = max(_raw_pct98 * 1.08, 3.6)  # show 98th percentile + headroom
            # Clip violin plot data at _violin_cap — outliers are shown as discrete scatter above
            f17_plot['Active_Site_RMSD'] = np.clip(f17_plot['Active_Site_RMSD'].values, 0, _violin_cap)
            # Compute % of data in each zone (across all tiers together, using unclipped values)
            total_n = len(f17_df)
            pct_green  = (f17_df['Active_Site_RMSD'] <= 1.0).sum() / total_n * 100
            pct_yellow = ((f17_df['Active_Site_RMSD'] > 1.0) & (f17_df['Active_Site_RMSD'] <= 2.0)).sum() / total_n * 100
            pct_red    = (f17_df['Active_Site_RMSD'] > 2.0).sum() / total_n * 100
            fig, ax = plt.subplots(figsize=(12, 7))
            # Coloured quality zones — distinct colours at higher alpha so they read apart
            ax.axhspan(0.0, 1.0, alpha=0.13, color='#009E73', zorder=0)   # green  = Excellent
            ax.axhspan(1.0, 2.0, alpha=0.13, color='#E69F00', zorder=0)   # amber  = Acceptable
            ax.axhspan(2.0, y_ceil + 0.1, alpha=0.13, color='#D55E00', zorder=0)  # red = Diverged
            # Zone boundary lines
            ax.axhline(y=1.0, color='#009E73', linestyle='--', alpha=0.7, linewidth=1.1)
            ax.axhline(y=2.0, color='#D55E00', linestyle='--', alpha=0.7, linewidth=1.1)
            # Violin plot — uses clipped data so bodies remain within the continuous axis
            sns.violinplot(data=f17_plot, x='degrader_tier', y='Active_Site_RMSD',
                           order=valid_tiers_f17, palette=TIER_PALETTE,
                           inner='quartile', cut=0, linewidth=1.1, ax=ax, zorder=3)
            # Very sparse strip for large-n tiers only
            _large_n_tiers17 = [t for t in valid_tiers_f17
                                 if len(f17_plot[f17_plot['degrader_tier'] == t]) > 30]
            if _large_n_tiers17:
                _f17_large = f17_plot[f17_plot['degrader_tier'].isin(_large_n_tiers17)]
                sns.stripplot(data=_f17_large, x='degrader_tier', y='Active_Site_RMSD',
                              order=valid_tiers_f17, color='black', alpha=0.04,
                              size=1.5, jitter=True, ax=ax, zorder=2)
            # For tiers with very few data points (≤30), draw individual dots prominently
            # so every value is visible — avoids them being lost inside the violin body
            for _i17, _tier17 in enumerate(valid_tiers_f17):
                _sub17 = f17_plot.loc[f17_plot['degrader_tier'] == _tier17, 'Active_Site_RMSD'].dropna()
                if len(_sub17) <= 30:
                    _col17 = TIER_PALETTE.get(_tier17, '#999')
                    _jit17 = np.random.default_rng(42).uniform(-0.12, 0.12, size=len(_sub17))
                    ax.scatter(_i17 + _jit17, _sub17.values,
                               color=_col17, edgecolors='black', s=55,
                               linewidths=0.8, alpha=0.88, zorder=7)
            # Outlier scatter: actual RMSD values above 3 Å (unclipped)
            for _ti, tier in enumerate(valid_tiers_f17):
                sub_f6 = f17_df[f17_df['degrader_tier'] == tier]
                _out = sub_f6[sub_f6['Active_Site_RMSD'] > _violin_cap]['Active_Site_RMSD'].values
                if len(_out) == 0:
                    continue
                _jit_f6 = np.random.default_rng(17 + _ti).uniform(-0.22, 0.22, size=len(_out))
                ax.scatter(_ti + _jit_f6, _out,
                           color=TIER_PALETTE.get(tier, '#999'), s=38, alpha=0.65,
                           marker='o', edgecolors='white', linewidths=0.4, zorder=4)
                if len(_out) >= 1:
                    _lbl_col = TIER_PALETTE.get(tier, '#999')
                    ax.text(_ti, y_ceil + 0.04, f'+{len(_out)}',
                            ha='center', va='bottom', fontsize=9.5,
                            color=_lbl_col, fontweight='bold',
                            zorder=30, clip_on=False,
                            bbox=dict(boxstyle='round,pad=0.18', fc='white',
                                      ec=_lbl_col, alpha=1.0, linewidth=1.0))

            # Capped axis, fixed ticks, and dashed break marker at 3 Å
            ax.axhline(3.0, color='#888', linestyle=':', linewidth=0.9, alpha=0.7)
            ax.set_ylim(0, y_ceil)
            _y_ticks_f6 = [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
            _y_extra_f6 = [t for t in np.arange(3.5, y_ceil, 0.5) if t <= y_ceil - 0.1]
            _y_ticks_f6 += _y_extra_f6
            ax.set_yticks(_y_ticks_f6)
            ax.set_yticklabels([f'{t:.1f}' for t in _y_ticks_f6])
            ax.tick_params(axis='y', labelsize=7)
            ax.set_xlabel("Degrader Tier", fontsize=11)
            ax.set_ylabel("Active Site RMSD (Å)  — lower = better", fontsize=10)
            ax.text(0.99, 0.99, f' ',
                    ha='right', va='top', fontsize=7, color='#888', style='italic',
                    transform=ax.transAxes)
            # Build two-line x-tick labels: tier name + n / median stats (using unclipped values)
            _xtick_labels17 = []
            for _tier17 in valid_tiers_f17:
                _sub17 = f17_df.loc[f17_df['degrader_tier'] == _tier17,
                                    'Active_Site_RMSD'].dropna()
                if len(_sub17) > 0:
                    _med17 = float(_sub17.median())
                    _xtick_labels17.append(
                        f'{_tier17}\nn={len(_sub17):,}  med={_med17:.2f}Å'
                    )
                else:
                    _xtick_labels17.append(_tier17)
            ax.set_xticks(range(len(valid_tiers_f17)))
            ax.set_xticklabels(_xtick_labels17, rotation=40, ha='right', fontsize=7.5)
            for _tick17, _tier17 in zip(ax.get_xticklabels(), valid_tiers_f17):
                _tick17.set_color(TIER_PALETTE.get(_tier17, 'black'))
                _tick17.set_fontweight('bold')

            # Per-zone count labels (90°-rotated) inside each RMSD zone band
            _zone_count_defs = [
                (0.0, 1.0,    '#009E73', '<1 Å'),
                (1.0, 2.0,    '#E69F00', '1–2 Å'),
                (2.0, 3.0,    '#D55E00', '2–3 Å'),
                (3.0, np.inf, '#990000', '>3 Å'),
            ]
            for _zlo, _zhi, _zcol, _zlbl in _zone_count_defs:
                _zmask = (f17_df['Active_Site_RMSD'] >= _zlo) & (f17_df['Active_Site_RMSD'] < _zhi)
                _zn = int(_zmask.sum())
                _zy = min(_zhi, y_ceil) * 0.98
                ax.text(len(valid_tiers_f17) - 0.3, _zy,
                        f'n={_zn:,}  {_zlbl}',
                        ha='center', va='top', fontsize=6.5, color=_zcol, fontweight='bold',
                        rotation=90, zorder=25, clip_on=False,
                        bbox=dict(boxstyle='round,pad=0.18', fc='white', ec=_zcol,
                                  alpha=1.0, linewidth=0.5))

            # Zone band legend (top-left) replaces inline text labels
            from matplotlib.patches import Patch as _P6z
            _zone_legend = [
                _P6z(facecolor='#D4EDDA', alpha=0.55, edgecolor='none', label=f'Excellent  (< 1.0 Å)  {pct_green:.0f}%'),
                _P6z(facecolor='#FFF3CD', alpha=0.55, edgecolor='none', label=f'Acceptable  (1.0–2.0 Å)  {pct_yellow:.0f}%'),
                _P6z(facecolor='#F8D7DA', alpha=0.55, edgecolor='none', label=f'Diverged  (> 2.0 Å)  {pct_red:.0f}%'),
            ]
            _tier_handles, _tier_labels = ax.get_legend_handles_labels()
            _combined_handles = _zone_legend + _tier_handles
            _combined_labels  = [h.get_label() for h in _zone_legend] + _tier_labels
            ax.legend(_combined_handles, _combined_labels,
                      loc='upper left', fontsize=8, framealpha=0.92)
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_06_ActiveSite_RMSD_by_Tier.png", dpi=300, bbox_inches='tight')
            plt.close()
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_06_ActiveSite_RMSD_by_Tier.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 06 skipped: {e}")

    # --- Figure 07: Feature Correlation Matrix (Spearman ρ) ---
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
    for _alias in ['tier_val', 'tier_value', 'Tier_Score', 'tier_score',
                   'Tier_Numeric', 'degrader_tier_numeric', 'Degrader_Tier_Numeric']:
        if _alias in df.columns and 'tier_numeric' not in df.columns:
            df['tier_numeric'] = pd.to_numeric(df[_alias], errors='coerce')
    # Last-resort: derive numeric rank from the categorical degrader_tier string
    if 'tier_numeric' not in df.columns and 'degrader_tier' in df.columns:
        _tier_rank_map = {t: i for i, t in enumerate(TIER_ORDER_LOGIC)}
        df['tier_numeric'] = df['degrader_tier'].map(_tier_rank_map)
    corr_cols = [c for c in _corr_col_labels if c in df.columns]
    if len(corr_cols) >= 2:
        from scipy import stats as _sp_stats
        from scipy.cluster import hierarchy as _sch
        import scipy.spatial.distance as _ssd

        # Fill each column with its median instead of dropping rows —
        # prevents tier_numeric losing all variance when NaN rows are removed.
        sub_f3 = df[corr_cols].copy()
        for _c in corr_cols:
            sub_f3[_c] = pd.to_numeric(sub_f3[_c], errors='coerce')
            _med_c = sub_f3[_c].median()
            sub_f3[_c] = sub_f3[_c].fillna(_med_c if not np.isnan(_med_c) else 0)
        sub_f3 = sub_f3.dropna()          # only drop rows still fully NaN
        corr_raw = sub_f3.corr(method='spearman')

        # Build p-value matrix for significance stars
        _n_f3 = len(sub_f3)
        p_mat = np.ones((len(corr_cols), len(corr_cols)))
        for _ii, _c1 in enumerate(corr_cols):
            for _jj, _c2 in enumerate(corr_cols):
                if _ii != _jj:
                    _, _pv = _sp_stats.spearmanr(sub_f3[_c1], sub_f3[_c2])
                    p_mat[_ii, _jj] = _pv
        p_df_f3 = pd.DataFrame(p_mat, index=corr_raw.index, columns=corr_raw.columns)

        # Benjamini-Hochberg FDR correction on all pairwise p-values
        _n_feat_f7 = len(corr_cols)
        _pv_flat_f7, _pv_idx_f7 = [], []
        for _ii7 in range(_n_feat_f7):
            for _jj7 in range(_n_feat_f7):
                if _ii7 != _jj7:
                    _pv_flat_f7.append(p_mat[_ii7, _jj7])
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
        p_adj_df_f7 = pd.DataFrame(p_adj_mat_f7, index=corr_raw.index, columns=corr_raw.columns)

        # Hierarchical clustering reorder (features grouped by similarity)
        _dist = _ssd.squareform(np.clip(1 - np.abs(corr_raw.values), 0, 2))
        _link = _sch.linkage(_dist, method='average')
        _ord  = _sch.leaves_list(_link)
        corr_raw  = corr_raw.iloc[_ord, _ord]
        p_df_f3   = p_df_f3.iloc[_ord, _ord]
        p_adj_df_f7 = p_adj_df_f7.iloc[_ord, _ord]

        # Rename for display
        _display_names = [_corr_col_labels[corr_cols[i]] for i in _ord]
        corr_raw.index   = _display_names
        corr_raw.columns = _display_names
        p_df_f3.index    = _display_names
        p_df_f3.columns  = _display_names
        p_adj_df_f7.index   = _display_names
        p_adj_df_f7.columns = _display_names

        # Build annotation: ρ value + significance stars
        def _sig_stars(p):
            return '***' if p < 0.001 else ('**' if p < 0.01 else ('*' if p < 0.05 else ''))

        annot_f3 = pd.DataFrame('', index=corr_raw.index, columns=corr_raw.columns)
        for _ri in corr_raw.index:
            for _ci in corr_raw.columns:
                _rho  = corr_raw.loc[_ri, _ci]
                _pv   = p_adj_df_f7.loc[_ri, _ci]
                if pd.isna(_rho):
                    annot_f3.loc[_ri, _ci] = '—'
                    continue
                _star = _sig_stars(_pv)
                annot_f3.loc[_ri, _ci] = f'{_rho:.2f}\n{_star}' if _star else f'{_rho:.2f}'

        # k=1 excludes diagonal from mask so self-correlation (1.00) is shown
        mask_f3 = np.triu(np.ones_like(corr_raw, dtype=bool), k=1)
        for _dn in corr_raw.index:
            annot_f3.loc[_dn, _dn] = '1.00'
        # Dendrogram removed — features remain reordered by hierarchical clustering
        # (see _ord above); the footnote states this. Single-axes heatmap only.
        fig_f3, ax_f3 = plt.subplots(figsize=(12, 11))
        sns.heatmap(corr_raw, mask=mask_f3, annot=annot_f3, fmt='',
                    cmap='coolwarm', vmin=-1, vmax=1, center=0,
                    square=True, linewidths=0.5,
                    annot_kws={'size': 8.5},
                    cbar_kws={'label': 'Spearman ρ  (−1 = perfect negative, +1 = perfect positive)',
                              'shrink': 0.75},
                    ax=ax_f3)
        # Footnote inside the plot — top-right corner, well clear of heatmap cells
        ax_f3.text(0.99, 0.99,
                   f'★ p<0.05  ★★ p<0.01  ★★★ p<0.001  (BH FDR-corrected)  |  Features reordered by hierarchical clustering  |  n = {_n_f3:,} complexes',
                   transform=ax_f3.transAxes, ha='right', va='top',
                   fontsize=7.5, color='#555', style='italic',
                   bbox=dict(boxstyle='round,pad=0.35', fc='white', ec='#ccc',
                             alpha=0.88, linewidth=0.7))
        # Labels are now short (≤7 chars) — no truncation needed
        ax_f3.set_xticklabels([l.get_text() for l in ax_f3.get_xticklabels()],
                              fontsize=9, rotation=45, ha='right')
        ax_f3.set_yticklabels([l.get_text() for l in ax_f3.get_yticklabels()],
                              fontsize=9, rotation=90, va='center', ha='right')

        # Colour significance stars: white on dark cells (|ρ|>0.5), dark on light cells
        for _txt_art in ax_f3.texts:
            _raw = _txt_art.get_text()
            _lines = _raw.split('\n')
            if len(_lines) == 2 and _lines[1] and all(c == '*' for c in _lines[1]):
                try:
                    _rho_val = float(_lines[0])
                    _txt_art.set_color('white' if abs(_rho_val) > 0.5 else '#1a1a1a')
                except ValueError:
                    pass
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_07_Feature_Correlations.png", dpi=300, bbox_inches='tight')
        plt.close()
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_07_Feature_Correlations.png').resolve()}")

    # --- Figure 08: Tier Quality Summary — Cleveland Dot Plot (multi-metric) ---
    # Each metric occupies a horizontal row; coloured dots show each tier's normalised
    # median score.  A thin grey range line connects the min-to-max dot for that metric,
    # immediately revealing which metrics discriminate tiers the most.
    # Raw median values annotate each dot for quantitative readability.
    _f08_metric_map = {
        'Boltz_Model_Confidence': 'AI Confidence', 'Mechanistic_Fingerprint_Score': 'Mech. Fingerprint',
        'Binding_Probability': 'Binding Prob.', 'identity_pct': 'Seq. Identity (%)',
        'SN2_Attack_Angle': 'SN2 Angle (°)', 'Active_Site_RMSD': 'RMSD (Å, inv.)',
        'SN2_Trajectory_Deviation_A': 'SN2 Traj. Dev. (Å)',
        'soft_catalytic_score': 'Soft Catalytic Score',
    }
    _f08_invert = {'Active_Site_RMSD', 'SN2_Trajectory_Deviation_A'}
    _f08_cols = [c for c in _f08_metric_map if c in df.columns and 'degrader_tier' in df.columns]
    if len(_f08_cols) >= 2:
        try:
            _diag08_parts = []
            for _t08d in existing_tiers:
                _s08d = df[df['degrader_tier'] == _t08d]
                _v08d = pd.to_numeric(
                    _s08d['Mechanistic_Fingerprint_Score'] if 'Mechanistic_Fingerprint_Score' in _s08d.columns
                    else pd.Series(dtype=float), errors='coerce').dropna()
                if len(_v08d):
                    _diag08_parts.append(f"{_t08d}:{float(_v08d.mean()):.2f}(n={len(_v08d)})")
                else:
                    _diag08_parts.append(f"{_t08d}:NA")
            pass  # Fig08 diag suppressed
            _f08_rows = []
            for tier in existing_tiers:
                sub = df[df['degrader_tier'] == tier]
                if sub.empty:
                    continue
                row = {'Tier': tier}
                for col in _f08_cols:
                    vals = pd.to_numeric(sub[col], errors='coerce').dropna()
                    row[_f08_metric_map[col]] = float(vals.mean()) if len(vals) else np.nan
                _f08_rows.append(row)
            if _f08_rows:
                _f08_df   = pd.DataFrame(_f08_rows).set_index('Tier')
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
                    meds = [float(np.median(rng_b.choice(vals, size=len(vals))))
                            for _ in range(n_boot)]
                    return float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))

                # Precompute per-tier bootstrap CIs for each metric (raw scale)
                _f08_ci = {}  # (tier, metric) -> (lo_norm, hi_norm)
                for col in _f08_cols:
                    disp = _f08_metric_map[col]
                    mn_raw = float(_f08_df[disp].min())
                    mx_raw = float(_f08_df[disp].max())
                    _inv = col in _f08_invert
                    for tier in tiers_f08:
                        sub_ci = df[df['degrader_tier'] == tier]
                        raw_vals = pd.to_numeric(sub_ci[col], errors='coerce').dropna().values
                        if len(raw_vals) >= 2:
                            lo_r, hi_r = _bsci(raw_vals)
                            if mx_raw > mn_raw:
                                lo_n = (lo_r - mn_raw) / (mx_raw - mn_raw)
                                hi_n = (hi_r - mn_raw) / (mx_raw - mn_raw)
                            else:
                                lo_n, hi_n = 0.5, 0.5
                            if _inv:
                                lo_n, hi_n = 1 - hi_n, 1 - lo_n
                            _f08_ci[(tier, disp)] = (float(np.clip(lo_n, 0, 1)),
                                                     float(np.clip(hi_n, 0, 1)))
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
                    ax08.xaxis.grid(True, color='#EBEBEB', linewidth=0.7, linestyle='-', zorder=0)
                    ax08.yaxis.grid(False)

                    for _rmi, metric in enumerate(_row_mets):
                        _base_y = _rmi * 1.0   # row spacing within subplot
                        norm_vals  = _f08_norm[metric].values
                        valid_mask = ~np.isnan(norm_vals)
                        _is_tied08 = _met_tied_map08.get(metric, False)

                        _dot_data08 = []
                        for ti, tier in enumerate(tiers_f08):
                            nv = float(_f08_norm.loc[tier, metric])
                            rv = float(_f08_df.loc[tier, metric])
                            if np.isnan(nv):
                                continue
                            _dot_data08.append((nv, rv, TIER_PALETTE.get(tier, '#999'), tier))
                        _dot_data08.sort(key=lambda d: d[0])
                        _nd08 = len(_dot_data08)
                        _all_nvs08 = [d[0] for d in _dot_data08]
                        _is_tied08 = bool(_nd08 > 1 and (max(_all_nvs08) - min(_all_nvs08)) < 0.02)

                        if _is_tied08:
                            _raw_tied08 = float(_dot_data08[0][1]) if _dot_data08 else float('nan')
                            if not any(m == metric for m, _ in _tied_metrics08):
                                _tied_metrics08.append((metric, _raw_tied08))

                        # Range line
                        if valid_mask.sum() >= 2:
                            if _is_tied08:
                                ax08.hlines(_base_y, 0.0, 1.0, colors='#D4A000',
                                            linewidth=1.0, linestyle='--', alpha=0.55, zorder=1)
                            else:
                                ax08.hlines(_base_y, float(np.nanmin(norm_vals)),
                                            float(np.nanmax(norm_vals)),
                                            colors='#CCCCCC', linewidth=2.5, zorder=1)

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
                                         edgecolors='black', linewidths=0.7)

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
                                ax08.text(_nv_c, _y_pos08, f'{rv:.2f}',
                                          ha='center', va='bottom', fontsize=6.5,
                                          color=col_f08, fontweight='bold', zorder=5,
                                          rotation=90,
                                          bbox=dict(boxstyle='round,pad=0.07', fc='white',
                                                    ec='none', alpha=0.75))

                    ax08.set_yticks([i * 1.0 for i in range(len(_row_mets))])
                    ax08.set_yticklabels(_row_mets, fontsize=10.5)
                    ax08.set_ylim(-0.5, len(_row_mets) - 0.5 + 0.55)
                    ax08.set_xlim(-0.08, 1.12)
                    ax08.spines[['top', 'right']].set_visible(False)

                for _ax08x in axes_f08:
                    _ax08x.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
                    _ax08x.set_xticklabels(['0\n(worst)', '0.25', '0.50', '0.75', '1.0\n(best)'],
                                           fontsize=8.5)
                axes_f08[-1].set_xlabel(
                    "Normalised Median Score per Tier  (0 = worst, 1 = best within metric)",
                    fontsize=11)

                from matplotlib.lines import Line2D as _L08
                _leg08 = [_L08([0],[0], marker='o', color='w',
                               markerfacecolor=TIER_PALETTE.get(t,'#999'),
                               markeredgecolor='black', markeredgewidth=0.8,
                               markersize=10, label=t) for t in tiers_f08]
                _leg08.append(_L08([0],[0], color='none', linewidth=0,
                                   label='★ RMSD metrics: inverted — lower RMSD = better score'))
                _leg08.append(_L08([0],[0], color='#999', linewidth=5.5, alpha=0.3,
                                   label='95% bootstrap CI on median (1 000 resamples)'))
                if _any_tied08:
                    _leg08.append(_L08([0],[0], color='#D4A000', linewidth=1.5,
                                       linestyle='--',
                                       label='⚡ Amber dashed = metric tied across all tiers'))
                    for _tm_name, _tm_val in _tied_metrics08:
                        _leg08.append(_L08([0],[0], color='none', linewidth=0,
                                           label=f'   {_tm_name}: all tiers = {_tm_val:.2f}'))
                axes_f08[0].legend(handles=_leg08, loc='upper left',
                                   bbox_to_anchor=(0.0, 1.0),
                                   ncol=max(1, (len(_leg08) + 1) // 2),
                                   fontsize=7.5, framealpha=0.92, fancybox=True)
                plt.tight_layout()
                plt.savefig(out_dir / "Figure_08_Tier_Quality_DotPlot.png", dpi=300, bbox_inches='tight')
                plt.close()
                reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_08_Tier_Quality_DotPlot.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 08 skipped: {e}")

    # --- Figure 09: Mechanistic State Cross-Tab (Halide Stabilisation × Carboxylate Clamp) ---
    # Column alias detection for both mechanistic columns
    for _hs_alias in [
        'Halide_Stabilisation', 'halide_stabilisation', 'halide_stabilization',
        'HalideStabilisation', 'halide_shield', 'Halide_Shield',
        'Aromatic_Shield', 'aromatic_shield', 'halide_interaction', 'Halide_Interaction',
        'has_halide_shield', 'Has_Halide_Shield', 'halide_clamp', 'Halide_Clamp',
    ]:
        if _hs_alias in df.columns and 'Halide_Stabilisation' not in df.columns:
            df['Halide_Stabilisation'] = df[_hs_alias]
    # Partial-match fallback: any column with 'halide' or 'stabilisa'/'shield'
    if 'Halide_Stabilisation' not in df.columns:
        _hs_fb = next((c for c in df.columns
                       if any(k in c.lower() for k in ('halide', 'stabilisa', 'shield'))), None)
        if _hs_fb:
            df['Halide_Stabilisation'] = df[_hs_fb]

    for _cc_alias in [
        'Carboxylate_Clamp', 'carboxylate_clamp', 'CarboxylateClamp',
        'carboxylate_clamping', 'Carboxylate_Clamping',
        'Asp_Clamp', 'asp_clamp', 'ASP_Clamp', 'asp_clamping',
        'carboxylate_contact', 'Carboxylate_Contact',
        'has_carboxylate_clamp', 'Has_Carboxylate_Clamp',
    ]:
        if _cc_alias in df.columns and 'Carboxylate_Clamp' not in df.columns:
            df['Carboxylate_Clamp'] = df[_cc_alias]
    # Partial-match fallback: any column with 'clamp' or 'carboxyl'
    if 'Carboxylate_Clamp' not in df.columns:
        _cc_fb = next((c for c in df.columns
                       if any(k in c.lower() for k in ('clamp', 'carboxyl', 'asp_clamp'))), None)
        if _cc_fb:
            df['Carboxylate_Clamp'] = df[_cc_fb]
    if 'Halide_Stabilisation' not in df.columns or 'Carboxylate_Clamp' not in df.columns:
        _mech_missing = []
        if 'Halide_Stabilisation' not in df.columns:
            _mech_missing.append('Halide_Stabilisation')
        if 'Carboxylate_Clamp' not in df.columns:
            _mech_missing.append('Carboxylate_Clamp')
        _mech_candidates = [c for c in df.columns
                            if any(k in c.lower() for k in ('halide', 'clamp', 'carboxyl', 'stabil'))]
        reporter.log(f"  ! Figure 09 skipped: columns not found: {_mech_missing}. "
                     f"Partial matches in CSV: {_mech_candidates}. "
                     f"Full column list: {list(df.columns)}")
    if 'Halide_Stabilisation' in df.columns and 'Carboxylate_Clamp' in df.columns and 'degrader_tier' in df.columns:
        try:
            mech_df = df.copy()
            mech_df['HS'] = mech_df['Halide_Stabilisation'].map(
                lambda x: 'Stabilised' if str(x).strip().lower() in ('true', '1', 'yes') else 'Unstabilised')
            mech_df['CC'] = mech_df['Carboxylate_Clamp'].map(
                lambda x: 'Clamped' if str(x).strip().lower() in ('true', '1', 'yes') else 'Unclamped')
            mech_df['Mech_State'] = mech_df['HS'] + '\n' + mech_df['CC']

            state_order = ['Stabilised\nClamped', 'Stabilised\nUnclamped',
                           'Unstabilised\nClamped', 'Unstabilised\nUnclamped']
            cross = pd.crosstab(mech_df['degrader_tier'], mech_df['Mech_State'])
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
                _chi2_title = (f"Mechanistic State Cross-Tab  (Halide Stabilisation × Carboxylate Clamp)"
                               f"  ·  χ²({_dof_09}) = {_chi2_09:.1f},  "
                               f"p {'< 0.001' if _p_09 < 0.001 else f'= {_p_09:.3f}'}")
            except Exception:
                _chi2_title = "Mechanistic State Cross-Tab  (Halide Stabilisation × Carboxylate Clamp)"

            fig19, ax19 = plt.subplots(figsize=(11, 7))
            sns.heatmap(cross_pct, annot=annot_labels, fmt='', cmap='YlOrRd', linewidths=0.5,
                        annot_kws={"size": 10, "weight": "bold"},
                        cbar_kws={'label': '% of tier'}, ax=ax19)
            ax19.text(0.99, 1.01, _chi2_title, transform=ax19.transAxes,
                      ha='right', va='bottom', fontsize=8.5, style='italic',
                      color='#333', clip_on=False,
                      bbox=dict(boxstyle='round,pad=0.25', fc='white',
                                ec='#ccc', alpha=0.95, linewidth=0.6))
            ax19.set_xlabel("Mechanistic State", fontsize=11)
            ax19.set_ylabel("Degrader Tier", fontsize=11)
            ax19.set_xticklabels(ax19.get_xticklabels(), rotation=15, ha='right', fontsize=9)
            # Colour y-tick labels by tier (same palette as all other figures)
            for tick in ax19.get_yticklabels():
                tier_name = tick.get_text()
                tick.set_color(TIER_PALETTE.get(tier_name, 'black'))
                tick.set_fontweight('bold')
                tick.set_rotation(0)
            plt.tight_layout()
            plt.savefig(out_dir / "Figure_09_Mech_State_CrossTab.png", dpi=300, bbox_inches='tight')
            plt.close()
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_09_Mech_State_CrossTab.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 09 skipped: {e}")

    # --- Figure 10: Mechanistic Fingerprint Score — Mean±CI dot plot ---
    if 'Mechanistic_Fingerprint_Score' in df.columns and 'degrader_tier' in df.columns:
        from scipy import stats as _scipy_stats
        f13_data = df.dropna(subset=['Mechanistic_Fingerprint_Score'])
        fig, ax = plt.subplots(figsize=(13, 8))
        ax.axhspan(0.88, 1.01, alpha=0.10, color='#009E73', zorder=0, label='_nolegend_')
        ax.axhspan(0.50, 0.88, alpha=0.09, color='#E69F00', zorder=0, label='_nolegend_')
        ax.axhspan(0.00, 0.50, alpha=0.08, color='#D55E00', zorder=0, label='_nolegend_')
        ax.text(0.99, 0.942, 'Strong zone  (≥0.88)', color='#007A50',
                fontsize=8, ha='right', va='center', fontweight='bold',
                style='italic', transform=ax.transAxes,
                bbox=dict(boxstyle='round,pad=0.15', fc='#D6F5EB', ec='#009E73',
                          alpha=0.85, linewidth=0.6))
        ax.text(0.99, 0.798, 'Moderate zone  (0.50–0.88)', color='#8A6000',
                fontsize=8, ha='right', va='center', fontweight='bold',
                style='italic', transform=ax.transAxes,
                bbox=dict(boxstyle='round,pad=0.15', fc='#FFF3CC', ec='#E69F00',
                          alpha=0.85, linewidth=0.6))
        ax.text(0.99, 0.435, 'Weak zone  (<0.50)', color='#A03000',
                fontsize=8, ha='right', va='center', fontweight='bold',
                style='italic', transform=ax.transAxes,
                bbox=dict(boxstyle='round,pad=0.15', fc='#FDECEA', ec='#D55E00',
                          alpha=0.85, linewidth=0.6))
        ax.axhline(y=0.88, color='#009E73', linestyle='--', linewidth=1.2, alpha=0.7)
        ax.axhline(y=0.50, color='#D55E00', linestyle=':', linewidth=1.0, alpha=0.7)
        sns.stripplot(data=f13_data, x='degrader_tier', y='Mechanistic_Fingerprint_Score',
                      order=existing_tiers, palette=TIER_PALETTE, alpha=0.15, size=3.0,
                      jitter=0.28, ax=ax, zorder=1)
        _f10_stats = {}
        for _t10 in existing_tiers:
            _v10 = f13_data.loc[f13_data['degrader_tier'] == _t10,
                                 'Mechanistic_Fingerprint_Score'].dropna()
            if len(_v10) >= 2:
                _f10_stats[_t10] = {
                    'mean': float(_v10.mean()),
                    'ci95': float(_scipy_stats.sem(_v10) * _scipy_stats.t.ppf(0.975, len(_v10) - 1)),
                    'med':  float(_v10.median()),
                    'n':    len(_v10),
                    'vals': _v10,
                }
        for i, tier in enumerate(existing_tiers):
            if tier not in _f10_stats:
                continue
            _st = _f10_stats[tier]
            mean_v, ci95, med_v, vals = _st['mean'], _st['ci95'], _st['med'], _st['vals']
            ax.plot([i, i], [mean_v - ci95, mean_v + ci95], color='black', linewidth=2.2, zorder=4)
            ax.scatter([i], [mean_v], color=TIER_PALETTE.get(tier, '#999'), s=180,
                       zorder=5, edgecolors='black', linewidths=1.2)
            pct_strong   = (vals >= 0.88).sum() / len(vals) * 100
            pct_moderate = ((vals >= 0.50) & (vals < 0.88)).sum() / len(vals) * 100
            pct_weak     = (vals < 0.50).sum() / len(vals) * 100
            if pct_strong > 0:
                ax.text(i, 0.91, f'{pct_strong:.0f}%', ha='center', va='bottom',
                        fontsize=7, color='#007A50', fontweight='bold', zorder=10,
                        bbox=dict(boxstyle='round,pad=0.10', fc='white', ec='#009E73',
                                  alpha=0.75, linewidth=0.5))
            if pct_moderate > 0:
                ax.text(i, 0.69, f'{pct_moderate:.0f}%', ha='center', va='center',
                        fontsize=7, color='#8A6000', fontweight='bold', zorder=10,
                        bbox=dict(boxstyle='round,pad=0.10', fc='white', ec='#E69F00',
                                  alpha=0.75, linewidth=0.5))
            if pct_weak > 0:
                ax.text(i, 0.25, f'{pct_weak:.0f}%', ha='center', va='center',
                        fontsize=7, color='#A03000', fontweight='bold', zorder=10,
                        bbox=dict(boxstyle='round,pad=0.10', fc='white', ec='#D55E00',
                                  alpha=0.75, linewidth=0.5))
        ax.set_ylim(-0.02, 1.06)
        ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.88, 0.9, 1.0])
        ax.set_xlabel("Degrader Tier", fontsize=11)
        ax.set_ylabel("Mechanistic Fingerprint Score  (fraction of 8 catalytic residues within 6 Å)", fontsize=10)
        ax.set_xticks(range(len(existing_tiers)))
        _xtlbl10 = []
        for _t10 in existing_tiers:
            if _t10 in _f10_stats:
                _s = _f10_stats[_t10]
                _lbl = f"{_t10}\nμ={_s['mean']:.3f}  ±{_s['ci95']:.3f}"
                if abs(_s['med'] - _s['mean']) > 0.005:
                    _lbl += f"\nmed={_s['med']:.3f}"
                _lbl += f"\nn={_s['n']:,}"
            else:
                _lbl = _t10
            _xtlbl10.append(_lbl)
        ax.set_xticklabels(_xtlbl10, rotation=40, ha='right', fontsize=7.5)
        for _tick13x, _tier13x in zip(ax.get_xticklabels(), existing_tiers):
            _tick13x.set_color(TIER_PALETTE.get(_tier13x, 'black'))
            _tick13x.set_fontweight('bold')
        ax.yaxis.grid(True, color='#DCDCDC', linewidth=0.6, alpha=0.70, zorder=0)
        ax.set_axisbelow(True)
        _means10 = [_f10_stats[t]['mean'] for t in existing_tiers if t in _f10_stats]
        if _means10 and (max(_means10) - min(_means10)) < 0.05:
            _mean10_all = float(np.mean(_means10))
            ax.axhline(y=_mean10_all, color='#D4A000', linewidth=1.5,
                       linestyle='--', alpha=0.65, zorder=2)
            ax.text(0.99, 0.52,
                    f'⚡ All tier means within {max(_means10)-min(_means10):.3f} '
                    f'— minimal inter-tier separation\n'
                    f'(MFS cluster near {_mean10_all:.3f}; see also Fig 08 tied-metric row)',
                    transform=ax.transAxes, ha='right', va='center', fontsize=8,
                    color='#9A6B00', fontweight='bold', style='italic',
                    bbox=dict(boxstyle='round,pad=0.15', fc='#FFF8E7', ec='#D4A000',
                              alpha=0.90, linewidth=0.8), zorder=8)
        from matplotlib.lines import Line2D as _L10
        _lh10 = [
            _L10([0],[0], marker='o', color='none', markerfacecolor='#888888',
                 markeredgecolor='none', markersize=5, alpha=0.35,
                 label='Individual complex'),
            _L10([0],[0], marker='D', color='none', markerfacecolor='#888888',
                 markeredgecolor='black', markersize=8, markeredgewidth=1.2,
                 label='Tier mean (◆)'),
            _L10([0],[0], color='black', linewidth=2.0, label='95% confidence interval'),
        ]
        ax.legend(handles=_lh10, loc='lower left', fontsize=8.5,
                  framealpha=0.92, fancybox=True,
                  ncol=3).set_zorder(20)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_10_Mechanistic_Fingerprint_by_Tier.png", dpi=300, bbox_inches='tight')
        plt.close()
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_10_Mechanistic_Fingerprint_by_Tier.png').resolve()}")

    # --- Figure 11: SN2 Attack Angle — ECDF by Tier ---
    if 'SN2_Attack_Angle' in df.columns and 'degrader_tier' in df.columns:
        plot_df = df[df['SN2_Attack_Angle'] > 0].copy()
        valid_tiers_f14 = [t for t in existing_tiers if t in plot_df['degrader_tier'].values]
        fig, ax = plt.subplots(figsize=(12, 7))
        # Shade tier quality zones + count labels inside each zone band
        _zone_defs11 = [(175, 180, '#009E73'), (165, 175, '#56B4E9'),
                        (155, 165, '#0072B2'), (145, 155, '#CC79A7')]
        for lo, hi, col in _zone_defs11:
            ax.axvspan(lo, hi, alpha=0.07, color=col)
            ax.text((lo + hi) / 2, 1.025, f'≥{lo}°', ha='center', va='bottom',
                    fontsize=7.5, color=col, fontweight='bold', transform=ax.get_xaxis_transform())
            _zn11 = int(((plot_df['SN2_Attack_Angle'] >= lo) & (plot_df['SN2_Attack_Angle'] < hi)).sum())
            ax.text((lo + hi) / 2, 0.10, f'n={_zn11:,}',
                    ha='center', va='center', fontsize=6.0, color=col, fontweight='bold',
                    rotation=90, zorder=10,
                    transform=ax.get_xaxis_transform(),
                    bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.68))
        # ECDF per tier — annotations alternate above/below the scatter dot to avoid overlap
        _f14_annot_side = 1   # +1 = above dot, -1 = below dot; toggle each tier
        for tier in valid_tiers_f14:
            tier_angles = np.sort(plot_df.loc[plot_df['degrader_tier'] == tier, 'SN2_Attack_Angle'].values)
            ecdf_y = np.arange(1, len(tier_angles) + 1) / len(tier_angles)
            col = TIER_PALETTE.get(tier, '#999')
            ax.step(tier_angles, ecdf_y, where='post', color=col, linewidth=2.5,
                    label=f'{tier}  (n={len(tier_angles):,})')
            _n_ecdf = len(tier_angles)
            _eps_dkw = np.sqrt(np.log(2 / 0.05) / (2 * _n_ecdf))
            _ecdf_y_lo = np.clip(ecdf_y - _eps_dkw, 0, 1)
            _ecdf_y_hi = np.clip(ecdf_y + _eps_dkw, 0, 1)
            ax.fill_between(tier_angles, _ecdf_y_lo, _ecdf_y_hi,
                            color=TIER_PALETTE.get(tier, '#999'), alpha=0.15, zorder=1)
            med = float(np.median(tier_angles))
            med_y = float(np.interp(med, tier_angles, ecdf_y))
            ax.scatter([med], [med_y], color=col, s=70, zorder=6, edgecolors='black', linewidths=0.8)
            # Alternate label side: odd tiers above (+14 pts), even tiers below (-14 pts)
            _dy14 = 14 * _f14_annot_side
            _va14 = 'bottom' if _f14_annot_side > 0 else 'top'
            ax.annotate(f'{med:.0f}°', (med, med_y),
                        xytext=(4, _dy14), textcoords='offset points',
                        fontsize=8, color=col, fontweight='bold', va=_va14,
                        bbox=dict(boxstyle='round,pad=0.15', fc='white', ec=col,
                                  alpha=0.80, linewidth=0.6),
                        arrowprops=dict(arrowstyle='->', color=col, lw=0.7,
                                        shrinkA=0, shrinkB=3))
            _f14_annot_side *= -1   # flip side for next tier
        # Vertical threshold lines
        for angle, tier_key in [(175, 'Perfect_A'), (165, 'Perfect_B'), (155, 'Best_A'), (145, 'Best_B')]:
            ax.axvline(x=angle, color=TIER_PALETTE.get(tier_key, '#999'), linestyle='--', alpha=0.65, linewidth=1.0)
        ax.set_xlim(30, 185)
        ax.set_ylim(0, 1.09)
        ax.set_xticks([30, 45, 90, 130, 145, 155, 165, 175, 180])
        ax.set_yticks([i/10 for i in range(0, 11)])
        ax.set_yticklabels([f'{i*10}%' for i in range(0, 11)])
        ax.set_xlabel("SN2 Attack Angle (°)  — 180° = ideal linear nucleophilic back-attack", fontsize=11)
        ax.set_ylabel("Cumulative Fraction of Complexes in Tier", fontsize=11)
        # upper left — ECDF lines fan right so top-left is always clear
        from matplotlib.lines import Line2D as _L11
        _hdr11 = _L11([0],[0], color='none', linewidth=0, label='Tier (●=median)')
        _h11, _l11 = ax.get_legend_handles_labels()
        _leg14 = ax.legend(handles=[_hdr11] + _h11, labels=['Tier (●=median)'] + _l11,
                           loc='upper left', bbox_to_anchor=(0.01, 0.99),
                           fontsize=6, framealpha=0.92, fancybox=True,
                           ncol=max(2, (len(_h11) + 1) // 2))
        _leg14.set_zorder(20)
        # Thin horizontal grid lines at every 10% ECDF level for easy reading
        ax.set_axisbelow(True)
        ax.yaxis.grid(True, color='#DCDCDC', linewidth=0.6, linestyle='-', alpha=0.85, zorder=0)
        ax.xaxis.grid(True, color='#EBEBEB', linewidth=0.5, linestyle=':', alpha=0.7, zorder=0)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_11_SN2_Angle_by_Tier.png", dpi=300, bbox_inches='tight')
        plt.close()
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_11_SN2_Angle_by_Tier.png').resolve()}")

    # --- Figure 12a: SN2 Mechanism Geometry Scatter ---
    if 'Dist_Nucleophile_ASP110' in df.columns and 'SN2_Attack_Angle' in df.columns:
        _n_raw9 = len(df)
        plot_df12 = df[
            (df['Dist_Nucleophile_ASP110'] > 0) &
            (df['Dist_Nucleophile_ASP110'] < 5.0) &
            (df['SN2_Attack_Angle'] >= 0)
        ].copy()
        n_excluded12 = _n_raw9 - len(plot_df12)
        plot_df12['Stabilised'] = plot_df12.get('halide_stabilisation_score',
                                                 pd.Series(0.0, index=plot_df12.index)) > 0.5
        plot_df12['Status'] = plot_df12['Stabilised'].map({True: 'Aromatic Shield', False: 'Unstabilised'})

        fig, ax = plt.subplots(figsize=(12, 8))
        markers12 = {'Aromatic Shield': 'o', 'Unstabilised': 'X'}
        plot_df12_sorted = plot_df12.copy()
        plot_df12_sorted['_tier_ord'] = plot_df12_sorted['degrader_tier'].map(
            {t: i for i, t in enumerate(reversed(existing_tiers))})
        plot_df12_sorted = plot_df12_sorted.sort_values('_tier_ord')
        sns.scatterplot(data=plot_df12_sorted, x='Dist_Nucleophile_ASP110', y='SN2_Attack_Angle',
                        hue='degrader_tier', hue_order=existing_tiers, style='Status',
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
                _td12 = plot_df12[plot_df12['degrader_tier'] == tier]
                if len(_td12) < 30:
                    continue
                _col12 = TIER_PALETTE.get(tier, '#999')
                try:
                    _kde12 = _gkde12a(np.vstack([
                        _td12['Dist_Nucleophile_ASP110'].values,
                        _td12['SN2_Attack_Angle'].values
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

        # Thresholds sourced DIRECTLY from CFG (not hardcoded) so the guide lines
        # exactly match the production tier gates in 02_Production.
        _nuc_cfg12 = CFG.TIER_NUC_DIST
        _ang_cfg12 = CFG.TIER_ANGLE_MIN
        ax.axvspan(0, _nuc_cfg12['Perfect_A'], alpha=0.06, color='#009E73', label='_nolegend_')
        ax.axhspan(_ang_cfg12['Perfect_A'], 180, alpha=0.06, color='#009E73', label='_nolegend_')

        for _idx12, tier_key in enumerate(['Perfect_A', 'Perfect_B', 'Best_A', 'Best_B']):
            dist = _nuc_cfg12[tier_key]
            col = TIER_PALETTE.get(tier_key, '#999')
            ax.axvline(x=dist, color=col, linestyle='--', alpha=0.65, linewidth=1.2)
            ax.text(dist, 186 if _idx12 % 2 == 0 else 181, f'{tier_key}\n≤{dist}Å', color=col,
                    fontsize=6.5, va='bottom', ha='center', alpha=0.90,
                    bbox=dict(boxstyle='round,pad=0.15', fc='white', alpha=0.75, ec=col,
                              linewidth=0.5))

        for tier_key in ['Perfect_A', 'Perfect_B', 'Best_A', 'Best_B']:
            angle = _ang_cfg12[tier_key]
            col = TIER_PALETTE.get(tier_key, '#999')
            ax.axhline(y=angle, color=col, linestyle=':', alpha=0.65, linewidth=1.2, zorder=1)
            ax.text(0.02, angle, f'≥{angle}° ({tier_key})', color=col,
                    fontsize=6.5, ha='left', va='bottom', alpha=0.90, zorder=6,
                    transform=ax.get_yaxis_transform(),
                    bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none',
                              alpha=0.70, linewidth=0))

        ax.set_xlim(0, 5.0)
        ax.set_ylim(0, 195)
        ax.set_yticks(range(0, 196, 20))
        ax.set_xlabel("ASP110 Nucleophile → Carbon Distance (Å)  — shorter = closer to reaction geometry",
                      fontsize=10)
        ax.set_ylabel("SN2 Attack Angle (°)  — 180° = perfect linear back-attack", fontsize=10)
        ax.text(0.01, 0.01,
                f'n = {len(plot_df12):,} complexes shown  |  {n_excluded12:,} excluded (dist = 0 or dist ≥ 5 Å)',
                transform=ax.transAxes, fontsize=7.5, color='#666', style='italic', va='bottom',
                bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#ccc', alpha=0.82))
        handles12, labels12 = ax.get_legend_handles_labels()
        clean_labels12, clean_handles12 = [], []
        for h12, l12 in zip(handles12, labels12):
            if l12 in ('degrader_tier', 'Status'):
                continue
            clean_handles12.append(h12)
            clean_labels12.append(l12)
        ax.legend(clean_handles12, clean_labels12,
                  loc='upper left', bbox_to_anchor=(0.01, 0.99),
                  fontsize=7, framealpha=0.92, fancybox=True,
                  ncol=max(2, (len(clean_handles12) + 1) // 2),
                  handlelength=1.2, handletextpad=0.4)
        ax.set_axisbelow(True)
        ax.grid(color='#E0E0E0', linewidth=0.5, alpha=0.7)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_12a_Mechanism_Geometry_Scatter.png",
                    dpi=300, bbox_inches='tight')
        plt.close()
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_12a_Mechanism_Geometry_Scatter.png').resolve()}")

    # Figure 12b: PA companion (mechanistic scatter overlay)
    if _pa_has_imgs:
        _fig_12b_pa_mechanistic(df, _pa, _imgs, out_dir, reporter)

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
    # --- Figure 13: Molecular Interaction Profile — Stacked bars (bond types) + F-engagement line ---
    # Single chart, dual y-axis: stacked bars per tier show the total interaction count and
    # its bond-type composition (left y-axis); overlaid connected dot-line shows the fluorine
    # engagement ratio per tier on the right y-axis, revealing whether richer interaction
    # profiles correlate with higher fluorine utilisation.
    _has_int   = bool(present_int_cols) and 'degrader_tier' in df.columns
    _has_f16   = ('interacting_fluorine_count' in df.columns and
                  'total_fluorine_count' in df.columns and 'degrader_tier' in df.columns)
    if _has_int or _has_f16:
        try:
            from scipy import stats as _sc_stats15
            # Pre-compute tier count for dynamic figure height
            n_tiers_f13 = len([t for t in existing_tiers if t in df['degrader_tier'].values]) if 'degrader_tier' in df.columns else 6
            fig, (ax, ax_hm) = plt.subplots(1, 2, figsize=(20, max(5.5, n_tiers_f13 * 0.95 + 2.5)),
                                             gridspec_kw={'width_ratios': [3, 2], 'wspace': 0.35})
            ax_r15 = ax.twinx()   # right y-axis for fluorine engagement line

            # Left axis: stacked bar chart — one bar per tier, stacked by bond type
            if _has_int:
                tier_int = df.groupby('degrader_tier')[[present_int_cols[k] for k in present_int_cols]].mean()
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
                    bars = ax.bar(bar_x, vals_pct, bottom=bar_bottom,
                                  color=bond_palette.get(bond_type, '#999'),
                                  edgecolor='white', linewidth=0.5, width=0.55,
                                  label=bond_type, zorder=2)
                    # Annotate segment with absolute mean count; skip tiny segments
                    for xi, vi, va, bi in zip(bar_x, vals_pct, vals_abs, bar_bottom):
                        if vi < 3.0:
                            continue
                        mid_y = bi + vi / 2
                        _fs = 7 if vi >= 8.0 else 5.5
                        ax.text(xi, mid_y, f'{va:.1f}',
                                ha='center', va='center', fontsize=_fs,
                                color='white', fontweight='bold', zorder=3)
                    bar_bottom = bar_bottom + vals_pct
                # Total mean count label on top (absolute, not the 100% bar height)
                for xi, tot_abs in zip(bar_x, _row_totals_abs.values):
                    ax.text(xi, 102.0, f'Σ={tot_abs:.0f}',
                            ha='center', va='bottom', fontsize=8, fontweight='bold', color='black', zorder=4)
                ax.set_xticks(bar_x)
                ax.set_xticklabels(tier_labels, rotation=35, ha='right', fontsize=9)
                ax.set_ylim(0, 115)
                ax.set_yticks([0, 25, 50, 75, 100])
                ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=9)
                # Left axis: steel-blue tick labels matching the blue gridlines
                ax.tick_params(axis='y', labelsize=9, labelcolor='#2C6FAC')
                ax.yaxis.label.set_color('#2C6FAC')
                ax.set_xlabel("Degrader Tier", fontsize=11)
                ax.set_ylabel("Interaction type proportion  (%)  —  annotation = mean count", fontsize=11)

            # Right axis: Fluorine Engagement Ratio connected line per tier
            if _has_f16:
                f16_df = df.copy()
                f16_df['FER'] = np.where(
                    f16_df['total_fluorine_count'] > 0,
                    f16_df['interacting_fluorine_count'] / f16_df['total_fluorine_count'],
                    np.nan
                )
                f16_df = f16_df.dropna(subset=['FER'])
                fer_xs, fer_meds, fer_lo_ci, fer_hi_ci = [], [], [], []
                for i, tier in enumerate(existing_tiers):
                    if tier not in f16_df['degrader_tier'].values:
                        continue
                    vals = f16_df.loc[f16_df['degrader_tier'] == tier, 'FER'].dropna()
                    if len(vals) < 2:
                        continue
                    med_v = float(vals.median())
                    ci95  = float(_sc_stats15.sem(vals) * _sc_stats15.t.ppf(0.975, len(vals) - 1))
                    fer_xs.append(i)
                    fer_meds.append(med_v)
                    fer_lo_ci.append(max(0, med_v - ci95))
                    fer_hi_ci.append(min(1, med_v + ci95))
                if fer_xs:
                    ax_r15.plot(fer_xs, fer_meds, color='#E69F00', linewidth=2.6,
                                marker='o', markersize=9, markeredgecolor='black',
                                markeredgewidth=0.9, zorder=8,
                                label='Fluorine Engagement Ratio  (median ± 95% CI)')
                    ax_r15.fill_between(fer_xs, fer_lo_ci, fer_hi_ci,
                                        color='#E69F00', alpha=0.22, zorder=7)
                    # Filled zone backgrounds on right axis
                    ax_r15.axhspan(0.75, 1.20, alpha=0.15, color='#009E73', zorder=0)
                    ax_r15.axhspan(0.50, 0.75, alpha=0.13, color='#0072B2', zorder=0)
                    ax_r15.axhspan(0.00, 0.50, alpha=0.12, color='#D55E00', zorder=0)
                    ax_r15.axhline(y=0.5, color='#0072B2', linestyle='--',
                                   alpha=0.70, linewidth=1.3, zorder=6)
                    ax_r15.axhline(y=1.0, color='#009E73', linestyle=':',
                                   alpha=0.75, linewidth=1.3, zorder=6)
                    # FER value labels: all placed in a fixed band above the right axis
                    # (y = 1.08–1.20) to clear the stacked bars.  Cycle through 3 y-levels
                    # and alternate x-offset ±0.15 so adjacent labels never touch.
                    for _fi15, (xi, yi) in enumerate(zip(fer_xs, fer_meds)):
                        _lbl_y15 = 1.08 + (_fi15 % 3) * 0.05   # 1.08 / 1.13 / 1.18
                        _x_off15 = 0.12 * (1 if _fi15 % 2 == 0 else -1)
                        ax_r15.annotate(f'{yi:.2f}',
                                        xy=(xi, yi), xycoords='data',
                                        xytext=(xi + _x_off15, _lbl_y15), textcoords='data',
                                        ha='center', va='bottom', fontsize=7.5,
                                        color='#A06000', fontweight='bold', zorder=10,
                                        arrowprops=dict(arrowstyle='->', color='#E69F00',
                                                        lw=0.5, mutation_scale=6,
                                                        shrinkA=0, shrinkB=2),
                                        bbox=dict(boxstyle='round,pad=0.09', fc='white',
                                                  ec='#E69F00', alpha=0.94, linewidth=0.5))
                ax_r15.set_ylim(0, 1.30)
                ax_r15.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
                ax_r15.set_yticklabels(['0%', '25%', '50%', '75%', '100%'], fontsize=9)
                ax_r15.set_ylabel("Fluorine Engagement Ratio  (interacting F / total F per ligand)",
                                  fontsize=11, rotation=270, labelpad=14, color='#A06000')
                ax_r15.tick_params(axis='y', labelcolor='#A06000', labelsize=9)
                # Rotated zone labels — left side of figure, written bottom-to-top
                # ylim=(0, 1.30): Low midpoint=0.25→0.192; Mid midpoint=0.625→0.481; High midpoint=0.875→0.673
                ax_r15.text(0.99, 0.192, 'Low engagement (<50%)',
                            ha='center', va='center', rotation=90,
                            fontsize=7.5, color='#D55E00', style='italic', fontweight='bold',
                            transform=ax_r15.transAxes,
                            bbox=dict(boxstyle='round,pad=0.10', fc='white', ec='#D55E00',
                                      alpha=0.80, linewidth=0.5))
                ax_r15.text(0.99, 0.481, '50% engagement',
                            ha='center', va='center', rotation=90,
                            fontsize=7.5, color='#0072B2', style='italic', fontweight='bold',
                            transform=ax_r15.transAxes,
                            bbox=dict(boxstyle='round,pad=0.10', fc='white', ec='#0072B2',
                                      alpha=0.80, linewidth=0.5))
                ax_r15.text(0.99, 0.673, 'Full engagement',
                            ha='center', va='center', rotation=90,
                            fontsize=7.5, color='#009E73', style='italic', fontweight='bold',
                            transform=ax_r15.transAxes,
                            bbox=dict(boxstyle='round,pad=0.10', fc='white', ec='#009E73',
                                      alpha=0.80, linewidth=0.5))

            # Dual-colour grid system: left = steel-blue (matches left tick labels),
            # right = amber (matches right tick labels). Both BEHIND bars via set_axisbelow.
            ax.set_axisbelow(True)
            ax.yaxis.grid(True, color='#2C6FAC', linewidth=0.55, linestyle='-',
                          alpha=0.20, zorder=0)
            ax.xaxis.grid(False)
            ax_r15.yaxis.grid(True, color='#A06000', linewidth=0.55, linestyle='--',
                              alpha=0.18, zorder=0)
            ax_r15.set_axisbelow(True)

            # Split legend: bond-type bars in lower-left (typically sparse area),
            # FER line separately via ax_r15 annotation to avoid overlap with bars.
            handles_bars, labels_bars = ax.get_legend_handles_labels()
            handles_line, labels_line = ax_r15.get_legend_handles_labels()
            _all_h15 = handles_bars + handles_line
            _all_l15 = labels_bars + labels_line
            # Single row across the top; if too many items, allow a second row (ncol = ceil/2)
            import math as _math15
            _ncol15 = len(_all_h15) if len(_all_h15) <= 8 else _math15.ceil(len(_all_h15) / 2)
            _leg15 = ax.legend(_all_h15, _all_l15,
                               loc='upper left', ncol=_ncol15,
                               fontsize=7.5, framealpha=0.93, fancybox=True,
                               borderpad=0.4, labelspacing=0.3,
                               handlelength=1.0, handletextpad=0.4,
                               columnspacing=0.8)
            _leg15.set_zorder(20)

            # Heatmap: tier × bond-type proportion (%) — only when interaction data exist
            if _has_int and 'tier_int' in dir():
                import seaborn as _sns13
                _hm_data = tier_int.T  # bond types as rows, tiers as columns
                # tier_int is already normalised to 100% (from earlier in the function)
                _sns13.heatmap(_hm_data, ax=ax_hm,
                               cmap='YlOrRd', vmin=0, vmax=100,
                               annot=True, fmt='.0f', annot_kws={'size': 8},
                               linewidths=0.4, linecolor='white',
                               cbar_kws={'label': 'Proportion (%)', 'shrink': 0.8})
                ax_hm.set_xlabel("Degrader Tier  —  Interaction Type Proportion (%)", fontsize=9.5)
                ax_hm.set_ylabel("Interaction Type", fontsize=10)
                ax_hm.tick_params(axis='x', labelrotation=45, labelsize=8.5)
                ax_hm.tick_params(axis='y', labelrotation=0,  labelsize=8)
            else:
                ax_hm.set_visible(False)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_13a_Molecular_Interaction_Profile.png",
                        dpi=300, bbox_inches='tight')
            plt.close()
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_13a_Molecular_Interaction_Profile.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 13 skipped: {e}")

    # Figure 13b: PA companion (interaction profile overlay)
    if _pa_has_imgs:
        _fig_13b_pa_interactions(df, _pa, _imgs, out_dir, reporter)

    # --- Figure 14: Fluorine Engagement Ratio by Tier (box + strip + median trend line) ---
    if 'interacting_fluorine_count' in df.columns and 'total_fluorine_count' in df.columns and 'degrader_tier' in df.columns:
        try:
            f16_df = df.copy()
            f16_df['FER'] = np.where(
                f16_df['total_fluorine_count'] > 0,
                f16_df['interacting_fluorine_count'] / f16_df['total_fluorine_count'],
                np.nan)
            f16_df = f16_df.dropna(subset=['FER'])
            valid_t16 = [t for t in existing_tiers if t in f16_df['degrader_tier'].values]
            fig, ax = plt.subplots(figsize=(12, 7))

            ax.axhspan(0.75, 1.01, alpha=0.07, color='#009E73', zorder=0)
            ax.axhspan(0.50, 0.75, alpha=0.06, color='#E69F00', zorder=0)
            ax.axhspan(0.00, 0.50, alpha=0.06, color='#D55E00', zorder=0)
            ax.axhline(y=0.75, color='#009E73', linestyle='--', alpha=0.6, linewidth=1.1)
            ax.axhline(y=0.50, color='#E69F00', linestyle=':', alpha=0.6, linewidth=1.1)
            ax.set_axisbelow(True)
            # KDE violin + IQR box overlay (no boxplot — keeps figure clean)
            import scipy.stats as _stats14
            for _ti14, tier14 in enumerate(valid_t16):
                _vals14 = f16_df.loc[f16_df['degrader_tier'] == tier14, 'FER'].dropna().values
                if len(_vals14) < 4:
                    continue
                from scipy.stats import gaussian_kde as _gkde14
                _kde14  = _gkde14(_vals14)
                _y14    = np.linspace(np.min(_vals14), np.max(_vals14), 200)
                _w14    = _kde14(_y14)
                _w14   /= _w14.max()
                _w14   *= 0.35   # half-width of violin
                ax.fill_betweenx(_y14, _ti14 - _w14, _ti14 + _w14,
                                 alpha=0.55, color=TIER_PALETTE.get(tier14, '#999'))
                # IQR box on top
                _q25, _q75 = np.percentile(_vals14, [25, 75])
                _med14 = np.median(_vals14)
                ax.vlines(_ti14, _q25, _q75, color='#333', linewidth=2.5, zorder=5)
                ax.scatter([_ti14], [_med14], color='white', s=30, zorder=6,
                           edgecolors='#333', linewidths=1.2)

            # Per-tier median statistics badge
            _med16_xs, _med16_ys = [], []
            for i, tier in enumerate(valid_t16):
                sub16 = f16_df.loc[f16_df['degrader_tier'] == tier, 'FER'].dropna()
                if len(sub16) == 0:
                    continue
                med = float(sub16.median())
                n16 = len(sub16)
                _med16_xs.append(i)
                _med16_ys.append(med)
            # FER connecting line across tier medians
            if len(_med16_xs) >= 2:
                ax.plot(_med16_xs, _med16_ys, color='#E69F00', linewidth=2.2,
                        linestyle='--', alpha=0.85, zorder=4,
                        marker='D', markersize=6, markeredgecolor='#333',
                        markeredgewidth=0.7, label='FER trend (median)')
            for i, tier in enumerate(valid_t16):
                sub16 = f16_df.loc[f16_df['degrader_tier'] == tier, 'FER'].dropna()
                if len(sub16) == 0:
                    continue
                med = float(sub16.median())
                n16 = len(sub16)
                ax.text(i, med + 0.022, f'med={med:.2f}\nn={n16:,}',
                        ha='center', va='bottom', fontsize=6.0, fontweight='bold',
                        zorder=7, color='#111',
                        bbox=dict(boxstyle='round,pad=0.2', fc='white', ec='#ccc',
                                  linewidth=0.6, alpha=0.85))

            ax.set_ylim(0, 1.08)
            ax.set_yticks([0, 0.25, 0.50, 0.75, 1.00])
            ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
            # Color y-tick labels to match zone backgrounds
            _ytick_colors16 = {
                '0%':   '#B84000',
                '25%':  '#B84000',
                '50%':  '#8A6000',
                '75%':  '#007A50',
                '100%': '#007A50',
            }
            for tick16 in ax.get_yticklabels():
                tick16.set_color(_ytick_colors16.get(tick16.get_text(), '#111'))
                tick16.set_fontweight('bold')
            # Start-y (axes fraction) = bottom border of each label's OWN zone, so
            # every label sits on its category's lower boundary line (va='bottom',
            # text rises into its zone). Low's zone bottom coincides with the axes
            # floor, which is expected — High/Moderate now match that anchoring.
            for (_zy14_start, _ztxt14, _zcol14, _zec14) in [
                (0.73, 'High engagement (≥75%)', '#007A50', '#009E73'),
                (0.50, 'Moderate (50–75%)',       '#8A6000', '#E69F00'),
                (0.10, 'Low engagement (<50%)',   '#B84000', '#D55E00'),
            ]:
                ax.text(0.985, _zy14_start, _ztxt14,
                        ha='center', va='bottom', clip_on=True,
                        fontsize=5.5, color=_zcol14, style='italic',
                        fontweight='bold', rotation=90, transform=ax.transAxes,
                        bbox=dict(boxstyle='round,pad=0.12', fc='white', ec=_zec14,
                                  alpha=0.80, linewidth=0.5))
            ax.set_xlabel("Degrader Tier", fontsize=11)
            ax.set_ylabel("Fluorine Engagement Ratio  (interacting F / total F)", fontsize=11)
            ax.set_xticks(range(len(valid_t16)))
            ax.set_xticklabels(valid_t16, rotation=35, ha='right', fontsize=9)
            for tick16x, tier16x in zip(ax.get_xticklabels(), valid_t16):
                tick16x.set_color(TIER_PALETTE.get(tier16x, 'black'))
                tick16x.set_fontweight('bold')

            # Legend — violin + IQR box description
            from matplotlib.patches import Patch as _P16
            from matplotlib.lines import Line2D as _L16
            _leg16_h = [
                _P16(facecolor='#888', alpha=0.55, label='KDE violin (per tier)'),
                _L16([0], [0], color='#333', linewidth=2.5, label='IQR (25–75%)'),
                _L16([0], [0], color='white', marker='o', markersize=6,
                     markeredgecolor='#333', markeredgewidth=1.2,
                     linestyle='None', label='Median'),
            ]
            _leg16 = ax.legend(handles=_leg16_h,
                               loc='lower left', bbox_to_anchor=(0.01, 0.01),
                               fontsize=6, framealpha=0.92, fancybox=True, ncol=3,
                               labelspacing=0.2, handlelength=0.8, handletextpad=0.3)
            _leg16.set_zorder(20)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_14_Fluorine_Engagement_by_Tier.png", dpi=300, bbox_inches='tight')
            plt.close()
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_14_Fluorine_Engagement_by_Tier.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 14 skipped: {e}")

    # --- Figure 15: Binding Energetics — Binding Probability (violin) + Product Inhibition (line) ---
    # Single chart, dual y-axis: violin plots show the distribution of Binding Probability
    # per tier on the left axis (zoomed to ~0.97-1.0 since all values cluster there);
    # a connected dot-line shows Product Inhibition Penalty per tier on the right axis.
    # Together they reveal whether the tiers that bind well also carry high inhibition risk.
    bind_col   = next((c for c in ["Binding_Probability", "Likelihood_Degrader_Score"] if c in df.columns), None)
    _has_bind  = bind_col is not None and 'degrader_tier' in df.columns
    _has_pip   = 'product_inhibition_penalty' in df.columns and 'degrader_tier' in df.columns
    if _has_bind or _has_pip:
        try:
            from scipy import stats as _sc_stats12
            fig, ax = plt.subplots(figsize=(13, 7))
            ax_r = ax.twinx()   # right y-axis for product inhibition line

            # Left axis: Binding Probability violin per tier
            if _has_bind:
                f12_data = df.dropna(subset=[bind_col])
                y_lo_bp = max(0.0, float(f12_data[bind_col].quantile(0.01)) - 0.003)
                y_hi_bp = min(1.0, float(f12_data[bind_col].max()) + 0.001)
                sns.violinplot(data=f12_data, x='degrader_tier', y=bind_col, order=existing_tiers,
                               palette=TIER_PALETTE, inner='box', linewidth=1.2, cut=0,
                               alpha=0.75, ax=ax)
                for i, tier in enumerate(existing_tiers):
                    med = f12_data.loc[f12_data['degrader_tier'] == tier, bind_col].median()
                    if not np.isnan(med):
                        # Red median line across the violin width
                        ax.hlines(med, i - 0.22, i + 0.22, colors='#CC0000',
                                  linewidth=2.0, zorder=7, linestyle='-')
                        ax.text(i + 0.25, med, f'  {med:.4f}', va='center', ha='left',
                                fontsize=7, color='#CC0000', fontweight='bold', zorder=8)
                ax.set_ylim(y_lo_bp, y_hi_bp)
                ax.set_ylabel("Binding Probability Score  (sigmoid of SN2 geometry)", fontsize=11)
            else:
                ax.set_visible(False)

            # Right axis: Product Inhibition Penalty connected dot-line per tier
            if _has_pip:
                f20_df = df.dropna(subset=['product_inhibition_penalty']).copy()
                valid_tiers_f20 = [t for t in existing_tiers if t in f20_df['degrader_tier'].values]
                pip_xs, pip_meds, pip_lo_ci, pip_hi_ci = [], [], [], []
                for tier in valid_tiers_f20:
                    x_pos = existing_tiers.index(tier) if tier in existing_tiers else 0
                    vals  = f20_df.loc[f20_df['degrader_tier'] == tier, 'product_inhibition_penalty'].dropna()
                    if len(vals) < 2:
                        continue
                    med_v = float(vals.median())
                    ci95  = float(_sc_stats12.sem(vals) * _sc_stats12.t.ppf(0.975, len(vals) - 1))
                    pip_xs.append(x_pos)
                    pip_meds.append(med_v)
                    pip_lo_ci.append(med_v - ci95)
                    pip_hi_ci.append(med_v + ci95)
                if pip_xs:
                    ax_r.plot(pip_xs, pip_meds, color='#D55E00', linewidth=2.5,
                              marker='D', markersize=9, markeredgecolor='black',
                              markeredgewidth=0.9, zorder=7, label='Product Inhibition Penalty\n(median ± 95 % CI)')
                    ax_r.fill_between(pip_xs, pip_lo_ci, pip_hi_ci, color='#D55E00',
                                      alpha=0.20, zorder=6)
                    # Median annotation
                    for xi, yi in zip(pip_xs, pip_meds):
                        ax_r.text(xi, yi + (max(pip_meds) - min(pip_meds)) * 0.04,
                                  f'{yi:.3f}', ha='center', va='bottom', fontsize=8,
                                  color='#D55E00', fontweight='bold', zorder=8)
                ax_r.set_ylabel("Product Inhibition Penalty Score  (higher = greater feedback risk)",
                                fontsize=11, rotation=270, labelpad=14, color='#D55E00')
                ax_r.tick_params(axis='y', labelcolor='#D55E00', labelsize=9)

            ax.set_xticks(range(len(existing_tiers)))
            ax.set_xticklabels(existing_tiers, rotation=40, ha='right', fontsize=9)
            for _tick12x, _tier12x in zip(ax.get_xticklabels(), existing_tiers):
                _tick12x.set_color(TIER_PALETTE.get(_tier12x, 'black'))
                _tick12x.set_fontweight('bold')
            # Left axis (binding prob.): steel-blue labels matching blue gridlines
            ax.tick_params(axis='y', labelsize=9, labelcolor='#2C6FAC')
            ax.yaxis.label.set_color('#2C6FAC')
            # Right axis (PIP): amber labels matching amber gridlines
            ax_r.tick_params(axis='y', labelsize=9, labelcolor='#A06000')
            ax_r.yaxis.label.set_color('#A06000')
            ax.set_xlabel("Degrader Tier", fontsize=11)

            # Dual-colour grids: left-axis = steel-blue; right-axis = amber
            ax.set_axisbelow(True)
            ax.yaxis.grid(True, color='#2C6FAC', linewidth=0.55, linestyle='-', alpha=0.20, zorder=0)
            ax.xaxis.grid(False)
            ax_r.yaxis.grid(True, color='#A06000', linewidth=0.55, linestyle='--', alpha=0.18, zorder=0)
            ax_r.set_axisbelow(True)

            # Legend — bottom-left (violins tend to be taller on right)
            if _has_pip and pip_xs:
                from matplotlib.lines import Line2D as _Line12
                _pip_handle = _Line12([0], [0], color='#D55E00', linewidth=2.5,
                                      marker='D', markersize=8,
                                      label='Product Inhibition Penalty  (median ± 95% CI)')
                _leg12 = ax_r.legend(handles=[_pip_handle], loc='lower left',
                                     fontsize=8.5, framealpha=0.92, fancybox=True)
                _leg12.set_zorder(20)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_15_Binding_Energetics.png", dpi=300, bbox_inches='tight')
            plt.close()
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_15_Binding_Energetics.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 15 skipped: {e}")

    # --- Figure 16: Product Inhibition Penalty by Tier ---
    if 'product_inhibition_penalty' in df.columns and 'degrader_tier' in df.columns:
        try:
            f20_df = df.dropna(subset=['product_inhibition_penalty']).copy()
            f20_df = f20_df[f20_df['product_inhibition_penalty'] > 0]
            valid_t20 = [t for t in existing_tiers if t in f20_df['degrader_tier'].values]
            if not f20_df.empty and valid_t20:
                fig, ax = plt.subplots(figsize=(12, 6))
                # Use data-driven zone thresholds so zones span the full data range
                _p_all = f20_df['product_inhibition_penalty']
                _z_lo  = float(_p_all.quantile(0.25))   # lower 25% = low risk
                _z_mid = float(_p_all.quantile(0.75))   # upper 25% = high risk
                _z_max = float(_p_all.max()) * 1.10
                ax.axhspan(0,      _z_lo,  alpha=0.18, color='#009E73', zorder=0)
                ax.axhspan(_z_lo,  _z_mid, alpha=0.16, color='#E69F00', zorder=0)
                ax.axhspan(_z_mid, _z_max, alpha=0.14, color='#D55E00', zorder=0)
                ax.axhline(y=_z_lo,  color='#009E73', linestyle='--', alpha=0.6, linewidth=1.1)
                ax.axhline(y=_z_mid, color='#D55E00', linestyle=':', alpha=0.6, linewidth=1.1)
                sns.boxplot(data=f20_df, x='degrader_tier', y='product_inhibition_penalty',
                            order=valid_t20, palette=TIER_PALETTE, linewidth=1.2,
                            flierprops=dict(marker='o', markersize=2, alpha=0.3), ax=ax)
                sns.stripplot(data=f20_df, x='degrader_tier', y='product_inhibition_penalty',
                              order=valid_t20, color='black', alpha=0.10, size=2.0, jitter=True, ax=ax)
                for i, tier in enumerate(valid_t20):
                    sub20 = f20_df.loc[f20_df['degrader_tier'] == tier, 'product_inhibition_penalty'].dropna()
                    if len(sub20) == 0:
                        continue
                    med = float(sub20.median())
                    _bg_col20 = TIER_PALETTE.get(tier, '#666')
                    # Adaptive text: white on dark backgrounds, near-black on light ones
                    _txt_col20 = _text_color(_bg_col20)
                    ax.text(i, med, f'{med:.1f}', ha='center', va='center',
                            fontsize=8.5, fontweight='bold', zorder=7,
                            color=_txt_col20,
                            bbox=dict(boxstyle='round,pad=0.15', fc=_bg_col20,
                                      ec='white', linewidth=0.6, alpha=0.88))
                ax.set_xlabel("Degrader Tier", fontsize=11)
                ax.set_ylabel("Product Inhibition Penalty  (higher = greater feedback risk)", fontsize=11)
                ax.set_xticks(range(len(valid_t20)))
                ax.set_xticklabels(valid_t20, rotation=35, ha='right', fontsize=9)
                for _tick20x, _tier20x in zip(ax.get_xticklabels(), valid_t20):
                    _tick20x.set_color(TIER_PALETTE.get(_tier20x, 'black'))
                    _tick20x.set_fontweight('bold')
                ax.set_axisbelow(True)
                # Dual-colour grid: left-axis (penalty) = steel-blue, standard gridlines
                ax.yaxis.grid(True, color='#2C6FAC', linewidth=0.55, linestyle='-',
                              alpha=0.22, zorder=0)
                ax.xaxis.grid(False)
                ax.tick_params(axis='y', labelcolor='#2C6FAC', labelsize=9)
                ax.yaxis.label.set_color('#2C6FAC')
                # Zone backgrounds as legend entries (patch legend, upper right)
                _ymax20 = _z_max
                ax.set_ylim(bottom=0, top=_ymax20)
                from matplotlib.patches import Patch as _Patch20
                _zone_handles20 = [
                    _Patch20(facecolor='#009E73', alpha=0.45, edgecolor='#007A50',
                             linewidth=0.8, label=f'Low risk  (< {_z_lo:.0f})'),
                    _Patch20(facecolor='#E69F00', alpha=0.45, edgecolor='#b07000',
                             linewidth=0.8, label=f'Moderate  ({_z_lo:.0f}–{_z_mid:.0f})'),
                    _Patch20(facecolor='#D55E00', alpha=0.45, edgecolor='#b03000',
                             linewidth=0.8, label=f'High risk  (> {_z_mid:.0f})'),
                ]
                _leg20 = ax.legend(handles=_zone_handles20, loc='upper right',
                                   ncol=3, fontsize=8, framealpha=0.92, fancybox=True,
                                   title='Inhibition Risk Zone', title_fontsize=7.5)
                _leg20.set_zorder(20)

                # ── Relative tier comparison annotations ──────────────────────
                # Labels show best/worst tier name + median; no absolute risk qualifier
                # (both tiers can be in the HIGH zone, so "lowest risk" would be misleading)
                _best_t20 = valid_t20[0]
                _best_med20 = float(f20_df.loc[f20_df['degrader_tier'] == _best_t20,
                                               'product_inhibition_penalty'].median())
                ax.annotate(f'Best tier ({_best_t20})\nmedian = {_best_med20:.0f}',
                            xy=(0, _best_med20),
                            xytext=(0, _best_med20 + _ymax20 * 0.22),
                            ha='center', va='bottom', fontsize=8, color='#007A50',
                            fontweight='bold', style='italic',
                            arrowprops=dict(arrowstyle='->', color='#007A50',
                                            lw=1.4, connectionstyle='arc3,rad=0.25'),
                            zorder=8)
                _worst_t20 = valid_t20[-1]
                _worst_med20 = float(f20_df.loc[f20_df['degrader_tier'] == _worst_t20,
                                                'product_inhibition_penalty'].median())
                _x_worst20 = len(valid_t20) - 1
                ax.annotate(f'Worst tier ({_worst_t20})\nmedian = {_worst_med20:.0f}',
                            xy=(_x_worst20, _worst_med20),
                            xytext=(_x_worst20, _worst_med20 + _ymax20 * 0.22),
                            ha='center', va='bottom', fontsize=8, color='#A03000',
                            fontweight='bold', style='italic',
                            arrowprops=dict(arrowstyle='->', color='#A03000',
                                            lw=1.4, connectionstyle='arc3,rad=-0.25'),
                            zorder=8)
                # Kruskal-Wallis significance note
                try:
                    from scipy.stats import kruskal as _kw20
                    _kw_groups = [f20_df.loc[f20_df['degrader_tier'] == t,
                                             'product_inhibition_penalty'].dropna().values
                                  for t in valid_t20 if len(f20_df[f20_df['degrader_tier']==t]) >= 3]
                    if len(_kw_groups) >= 2:
                        _, _kw_p = _kw20(*_kw_groups)
                        _pstr = ('p<0.001' if _kw_p < 0.001 else
                                 f'p={_kw_p:.3f}')
                        ax.text(0.01, 0.99, f'Kruskal–Wallis {_pstr}  (across tiers)',
                                transform=ax.transAxes, ha='left', va='top',
                                fontsize=8, color='#333', style='italic')
                except Exception:
                    pass

                # Per-tier coverage: fraction of tier's total rows that have this metric
                # product_inhibition_penalty only present when interaction_density > 1.5
                # → strong selection bias; must annotate to prevent misinterpretation
                _n_total = len(df)
                _n_metric = len(f20_df)
                _global_cov = 100 * _n_metric / max(_n_total, 1)
                _cov16_parts = []
                for i, tier in enumerate(valid_t20):
                    _n_tier_all = int((df['degrader_tier'] == tier).sum())
                    _n_tier_pip = int(((f20_df['degrader_tier'] == tier)).sum())
                    _cov_pct = 100 * _n_tier_pip / max(_n_tier_all, 1)
                    ax.text(i, -_ymax20 * 0.055, f'{_n_tier_pip}/{_n_tier_all}\n({_cov_pct:.0f}%)',
                            ha='center', va='top', fontsize=6.5, color='#555',
                            style='italic',
                            bbox=dict(boxstyle='round,pad=0.08', fc='#FFF8E1',
                                      ec='#E69F00', alpha=0.80, linewidth=0.5))
                    _cov16_parts.append(f"{tier}:{_n_tier_pip}/{_n_tier_all}({_cov_pct:.0f}%)")
                pass  # Fig16 cov suppressed
                ax.set_ylim(bottom=-_ymax20 * 0.14, top=_ymax20)

                plt.tight_layout()
                plt.savefig(out_dir / "Figure_16_Product_Inhibition_by_Tier.png", dpi=300, bbox_inches='tight')
                plt.close()
                reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_16_Product_Inhibition_by_Tier.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 16 skipped: {e}")

    # --- Figure 17a: Chemical Space UMAP Manifold ---
    if 'UMAP_X' in df.columns:
        fig, ax = plt.subplots(figsize=(10, 7))
        # Local colour overrides so Perfect_B and Best_A are visually distinct
        # even if TIER_PALETTE assigns them similar shades.
        _f05_palette = dict(TIER_PALETTE)
        _f05_palette['Perfect_B'] = '#C0392B'   # deep crimson — distinct from Perfect_A gold
        _f05_palette['Best_A']    = '#8E44AD'   # purple — distinct from Best_B blue
        # Draw largest tiers first (background) → smallest tiers last (foreground)
        tiers_by_size = sorted(existing_tiers,
                               key=lambda t: df[df['degrader_tier'] == t].shape[0],
                               reverse=True)
        for tier in tiers_by_size:
            sub = df[df['degrader_tier'] == tier]
            if sub.empty:
                continue
            is_pa  = tier == 'Perfect_A'
            is_top = tier in ('Perfect_B', 'Best_A')
            ax.scatter(sub['UMAP_X'], sub['UMAP_Y'],
                       c=_f05_palette.get(tier, '#999'),
                       alpha=0.95 if is_pa else (0.75 if is_top else 0.45),
                       s=85 if is_pa else (60 if is_top else 35),
                       linewidths=0.5 if (is_pa or is_top) else 0,
                       edgecolors='black' if (is_pa or is_top) else 'none',
                       label=f'{tier}  (n={len(sub):,})',
                       zorder=12 if is_pa else (5 if is_top else 2))
        gems = df[df['Conflict_Category'] == 'Hidden Gem'] if 'Conflict_Category' in df.columns else pd.DataFrame()
        if not gems.empty:
            ax.scatter(gems['UMAP_X'], gems['UMAP_Y'],
                       color=CONFLICT_PALETTE.get("Hidden Gem", '#CC79A7'),
                       marker='*', s=160, label='Hidden Gems', edgecolor='k', zorder=6)

        # KDE contour outlines per tier — top tiers only, overlaid on scatter
        try:
            from scipy.stats import gaussian_kde as _gkde05
            _top_tiers_05 = [t for t in existing_tiers[:4] if len(df[df['degrader_tier'] == t]) >= 15]
            _cont_alphas05 = [0.80, 0.68, 0.56, 0.45]
            for _ti05, _tier05 in enumerate(_top_tiers_05):
                _sub05 = df[df['degrader_tier'] == _tier05].dropna(subset=['UMAP_X', 'UMAP_Y'])
                if len(_sub05) < 15:
                    continue
                _col05 = _f05_palette.get(_tier05, '#999')
                _xy05 = np.vstack([_sub05['UMAP_X'], _sub05['UMAP_Y']])
                try:
                    _kde05 = _gkde05(_xy05)
                    _xx05 = np.linspace(float(df['UMAP_X'].min()), float(df['UMAP_X'].max()), 100)
                    _yy05 = np.linspace(float(df['UMAP_Y'].min()), float(df['UMAP_Y'].max()), 100)
                    _XX05, _YY05 = np.meshgrid(_xx05, _yy05)
                    _ZZ05 = _kde05(np.vstack([_XX05.ravel(), _YY05.ravel()])).reshape(_XX05.shape)
                    ax.contourf(_XX05, _YY05, _ZZ05, levels=2, colors=[_col05], alpha=0.07, zorder=1)
                    ax.contour(_XX05, _YY05, _ZZ05, levels=3, colors=[_col05],
                               alpha=_cont_alphas05[_ti05], linewidths=1.1, zorder=3)
                except Exception:
                    pass
        except ImportError:
            pass

        # Centroid markers + labels placed in the NEAREST EMPTY GRID CELL to avoid
        # overlapping data points.  A 20×20 occupancy grid flags cells that contain
        # UMAP points; label candidates shift outward until an empty cell is found.
        _ux5 = df['UMAP_X'].values
        _uy5 = df['UMAP_Y'].values
        _x_min5, _x_max5 = float(_ux5.min()), float(_ux5.max())
        _y_min5, _y_max5 = float(_uy5.min()), float(_uy5.max())
        _x_span5 = _x_max5 - _x_min5
        _y_span5 = _y_max5 - _y_min5
        _GRID = 24
        _occ5 = np.zeros((_GRID, _GRID), dtype=int)
        _xi_idx = np.clip((((_ux5 - _x_min5) / (_x_span5 + 1e-9)) * _GRID).astype(int), 0, _GRID - 1)
        _yi_idx = np.clip((((_uy5 - _y_min5) / (_y_span5 + 1e-9)) * _GRID).astype(int), 0, _GRID - 1)
        for _xi5, _yi5 in zip(_xi_idx, _yi_idx):
            _occ5[_yi5, _xi5] += 1

        def _find_empty_cell5(cx, cy, n_tries=30):
            """Return (lx, ly) near centroid (cx,cy) in an empty grid cell."""
            _ci = int(np.clip((cx - _x_min5) / (_x_span5 + 1e-9) * _GRID, 0, _GRID - 1))
            _ri = int(np.clip((cy - _y_min5) / (_y_span5 + 1e-9) * _GRID, 0, _GRID - 1))
            for _radius in range(1, n_tries):
                for _dc in range(-_radius, _radius + 1):
                    for _dr in range(-_radius, _radius + 1):
                        if abs(_dc) != _radius and abs(_dr) != _radius:
                            continue
                        _nc, _nr = _ci + _dc, _ri + _dr
                        if 0 <= _nc < _GRID and 0 <= _nr < _GRID and _occ5[_nr, _nc] == 0:
                            _lx = _x_min5 + (_nc + 0.5) / _GRID * _x_span5
                            _ly = _y_min5 + (_nr + 0.5) / _GRID * _y_span5
                            return _lx, _ly
            # Fallback: just above the centroid
            return cx, cy + _y_span5 * 0.08

        for tier in existing_tiers:
            sub_t = df[df['degrader_tier'] == tier]
            if len(sub_t) < 3:
                continue
            cx = float(sub_t['UMAP_X'].median())
            cy = float(sub_t['UMAP_Y'].median())
            col_t = _f05_palette.get(tier, '#999')
            ax.scatter(cx, cy, marker='X', s=85, color=col_t,
                       edgecolors='black', linewidths=1.0, zorder=8)
            _lbl_x, _lbl_y = _find_empty_cell5(cx, cy)
            ax.annotate(tier, xy=(cx, cy),
                        xytext=(_lbl_x, _lbl_y),
                        ha='center', va='center', fontsize=7.5,
                        color=col_t, fontweight='normal', zorder=9,
                        arrowprops=dict(arrowstyle='->', color=col_t,
                                        lw=0.9, shrinkA=0, shrinkB=4))

        ax.set_xlabel("UMAP Dimension 1  (distances reflect chemical similarity)", fontsize=10)
        ax.set_ylabel("UMAP Dimension 2", fontsize=10)
        ax.grid(False)
        # Rebuild handles in tier order (Perfect_A → worst), not size order
        _h5, _l5 = [], []
        from matplotlib.lines import Line2D as _L5
        for t in existing_tiers:
            sub_t = df[df['degrader_tier'] == t]
            if sub_t.empty:
                continue
            _h5.append(_L5([0], [0], marker='o', color='w',
                           markerfacecolor=_f05_palette.get(t, '#999'),
                           markeredgecolor='black', markeredgewidth=0.5,
                           markersize=8, label=f'{t}  (n={len(sub_t):,})'))
            _l5.append(f'{t}  (n={len(sub_t):,})')
        if not df[df['Conflict_Category'] == 'Hidden Gem'].empty if 'Conflict_Category' in df.columns else False:
            from matplotlib.lines import Line2D as _L5s
            _h5.append(_L5s([0], [0], marker='*', color='w',
                            markerfacecolor=CONFLICT_PALETTE.get("Hidden Gem", '#CC79A7'),
                            markeredgecolor='k', markersize=10, label='★ Hidden Gem (atypical)'))
        # Add centroid marker explanation
        _h5.append(_L5([0], [0], marker='X', color='w', markerfacecolor='#888',
                       markeredgecolor='black', markeredgewidth=0.8,
                       markersize=8, label='✕ Tier centroid (median UMAP)'))
        ax.legend(handles=_h5, loc='lower left',
                  ncol=2, fontsize=4.25, framealpha=0.92, fancybox=True)
        plt.tight_layout()
        plt.savefig(out_dir / "Figure_17a_Chemical_Space_Map.png", dpi=300, bbox_inches='tight')
        plt.close()
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_17a_Chemical_Space_Map.png').resolve()}")

    # Figure 17b: PA companion (chemical space landscape)
    if _pa_has_imgs:
        _fig_17b_pa_landscape(df, _pa, _imgs, out_dir, reporter)

    # --- Figure 18: Candidate Fingerprint: 18a (top-5 hits) + 18b (one per tier) ---
    _radar_labels = {
        "Boltz_Model_Confidence":   "AI Conf.",
        "iptm":                     "ipTM",
        "Interaction_Density_Norm": "Int.Den",
        "mean_plddt":               "pLDDT",
        "Binding_Probability":      "Bind.Prob",
        "SN2_Attack_Angle":         "SN2(°)",
        "Mechanistic_Fingerprint_Score": "Mech.FP",
        "Active_Site_RMSD":         "RMSD(Å)",
        "identity_pct":             "Seq.ID%",
        "SN2_Trajectory_Deviation_A": "Traj.Dev",
        "soft_catalytic_score":     "Soft.Cat",
    }

    def _draw_radar(rows_df, baseline_df, metrics, fig_path, subtitle,
                    label_mode='tier'):
        # label_mode='rank'  → "Fluoroacetate (Rank 1)"
        # label_mode='tier'  → "Fluoroacetate (Perfect_A)"
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
        fig_r, ax_r = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))
        colours_r = ["#057759", "#0BF1E2", "#E69F00", "#CC79A7", "#0072B2",
                     "#56B4E9", "#F0E442", "#009E73", "#D55E00", "#CC79A7"]
        lname_col = next((c for c in ['Ligand_Name', 'ligand'] if c in rows_df.columns), None)

        _ring_levels = [0.2, 0.4, 0.6, 0.8, 1.0]
        _ring_angles = np.linspace(0, 2 * np.pi, 300)
        for _rl in _ring_levels:
            ax_r.plot(_ring_angles, [_rl] * 300,
                      color='#CC2222', linewidth=0.55, alpha=0.45, zorder=1,
                      linestyle='-', solid_capstyle='round')

        import re as _re18
        for i in range(len(rows_df)):
            vals  = _scaler_vals.iloc[i].values.flatten().tolist() + [_scaler_vals.iloc[i].values[0]]
            _raw  = str(rows_df.iloc[i][lname_col]) if lname_col else f"Hit {i+1}"
            # Strip leading "26_" style number prefix
            _clean = _re18.sub(r'^\d+_', '', _raw)
            tier  = str(rows_df.iloc[i].get('degrader_tier', ''))
            if label_mode == 'rank':
                _rank_val = int(rows_df.iloc[i].get('Scientific_Rank', i + 1))
                lname_lbl = f"{_clean}  (Rank {_rank_val})"
            else:
                lname_lbl = f"{_clean}  ({tier})"
            c = colours_r[i % len(colours_r)]
            ax_r.plot(angles, vals, linewidth=2.2, linestyle='solid',
                      label=lname_lbl, color=c, zorder=3)
            ax_r.fill(angles, vals, color=c, alpha=0.08, zorder=2)
            ax_r.scatter(angles[:-1], vals[:-1],
                         color=c, s=55, zorder=5, edgecolors='white', linewidths=0.8)

        if not baseline_df.empty:
            bvals = _scaler_vals.iloc[-1].values.flatten().tolist() + [_scaler_vals.iloc[-1].values[0]]
            ax_r.plot(angles, bvals, linewidth=2, linestyle='--', color='#999',
                      label='Decoy avg.  (reference)', zorder=3)
            ax_r.fill(angles, bvals, color='#999', alpha=0.05, zorder=2)
            ax_r.scatter(angles[:-1], bvals[:-1],
                         color='#999', s=45, zorder=5, edgecolors='white', linewidths=0.8)

        _radar_fs = 8
        ax_r.set_xticks(angles[:-1])
        ax_r.set_xticklabels([])   # clear default labels; draw manually below
        ax_r.tick_params(axis='x', pad=18)
        # Draw spoke labels facing outward from the centre
        for _ri_s, (_ang_s, _lbl_s) in enumerate(zip(angles[:-1], spoke_labels)):
            _ang_deg_s = np.degrees(_ang_s)
            if 85 < _ang_deg_s < 95 or 265 < _ang_deg_s < 275:
                _ha_s = 'center'
            elif _ang_deg_s < 180:
                _ha_s = 'left'
            else:
                _ha_s = 'right'
            _rot_s = _ang_deg_s - 90
            if _ang_deg_s > 180:
                _rot_s += 180
            ax_r.text(_ang_s, ax_r.get_ylim()[1] * 1.15, _lbl_s,
                      ha=_ha_s, va='center', fontsize=_radar_fs,
                      rotation=_rot_s, rotation_mode='anchor')
        ax_r.set_yticks(_ring_levels)
        ax_r.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"],
                             color="#CC2222", size=7.5, fontweight='bold')
        ax_r.set_rlabel_position(45)
        ax_r.set_ylim(0, 1.10)
        ax_r.yaxis.grid(False)
        ax_r.xaxis.grid(True, color='#CCCCCC', linewidth=0.6, alpha=0.7)

        n_lines = len(rows_df) + (0 if baseline_df.empty else 1)
        ax_r.legend(loc='lower center', bbox_to_anchor=(0.5, -0.18),
                    ncol=min(3, n_lines), fontsize=_radar_fs,
                    framealpha=0.88, fancybox=True)
        plt.tight_layout(rect=[0, 0.15, 1, 1])
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close()

    try:
        metrics_r = [m for m in _radar_labels if m in df.columns]
        if metrics_r:
            worst_tier = existing_tiers[-1] if existing_tiers else None
            worst_avg  = (df[df['degrader_tier'] == worst_tier].mean(numeric_only=True).to_frame().T
                          if worst_tier else pd.DataFrame())

            # --- 18a: All Perfect_A complexes sorted by Scientific_Rank ---
            _pa18a = df[df['degrader_tier'] == 'Perfect_A'].copy()
            _sort_col18a = 'Scientific_Rank' if 'Scientific_Rank' in _pa18a.columns \
                           else 'Boltz_Model_Confidence'
            _asc18a = _sort_col18a == 'Scientific_Rank'
            _pa18a  = _pa18a.sort_values(_sort_col18a, ascending=_asc18a)
            if _pa18a.empty:
                _pa18a = df.sort_values('Boltz_Model_Confidence', ascending=False).head(8)
            _draw_radar(_pa18a, worst_avg, metrics_r,
                        out_dir / "Figure_18a_Fingerprint_TopHits.png",
                        "All Perfect_A complexes vs. worst-tier baseline",
                        label_mode='rank')
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_18a_Fingerprint_TopHits.png').resolve()}")

            # --- 18b: One representative per tier (best Boltz confidence per tier) ---
            reps_b = []
            for t in existing_tiers:
                sub_t = df[df['degrader_tier'] == t]
                if not sub_t.empty:
                    reps_b.append(sub_t.sort_values('Boltz_Model_Confidence', ascending=False).iloc[[0]])
            rep_df = pd.concat(reps_b) if reps_b else df.head(len(existing_tiers))
            _draw_radar(rep_df, pd.DataFrame(), metrics_r,
                        out_dir / "Figure_18b_Fingerprint_TierReps.png",
                        "One best representative per tier (Perfect_A → worst)",
                        label_mode='tier')
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_18b_Fingerprint_TierReps.png').resolve()}")
    except Exception as e:
        reporter.log(f"  ! Radar chart skipped: {e}")

    # --- Figure 19: Dual-Objective Tier Assessment — SN2 Geometry vs AI Confidence ---
    # Replaces the broken Pareto-rank approach (Binding_Probability was ~constant 0.99).
    # New design: for each tier, compute % meeting SN2 ≥ 165° (geometry) and
    # % meeting Boltz_Model_Confidence ≥ 0.75 (AI) — showing where each tier sits
    # in the biologically meaningful substrate-vs-inhibitor space.
    if "SN2_Attack_Angle" in df.columns and "Boltz_Model_Confidence" in df.columns:
        try:
            from matplotlib.lines import Line2D as _L19
            import matplotlib.patches as _mp19

            _d19 = df[['SN2_Attack_Angle', 'Boltz_Model_Confidence', 'degrader_tier']].copy()
            _d19['SN2_Attack_Angle']      = pd.to_numeric(_d19['SN2_Attack_Angle'],      errors='coerce')
            _d19['Boltz_Model_Confidence'] = pd.to_numeric(_d19['Boltz_Model_Confidence'], errors='coerce')
            _d19 = _d19.dropna()
            _tiers19 = [t for t in TIER_ORDER_LOGIC if t in _d19['degrader_tier'].unique()]

            # Per-tier metrics
            _rows19 = []
            for _t19 in _tiers19:
                _sub = _d19[_d19['degrader_tier'] == _t19]
                _n   = len(_sub)
                _pct_sn2  = float((_sub['SN2_Attack_Angle'] >= CFG.SUBSTRATE_ANGLE_MIN).mean() * 100)
                _pct_conf = float((_sub['Boltz_Model_Confidence'] >= CFG.SUBSTRATE_CONF_MIN).mean() * 100)
                _pct_both = float(((_sub['SN2_Attack_Angle'] >= CFG.SUBSTRATE_ANGLE_MIN) & (_sub['Boltz_Model_Confidence'] >= CFG.SUBSTRATE_CONF_MIN)).mean() * 100)
                _pct_inh  = float((_sub['SN2_Attack_Angle'] < CFG.INHIBITOR_ANGLE_MAX).mean() * 100)
                _med_sn2  = float(_sub['SN2_Attack_Angle'].median())
                _med_conf = float(_sub['Boltz_Model_Confidence'].median())
                _rows19.append({'tier': _t19, 'n': _n,
                                'pct_substrate': _pct_both,
                                'pct_sn2_ok': _pct_sn2,
                                'pct_conf_ok': _pct_conf,
                                'pct_inhibitor': _pct_inh,
                                'med_sn2': _med_sn2,
                                'med_conf': _med_conf})
            _df19 = pd.DataFrame(_rows19)

            # Figure 19 split into two standalone panels: 19a (tier % success bars)
            # and 19b (Conf × SN2 median scatter). Each saved as its own PNG.
            fig19a, ax19L = plt.subplots(figsize=(7.5, 7))
            fig19b, ax19R = plt.subplots(figsize=(7.5, 7))

            # ── Panel 19a: tier % success (substrate vs inhibitor) ───────────
            _y19 = np.arange(len(_df19))
            _bar_h = 0.25
            # Master container bar per tier — grey track, outline coloured to match
            # the tier's y-axis label colour. Shorter height leaves clear space
            # between adjacent tier groups; outline 70% of previous thickness.
            from matplotlib.colors import to_rgba as _to_rgba19
            _master_ec19 = [TIER_PALETTE.get(t, '#BDBDBD') for t in _df19['tier'].tolist()]
            # Faint tint of each tier's own colour as the track fill (alpha baked in
            # so the coloured outline stays vivid) — cohesive, not flat grey.
            _master_fc19 = [_to_rgba19(c, 0.13) for c in _master_ec19]
            ax19L.barh(_y19, 105, height=0.82, color=_master_fc19,
                       edgecolor=_master_ec19, linewidth=1.4, zorder=2)
            ax19L.barh(_y19 + _bar_h,  _df19['pct_substrate'].values,  height=_bar_h,
                       color='#009E73', alpha=0.85, label='Substrate (SN2≥165° + Conf≥0.75)', zorder=3)
            ax19L.barh(_y19,            _df19['pct_sn2_ok'].values,     height=_bar_h,
                       color='#56B4E9', alpha=0.75, label='SN2≥165° (any conf)', zorder=3)
            ax19L.barh(_y19 - _bar_h,  _df19['pct_inhibitor'].values,  height=_bar_h,
                       color='#D55E00', alpha=0.80, label='Potential inhibitor (SN2<145°)', zorder=3)

            # Value labels — inside (centred, white) when the bar is wide enough,
            # else just outside the bar's right end (dark text).
            for _i19, row19 in _df19.iterrows():
                _yi = int(_i19)
                for _pct_val, _yo in [(row19['pct_substrate'], _bar_h),
                                      (row19['pct_sn2_ok'], 0),
                                      (row19['pct_inhibitor'], -_bar_h)]:
                    if _pct_val < 1.0:
                        continue
                    if _pct_val >= 12:
                        ax19L.text(_pct_val / 2, _yi + _yo, f'{_pct_val:.0f}%',
                                   va='center', ha='center', fontsize=6,
                                   color='white', fontweight='bold', zorder=5)
                    else:
                        ax19L.text(_pct_val + 1.0, _yi + _yo, f'{_pct_val:.0f}%',
                                   va='center', ha='left', fontsize=6,
                                   color='#333', zorder=5)

            ax19L.set_yticks(_y19)
            ax19L.set_yticklabels(
                [f'{r["tier"]}  (n={r["n"]:,})' for _, r in _df19.iterrows()],
                fontsize=7
            )
            for tick, tier in zip(ax19L.get_yticklabels(), _df19['tier'].tolist()):
                tick.set_color(TIER_PALETTE.get(tier, 'black'))
                tick.set_fontweight('bold')
            ax19L.set_xlim(0, 105)
            ax19L.axvline(25, color='#999', lw=0.7, ls=':', zorder=2.4)
            ax19L.axvline(50, color='#999', lw=0.7, ls=':', zorder=2.4)
            ax19L.axvline(75, color='#999', lw=0.7, ls=':', zorder=2.4)
            ax19L.set_xlabel('Percentage of complexes (%)', fontsize=9)
            ax19L.tick_params(axis='x', labelsize=8)
            # Legend: single row, 70% font (7.5→5.25), inside bottom-right
            ax19L.legend(loc='lower right', ncol=3, fontsize=5.25,
                         framealpha=0.92, fancybox=True, borderpad=0.4,
                         handlelength=1.0, handletextpad=0.3, columnspacing=0.7)
            ax19L.invert_yaxis()

            # ── Right panel: scatter of tier medians in Conf × SN2 space ─────
            # Background hexbin of all complexes (shows density landscape)
            _hb19 = ax19R.hexbin(_d19['Boltz_Model_Confidence'], _d19['SN2_Attack_Angle'],
                                  gridsize=40, cmap='Greys', mincnt=1, alpha=0.4, linewidths=0.2, zorder=1)

            # Substrate / inhibitor zone shading — bounds from CFG
            ax19R.axhspan(CFG.SUBSTRATE_ANGLE_MIN, 185, color='#009E73', alpha=0.10, zorder=0)
            ax19R.axhspan(0, CFG.INHIBITOR_ANGLE_MAX, color='#D55E00', alpha=0.08, zorder=0)
            ax19R.axhline(CFG.SUBSTRATE_ANGLE_MIN, color='#009E73', lw=1.2, ls='--', alpha=0.7, zorder=2)
            ax19R.axhline(CFG.INHIBITOR_ANGLE_MAX, color='#D55E00', lw=1.2, ls='--', alpha=0.7, zorder=2)
            ax19R.axvline(CFG.SUBSTRATE_CONF_MIN, color='#0072B2', lw=1.2, ls='--', alpha=0.7, zorder=2)

            # Tier median diamonds with IQR error bars — collect positions first
            _tier_pts19 = []   # (tier, cx, cy, colour)
            for _t19 in _tiers19:
                _sub = _d19[_d19['degrader_tier'] == _t19]
                if _sub.empty:
                    continue
                _cx = float(_sub['Boltz_Model_Confidence'].median())
                _cy = float(_sub['SN2_Attack_Angle'].median())
                _ex_lo = _cx - float(_sub['Boltz_Model_Confidence'].quantile(0.25))
                _ex_hi = float(_sub['Boltz_Model_Confidence'].quantile(0.75)) - _cx
                _ey_lo = _cy - float(_sub['SN2_Attack_Angle'].quantile(0.25))
                _ey_hi = float(_sub['SN2_Attack_Angle'].quantile(0.75)) - _cy
                _col19 = TIER_PALETTE.get(_t19, '#999')
                ax19R.errorbar(_cx, _cy,
                               xerr=[[_ex_lo], [_ex_hi]],
                               yerr=[[_ey_lo], [_ey_hi]],
                               fmt='D', color=_col19, ms=11,
                               ecolor=_col19, elinewidth=1.5, capsize=3,
                               markeredgecolor='black', markeredgewidth=0.8,
                               zorder=6, label=f'{_t19}  (med SN2={_cy:.0f}°, conf={_cx:.3f})')
                _tier_pts19.append((_t19, _cx, _cy, _col19))

            ax19R.set_xlim(0.72, 1.01)
            ax19R.set_ylim(-5, 185)

            # Tier name labels parked in the empty left column, evenly spread on y so
            # none overlap; gentle curved arrows point from each label to its diamond.
            _pts_sorted19 = sorted(_tier_pts19, key=lambda p: p[2], reverse=True)
            _n_lab19 = len(_pts_sorted19)
            _lab_ys19 = list(np.linspace(160, 25, _n_lab19)) if _n_lab19 > 1 else [90]
            for (_t19, _cx, _cy, _col19), _laby in zip(_pts_sorted19, _lab_ys19):
                ax19R.annotate(_t19, xy=(_cx, _cy),
                               xytext=(0.83, _laby),
                               ha='right', va='center',
                               fontsize=6.5, color=_col19, fontweight='bold',
                               arrowprops=dict(arrowstyle='-|>', color=_col19, lw=1.3,
                                               mutation_scale=15, shrinkA=3, shrinkB=5,
                                               connectionstyle='arc3,rad=0.12'),
                               bbox=dict(boxstyle='round,pad=0.18', fc='white',
                                         ec=_col19, alpha=0.92, linewidth=0.8),
                               zorder=8)

            ax19R.set_xlabel('Boltz Model Confidence', fontsize=9)
            ax19R.set_ylabel('SN2 Attack Angle (°)', fontsize=9)
            ax19R.tick_params(axis='both', labelsize=8)

            # Zone / threshold labels — all inside the axes, single line each
            ax19R.text(0.87, 181, 'Substrate zone', color='#007A52', fontsize=7,
                       style='italic', ha='center', va='center')
            ax19R.text(0.87, 8, 'Inhibitor zone', color='#C04000', fontsize=7,
                       style='italic', ha='center', va='center')
            ax19R.text(0.745, 95, 'Conf threshold = 0.75', color='#0072B2', fontsize=6.5,
                       ha='center', va='center', style='italic', rotation=90)

            # Legend: 3 columns, shrunk icons, fixed outside above the top-left corner
            ax19R.legend(loc='lower left', bbox_to_anchor=(0.0, 1.02),
                         ncol=3, fontsize=5.25, framealpha=0.90, fancybox=True,
                         markerscale=0.45, handlelength=1.0, handletextpad=0.3,
                         columnspacing=0.8, labelspacing=0.3, borderpad=0.3)

            fig19a.savefig(out_dir / "Figure_19a_Tier_Success_Rates.png",
                           dpi=300, bbox_inches='tight')
            fig19b.savefig(out_dir / "Figure_19b_Conf_SN2_Landscape.png",
                           dpi=300, bbox_inches='tight')
            plt.close(fig19a)
            plt.close(fig19b)
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_19a_Tier_Success_Rates.png').resolve()}")
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_19b_Conf_SN2_Landscape.png').resolve()}")
        except Exception as e:
            reporter.log(f"  ! Figure 19 skipped: {e}")

    # --- Figure 20: Conflict Composition — 100% Stacked Horizontal Bar per Tier ---
    # One bar per tier (best → worst on y-axis).  Left portion = favourable categories
    # (Consensus High, Hidden Gem); right portion = concerning categories.
    # A thin vertical divider at the "favourable total" position separates the two sides.
    # Segment labels show % only for segments ≥ 5%.  Y-axis labels include tier n-count.
    if 'Conflict_Category' in df.columns and 'degrader_tier' in df.columns:
        cat_colours_f2 = {
            'Consensus High': '#009E73', 'Hidden Gem': '#CC79A7',
            'Decoy': '#D55E00', 'Consensus Low': '#777777', 'Ambiguous': '#E69F00'}
        # Canonical plotting order: favourable → concerning, left → right
        cat_order_f2 = ['Consensus High', 'Hidden Gem', 'Ambiguous', 'Consensus Low', 'Decoy']
        present_cats_f2 = [c for c in cat_order_f2 if c in df['Conflict_Category'].unique()]

        ct_f2 = pd.crosstab(df['degrader_tier'], df['Conflict_Category'])
        tier_rows_f2 = [t for t in existing_tiers if t in ct_f2.index]
        ct_f2 = ct_f2.reindex(index=tier_rows_f2, columns=present_cats_f2, fill_value=0)
        tier_n_f2   = ct_f2.sum(axis=1)
        # Grand total = all rows with a valid degrader_tier (not filtered by Conflict_Category)
        _total_rows_f2 = len(df)
        grand_total    = int(df['degrader_tier'].notna().sum())
        pct_f2         = ct_f2.div(tier_n_f2, axis=0).fillna(0) * 100
        # Detect control cases — generic mask first, then specific 3R3U / DeHa4
        _ctrl_mask = pd.Series(False, index=df.index)
        for _ctrl_col in ['is_control', 'control', 'Control', 'is_Control']:
            if _ctrl_col in df.columns:
                _ctrl_mask = df[_ctrl_col].astype(bool)
                break
        if not _ctrl_mask.any() and 'job_name' in df.columns:
            _ctrl_mask = df['job_name'].str.lower().str.contains('control|ctrl|deha4_ref|3r3u', na=False)
        _n_controls = int(_ctrl_mask.sum())
        # Identify the two reference controls by their exact job-name patterns:
        #   3R3U  control: job_name contains '3r3u'  AND 'control'
        #   DeHa4 control: job_name contains 'deha4' AND 'control'
        # E.g.: "0000000_4_3R3U_Control_26_Fluoroacetate"
        #       "0000000_02_DeHa4_Control_26_Fluoroacetate"
        _ctrl_3R3U_mask  = pd.Series(False, index=df.index)
        _ctrl_deha4_mask = pd.Series(False, index=df.index)
        if 'job_name' in df.columns:
            _jn_lower = df['job_name'].str.lower()
            _ctrl_3R3U_mask = (
                _jn_lower.str.contains('3r3u', na=False) &
                _jn_lower.str.contains('control|ctrl', na=False)
            )
            _ctrl_deha4_mask = (
                _jn_lower.str.contains('deha4', na=False) &
                _jn_lower.str.contains('control|ctrl', na=False)
            )
            # Fallback: if 'control' keyword absent, match on 3r3u / deha4 alone
            if not _ctrl_3R3U_mask.any():
                _ctrl_3R3U_mask = _jn_lower.str.contains('3r3u', na=False)
            if not _ctrl_deha4_mask.any():
                _ctrl_deha4_mask = _jn_lower.str.contains('deha4', na=False)
        _n_3R3U  = int(_ctrl_3R3U_mask.sum())
        _n_deha4 = int(_ctrl_deha4_mask.sum())

        # Y positions (top tier at top)
        tier_y_f2   = np.arange(len(tier_rows_f2) - 1, -1, -1)
        bar_h       = 0.60

        # Height: fit bars tightly — no whitespace padding below lowest tier
        _fh_f2 = max(4.0, len(tier_rows_f2) * 0.65 + 1.8)
        fig_f2, ax_f2 = plt.subplots(figsize=(13, _fh_f2))
        ax_f2.set_facecolor('#FAFAFA')

        left_f2 = np.zeros(len(tier_rows_f2))
        for cat in present_cats_f2:
            vals   = pct_f2[cat].values
            counts = ct_f2[cat].values
            ax_f2.barh(tier_y_f2, vals, left=left_f2, height=bar_h,
                       color=cat_colours_f2.get(cat, '#999'),
                       edgecolor='white', linewidth=0.7, zorder=2)
            for yi, vi, li, ni in zip(tier_y_f2, vals, left_f2, counts):
                if vi < 4:
                    continue
                txt_col = 'white' if vi > 12 else 'black'
                if vi >= 9:
                    lbl = f'{vi:.0f}%\n(n={int(ni):,})'
                    fs  = 7.0
                else:
                    lbl = f'{vi:.0f}%'
                    fs  = 7.5
                ax_f2.text(li + vi / 2, yi, lbl,
                           ha='center', va='center', fontsize=fs,
                           color=txt_col, fontweight='bold', zorder=4,
                           linespacing=1.1)
            left_f2 = left_f2 + vals

        # Vertical divider at boundary between favourable and concerning categories
        fav_cats_f2 = [c for c in ['Consensus High', 'Hidden Gem'] if c in present_cats_f2]
        fav_pct_f2  = pct_f2[fav_cats_f2].sum(axis=1).values if fav_cats_f2 else None
        if fav_pct_f2 is not None:
            for yi, xv in zip(tier_y_f2, fav_pct_f2):
                ax_f2.plot([xv, xv], [yi - bar_h / 2, yi + bar_h / 2],
                           color='white', linewidth=2.2, zorder=5)

        # Reference line at 50%
        ax_f2.axvline(50, color='#555555', linewidth=1.1, linestyle='--', alpha=0.5, zorder=1)
        ax_f2.text(50.5, len(tier_rows_f2) - 0.5, '50%',
                   va='top', ha='left', fontsize=8, color='#555')

        # Y-axis labels: tier name + n-count in matching colour
        ax_f2.set_yticks(tier_y_f2)
        ax_f2.set_yticklabels(
            [f'{t}  (n={tier_n_f2.loc[t]:,})' for t in tier_rows_f2],
            fontsize=10)
        for tick, tier in zip(ax_f2.get_yticklabels(), tier_rows_f2):
            tick.set_color(TIER_PALETTE.get(tier, 'black'))
            tick.set_fontweight('bold')

        ax_f2.set_xlim(0, 100)
        # Clamp y-range to exactly the number of tiers — prevents blank canvas below
        ax_f2.set_ylim(-0.55, len(tier_rows_f2) - 0.45)
        ax_f2.set_xlabel("Percentage of complexes within tier  (%)", fontsize=11)
        ax_f2.set_xticks([0, 25, 50, 75, 100])
        ax_f2.set_xticklabels(['0%', '25%', '50%', '75%', '100%'])
        ax_f2.set_axisbelow(True)
        ax_f2.xaxis.grid(True, alpha=0.30, linestyle=':', color='#888', zorder=0)
        #ax_f2.set_title("Conflict-Category Composition by Degrader Tier", fontsize=13, pad=10)

        # Top-right badge: total jobs | tiered count — lifted above the axes frame
        ax_f2.text(0.99, 1.06,
                   f'Total jobs: {_total_rows_f2:,}   |   Tiered: {grand_total:,}',
                   transform=ax_f2.transAxes, ha='right', va='bottom',
                   fontsize=8.5, color='#333', fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.18', fc='white', ec='#ccc',
                             alpha=0.88, linewidth=0.5))

        # Section labels — placed INSIDE the chart at the top of the bars so they do not
        # collide with the "Total jobs | Tiered" badge that floats above the axes frame.
        # get_xaxis_transform(): x = data coords (0–100 %), y = axes fraction [0,1].
        if fav_cats_f2:
            mid_fav = np.mean(fav_pct_f2) / 2
            mid_con = np.mean(fav_pct_f2) + (100 - np.mean(fav_pct_f2)) / 2
            ax_f2.text(mid_fav, 1.015, '▶ Favourable',
                       ha='center', va='bottom', fontsize=9.5,
                       color='#009E73', fontweight='bold',
                       transform=ax_f2.get_xaxis_transform(),
                       bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='#009E73',
                                 alpha=0.80, linewidth=0.5))
            ax_f2.text(mid_con, 1.015, 'Concerning ◀',
                       ha='center', va='bottom', fontsize=9.5,
                       color='#D55E00', fontweight='bold',
                       transform=ax_f2.get_xaxis_transform(),
                       bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='#D55E00',
                                 alpha=0.80, linewidth=0.5))

        # Arrow annotations for small segments (< 4%) that couldn't fit inline text
        left_f2_tmp = np.zeros(len(tier_rows_f2))
        for cat in present_cats_f2:
            vals   = pct_f2[cat].values
            counts = ct_f2[cat].values
            for yi, vi, li, ni in zip(tier_y_f2, vals, left_f2_tmp, counts):
                if vi > 0 and vi < 4 and ni > 0:
                    mid_x = li + vi / 2
                    col_arr = cat_colours_f2.get(cat, '#999')
                    ax_f2.annotate(f'n={int(ni):,}',
                                   xy=(mid_x, yi),
                                   xytext=(mid_x, yi + 0.38),
                                   ha='center', va='bottom', fontsize=6.5,
                                   color=col_arr, fontweight='bold',
                                   arrowprops=dict(arrowstyle='->', color=col_arr,
                                                   lw=0.7, shrinkB=2),
                                   bbox=dict(boxstyle='round,pad=0.10', fc='white',
                                             ec=col_arr, alpha=0.85, linewidth=0.5),
                                   zorder=6)
            left_f2_tmp = left_f2_tmp + vals

        # Highlight each control complex individually — labelled protein+ligand
        # (e.g. "DeHa4+TFA"). 6 controls: {3R3U, DeHa4} × {TFA, Fluoroacetate,
        # Difluoroacetate}. Stacked vertically per tier row so they never overlap.
        def _lig_short_f20(_idx):
            _ln = ''
            if 'Ligand_Name' in df.columns:
                _ln = str(df.at[_idx, 'Ligand_Name'])
            if (not _ln or _ln.lower() == 'nan') and 'job_name' in df.columns:
                _jn = str(df.at[_idx, 'job_name'])
                _ln = _jn.split('Control_', 1)[1] if 'Control_' in _jn else _jn
            # strip a leading numeric index like "25_TFA" → "TFA"
            if '_' in _ln and _ln.split('_', 1)[0].isdigit():
                _ln = _ln.split('_', 1)[1]
            return _ln.replace('_', ' ').strip()

        _ctrl_specs_f20 = [(_ctrl_3R3U_mask, '3R3U', '#8B0000'),
                           (_ctrl_deha4_mask, 'DeHa4', '#00408B')]
        _ctrl_by_tier_f20 = {}   # tier -> [(label, colour), ...]
        for _cmask, _cprot, _ccol in _ctrl_specs_f20:
            if not (_cmask.any() and 'job_name' in df.columns):
                continue
            for _cidx in df.index[_cmask]:
                _ctier = df.at[_cidx, 'degrader_tier']
                if _ctier not in tier_rows_f2:
                    continue
                _ctrl_by_tier_f20.setdefault(_ctier, []).append(
                    (f'{_cprot}+{_lig_short_f20(_cidx)}', _ccol))
        for _ctier, _items in _ctrl_by_tier_f20.items():
            _cy = tier_y_f2[tier_rows_f2.index(_ctier)]
            _k  = len(_items)
            _offs = np.linspace(0.34, -0.34, _k) if _k > 1 else [0.0]
            for (_clbl, _ccol), _off in zip(_items, _offs):
                ax_f2.annotate(f'◀ {_clbl}', xy=(100, _cy),
                               xytext=(99, _cy + _off),
                               ha='right', va='center', fontsize=6.8,
                               color=_ccol, fontweight='bold',
                               arrowprops=dict(arrowstyle='->', color=_ccol, lw=0.9,
                                               shrinkA=1, shrinkB=2),
                               bbox=dict(boxstyle='round,pad=0.12', fc='white',
                                         ec=_ccol, alpha=0.90, linewidth=0.5),
                               zorder=7)

        from matplotlib.patches import Patch as _P02
        _leg02 = [_P02(facecolor=cat_colours_f2.get(c, '#999'), edgecolor='white',
                       linewidth=0.5, label=c) for c in present_cats_f2]
        # Single-row legend below chart — ncol = number of categories so all fit on one line
        _leg02_obj = ax_f2.legend(handles=_leg02,
                                  loc='upper center', bbox_to_anchor=(0.5, -0.10),
                                  ncol=len(present_cats_f2),
                                  fontsize=8.5, framealpha=0.92, fancybox=True,
                                  title='Conflict Category', title_fontsize=8)
        _leg02_obj.set_zorder(20)

        fig_f2.tight_layout(rect=[0, 0.07, 1, 1])
        plt.savefig(out_dir / "Figure_20_Conflict_Composition.png", dpi=300, bbox_inches='tight')
        plt.close()
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_20_Conflict_Composition.png').resolve()}")

    # --- Figure 21: Hidden Gems — Slope Chart (SN2 Angle vs AI Confidence) ---
    # Each Hidden Gem complex is a lollipop row showing its absolute SN2 attack angle
    # (physics metric, left dot) and Boltz Model Confidence (AI metric, right dot).
    # A large gap means strong mechanistic geometry but undervalued by AI — the
    # "hidden gem" story in a single visual stroke.
    _phys_col21_abs = next((c for c in ['SN2_Attack_Angle', 'Scientific_Rank', 'Pareto_Rank'] if c in df.columns), None)
    _ai_col21_abs   = next((c for c in ['Boltz_Model_Confidence', 'Ensemble_Data_Rank'] if c in df.columns), None)
    if 'Conflict_Category' in df.columns and _phys_col21_abs and _ai_col21_abs:
        try:
            gems = df[df['Conflict_Category'] == 'Hidden Gem'].copy()
            if not gems.empty:
                # Use absolute values for physics and AI metrics
                gems['_phys_pct21'] = pd.to_numeric(gems[_phys_col21_abs], errors='coerce')
                gems['_ai_pct21']   = pd.to_numeric(gems[_ai_col21_abs],   errors='coerce')
                gems = gems.dropna(subset=['_phys_pct21', '_ai_pct21'])
                # Normalise to 0–100 scale for consistent axis display
                _phys_min21, _phys_max21 = float(gems['_phys_pct21'].min()), float(gems['_phys_pct21'].max())
                _ai_min21,   _ai_max21   = float(gems['_ai_pct21'].min()),   float(gems['_ai_pct21'].max())
                _phys_rng21 = _phys_max21 - _phys_min21 if _phys_max21 > _phys_min21 else 1.0
                _ai_rng21   = _ai_max21   - _ai_min21   if _ai_max21   > _ai_min21   else 1.0
                gems['_phys_pct21'] = (gems['_phys_pct21'] - _phys_min21) / _phys_rng21 * 100
                gems['_ai_pct21']   = (gems['_ai_pct21']   - _ai_min21)   / _ai_rng21   * 100
                gems['_gap21']      = gems['_phys_pct21'] - gems['_ai_pct21']

                lname_col21 = next((c for c in ['complex_id', 'Ligand_Name', 'ligand',
                                                'Protein_Name', 'protein']
                                    if c in gems.columns), None)

                # ── Horizontal dual-dot lollipop chart ──────────────────────────
                # Each Hidden Gem is a row.  Two coloured dots show its Physics rank
                # (blue circle) and AI rank (orange diamond) on a shared x-axis (0–100%).
                # A grey connecting line spans between the two dots; width ∝ gap.
                # Sorted by AI rank (best at top) so the "most overlooked" gems sit
                # at the bottom where they are easy to read.
                from matplotlib.lines import Line2D as _L21
                gems_s21 = gems.sort_values('_ai_pct21', ascending=True)  # best at top after invert
                _ng21 = len(gems_s21)
                _y_pos21 = np.arange(_ng21)

                # Dynamic x-axis range so clustered low-percentile gems fill the axis
                _all_pct21 = (list(gems_s21['_phys_pct21'].values) +
                              list(gems_s21['_ai_pct21'].values))
                _x_data_max21 = max(_all_pct21) if _all_pct21 else 5.0
                _x_hi21 = float(np.clip(_x_data_max21 * 4, 20, 100))

                fig21, ax21 = plt.subplots(figsize=(12, max(4.5, _ng21 * 0.70 + 2.2)))
                ax21.set_facecolor('#FAFAFA')
                ax21.set_axisbelow(True)
                ax21.xaxis.grid(True, color='#DCDCDC', linewidth=0.7, zorder=0)

                # Background zone bands — tertile reference lines
                ax21.axvspan(0,  33, alpha=0.04, color='#D55E00', zorder=0)  # lower third
                ax21.axvspan(33, 67, alpha=0.03, color='#E69F00', zorder=0)  # middle
                ax21.axvspan(67, 100, alpha=0.04, color='#009E73', zorder=0)  # top third
                for _xv21 in [33, 67]:
                    ax21.axvline(_xv21, color='#CCCCCC', linewidth=0.8,
                                 linestyle='--', zorder=1)

                for _gi21, (_, row21) in enumerate(gems_s21.iterrows()):
                    _phys21 = float(row21['_phys_pct21'])
                    _ai21   = float(row21['_ai_pct21'])
                    _gap21  = float(row21['_gap21'])
                    _tier21 = (row21.get('degrader_tier', '')
                               if 'degrader_tier' in gems.columns else '')
                    _col21  = TIER_PALETTE.get(_tier21, '#AA3377')
                    _y21    = _gi21

                    # Connecting line — width proportional to |gap|
                    _lw21 = np.clip(0.8 + abs(_gap21) / 25, 0.8, 3.5)
                    _lo_x21, _hi_x21 = sorted([_phys21, _ai21])
                    ax21.hlines(_y21, _lo_x21, _hi_x21,
                                colors='#CCCCCC', linewidth=_lw21,
                                alpha=0.75, zorder=2)

                    # Physics dot (blue circle, red outline)
                    ax21.scatter(_phys21, _y21, color='#0072B2', s=130, zorder=4,
                                 marker='o', edgecolors='red', linewidths=1.2)
                    ax21.text(_phys21, _y21 + 0.34, f'{_phys21:.0f}%',
                              va='bottom', ha='center', fontsize=5.5, color='#0072B2',
                              fontweight='bold', zorder=5)
                    # AI dot (golden diamond, no border)
                    ax21.scatter(_ai21, _y21, color='#E69F00', s=130, zorder=4,
                                 marker='D', edgecolors='none')
                    ax21.text(_ai21, _y21 - 0.34, f'{_ai21:.0f}%',
                              va='top', ha='center', fontsize=5.5, color='#E69F00',
                              fontweight='bold', zorder=5)

                    # Value labels just outside the right dot
                    _rightmost21 = max(_phys21, _ai21)
                    ax21.text(_rightmost21 + 1.2, _y21,
                              f'Δ={abs(_gap21):.1f}%',
                              va='center', ha='left', fontsize=7.5,
                              color='#555', fontweight='bold')

                # Y-axis: gem names (reversed so top row = best AI rank)
                _gem_labels21 = []
                for _gi_lbl21, (_, row21) in enumerate(gems_s21.iterrows()):
                    if lname_col21:
                        _gem_labels21.append(str(row21[lname_col21])[:28])
                    else:
                        _gem_labels21.append(f'Hidden Gem {_gi_lbl21 + 1}')
                ax21.set_yticks(_y_pos21)
                ax21.set_yticklabels(_gem_labels21, fontsize=8.5)
                for _tick21, (_, row21) in zip(ax21.get_yticklabels(), gems_s21.iterrows()):
                    _t21 = (row21.get('degrader_tier', '')
                            if 'degrader_tier' in gems.columns else '')
                    _tick21.set_color(TIER_PALETTE.get(_t21, '#333'))
                    _tick21.set_fontweight('bold')

                ax21.set_xlim(-1, _x_hi21 + 2)
                _x21_phys_lbl = ('SN2 Attack Angle' if _phys_col21_abs == 'SN2_Attack_Angle'
                                 else 'Physics Rank')
                _x21_ai_lbl   = ('Boltz Confidence' if _ai_col21_abs == 'Boltz_Model_Confidence'
                                 else 'AI Rank')
                ax21.set_xlabel(
                    f'Normalised Score  ({_x21_phys_lbl} / {_x21_ai_lbl})  '
                    '(0 = lowest · 100 = highest)',
                    fontsize=10)
                _all_tick_vals21 = [0, 5, 10, 15, 20, 25, 33, 50, 67, 80, 90, 100]
                _tick_lbl_map21  = {0: '0%', 5: '5%', 10: '10%', 15: '15%',
                                    20: '20%', 25: '25%', 33: '33rd', 50: '50%',
                                    67: '67th', 80: '80%', 90: '90%', 100: '100%'}
                _vis_ticks21 = [t for t in _all_tick_vals21 if t <= _x_hi21 + 1]
                ax21.set_xticks(_vis_ticks21)
                ax21.set_xticklabels([_tick_lbl_map21[t] for t in _vis_ticks21], fontsize=8.5)
                ax21.set_ylim(-0.6, _ng21 - 0.4)

                # Zone labels at top — only render those within the visible x-range
                _zone_lbls21 = [
                    (16.5,  'Lower third', '#D55E00', 33),
                    (50.0,  'Middle third', '#8A6000', 67),
                    (83.5,  'Top third',   '#009E73', 100),
                ]
                for _zx21, _ztxt21, _zcol21, _zmax21 in _zone_lbls21:
                    if _zx21 <= _x_hi21 + 1:
                        ax21.text(_zx21, 0.97, _ztxt21, ha='center', va='top',
                                  fontsize=7.5, color=_zcol21, style='italic',
                                  transform=ax21.get_xaxis_transform())

                # Legend
                _lh21 = [
                    _L21([0],[0], marker='o', color='w', markerfacecolor='#0072B2',
                         markeredgecolor='red', markeredgewidth=1.2,
                         markersize=9, label=f'Physics score ({_x21_phys_lbl})'),
                    _L21([0],[0], marker='D', color='w', markerfacecolor='#E69F00',
                         markeredgecolor='none',
                         markersize=9, label=f'AI score ({_x21_ai_lbl})'),
                ]
                # Add tier entries if multiple tiers
                _tier_in_gems = [t for t in existing_tiers
                                 if 'degrader_tier' in gems.columns and
                                 t in gems['degrader_tier'].values]
                _lh21 += [_L21([0],[0], color=TIER_PALETTE.get(t,'#999'), linewidth=3,
                                label=f'{t}  (n={len(gems[gems["degrader_tier"]==t]):,})')
                          for t in _tier_in_gems]
                ax21.legend(handles=_lh21, loc='lower right', fontsize=8,
                            framealpha=0.92, fancybox=True, ncol=2,
                            title='Marker type  |  Tier (line colour)',
                            title_fontsize=7.5)

                # Badge: summary stats
                n_gems21 = len(gems)
                med_gap21 = float(gems['_gap21'].median())
                ax21.text(0.99, 0.99,
                          f'Hidden Gems: {n_gems21:,}\nMedian |Δ|: {abs(med_gap21):.1f}%',
                          transform=ax21.transAxes, va='top', ha='right',
                          fontsize=8.5, color='#333',
                          bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='#ccc',
                                    alpha=0.88, linewidth=0.5))

                plt.tight_layout()
                df.drop(columns=['_phys_pct21', '_ai_pct21'], errors='ignore', inplace=True)
                plt.savefig(out_dir / "Figure_21_Hidden_Gems_DeepDive.png",
                            dpi=300, bbox_inches='tight')
                plt.close(fig21)
                reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_21_Hidden_Gems_DeepDive.png').resolve()}")
            else:
                reporter.log("  ! Figure 21 skipped: No Hidden Gems found in dataset.")
        except Exception as e:
            reporter.log(f"  ! Figure 21 skipped: {e}")

    # --- Figure 22: Category Overlap — Euler / Venn Diagram ---
    try:
        _phys_col22 = next((c for c in ['Pareto_Rank', 'Ensemble_Data_Rank'] if c in df.columns), None)
        _ai_col22   = next((c for c in ['Boltz_Model_Confidence', 'iptm', 'ptm'] if c in df.columns), None)
        _tier_strong22 = {'Perfect_A', 'Perfect_B', 'Best_A'}

        if _phys_col22 and _ai_col22 and 'degrader_tier' in df.columns:
            _n22 = len(df)
            _phys_thresh22 = float(df[_phys_col22].quantile(0.33))
            _set_A = set(df.index[df[_phys_col22] <= _phys_thresh22])
            _ai_thresh22   = float(df[_ai_col22].quantile(0.67))
            _set_B = set(df.index[pd.to_numeric(df[_ai_col22], errors='coerce') >= _ai_thresh22])
            _set_C = set(df.index[df['degrader_tier'].isin(_tier_strong22)])

            _ctrl_3R3U_idx22  = None
            _ctrl_deha4_idx22 = None
            if 'job_name' in df.columns:
                _jn22 = df['job_name'].str.lower()
                _3r3u_m22  = (_jn22.str.contains('3r3u',  na=False) &
                              _jn22.str.contains('control|ctrl', na=False))
                _deha4_m22 = (_jn22.str.contains('deha4', na=False) &
                              _jn22.str.contains('control|ctrl', na=False))
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

            fig22, ax22 = plt.subplots(figsize=(11, 9))
            ax22.set_aspect('equal')
            ax22.set_xlim(0.5, 9.5)
            ax22.set_ylim(0.4, 8.7)
            ax22.axis('off')

            _r22 = 2.65
            _cx_A, _cy_A = 3.5, 5.8
            _cx_B, _cy_B = 6.5, 5.8
            _cx_C, _cy_C = 5.0, 3.4
            for _cx22, _cy22, _col22 in [
                    (_cx_A, _cy_A, '#0072B2'),
                    (_cx_B, _cy_B, '#D55E00'),
                    (_cx_C, _cy_C, '#009E73')]:
                ax22.add_patch(_Circ22((_cx22, _cy22), _r22,
                                       facecolor=_col22, alpha=0.20,
                                       edgecolor=_col22, linewidth=2.2, zorder=2))

            _abc_cx, _abc_cy, _r_abc = 5.0, 5.10, 0.75
            ax22.add_patch(_Circ22((_abc_cx, _abc_cy), _r_abc,
                                   facecolor='#AAAAAA', alpha=0.15,
                                   edgecolor='#444', linewidth=2.0,
                                   linestyle=':', zorder=4))

            _set_PA22 = set(df.index[df['degrader_tier'] == 'Perfect_A'])
            _n_PA22   = len(_set_PA22)
            _pa_cx22, _pa_cy22, _r_PA22, _pa_col22 = _cx_C, _cy_C, 0.50, '#CC79A7'
            if _n_PA22 > 0:
                _fA = len(_set_PA22 & _set_A) / _n_PA22
                _fB = len(_set_PA22 & _set_B) / _n_PA22
                _fC = len(_set_PA22 & _set_C) / _n_PA22
                _ww = _fA + _fB + _fC
                if _ww > 0:
                    _pa_cx22 = (_fA * _cx_A + _fB * _cx_B + _fC * _cx_C) / _ww
                    _pa_cy22 = (_fA * _cy_A + _fB * _cy_B + _fC * _cy_C) / _ww
                ax22.add_patch(_Circ22((_pa_cx22, _pa_cy22), _r_PA22,
                                       facecolor=_pa_col22, alpha=0.28,
                                       edgecolor=_pa_col22, linewidth=2.5,
                                       linestyle='--', zorder=5))
                ax22.annotate(f'Perfect_A\n(n={_n_PA22:,})',
                              xy=(_pa_cx22 - _r_PA22 * 0.7, _pa_cy22 + _r_PA22 * 0.7),
                              xytext=(_pa_cx22 - 2.2, _pa_cy22 + 1.4),
                              ha='center', va='bottom', fontsize=8,
                              color=_pa_col22, fontweight='bold', zorder=9,
                              bbox=dict(boxstyle='round,pad=0.20', fc='#F9E4F2',
                                        ec=_pa_col22, alpha=0.93, linewidth=1.3),
                              arrowprops=dict(arrowstyle='->', color=_pa_col22,
                                              lw=1.1, shrinkA=0, shrinkB=3))

            # Primary set labels — uniform size, bold, distinct filled-and-outlined box
            _main_lbl_kw22 = dict(ha='center', va='center', fontsize=11,
                                  fontweight='bold', linespacing=1.25, zorder=6)
            ax22.text(2.15, 7.25, 'Physics-Strong', color='#0072B2', **_main_lbl_kw22,
                      bbox=dict(boxstyle='round,pad=0.40', fc='#D6EAF8', ec='#0072B2',
                                alpha=0.95, linewidth=1.8))
            ax22.text(7.45, 7.25, 'AI-Strong', color='#D55E00', **_main_lbl_kw22,
                      bbox=dict(boxstyle='round,pad=0.40', fc='#FCE4D0', ec='#D55E00',
                                alpha=0.95, linewidth=1.8))
            ax22.text(5.00, 1.55, 'Tier-Strong', color='#009E73', **_main_lbl_kw22,
                      bbox=dict(boxstyle='round,pad=0.40', fc='#D4EFDF', ec='#009E73',
                                alpha=0.95, linewidth=1.8))

            _region_data22 = [
                (1.85, 5.80, _only_A,    '#0072B2', 'Physics only'),
                (8.15, 5.80, _only_B,    '#D55E00', 'AI only'),
                (5.00, 2.20, _only_C,    '#009E73', 'Tier only'),
                (5.00, 6.90, _AB_not_C,  '#1A4080', 'Physics ∩ AI'),
                (2.80, 4.40, _AC_not_B,  '#005840', 'Physics ∩ Tier'),
                (7.20, 4.40, _BC_not_A,  '#803000', 'AI ∩ Tier'),
                (5.00, 5.30, _ABC,       '#222',    'All three'),
            ]
            # Second-layer region labels — lighter square boxes, name normal / count bold,
            # visually subordinate to the primary set labels above.
            for _rx22, _ry22, _rn22, _rc22, _rl22 in _region_data22:
                _pct22 = _rn22 / _n22 * 100 if _n22 else 0
                ax22.text(_rx22, _ry22,
                          f'{_rl22}\n' + r'$\bf{' + f'{_rn22:,}' + r'}$' + f'  ({_pct22:.1f}%)',
                          ha='center', va='center', fontsize=7,
                          color=_rc22, linespacing=1.35,
                          bbox=dict(boxstyle='square,pad=0.32', fc='#FBFBFB', ec=_rc22,
                                    alpha=0.92, linewidth=0.9), zorder=5)

            _ctrl_items22 = [
                (_ctrl_3R3U_idx22,  '3R3U (control)',  '#8B0000'),
                (_ctrl_deha4_idx22, 'DeHa4 (control)', '#00408B'),
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
                ax22.scatter([_mx22], [_my22], marker='*', s=280,
                             color=_ccol22, edgecolors='black',
                             linewidths=0.6, zorder=16, clip_on=False)
                _ctrl_legend_h22.append(
                    _L22([0], [0], marker='*', color='w',
                         markerfacecolor=_ccol22, markeredgecolor='black',
                         markersize=11, markeredgewidth=0.6,
                         label=f'★  {_clbl22}'))

            _leg22 = [
                _Patch22(facecolor='#0072B2', alpha=0.45, edgecolor='#0072B2',
                         linewidth=1.2, label='Physics-Strong  (top-33% Pareto rank)'),
                _Patch22(facecolor='#D55E00', alpha=0.45, edgecolor='#D55E00',
                         linewidth=1.2, label='AI-Strong  (top-33% confidence)'),
                _Patch22(facecolor='#009E73', alpha=0.45, edgecolor='#009E73',
                         linewidth=1.2, label='Tier-Strong  (Perfect_A/B, Best_A)'),
            ]
            if _n_PA22 > 0:
                _leg22.append(_L22([0], [0], color='#CC79A7', linewidth=2.2,
                                   linestyle='--',
                                   label=f'Perfect_A overlay  (n={_n_PA22:,})'))
            _leg22.append(_L22([0], [0], color='#444', linewidth=1.8,
                               linestyle=':', marker='o', markersize=12,
                               markerfacecolor='none', markeredgecolor='#444',
                               markeredgewidth=1.5,
                               label=f'All-three intersection  ({_ABC:,})'))
            _leg22.extend(_ctrl_legend_h22)
            # Legend pulled close under the diagram (footnote clutter removed)
            ax22.legend(handles=_leg22,
                        loc='upper center', bbox_to_anchor=(0.5, 0.02),
                        ncol=3, fontsize=8, framealpha=0.90, fancybox=True,
                        borderpad=0.5, columnspacing=1.2, handletextpad=0.5)

            plt.tight_layout()
            plt.savefig(out_dir / "Figure_22_Category_Overlap_Euler.png",
                        dpi=300, bbox_inches='tight')
            plt.close()
            reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_22_Category_Overlap_Euler.png').resolve()}")
        else:
            reporter.log("  ! Figure 22 skipped: requires Pareto_Rank, confidence column, and degrader_tier")
    except Exception as e:
        reporter.log(f"  ! Figure 22 skipped: {e}")

    _fig23_multitarget(df, out_dir, reporter)
    _fig24_sankey(df, out_dir, reporter)
    _fig25_pfas_size(df, out_dir, reporter)


# ===============================================================================
# SECTION 4C: Publication Assembly Figures (23 / 24)
# ===============================================================================

def _load_pfas_smiles() -> dict:
    """Load PFAS SMILES from D_PFAS27_Tue.smi in the script directory."""
    path = Path(__file__).resolve().parent / "D_PFAS27_Tue.smi"
    if not path.exists():
        return {}
    smiles: dict = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith('#'):
            continue
        parts = line.strip().split()
        if len(parts) >= 2:
            smiles[parts[1]] = parts[0]
    return smiles


def _render_mol_image(smiles_str: str, size=(180, 150)):
    """Render SMILES to a matplotlib-compatible image array via RDKit Cairo."""
    try:
        import io as _io
        from rdkit import Chem
        from rdkit.Chem.Draw import rdMolDraw2D
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return None
        rdMolDraw2D.PrepareMolForDrawing(mol)
        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()
        opts.padding = 0.08
        opts.bondLineWidth = 1.8
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        buf = _io.BytesIO(drawer.GetDrawingText())
        buf.seek(0)
        return plt.imread(buf)
    except Exception:
        return None


def _fig23_multitarget(df: pd.DataFrame, out_dir: Path, reporter) -> None:
    """Figure 23: Top 25 multi-target proteins (stacked bar)."""
    try:
        import matplotlib.patches as _mp
        from matplotlib.offsetbox import AnnotationBbox, OffsetImage

        d = _utils_mod.standardise_dataframe_tiers(df.copy(), CFG)
        d['tier'] = d[CFG.COL_TIER].astype(str)

        tier_order = [t for t in CFG.TIER_ORDER if t in TIER_PALETTE]
        tier_colors = {t: TIER_PALETTE.get(t, '#999') for t in tier_order}
        d['tier_rank'] = d['tier'].map(CFG.TIER_RANK).fillna(0)
        d['prot'] = d[CFG.COL_PROT] if CFG.COL_PROT in d.columns else d.get('Protein_Name', d.get('protein', ''))
        d['lig'] = d[CFG.COL_LIG] if CFG.COL_LIG in d.columns else d.get('Ligand_Name', d.get('ligand', ''))
        d['Ensemble_Score'] = pd.to_numeric(
            d.get('Ensemble_Score', 0), errors='coerce'
        ).fillna(0)

        best_pairs = (
            d.sort_values(['prot', 'lig', 'tier_rank', 'Ensemble_Score'],
                          ascending=[True, True, False, False])
             .groupby(['prot', 'lig'], as_index=False)
             .first()
        )
        protein_summary = best_pairs.groupby('prot').agg(
            total_pfases=('lig', 'nunique')
        ).reset_index()
        tier_counts_piv = (
            best_pairs.pivot_table(
                index='prot', columns='tier', values='lig',
                aggfunc='nunique', fill_value=0
            ).reindex(columns=tier_order, fill_value=0)
        )
        protein_summary = protein_summary.join(tier_counts_piv, on='prot').fillna(0)
        protein_summary = protein_summary.astype({t: int for t in tier_order})
        # Quality-weighted sort: sum of (tier_rank × count) then total as tiebreaker
        _tier_wt23 = {'Perfect_A': 6, 'Perfect_B': 5, 'Best_A': 4,
                      'Best_B': 3, 'Good': 2, 'Decoy': 1}
        protein_summary['_qscore'] = sum(
            protein_summary[t] * _tier_wt23.get(t, 0) for t in tier_order
        )
        protein_summary = protein_summary.sort_values(
            ['_qscore', 'total_pfases', 'Perfect_A', 'Perfect_B'],
            ascending=False
        )
        top_proteins = protein_summary.head(25).set_index('prot')

        pfas_order = [
            '25_TFA', '26_Fluoroacetate', '27_Difluoroacetate', '7_PFBA', '8_PFPeA',
            '6_PFHxA', '13_PFHpA', '1_PFOA', '4_PFNA', '10_PFDA', '12_PFUnDA',
            '11_PFDoDA', '14_PFTrDA', '17_PFTeDA', '18_PFHxDA', '19_PFODA',
            '5_PFBS', '9_PFPeS', '3_PFHxS', '15_PFHpS', '2_PFOS', '16_PFDS',
            '20_GenX', '21_ADONA', '24_C6O4', '22_6-2-FTOH', '23_8-2-FTOH',
        ]
        pfas_short  = {p: p.split('_', 1)[1] if '_' in p else p for p in pfas_order}
        smiles_map  = _load_pfas_smiles()
        pfas_images = {
            p: _render_mol_image(smiles_map.get(p, ''), size=(180, 140))
            for p in pfas_order
        }
        ligand_tier_counts = (
            best_pairs.groupby(['lig', 'tier'])
            .agg(protein_count=('prot', 'nunique'))
            .reset_index()
        )
        ligand_tier_matrix = (
            ligand_tier_counts.pivot(index='lig', columns='tier', values='protein_count')
            .reindex(pfas_order, fill_value=0)[tier_order].fillna(0)
        )
        max_total = max(int(ligand_tier_matrix.sum(axis=1).max()), 1)

        # ── Figure 23a: stacked bar — top-25 proteins ────────────────────────
        fig23a, ax = plt.subplots(figsize=(16, 9))
        bar_data = top_proteins[tier_order]
        bar_data.plot(kind='barh', stacked=True, ax=ax,
                      color=[tier_colors[t] for t in tier_order],
                      edgecolor='none', legend=False)
        ax.invert_yaxis()
        ax.set_xlabel(
            'Number of unique PFAS ligands degraded  (best tier per protein–ligand pair)',
            fontsize=10)
        ax.set_ylabel('Protein  (top 25, ranked by quality-weighted degradation breadth)',
                      fontsize=10)
        ax.tick_params(axis='y', labelsize=9)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.xaxis.grid(True, color='#DDDDDD', linewidth=0.5, zorder=0)
        ax.set_axisbelow(True)
        for i, prot in enumerate(bar_data.index):
            total  = int(bar_data.loc[prot].sum())
            pa_cnt = int(bar_data.loc[prot, 'Perfect_A']) \
                     if 'Perfect_A' in bar_data.columns else 0
            _lbl = f'{total}' if pa_cnt == 0 else f'{total}  ★×{pa_cnt}'
            ax.text(total + 0.15, i, _lbl, va='center', fontsize=8.5,
                    color=tier_colors['Perfect_A'] if pa_cnt > 0 else '#444',
                    fontweight='bold' if pa_cnt > 0 else 'normal')
        legend_h23 = [_mp.Patch(color=tier_colors[t], label=t.replace('_', ' '))
                      for t in tier_order]
        fig23a.legend(handles=legend_h23, loc='lower center', ncol=6, frameon=False,
                      fontsize=9.5,
                      title='Degrader tier  (★×N = Perfect_A count per protein)',
                      title_fontsize=9)
        # Figure title removed per request.
        plt.tight_layout(rect=[0, 0.07, 1, 1])
        plt.savefig(out_dir / "Figure_23_Top25_Multitarget_Proteins.png",
                    dpi=300, bbox_inches='tight')
        plt.close(fig23a)
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_23_Top25_Multitarget_Proteins.png').resolve()}")

    except Exception as e:
        reporter.log(f"  ! Figure 23 skipped: {e}")


def _fig24_sankey(df: pd.DataFrame, out_dir: Path, reporter) -> None:
    """
    Figure 24: Extended Sankey — multiple geometric constraints → final tier.
    Layout (up to 8 columns): All Complexes (thin) → Nuc Dist → Clamp ARG111 →
    Clamp ARG114 → Acid ASP134 → Base HIS277 → SN2 Angle → Final Tier.
    Columns absent from the data are silently dropped.
    """
    try:
        from matplotlib.path import Path as _MplPath
        from matplotlib.patches import PathPatch as _PathPatch
        import matplotlib.patches as _mp
        import matplotlib.colors as _mcol24
        import textwrap as _tw24

        d = _utils_mod.standardise_dataframe_tiers(df.copy(), CFG)
        # Drop any stale computed columns so re-assignment never raises "already exists"
        for _stale24 in ['nuc_cat', 'clamp_cat', 'stab_cat', 'triad_cat',
                         'ang_cat', 'mech_cat', 'tier_cat']:
            if _stale24 in d.columns:
                d.drop(columns=[_stale24], inplace=True)

        # ── Binning — each column mirrors a REAL production tier gate, with cut
        #    points taken straight from CFG (not hand-picked). The five redundant
        #    per-residue ligand distances are replaced by the actual rule checks:
        #    nucleophile distance · carboxylate clamp · halide stabilisation ·
        #    catalytic-triad geometry · SN2 angle · mechanistic fingerprint. ─────
        def _num24(col, default=np.inf):
            return pd.to_numeric(d[col], errors='coerce') if col in d.columns \
                else pd.Series(default, index=d.index)

        # 1) Nucleophile distance — cuts = CFG.TIER_NUC_DIST gate values
        _nd = CFG.TIER_NUC_DIST
        _nuc_bins_24 = [-np.inf, _nd['Perfect_A'], _nd['Best_A'], _nd['Best_B'], _nd['Good'], np.inf]
        _nuc_lbls_24 = [f'≤{_nd["Perfect_A"]:.1f}Å',
                        f'{_nd["Perfect_A"]:.1f}–{_nd["Best_A"]:.1f}Å',
                        f'{_nd["Best_A"]:.1f}–{_nd["Best_B"]:.1f}Å',
                        f'{_nd["Best_B"]:.1f}–{_nd["Good"]:.1f}Å',
                        f'>{_nd["Good"]:.1f}Å']
        d['nuc_cat'] = (pd.cut(_num24('Dist_Nucleophile_ASP110'),
                               bins=_nuc_bins_24, labels=_nuc_lbls_24)
                        .astype(str).fillna('Unknown')) if 'Dist_Nucleophile_ASP110' in d.columns else 'Unknown'

        # 2) Carboxylate clamp (ARG111 OR ARG114 ≤5 Å) — production clamp_ok
        d['clamp_cat'] = np.where(_num24('carboxylate_clamp_integrity', 0).fillna(0) >= 0.5,
                                  'Clamp intact', 'Clamp broken')
        # 3) Halide stabilisation (Trp/Tyr/polar sidechain ≤5.5 Å)
        d['stab_cat'] = np.where(_num24('halide_stabilisation_score', 0).fillna(0) >= 0.5,
                                 'Stabilised', 'Unstabilised')
        # 4) Catalytic-triad geometry (internal Nuc–Base & Base–Acid) — CFG cuts
        _nb24, _ba24 = _num24('dist_nuc_base_internal'), _num24('dist_base_acid_internal')
        _nbm, _bam = CFG.TIER_NB_MAX, CFG.TIER_BA_MAX
        _triad_lbls_24 = ['Triad tight', 'Triad moderate', 'Triad loose']
        d['triad_cat'] = np.where((_nb24 <= _nbm['Perfect_A']) & (_ba24 <= _bam['Perfect_A']), 'Triad tight',
                          np.where((_nb24 <= _nbm['Best_A']) & (_ba24 <= _bam['Best_A']), 'Triad moderate',
                                   'Triad loose'))
        # 5) SN2 attack angle — cuts = CFG.TIER_ANGLE_MIN
        _am = CFG.TIER_ANGLE_MIN
        _ang_ord_24 = [f'<{_am["Best_B"]:.0f}°',
                       f'{_am["Best_B"]:.0f}–{_am["Best_A"]:.0f}°',
                       f'{_am["Best_A"]:.0f}–{_am["Perfect_B"]:.0f}°',
                       f'{_am["Perfect_B"]:.0f}–{_am["Perfect_A"]:.0f}°',
                       f'≥{_am["Perfect_A"]:.0f}°']
        d['ang_cat'] = (pd.cut(_num24('SN2_Attack_Angle', -1),
                               bins=[-np.inf, _am['Best_B'], _am['Best_A'], _am['Perfect_B'], _am['Perfect_A'], np.inf],
                               labels=_ang_ord_24)
                        .astype(str).fillna('Unknown')) if 'SN2_Attack_Angle' in d.columns else 'Unknown'
        # 6) Mechanistic fingerprint score — cuts = CFG.TIER_MECH_MIN
        _mm = CFG.TIER_MECH_MIN
        _mech_lbls_24 = [f'<{_mm["Best_A"]:.1f}',
                         f'{_mm["Best_A"]:.1f}–{_mm["Perfect_B"]:.1f}',
                         f'{_mm["Perfect_B"]:.1f}–{_mm["Perfect_A"]:.1f}',
                         f'≥{_mm["Perfect_A"]:.1f}']
        d['mech_cat'] = (pd.cut(_num24('mechanistic_score', np.nan),
                                bins=[-np.inf, _mm['Best_A'], _mm['Perfect_B'], _mm['Perfect_A'], np.inf],
                                labels=_mech_lbls_24)
                         .astype(str).fillna('Unknown')) if 'mechanistic_score' in d.columns else 'Unknown'

        # Order mirrors 02_Production tier cascade: Perfect_A → … → Good → Poor → Decoy
        _tier_ord_24 = CFG.TIER_ORDER + ['Other']
        d['tier_cat'] = d[CFG.COL_TIER].astype(str)
        d.loc[~d['tier_cat'].isin(_tier_ord_24[:-1]), 'tier_cat'] = 'Other'

        # ── Column definitions (best category first → drawn on top) ────────────
        _all_intermed_24 = [
            ('nuc_cat',   _nuc_lbls_24,                'Nuc Distance',         '(Å · ASP110)',          '#1B4D2E'),
            ('clamp_cat', ['Clamp intact', 'Clamp broken'], 'Carboxylate Clamp', '(ARG111 / ARG114)',   '#2C4A1E'),
            ('stab_cat',  ['Stabilised', 'Unstabilised'],   'Halide Stabilisation', '(Trp / Tyr / polar)', '#1E4A3A'),
            ('triad_cat', _triad_lbls_24,              'Triad Geometry',       '(Nuc–Base · Base–Acid)', '#1E3A4A'),
            ('ang_cat',   list(reversed(_ang_ord_24)), 'SN2 Angle',            '(degrees)',             '#1B3A5E'),
            ('mech_cat',  list(reversed(_mech_lbls_24)), 'Mech Fingerprint',   '(score 0–1)',           '#3A1E4A'),
        ]
        _valid_im_24 = [(c, o, t, s, b) for c, o, t, s, b in _all_intermed_24
                        if (d[c] != 'Unknown').any()]

        _col_seq_24  = ['all'] + [c for c, *_ in _valid_im_24] + ['tier_cat']
        _ord_seq_24  = [None]  + [o for _, o, *_ in _valid_im_24] + [_tier_ord_24]
        _ttl_seq_24  = ['All Complexes'] + [t for _, _, t, *_ in _valid_im_24] + ['Final Tier']
        _sub_seq_24  = ['(starting pool)'] + [s for _, _, _, s, *_ in _valid_im_24] + ['(degradation tier)']
        _bg_seq_24   = ['#2E4053'] + [b for _, _, _, _, b in _valid_im_24] + ['#2E1B5E']
        _n_cols_24   = len(_col_seq_24)

        # ── Colours — sourced from CFG § 9.10 (best=green … worst=red) ─────────
        _colmap_24 = {
            'nuc_cat':   dict(zip(_nuc_lbls_24, CFG.SANKEY_GRAD5)),
            'clamp_cat': dict(CFG.SANKEY_CLAMP_COLOUR),
            'stab_cat':  dict(CFG.SANKEY_STAB_COLOUR),
            'triad_cat': dict(CFG.SANKEY_TRIAD_COLOUR),
            # angle labels ascend (<… → ≥…); colour worst→best so ≥175° is green
            'ang_cat':   dict(zip(_ang_ord_24, list(reversed(CFG.SANKEY_GRAD5)))),
            'mech_cat':  dict(zip(_mech_lbls_24, CFG.SANKEY_MECH_GRAD)),
        }
        _tier_clr_24  = {t: TIER_PALETTE.get(t, '#999')
                         for t in ['Perfect_A','Perfect_B','Best_A','Best_B','Good','Poor','Decoy','Other']}
        _tier_alp_24  = {
            'Perfect_A': 0.72, 'Perfect_B': 0.68, 'Best_A': 0.62,
            'Best_B': 0.57, 'Good': 0.48, 'Poor': 0.42, 'Decoy': 0.36, 'Other': 0.28,
        }

        def _clr_24(cat_col, lbl):
            if cat_col == 'tier_cat':
                return TIER_PALETTE.get(lbl, _tier_clr_24.get(lbl, '#CCC'))
            return _colmap_24.get(cat_col, {}).get(lbl, '#CCC')

        def _dark_24(cat_col, lbl):   # text colour now auto-picked by luminance
            return False

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
            # Stack TOP-DOWN so the first label in `order` sits at the top of the
            # column (best category on top — Perfect_A, ≤3.0Å, ≥175°, …).
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
            # Tuck both ends a hair INTO the node boxes so the junction is flush —
            # the boxes (higher zorder) cover the overlap, leaving no gap/overflow.
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
                                      edgecolor='none', zorder=zorder))
            _pct = max(_ha, _hb) / _BH_24 * 100
            if _ha > 0.033 and _hb > 0.033 and _pct >= 5.0:
                ax24.text((x0 + bw_s + x1) / 2,
                          ((y0a + _ha / 2) + (y0b + _hb / 2)) / 2,
                          f'{_pct:.0f}%',
                          ha='center', va='center', fontsize=14, fontweight='bold',
                          color='white', zorder=zorder + 1,
                          bbox=dict(boxstyle='round,pad=0.05', fc=color, ec='none', alpha=0.70))

        # ── Node box ──────────────────────────────────────────────────────────
        def _dbox_24(x, bw, y0, y1, label, count, box_color, lbl_dark=False):
            _h = max(y1 - y0, 0.001)
            # Slightly blunt corners (pad=0 → no halo; small rounding softens edges)
            ax24.add_patch(_mp.FancyBboxPatch(
                (x, y0), bw, _h,
                boxstyle='round,pad=0,rounding_size=0.004', mutation_aspect=1.0,
                facecolor=box_color, edgecolor='white',
                linewidth=1.0, alpha=0.96, zorder=5))
            if _h < 0.006:
                return
            # Auto-contrast: white text on dark boxes, dark text on light boxes,
            # chosen from the box's own luminance (fixes dark text on dark fills).
            _rr, _gg, _bb = _mcol24.to_rgb(box_color)
            _lum = 0.299 * _rr + 0.587 * _gg + 0.114 * _bb
            lbl_col = 'white' if _lum < 0.6 else '#111111'
            # Show the actual COUNT (not %) — a rare class like Perfect_A (n=8) must
            # read "8", never a rounded "0.0%".
            _cnt_s = f'{int(count):,}'
            # Single line "label  count" (wrap only if long). Font scales with box
            # height so the value FILLS the otherwise-empty top/bottom space; thin
            # boxes keep a readable minimum and simply overlay the box ("on top").
            _txt = _tw24.fill(f'{label}  {_cnt_s}', width=15)
            _fs  = float(np.clip(_h * 520, 9.0, 16.0))
            ax24.text(x + bw / 2, y0 + _h / 2, _txt,
                      ha='center', va='center', fontsize=_fs, fontweight='bold',
                      color=lbl_col, zorder=10, linespacing=1.0)

        # ── Figure ────────────────────────────────────────────────────────────
        fig24, ax24 = plt.subplots(figsize=(30, 16))
        ax24.set_xlim(0, 1); ax24.set_ylim(0, 1.05); ax24.axis('off')
        fig24.patch.set_facecolor('white')

        # ── Column 0: thin "All Complexes" bar (90°-rotated label) ────────────
        ax24.add_patch(_mp.FancyBboxPatch(
            (_xs_24[0], _BY0_24), _bw0_24, _BH_24,
            boxstyle='round,pad=0,rounding_size=0.004', mutation_aspect=1.0,
            facecolor='#2E4053', edgecolor='white', linewidth=1.0, alpha=0.96, zorder=5))
        ax24.text(_xs_24[0] + _bw0_24 / 2, _BY0_24 + _BH_24 / 2,
                  f'All Complexes  {_total_j_24:,}',
                  ha='center', va='center', fontsize=15, fontweight='bold',
                  color='white', zorder=6, rotation=90)

        # ── Node boxes for columns 1..n_cols-1 ───────────────────────────────
        for _ci24, _cc24 in enumerate(_col_seq_24[1:], start=1):
            for lbl, (y0, y1) in _cpos_24.get(_cc24, {}).items():
                _dbox_24(_xs_24[_ci24], _bw_24, y0, y1, lbl,
                         int(_ccnt_24.get(_cc24, {}).get(lbl, 0)),
                         _clr_24(_cc24, lbl), _dark_24(_cc24, lbl))

        # ── Column header banners — extra gap so they clear the node content ──
        _hdr_gap_24 = 0.078   # wider gap → more separation between content and headers
        _hdr_y_24   = _BY0_24 + _BH_24 + _hdr_gap_24
        _hbox_24    = []      # (left, right) x-extent of each header banner
        for _ci24, (_htxt, _hsub, _hbg) in enumerate(zip(_ttl_seq_24, _sub_seq_24, _bg_seq_24)):
            if _ci24 == 0:
                # thin source bar → give its header a wider banner with horizontal,
                # wrapped text so the full title shows (no rotation/clipping)
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
                boxstyle='round,pad=0,rounding_size=0.005', mutation_aspect=1.0,
                facecolor=_hbg, edgecolor='none', alpha=0.95, zorder=7))
            ax24.text(_hx + _cbw / 2, _hdr_y_24 + 0.038, _title,
                      ha='center', va='center', fontsize=14, fontweight='bold',
                      color='white', zorder=8, linespacing=1.0)
            ax24.text(_hx + _cbw / 2, _hdr_y_24 + 0.012, _hsub,
                      ha='center', va='center', fontsize=12, color='#CFE0EA', zorder=8)

        # ── Flow arrows between headers — anchored exactly at banner edges ────
        _ya24 = _hdr_y_24 + 0.030
        for _ci24 in range(_n_cols_24 - 1):
            # Tail flush on this banner's right edge, head flush on next banner's
            # left edge — arrow touches both boxes (no gap, matches box shape).
            _x_start = _hbox_24[_ci24][1]
            _x_end   = _hbox_24[_ci24 + 1][0]
            if _x_end <= _x_start:
                continue
            ax24.annotate('', xy=(_x_end, _ya24), xytext=(_x_start, _ya24),
                          arrowprops=dict(arrowstyle='->', color='#888', lw=1.4,
                                          mutation_scale=12, shrinkA=0, shrinkB=0),
                          zorder=9)

        # ── Ribbons col 0 → col 1 ─────────────────────────────────────────────
        _fc24    = _col_seq_24[1]
        _fo24    = _ord_seq_24[1]
        _flow01  = d.groupby([_fc24, 'tier_cat']).size().reset_index(name='count')
        _flow01['_nr'] = _flow01[_fc24].apply(lambda v: _fo24.index(v) if v in _fo24 else 999)
        _flow01['_tr'] = _flow01['tier_cat'].map(_rank_m_24).fillna(-99)
        _flow01  = _flow01.sort_values(['_nr', '_tr'])
        _used_all_24 = _BY0_24
        _dest01_24   = {k: v[0] for k, v in _cpos_24.get(_fc24, {}).items()}
        for _, r01 in _flow01.iterrows():
            if r01['count'] <= 0: continue
            nc, tc = r01[_fc24], r01['tier_cat']
            if nc not in _cpos_24.get(_fc24, {}): continue
            _h_d = _cpos_24[_fc24][nc][1] - _cpos_24[_fc24][nc][0]
            _rhs = _BH_24 * r01['count'] / _total_j_24
            _rht = _h_d  * r01['count'] / max(int(_ccnt_24.get(_fc24, {}).get(nc, 1)), 1)
            y0a, y0b = _used_all_24, _dest01_24.get(nc, 0)
            _used_all_24 = y0a + _rhs
            _dest01_24[nc] = y0b + _rht
            _brib_24(_xs_24[0], _bw0_24, _xs_24[1],
                     y0a, y0a + _rhs, y0b, y0b + _rht,
                     TIER_PALETTE.get(tc, _tier_clr_24.get(tc, '#999')),
                     _tier_alp_24.get(tc, 0.40), zorder=2)

        # ── Ribbons col i → col i+1 (i = 1 .. n_cols-2) ─────────────────────
        for _ri in range(1, _n_cols_24 - 1):
            _sc24 = _col_seq_24[_ri]
            _tc24 = _col_seq_24[_ri + 1]
            _sp24 = _cpos_24.get(_sc24, {})
            _tp24 = _cpos_24.get(_tc24, {})
            _sn24 = _ccnt_24.get(_sc24, {})
            _tn24 = _ccnt_24.get(_tc24, {})
            _grp_cols_24 = [_sc24, _tc24] if _tc24 == 'tier_cat' else [_sc24, _tc24, 'tier_cat']
            _frm  = (d.groupby(_grp_cols_24)
                     .size().reset_index(name='count'))
            if 'tier_cat' not in _frm.columns:
                _frm['tier_cat'] = _frm[_tc24]
            _frm['_tr'] = _frm['tier_cat'].map(_rank_m_24).fillna(-99)
            _frm = _frm.sort_values('_tr')
            _so24 = {k: v[0] for k, v in _sp24.items()}
            _to24 = {k: v[0] for k, v in _tp24.items()}
            for _, rr in _frm.iterrows():
                if rr['count'] <= 0: continue
                s0, t0, tc = rr[_sc24], rr[_tc24], rr['tier_cat']
                if s0 not in _sp24 or t0 not in _tp24: continue
                _hs = _sp24[s0][1] - _sp24[s0][0]
                _ht = _tp24[t0][1] - _tp24[t0][0]
                _rhs = _hs * rr['count'] / max(int(_sn24.get(s0, 1)), 1)
                _rht = _ht * rr['count'] / max(int(_tn24.get(t0, 1)), 1)
                y0a, y0b = _so24.get(s0, 0), _to24.get(t0, 0)
                _so24[s0] = y0a + _rhs
                _to24[t0] = y0b + _rht
                _brib_24(_xs_24[_ri], _bw_24, _xs_24[_ri + 1],
                         y0a, y0a + _rhs, y0b, y0b + _rht,
                         TIER_PALETTE.get(tc, _tier_clr_24.get(tc, '#999')),
                         _tier_alp_24.get(tc, 0.40), zorder=2 + _ri)

        # Legend removed — Final Tier column already carries the tier names/colours.
        # Figure title removed per request.
        plt.tight_layout(pad=0.3)
        plt.savefig(out_dir / "Figure_24_Sankey_Workflow.png",
                    dpi=300, bbox_inches='tight', facecolor='white')
        plt.close(fig24)
        reporter.log(f"  ✔ Saved: {(out_dir / 'Figure_24_Sankey_Workflow.png').resolve()}")

    except Exception as e:
        reporter.log(f"  ! Figure 24 skipped: {e}")


def _fig25_pfas_size(df: pd.DataFrame, out_dir: Path, reporter) -> None:
    """Figure 25a–25e — PFAS chain-length selectivity saved as five separate PNGs."""
    try:
        _req25 = ['total_fluorine_count', 'SN2_Attack_Angle', 'Boltz_Model_Confidence', 'degrader_tier']
        _miss25 = [c for c in _req25 if c not in df.columns]
        if _miss25:
            reporter.log(f"  ! Figure 25 skipped: missing columns {_miss25}")
            return

        _d25 = df[_req25].copy()
        for _c25 in _req25[:-1]:
            _d25[_c25] = pd.to_numeric(_d25[_c25], errors='coerce')
        _d25 = _d25.dropna(subset=_req25)
        _d25 = _d25[_d25['degrader_tier'].isin(TIER_ORDER_LOGIC)]
        if len(_d25) < 30:
            reporter.log("  ! Figure 25 skipped: insufficient data after filtering")
            return

        # PFAS chain-length bins (F-count as proxy: C4≈9F, C6≈13F, C8≈17F, C10≈21F)
        _f25_bins   = [0, 8, 13, 18, 24, np.inf]
        _f25_labels = [
            'Very Short\n(≤C4, ≤8F)',
            'Short-chain\n(C5–6, 9–13F)',
            'PFOA/PFOS\n(C7–9, 14–18F)',
            'Long-chain\n(C10–12, 19–24F)',
            'Ultra-long\n(≥C13, ≥25F)',
        ]
        _f25_labels_clean = ['≤C4 (≤8F)', 'C5–6 (9–13F)', 'C7–9 (14–18F)',
                             'C10–12 (19–24F)', '≥C13 (≥25F)']
        _d25['_size_bin'] = pd.cut(
            _d25['total_fluorine_count'], bins=_f25_bins,
            labels=_f25_labels, right=True
        )
        _d25 = _d25.dropna(subset=['_size_bin'])

        # Mechanistic outcome classification
        _ang25  = _d25['SN2_Attack_Angle']
        _conf25 = _d25['Boltz_Model_Confidence']
        _SA25, _IA25, _SC25 = CFG.SUBSTRATE_ANGLE_MIN, CFG.INHIBITOR_ANGLE_MAX, CFG.SUBSTRATE_CONF_MIN
        _d25['_outcome'] = np.where(
            (_ang25 >= _SA25) & (_conf25 >= _SC25), 'Substrate',
            np.where(
                (_ang25 < _IA25) & (_conf25 >= _SC25), 'Potential Inhibitor',
                np.where(
                    (_ang25 >= _SA25) & (_conf25 < _SC25), 'Reactive (low conf)',
                    np.where((_ang25 >= _IA25) & (_ang25 < _SA25), 'Borderline', 'Non-reactive')
                )
            )
        )
        _oc25_order  = ['Substrate', 'Borderline', 'Reactive (low conf)',
                        'Non-reactive', 'Potential Inhibitor']
        _oc25_colors = dict(CFG.OUTCOME_COLOUR)   # sourced from CFG § 9.9

        # ── Panel A (25a): hexbin landscape + tier scatter + rolling median ──
        fig25a, axA = plt.subplots(figsize=(14, 8))
        axA.axhspan(CFG.SUBSTRATE_ANGLE_MIN, 185, color='#009E73', alpha=0.08, zorder=0)
        axA.axhspan(90, CFG.INHIBITOR_ANGLE_MAX, color='#D55E00', alpha=0.08, zorder=0)
        axA.axhline(CFG.SUBSTRATE_ANGLE_MIN, color='#009E73', lw=1.2, ls='--', alpha=0.7, zorder=2)
        axA.axhline(CFG.INHIBITOR_ANGLE_MAX, color='#D55E00', lw=1.2, ls='--', alpha=0.7, zorder=2)
        _hb25 = axA.hexbin(
            _d25['total_fluorine_count'], _d25['SN2_Attack_Angle'],
            gridsize=50, mincnt=1, cmap='YlOrRd', alpha=0.65, zorder=1, linewidths=0.2
        )
        _cb25 = fig25a.colorbar(_hb25, ax=axA, pad=0.01, aspect=30, shrink=0.85)
        _cb25.set_label('Count per hex', fontsize=9)
        _cb25.ax.tick_params(labelsize=8)
        _samp_A25 = _d25.sample(min(3000, len(_d25)), random_state=42)
        for _t25a in sorted(TIER_ORDER_LOGIC, key=lambda t: t == 'Perfect_A'):
            _is_pa25a = (_t25a == 'Perfect_A')
            # Perfect_A is crucial and rare → plot ALL of its points (never sampled)
            _src25a = _d25 if _is_pa25a else _samp_A25
            _ts25 = _src25a[_src25a['degrader_tier'] == _t25a]
            if _ts25.empty:
                continue
            axA.scatter(
                _ts25['total_fluorine_count'], _ts25['SN2_Attack_Angle'],
                color=TIER_PALETTE.get(_t25a, '#999'),
                alpha=0.85 if _is_pa25a else 0.45,
                s=60 if _is_pa25a else 14,
                edgecolors='black' if _is_pa25a else 'none',
                linewidths=0.8 if _is_pa25a else 0,
                zorder=5 if _is_pa25a else 3, label=_t25a
            )
        _roll25 = _d25[['total_fluorine_count', 'SN2_Attack_Angle']].sort_values('total_fluorine_count')
        if len(_roll25) >= 20:
            _win25 = max(20, len(_roll25) // 40)
            _roll25['_med'] = _roll25['SN2_Attack_Angle'].rolling(
                window=_win25, center=True, min_periods=5).median()
            axA.plot(_roll25['total_fluorine_count'], _roll25['_med'],
                     color='black', lw=2.2, zorder=6, label='Rolling median SN2')
        axA.text(0.99, 0.92, 'SUBSTRATE ZONE  (SN2 ≥ 165°)',
                 transform=axA.transAxes, ha='right', va='top', fontsize=9.5,
                 color='#007A52', fontweight='bold',
                 bbox=dict(facecolor='white', alpha=0.75, edgecolor='none'))
        axA.text(0.99, 0.10, 'POTENTIAL INHIBITOR ZONE  (SN2 < 145°)',
                 transform=axA.transAxes, ha='right', va='bottom', fontsize=9.5,
                 color='#C04000', fontweight='bold',
                 bbox=dict(facecolor='white', alpha=0.75, edgecolor='none'))
        for _fb25 in [8, 13, 18, 24]:
            axA.axvline(_fb25, color='#444', lw=0.8, ls=':', alpha=0.5, zorder=2)
        _bin_xmids25   = [4, 10.5, 15.5, 21, 28]
        _bin_toplbls25 = ['Very Short', 'Short chain', 'PFOA/PFOS', 'Long chain', 'Ultra-long']
        axA.set_ylim(85, 186)
        _fc_max25 = int(_d25['total_fluorine_count'].max()) + 2
        axA.set_xticks(range(0, _fc_max25 + 1, 2))
        for _bx25, _bl25 in zip(_bin_xmids25, _bin_toplbls25):
            axA.text(_bx25, 184, _bl25, ha='center', va='top',
                     fontsize=7.5, color='#444', style='italic')
        axA.set_xlabel('Total fluorine count  (proxy for carbon chain length)', fontsize=11)
        axA.set_ylabel('SN2 Attack Angle (°)', fontsize=11)
        axA.legend(loc='upper left', bbox_to_anchor=(0.0, 1.05),
                   ncol=len(TIER_ORDER_LOGIC) + 1, fontsize=7,
                   framealpha=0.92, fancybox=True, handlelength=0.8, borderpad=0.3)
        plt.tight_layout()
        _out25a = out_dir / "Figure_25a_PFAS_Size_Hexbin_Landscape.png"
        fig25a.savefig(_out25a, dpi=300, bbox_inches='tight')
        plt.close(fig25a)
        reporter.log(f"  ✔ Saved: {_out25a.resolve()}")

        # ── Figure 25b: ONE figure, two stories per chain-length bin — a left
        #    stacked bar (mechanistic outcome) and a right stacked bar (degrader
        #    tier) placed side by side within each bin group. ──────────────────
        _oc_ct25 = (
            _d25.groupby(['_size_bin', '_outcome'], observed=True)
            .size().unstack(fill_value=0)
        )
        _oc_ct25 = _oc_ct25.reindex(columns=[c for c in _oc25_order if c in _oc_ct25.columns])
        _oc_pct25 = _oc_ct25.div(_oc_ct25.sum(axis=1), axis=0) * 100

        _tier_ct25 = (
            _d25.groupby(['_size_bin', 'degrader_tier'], observed=True)
            .size().unstack(fill_value=0)
        )
        _tier_ct25 = _tier_ct25.reindex(columns=[t for t in TIER_ORDER_LOGIC if t in _tier_ct25.columns])
        _tier_pct25 = _tier_ct25.div(_tier_ct25.sum(axis=1), axis=0) * 100

        _xlabels_25 = [_f25_labels_clean[_f25_labels.index(lb)] if lb in _f25_labels else str(lb)
                       for lb in _oc_pct25.index.tolist()]

        # Outcome palette + per-bin colours — sourced from CFG (§ 9.9, § 9.9b).
        _oc25_colors_b = dict(CFG.OUTCOME_COLOUR)
        _bin_palette_25 = list(CFG.PFAS_SIZE_BIN_COLOUR)
        _bin_cols_25 = [_bin_palette_25[i % len(_bin_palette_25)] for i in range(len(_oc_pct25))]

        fig25b, axB = plt.subplots(figsize=(15, 8))
        _xB, _wB, _offB = np.arange(len(_oc_pct25)), 0.27, 0.19   # bars 75% of prior width

        # Master container bar per x-category (like Fig 19a) — faint bin-coloured
        # fill + matching coloured outline, holding both child stacked bars.
        from matplotlib.colors import to_rgba as _to_rgba25
        for _xi, _bc in zip(_xB, _bin_cols_25):
            axB.bar(_xi, 104, 0.74, bottom=0, color=_to_rgba25(_bc, 0.10),
                    edgecolor=_bc, linewidth=1.4, zorder=2)

        def _stack_at_25(xpos, pct_df, ct_df, colour_fn, side):
            # Inside label when the segment is tall enough (≥4%); otherwise a small
            # leader-arrow points to the thin segment with the value placed in the
            # gap beside the bar (left bar → labels left, right bar → labels right).
            _bot = np.zeros(len(pct_df))
            _dx  = -0.30 if side == 'left' else 0.30
            _ha  = 'right' if side == 'left' else 'left'
            _esign = -1 if side == 'left' else 1
            for _col in pct_df.columns:
                _vals = pct_df[_col].values
                _ccol = colour_fn(_col)
                axB.bar(xpos, _vals, _wB, bottom=_bot, color=_ccol,
                        edgecolor='white', linewidth=0.6, zorder=3)
                for _xi, (_v, _b) in enumerate(zip(_vals, _bot)):
                    if _v <= 0:
                        continue
                    _ymid = _b + _v / 2
                    if _v >= 4:
                        axB.text(xpos[_xi], _ymid, f'{_v:.0f}%', ha='center',
                                 va='center', fontsize=8.5, color='white',
                                 fontweight='bold', zorder=4)
                    elif _v >= 0.8:
                        axB.annotate(f'{_v:.0f}%',
                                     xy=(xpos[_xi] + _esign * _wB / 2, _ymid),
                                     xytext=(xpos[_xi] + _dx, _ymid),
                                     ha=_ha, va='center', fontsize=8.5,
                                     color=_ccol, fontweight='bold', zorder=7,
                                     arrowprops=dict(arrowstyle='-', color=_ccol,
                                                     lw=0.5, shrinkA=1, shrinkB=1))
                _bot = _bot + _vals
            for _xi, _n in enumerate(ct_df.sum(axis=1).values):
                axB.text(xpos[_xi], 101.5, f'n={_n:,}', ha='center', va='bottom',
                         fontsize=6.3, color=_bin_cols_25[_xi], fontweight='bold')

        _stack_at_25(_xB - _offB, _oc_pct25, _oc_ct25, lambda c: _oc25_colors_b.get(c, '#999'), 'left')
        _stack_at_25(_xB + _offB, _tier_pct25, _tier_ct25, lambda t: TIER_PALETTE.get(t, '#999'), 'right')
        # Mini-tags under each paired bar so the two stories are unambiguous
        for _xi, _bc in zip(_xB, _bin_cols_25):
            axB.text(_xi - _offB, -1.5, 'outcome', ha='center', va='top',
                     fontsize=8.5, color=_bc, style='italic', fontweight='bold')
            axB.text(_xi + _offB, -1.5, 'tier', ha='center', va='top',
                     fontsize=8.5, color=_bc, style='italic', fontweight='bold')

        axB.set_ylim(0, 108)
        axB.set_xlim(-0.6, len(_oc_pct25) - 0.4)
        axB.set_xticks(_xB)
        axB.set_xticklabels(_xlabels_25, fontsize=9)
        for _tk, _bc in zip(axB.get_xticklabels(), _bin_cols_25):
            _tk.set_color(_bc)
            _tk.set_fontweight('bold')
        axB.tick_params(axis='x', pad=16)   # room for the outcome/tier mini-tags
        axB.set_ylabel('Percentage of complexes (%)', fontsize=10)
        axB.set_xlabel('PFAS chain-length bin', fontsize=10)
        axB.spines['top'].set_visible(False)
        axB.spines['right'].set_visible(False)
        axB.yaxis.grid(True, color='#ECECEC', linewidth=0.7, zorder=0)
        axB.set_axisbelow(True)

        from matplotlib.patches import Patch as _P25b
        _oc_h25 = [_P25b(facecolor=_oc25_colors_b.get(c, '#999'), edgecolor='white',
                         linewidth=0.5, label=c) for c in _oc_pct25.columns]
        _ti_h25 = [_P25b(facecolor=TIER_PALETTE.get(t, '#999'), edgecolor='white',
                         linewidth=0.5, label=t) for t in _tier_pct25.columns]
        # Single merged legend in ONE row, top-left — plain (no title), normal font.
        _all_h25 = _oc_h25 + _ti_h25
        axB.legend(handles=_all_h25, labels=[h.get_label() for h in _all_h25],
                   loc='lower left', bbox_to_anchor=(0.0, 1.005),
                   ncol=len(_all_h25), fontsize=8.5,
                   framealpha=0.92, fancybox=True, handlelength=1.0,
                   handletextpad=0.4, columnspacing=0.8, borderpad=0.4)
        plt.tight_layout()
        _out25b = out_dir / "Figure_25b_PFAS_Size_Composition_Merged.png"
        fig25b.savefig(_out25b, dpi=300, bbox_inches='tight')
        plt.close(fig25b)
        reporter.log(f"  ✔ Saved: {_out25b.resolve()}")

        # ── Figure 25c: confidence distribution boxplot per chain-length bin ──
        # (Tier-composition panel moved into 25b; 25c is now a single box plot.)
        fig25c, axC = plt.subplots(figsize=(12, 6))
        _c25_plot_df = _d25[['_size_bin', 'Boltz_Model_Confidence']].copy()
        _c25_order   = [lb for lb in _f25_labels if lb in _d25['_size_bin'].cat.categories]
        _c25_clean   = [_f25_labels_clean[_f25_labels.index(lb)] if lb in _f25_labels else str(lb)
                        for lb in _c25_order]
        _c25_plot_df['_bin_clean'] = _c25_plot_df['_size_bin'].map(
            {lb: cl for lb, cl in zip(_c25_order, _c25_clean)}
        )
        sns.boxplot(data=_c25_plot_df, x='_bin_clean', y='Boltz_Model_Confidence',
                    order=_c25_clean, ax=axC,
                    color='#56B4E9', width=0.55, linewidth=1.2,
                    fliersize=0, showfliers=False, zorder=3)
        # Stripplot — sample max 300 per bin to avoid overplotting
        _strip_samp = pd.concat([
            _c25_plot_df[_c25_plot_df['_bin_clean'] == cl].sample(min(300, int((_c25_plot_df['_bin_clean'] == cl).sum())), random_state=42)
            for cl in _c25_clean if ((_c25_plot_df['_bin_clean'] == cl)).any()
        ])
        sns.stripplot(data=_strip_samp, x='_bin_clean', y='Boltz_Model_Confidence',
                      order=_c25_clean, ax=axC,
                      color='#0072B2', alpha=0.25, size=2.5, jitter=True, zorder=2)
        # Median trend line
        _meds_C25 = [float(_c25_plot_df[_c25_plot_df['_bin_clean'] == cl]['Boltz_Model_Confidence'].median())
                     for cl in _c25_clean]
        axC.plot(np.arange(len(_meds_C25)), _meds_C25,
                 color='#D55E00', lw=2.0, marker='D', markersize=6,
                 zorder=6, label='Median confidence')
        # The 0.75 substrate-confidence threshold lies below the whole data range,
        # so it is noted in the legend instead of drawn as an off-axis line.
        from matplotlib.lines import Line2D as _L25c
        _thr_handle25 = _L25c([0], [0], color='#009E73', lw=1.2, ls='--',
                              label='Substrate conf. threshold = 0.75  (all bins above)')
        axC.set_ylim(0.86, 1.005)   # tight to the data → minimal white space
        axC.set_ylabel('Boltz Model Confidence', fontsize=9.5)
        axC.set_xlabel('PFAS chain-length bin', fontsize=10)
        axC.tick_params(axis='x', labelsize=8.5)
        axC.spines['top'].set_visible(False)
        axC.spines['right'].set_visible(False)
        _med_h25, _med_l25 = axC.get_legend_handles_labels()
        axC.legend(_med_h25 + [_thr_handle25], _med_l25 + [_thr_handle25.get_label()],
                   loc='lower left', fontsize=7.5, framealpha=0.90)
        plt.tight_layout()
        _out25c = out_dir / "Figure_25c_PFAS_Size_Confidence_Boxplot.png"
        fig25c.savefig(_out25c, dpi=300, bbox_inches='tight')
        plt.close(fig25c)
        reporter.log(f"  ✔ Saved: {_out25c.resolve()}")

    except Exception as e:
        reporter.log(f"  ! Figure 25 skipped: {e}")


# ===============================================================================
# SECTION 5: FIGURE DESCRIPTIONS & REPORTING
# ===============================================================================

def write_figure_descriptions(out_dir: Path):
    """
    Writes a human-readable log file describing every figure produced by this pipeline.

    Figures are listed in scientific narrative order:
      Part 1  (01–03): Dataset Overview
      Part 2  (04–05): AI Prediction Quality
      Part 3  (06–08): Structural Validation
      Part 4  (09–12): Mechanistic Analysis
      Part 5  (13–14): Ligand Interactions
      Part 6  (15–16): Binding Energetics
      Part 7  (17):    Chemical Space
      Part 8  (18–22): Multi-metric Synthesis
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
        "Figure 01 — Figure_01_Tier_Distribution.png",
        "  Title   : Distribution of Catalytic Tiers",
        "  Type    : Vertical bar chart — count + percentage annotations inside bars",
        "  X-axis  : Tier name (ordered best→worst: Perfect_A, Perfect_B, Best_A, Best_B, Good, Poor)",
        "  Y-axis  : Count of protein–ligand complexes",
        "  Inset   : Pie chart showing Boltz-2 diffusion model selection frequency (model_0–model_4)",
        "  Look for: The fraction of candidates reaching each quality level.",
        "            Inset reveals whether one Boltz model dominates (~model_0 ≈ 70%).",
        "",
        "-" * 80,
        "Figure 02 — Figure_02_Alignment_Grades.png",
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
        "Figure 03 — Figure_03_Tier_Grade_Distribution.png",
        "  Title   : Tier x Alignment Grade — Stacked 100% Horizontal Bar Chart",
        "  Type    : 100% stacked horizontal bars (one row per tier)",
        "  Y-axis  : Degrader tiers (each row normalised to 100%)",
        "  X-axis  : Grade bands (A [>=90%] to I [<20%]); 10% intervals with pastel background stripes",
        "  Colour  : Grade-specific green-to-brown ramp (same palette as Figure 02)",
        "  Annotations: Count + % for segments >=2.5%",
        "  Look for: Do higher-tier proteins cluster at higher sequence identity grades?",
        "            If Perfect/Best rows are dominated by Grade A-C, conservation tracks quality.",
        "",
        "=" * 80,
        "PART 2 — AI PREDICTION QUALITY",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 04a — Figure_04a_AI_Quality_Assessment.png",
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
        "Figure 04b — Figure_04b_PA_AI_Quality_Space.png",
        "  Title   : Perfect_A AI Quality Space (hexbin density + structure thumbnails)",
        "  Type    : Hexbin density scatter of pTM vs. ipTM; inset PyMOL structure thumbnails",
        "  Axes    : X = pTM (global fold confidence, 0-1); Y = ipTM (interface confidence, 0-1)",
        "  Stars   : Mark Perfect_A PFAS positions in the pTM/ipTM confidence space",
        "  Insets  : PyMOL-rendered protein-ligand structures for each Perfect_A representative",
        "  Look for: Perfect_A stars in the top-right quadrant (high pTM AND high ipTM).",
        "            Structure thumbnails show the active-site geometry at a glance.",
        "",
        "-" * 80,
        "Figure 05 — Figure_05_pTM_vs_ipTM_by_Tier.png",
        "  Title   : pTM vs. ipTM — 2-D Scatter by Tier",
        "  Type    : Scatter with diagonal reference line (pTM = ipTM) + coloured deviation zones",
        "  X-axis  : pTM — overall protein fold confidence (0-1); 1.0 = perfect fold prediction",
        "  Y-axis  : ipTM — interface confidence (0-1); 1.0 = perfect interface prediction",
        "  Diagonal: pTM = ipTM reference; +/-0.03 band; above diagonal = preferred (interface > fold)",
        "  Colour  : Catalytic tier; diamond markers show per-tier median positions",
        "  Zones   : Green band above diagonal = interface-dominant (preferred); red band below = fold-dominant",
        "  Look for: Best-tier points cluster above the diagonal and in the top-right corner.",
        "            Decoy complexes often fall below the diagonal (fold confident, interface uncertain).",
        "",
        "=" * 80,
        "PART 3 — STRUCTURAL VALIDATION",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 06 — Figure_06_ActiveSite_RMSD_by_Tier.png",
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
        "Figure 07 — Figure_07_Feature_Correlations.png",
        "  Title   : Feature Correlation Matrix (Spearman rho) with Significance & Clustering",
        "  Type    : Triangular heatmap (lower triangle only) with hierarchical clustering",
        "  Cells   : Spearman rank correlation coefficient (-1 to +1) + significance stars",
        "  Stars   : *** p<0.001   ** p<0.01   * p<0.05   (no star = not significant)",
        "  Order   : Features reordered by hierarchical clustering (similar features grouped together)",
        "  Look for: Strong red cells = features that rise and fall together (redundant or causal).",
        "            Strong blue = features that are inversely related.",
        "            AI Confidence + Binding Prob. cluster together — same underlying signal.",
        "            Pareto Rank inverts against physics/AI scores (lower rank = better = higher score).",
        "",
        "-" * 80,
        "Figure 08 — Figure_08_Tier_Quality_DotPlot.png",
        "  Title   : Tier Quality Summary — Normalised Multi-metric Cleveland Dot Plot",
        "  Type    : Horizontal Cleveland dot plot; one metric per row, tiers as coloured dots",
        "  Y-axis  : 6 normalised metrics (0=worst, 1=best within each metric):",
        "            AI Confidence, Mech. Fingerprint, Binding Prob., Seq. Identity, SN2 Angle, RMSD (inverted)",
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
        "Figure 09 — Figure_09_Mech_State_CrossTab.png",
        "  Title   : Mechanistic State Cross-Tab (Halide Stabilisation x Carboxylate Clamp by Tier)",
        "  Type    : Heatmap (count + row %) — rows = tiers; columns = mechanistic state combinations",
        "  States  : Stabilised/Unstabilised (TRP/TYR aromatic shield) x Clamped/Unclamped (ASP carboxylate)",
        "  Look for: Virtually all complexes across all tiers show Stabilised+Clamped (98-100%).",
        "            This confirms the DEHA4 active-site architecture is robustly conserved.",
        "            Any tier deviating from 100% warrants structural inspection.",
        "",
        "-" * 80,
        "Figure 10 — Figure_10_Mechanistic_Fingerprint_by_Tier.png",
        "  Title   : Mechanistic Fingerprint Score by Tier",
        "  Type    : Mean +/- 95% CI dot plot with jittered individual data points behind",
        "  X-axis  : Catalytic tier (ordered best to worst)",
        "  Y-axis  : Mechanistic Fingerprint Score (0-1; fraction of 8 catalytic residues within 6 Ang of ligand)",
        "  Zones   : Green (>=0.88 = strong), Yellow (0.50-0.88 = moderate), Red (<0.50 = weak)",
        "  Annotations: Large dot = mean; vertical bar = 95% CI; % label = fraction in strong zone",
        "  Look for: Most complexes score >=0.88 across all tiers — active-site architecture conserved.",
        "            Poor tier has the widest spread (some missing residues in contact).",
        "",
        "-" * 80,
        "Figure 11 — Figure_11_SN2_Angle_by_Tier.png",
        "  Title   : SN2 Attack Angle — Empirical Cumulative Distribution (ECDF) by Tier",
        "  Type    : ECDF step plot per tier",
        "  X-axis  : SN2 attack angle (degrees); 180 deg = ideal linear nucleophilic back-attack",
        "            Vertical dashed lines = tier-specific angle thresholds (145, 155, 165, 175 deg)",
        "  Y-axis  : Cumulative fraction of complexes in each tier with angle <= X (read as %)",
        "  Dots    : Median angle per tier",
        "  How to read: Y% of complexes in this tier have SN2 angle <= X degrees",
        "  Look for: Perfect/Best tiers: ECDF curves shifted right (most complexes at high angles).",
        "            Good/Poor: curves shifted left. Tier thresholds show what % meet each criterion.",
        "",
        "-" * 80,
        "Figure 12a — Figure_12a_Mechanism_Geometry_Scatter.png",
        "  Title   : SN2 Attack Angle vs. Nucleophile-Ligand Distance",
        "  Type    : Scatter plot with tier threshold lines and ideal-zone shading",
        "  X-axis  : ASP110 nucleophile to electrophilic carbon distance (Ang); shorter = reaction-ready",
        "            Clipped to <=5 Ang to remove outliers with no catalytic contact.",
        "  Y-axis  : SN2 attack angle (degrees); 180 = perfect linear back-attack geometry",
        "  Shape   : Circle = halide-stabilised (TRP/TYR aromatic shield present); X = not stabilised",
        "  Dashed lines: Vertical = distance thresholds per tier; Horizontal = angle thresholds.",
        "  Shading : Green zone = ideal geometry (short distance AND high angle).",
        "  Look for: Perfect/Best tier points cluster top-left (short distance + high angle).",
        "            Good/Poor tier points scatter widely — geometry less constrained.",
        "",
        "-" * 80,
        "Figure 12b — Figure_12b_PA_Mechanistic_Quality_Space.png",
        "  Title   : Perfect_A Mechanistic Quality Space (hexbin density + structure thumbnails)",
        "  Type    : Hexbin density scatter; inset PyMOL protein-ligand structure thumbnails",
        "  Axes    : X = Dist_ASP110 (Ang); Y = SN2_Attack_Angle (degrees)",
        "  Stars   : Highlight Perfect_A PFAS complexes within the mechanistic space",
        "  Insets  : PyMOL-rendered active-site views for each Perfect_A representative protein",
        "  Look for: Perfect_A stars clustered in the ideal mechanistic zone (top-left corner).",
        "            Insets confirm the geometry seen in data is reflected in the 3-D structure.",
        "",
        "=" * 80,
        "PART 5 — LIGAND INTERACTIONS",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 13a — Figure_13a_Molecular_Interaction_Profile.png",
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
        "            Panel B: Perfect/Best tiers engage a higher fraction of ligand fluorine.",
        "",
        "-" * 80,
        "Figure 13b — Figure_13b_PA_Interaction_Quality_Space.png",
        "  Title   : Perfect_A Interaction Quality Space (hexbin density + structure thumbnails)",
        "  Type    : Hexbin density scatter of Binding_Prob vs. ipTM; inset PyMOL structures",
        "  Axes    : X = Binding Probability; Y = ipTM (interface confidence)",
        "  Stars   : Mark Perfect_A PFAS positions in the binding-confidence 2-D space",
        "  Insets  : PyMOL-rendered protein-ligand structures for each Perfect_A representative",
        "  Look for: Perfect_A stars clustered top-right (high binding probability + high interface confidence).",
        "",
        "-" * 80,
        "Figure 14 — Figure_14_Fluorine_Engagement_by_Tier.png",
        "  Title   : Fluorine Engagement Ratio by Tier (box + strip + median trend line)",
        "  Type    : Dual-axis (box + strip plot on left; median trend line on right)",
        "  Left axis : Box + strip of Fluorine Engagement Ratio (FER) per tier",
        "              Quality zones: Green >=0.75 (high), Yellow 0.50-0.75 (moderate), Red <0.50 (low)",
        "  Right axis: Median FER trend line per tier (blue diamonds with 95% CI band)",
        "  Annotations: med=X, n=Y badge above each box",
        "  Look for: Do higher tiers engage a greater fraction of ligand fluorine atoms?",
        "            A rising trend line from Poor to Perfect confirms tier-quality tracks F-engagement.",
        "",
        "=" * 80,
        "PART 6 — BINDING ENERGETICS",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 15 — Figure_15_Binding_Energetics.png",
        "  Title   : Binding Energetics — Binding Probability + Product Inhibition Penalty  [2-panel]",
        "  Panel A : Violin plot of Binding Probability Score per tier",
        "            Y-axis intentionally zoomed (all values >0.97) — score is essentially a pass/fail threshold",
        "            Median annotated inside each violin; narrow violins = genuine tier similarity",
        "  Panel B : Connected dot-line of Product Inhibition Penalty Score per tier (right axis)",
        "            Higher penalty = greater risk of product feedback inhibition after C-F cleavage",
        "            Shaded band = 95% CI around median",
        "  Look for: Panel A: Poor tier may show a marginally broader/lower distribution.",
        "            Panel B: Perfect/Best tiers may show higher penalty — tighter binding means product",
        "            also binds more tightly, increasing feedback risk.",
        "",
        "-" * 80,
        "Figure 16 — Figure_16_Product_Inhibition_by_Tier.png",
        "  Title   : Product Inhibition Penalty by Tier",
        "  Type    : Box + strip plot with quality zone backgrounds + directional arrows",
        "  X-axis  : Catalytic tier (ordered best to worst)",
        "  Y-axis  : product_inhibition_penalty score (higher = greater feedback inhibition risk)",
        "  Zones   : Green (0-5 low), Yellow (5-15 moderate), Red (>15 high)",
        "  Annotations: Median badge per box; directional arrows on first/last tier boxes",
        "  Stats   : Kruskal-Wallis H-test p-value annotated",
        "  Look for: Top-tier complexes with high penalties warrant attention — cleaved product",
        "            may re-inhibit the enzyme, reducing catalytic turnover.",
        "  Coverage note: ⚠  product_inhibition_penalty is computed only when interaction_density > 1.5.",
        "            Omission is non-random (low-density complexes excluded). Per-box labels show",
        "            'complexes with data / tier total (coverage %)'. Interpret tier differences",
        "            with caution — missing data are structurally biased, not random.",
        "",
        "=" * 80,
        "PART 7 — CHEMICAL SPACE",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 17a — Figure_17a_Chemical_Space_Map.png",
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
        "Figure 17b — Figure_17b_PA_Chemical_Space_Landscape.png",
        "  Title   : Perfect_A Chemical Space Landscape (KDE density + structure thumbnails)",
        "  Type    : KDE density contour overlay on UMAP scatter; inset PyMOL structure thumbnails",
        "  Axes    : UMAP dimensions 1 & 2",
        "  Stars   : Gold/green star markers indicate Perfect_A PFAS positions",
        "  Insets  : PyMOL-rendered protein-ligand structures for each Perfect_A representative",
        "  Look for: Density ridgelines isolating Perfect_A from lower-tier complexes.",
        "            Structure thumbnails reveal active-site geometry at a glance.",
        "            If Perfect_A forms a tight cluster, they share structural/chemical features.",
        "",
        "=" * 80,
        "PART 8 — MULTI-METRIC SYNTHESIS",
        "=" * 80,
        "",
        "-" * 80,
        "Figure 18a — Figure_18a_Fingerprint_TopHits.png",
        "  Title   : Candidate Fingerprint — Top-5 Hits vs. Poor-Tier Baseline  [Radar]",
        "  Type    : Radar / spider chart with normalised metric spokes (0-1 each)",
        "  Spokes  : 5 normalised metrics: AI Confidence, ipTM, Interaction Density, mean pLDDT, Binding Prob.",
        "  Lines   : Top-5 best-tier hits (coloured lines) vs. Poor-tier average baseline (red dashed)",
        "  Look for: Hits forming large polygons outperform on multiple axes simultaneously.",
        "            Spokes where hits touch 1.0 = this metric is at its best possible value.",
        "            Spokes where baseline and hits overlap = metric does not discriminate.",
        "",
        "-" * 80,
        "Figure 18b — Figure_18b_Fingerprint_TierReps.png",
        "  Title   : Candidate Fingerprint — One Representative per Tier  [Radar]",
        "  Type    : Radar / spider chart; one line per tier (best Boltz confidence per tier selected)",
        "  Spokes  : Same 5 normalised metrics as Figure 18a",
        "  Look for: Which tier consistently achieves the largest polygon (best on all metrics)?",
        "            Which metrics most differentiate tiers (largest gap between lines)?",
        "",
        "-" * 80,
        "Figure 19a — Figure_19a_Tier_Success_Rates.png",
        "  Title   : Tier success rates — substrate vs inhibitor geometry",
        "  Type    : Grouped horizontal bar chart (one tier per row)",
        "  X-axis  : Percentage of complexes (%)",
        "  Bars    : Substrate (SN2≥165° + Conf≥0.75); SN2≥165° (any conf); Potential inhibitor (SN2<145°)",
        "  Look for: Top tiers carry higher substrate-geometry pass rates and lower inhibitor fractions.",
        "",
        "Figure 19b — Figure_19b_Conf_SN2_Landscape.png",
        "  Title   : Tier medians in Confidence × SN2-angle space",
        "  Type    : Scatter of per-tier median diamonds (IQR error bars) over a hexbin density background",
        "  X-axis  : Boltz Model Confidence",
        "  Y-axis  : SN2 Attack Angle (°)",
        "  Zones   : Substrate zone (SN2≥165°, green); Inhibitor zone (SN2<145°, orange); Conf threshold 0.75",
        "  Diamonds: Per-tier median position, coloured by tier, annotated with median SN2 + confidence.",
        "  Look for: Top tiers sit inside the substrate zone, right of the 0.75 confidence threshold.",
        "",
        "-" * 80,
        "Figure 20 — Figure_20_Conflict_Composition.png",
        "  Title   : Conflict Category x Tier Composition — Stacked Bar",
        "  Type    : Side-by-side stacked bar (absolute counts left, % composition right)",
        "  X-axis  : Conflict category (Consensus High, Hidden Gem, Ambiguous, Consensus Low, Decoy)",
        "  Y-axis left : Absolute number of complexes in each category x tier combination",
        "  Y-axis right: Percentage composition of each conflict category by tier (sums to 100%)",
        "  Colour  : Degrader tier (same palette as all other figures)",
        "  Definitions:",
        "    Consensus High  = top 50% on both physics AND AI confidence",
        "    Decoy           = high AI confidence but weak physics geometry (false positive risk)",
        "    Hidden Gem      = strong physics but low AI confidence (under-estimated by AI)",
        "    Ambiguous       = intermediate on both axes",
        "    Consensus Low   = bottom 50% on both axes",
        "  Look for: Consensus High dominated by Perfect/Best tiers = metrics agree.",
        "            Decoys concentrated in Poor/Good = AI over-confident on weaker candidates.",
        "            The % panel shows tier composition normalised — compare categories fairly.",
        "",
        "-" * 80,
        "Figure 21 — Figure_21_Hidden_Gems_DeepDive.png",
        "  Title   : Hidden Gems — Physics-Good / AI-Missed Conflicts",
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
        "Figure 22 — Figure_22_Category_Overlap_Euler.png",
        "  Title   : Category Overlap — Euler / Venn Diagram",
        "  Type    : Manual circle patches (3 main circles + Perfect_A dashed overlay)",
        "  Sets    :",
        "    A (Physics-Strong)  : Pareto_Rank <= 33rd percentile (top 33% physics geometry)",
        "    B (AI-Strong)       : Boltz_Model_Confidence / ipTM / pTM >= 67th percentile (top 33% AI)",
        "    C (Tier-Strong)     : degrader_tier in {Perfect_A, Perfect_B, Best_A}",
        "  Overlays:",
        "    Dashed pink circle  : Perfect_A — weighted to their position across A/B/C",
        "    Dotted grey circle  : All-three intersection highlight",
        "    Star markers        : 3R3U and DeHa4 reference controls (positioned in correct region)",
        "  Annotations: 7 region count boxes (count + row %); 4-column legend",
        "  Look for: How many top-tier candidates excel across all three metrics simultaneously?",
        "            The A+B+C intersection = Consensus Best on physics, AI, AND catalytic tier.",
        "            Reference controls (3R3U, DeHa4) should land in the A+B+C region.",
        "",
        "-" * 80,
        "Figure 23 — Figure_23_Top25_Multitarget_Proteins.png",
        "  Type    : Horizontal stacked bar chart (one bar per protein, stacked by tier)",
        "  X-axis  : Number of unique PFAS ligands degraded (best tier per protein–ligand pair)",
        "  Y-axis  : Protein (top 25, ranked by quality-weighted degradation breadth)",
        "  Colours : 6-tier palette (Perfect_A green → Decoy grey)",
        "  Look for: Wide bars = broad-spectrum PFAS degraders. Deep-green left stack = top-tier activity across many ligands.",
        "",
        "-" * 80,
        "Figure 23b — Figure_23b_Ligand_Network.png",
        "  Title   : PFAS Ligand Network — Tier-Coloured Protein Interactions",
        "  Type    : Radial spoke network; 27 ligands on periphery, tier-coloured flows from centre",
        "  Spokes  : Width proportional to number of proteins at that tier for the ligand",
        "  Colours : Same 6-tier palette; flow colour = tier of the protein–ligand pairing",
        "  Labels  : PFAS short name + total protein count outside each spoke cluster",
        "  Look for: Short-chain / polar PFAS (TFA, Fluoroacetate) have densest green spokes.",
        "            Long-chain PFAS show fewer Perfect_A hits — used to guide substrate scope claims.",
        "",
        "-" * 80,
        "Figure 24 — Figure_24_Sankey_Workflow.png",
        "  Title   : Sankey Pathway — full FAcD tier-decision rules → Final Tier",
        "  Type    : Flat Sankey (ribbon polygons); node boxes labelled with COUNTS",
        "  Columns : Nuc Distance (CFG cuts) · Carboxylate Clamp (ARG111/114) ·",
        "            Halide Stabilisation · Triad Geometry (Nuc–Base · Base–Acid) ·",
        "            SN2 Angle · Mechanistic Fingerprint · Final Tier",
        "  Order   : best category on top in every column (Perfect_A, ≤3.0 Å, ≥175°, ≥0.9)",
        "  Bins    : all cut points sourced from CFG (TIER_NUC_DIST / TIER_ANGLE_MIN /",
        "            TIER_NB_MAX / TIER_BA_MAX / TIER_MECH_MIN) — match 02_Production gates",
        "  Ribbons : Width proportional to complex count; coloured by destination tier",
        "  Look for: Perfect_A requires ALL rules to pass (tight nuc + clamp + stabilisation",
        "            + tight triad + ≥175° + mech ≥0.9); the narrowing chain shows the attrition.",
        "",
        "-" * 80,
        "Figure 25a — Figure_25a_PFAS_Size_Hexbin_Landscape.png",
        "  Title   : PFAS Chain-Length vs SN2 Geometry — Substrate Preference and Inhibition Risk",
        "  Type    : Hexbin density (F-count × SN2 angle) with tier scatter overlay and rolling median",
        "  X-axis  : Total fluorine count (proxy for carbon chain length); ticks every 2 units",
        "  Y-axis  : SN2 Attack Angle (°); green zone ≥165° = substrate; red zone <145° = inhibitor",
        "  Look for: Rolling median SN2 trend across chain lengths; tier scatter reveals Perfect_A",
        "            distribution relative to substrate zone.",
        "",
        "-" * 80,
        "Figure 25b — Figure_25b_PFAS_Size_Composition_Merged.png",
        "  Title   : Chain-length Composition — Mechanistic Outcome + Degrader Tier (merged)",
        "  Type    : Single panel; per chain-length bin a LEFT stacked bar (outcome) and",
        "            a RIGHT stacked bar (degrader tier) side by side.",
        "  Left bar: Substrate (green), Borderline (amber), Reactive low conf (blue),",
        "            Non-reactive (grey), Potential Inhibitor (orange-red)",
        "  Right bar: Degrader-tier composition (Perfect_A → Decoy palette)",
        "  Look for: Substrate % peaks at short/medium chains; inhibition risk rises with chain length.",
        "",
        "-" * 80,
        "Figure 25c — Figure_25c_PFAS_Size_Confidence_Boxplot.png",
        "  Title   : AI Binding Confidence by Chain Length",
        "  Type    : Violin plot of Boltz Model Confidence per PFAS size bin with median trend line",
        "  Look for: Confidence remains high even for geometrically unfeasible long-chain complexes",
        "            — the core AI/physics disconnect for PFAS remediation.",
        "",
        "=" * 80,
    ]
    desc_path = out_dir / "05_Figure_Descriptions.txt"
    desc_path.write_text("\n".join(lines), encoding='utf-8')
    return desc_path


# ===============================================================================
# SECTION 6: MAIN EXECUTION
# ===============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Boltz-2 Master Validation & Phylogeny Framework",
        usage="%(prog)s <run_folder>  (e.g. Boltz-2_Run_20260309T085406Z)"
    )
    parser.add_argument("run", help="Name of the Boltz-2 run folder (e.g. Boltz-2_Run_20260309T085406Z)")
    args = parser.parse_args()

    root_dir = Path.cwd()
    run_path = root_dir / args.run
    if not run_path.exists():
        print(f"Error: run folder not found: {run_path}")
        sys.exit(1)
    
    _utils_mod.print_script_banner(
        "03_Validation_Figures_FAcDs.py",
        "Multi-Objective Ranking  ·  Pareto Frontiers  ·  Publication Figures",
    )
    print(f"  Run Name : {run_path.name}", flush=True)

    prod_dir = run_path / "1_Boltz2_Production"
    out_dir = run_path / "3_Validation_Figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    np.random.seed(42)

    reporter = ReportManager(out_dir)

    try:
        # -------------------------------------------------------------------------------
        # Pipeline Execution Phase
        # -------------------------------------------------------------------------------
        df, features = load_and_prep_data(prod_dir, reporter)
        df = perform_advanced_ranking(df, features, out_dir, reporter)
        df = analyse_conflicts(df, out_dir, reporter)

        final_csv = out_dir / "03_Final_Validated_Master.csv"
        df.to_csv(final_csv, index=False)
        reporter.log(f"Final Validated Dataset Saved: {final_csv.resolve()}")

        generate_comprehensive_figures(df, features, out_dir, reporter)

        _thumb_dir_cleanup = out_dir / "_pa_thumbnails"
        if _thumb_dir_cleanup.exists():
            shutil.rmtree(_thumb_dir_cleanup)

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
