"""Read-only loaders for the FAcDs data that feed the 20 figures.

Nothing here writes outside TEST/. The big source CSVs live under the run directory and are only
READ; small derived tables are cached under TEST/data_cache/ so repeated runs stay fast. Every
figure draws from real project data — the ranked candidate table, the per-model QC variance table,
the PCA loadings and the statistical-test summary — never from simulated values.
"""
from __future__ import annotations
import functools
from pathlib import Path
import numpy as np
import pandas as pd

_ROOT = Path("/mnt/wdpassport/FAcDs/02_Data/Boltz-2_Run_20260309T085406Z")
RANKED = _ROOT / "1_Boltz2_Production" / "6_Boltz2_FAcDs_Ranked_20260716_225324.csv"
ANALYSIS = _ROOT / "3_Validation_Figures" / "01_Analysis_Data"
VARIANCE = ANALYSIS / "07_Boltz2_MultiModel_QC_Variance.csv"
PCA_LOAD = ANALYSIS / "02_PCA_Loadings.csv"
STATS = ANALYSIS / "06_Statistical_Tests.csv"

TIER_ORDER = ["Tier_1A", "Tier_1B", "Tier_2A", "Tier_2B", "Tier_3", "Tier_4", "Tier_5_Decoy"]

# Short ligand labels — most PFAS names are already their standard abbreviation; only the two
# short-chain acids spell out, so map just those (Fluoroacetate → FA, Difluoroacetate → DFA).
LIGAND_SHORT = {"Fluoroacetate": "FA", "Difluoroacetate": "DFA"}

# Compact display labels for the metric columns, so heatmap/axis ticks stay inside their panel.
METRIC_LABELS = {
    "iptm": "iPTM", "ptm": "pTM", "mean_plddt": "pLDDT", "active_site_plddt": "AS pLDDT",
    "Boltz_Model_Confidence": "Confidence", "Binding_Probability_Score": "Binding P",
    "custom_affinity_score": "Affinity", "interaction_density": "Int. density",
    "competence_score": "Competence", "feasibility_factor": "Feasibility",
    "mechanistic_score_effective": "Mech. score", "mechanistic_score": "Mech. score (raw)",
    "soft_catalytic_score": "Soft cat.", "catalytic_constellation_score": "Cat. constel.",
    "halide_stabilisation_score": "Halide stab.", "carboxylate_clamp_integrity": "Cbx clamp",
    "pocket_occupancy": "Pocket occ.", "num_interactions": "N interact.",
    "count_hydrogen_bond": "H-bonds", "count_salt_bridge": "Salt bridges",
    "count_hydrophobic": "Hydrophobic", "count_pi_stacking": "π-stack",
    "SN2_Attack_Angle": "S$_N$2 angle", "scissile_cf_bde": "C–F BDE",
    "sn2_backside_occlusion": "Backside occ.", "Interaction_Density_Norm": "Int. density",
    "Chemical_Affinity_Score": "Affinity", "Binding_Probability": "Binding P",
}


def mlabel(col: str) -> str:
    """Short display label for a metric column (falls back to a tidied raw name)."""
    return METRIC_LABELS.get(col, col.replace("_", " "))

# The subset of ranked columns the figures actually use — keeps 58k×138 down to a light frame.
_USECOLS = [
    "Scientific_Rank", "job_name", "Protein_Name", "Ligand_Name", "degrader_tier", "is_degrader",
    "Boltz_Model_Confidence", "iptm", "ptm", "mean_plddt", "Binding_Probability_Score",
    "custom_affinity_score", "interaction_density", "SN2_Attack_Angle", "scissile_cf_bde",
    "sn2_backside_occlusion", "beta_f_count", "feasibility_factor", "competence_score",
    "mechanistic_score", "mechanistic_score_effective", "soft_catalytic_score",
    "catalytic_constellation_score", "active_site_plddt", "num_interactions",
    "count_hydrogen_bond", "count_hydrophobic", "count_salt_bridge", "count_fluorine_contact",
    "count_fluorine_polar", "count_fluorous_hydrophobic", "count_pi_stacking",
    "Dist_Nucleophile", "Dist_Base", "Dist_Acid", "Dist_Clamp1", "Dist_Clamp2",
    "halide_stabilisation_score", "carboxylate_clamp_integrity", "pocket_occupancy", "fit_ratio",
    "ligand_volume", "active_site_volume", "SN2_Trajectory_Deviation_A", "MD_Selected", "MD_Rank",
]


@functools.lru_cache(maxsize=1)
def ranked() -> pd.DataFrame:
    """The per-candidate ranked table (58k rows), light column subset, tier as an ordered category."""
    cols = pd.read_csv(RANKED, nrows=0).columns
    use = [c for c in _USECOLS if c in cols]
    df = pd.read_csv(RANKED, usecols=use, low_memory=False)
    if "degrader_tier" in df:
        df["degrader_tier"] = pd.Categorical(df["degrader_tier"], categories=TIER_ORDER, ordered=True)
    df["ligand_short"] = (df["Ligand_Name"].astype(str).str.replace(r"^\d+_", "", regex=True)
                          .replace(LIGAND_SHORT))
    return df


@functools.lru_cache(maxsize=1)
def variance() -> pd.DataFrame:
    """Per-model QC replicates (ptm/iptm/ligand_iptm/confidence/sn2 geometry) — the spread source."""
    df = pd.read_csv(VARIANCE, low_memory=False)
    return df


@functools.lru_cache(maxsize=1)
def pca_loadings() -> pd.DataFrame:
    df = pd.read_csv(PCA_LOAD, index_col=0)
    df.index.name = "feature"
    return df


@functools.lru_cache(maxsize=1)
def stats() -> pd.DataFrame:
    return pd.read_csv(STATS, low_memory=False)


# ---- shared derived views -------------------------------------------------------------------
NUMERIC_METRICS = [
    "iptm", "ptm", "mean_plddt", "Boltz_Model_Confidence", "Binding_Probability_Score",
    "custom_affinity_score", "interaction_density", "competence_score", "feasibility_factor",
    "mechanistic_score_effective", "soft_catalytic_score", "catalytic_constellation_score",
    "halide_stabilisation_score", "carboxylate_clamp_integrity", "pocket_occupancy",
    "num_interactions", "count_hydrogen_bond", "count_salt_bridge", "SN2_Attack_Angle",
    "scissile_cf_bde", "sn2_backside_occlusion",
]


def metric_matrix() -> pd.DataFrame:
    """Numeric metric columns only (for correlation / PCA-style figures), NaNs dropped column-wise."""
    df = ranked()
    cols = [c for c in NUMERIC_METRICS if c in df.columns]
    m = df[cols].apply(pd.to_numeric, errors="coerce")
    return m.loc[:, m.notna().sum() > 100]


def per_ligand(metric: str) -> pd.DataFrame:
    """mean / sd / n / se of `metric` per ligand, with the modal tier and degrader flag attached."""
    df = ranked()
    g = df.groupby("ligand_short", observed=True)[metric].agg(["mean", "std", "count"])
    g["se"] = g["std"] / np.sqrt(g["count"].clip(lower=1))
    tier = df.groupby("ligand_short", observed=True)["degrader_tier"].agg(
        lambda s: s.mode().iloc[0] if len(s.mode()) else np.nan)
    deg = df.groupby("ligand_short", observed=True)["is_degrader"].mean() >= 0.5
    g["tier"] = tier
    g["is_degrader"] = deg
    return g.dropna(subset=["mean"]).sort_values("mean", ascending=False)
