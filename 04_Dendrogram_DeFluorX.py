#!/usr/bin/env python3

"""
===============================================================================
DeFluorX Pipeline  |  Step 04  |  Phylogenetic Analysis & Deployment
===============================================================================
Computes UPGMA dendrograms from alignment-free K-mer (k=3) cosine sequence
distances and deploys interactive D3.js HTML visualisations with per-tier
isolations.

Scientific scope - read before citing
─────────────────────────────────────
This is a fast sequence-similarity dendrogram for clustering and visualisation,
NOT a substitution-model dendrogram. K-mer cosine + UPGMA assumes a constant
evolutionary rate (ultrametricity) and applies no substitution model, no
indel/back-mutation handling, and no branch-support estimation. Treat the tree as
a similarity map of the screened candidates - do not draw evolutionary-rate or
ancestry conclusions from branch lengths. For publication-grade phylogenetics,
build an MSA (MAFFT / Clustal-Ω) and a maximum-likelihood or Bayesian tree
(IQ-TREE / RAxML / MrBayes) with bootstrap or posterior support, then overlay the
tier classification produced here.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 20 August 2026 <─────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 04_Dendrogram_DeFluorX.py
  Role          : Phylogenetic analysis and interactive tree visualisation.
  Imports from  : 00_02_Project_Utils_DeFluorX.py  (console_info / console_separator)
  Reads         : <Run>/3_Validation_Figures/01_Analysis_Data/03_Figure_Enriched_Dataset.csv
                  <Run>/1_Boltz2_Production/1_Input_Data/*.fasta
  Writes        : <Run>/4_Dendrogram/01_Global_Master_Dendrogram.tree
                  <Run>/4_Dendrogram/02_Global_Master_Matrix_Data.csv
                  <Run>/4_Dendrogram/03_Global_Master_Interactive_App.html
                  <Run>/4_Dendrogram/04_Global_Master_Delivery_Suite.zip
                  <Run>/4_Dendrogram/05_Tiers/<Tier>_*  (per-tier tree + HTML)
                  <Run>/4_Dendrogram/00_Dendrogram.log
  Upstream      : 03_Validation_Figures_DeFluorX.py → writes 03_Figure_Enriched_Dataset.csv
  Downstream    : None (terminal analysis step)
───────────────────────────────────────────────────────────────────────────────

── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Matrix Scaling: Distance matrix calculation scales O(N^2) in both time and
     memory; datasets exceeding 5k proteins may require HPC nodes with 64GB+ RAM.
  2. K-mer Sensitivity: Tree topology is derived from k-mer (k=3) frequencies;
     may differ from traditional MSA-based maximum-likelihood trees.
  3. Interactive Delivery: HTML apps rely on D3.js (loaded via CDN); local
     viewing requires an active internet connection or offline D3 assets.
  4. Zip Bottleneck: Creating the delivery suite for runs with thousands of
     small files can be slow on network-mounted filesystems.
───────────────────────────────────────────────────────────────────────────────

Usage:
    conda activate PFAS
    python 04_Dendrogram_DeFluorX.py Boltz-2_Run_20260309T085406Z
───────────────────────────────────────────────────────────────────────────────

-------------------------------------------------------------------------------
Scientific References:
    1. Hierarchical clustering & pairwise distances (SciPy):
       - Virtanen, P. et al. (2020) SciPy 1.0. Nature Methods 17:261–272.
       - DOI: https://doi.org/10.1038/s41592-019-0686-2
    2. UPGMA agglomerative clustering:
       - Sokal, R.R. & Michener, C.D. (1958) A statistical method for evaluating
         systematic relationships. Univ Kansas Sci Bull 38:1409–1438.
    3. Sequence parsing (Biopython SeqIO):
       - Cock, P.J.A. et al. (2009) Biopython. Bioinformatics 25:1422–1423.
       - DOI: https://doi.org/10.1093/bioinformatics/btp163
    4. Interactive visualisation (D3.js):
       - Bostock, M., Ogievetsky, V. & Heer, J. (2011) D3: Data-Driven Documents.
         IEEE Trans Vis Comput Graph 17:2301–2309. DOI: https://doi.org/10.1109/TVCG.2011.185
    5. Numerics: Harris, C.R. et al. (2020) Array programming with NumPy. Nature 585:357–362.
       - DOI: https://doi.org/10.1038/s41586-020-2649-2
-------------------------------------------------------------------------------
"""

# =============================================================================
# SECTION 1: IMPORTS & CONFIGURATION
# =============================================================================

# -----------------------------------------------------------------------------
# Step 1.1: Standard Library Imports
# -----------------------------------------------------------------------------
import sys
import argparse
import json
import zipfile
import re
import traceback
from pathlib import Path
from collections import Counter

# -----------------------------------------------------------------------------
# Step 1.2: Scientific Stack Imports
# -----------------------------------------------------------------------------
'''
CPU usage cap (total cores − 2; mirrors CFG.PREP_CPU_RESERVE). Reserve 2 cores
for OS/desktop stability by limiting the thread-pool maths libraries (BLAS /
MKL / OpenMP / NumExpr). Must precede numpy/scipy import to take effect;
setdefault() preserves any value exported by the caller or pipeline runner.
'''
import os as _os
_CPU_CAP = str(max(1, (_os.cpu_count() or 4) - 2))
for _tv in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_tv, _CPU_CAP)

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, to_tree
from Bio import SeqIO

# -----------------------------------------------------------------------------
# Step 1.3: Pipeline modules (00_01) via importlib
# -----------------------------------------------------------------------------
import importlib.util as _ilu
# Consolidated top-level imports; any optional/heavy dependency stays local to its caller.
import time as _time

def _load_module(name: str, path: Path):
    if not path.exists():
        raise FileNotFoundError(f"Required module not found: {path}")
    spec = _ilu.spec_from_file_location(name, str(path))
    mod  = _ilu.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

_utils_mod      = _load_module("ProjectUtils", Path(__file__).resolve().parent / "00_02_Project_Utils_DeFluorX.py")
_cfg_mod        = _load_module("ProjectConfig", Path(__file__).resolve().parent / "00_01_Project_Config_DeFluorX.py")
CFG             = _cfg_mod.CFG()
_console_info   = _utils_mod.console_info
_console_sep    = _utils_mod.console_separator
_setup_logging  = _utils_mod.setup_logging

# -----------------------------------------------------------------------------
# Step 1.4: Global Configuration
# -----------------------------------------------------------------------------
SEPARATOR = "-" * 80

# Fallback column mapping for compatibility across slightly different CSV schemas
COL_MAP = CFG.VIS_PHYLO_COLUMN_MAP


# =============================================================================
# SECTION 2: LOGGING & UTILITIES
# =============================================================================

logger = None  # Initialised in main()

def console_info(msg: str) -> None:
    _console_info(msg, logger)

def console_separator() -> None:
    _console_sep(logger, heavy=True)


ReportManager = _utils_mod.ReportManager   # shared logger (00_02)


def _make_reporter(out_dir: Path):
    """Construct the shared ReportManager with this step's log path/header/logger."""
    return ReportManager(
        out_dir / "00_Dendrogram.log",
        "BOLTZ-2 PHYLOGENY PIPELINE REPORT",
        # Non-logging log_fn: ReportManager.log already appends every line to the log file itself, so a
        # logger-bound log_fn (console_info → logger → same file) would write each line twice. print keeps
        # it on the console; the [LOG] append keeps it in the file, once.
        separator=SEPARATOR, rule_width=79, log_fn=print)


def clean_id(name: str) -> str:
    """Normalises protein IDs for fuzzy matching between FASTA headers and CSV rows."""
    if not isinstance(name, str): return "UNKNOWN"
    return re.sub(r"[^a-zA-Z0-9]", "", name).lower()


# =============================================================================
# SECTION 3: BIO-MATHEMATICS
# =============================================================================

# -----------------------------------------------------------------------------
# Step 3.1: K-mer Profiling
# -----------------------------------------------------------------------------

def get_kmer_counts(seq: str, k: int = CFG.DENDRO_KMER_SIZE) -> dict:
    """Generates K-mer frequency profile for a given amino acid sequence."""
    seq = seq.upper()
    return dict(Counter(seq[i:i+k] for i in range(len(seq) - k + 1)))


'''
NB: pairwise cosine distances are computed in bulk via scipy.spatial.distance.pdist
(metric='cosine') inside generate_upgma_newick - one vectorised BLAS call, so no
scalar per-pair distance helper is needed.
'''


# -----------------------------------------------------------------------------
# Step 3.2: UPGMA Clustering
# -----------------------------------------------------------------------------

def generate_upgma_newick(sequences: dict) -> tuple[str, list]:
    """Computes pairwise cosine distance matrix and performs UPGMA clustering to Newick."""
    labels = list(sequences.keys())
    n = len(labels)

    if n == 0: return "();", []
    if n == 1: return f"({labels[0]}:0.0);", labels

    profiles = [get_kmer_counts(sequences[l]) for l in labels]

    '''
    Build dense k-mer matrix and compute all pairwise cosine distances in one
    vectorised BLAS call (pdist) instead of an O(N²) Python loop.
    '''
    all_kmers      = sorted(set(k for p in profiles for k in p))
    kmer_mat       = np.array([[p.get(k, 0) for k in all_kmers] for p in profiles],
                               dtype=np.float64)
    condensed_dist = pdist(kmer_mat, metric=CFG.DENDRO_DISTANCE_METRIC)
    condensed_dist = np.nan_to_num(condensed_dist, nan=1.0)  # safety: zero-vector → max dist
    Z              = linkage(condensed_dist, method=CFG.DENDRO_LINKAGE_METHOD)
    tree_node      = to_tree(Z, rd=False)

    def _nwk_safe(lbl) -> str:
        # Newick metacharacters ( ) , : ; [ ] ' " and whitespace corrupt the tree
        # string and break downstream parsers; collapse any run of them to '_'.
        return re.sub(r"[\s(),:;\[\]'\"]+", "_", str(lbl))

    def build_newick(node, parentdist):
        '''
        Clamp branch length at 0: a non-monotonic linkage can give node.dist >
        parentdist, which would emit a negative branch length (rejected by most
        tree parsers). max(0.0, …) keeps the Newick valid.
        UPGMA is ultrametric: node height = half the cophenetic (merge) distance
        scipy stores in node.dist, so each branch length is (parent−node)/2. Without
        the /2 every leaf-to-leaf path is exactly 2× the true cosine distance.
        '''
        if node.is_leaf():
            return f"{_nwk_safe(labels[node.id])}:{max(0.0, (parentdist - node.dist) / 2.0):.4f}"
        left_str  = build_newick(node.left,  node.dist)
        right_str = build_newick(node.right, node.dist)
        return f"({left_str},{right_str}):{max(0.0, (parentdist - node.dist) / 2.0):.4f}"

    left      = build_newick(tree_node.left,  tree_node.dist)
    right     = build_newick(tree_node.right, tree_node.dist)
    newick_str = f"({left},{right});"

    return newick_str, labels


# =============================================================================
# SECTION 4: PHYLOGENY DEPLOYMENT
# =============================================================================

# -----------------------------------------------------------------------------
# Step 4.1: Suite Packaging
# -----------------------------------------------------------------------------

def package_deployment(out_dir: Path, prefix: str, nwk_str: str,
                        csv_df: pd.DataFrame, labels: list,
                        reporter: ReportManager):
    """Writes Newick, CSV, interactive HTML app, and bundles them into a ZIP delivery suite."""
    out_dir.mkdir(parents=True, exist_ok=True)

    nwk_path  = out_dir / f"01_{prefix}_Dendrogram.tree"
    csv_path  = out_dir / f"02_{prefix}_Matrix_Data.csv"
    html_path = out_dir / f"03_{prefix}_Interactive_App.html"
    zip_path  = out_dir / f"04_{prefix}_Delivery_Suite.zip"

    nwk_path.write_text(nwk_str)

    csv_str = csv_df.to_csv(index=False)
    csv_path.write_text(csv_str)

    '''
    Deduplicate to one row per protein for the embedded HTML state.
    The full CSV (all ligand×protein rows) can be multi-MB and slow browsers.
    '''
    prot_col = next((c for c in CFG.VIS_PHYLO_COLUMN_MAP["Protein_Name"] if c in csv_df.columns), None)
    if prot_col and len(csv_df) > len(labels):
        embed_df = csv_df.drop_duplicates(subset=[prot_col], keep="first")
    else:
        embed_df = csv_df
    state_json   = json.dumps({"tree": nwk_str, "csv": embed_df.to_csv(index=False)})
    html_content = HTML_APP_TEMPLATE.replace(
        "<!-- INJECT_EMBEDDED_STATE -->",
        f"<script>window.EMBEDDED_STATE = {state_json};</script>"
    ).replace(
        "{{ JSON_BASE_COLORS }}", json.dumps(CFG.TIER_COLOUR)
    ).replace(
        "{{ JSON_TIER_ORDER }}",  json.dumps(CFG.TIER_ORDER)
    ).replace(
        "{{ JSON_TIER_ALPHAS }}", json.dumps(CFG.VIS_TIER_ALPHAS)
    ).replace(
        "{{ JSON_LIGAND_COLORS }}", json.dumps(list(CFG.VIS_LIGAND_SERIES))
    ).replace(
        "{{ JSON_LIGAND_SHORT }}",  json.dumps(CFG.VIS_LIGAND_SHORT)
    ).replace(
        "{{ JSON_CLADE_COLORS }}",  json.dumps(list(CFG.VIS_CLADE_SERIES))
    ).replace(
        '<div id="upload-section">',
        '<div id="upload-section" style="display:none;">'
    )
    html_path.write_text(html_content, encoding="utf-8")

    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.write(nwk_path,  nwk_path.name)
        zf.write(csv_path,  csv_path.name)
        zf.write(html_path, html_path.name)

    _lig_col = next((c for c in CFG.VIS_PHYLO_COLUMN_MAP["Ligand_Name"] if c in csv_df.columns), None)
    _n_lig = int(csv_df[_lig_col].nunique()) if _lig_col else 0
    _lig_txt = f" × {_n_lig} ligands" if _n_lig else ""
    reporter.log(f"  ✔ Suite Generated: {prefix} ({len(labels)} proteins{_lig_txt}) -> {out_dir.resolve()}")


# -----------------------------------------------------------------------------
# Step 4.2: Dendrogram Engine
# -----------------------------------------------------------------------------

def generate_phylogenies(df: pd.DataFrame, prod_dir: Path,
                          out_dir: Path, reporter: ReportManager):
    """Drives global and per-tier phylogenetic tree generation."""
    with open(reporter.path, "a") as f: f.write("\n[LOG] Dendrogram Engine Started\n")
    console_separator()

    input_dir = next(
        (prod_dir / n for n in ["1_Input_Data", "1_Input_FASTA_and_SMILES"]
         if (prod_dir / n).exists()),
        prod_dir / "1_Input_Data"
    )
    tier_dir = out_dir / "05_Tiers"
    out_dir.mkdir(parents=True, exist_ok=True)

    prot_col = next((c for c in COL_MAP["Protein_Name"] if c in df.columns), None)
    tier_col = next((c for c in COL_MAP["Tier"]         if c in df.columns), None)

    if not prot_col:
        reporter.log("  ! Critical Error: Could not identify Protein column in CSV. Skipping.")
        return

    # Sorted: max() returns the FIRST maximal element, so two FASTAs with the same record
    # count would otherwise be separated by filesystem order.
    fasta_files = sorted(input_dir.glob("*.fasta")) + sorted(input_dir.glob("*.fa"))
    if not fasta_files:
        reporter.log("  ! Critical Error: No FASTA file found in input directory. Skipping.")
        return

    '''
    The input directory can hold several FASTAs (e.g. the merged panel + a
    single-sequence reference such as DeHa4_Ref.fasta). Selecting the first
    glob hit is order-dependent and may pick the 1-sequence reference, leaving
    nothing to match the CSV. Choose the full panel: a name flagged as the
    merged input wins, otherwise the file with the most records.
    '''
    def _record_count(fp):
        try:
            with open(fp) as fh:
                return sum(1 for ln in fh if ln.startswith(">"))
        except OSError:
            return 0
    merged_hint = [f for f in fasta_files if any(t in f.name.lower() for t in ("merged", "_inp_", "inp_"))]
    fasta_path = max(merged_hint or fasta_files, key=_record_count)
    if _record_count(fasta_path) == 0:      # every hinted file is empty → fall back to any non-empty FASTA
        fasta_path = max(fasta_files, key=_record_count)

    fasta_dict = {}
    for r in SeqIO.parse(str(fasta_path), "fasta"):
        raw_id    = r.id.split("|")[0].strip()
        clean_key = re.sub(r"^\d+_", "", raw_id, count=1)
        fasta_dict[clean_key] = str(r.seq)

    # -----------------------------------------------------------------------------
    # Phase 1: Global Master Dendrogram
    # -----------------------------------------------------------------------------
    with open(reporter.path, "a") as f: f.write("\n[LOG] Phase 1: Generating Global Master Dendrogram\n")

    valid_csv_prots = set(df[prot_col].astype(str).unique())
    # Pre-compute the cleaned CSV ids once (O(M)) so the membership test below is O(1)
    # per FASTA key rather than an O(N·M) regex recomputation.
    _valid_clean_ids = {clean_id(p) for p in valid_csv_prots}
    global_seqs = {
        k: v for k, v in fasta_dict.items()
        if clean_id(k) in _valid_clean_ids
    }

    '''
    Reporting only (does NOT change the matching): name the CSV proteins with no FASTA sequence, so
    their absence from the tree is logged rather than silent. The controls (3R3U_Control, DeHa4_Control)
    are absent by construction - their sequences are not in the merged FASTA.
    '''
    _fasta_clean = {clean_id(k) for k in fasta_dict}
    _unmatched = sorted(p for p in valid_csv_prots if clean_id(p) not in _fasta_clean)
    if _unmatched:
        reporter.log(f"  ℹ {len(_unmatched)} CSV protein(s) not in the merged FASTA, so absent from the "
                     f"tree: {', '.join(_unmatched[:5])}{' …' if len(_unmatched) > 5 else ''}")

    if not global_seqs:
        reporter.log("  ! Name mismatch between FASTA and CSV - no sequences matched. Skipping global dendrogram.")
        return

    nwk_str, labels = generate_upgma_newick(global_seqs)
    package_deployment(out_dir, "Global_Master", nwk_str, df, labels, reporter)

    # -----------------------------------------------------------------------------
    # Phase 2: Tier-Specific Phylogenies
    # -----------------------------------------------------------------------------
    _TIER_RANK_ORDER = CFG.TIER_ORDER

    if tier_col:
        df = _utils_mod.standardise_dataframe_tiers(df, CFG)
        with open(reporter.path, "a") as f: f.write("\n[LOG] Phase 2: Generating Isolated Tier Phylogenies\n")
        # Every classified tier gets its own phylogeny, the decoy bucket included (it is the largest and
        # its clustering shows what was rejected); only genuinely unclassified rows ("Unknown") are dropped.
        _raw_tiers = [t for t in df[tier_col].dropna().unique() if t != "Unknown"]
        def _tier_key(t):
            try: return _TIER_RANK_ORDER.index(t)
            except ValueError: return len(_TIER_RANK_ORDER)
        _sorted_tiers = sorted(_raw_tiers, key=_tier_key)
        for _tier_idx, tier in enumerate(_sorted_tiers, 1):
            tier_subfolder = f"{_tier_idx:02d}_{tier}"

            tier_df   = df[df[tier_col] == tier].copy()
            if tier_df.empty: continue

            tier_prots = set(tier_df[prot_col].astype(str).unique())
            _tier_clean_ids = {clean_id(p) for p in tier_prots}
            tier_seqs  = {
                k: v for k, v in global_seqs.items()
                if clean_id(k) in _tier_clean_ids
            }

            if len(tier_seqs) < 2:
                reporter.log(f"  - Skipping {tier}: requires at least 2 proteins for a tree.")
                continue

            nwk_str, labels = generate_upgma_newick(tier_seqs)
            package_deployment(tier_dir / tier_subfolder, tier, nwk_str, tier_df, labels, reporter)
    else:
        reporter.log("  ! Skipping Phase 2: No tier column found in CSV.")

    with open(reporter.path, "a") as f: f.write("\n[LOG] ✔ Dendrogram Pipeline Completed Successfully\n")


# =============================================================================
# SECTION 5: MAIN EXECUTION
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Boltz-2 Phylogenetic Analysis Pipeline")
    parser.add_argument("run", nargs="?", help="Name of the Run Folder (e.g., Boltz-2_Run_2026...)")
    args = parser.parse_args()

    root_dir = Path.cwd()

    if args.run:
        run_path = root_dir / args.run
        if not run_path.exists():
            print(f"Error: Specified run folder not found: {run_path}")
            sys.exit(1)
    else:
        runs = sorted(
            [d for d in root_dir.iterdir() if d.is_dir() and "Boltz-2_Run_" in d.name],
            key=lambda x: (x.stat().st_mtime, x.name)   # name breaks an mtime tie deterministically
        )
        if not runs:
            print("Error: No 'Boltz-2_Run_*' folders found in current directory.")
            sys.exit(1)
        run_path = runs[-1]

    _utils_mod.print_script_banner(
        "04_Dendrogram_DeFluorX.py",
        "Phylogenetic Tree Construction  ·  Sequence Clustering  ·  Taxonomic Analysis",
    )
    print(f"  Run Name : {run_path.name}", flush=True)

    prod_dir = run_path / "1_Boltz2_Production"
    val_dir  = run_path / "3_Validation_Figures"
    out_dir  = run_path / "4_Dendrogram"
    out_dir.mkdir(parents=True, exist_ok=True)

    '''
    Build the reporter FIRST: its __init__ opens the log "w" (writes the header), so doing it after
    _setup_logging would truncate the handler's first lines (e.g. "Loaded: <csv>"). Header first, then
    the logging handler appends.
    '''
    reporter = _make_reporter(out_dir)

    global logger
    logger = _setup_logging(out_dir / "00_Dendrogram.log", "04_Dendrogram")

    # Locate the figure-enriched dataset produced by 03_Validation_Figures_DeFluorX.py
    # (written under 01_Analysis_Data; the root and rglob lookups cover a non-default out_dir layout).
    csv_candidates = (
        sorted(val_dir.glob(f"01_Analysis_Data/{CFG.FILE_VALIDATED_MASTER}")) or
        sorted(val_dir.glob(CFG.FILE_VALIDATED_MASTER)) or
        sorted(val_dir.rglob(CFG.FILE_VALIDATED_MASTER))
    )
    if not csv_candidates:
        print(f"Error: No {CFG.FILE_VALIDATED_MASTER} found in {val_dir.resolve()}")
        print("       Run 03_Validation_Figures_DeFluorX.py first to generate it.")
        sys.exit(1)

    csv_path = csv_candidates[0]
    df = pd.read_csv(csv_path, low_memory=False)
    console_info(f"Loaded: {csv_path.name} ({len(df):,} rows)")

    try:
        generate_phylogenies(df, prod_dir, out_dir, reporter)
    except Exception as e:
        print(f"\n[CRITICAL ERROR] Dendrogram Pipeline Failed: {e}")
        traceback.print_exc()
        sys.exit(1)


# =============================================================================
# SECTION 6: HTML APPLICATION TEMPLATE (EMBEDDED D3.JS ENGINE)
# =============================================================================

HTML_APP_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>PFAS Degradation Dendrogram Explorer</title>
    <script src="https://cdn.tailwindcss.com"></script>
    <script src="https://d3js.org/d3.v7.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/PapaParse/5.3.2/papaparse.min.js"></script>
    <script src="https://unpkg.com/@phosphor-icons/web"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/jszip/3.10.1/jszip.min.js"></script>

    <style>
        body { font-family: 'Inter', sans-serif; overflow: hidden; background-color: #f8fafc; }
        .node circle { cursor: pointer; stroke-width: 1.5px; transition: stroke-width 0.2s, stroke 0.2s; }
        .node circle:hover { stroke: #1e293b; stroke-width: 3px; }
        .node text { font-size: 11px; font-family: sans-serif; pointer-events: none; }
        .link { fill: none; stroke-width: 1.5px; opacity: 0.9; cursor: pointer; transition: stroke-width 0.2s; }
        .link:hover { stroke-width: 3.5px; }

        .dist-label-bg { font-size: 9px; fill: none; stroke: #f8fafc; stroke-width: 3px; stroke-linejoin: round; pointer-events: none; font-family: monospace; }
        .dist-label-fg { font-size: 9px; fill: #64748b; pointer-events: none; font-family: monospace; }

        .heatmap-cell-rect { cursor: pointer; transition: stroke 0.2s, opacity 0.3s, x 0.5s, y 0.5s, width 0.5s, height 0.5s; }
        .heatmap-cell-rect:hover { stroke: #1e293b; stroke-width: 2px; }

        ::-webkit-scrollbar { width: 6px; height: 6px; }
        ::-webkit-scrollbar-track { background: transparent; }
        ::-webkit-scrollbar-thumb { background: #cbd5e1; border-radius: 4px; }
        ::-webkit-scrollbar-thumb:hover { background: #94a3b8; }

        #sidebar-wrapper { transition: width 0.3s ease-in-out; }
        .sidebar-collapsed { width: 0 !important; border-right: none !important; }
        .sidebar-collapsed #sidebar-content { opacity: 0; pointer-events: none; }

        .tooltip {
            position: absolute;
            text-align: left;
            padding: 12px;
            font-size: 12px;
            background: rgba(15, 23, 42, 0.95);
            color: #f8fafc;
            border-radius: 8px;
            pointer-events: none;
            opacity: 0;
            transition: opacity 0.2s;
            box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
            z-index: 50;
            max-width: 280px;
        }
    </style>
    <!-- INJECT_EMBEDDED_STATE -->
</head>
<body class="flex h-screen w-full text-slate-800">

    <!-- Smart Sliding Sidebar -->
    <div id="sidebar-wrapper" class="relative h-full z-30 w-80 shrink-0 bg-white/95 backdrop-blur shadow-[4px_0_24px_rgba(0,0,0,0.05)] border-r border-slate-200">

        <div id="sidebar-content" class="w-80 h-full flex flex-col overflow-y-auto transition-opacity duration-300">
            <div class="p-6 border-b border-slate-200 bg-white shrink-0">
                <h1 class="text-xl font-bold flex items-center gap-2 text-indigo-700 whitespace-nowrap">
                    <i class="ph ph-tree-structure text-2xl"></i>
                    Dendrogram Explorer
                </h1>
            </div>

            <div class="p-6 flex-1 flex flex-col gap-6">

                <div id="upload-section">
                    <div class="space-y-4 pt-2 border-slate-200">
                        <h2 class="text-sm font-semibold uppercase tracking-wider text-slate-400">Load Data</h2>
                        <div>
                            <label class="block text-xs font-medium text-slate-600 mb-1">Upload Tree (.tree) OR Sequences (.fasta)</label>
                            <input type="file" id="tree-upload" accept=".tree,.txt,.fasta,.fa,.faa" class="block w-full text-xs text-slate-500 file:mr-3 file:py-1.5 file:px-3 file:rounded-md file:border-0 file:text-xs file:font-semibold file:bg-slate-100 file:text-slate-700 hover:file:bg-slate-200 cursor-pointer border border-slate-200 rounded-md">
                        </div>
                        <div>
                            <label class="block text-xs font-medium text-slate-600 mb-1">Upload Boltz Master CSV</label>
                            <input type="file" id="csv-upload" accept=".csv" class="block w-full text-xs text-slate-500 file:mr-3 file:py-1.5 file:px-3 file:rounded-md file:border-0 file:text-xs file:font-semibold file:bg-slate-100 file:text-slate-700 hover:file:bg-slate-200 cursor-pointer border border-slate-200 rounded-md">
                        </div>
                        <button id="btn-plot-tree" class="w-full py-2.5 bg-emerald-500 hover:bg-emerald-600 text-white text-xs font-bold uppercase tracking-wide rounded-md transition-colors flex items-center justify-center gap-2 shadow-md">
                            <i class="ph ph-tree text-lg"></i> Plot Tree & Matrix
                        </button>
                    </div>
                </div>

                <!-- 2. Global Stats -->
                <div class="space-y-3 pt-4 border-t border-slate-200">
                    <h2 class="text-sm font-semibold uppercase tracking-wider text-slate-400">Matrix Statistics</h2>
                    <div class="grid grid-cols-2 gap-2 mb-2">
                        <div class="bg-slate-50 p-2 rounded border border-slate-100 text-center">
                            <span class="block text-[10px] text-slate-400 uppercase">Tree Proteins</span>
                            <span id="stat-total-proteins" class="block text-base font-bold text-slate-700">-</span>
                        </div>
                        <div class="bg-slate-50 p-2 rounded border border-slate-100 text-center">
                            <span class="block text-[10px] text-slate-400 uppercase">Visible Proteins</span>
                            <span id="stat-visible-proteins" class="block text-base font-bold text-slate-700">-</span>
                        </div>
                        <div class="bg-slate-50 p-2 rounded border border-slate-100 text-center col-span-2">
                            <span class="block text-[10px] text-slate-400 uppercase">Active Ligands</span>
                            <span id="stat-active-ligands" class="block text-base font-bold text-slate-700">-</span>
                        </div>
                    </div>
                    <div id="tier-breakdown"></div>
                </div>

                <!-- 3. Dynamic Ligand Filter -->
                <div class="space-y-2 pt-4 border-t border-slate-200">
                    <div class="flex items-center justify-between">
                        <h2 class="text-sm font-semibold uppercase tracking-wider text-slate-400">Ligands</h2>
                        <div class="flex gap-2 text-[10px] font-medium">
                            <button id="ligand-all" class="text-indigo-600 hover:text-indigo-800 transition">All</button>
                            <button id="ligand-none" class="text-slate-500 hover:text-slate-700 transition">None</button>
                        </div>
                    </div>
                    <div id="ligand-container" class="space-y-1 text-xs max-h-40 overflow-y-auto pr-2 bg-slate-50/50 rounded p-1 border border-slate-100">
                        <span class="text-slate-400 italic px-1">Awaiting Data...</span>
                    </div>
                </div>

                <!-- 4. Dynamic Tier Filter -->
                <div class="space-y-2 pt-4 border-t border-slate-200">
                    <div class="flex items-center justify-between">
                        <h2 class="text-sm font-semibold uppercase tracking-wider text-slate-400">Tiers (Colors)</h2>
                        <div class="flex gap-2 text-[10px] font-medium">
                            <button id="tier-all" class="text-indigo-600 hover:text-indigo-800 transition">All</button>
                            <button id="tier-none" class="text-slate-500 hover:text-slate-700 transition">None</button>
                        </div>
                    </div>
                    <div id="tier-container" class="space-y-1 text-xs max-h-40 overflow-y-auto pr-2 bg-slate-50/50 rounded p-1 border border-slate-100">
                        <span class="text-slate-400 italic px-1">Awaiting Data...</span>
                    </div>
                </div>

                <!-- 5. Export Module -->
                <div id="export-module" class="space-y-3 pt-4 border-t border-slate-200 pb-6">
                    <h2 class="text-sm font-semibold uppercase tracking-wider text-slate-400">Export Matrix</h2>
                    <div class="grid grid-cols-2 gap-2">
                        <button id="btn-export-svg" class="col-span-2 py-1.5 bg-slate-100 hover:bg-slate-200 text-slate-700 text-xs font-semibold rounded transition flex items-center justify-center gap-1 border border-slate-300">
                            <i class="ph ph-vector-two"></i> Vector SVG
                        </button>
                    </div>
                </div>

            </div>
        </div>

        <button id="btn-toggle-sidebar" class="absolute top-1/2 -right-8 w-8 h-16 bg-white border border-slate-200 border-l-0 rounded-r-lg shadow-[4px_0_10px_rgba(0,0,0,0.05)] flex items-center justify-center text-slate-500 hover:text-indigo-600 transition-colors z-40" style="transform: translateY(-50%);">
            <i id="sidebar-toggle-icon" class="ph ph-caret-left text-xl"></i>
        </button>
    </div>

    <!-- Main Visualization Area -->
    <div class="flex-1 relative bg-slate-50" id="canvas-container">
        <div class="absolute bottom-4 left-4 z-20 flex flex-col bg-white/90 backdrop-blur rounded-md shadow-sm border border-slate-200 p-3 gap-2 text-xs text-slate-500 font-medium">
            <span class="flex items-center gap-2"><i class="ph ph-cursor-click text-indigo-500"></i> Free-drag nodes in all directions</span>
            <span class="flex items-center gap-2"><i class="ph ph-mouse-right-click text-indigo-500"></i> Right-click branches to paint clade colours</span>
            <span class="flex items-center gap-2"><i class="ph ph-arrow-u-up-left text-indigo-500"></i> <b>Ctrl + Z</b> to undo any action</span>
        </div>

        <div class="absolute top-4 right-4 z-20 flex bg-white rounded-md shadow-sm border border-slate-200 p-1 gap-1">
            <button id="btn-rotate" class="p-2 hover:bg-slate-100 rounded text-slate-600 transition" title="Rotate Tree (Horizontal ↔ Vertical)"><i class="ph ph-arrows-clockwise text-lg"></i></button>
            <div class="w-px bg-slate-200 mx-1 my-1"></div>
            <button id="btn-circular" class="p-2 hover:bg-slate-100 rounded text-slate-600 transition" title="Circular / Radial Tree Layout"><i id="icon-circular" class="ph ph-circles-three-plus text-lg"></i></button>
            <div class="w-px bg-slate-200 mx-1 my-1"></div>
            <button id="btn-toggle-details" class="p-2 hover:bg-slate-100 rounded text-slate-600 transition" title="Toggle Full Details Mode"><i id="icon-toggle-details" class="ph ph-magnifying-glass-plus text-lg"></i></button>
            <div class="w-px bg-slate-200 mx-1 my-1"></div>
            <button id="btn-align-labels" class="p-2 hover:bg-slate-100 rounded text-slate-600 transition" title="Align labels →"><i id="icon-align-labels" class="ph ph-align-right text-lg"></i></button>
            <div class="w-px bg-slate-200 mx-1 my-1"></div>
            <button id="btn-reset" class="p-2 hover:bg-slate-100 rounded text-slate-600 transition" title="Reset Zoom"><i class="ph ph-corners-out text-lg"></i></button>
        </div>

        <div id="loading-overlay" class="absolute inset-0 bg-slate-50/80 backdrop-blur-sm z-30 flex flex-col items-center justify-center hidden">
            <i class="ph ph-spinner-gap animate-spin text-4xl text-indigo-600 mb-4"></i>
            <p class="text-sm font-semibold text-slate-600" id="loading-text">Processing Data...</p>
        </div>

        <svg id="dendrogram-svg" class="w-full h-full cursor-move"></svg>
    </div>

    <!-- Tooltip & Color Palette -->
    <div id="tooltip" class="tooltip"></div>
    <div id="colour-palette" class="absolute hidden bg-white shadow-[0_4px_15px_rgba(0,0,0,0.1)] border border-slate-200 rounded-md p-2 z-50 flex flex-wrap w-36 gap-1.5 cursor-pointer"></div>

    <script>
        function parseNewick(a) {
            let e = [], r = {}, s = a.split(/(;|\(|\)|,)/), t = 0;
            for (; t < s.length; t++) {
                let n = s[t];
                if ("(" === n) { let c = {}; r.children = [c]; e.push(r); r = c; }
                else if ("," === n) { let c = {}; e[e.length - 1].children.push(c); r = c; }
                else if (")" === n) { r = e.pop(); }
                else if (";" === n) break;
                else if (n.trim() !== "") {
                    let h = n.split(":");
                    r.name = h[0].trim();
                    if (h.length > 1) r.length = parseFloat(h[1]);
                }
            }
            return r;
        }

        const BASE_TIER_COLORS = {{ JSON_BASE_COLORS }};
        const LIGAND_COLORS = {{ JSON_LIGAND_COLORS }};
        // CFG.VIS_LIGAND_SHORT (single source of truth) - the ONE ligand shortener; display-only,
        // the raw Ligand_Name stays the key for filtering / colouring / matching.
        const LIGAND_SHORT = {{ JSON_LIGAND_SHORT }};
        function ligShort(n){ if(n==null) return n; const s=String(n).replace(/^\d+_/,''); return LIGAND_SHORT[s.toLowerCase()] || s; }
        const CLADE_COLORS = [...LIGAND_COLORS, ...{{ JSON_CLADE_COLORS }}];
        const TIER_RANKING = {{ JSON_TIER_ORDER }};
        const CELL_PADDING = 2;

        let state = {
            treeData: null, rawNwk: "", rawCsv: "", csvData: {},
            ligands: new Set(), activeLigands: new Set(),
            tiers: new Set(), activeTiers: new Set(),
            root: null, isDetailedView: false, isExpanded: true,
            orientation: 'horizontal',   // 'horizontal' | 'vertical' | 'circular'
            alignLabelsRight: false
        };

        let historyStack = [];

        const svg = d3.select("#dendrogram-svg");
        const container = document.getElementById('canvas-container');
        const tooltip = d3.select("#tooltip");
        const colorPalette = document.getElementById('colour-palette');
        const loadingOverlay = document.getElementById('loading-overlay');
        const loadingText = document.getElementById('loading-text');

        let contextNode = null;

        colorPalette.innerHTML = '';
        CLADE_COLORS.forEach(c => {
            let div = document.createElement('div');
            div.className = "w-6 h-6 rounded cursor-pointer hover:scale-110 transition-transform shadow-sm";
            div.style.backgroundColor = c;
            div.onclick = () => {
                if(contextNode) {
                    saveToHistory();
                    contextNode.each(child => child.data.cladeColor = c);
                    updateTree(state.root);
                }
                colorPalette.classList.add('hidden');
            };
            colorPalette.appendChild(div);
        });

        document.addEventListener('click', () => colorPalette.classList.add('hidden'));
        document.addEventListener('keydown', (e) => { if (e.ctrlKey && e.key.toLowerCase() === 'z') { e.preventDefault(); undo(); } });

        const zoom = d3.zoom().scaleExtent([0.002, 8]).on("zoom", (e) => g.attr("transform", e.transform));
        svg.call(zoom);
        const g = svg.append("g");

        const gLines = g.append("g").attr("class", "layer-lines");
        const gLinks = g.append("g").attr("class", "layer-links");
        const gNodes = g.append("g").attr("class", "layer-nodes");
        const gHeaders = g.append("g").attr("class", "layer-headers");

        function normalizeProtName(name) {
            if (!name) return "";
            return String(name).replace(/[^a-zA-Z0-9]/g, '').toLowerCase();
        }

        function showLoading(text) { loadingText.innerText = text; loadingOverlay.classList.remove('hidden'); }
        function hideLoading() { loadingOverlay.classList.add('hidden'); }

        function getColorForTier(tier) {
            if (BASE_TIER_COLORS[tier]) return BASE_TIER_COLORS[tier];
            let hash = 0;
            for (let i = 0; i < tier.length; i++) hash = tier.charCodeAt(i) + ((hash << 5) - hash);
            const c = (hash & 0x00FFFFFF).toString(16).toUpperCase();
            return "#" + "00000".substring(0, 6 - c.length) + c;
        }

        function getContrastColor(tier) {
            if (tier === TIER_RANKING[TIER_RANKING.length-1] || tier === "Unknown") return "#334155";
            return "#ffffff";
        }

        function saveToHistory() {
            let visState = {};
            state.root.each(d => {
                visState[d.data.uid] = { offsetX: d.data.offsetX, offsetY: d.data.offsetY, cladeColor: d.data.cladeColor, collapsed: !!d._children };
            });
            historyStack.push(JSON.stringify(visState));
            if(historyStack.length > 50) historyStack.shift();
        }

        function undo() {
            if(historyStack.length > 0) {
                let visState = JSON.parse(historyStack.pop());
                state.root.each(d => {
                    const s = visState[d.data.uid];
                    if(s) {
                        d.data.offsetX = s.offsetX; d.data.offsetY = s.offsetY; d.data.cladeColor = s.cladeColor;
                        if (s.collapsed && d.children) { d._children = d.children; d.children = null; }
                        else if (!s.collapsed && d._children) { d.children = d._children; d._children = null; }
                    }
                });
                updateTree(state.root);
            }
        }

        function triggerDownload(url, filename) {
            const a = document.createElement('a');
            a.href = url; a.download = filename;
            document.body.appendChild(a); a.click(); document.body.removeChild(a);
        }

        // Serialise the live dendrogram/matrix SVG (with print styles, tier legend and a padded
        // viewBox) into a standalone, self-contained SVG string - TRUE VECTOR, no canvas raster step.
        function buildExportSvgData() {
                const svgNode = document.getElementById("dendrogram-svg");
                const gNode = svgNode.querySelector("g");
                const originalViewBox = svgNode.getAttribute("viewBox");
                const originalTransform = gNode.getAttribute("transform");

                svgNode.setAttribute("xmlns", "http://www.w3.org/2000/svg");
                const styleElement = document.createElementNS("http://www.w3.org/2000/svg", "style");
                styleElement.textContent = `
                    text { font-family: sans-serif; }
                    .link { fill: none; stroke-width: 1.5px; opacity: 0.9; }
                    .dist-label-bg { font-size: 9px; fill: none; stroke: #f8fafc; stroke-width: 3px; stroke-linejoin: round; }
                    .dist-label-fg { font-size: 9px; fill: #64748b; }
                    .heatmap-header { font-size: 10px; font-weight: bold; }
                    .matrix-col-line, .matrix-cap-line { stroke-width: 1.5px; stroke-dasharray: 4 4; opacity: 0.9; }
                    .c-t1-fg { font-size: 10px; font-weight: bold; }
                    .c-t2-fg { font-size: 9px; }
                `;
                svgNode.insertBefore(styleElement, svgNode.firstChild);

                const legendGroup = document.createElementNS("http://www.w3.org/2000/svg", "g");
                let currentX = 0;
                Array.from(state.activeTiers).sort((a,b) => TIER_RANKING.indexOf(a) - TIER_RANKING.indexOf(b)).forEach(tier => {
                    const rect = document.createElementNS("http://www.w3.org/2000/svg", "rect");
                    rect.setAttribute("x", currentX); rect.setAttribute("y", 0); rect.setAttribute("width", 14); rect.setAttribute("height", 14); rect.setAttribute("rx", 3); rect.setAttribute("fill", getColorForTier(tier));
                    legendGroup.appendChild(rect);
                    const text = document.createElementNS("http://www.w3.org/2000/svg", "text");
                    text.setAttribute("x", currentX + 22); text.setAttribute("y", 11); text.setAttribute("font-size", "14px"); text.setAttribute("font-weight", "bold"); text.setAttribute("fill", "#334155"); text.textContent = tier.replace('_', ' ');
                    legendGroup.appendChild(text);
                    currentX += 120;
                });
                gNode.appendChild(legendGroup);

                const bbox = gNode.getBBox();
                legendGroup.setAttribute("transform", `translate(${bbox.x}, ${bbox.y + bbox.height + 60})`);
                const finalBbox = gNode.getBBox();

                const padding = 150;
                const finalWidth = finalBbox.width + padding*2;
                const finalHeight = finalBbox.height + padding*2;

                svgNode.setAttribute("viewBox", `${finalBbox.x - padding} ${finalBbox.y - padding} ${finalWidth} ${finalHeight}`);
                const oldWidth = svgNode.getAttribute("width");
                const oldHeight = svgNode.getAttribute("height");
                svgNode.setAttribute("width", finalWidth);
                svgNode.setAttribute("height", finalHeight);
                gNode.removeAttribute("transform");

                const svgData = new XMLSerializer().serializeToString(svgNode);

                gNode.removeChild(legendGroup);
                svgNode.removeChild(styleElement);
                if(oldWidth) svgNode.setAttribute("width", oldWidth); else svgNode.removeAttribute("width");
                if(oldHeight) svgNode.setAttribute("height", oldHeight); else svgNode.removeAttribute("height");
                if(originalViewBox) svgNode.setAttribute("viewBox", originalViewBox); else svgNode.removeAttribute("viewBox");
                gNode.setAttribute("transform", originalTransform);

                return { svgData, finalWidth, finalHeight };
        }

        document.getElementById('btn-export-svg').addEventListener('click', () => {
            if (!state.root) return alert("Please plot a tree first.");
            showLoading("Generating Vector SVG...");
            try {
                const { svgData } = buildExportSvgData();
                // Prepend the XML declaration so the file opens as a standalone SVG in any viewer/editor.
                const doc = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n' + svgData;
                const url = URL.createObjectURL(new Blob([doc], {type: 'image/svg+xml;charset=utf-8'}));
                triggerDownload(url, "PFAS_Matrix_Export.svg");
                URL.revokeObjectURL(url);
            } finally {
                hideLoading();
            }
        });

        // View Controls
        document.getElementById('btn-toggle-sidebar').addEventListener('click', () => {
            const sidebar = document.getElementById('sidebar-wrapper');
            const icon = document.getElementById('sidebar-toggle-icon');
            sidebar.classList.toggle('sidebar-collapsed');
            icon.classList.replace(sidebar.classList.contains('sidebar-collapsed') ? 'ph-caret-left' : 'ph-caret-right', sidebar.classList.contains('sidebar-collapsed') ? 'ph-caret-right' : 'ph-caret-left');
            setTimeout(() => document.getElementById('btn-reset').click(), 300);
        });

        document.getElementById('btn-rotate').addEventListener('click', () => {
            if(!state.root) return;
            state.orientation = state.orientation === 'horizontal' ? 'vertical' : 'horizontal';
            updateTree(state.root);
            setTimeout(() => document.getElementById('btn-reset').click(), 100);
        });

        document.getElementById('btn-toggle-details').addEventListener('click', () => {
            if(!state.root) return;
            state.isDetailedView = !state.isDetailedView;
            document.getElementById('icon-toggle-details').className = state.isDetailedView ? "ph ph-magnifying-glass-minus text-lg" : "ph ph-magnifying-glass-plus text-lg";
            updateTree(state.root);
            setTimeout(() => document.getElementById('btn-reset').click(), 100);
        });

        // Circular / radial layout toggle - fully independent single-click
        document.getElementById('btn-circular').addEventListener('click', () => {
            if (!state.root) return;
            const goingCircular = state.orientation !== 'circular';
            state.orientation = goingCircular ? 'circular' : 'horizontal';
            // Disable align-labels mode when entering circular (incompatible)
            if (goingCircular && state.alignLabelsRight) {
                state.alignLabelsRight = false;
                document.getElementById('btn-align-labels').style.background = '';
                document.getElementById('icon-align-labels').style.color = '';
            }
            const btn  = document.getElementById('btn-circular');
            const icon = document.getElementById('icon-circular');
            btn.style.background  = goingCircular ? '#eef2ff' : '';
            icon.style.color      = goingCircular ? '#6366f1' : '';
            updateTree(state.root);
            setTimeout(() => document.getElementById('btn-reset').click(), 150);
        });

        // Align-labels-right toggle - moves all leaf labels flush to the matrix edge
        document.getElementById('btn-align-labels').addEventListener('click', () => {
            if (!state.root) return;
            state.alignLabelsRight = !state.alignLabelsRight;
            const btn  = document.getElementById('btn-align-labels');
            const icon = document.getElementById('icon-align-labels');
            btn.style.background  = state.alignLabelsRight ? '#eef2ff' : '';
            icon.style.color      = state.alignLabelsRight ? '#6366f1' : '';
            // Reset layout so depth-step recalculates for the new mode
            state.root.each(d => { d.manualY = undefined; });
            updateTree(state.root);
            setTimeout(() => document.getElementById('btn-reset').click(), 100);
        });

        document.getElementById('btn-reset').addEventListener('click', () => {
            if (!state.root) return;
            // Circular layout: centre the radial tree at canvas centre
            if (state.orientation === 'circular') {
                svg.transition().duration(750).call(zoom.transform,
                    d3.zoomIdentity.translate(container.clientWidth / 2, container.clientHeight / 2).scale(0.9));
                return;
            }
            const DIM_LIG = state.isDetailedView ? 140 : 16;
            const isV = state.orientation === 'vertical';
            const MATRIX_GAP = (state.alignLabelsRight && !isV) ? 190 : 60;
            const matrixDim = state.activeLigands.size * DIM_LIG;

            const leaves = state.root.leaves();
            const maxDepthY = d3.max(leaves, d => d.finalY) || 0;
            const minLeafX = d3.min(leaves, d => d.finalX) || 0;
            const maxLeafX = d3.max(leaves, d => d.finalX) || 0;

            if (isV) {
                const totalHeight = maxDepthY + 60 + matrixDim + 80;
                const totalWidth = (maxLeafX - minLeafX) + 200;
                const scale = Math.min(1.2, container.clientWidth / totalWidth, container.clientHeight / totalHeight) * 0.9;
                const cx = (minLeafX + maxLeafX) / 2;
                const cy = totalHeight / 2;
                const tx = container.clientWidth / 2 - cx * scale;
                const ty = container.clientHeight / 2 - cy * scale + 50;
                svg.transition().duration(750).call(zoom.transform, d3.zoomIdentity.translate(tx, ty).scale(scale));
            } else {
                const totalWidth = maxDepthY + MATRIX_GAP + matrixDim + 80;
                const totalHeight = (maxLeafX - minLeafX) + 200;
                const scale = Math.min(1.2, container.clientWidth / totalWidth, container.clientHeight / totalHeight) * 0.9;
                const cx = totalWidth / 2;
                const cy = (minLeafX + maxLeafX) / 2;
                const tx = container.clientWidth / 2 - cx * scale + 50;
                const ty = container.clientHeight / 2 - cy * scale;
                svg.transition().duration(750).call(zoom.transform, d3.zoomIdentity.translate(tx, ty).scale(scale));
            }
        });

        // Filters
        document.getElementById('ligand-all').addEventListener('click', () => { document.querySelectorAll('.ligand-checkbox').forEach(cb => { cb.checked = true; state.activeLigands.add(cb.value); }); updateColorsAndStats(); });
        document.getElementById('ligand-none').addEventListener('click', () => { document.querySelectorAll('.ligand-checkbox').forEach(cb => cb.checked = false); state.activeLigands.clear(); updateColorsAndStats(); });
        document.getElementById('tier-all').addEventListener('click', () => { document.querySelectorAll('.tier-checkbox').forEach(cb => { cb.checked = true; state.activeTiers.add(cb.value); }); updateColorsAndStats(); });
        document.getElementById('tier-none').addEventListener('click', () => { document.querySelectorAll('.tier-checkbox').forEach(cb => cb.checked = false); state.activeTiers.clear(); updateColorsAndStats(); });

        document.getElementById('btn-plot-tree').addEventListener('click', () => {
            const treeFile = document.getElementById('tree-upload').files[0];
            const csvFile = document.getElementById('csv-upload').files[0];

            if (!treeFile || !csvFile) { alert("Please upload both Tree and CSV."); return; }
            showLoading("Parsing Files...");

            const readerT = new FileReader();
            readerT.onload = e => {
                state.rawNwk = e.target.result; state.treeData = parseNewick(state.rawNwk);
                const readerC = new FileReader();
                readerC.onload = e2 => {
                    state.rawCsv = e2.target.result;
                    state.ligands.clear(); state.tiers.clear(); state.csvData = {};
                    Papa.parse(state.rawCsv, {
                        header: true, dynamicTyping: true, skipEmptyLines: true, transformHeader: h => h.trim(),
                        chunk: function(results) { processCsvChunk(results.data); results.data = []; },
                        complete: function() { finalizeCsvProcessing(); initializeTree(); hideLoading(); }
                    });
                };
                readerC.readAsText(csvFile);
            };
            readerT.readAsText(treeFile);
        });

        function processCsvChunk(data) {
            data.forEach(row => {
                const protRaw = row['Protein_Name'] || row.protein || row.protein_id || (row.job_name ? row.job_name.split('_')[1] : null);
                const ligRaw = row['Ligand_Name'] || row.ligand || (row.job_name ? row.job_name.split('_')[2] : null);
                const tier = row.Degrader_Tier || row.degrader_tier || row.Tier || TIER_RANKING[TIER_RANKING.length-1];
                const score = row.ActiveSite_Conservation_Score || row.Binding_Probability_Score || 0;

                if (!protRaw || !ligRaw) return;
                const normProt = normalizeProtName(protRaw);
                const lig = String(ligRaw).trim();

                if (!state.csvData[normProt]) state.csvData[normProt] = {};
                state.csvData[normProt][lig] = { tier, score, originalProtName: protRaw };

                state.ligands.add(lig); state.tiers.add(tier);
            });
        }

        function finalizeCsvProcessing() {
            const sortedLigands = Array.from(state.ligands).sort((a, b) => a.localeCompare(b, undefined, {numeric: true, sensitivity: 'base'}));
            state.activeLigands = new Set(sortedLigands);

            const sortedTiers = Array.from(state.tiers).sort((a, b) => {
                const rankA = TIER_RANKING.indexOf(a) !== -1 ? TIER_RANKING.indexOf(a) : 99;
                const rankB = TIER_RANKING.indexOf(b) !== -1 ? TIER_RANKING.indexOf(b) : 99;
                if (rankA !== rankB) return rankA - rankB; return a.localeCompare(b);
            });
            state.activeTiers = new Set(sortedTiers);
            renderFilterUI(sortedLigands, sortedTiers);
        }

        function renderFilterUI(ligands, tiers) {
            const ligContainer = document.getElementById('ligand-container');
            ligContainer.innerHTML = ligands.map(lig => `
                <label class="flex items-center gap-2 cursor-pointer hover:bg-slate-100 p-1.5 px-2 rounded transition">
                    <input type="checkbox" class="ligand-checkbox w-3.5 h-3.5 rounded border-slate-300 text-indigo-600" value="${lig}" checked>
                    <span class="text-slate-600 font-medium select-none truncate" title="${lig}">${ligShort(lig)}</span>
                </label>
            `).join('');
            document.querySelectorAll('.ligand-checkbox').forEach(cb => cb.addEventListener('change', (e) => { if (e.target.checked) state.activeLigands.add(e.target.value); else state.activeLigands.delete(e.target.value); updateColorsAndStats(); }));

            const tierContainer = document.getElementById('tier-container');
            tierContainer.innerHTML = tiers.map(tier => `
                <label class="flex items-center gap-2 cursor-pointer hover:bg-slate-100 p-1.5 px-2 rounded transition">
                    <input type="checkbox" class="tier-checkbox w-3.5 h-3.5 rounded border-slate-300 text-indigo-600" value="${tier}" checked>
                    <div class="w-3 h-3 rounded-full shadow-sm" style="background-color: ${getColorForTier(tier)}"></div>
                    <span class="text-slate-600 font-medium select-none truncate" title="${tier}">${tier.replace('_', ' ')}</span>
                </label>
            `).join('');
            document.querySelectorAll('.tier-checkbox').forEach(cb => cb.addEventListener('change', (e) => { if (e.target.checked) state.activeTiers.add(e.target.value); else state.activeTiers.delete(e.target.value); updateColorsAndStats(); }));
        }

        function getBestActiveMatch(rawProteinName) {
            const normName = normalizeProtName(rawProteinName);
            if (!state.csvData[normName]) return { tier: "Unknown", ligand: null, score: 0, originalProtName: rawProteinName };

            let bestRank = 999; let bestTier = TIER_RANKING[TIER_RANKING.length-1]; let bestLigand = null; let bestScore = -999;
            const pData = state.csvData[normName];

            for (const lig of state.activeLigands) {
                if (pData[lig] && state.activeTiers.has(pData[lig].tier)) {
                    const ligTier = pData[lig].tier;
                    const effectiveRank = TIER_RANKING.indexOf(ligTier) !== -1 ? TIER_RANKING.indexOf(ligTier) : 99;
                    if (effectiveRank < bestRank || (effectiveRank === bestRank && pData[lig].score > bestScore)) {
                        bestRank = effectiveRank; bestTier = ligTier; bestLigand = lig; bestScore = pData[lig].score;
                    }
                }
            }
            if (bestLigand === null) return { tier: TIER_RANKING[TIER_RANKING.length-1], ligand: null, score: 0, originalProtName: rawProteinName };
            return { tier: bestTier, ligand: bestLigand, score: bestScore, originalProtName: pData[bestLigand].originalProtName };
        }

        function updateColorsAndStats() {
            if (!state.root) return;

            let totalProteins = state.root.leaves().length;
            let visibleProteins = 0; let tierCounts = {};
            TIER_RANKING.forEach(t => tierCounts[t] = 0);

            gNodes.selectAll("g.node").transition().duration(500).style("opacity", d => (d.children || d._children || getBestActiveMatch(d.data.name).ligand !== null) ? 1 : 0.1);
            gLinks.selectAll("g.dist-group").transition().duration(500).style("opacity", d => (d.target.children || d.target._children || getBestActiveMatch(d.target.data.name).ligand !== null) ? 1 : 0.1);
            gNodes.selectAll("circle.node-circle").transition().duration(500).style("fill", d => d._children ? (d.data.cladeColor || "#94a3b8") : (getBestActiveMatch(d.data.name).tier === "Unknown" ? "#e2e8f0" : getColorForTier(getBestActiveMatch(d.data.name).tier)));

            updateTree(state.root);

            state.root.leaves().forEach(leaf => {
                if (leaf.children || leaf._children) return;
                const match = getBestActiveMatch(leaf.data.name);
                if (match.ligand !== null || match.tier === "Unknown") visibleProteins++;

                const pData = state.csvData[normalizeProtName(leaf.data.name)];
                if (pData) {
                    state.activeLigands.forEach(lig => {
                        if (pData[lig] && state.activeTiers.has(pData[lig].tier) && pData[lig].tier !== "Unknown" && pData[lig].tier !== TIER_RANKING[TIER_RANKING.length-1]) {
                            tierCounts[pData[lig].tier] = (tierCounts[pData[lig].tier] || 0) + 1;
                        }
                    });
                }
            });

            document.getElementById('stat-total-proteins').innerText = totalProteins;
            document.getElementById('stat-visible-proteins').innerText = visibleProteins;
            document.getElementById('stat-active-ligands').innerText = `${state.activeLigands.size} / ${state.ligands.size}`;

            let breakdownHtml = `<div class="grid grid-cols-2 gap-1 mt-2 text-[10px]">`;
            TIER_RANKING.forEach(t => {
                if(t === TIER_RANKING[TIER_RANKING.length-1]) return;
                breakdownHtml += `<div class="flex justify-between items-center bg-white px-2 py-1.5 rounded border border-slate-100 shadow-sm"><span style="color:${getColorForTier(t)}" class="font-bold">${t.replace('_',' ')}</span> <span class="text-slate-700 font-bold">${tierCounts[t] || 0}</span></div>`;
            });
            document.getElementById('tier-breakdown').innerHTML = breakdownHtml + `</div>`;
        }

        function initializeTree() {
            if (!state.treeData) return;
            // Reset layout-mode toggles on new tree load
            state.orientation     = 'horizontal';
            state.alignLabelsRight = false;
            ['btn-circular','btn-align-labels'].forEach(id => {
                const b = document.getElementById(id);
                if (b) b.style.background = '';
            });
            ['icon-circular','icon-align-labels'].forEach(id => {
                const ic = document.getElementById(id);
                if (ic) ic.style.color = '';
            });
            state.root = d3.hierarchy(state.treeData, d => d.children);
            let i = 0;
            state.root.each(d => { d.data.uid = "node_" + i++; d.data.cladeColor = '#cbd5e1'; d.manualY = undefined; d.data.offsetX = 0; d.data.offsetY = 0; d._angle = 0; d._radius = 0; d.cx = 0; d.cy = 0; });
            if (state.root.children) state.root.children.forEach((child, idx) => child.each(d => d.data.cladeColor = CLADE_COLORS[idx % CLADE_COLORS.length]));
            state.root.x0 = 0; state.root.y0 = 0;
            updateTree(state.root);
            setTimeout(() => document.getElementById('btn-reset').click(), 100);
            updateColorsAndStats();
        }

        // ─── Circular / radial tree layout ────────────────────────────────────
        function updateTreeCircular(source, skipLayout) {
            const transDuration = skipLayout ? 0 : 500;
            const leaves        = state.root.leaves();
            const nLeaves       = leaves.length;
            if (nLeaves === 0) return;

            const visibleLigands = Array.from(state.activeLigands).sort(
                (a, b) => a.localeCompare(b, undefined, {numeric: true, sensitivity: 'base'}));
            const nLig = visibleLigands.length;

            const RING_W    = 18;                      // radial thickness of each ligand ring
            const RING_GAP  = 4;                       // gap between rings
            const RING_STEP = RING_W + RING_GAP;
            const R         = Math.max(250, nLeaves * 30);   // tree branch radius
            const LABEL_GAP = 105;                     // leaf-label text space before rings
            const RINGS_START = R + LABEL_GAP;         // first ring inner edge

            if (!skipLayout) {
                const treeLayout = d3.tree()
                    .size([2 * Math.PI * 0.97, R])     // slight gap at top for readability
                    .separation((a, b) => 1.2);
                treeLayout(state.root);
                // Offset so gap sits at 12 o'clock
                state.root.each(d => { d.x += Math.PI * 0.015; });
            }
            const nodes = state.root.descendants();
            const links = state.root.links();

            nodes.forEach(d => {
                if (!skipLayout) {
                    d._angle  = d.x;
                    d._radius = d.y;
                    d.cx = d.y * Math.cos(d.x - Math.PI / 2);
                    d.cy = d.y * Math.sin(d.x - Math.PI / 2);
                }
                d.finalX = d.cx + (d.data.offsetX || 0);
                d.finalY = d.cy + (d.data.offsetY || 0);
            });

            // Hide rectangular matrix elements
            gLines.selectAll("rect.matrix-col-bg").style("opacity", 0);
            gLines.selectAll("line.matrix-col-line").style("opacity", 0);
            gHeaders.selectAll("text.heatmap-header").style("opacity", 0);

            // ── Concentric ring track backgrounds ─────────────────────────────
            const ringTrackData = visibleLigands.map((lig, i) => ({lig, i}));
            const ringTracks = gLines.selectAll("circle.ring-track").data(ringTrackData, d => d.lig);
            ringTracks.enter().append("circle").attr("class", "ring-track")
                .style("pointer-events", "none")
                .merge(ringTracks)
                .attr("r",   d => RINGS_START + d.i * RING_STEP + RING_W / 2)
                .attr("fill", "none")
                .style("stroke", "#d1d5db")
                .style("stroke-opacity", 0.55)
                .attr("stroke-width", RING_W);
            ringTracks.exit().remove();

            // ── Arc cells (leaf × ligand) ─────────────────────────────────────
            const angSlot = (2 * Math.PI * 0.97 / nLeaves) * 0.92;
            const arcGen  = d3.arc();
            const arcData = [];
            leaves.forEach(leaf => {
                const pData = state.csvData[normalizeProtName(leaf.data.name)] || {};
                visibleLigands.forEach((lig, ligIdx) => {
                    const lData = pData[lig] || {tier: TIER_RANKING[TIER_RANKING.length-1], score: 0};
                    arcData.push({
                        uid: leaf.data.uid + '__' + lig,
                        leaf, lig, ligIdx,
                        tier: lData.tier, score: lData.score,
                        protName: lData.originalProtName || leaf.data.name
                    });
                });
            });

            const arcCells = gLines.selectAll("path.ring-cell").data(arcData, d => d.uid);
            arcCells.enter().append("path").attr("class", "ring-cell")
                .style("cursor", "pointer")
                .on("mouseover", function(event, d) {
                    tooltip.html(
                        `<strong class="text-indigo-300 border-b border-slate-700 pb-1 block mb-2">${d.protName}</strong>` +
                        `<div class="grid grid-cols-2 gap-x-4 gap-y-1">` +
                        `<span class="text-slate-400">Ligand:</span><span>${ligShort(d.lig)}</span>` +
                        `<span class="text-slate-400">Tier:</span>` +
                        `<span style="color:${getColorForTier(d.tier)};font-weight:bold">${d.tier.replace('_',' ')}</span>` +
                        `<span class="text-slate-400">Score:</span><span>${parseFloat(d.score||0).toFixed(2)}</span>` +
                        `</div>`
                    ).style("opacity", 1);
                })
                .on("mousemove", (event) => tooltip.style("left",(event.pageX+15)+"px").style("top",(event.pageY-15)+"px"))
                .on("mouseout",  () => tooltip.style("opacity", 0))
                .merge(arcCells)
                .transition().duration(transDuration)
                .style("opacity", d => {
                    if (!state.activeTiers.has(d.tier)) return 0.06;
                    const tierOpacity = {{ JSON_TIER_ALPHAS }};
                    return tierOpacity[d.tier] !== undefined ? tierOpacity[d.tier] : 0.5;
                })
                .attr("fill", d => LIGAND_COLORS[d.ligIdx % LIGAND_COLORS.length])
                .attr("stroke", "none")
                .attr("d", d => arcGen({
                    innerRadius: RINGS_START + d.ligIdx * RING_STEP,
                    outerRadius: RINGS_START + d.ligIdx * RING_STEP + RING_W,
                    startAngle:  d.leaf._angle - angSlot / 2,
                    endAngle:    d.leaf._angle + angSlot / 2
                }));
            arcCells.exit().remove();

            // ── Ring labels: curved textPath along each ring arc at 12 o'clock ──
            // Invisible guide arcs inside gHeaders define the curved path for textPath.
            // Coordinates are in gHeaders local space (tree-centre = origin).
            const guideArcs = gHeaders.selectAll("path.ring-arc-guide").data(ringTrackData, d => d.lig);
            guideArcs.enter().append("path").attr("class", "ring-arc-guide")
                .attr("fill", "none").attr("stroke", "none").style("pointer-events", "none")
                .merge(guideArcs)
                .attr("id", d => `ring-arc-guide-${d.i}`)
                .attr("d", d => {
                    const r  = RINGS_START + d.i * RING_STEP + RING_W / 2;  // on the ring band
                    const a1 = -Math.PI / 2 - Math.PI / 5;   // ±36° around 12 o'clock
                    const a2 = -Math.PI / 2 + Math.PI / 5;
                    return `M ${(r * Math.cos(a1)).toFixed(2)} ${(r * Math.sin(a1)).toFixed(2)} ` +
                           `A ${r.toFixed(2)} ${r.toFixed(2)} 0 0 1 ` +
                           `${(r * Math.cos(a2)).toFixed(2)} ${(r * Math.sin(a2)).toFixed(2)}`;
                });
            guideArcs.exit().remove();

            const ringLabels = gHeaders.selectAll("text.ring-label").data(ringTrackData, d => d.lig);
            const rlEnter = ringLabels.enter().append("text").attr("class", "ring-label")
                .style("font-size", "10px").style("font-family", "sans-serif")
                .style("font-weight", "bold").style("cursor", "pointer")
                .on("click", function(event, d) {
                    if (state.activeLigands.has(d.lig)) state.activeLigands.delete(d.lig);
                    else state.activeLigands.add(d.lig);
                    document.querySelectorAll('.ligand-checkbox')
                        .forEach(cb => { if (cb.value === d.lig) cb.checked = state.activeLigands.has(d.lig); });
                    updateColorsAndStats();
                });
            rlEnter.append("textPath")
                .attr("startOffset", "50%").attr("text-anchor", "middle");
            ringLabels.merge(rlEnter)
                .style("opacity", d => state.activeLigands.has(d.lig) ? 1 : 0.3)
                .style("fill", (d, i) => LIGAND_COLORS[i % LIGAND_COLORS.length])
                .style("paint-order", "stroke fill")
                .style("-webkit-text-stroke", "2px white")
                .select("textPath")
                .attr("href", d => `#ring-arc-guide-${d.i}`)
                .text(d => ligShort(d.lig));
            ringLabels.exit().remove();

            // ── Nodes ─────────────────────────────────────────────────────────
            const node = gNodes.selectAll("g.node").data(nodes, d => d.data.uid);
            const nodeEnter = node.enter().append("g").attr("class", "node")
                .attr("transform", () => `translate(${source.finalX || 0},${source.finalY || 0})`)
                .on("contextmenu", function(event, d) {
                    event.preventDefault(); event.stopPropagation(); contextNode = d;
                    colorPalette.style.left = event.pageX + "px"; colorPalette.style.top = event.pageY + "px";
                    colorPalette.classList.remove('hidden');
                })
                .on("mouseover", (event, d) => {
                    if (!d.data.name) return;
                    const match = getBestActiveMatch(d.data.name);
                    let tipHtml = `<strong class="text-indigo-300 border-b border-slate-700 pb-1 block mb-2">${d.data.name}</strong>`;
                    if (match.ligand) tipHtml += `<div class="grid grid-cols-2 gap-x-4 gap-y-1"><span class="text-slate-400">Best ligand:</span><span>${ligShort(match.ligand)}</span><span class="text-slate-400">Tier:</span><span style="color:${getColorForTier(match.tier)};font-weight:bold">${match.tier}</span></div>`;
                    else tipHtml += `<span class="text-slate-400 italic">Internal node</span>`;
                    tooltip.html(tipHtml).style("opacity", 1);
                })
                .on("mousemove", (event) => tooltip.style("left",(event.pageX+15)+"px").style("top",(event.pageY-15)+"px"))
                .on("mouseout",  () => tooltip.style("opacity", 0));

            nodeEnter.call(d3.drag()
                .on("start", () => saveToHistory())
                .on("drag",  (event, d) => {
                    d.data.offsetX = (d.data.offsetX||0) + event.dx;
                    d.data.offsetY = (d.data.offsetY||0) + event.dy;
                    updateTreeCircular(state.root, true);
                }));

            nodeEnter.append("circle").attr("class", "node-circle").attr("r", 1e-6)
                .on("click", (event, d) => {
                    if (event.defaultPrevented) return;
                    saveToHistory();
                    if (d.children) { d._children = d.children; d.children = null; }
                    else { d.children = d._children; d._children = null; }
                    updateTreeCircular(d);
                });
            nodeEnter.append("line").attr("class", "leaf-leader")
                .attr("x1", 0).attr("y1", 0).attr("x2", 0).attr("y2", 0)
                .style("pointer-events", "none").style("stroke-width", "1px").style("opacity", 0);
            nodeEnter.append("text").attr("class", "node-text").attr("dy", ".35em")
                .style("fill-opacity", 1e-6).style("fill", "#334155");

            const nodeUpdate = nodeEnter.merge(node).transition().duration(transDuration)
                .attr("transform", d => `translate(${d.finalX},${d.finalY})`);

            nodeUpdate.select("circle")
                .attr("r",  d => d.children || d._children ? 5 : (d.data.name ? 4 : 3))
                .attr("cx", 0)
                .style("fill",         d => d._children ? (d.data.cladeColor || "#94a3b8") : getColorForTier(getBestActiveMatch(d.data.name).tier))
                .style("stroke",       d => d._children ? d.data.cladeColor : "#fff")
                .style("stroke-width", "1.5px");

            nodeUpdate.select("line.leaf-leader")
                .attr("transform", d => {
                    if (d.children || d._children || !d.data.name) return "rotate(0)";
                    const angleDeg = (d._angle - Math.PI / 2) * 180 / Math.PI;
                    const flip = d._angle > Math.PI;
                    return `rotate(${flip ? angleDeg + 180 : angleDeg})`;
                })
                .attr("x1", d => (d.children || d._children || !d.data.name) ? 0 : (d._angle > Math.PI ? -5 : 5))
                .attr("x2", d => (d.children || d._children || !d.data.name) ? 0 : (d._angle > Math.PI ? -11 : 11))
                .attr("y1", 0).attr("y2", 0)
                .style("stroke", d => {
                    if (d.children || d._children) return "none";
                    const m = getBestActiveMatch(d.data.name);
                    const li = visibleLigands.indexOf(m.ligand);
                    return li >= 0 ? LIGAND_COLORS[li % LIGAND_COLORS.length] : "#94a3b8";
                })
                .style("stroke-dasharray", d => (d.children || d._children || !d.data.name) ? "none" : "3 2")
                .style("opacity", d => (d.children || d._children || !d.data.name) ? 0 : 0.55);

            nodeUpdate.select("text")
                .attr("transform", d => {
                    const angleDeg = (d._angle - Math.PI / 2) * 180 / Math.PI;
                    const flip = d._angle > Math.PI;
                    return `rotate(${flip ? angleDeg + 180 : angleDeg})`;
                })
                .attr("x", d => (d.children || d._children) ? 0 : (d._angle > Math.PI ? -14 : 14))
                .attr("y",           d => (d.children || d._children) ? -12 : 0)
                .attr("text-anchor", d => (d.children || d._children) ? "middle" : (d._angle > Math.PI ? "end" : "start"))
                .style("font-size",  d => (d.children || d._children) ? "9px" : "10px")
                .style("paint-order", "stroke fill")
                .style("stroke", "#ffffff").style("stroke-width", "3px").style("stroke-linejoin", "round")
                .style("fill", d => {
                    if (d.children || d._children || !d.data.name) return "#334155";
                    const m = getBestActiveMatch(d.data.name);
                    const li = visibleLigands.indexOf(m.ligand);
                    return li >= 0 ? LIGAND_COLORS[li % LIGAND_COLORS.length] : "#334155";
                })
                .text(d => { const n = d.data.name || ''; return n.length > 16 ? n.substring(0, 15) + '…' : n; })
                .style("fill-opacity", 1);

            // Remove rectangular heatmap cells and any align-stem lines from prior rectangular render
            nodeUpdate.selection().selectAll("g.heatmap-cell-group").remove();
            nodeUpdate.selection().selectAll("line.align-stem").style("opacity", 0);

            // ── Links (radial elbow: arc at parent radius + radial to child) ──
            function circularDiag(lnk) {
                const s = lnk.source, t = lnk.target;
                const sa = s._angle - Math.PI / 2;
                const da = t._angle - Math.PI / 2;
                const sr = s._radius;
                const midX = sr * Math.cos(da), midY = sr * Math.sin(da);
                const dAngle = da - sa;
                const largeArc = Math.abs(dAngle) > Math.PI ? 1 : 0;
                const sweep    = dAngle > 0 ? 1 : 0;
                return `M ${s.finalX} ${s.finalY}` +
                       ` A ${sr} ${sr} 0 ${largeArc} ${sweep} ${midX} ${midY}` +
                       ` L ${t.finalX} ${t.finalY}`;
            }

            const link = gLinks.selectAll("path.link").data(links, d => d.target.data.uid);
            link.enter().insert("path", "g").attr("class", "link")
                .attr("stroke", d => d.target.data.cladeColor || '#cbd5e1')
                .attr("d", "M 0 0")
                .on("contextmenu", function(event, d) {
                    event.preventDefault(); event.stopPropagation(); contextNode = d.target;
                    colorPalette.style.left = event.pageX + "px"; colorPalette.style.top = event.pageY + "px";
                    colorPalette.classList.remove('hidden');
                })
                .merge(link).transition().duration(transDuration)
                .attr("stroke", d => d.target.data.cladeColor || '#cbd5e1')
                .attr("d", circularDiag);
            link.exit().transition().duration(transDuration).attr("d", "M 0 0").remove();

            const distGroupC = gLinks.selectAll("g.dist-group").data(links, d => d.target.data.uid);
            const distGroupCEnter = distGroupC.enter().insert("g", "g").attr("class", "dist-group")
                .attr("transform", "translate(0,0)").style("opacity", 1e-6);
            distGroupCEnter.append("text").attr("class", "dist-label-bg").attr("text-anchor", "middle")
                .text(d => (d.target.data.length !== undefined && d.target.data.length !== null) ? d.target.data.length.toFixed(3) : "");
            distGroupCEnter.append("text").attr("class", "dist-label-fg").attr("text-anchor", "middle")
                .text(d => (d.target.data.length !== undefined && d.target.data.length !== null) ? d.target.data.length.toFixed(3) : "");
            distGroupCEnter.merge(distGroupC).transition().duration(transDuration)
                .attr("transform", d => {
                    const da = d.target._angle - Math.PI / 2;
                    const sr = d.source._radius;
                    const mx = (sr * Math.cos(da) + d.target.finalX) / 2;
                    const my = (sr * Math.sin(da) + d.target.finalY) / 2;
                    return `translate(${mx}, ${my})`;
                })
                .style("opacity", 1);
            distGroupC.exit().transition().duration(transDuration).style("opacity", 1e-6).remove();

            node.exit().transition().duration(transDuration).style("opacity", 1e-6).remove();

            nodes.forEach(d => { d.x0 = d.finalX; d.y0 = d.finalY; });
        }

        function updateTree(source, skipLayout = false) {
            if (state.orientation === 'circular') { updateTreeCircular(source, skipLayout); return; }
            // Restore rectangular matrix elements; clear any circular-mode ring elements
            gLines.selectAll("rect.matrix-col-bg").style("opacity", null);
            gLines.selectAll("line.matrix-col-line").style("opacity", null);
            gHeaders.selectAll("text.heatmap-header").style("opacity", null);
            gLines.selectAll("circle.ring-track").remove();
            gLines.selectAll("path.ring-cell").remove();
            gHeaders.selectAll("text.ring-label").remove();
            gHeaders.selectAll("path.ring-arc-guide").remove();

            const isV = state.orientation === 'vertical';
            const NODE_SPACING = state.isDetailedView ? 70 : 45;
            const DIM_LIG = state.isDetailedView ? 180 : 28;
            const DIM_PROT = state.isDetailedView ? 40 : 14;
            const transDuration = skipLayout ? 0 : 500;
            const depthStep = state.alignLabelsRight ? 80 : 140;
            const MATRIX_GAP = (state.alignLabelsRight && !isV) ? 190 : 60;

            if (!skipLayout) { const treeLayout = d3.tree().nodeSize([NODE_SPACING, depthStep]); treeLayout(state.root); }

            const nodes = state.root.descendants();
            const links = state.root.links();
            const leaves = state.root.leaves();

            nodes.forEach(d => {
                if (!skipLayout) { if (d.manualY === undefined) d.manualY = d.y; else d.y = d.manualY; }
                d.finalX = d.x + (d.data.offsetX || 0); d.finalY = d.y + (d.data.offsetY || 0);
            });

            const maxDepthY = d3.max(leaves, d => d.finalY) || 0;
            const minLeafX = d3.min(leaves, d => d.finalX) || 0;
            const maxLeafX = d3.max(leaves, d => d.finalX) || 0;
            const visibleLigands = Array.from(state.activeLigands).sort((a, b) => a.localeCompare(b, undefined, {numeric: true, sensitivity: 'base'}));

            // GRID LINES
            const colBgs = gLines.selectAll("rect.matrix-col-bg").data(visibleLigands, d => d);
            colBgs.enter().append("rect").attr("class", "matrix-col-bg").style("fill", (d, i) => i % 2 === 0 ? "#00000000" : "#0f172a08")
                .merge(colBgs).transition().duration(transDuration)
                .attr("x", (d, i) => isV ? minLeafX - NODE_SPACING : maxDepthY + MATRIX_GAP + (i * DIM_LIG) - CELL_PADDING/2)
                .attr("y", (d, i) => isV ? maxDepthY + MATRIX_GAP + (i * DIM_LIG) - CELL_PADDING/2 : minLeafX - NODE_SPACING)
                .attr("width", isV ? maxLeafX - minLeafX + 2*NODE_SPACING : DIM_LIG)
                .attr("height", isV ? DIM_LIG : maxLeafX - minLeafX + 2*NODE_SPACING)
                .style("opacity", state.isDetailedView ? 1 : 0);
            colBgs.exit().remove();

            const colLines = gLines.selectAll("line.matrix-col-line").data(visibleLigands, d => d);
            colLines.enter().append("line").attr("class", "matrix-col-line").style("stroke-dasharray", "4 4").style("stroke-width", "1.5px").style("opacity", 0.9)
                .merge(colLines).transition().duration(transDuration)
                .attr("x1", (d, i) => isV ? minLeafX - NODE_SPACING : maxDepthY + MATRIX_GAP + (i * DIM_LIG) - CELL_PADDING/2)
                .attr("y1", (d, i) => isV ? maxDepthY + MATRIX_GAP + (i * DIM_LIG) - CELL_PADDING/2 : minLeafX - NODE_SPACING)
                .attr("x2", (d, i) => isV ? maxLeafX + NODE_SPACING : maxDepthY + MATRIX_GAP + (i * DIM_LIG) - CELL_PADDING/2)
                .attr("y2", (d, i) => isV ? maxDepthY + MATRIX_GAP + (i * DIM_LIG) - CELL_PADDING/2 : maxLeafX + NODE_SPACING)
                .style("stroke", (d, i) => LIGAND_COLORS[i % LIGAND_COLORS.length]);
            colLines.exit().remove();

            // HEADERS
            const headers = gHeaders.selectAll("text.heatmap-header").data(visibleLigands, d => d);
            headers.enter().append("text").attr("class", "heatmap-header").style("font-size", "10px").style("font-family", "sans-serif").style("font-weight", "bold").style("cursor", "pointer")
                .on("click", function(event, d) {
                    if(state.activeLigands.has(d)) state.activeLigands.delete(d); else state.activeLigands.add(d);
                    document.querySelectorAll('.ligand-checkbox').forEach(cb => { if(cb.value === d) cb.checked = state.activeLigands.has(d); });
                    updateColorsAndStats();
                })
                .merge(headers).transition().duration(transDuration).text(d => ligShort(d))
                .attr("transform", (d, i) => isV ? `translate(${maxLeafX + NODE_SPACING + 10}, ${maxDepthY + MATRIX_GAP + (i * DIM_LIG) + DIM_LIG/2 + 3}) rotate(0)` : `translate(${maxDepthY + MATRIX_GAP + (i * DIM_LIG) + DIM_LIG/2 + 3}, ${minLeafX - NODE_SPACING - 10}) rotate(-90)`)
                .style("fill", (d, i) => LIGAND_COLORS[i % LIGAND_COLORS.length]).style("text-anchor", "start").style("dominant-baseline", "middle").style("opacity", d => state.activeLigands.has(d) ? 1 : 0.3);
            headers.exit().remove();

            // NODES
            const node = gNodes.selectAll("g.node").data(nodes, d => d.data.uid);
            const nodeEnter = node.enter().append("g").attr("class", "node")
                .attr("transform", d => `translate(${isV ? source.finalX : source.finalY},${isV ? source.finalY : source.finalX})`)
                .on("contextmenu", function(event, d) { event.preventDefault(); event.stopPropagation(); contextNode = d; colorPalette.style.left = event.pageX + "px"; colorPalette.style.top = event.pageY + "px"; colorPalette.classList.remove('hidden'); })
                .on("mouseover", (event, d) => {
                    if (!d.data.name) return;
                    const match = getBestActiveMatch(d.data.name);
                    let tipHtml = `<strong class="text-indigo-300 border-b border-slate-700 pb-1 block mb-2">${d.data.name}</strong>`;
                    if (match.ligand) tipHtml += `<div class="grid grid-cols-2 gap-x-4 gap-y-1"><span class="text-slate-400">Ligand:</span> <span>${ligShort(match.ligand)}</span><span class="text-slate-400">Tier:</span> <span style="color:${getColorForTier(match.tier)}; font-weight:bold">${match.tier}</span></div>`;
                    else tipHtml += `<span class="text-slate-400 italic">Internal Node / No Filtered Match</span>`;
                    tooltip.html(tipHtml).style("opacity", 1);
                }).on("mousemove", (event) => tooltip.style("left", (event.pageX + 15) + "px").style("top", (event.pageY - 15) + "px")).on("mouseout", () => tooltip.style("opacity", 0));

            nodeEnter.call(d3.drag().on("start", () => saveToHistory()).on("drag", (event, d) => {
                if (isV) { d.data.offsetX = (d.data.offsetX || 0) + event.dx; d.data.offsetY = (d.data.offsetY || 0) + event.dy; }
                else { d.data.offsetY = (d.data.offsetY || 0) + event.dx; d.data.offsetX = (d.data.offsetX || 0) + event.dy; }
                updateTree(state.root, true);
            }));

            nodeEnter.append("circle").attr("class", "node-circle").attr("r", 1e-6)
                .on("click", (event, d) => { if (event.defaultPrevented) return; saveToHistory(); if (d.children) { d._children = d.children; d.children = null; } else { d.children = d._children; d._children = null; } updateTree(d); });
            nodeEnter.append("line").attr("class", "align-stem")
                .style("pointer-events", "none").style("opacity", 0);
            nodeEnter.append("text").attr("class", "node-text").attr("dy", ".35em").style("fill-opacity", 1e-6).style("fill", "#334155");

            const nodeUpdate = nodeEnter.merge(node).transition().duration(transDuration)
                .attr("transform", d => `translate(${isV ? d.finalX : d.finalY},${isV ? d.finalY : d.finalX})`);
            nodeUpdate.select("circle")
                .attr("r", d => d.data.name ? 6 : 4)
                .attr("cx", d => {
                    const isLeaf = !(d.children || d._children);
                    return (!isV && isLeaf && state.alignLabelsRight) ? maxDepthY - d.finalY : 0;
                })
                .style("fill", d => d._children ? (d.data.cladeColor || "#94a3b8") : (getBestActiveMatch(d.data.name).tier === "Unknown" ? "#e2e8f0" : getColorForTier(getBestActiveMatch(d.data.name).tier)))
                .style("stroke", d => d._children ? d.data.cladeColor : "#fff");
            // alignLabelsRight: horizontal leaf labels flush to matrix edge (maxDepthY+120)
            nodeUpdate.select("text")
                .attr("transform", d => isV ? "rotate(90)" : "rotate(0)")
                .attr("x", d => {
                    const isLeaf = !(d.children || d._children);
                    // Branch extends to maxDepthY; label sits just past it (independent of node position)
                    if (!isV && isLeaf && state.alignLabelsRight) return maxDepthY + 20 - d.finalY;
                    return isV ? 0 : (isLeaf ? 10 : -10);
                })
                .attr("y", d => isV ? ((d.children || d._children) ? -15 : 15) : 0)
                .attr("dx", d => (isV && !(d.children || d._children)) ? "10px" : "0px")
                .attr("text-anchor", d => {
                    const isLeaf = !(d.children || d._children);
                    if (!isV && isLeaf && state.alignLabelsRight) return "start";
                    return isV ? ((d.children || d._children) ? "middle" : "start") : ((d.children || d._children) ? "end" : "start");
                })
                .style("paint-order", "stroke fill")
                .style("stroke", "#ffffff").style("stroke-width", "3px").style("stroke-linejoin", "round")
                .text(d => { const n = d.data.name || ''; return n.length > 20 ? n.substring(0, 19) + '…' : n; })
                .style("fill-opacity", 1);

            // align-stem hidden: the extended bezier link IS the alignment line now
            nodeUpdate.selection().select("line.align-stem").style("opacity", 0);

            // MATRIX CELLS
            const cellGroups = nodeUpdate.selection().selectAll("g.heatmap-cell-group").data(d => {
                if (d.children || d._children) return [];
                const pData = state.csvData[normalizeProtName(d.data.name)] || {};
                return visibleLigands.map(lig => {
                    const lData = pData[lig] || {tier: TIER_RANKING[TIER_RANKING.length-1], score: 0, originalProtName: d.data.name};
                    return { lig: lig, tier: lData.tier, score: lData.score, parentY: d.finalY, protName: lData.originalProtName || d.data.name };
                });
            }, d => d.lig);

            const cellEnter = cellGroups.enter().append("g").attr("class", "heatmap-cell-group")
                .on("mouseover", function(event, d) {
                    tooltip.html(`<strong class="text-indigo-300 border-b border-slate-700 pb-1 block mb-2">${d.protName}</strong><div class="grid grid-cols-2 gap-x-4 gap-y-1"><span class="text-slate-400">Target Ligand:</span> <span>${ligShort(d.lig)}</span><span class="text-slate-400">Tier:</span> <span style="color:${getColorForTier(d.tier)}; font-weight:bold">${d.tier.replace('_', ' ')}</span></div>`).style("opacity", 1);
                    d3.select(this).select("rect").style("stroke", "#1e293b").style("stroke-width", "2px");
                }).on("mousemove", (event) => tooltip.style("left", (event.pageX + 15) + "px").style("top", (event.pageY - 15) + "px")).on("mouseout", function() { tooltip.style("opacity", 0); d3.select(this).select("rect").style("stroke", "none"); });

            cellEnter.append("rect").attr("class", "heatmap-cell-rect").attr("rx", 2);
            cellEnter.append("text").attr("class", "c-t1-fg").style("font-family", "sans-serif").style("font-weight", "bold").style("pointer-events", "none");
            cellEnter.append("text").attr("class", "c-t2-fg").style("font-family", "sans-serif").style("pointer-events", "none");

            const cellUpdate = cellEnter.merge(cellGroups).transition().duration(transDuration)
                .attr("transform", (d, i) => isV ? `translate(0, ${(maxDepthY - d.parentY) + MATRIX_GAP + (i * DIM_LIG)})` : `translate(${(maxDepthY - d.parentY) + MATRIX_GAP + (i * DIM_LIG)}, 0)`)
                .style("opacity", d => state.activeTiers.has(d.tier) ? 1 : 0.1);

            cellUpdate.select("rect").attr("width", isV ? DIM_PROT - CELL_PADDING : DIM_LIG - CELL_PADDING).attr("height", isV ? DIM_LIG - CELL_PADDING : DIM_PROT - CELL_PADDING).attr("x", isV ? -(DIM_PROT - CELL_PADDING)/2 : 0).attr("y", isV ? 0 : -(DIM_PROT - CELL_PADDING)/2).style("fill", d => state.activeTiers.has(d.tier) ? getColorForTier(d.tier) : "transparent");
            cellUpdate.selectAll(".c-t1-fg").attr("transform", isV ? "rotate(90)" : "rotate(0)").attr("x", isV ? 6 : 4).attr("y", isV ? 8 : -2).style("font-size", "10px").style("fill", d => getContrastColor(d.tier)).style("opacity", state.isDetailedView ? 1 : 0).text(d => `${(d.protName || '').substring(0,10)} | ${d.tier.replace('_',' ')}`);
            cellUpdate.selectAll(".c-t2-fg").attr("transform", isV ? "rotate(90)" : "rotate(0)").attr("x", isV ? 6 : 4).attr("y", isV ? -4 : 10).style("font-size", "9px").style("fill", d => getContrastColor(d.tier)).style("opacity", state.isDetailedView ? 1 : 0).text(d => `${(d.lig || '').substring(0,12)} | ${parseFloat(d.score).toFixed(1)}`);
            cellGroups.exit().remove();

            // LINKS
            // diagonal is defined here so it can access maxDepthY for alignLabelsRight extension.
            // When alignLabelsRight is on, leaf branch bezier extends to maxDepthY (no dashed stem needed).
            function diagonal(s, d, isVert) {
                if (isVert) {
                    return `M ${s.finalX} ${s.finalY} C ${s.finalX} ${(s.finalY + d.finalY) / 2}, ${d.finalX} ${(s.finalY + d.finalY) / 2}, ${d.finalX} ${d.finalY}`;
                }
                const isLeaf = !(d.children || d._children);
                const ty = (state.alignLabelsRight && !isVert && isLeaf) ? maxDepthY : d.finalY;
                return `M ${s.finalY} ${s.finalX} C ${(s.finalY + ty) / 2} ${s.finalX}, ${(s.finalY + ty) / 2} ${d.finalX}, ${ty} ${d.finalX}`;
            }
            const link = gLinks.selectAll("path.link").data(links, d => d.target.data.uid);
            const linkEnter = link.enter().insert("path", "g").attr("class", "link").attr("stroke", d => d.target.data.cladeColor || '#cbd5e1').attr("d", d => diagonal({finalX: source.finalX || 0, finalY: source.finalY || 0, children: null, _children: null}, {finalX: source.finalX || 0, finalY: source.finalY || 0, children: null, _children: null}, isV))
                .on("contextmenu", function(event, d) { event.preventDefault(); event.stopPropagation(); contextNode = d.target; colorPalette.style.left = event.pageX + "px"; colorPalette.style.top = event.pageY + "px"; colorPalette.classList.remove('hidden'); });
            linkEnter.merge(link).transition().duration(transDuration).attr("stroke", d => d.target.data.cladeColor || '#cbd5e1').attr("d", d => diagonal(d.source, d.target, isV));
            link.exit().transition().duration(transDuration).attr("d", () => `M ${isV ? source.finalX : source.finalY} ${isV ? source.finalY : source.finalX} L ${isV ? source.finalX : source.finalY} ${isV ? source.finalY : source.finalX}`).remove();

            // DIST LABELS
            const distGroup = gLinks.selectAll("g.dist-group").data(links, d => d.target.data.uid);
            const distGroupEnter = distGroup.enter().insert("g", "g").attr("class", "dist-group").attr("transform", d => `translate(${isV ? source.finalX : source.finalY}, ${isV ? source.finalY : source.finalX}) rotate(${isV ? 90 : 0})`).style("opacity", 1e-6);
            distGroupEnter.append("text").attr("class", "dist-label-bg").attr("text-anchor", "middle").text(d => (d.target.data.length !== undefined && d.target.data.length !== null) ? d.target.data.length.toFixed(3) : "");
            distGroupEnter.append("text").attr("class", "dist-label-fg").attr("text-anchor", "middle").text(d => (d.target.data.length !== undefined && d.target.data.length !== null) ? d.target.data.length.toFixed(3) : "");
            distGroupEnter.merge(distGroup).transition().duration(transDuration)
                .attr("transform", d => isV ? `translate(${(d.source.finalX + d.target.finalX)/2 + 4}, ${(d.source.finalY + d.target.finalY)/2}) rotate(90)` : `translate(${(d.source.finalY + d.target.finalY)/2}, ${(d.source.finalX + d.target.finalX)/2 - 5}) rotate(0)`)
                .style("opacity", d => (d.target.children || d.target._children || getBestActiveMatch(d.target.data.name).ligand !== null) ? 1 : 0.1);
            distGroup.exit().transition().duration(transDuration).style("opacity", 1e-6).remove();
            node.exit().transition().duration(transDuration).attr("transform", d => `translate(${isV ? source.finalX : source.finalY},${isV ? source.finalY : source.finalX})`).style("opacity", 1e-6).remove();

            nodes.forEach(d => { d.x0 = d.finalX; d.y0 = d.finalY; });
        }

        // diagonal is now defined inside updateTree (needs maxDepthY for alignLabelsRight extension)

        // --- App Bootloader (Restores Embedded State Automatically) ---
        window.onload = function() {
            if (window.EMBEDDED_STATE) {
                document.getElementById('upload-section').style.display = 'none';
                state.rawNwk = window.EMBEDDED_STATE.tree || window.EMBEDDED_STATE.nwk;
                state.treeData = parseNewick(state.rawNwk);
                state.rawCsv = window.EMBEDDED_STATE.csv;

                Papa.parse(state.rawCsv, {
                    header: true, dynamicTyping: true, skipEmptyLines: true, transformHeader: h => h.trim(),
                    chunk: function(results) { processCsvChunk(results.data); results.data = []; },
                    complete: function() { finalizeCsvProcessing(); initializeTree(); hideLoading(); }
                });
            }
        };
    </script>
</body>
</html>
"""

if __name__ == "__main__":
    _t0 = _time.perf_counter()
    main()
    _utils_mod.print_elapsed(_t0, "04_Dendrogram_DeFluorX.py")

