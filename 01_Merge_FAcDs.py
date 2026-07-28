#!/usr/bin/env python3
"""
===============================================================================
FAcDs Pipeline  |  Step 01  |  FASTA Sequence Merge & Deduplication
===============================================================================
Constructs a high-quality, non-redundant FASTA dataset by merging multiple
protein sequence sources using a strict, master-guided deduplication strategy.
Produces a merged FASTA, a detailed log, and a high-resolution QC dashboard.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 30 July 2026 <─────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Script        : 01_Merge_FAcDs.py
  Role          : Sequence merger and pre-processing pipeline wrapper.
  Imports from  : 00_02_Project_Utils_FAcDs.py  (clean_spines, print_script_banner,
                  print_elapsed, SEPARATOR_HEAVY/LIGHT/DASH)
                  00_01_Project_Config_FAcDs.py  (CFG - PREP_AMBIGUOUS_AA QC, CPU reserve)
  Reads         : User-supplied *.fasta files (master + secondary)
  Writes        : <output>.fasta   - merged, deduplicated sequence set
                  <output>.log     - inclusion/exclusion statistics
                  <output>.png     - QC dashboard (throughput + KDE)
  Upstream      : None (standalone data-curation step)
  Downstream    : 02_Production_FAcDs.py → consumes the merged FASTA as Boltz-2 input
───────────────────────────────────────────────────────────────────────────────

── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Sequence Identity Threshold: Assumes exact sequence matches for
     deduplication; does not handle fuzzy matching or SNP variants.
  2. FASTA Formatting: Highly sensitive to header formatting (expects UniProt/
     NCBI standard pipes).
  3. Memory: Input is streamed record-by-record (SeqIO.parse) and only sequence
     hashes + IDs are retained for deduplication, so the footprint scales with the
     UNIQUE-sequence count, not the raw file size - large FASTAs stream fine.
───────────────────────────────────────────────────────────────────────────────

Usage:
    python 01_Merge_FAcDs.py --master A_Labelled_15-Seq.fasta --secondary B_Downloaded-Blast_Uniprot_NCBI.fasta --output C_INP_Merged_for_Boltz-2.fasta

    (bash multi-line - use a single backslash, not double \\):
    python 01_Merge_FAcDs.py --master A_Labelled_15-Seq.fasta \
                       --secondary B_Downloaded-Blast_Uniprot_NCBI.fasta \
                       --output C_INP_Merged_for_Boltz-2.fasta

Purpose:
    Intended for curating aligned enzyme sequence sets prior to structure
    prediction, modelling, or phylogenetic analysis.

-------------------------------------------------------------------------------
Input FASTA Files:
    1. File 1 (MASTER):
       - Primary, trusted aligned sequence set.
       - No length filtering applied; basic QC (internal-stop / ambiguous-residue
         checks) and exact-sequence deduplication still apply.

    2. File 2 (Secondary):
       - Additional curated sequences.
       - Exact duplicates of File 1 are removed; remaining sequences are
         subjected to strict length filtering and multi-level deduplication.

-------------------------------------------------------------------------------
Output:
    - A single merged, non-redundant FASTA file.
    - A detailed console log (.log) summarising inclusion/exclusion statistics.
    - A High-Resolution PNG Visualisation Dashboard.
-------------------------------------------------------------------------------
Scientific References:
    1. Sequence parsing (Biopython SeqIO):
       - Cock, P.J.A. et al. (2009) Biopython. Bioinformatics 25:1422–1423.
       - DOI: https://doi.org/10.1093/bioinformatics/btp163
    2. Numerical & statistical tooling:
       - Harris, C.R. et al. (2020) Array programming with NumPy. Nature 585:357–362.
       - DOI: https://doi.org/10.1038/s41586-020-2649-2
       - Virtanen, P. et al. (2020) SciPy 1.0. Nature Methods 17:261–272.
       - DOI: https://doi.org/10.1038/s41592-019-0686-2
    3. Plotting (QC dashboard):
       - Hunter, J.D. (2007) Matplotlib. Comput Sci Eng 9:90–95.
       - DOI: https://doi.org/10.1109/MCSE.2007.55
-------------------------------------------------------------------------------
"""

# =============================================================================
# SECTION 1: CONFIGURATION & IMPORTS
# =============================================================================

# -------------------------------------------------------------------------------
# Step 1.1: Standard Library Imports
# -------------------------------------------------------------------------------
import argparse
import sys
import logging
import traceback
import re
from pathlib import Path
from typing import Set, Tuple, Dict
import hashlib

# -------------------------------------------------------------------------------
# Step 1.2: Scientific Stack (Matplotlib configured for headless/HPC servers)
# -------------------------------------------------------------------------------
"""
CPU usage cap (total cores − 2; mirrors CFG.PREP_CPU_RESERVE). Reserve 2 cores
for OS/desktop stability by limiting the thread-pool maths libraries (BLAS / MKL
/ OpenMP / NumExpr). Must precede numpy import to take effect; setdefault()
preserves any value exported by the caller or pipeline runner.
"""
import os as _os
_CPU_CAP = str(max(1, (_os.cpu_count() or 4) - 2))
for _tv in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
    _os.environ.setdefault(_tv, _CPU_CAP)

import matplotlib
matplotlib.use("Agg")  # Critical: Must be set before importing pyplot
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import skew, gaussian_kde

# -------------------------------------------------------------------------------
# Step 1.3: Biopython
# -------------------------------------------------------------------------------
from Bio import SeqIO

# -------------------------------------------------------------------------------
# Step 1.4: Pipeline utilities (00_02) + config (00_01) via importlib
# -------------------------------------------------------------------------------
"""
CFG supplies only PREP_AMBIGUOUS_AA for QC (no geometric thresholds used here).
"""
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

_utils_mod      = _load_module("ProjectUtils", Path(__file__).resolve().parent / "00_02_Project_Utils_FAcDs.py")
_cfg_mod        = _load_module("ProjectConfig", Path(__file__).resolve().parent / "00_01_Project_Config_FAcDs.py")
CFG             = _cfg_mod.CFG()
clean_spines    = _utils_mod.clean_spines

# -------------------------------------------------------------------------------
# Step 1.5: Global Constants
# -------------------------------------------------------------------------------
AMBIGUOUS_AA = CFG.PREP_AMBIGUOUS_AA  # Non-standard amino acids for QC

# =============================================================================
# SECTION 2: LOGGING INFRASTRUCTURE
# =============================================================================

def setup_logger(output_path: Path) -> logging.Logger:
    """
    Initialises a dual-handler logger (File + Console).
    Log is saved beside the output FASTA as 00_Merge.log - the same 00_<StepName>.log naming every
    step uses, so the log is instantly identifiable across the pipeline.
    """
    log_file = output_path.parent / "00_Merge.log"
    logger = logging.getLogger("FASTA_Merger")
    logger.setLevel(logging.INFO)

    # Prevent duplicate logs in interactive environments (Jupyter/IPython)
    if logger.hasHandlers():
        logger.handlers.clear()

    # Handlers
    file_handler = logging.FileHandler(log_file, mode="w")
    stream_handler = logging.StreamHandler(sys.stdout)

    # Format
    formatter = logging.Formatter("%(message)s")
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


# =============================================================================
# SECTION 3: CORE LOGIC & VALIDATORS
# =============================================================================

def clean_sequence_str(seq_obj) -> str:
    """
    Normalisation:
    1. Convert to String
    2. Upper Case
    3. Remove Gaps (-)
    4. Remove Trailing Stop Codons (*)
    """
    seq = str(seq_obj).upper()
    seq = seq.replace("-", "")
    # Strip a single trailing stop only; a double stop ("**") intentionally leaves an
    # internal "*" so is_valid_protein rejects the sequence as an internal stop codon.
    if seq.endswith("*"):
        seq = seq[:-1]
    return seq

def is_valid_protein(seq_str: str, max_ambiguous_percent: float = CFG.PREP_MAX_AMBIGUOUS_PCT) -> Tuple[bool, str]:
    """
    Quality Control (QC):
    Checks for internal stops and excessive ambiguous residues.
    """
    if not seq_str:
        return False, "Empty sequence"
    if "*" in seq_str:
        return False, "Internal stop codon (*)"

    len_seq = len(seq_str)
    ambiguous_count = sum(1 for aa in seq_str if aa in AMBIGUOUS_AA)

    if (ambiguous_count / len_seq) * 100 > max_ambiguous_percent:
        return False, f"Too many ambiguous residues ({ambiguous_count})"

    return True, "Valid"

def clean_header(description: str) -> str:
    """
    Extracts the primary accession ID and standardises the format.
    Removes parentheticals like '(2)'.
    """
    if not description or not description.strip():
        return "Unknown_Seq"

    # 1. Splitting by space isolates the primary ID from the trailing metadata
    raw_id = description.split()[0]

    """
    2. Remove any bracketed numbers like (2), (3) etc.
    e.g. 'GOI1_(2)' becomes 'GOI1_'
    """
    clean_id = re.sub(r"\(\d+\)", "", raw_id)

    return clean_id


# =============================================================================
# SECTION 4: PIPELINE ENGINE
# =============================================================================

def process_and_write(
    input_path: Path,
    output_handle,
    label: str,
    seen_sequences: Set[str],
    seen_ids: Set[str],
    logger: logging.Logger,
    min_len: int = 0,
    max_len: int = float("inf"),
    filter_length: bool = False,
    keep_gaps: bool = False
) -> Dict:
    """
    The Core Processor:
    Streams a FASTA file -> Cleans -> Filters -> Deduplicates -> Writes.
    Returns a dictionary of statistics for the visualisation dashboard.
    """

    stats = {
        "total": 0,
        "kept": 0,
        "dupes": 0,
        "length_fail": 0,
        "quality_fail": 0,
        "sum_len": 0,
        "min_len": float("inf"),
        "max_len": 0,
        "retained_lengths": []
    }

    logger.info(f"\n... Processing {label} file: {input_path}")

    try:
        iterator = SeqIO.parse(input_path, "fasta")
    except Exception as e:
        logger.error(f"Error reading {input_path}: {e}")
        return stats

    for rec in iterator:
        stats["total"] += 1

        # -------------------------------------------------------------------------------
        # Step 4.1: Normalisation
        # -------------------------------------------------------------------------------
        clean_seq = clean_sequence_str(rec.seq)

        # -------------------------------------------------------------------------------
        # Step 4.2: Quality Control
        # -------------------------------------------------------------------------------
        valid, _reason = is_valid_protein(clean_seq)
        if not valid:
            stats["quality_fail"] += 1
            logger.debug(f"  Quality fail [{rec.id}]: {_reason}")   # keep the cause (stop codon vs ambiguous)
            continue

        # -------------------------------------------------------------------------------
        # Step 4.3: Length Filter (Optional)
        # -------------------------------------------------------------------------------
        seq_len = len(clean_seq)
        if filter_length:
            if not (min_len <= seq_len <= max_len):
                stats["length_fail"] += 1
                continue

        # -------------------------------------------------------------------------------
        # Step 4.4: Deduplication (Strict Identity)
        # -------------------------------------------------------------------------------
        # Hash the sequence to save RAM before checking the 'seen' set
        clean_hash = hashlib.sha256(clean_seq.encode("utf-8")).hexdigest()

        if clean_hash in seen_sequences:
            stats["dupes"] += 1
            continue

        # -------------------------------------------------------------------------------
        # Step 4.5: Header Preprocessing (Cleaning the name & Collision check)
        # -------------------------------------------------------------------------------
        clean_id = clean_header(rec.description)

        """
        Disambiguate identical headers on distinct sequences: two different sequences that
        carry the same cleaned name would otherwise collide on ID, so a numeric suffix is
        appended to keep every record's identifier unique.
        """
        original_clean_id = clean_id
        counter = 1
        while clean_id in seen_ids:
            # If the name exists, add a _1, _2, etc. until it is unique
            clean_id = f"{original_clean_id}_{counter}"
            counter += 1

        # -------------------------------------------------------------------------------
        # Step 4.6: Update Database
        # -------------------------------------------------------------------------------
        seen_ids.add(clean_id)
        seen_sequences.add(clean_hash)

        # Update Statistics
        stats["sum_len"] += seq_len
        if seq_len < stats["min_len"]: stats["min_len"] = seq_len
        if seq_len > stats["max_len"]: stats["max_len"] = seq_len
        stats["retained_lengths"].append(seq_len)


        # -------------------------------------------------------------------------------
        # Step 4.7: Write to Stream (Single line per sequence)
        # -------------------------------------------------------------------------------
        if keep_gaps:
            # Preserve gap characters, but still upper-case and strip a trailing stop
            # so the emitted record matches the QC that ran on clean_seq. Tolerate a
            # stop followed by trailing gaps (e.g. "…*-"), which a plain endswith("*")
            # would miss and leak an internal stop into the final FASTA.
            seq_to_write = re.sub(r"\*(-*)$", r"\1", str(rec.seq).upper())
        else:
            seq_to_write = clean_seq

        # Generate the sequential ID based on the current number of retained sequences
        seq_number = len(seen_ids)
        formatted_header = f"{seq_number}_{clean_id}"

        # Write header and the entire sequence on a single line
        output_handle.write(f">{formatted_header}\n{seq_to_write}\n")

        stats["kept"] += 1

    return stats


# =============================================================================
# SECTION 5: VISUALISATION MODULE (DASHBOARD)
# =============================================================================

def apply_clean_spines(ax):
    """Delegate to 00_02_Project_Utils.clean_spines (mandatory import)."""
    clean_spines(ax)

def generate_plots(s1: Dict, s2: Dict, output_path: Path, logger: logging.Logger):
    """
    Generates a high-resolution dashboard.
    Features smart label placement and colourful KDE fillings.
    """

    # -------------------------------------------------------------------------------
    # Step 5.1: Palette Definition (Unified Consistency)
    # -------------------------------------------------------------------------------
    PALETTE = CFG.MERGE_QC_COLOUR       # single source (00_01 §visual palettes)

    # Configure Matplotlib fonts
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]

    # Setup Data
    data_map = [
        {"label": "Master",    "stats": s1, "colour": PALETTE["Master"]},
        {"label": "Secondary", "stats": s2, "colour": PALETTE["Secondary"]}
    ]

    # Calculate Totals
    total_input = sum(d["stats"]["total"] for d in data_map)
    total_kept = sum(d["stats"]["kept"] for d in data_map)
    all_lengths = []
    for d in data_map:
        all_lengths.extend(d["stats"]["retained_lengths"])

    # -------------------------------------------------------------------------------
    # Step 5.2: Canvas Setup (Dense Collage)
    # -------------------------------------------------------------------------------
    fig = plt.figure(figsize=(16, 9), facecolor=PALETTE["Bg"])

    # Layout: Top row (Bar + Violin), Bottom row (KDE)
    gs = fig.add_gridspec(2, 2, height_ratios=[0.65, 1.35], width_ratios=[1.3, 0.7], hspace=0.12, wspace=0.1)

    ax1 = fig.add_subplot(gs[0, 0]) # Top Left: Throughput
    ax2 = fig.add_subplot(gs[0, 1]) # Top Right: Violins
    ax3 = fig.add_subplot(gs[1, :]) # Bottom: KDE

    # -------------------------------------------------------------------------------
    # Step 5.3: Subplot 1 - Pipeline Throughput (Bar Chart)
    # -------------------------------------------------------------------------------
    bar_labels = ["Master", "Secondary", "FINAL\nDATASET"]
    bar_inputs = [s1["total"], s2["total"], total_input]
    bar_kept   = [s1["kept"], s2["kept"], total_kept]
    bar_colours = [PALETTE["Master"], PALETTE["Secondary"], PALETTE["Total"]]

    y_pos = np.arange(len(bar_labels))
    height = 0.6

    # Input (Light)
    ax1.barh(y_pos, bar_inputs, height, color=bar_colours, alpha=0.2,
             edgecolor="none", label="Total Input")

    # Retained (Solid)
    ax1.barh(y_pos, bar_kept, height, color=bar_colours, alpha=1.0,
             edgecolor="black", linewidth=0.8, label="Final Retained")

    ax1.set_yticks(y_pos)
    ax1.set_yticklabels(bar_labels, fontweight="bold", fontsize=CFG.VIS_FONT_AXIS_LABEL)
    ax1.invert_yaxis()
    ax1.set_xlabel("Number of Sequences", fontweight="bold", fontsize=CFG.VIS_FONT_TICK)
    # No Title for density

    # --- Smart Annotation Logic ---
    max_val = max(bar_inputs) if bar_inputs else 1
    """
    Estimation: approximate text width relative to axis (approx 15-20%)
    This prevents cramming text into small bars
    """
    width_threshold = max_val * 0.18

    for i, (inp, kp) in enumerate(zip(bar_inputs, bar_kept)):
        label_text = f"{kp:,} / {inp:,}" if i < 2 else f"TOTAL: {kp:,}"

        # Priority 1: Inside the Kept Bar (Contrasting White Text)
        if kp > width_threshold:
             ax1.text(kp - (max_val*0.02), i, label_text, va="center", ha="right",
                      fontsize=CFG.VIS_FONT_TICK, fontweight="bold", color="white")

        # Priority 2: Inside the Input Bar (Dark Text)
        elif inp > (kp + width_threshold):
             ax1.text(kp + (max_val*0.02), i, label_text, va="center", ha="left",
                      fontsize=CFG.VIS_FONT_TICK, fontweight="bold", color=PALETTE["Total"])

        # Priority 3: Outside (Dark Text)
        else:
             ax1.text(inp + (max_val*0.02), i, label_text, va="center", ha="left",
                      fontsize=CFG.VIS_FONT_TICK, fontweight="bold", color=PALETTE["Total"])

    apply_clean_spines(ax1)
    ax1.spines["left"].set_visible(False)
    ax1.tick_params(axis="y", length=0)

    ax1.legend(loc="upper right", frameon=True, fontsize=CFG.VIS_FONT_LEGEND, fancybox=True, framealpha=0.9)

    # -------------------------------------------------------------------------------
    # Step 5.4: Subplot 2 - Length Heterogeneity (Violin Plot)
    # -------------------------------------------------------------------------------
    violin_data = [d["stats"]["retained_lengths"] for d in data_map]
    safe_violin_data = [d if len(d) > 0 else [0] for d in violin_data]

    parts = ax2.violinplot(safe_violin_data, positions=[0, 1], vert=True, showmeans=False, showextrema=False, widths=0.75)

    for i, pc in enumerate(parts["bodies"]):
        pc.set_facecolor(data_map[i]["colour"])
        pc.set_alpha(0.7)
        pc.set_edgecolor("black")
        pc.set_linewidth(0.5)

    for i, d in enumerate(safe_violin_data):
        ax2.boxplot(d, positions=[i], widths=0.1, showfliers=False, patch_artist=True,
                    boxprops=dict(facecolor="white", alpha=0.6, linewidth=0.8),
                    whiskerprops=dict(linewidth=0.8), capprops=dict(linewidth=0.8),
                    medianprops=dict(color="black", linewidth=1.5))

    ax2.set_xticks([0, 1])
    ax2.set_xticklabels(["Master", "Secondary"], fontweight="bold", fontsize=CFG.VIS_FONT_TICK)
    ax2.set_ylabel("Length (AA)", fontweight="bold", fontsize=CFG.VIS_FONT_TICK)

    for i, tick in enumerate(ax2.get_xticklabels()):
        tick.set_color(data_map[i]["colour"])

    apply_clean_spines(ax2)
    ax2.grid(axis="y", linestyle=":", alpha=0.5)

    # -------------------------------------------------------------------------------
    # Step 5.5: Subplot 3 - Consolidated Architecture (KDE Plot)
    # -------------------------------------------------------------------------------
    ax3.set_xlabel("Sequence Length (Residues)", fontweight="bold", fontsize=CFG.VIS_FONT_AXIS_LABEL)
    ax3.set_ylabel("Density", fontweight="bold", fontsize=CFG.VIS_FONT_AXIS_LABEL)

    if len(all_lengths) > 5:
        min_x, max_x = min(all_lengths), max(all_lengths)
        pad_x = (max_x - min_x) * 0.15
        x_grid = np.linspace(max(0, min_x - pad_x), max_x + pad_x, 500)

        # 1. Plot Individual Sources (Colourful Fills)
        for d in data_map:
            data = d["stats"]["retained_lengths"]
            if len(data) > 2 and np.std(data) > 0:
                try:
                    kde = gaussian_kde(data)
                    y_grid = kde(x_grid)
                    # Increased alpha for better visibility
                    ax3.fill_between(x_grid, y_grid, color=d["colour"], alpha=0.15)
                    ax3.plot(x_grid, y_grid, color=d["colour"], linestyle="--", linewidth=1.5, alpha=0.9, label=f"{d['label']} (n={len(data)})")
                except Exception as e:
                    logger.debug(f"Skipped KDE plot for {d['label']} due to math error: {e}")

        # 2. Plot Combined Final (Robust Grey Fill)
        if len(all_lengths) > 2:
            try:
                kde_total = gaussian_kde(all_lengths)
                y_total = kde_total(x_grid)
                # Stronger grey fill
                ax3.fill_between(x_grid, y_total, color=PALETTE["grey_fill"], alpha=0.3)
                ax3.plot(x_grid, y_total, color=PALETTE["Total"], linewidth=2.5, label=f"Combined Final (N={len(all_lengths)})")

                mean_val = np.mean(all_lengths)
                ax3.axvline(mean_val, color=PALETTE["Total"], linestyle="-", linewidth=0.8, alpha=0.6)
                ax3.text(mean_val, max(y_total)*1.02, f"Mean: {mean_val:.1f}", ha="center", fontsize=CFG.VIS_FONT_TICK, color=PALETTE["Total"])

                # Comprehensive Statistics Box
                stats_text = (
                    f"DATASET STATISTICS\n"
                    f"------------------\n"
                    f"Total N  : {len(all_lengths):,}\n"
                    f"Mean     : {mean_val:.2f}\n"
                    f"Median   : {np.median(all_lengths):.1f}\n"
                    f"Range    : {min_x}-{max_x}\n"
                    f"Std Dev  : {np.std(all_lengths):.2f}\n"
                    f"Skewness : {skew(all_lengths):.2f}"
                )
                ax3.text(0.98, 0.95, stats_text, transform=ax3.transAxes, va="top", ha="right",
                         fontsize=CFG.VIS_FONT_AXIS_LABEL, fontfamily="monospace",
                         bbox=dict(facecolor="white", edgecolor=PALETTE["box_edge"], boxstyle="round,pad=0.6", alpha=0.95))

            except Exception as e:
                logger.debug(f"Skipped KDE plot for combined data due to math error: {e}")

        ax3.legend(loc="upper right", bbox_to_anchor=(0.82, 1.0), frameon=False, fontsize=CFG.VIS_FONT_AXIS_LABEL)
    else:
        ax3.text(0.5, 0.5, "Insufficient data for Density Plot", ha="center", transform=ax3.transAxes)

    apply_clean_spines(ax3)
    ax3.spines["left"].set_visible(False)
    ax3.set_yticks([])

    # Save
    plt.savefig(output_path.with_suffix(".png"), dpi=CFG.VIS_FIGURE_DPI, bbox_inches="tight")
    plt.close()
    return output_path.with_suffix(".png")


# =============================================================================
# SECTION 6: MAIN EXECUTION
# =============================================================================

def main():
    # -------------------------------------------------------------------------------
    # Step 6.1: Parse Arguments & Setup
    # -------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(description="Production Grade FASTA Merge (2 Files)")
    parser.add_argument("--master", required=True, help="File 1 (MASTER): Trusted.")
    parser.add_argument("--secondary", required=True, help="File 2 (Secondary): Deduped & Filtered.")
    parser.add_argument("--output", default="C_INP_Merged_for_Boltz-2.fasta", help="Output filename")
    parser.add_argument("--min-len", type=int, default=CFG.MERGE_SECONDARY_LEN_MIN, help="Min length for Secondary file")
    parser.add_argument("--max-len", type=int, default=CFG.MERGE_SECONDARY_LEN_MAX, help="Max length for Secondary file")
    parser.add_argument("--keep-gaps", action="store_true", help="Preserves '-' in output.")

    args = parser.parse_args()

    f1_path = Path(args.master)
    f2_path = Path(args.secondary)
    out_path = Path(args.output)

    logger = setup_logger(out_path)

    if not all(p.exists() for p in [f1_path, f2_path]):
        logger.error("Error: One or more input files do not exist.")
        sys.exit(1)

    _utils_mod.print_script_banner(
        "01_Merge_FAcDs.py",
        "FASTA Deduplication & Merge  ·  Length Filtering  ·  Sequence Standardisation",
    )
    logger.info(f"  Master    : {f1_path.name}")
    logger.info(f"  Secondary : {f2_path.name}  (Len: {args.min_len}–{args.max_len} aa)")
    logger.info(f"  Output    : {out_path.name}")

    seen_sequences = set()
    seen_ids = set()

    with open(out_path, "w") as out_handle:

        # -------------------------------------------------------------------------------
        # Step 6.2: Process Master File (Highest Priority)
        # -------------------------------------------------------------------------------
        # Rules: No filtering, adds to 'seen' database first.
        stats_master = process_and_write(
            f1_path, out_handle, "Master",
            seen_sequences, seen_ids, logger,
            keep_gaps=args.keep_gaps
        )

        # -------------------------------------------------------------------------------
        # Step 6.3: Process Secondary File
        # -------------------------------------------------------------------------------
        # Rules: Length filtering + Deduplicates against Master.
        stats_secondary = process_and_write(
            f2_path, out_handle, "Secondary",
            seen_sequences, seen_ids, logger,
            min_len=args.min_len, max_len=args.max_len,
            filter_length=True, keep_gaps=args.keep_gaps
        )

    # -------------------------------------------------------------------------------
    # Step 6.4: Generate Visual Report (QC Dashboard)
    # -------------------------------------------------------------------------------
    logger.info("\nGenerating Visual Report...")
    try:
        plot_path = generate_plots(stats_master, stats_secondary, out_path, logger)
    except Exception as e:
        logger.error(f"Visualisation failed: {e}")
        traceback.print_exc()
        plot_path = "FAILED"

    # -------------------------------------------------------------------------------
    # Step 6.5: Final Reporting & Shutdown
    # -------------------------------------------------------------------------------
    logger.info(_utils_mod.SEPARATOR_HEAVY)
    logger.info("FINAL REPORT")
    logger.info(_utils_mod.SEPARATOR_HEAVY)

    def log_stage(name, s):
        logger.info(f"\n[{name}]")
        logger.info(f"  Input Sequences       : {s['total']}")
        logger.info(f"  Dropped (Quality)     : {s['quality_fail']}")
        if s["length_fail"] > 0:
            logger.info(f"  Dropped (Length)      : {s['length_fail']}")
        logger.info(f"  Dropped (Duplicate)   : {s['dupes']}")
        logger.info(f"  RETAINED              : {s['kept']}")

    log_stage("File 1: Master", stats_master)
    log_stage("File 2: Secondary", stats_secondary)

    # Aggregate stats
    total_kept = len(seen_ids)
    all_lens = stats_master["retained_lengths"] + stats_secondary["retained_lengths"]

    final_mean = round(sum(all_lens) / total_kept, 2) if total_kept > 0 else 0
    min_l = min(all_lens) if all_lens else 0
    max_l = max(all_lens) if all_lens else 0

    logger.info("\n[FINAL DATASET]")
    logger.info(f"  Total Sequences       : {total_kept}")
    logger.info(f"  Unique Sequences      : {len(seen_sequences)}")
    _len_note = "  (master length-exempt)" if min_l < args.min_len else ""
    logger.info(f"  Length Range          : {min_l} - {max_l}{_len_note}")
    logger.info(f"  Mean Length           : {final_mean}")
    logger.info("\nFiles Saved:")
    logger.info(f"  1. FASTA : {out_path.resolve()}")
    logger.info(f"  2. LOG   : {(out_path.parent / '00_Merge.log').resolve()}")
    if plot_path != "FAILED":
        logger.info(f"  3. PLOTS : {plot_path.resolve()}")

    logger.info(_utils_mod.SEPARATOR_DASH)
    logger.info("✔ All Files Saved Successfully")

    logger.info(_utils_mod.SEPARATOR_LIGHT)

if __name__ == "__main__":
    _t0 = _time.perf_counter()
    main()
    _utils_mod.print_elapsed(_t0, "01_Merge_FAcDs.py")
