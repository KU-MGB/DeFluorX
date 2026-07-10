"""
===============================================================================
FAcDs Pipeline  |  MODULE 00_02  |  Shared Utilities
===============================================================================
Canonical source for console styling, logging infrastructure, matplotlib
spine helpers, MIC vector arithmetic, and geometric angle/dihedral functions.
All downstream scripts import from here — never duplicate these definitions.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 05 July 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Module        : 00_02_Project_Utils_FAcDs.py
  Role          : Shared utility library; no executable entry point.
  Imported by   : 01_Merge_FAcDs.py, 02_Production_FAcDs.py, 03_Validation_Figures_FAcDs.py,
                  04_Dendrogram_FAcDs.py, 05_TopN_and_PDB_Preparation_FAcDs.py,
                  05_TopN_and_PDB_Preparation_FAcDs.py, 06_SID_Prime-MMGBSA_FAcDs.py,
                  07_MD_Thermodynamics_QMMM_Engine_FAcDs.py
                  (also referenced by 00_01 for a ConsoleColours drift check)
  Reads         : (none — pure utility module)
  Writes        : (none — pure utility module)
───────────────────────────────────────────────────────────────────────────────

-------------------------------------------------------------------------------
Scientific References:
    1. Numerical arrays & vector geometry:
       - Harris, C.R. et al. (2020) Array programming with NumPy. Nature 585:357–362.
       - DOI: https://doi.org/10.1038/s41586-020-2649-2
    2. Plotting helpers (matplotlib spine/style utilities):
       - Hunter, J.D. (2007) Matplotlib. Comput Sci Eng 9:90–95.
       - DOI: https://doi.org/10.1109/MCSE.2007.55
    3. Structure / data / cheminformatics helpers (lazily imported on demand):
       - Gemmi: Wojdyr, M. (2022) J Open Source Softw 7:4200. DOI: https://doi.org/10.21105/joss.04200
       - pandas: McKinney, W. (2010) Data Structures for Statistical Computing in Python. Proc 9th Python in Science Conf 56–61. DOI: https://doi.org/10.25080/Majora-92bf1922-00a
       - RDKit: Landrum, G. (2006) RDKit: Open-source cheminformatics. https://www.rdkit.org
    Note: the angle/dihedral helpers implement standard vector geometry; metric
    definitions and their primary literature live in 00_01_Project_Config_FAcDs.py.
-------------------------------------------------------------------------------
"""

from __future__ import annotations

import logging
import os
import re as _re
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np

_ANSI_ESCAPE_RE = _re.compile(r"\033\[[0-9;]*[mKABCDEFGHJKSTfhilmnprsu]")


# =============================================================================
# SECTION 1: CONSOLE STYLING
# =============================================================================

class ConsoleColours:
    """ANSI terminal colour codes for pipeline console output."""
    OKGREEN = "\033[92m"   # green   — success / pass
    WARNING = "\033[93m"   # yellow  — caution
    FAIL    = "\033[91m"   # red     — error / fail
    OKBLUE  = "\033[94m"   # blue    — information
    MAGENTA = "\033[95m"   # magenta — script banners
    BOLD    = "\033[1m"    # bold    — section headers
    ENDC    = "\033[0m"    # reset   — end all formatting


# Horizontal separators — choose the weight that matches visual hierarchy.
SEPARATOR_HEAVY = "═" * 80   # major section boundary  (═══)
SEPARATOR_LIGHT = "─" * 80   # step / subsection       (───)
SEPARATOR_DASH  = "-" * 80   # info line / minor break  (---)


def safe_name(s: str) -> str:
    """Sanitises strings for secure usage as filenames, thereby preventing path-injection vulnerabilities."""
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)[:200]


# =============================================================================
# SECTION 2: LOGGING INFRASTRUCTURE
# =============================================================================

def setup_logging(
    log_path: Path | str,
    logger_name: str = "pfas_pipeline",
    mode: str = "a",
    timestamp: bool = True,
) -> logging.Logger:
    """
    Canonical file logger for all pipeline scripts.

    Creates a file-only handler; console output is handled by console_* below.
    Returns the Logger instance so each script can store it as a module global
    and pass it into console_* calls.

    Parameters
    ----------
    log_path    : Destination log file (parent directory is created if absent).
    logger_name : Unique name for this logger (prevents handler bleed between
                  scripts when both are imported in the same Python session).
    mode        : 'a' to append, 'w' to overwrite on each run.
    timestamp   : True  → "%(asctime)s | %(levelname)s | %(message)s"
                  False → "%(message)s"  (plain, for delivery manifests)
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log = logging.getLogger(logger_name)
    log.setLevel(logging.INFO)
    if log.handlers:
        log.handlers.clear()
    fmt = (
        "%(asctime)s | %(levelname)s | %(message)s"
        if timestamp
        else "%(message)s"
    )
    fh = logging.FileHandler(str(log_path), mode=mode, encoding="utf-8")
    fh.setFormatter(logging.Formatter(fmt, "%Y-%m-%d %H:%M:%S"))
    log.addHandler(fh)
    return log


# =============================================================================
# SECTION 3: CONSOLE OUTPUT FUNCTIONS
# =============================================================================
"""
All functions accept an optional `logger` parameter. Pass the module-level
logger from the calling script so output goes to both terminal and log file.
When `logger=None`, output is terminal-only (useful for standalone testing).
"""

def _strip_ansi(s: str) -> str:
    """Remove all ANSI/VT100 escape sequences from a string."""
    return _ANSI_ESCAPE_RE.sub("", s)


def console_title(msg: str, logger: logging.Logger | None = None) -> None:
    """Bold section header — printed and optionally written to log file."""
    print(f"\n{ConsoleColours.BOLD}{msg}{ConsoleColours.ENDC}", flush=True)
    if logger:
        logger.info(_strip_ansi(msg))


def console_info(msg: str, logger: logging.Logger | None = None) -> None:
    """Two-space-indented info line — printed and optionally written to log file."""
    # Ensure every line of a multi-line message is indented by two spaces.
    for line in str(msg).split("\n"):
        print(f"  {line}", flush=True)
    if logger:
        logger.info(_strip_ansi(msg))


def console_separator(
    logger: logging.Logger | None = None,
    heavy: bool = False,
) -> None:
    """Horizontal rule — printed and optionally written to log file.

    Parameters
    ----------
    logger : Logger instance to record output (optional).
    heavy  : True → SEPARATOR_HEAVY (═══), False → SEPARATOR_DASH (---).
    """
    sep = SEPARATOR_HEAVY if heavy else SEPARATOR_DASH
    print(sep, flush=True)
    if logger:
        logger.info(sep)


def print_script_banner(script_name: str, subtitle: str = "") -> None:
    """Print a coloured startup banner identifying which script is running.

    Parameters
    ----------
    script_name : Filename of the running script (e.g. "02_Production_FAcDs.py").
    subtitle    : One-line description shown below the script name.
    """
    _now = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"\n{SEPARATOR_HEAVY}", flush=True)
    print(
        f"  {ConsoleColours.MAGENTA}{ConsoleColours.BOLD}▶  {script_name}"
        f"{ConsoleColours.ENDC}  │  FAcDs Pipeline",
        flush=True,
    )
    if subtitle:
        print(f"  {subtitle}", flush=True)
    print(f"  Started : {_now}", flush=True)
    print(f"{SEPARATOR_HEAVY}\n", flush=True)


def print_elapsed(t0: float, script_name: str) -> None:
    """Print total pipeline-script elapsed time bookended by heavy separators.

    Parameters
    ----------
    t0          : ``time.perf_counter()`` value recorded before ``main()`` was called.
    script_name : Filename displayed in the output line (e.g. ``"02_Production_FAcDs.py"``).
    """
    _el = time.perf_counter() - t0
    _h, _rem = divmod(int(_el), 3600)
    _m, _s   = divmod(_rem, 60)
    _fmt = (f"{_h}h {_m:02d}m {_s:02d}s" if _h else
            f"{_m}m {_s:02d}s"            if _m else
            f"{_s}s")
    print(f"\n{SEPARATOR_HEAVY}", flush=True)
    print(
        f"  {ConsoleColours.OKGREEN}✔  {script_name}  —  Pipeline Phase Complete"
        f"  │  Total Elapsed: {_fmt}{ConsoleColours.ENDC}",
        flush=True,
    )
    print(f"{SEPARATOR_HEAVY}\n", flush=True)


class ReportManager:
    """Simultaneous console + file logger shared by the pipeline steps.

    The log path, banner header, section separator and rule width are supplied
    per step, and ``log_fn`` lets each caller route console output through its
    own logger-bound ``console_info`` so file logging is preserved. ``section``
    prints a bold title + separator and appends a divider to the log file.
    """
    def __init__(self, log_path, header, separator=SEPARATOR_LIGHT,
                 rule_width=80, log_fn=None):
        from pathlib import Path as _Path
        from datetime import datetime as _dt
        self.path = _Path(log_path)
        self.separator = separator
        self._log_fn = log_fn if log_fn is not None else console_info
        with open(self.path, "w") as f:
            f.write(header + "\n")
            f.write(f"Date: {_dt.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("=" * rule_width + "\n\n")

    def log(self, text: str):
        self._log_fn(text)
        with open(self.path, "a") as f:
            f.write(f"[LOG] {text}\n")

    def section(self, title: str):
        print(f"\n{ConsoleColours.BOLD}{title}{ConsoleColours.ENDC}", flush=True)
        print(self.separator, flush=True)
        with open(self.path, "a") as f:
            f.write(f"\n--- {title} ---\n")


# =============================================================================
# SECTION 4: MATPLOTLIB UTILITIES
# =============================================================================

def clean_spines(ax) -> None:
    """
    Academic-style axes: remove top/right spines, thin the remaining borders.

    Canonical replacement for any ``apply_clean_spines`` defined locally in
    individual scripts — import and call this function instead.
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(direction="out", length=4, width=0.8, labelsize=10)
    ax.grid(axis="x", linestyle="--", alpha=0.3)


# --- Central Ramachandran Plotting Helpers ------------------------------------

def _rama_get_atom_pos(res, name: str):
    """Find atom coordinates in a Gemmi residue structure."""
    atom = res.find_atom(name, "*")
    return atom.pos if atom else None


def compute_ramachandran_angles(st) -> list[tuple[str, int, float, float]]:
    """Extract (resname, resnum, phi, psi) for every residue that has both angles.

    Backbone dihedrals are computed only across genuine peptide bonds: the i-1→i and
    i→i+1 neighbours must be sequential in seqid AND covalently bonded (C–N ≤ 1.5 Å),
    so a modelled chain break, gap or insertion does not emit a spurious dihedral. All
    polymer chains in the model are processed.
    """
    import math
    import gemmi
    _PEPTIDE_CN_MAX = 1.5   # Å  C(i-1)–N(i) upper bound for a real peptide bond
    angles = []
    for chain in st[0]:
        residues = [r for r in chain if r.entity_type != gemmi.EntityType.NonPolymer
                    and r.entity_type != gemmi.EntityType.Water]
        for i, res in enumerate(residues):
            N  = _rama_get_atom_pos(res, "N")
            CA = _rama_get_atom_pos(res, "CA")
            C  = _rama_get_atom_pos(res, "C")
            if not (N and CA and C):
                continue
            phi = psi = None
            if i > 0:
                prev = residues[i - 1]
                C_prev = _rama_get_atom_pos(prev, "C")
                if C_prev and (int(res.seqid.num) - int(prev.seqid.num) == 1) \
                        and C_prev.dist(N) <= _PEPTIDE_CN_MAX:
                    phi = math.degrees(gemmi.calculate_dihedral(C_prev, N, CA, C))
            if i < len(residues) - 1:
                nxt = residues[i + 1]
                N_next = _rama_get_atom_pos(nxt, "N")
                if N_next and (int(nxt.seqid.num) - int(res.seqid.num) == 1) \
                        and C.dist(N_next) <= _PEPTIDE_CN_MAX:
                    psi = math.degrees(gemmi.calculate_dihedral(N, CA, C, N_next))
            if phi is not None and psi is not None:
                angles.append((res.name, int(res.seqid.num), phi, psi))
    return angles


def _rama_classify(phi: float, psi: float) -> str:
    """Classify a phi/psi pair as Favored, Allowed, or Outlier."""
    phi = ((phi + 180) % 360) - 180
    psi = ((psi + 180) % 360) - 180
    favored = (
        (-165 <= phi <= -30 and -70  <= psi <=  50) or
        (-180 <= phi <= -50 and (110 <= psi <= 180 or -180 <= psi <= -155)) or
        (  30 <= phi <=  90 and -25  <= psi <=  80)
    )
    if favored: return "Favored"
    allowed = (
        (-180 <= phi <=   0 and -100 <= psi <=  80) or
        (-180 <= phi <= -30 and   80 <= psi <= 180) or
        (-180 <= phi <= -30 and -180 <= psi <= -130) or
        (   0 <= phi <= 130 and  -50 <= psi <= 100)
    )
    return "Allowed" if allowed else "Outlier"


def _rama_stats(angles: list[tuple]) -> dict[str, Any]:
    """Calculate percentages of residues in favored, allowed, and outlier regions."""
    total = len(angles)
    counts = {"Favored": 0, "Allowed": 0, "Outlier": 0}
    for _, _, phi, psi in angles:
        counts[_rama_classify(phi, psi)] += 1
    pct = {k: (v / total * 100 if total else 0.0) for k, v in counts.items()}
    return {"total": total, "counts": counts, "pct": pct}


def _draw_rama_background(ax) -> None:
    """Draw the standard alpha/beta/L region backgrounds on a Ramachandran axes."""
    import matplotlib.pyplot as plt
    fav_c = "#dcedc8" # light green
    all_c = "#fff9c4" # light yellow

    alpha_fav  = plt.Polygon([(-165,-70),(-30,-70),(-30,50),(-165,50)], closed=True, fc=fav_c, ec="none", zorder=0)
    beta_fav1  = plt.Polygon([(-180,110),(-50,110),(-50,180),(-180,180)], closed=True, fc=fav_c, ec="none", zorder=0)
    beta_fav2  = plt.Polygon([(-180,-180),(-50,-180),(-50,-155),(-180,-155)], closed=True, fc=fav_c, ec="none", zorder=0)
    lhand_fav  = plt.Polygon([(30,-25),(90,-25),(90,80),(30,80)], closed=True, fc=fav_c, ec="none", zorder=0)

    allowed1   = plt.Polygon([(-180,-100),(0,-100),(0,80),(-180,80)], closed=True, fc=all_c, ec="none", zorder=0)
    allowed2   = plt.Polygon([(-180,80),(-30,80),(-30,180),(-180,180)], closed=True, fc=all_c, ec="none", zorder=0)
    allowed3   = plt.Polygon([(-180,-180),(-30,-180),(-30,-130),(-180,-130)], closed=True, fc=all_c, ec="none", zorder=0)
    allowed4   = plt.Polygon([(0,-50),(130,-50),(130,100),(0,100)], closed=True, fc=all_c, ec="none", zorder=0)

    for patch in [allowed1, allowed2, allowed3, allowed4]:
        ax.add_patch(patch)
    for patch in [alpha_fav, beta_fav1, beta_fav2, lhand_fav]:
        ax.add_patch(patch)

    ax.axhline(0, color="#bdbdbd", lw=1, zorder=1)
    ax.axvline(0, color="#bdbdbd", lw=1, zorder=1)


def save_ramachandran_comparison(angles_ref: list[tuple], angles_con: list[tuple],
                                 label_ref: str, label_con: str, out_path: Path | str,
                                 critical_res: dict[int, tuple[str, str]] = None,
                                 dpi: int = 300) -> None:
    """Save a side-by-side comparison Ramachandran PNG."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, angles, label, _colour in [
        (axes[0], angles_ref, label_ref, "#6a1b9a"),
        (axes[1], angles_con, label_con, "#1565c0"),
    ]:
        _draw_rama_background(ax)
        if angles:
            for marker_type, m_shape in [("General", "o"), ("Glycine", "^"), ("Proline", "s")]:
                m_phis, m_psis, m_cols = [], [], []
                for resname, resnum, phi, psi in angles:
                    if (marker_type == "Glycine" and resname == "GLY") or \
                       (marker_type == "Proline" and resname == "PRO") or \
                       (marker_type == "General" and resname not in ["GLY", "PRO"]):
                        m_phis.append(phi)
                        m_psis.append(psi)
                        classification = _rama_classify(phi, psi)
                        if classification == "Favored":
                            m_cols.append("#2e7d32")
                        elif classification == "Allowed":
                            m_cols.append("#f57f17")
                        else:
                            m_cols.append("#c62828")
                if m_phis:
                    ax.scatter(m_phis, m_psis, c=m_cols, marker=m_shape, s=20,
                               alpha=0.85, linewidths=0.4, edgecolors="white", zorder=3)

            if critical_res:
                triad_styles = {
                    "Acid": {"marker": "D", "fc": "#d81b60", "ec": "#880e4f", "s": 80},
                    "Base": {"marker": "p", "fc": "#1e88e5", "ec": "#0d47a1", "s": 100},
                    "Nuc":  {"marker": "*", "fc": "#00bcd4", "ec": "#006064", "s": 180}
                }
                for resname, resnum, phi, psi in angles:
                    if resnum in critical_res and critical_res[resnum][0] == resname:
                        role = critical_res[resnum][1]
                        style = triad_styles.get(role, {"marker": "X", "fc": "#9c27b0", "ec": "#4a148c", "s": 100})
                        ax.scatter([phi], [psi], c=style["fc"], marker=style["marker"], s=style["s"],
                                   alpha=1.0, linewidths=0.8, edgecolors=style["ec"], zorder=5)

        st = _rama_stats(angles)
        fav_pct = f"({st['pct']['Favored']:.1f}%)"
        allowed_pct = f"({st['pct']['Allowed']:.1f}%)"
        out_pct = f"({st['pct']['Outlier']:.1f}%)"
        tot_pct = "(100.0%)"

        txt = (f"Favoured  {st['counts']['Favored']:<4} {fav_pct:<8}\n"
               f"Allowed   {st['counts']['Allowed']:<4} {allowed_pct:<8}\n"
               f"Outlier   {st['counts']['Outlier']:<4} {out_pct:<8}\n"
               f"Total     {st['total']:<4} {tot_pct:<8}")
        ax.text(0.97, 0.97, txt, transform=ax.transAxes, fontsize=9,
                va="top", ha="right", multialignment="left", family="monospace",
                bbox=dict(fc="#ffffff", alpha=0.10, ec="#bdbdbd", boxstyle="round,pad=0.4"))

        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_aspect("equal")
        ax.set_xlabel("φ (phi) °", fontsize=11, fontweight="500")
        ax.set_ylabel("ψ (psi) °", fontsize=11, fontweight="500")
        ax.set_title(label, fontsize=12, fontweight="bold", pad=10)
        ax.set_xticks(range(-180, 181, 60))
        ax.set_yticks(range(-180, 181, 60))
        ax.tick_params(labelsize=9)
        ax.grid(True, linestyle=":", alpha=0.6, color="#9e9e9e", zorder=1)

    favoured_p = mpatches.Patch(color="#2e7d32", label="Favoured")
    allowed_p  = mpatches.Patch(color="#f57f17", label="Allowed")
    outlier_p  = mpatches.Patch(color="#c62828", label="Outlier")

    gen_m = mlines.Line2D([], [], color="none", marker="o", markerfacecolor="gray", markeredgecolor="white", markersize=7, label="General")
    gly_m = mlines.Line2D([], [], marker="^", color="none", markerfacecolor="gray", markeredgecolor="white", markersize=7, label="Glycine")
    pro_m = mlines.Line2D([], [], marker="s", color="none", markerfacecolor="gray", markeredgecolor="white", markersize=7, label="Proline")
    crit_handles = []
    if critical_res:
        triad_styles = {
            "Acid": {"marker": "D", "fc": "#d81b60", "ec": "#880e4f"},
            "Base": {"marker": "p", "fc": "#1e88e5", "ec": "#0d47a1"},
            "Nuc":  {"marker": "*", "fc": "#00bcd4", "ec": "#006064"}
        }
        for resnum, (resname, role) in sorted(critical_res.items(), key=lambda x: x[1][1]):
            style = triad_styles.get(role, {"marker": "X", "fc": "#9c27b0", "ec": "#4a148c"})
            handle = mlines.Line2D([], [], color="none", marker=style["marker"],
                                   markerfacecolor=style["fc"], markeredgecolor=style["ec"],
                                   markersize=9, label=f"{role}: {resname}{resnum}")
            crit_handles.append(handle)

    all_handles = [favoured_p, allowed_p, outlier_p, gen_m, gly_m, pro_m] + crit_handles

    fig.legend(handles=all_handles,
               loc="upper center", ncol=len(all_handles), fontsize=9, frameon=True,
               bbox_to_anchor=(0.5, -0.005), columnspacing=0.8, handletextpad=0.4)

    fig.suptitle("Ramachandran Comparison", fontsize=15, fontweight="bold", y=1.05)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def save_ramachandran_plot(angles: list[tuple], title: str, out_path: Path | str, dpi: int = 300) -> None:
    """Save a single-structure Ramachandran plot PNG."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 6))
    _draw_rama_background(ax)
    if angles:
        phis = [a[2] for a in angles]
        psis = [a[3] for a in angles]
        _rama_cls = [_rama_classify(p, s) for p, s in zip(phis, psis)]
        colours = ["#1b5e20" if c == "Favored"
                   else ("#f9a825" if c == "Allowed" else "#c62828")
                   for c in _rama_cls]
        ax.scatter(phis, psis, c=colours, s=14, alpha=0.75, linewidths=0, zorder=3)
    stats = _rama_stats(angles)
    legend_txt = (f"Favored  {stats['pct']['Favored']:.1f}%  ({stats['counts']['Favored']})\n"
                  f"Allowed  {stats['pct']['Allowed']:.1f}%  ({stats['counts']['Allowed']})\n"
                  f"Outlier   {stats['pct']['Outlier']:.1f}%  ({stats['counts']['Outlier']})\n"
                  f"Total: {stats['total']} residues")
    ax.text(0.98, 0.98, legend_txt, transform=ax.transAxes, fontsize=8,
            va="top", ha="right", family="monospace",
            bbox=dict(fc="white", alpha=0.7, ec="#cccccc", boxstyle="round,pad=0.3"))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("φ (phi) °", fontsize=11)
    ax.set_ylabel("ψ (psi) °", fontsize=11)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


# =============================================================================
# SECTION 5: GEOMETRIC & BIO-MATHEMATICAL UTILITIES
# =============================================================================

# -------------------------------------------------------------------------------
# Step 5.1: Minimum Image Convention (MIC) Vector
# -------------------------------------------------------------------------------

def get_mic_vector(pos1, pos2, box):
    """
    MIC displacement vector from pos2 → pos1, supporting orthorhombic and
    triclinic periodic boxes.

    Parameters
    ----------
    pos1 : numpy array of Cartesian coordinates (Å).
    pos2 : numpy array of Cartesian coordinates (Å).
    box  : 3×3 box matrix (triclinic) or None (vacuum / already-unwrapped).
    """
    # Coerce inputs to numpy arrays if they are gemmi.Position or list/tuple
    if hasattr(pos1, "x") and hasattr(pos1, "y") and hasattr(pos1, "z"):
        pos1 = np.array([pos1.x, pos1.y, pos1.z], dtype=float)
    else:
        pos1 = np.asarray(pos1, dtype=float)

    if hasattr(pos2, "x") and hasattr(pos2, "y") and hasattr(pos2, "z"):
        pos2 = np.array([pos2.x, pos2.y, pos2.z], dtype=float)
    else:
        pos2 = np.asarray(pos2, dtype=float)

    vec = pos1 - pos2
    if box is not None:
        try:
            b3       = _ensure_box_3x3(box)
            inv_box  = np.linalg.inv(b3)
            frac_vec = np.dot(vec, inv_box)
            frac_vec = frac_vec - np.round(frac_vec)
            vec      = np.dot(frac_vec, b3)
        except (np.linalg.LinAlgError, AttributeError, ValueError):
            diag = _box_diag(box)
            vec  = vec - diag * np.round(vec / diag)
    return vec


# -------------------------------------------------------------------------------
# Step 5.2: Distance
# -------------------------------------------------------------------------------

def distance(pos1, pos2, box=None) -> float:
    """PBC-corrected or Euclidean distance between two Cartesian coordinates (Å)."""
    return float(np.linalg.norm(get_mic_vector(np.asarray(pos1), np.asarray(pos2), box)))


def calculate_min_distance(
    frame_pos_a: np.ndarray,
    frame_pos_b: np.ndarray,
    box=None,
) -> float:
    """Minimum PBC-corrected distance between two sets of Cartesian positions."""
    return float(mic_dists_2d(np.asarray(frame_pos_a), np.asarray(frame_pos_b), box).min())


def _ensure_box_3x3(box) -> np.ndarray:
    """
    Coerce a periodic box to a (3,3) float64 matrix.
    Handles: (3,3) arrays, flat (9,) arrays (Schrödinger frame.box format),
    and (3,) diagonal vectors.
    """
    b = np.asarray(box, dtype=float)
    if b.shape == (3, 3):
        return b
    if b.size == 9:
        return b.reshape(3, 3)
    if b.size == 3:
        return np.diag(b)
    raise ValueError(f"Unrecognised box shape: {b.shape}")


def _box_diag(box) -> np.ndarray:
    """Orthorhombic fallback box lengths as a safe (3,) diagonal.

    Used when the full triclinic minimum-image transform fails (singular or
    non-invertible cell). Extracts the diagonal from a (3,3) / flat (9,) matrix
    or a (3,) vector, and replaces any zero length with 1e-6 to avoid a
    divide-by-zero in the orthorhombic wrap.
    """
    b = np.asarray(box, dtype=float)
    diag = (np.diag(b) if b.ndim == 2
            else b if b.size == 3
            else np.array([b.flat[0], b.flat[4], b.flat[8]]))
    return np.where(diag == 0, 1e-6, diag)


def mic_dists_2d(pos_a: np.ndarray, pos_b: np.ndarray, box) -> np.ndarray:
    """
    PBC-corrected pairwise distance matrix.
    pos_a (M, 3), pos_b (N, 3) → (M, N) float64 distance matrix.
    Box inverse is computed once for the batch (vs. once per pair in get_mic_vector).
    """
    vecs = pos_a[:, None, :] - pos_b[None, :, :]   # (M, N, 3)
    if box is not None:
        try:
            b3    = _ensure_box_3x3(box)
            inv_b = np.linalg.inv(b3)
            mn    = vecs.shape[0] * vecs.shape[1]
            frac  = vecs.reshape(mn, 3) @ inv_b
            frac -= np.round(frac)
            vecs  = (frac @ b3).reshape(vecs.shape)
        except (np.linalg.LinAlgError, AttributeError, ValueError):
            diag = _box_diag(box)
            vecs = vecs - diag * np.round(vecs / diag)
    return np.linalg.norm(vecs, axis=2)


# -------------------------------------------------------------------------------
# Step 5.3: Angle
# -------------------------------------------------------------------------------

def calculate_angle(p1, p2, p3, box=None) -> float:
    """
    Angle in degrees at vertex p2 (p1–p2–p3), PBC-corrected.

    Returns 0.0 if either leg has zero length.
    """
    v1 = get_mic_vector(p1, p2, box)
    v2 = get_mic_vector(p3, p2, box)
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    if n1 < 1e-6 or n2 < 1e-6:
        return 0.0
    cos_theta = np.dot(v1, v2) / (n1 * n2)
    return np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))


# -------------------------------------------------------------------------------
# Step 5.4: Dihedral
# -------------------------------------------------------------------------------

def calculate_dihedral(p1, p2, p3, p4, box=None) -> float:
    """
    Proper dihedral angle in degrees for p1–p2–p3–p4, PBC-corrected.
    Uses the Gram–Schmidt projection (numerically stable for near-linear bonds).
    """
    b0 = -get_mic_vector(p2, p1, box)
    b1 =  get_mic_vector(p3, p2, box)
    b2 =  get_mic_vector(p4, p3, box)
    b1_len = np.linalg.norm(b1)
    if b1_len < 1e-6:          # coincident central atoms (p2 ≡ p3) — dihedral undefined
        return 0.0
    b1 /= b1_len
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


# -------------------------------------------------------------------------------
# Step 5.5: Improper Dihedral (Walden Inversion / TS Flattening)
# -------------------------------------------------------------------------------

def calculate_improper_dihedral(p1, p2, p3, p4, box=None) -> float:
    """
    Out-of-plane improper dihedral at centre p1 — used to detect transition-state
    (Walden inversion) flattening of the sp3 electrophilic carbon.
    Ref (Walden inversion): Walden, P. (1896) Ber. Dtsch. Chem. Ges. 29:133–138.
    DOI: https://doi.org/10.1002/cber.18960290127

    Returns the angle between the p4 vector and the normal of the p2–p1–p3 plane
    (degrees from planarity).  Values near 0° indicate a near-planar TS geometry.

    Note: the reference plane is built from (p1, p2, p3), so for an asymmetric
    centre the returned angle has a mild dependence on which equatorial atom is
    passed as p4. The variation is small relative to the Walden TS gate
    (±WALDEN_IMPROPER_MAX) and does not change the planar/pyramidal verdict; it
    is reported as a planarity heuristic, not an exact symmetric out-of-plane
    distance.
    """
    v1 = get_mic_vector(p2, p1, box)
    v2 = get_mic_vector(p3, p1, box)
    v3 = get_mic_vector(p4, p1, box)
    n = np.cross(v1, v2)
    n_len = np.linalg.norm(n)
    if n_len < 1e-6:
        return 0.0
    n /= n_len
    v3_len = np.linalg.norm(v3)
    if v3_len < 1e-6:          # p4 coincides with the centre p1 — undefined out-of-plane angle
        return 0.0
    cos_theta = np.dot(v3, n) / v3_len
    return 90.0 - np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))


# -------------------------------------------------------------------------------
# Step 5.6: Bürgi–Dunitz Angle
# -------------------------------------------------------------------------------

def calculate_burgi_dunitz(nuc_pos, c_pos, o_pos) -> float:
    """
    Bürgi–Dunitz angle (Nucleophile–Carbon–Oxygen) in degrees; ideal ≈ 107° for
    nucleophilic addition at an sp² carbonyl carbon.

    AUXILIARY / REFERENCE metric only — NOT a tier gate. In the FAcD pipeline the
    catalytic event is SN2 at the sp³ α-carbon (the backside O–C–F attack angle,
    ideal 180°, is the gating geometry). This BD angle is reported solely as a
    nucleophile-vs-substrate-carboxylate pre-organisation descriptor and MUST be
    computed on the substrate carbonyl carbon (the –COO⁻ carbon), never on the
    sp³ SN2 carbon. Computing Nu–C(sp³)–F and calling it Bürgi–Dunitz is incorrect.
    """
    return calculate_angle(nuc_pos, c_pos, o_pos)


# -------------------------------------------------------------------------------
# Step 5.6b: Oδ Orientation Fallback for Smart-Lock nucleophile search
# -------------------------------------------------------------------------------

def find_nucleophile_od_fallback(
    asp_candidates: list,
    ca_pos: "np.ndarray",
    f_pos: "np.ndarray",
    max_dist: float,
    min_angle: float,
) -> "tuple | None":
    """
    Geometry-based Oδ fallback nucleophile search for the Smart-Lock algorithm.

    Scans ASP/ASH candidates for any Oδ (OD1 or OD2) within `max_dist` Å of
    substrate Cα whose approach vector subtends ≥ `min_angle`° at Cα relative
    to the leaving-fluorine arm.  The candidate closest to Cα that satisfies
    both criteria is returned.

    Parameters
    ----------
    asp_candidates : list of (resnum, od1_coords_np, od2_coords_np)
        Each element is a 3-tuple: integer residue number plus two (3,) NumPy
        arrays for OD1 and OD2 coordinates.  Pass None for a missing oxygen.
    ca_pos : np.ndarray, shape (3,)
        Coordinates of the substrate electrophilic carbon (Cα / warhead C).
    f_pos : np.ndarray, shape (3,)
        Coordinates of the departing fluorine atom.
    max_dist : float
        Maximum Oδ–Cα distance in Å (e.g. CFG.SMART_LOCK_OD_FALLBACK_DIST).
    min_angle : float
        Minimum Oδ–Cα–F angle in degrees (e.g. CFG.SMART_LOCK_OD_FALLBACK_ANGLE).

    Returns
    -------
    tuple (resnum: int, od_coords: np.ndarray) or None
        The residue number and coordinates of the best qualifying Oδ, or None
        if no candidate satisfies the geometric criteria.
    """
    best: "tuple | None" = None
    best_dist = float("inf")
    vec_c_f   = np.asarray(f_pos, dtype=float) - np.asarray(ca_pos, dtype=float)
    norm_cf   = float(np.linalg.norm(vec_c_f))
    if norm_cf < 1e-6:
        return None
    for resnum, od1, od2 in asp_candidates:
        for od_pos in (od1, od2):
            if od_pos is None:
                continue
            od_arr   = np.asarray(od_pos, dtype=float)
            ca_arr   = np.asarray(ca_pos, dtype=float)
            vec_c_od = od_arr - ca_arr
            d        = float(np.linalg.norm(vec_c_od))
            if d > max_dist:
                continue
            norm_od = float(np.linalg.norm(vec_c_od))
            if norm_od < 1e-6:
                continue
            cos_a = float(np.dot(vec_c_od, vec_c_f) / (norm_od * norm_cf))
            ang   = float(np.degrees(np.arccos(np.clip(cos_a, -1.0, 1.0))))
            if ang >= min_angle and d < best_dist:
                best_dist = d
                best      = (resnum, od_arr)
    return best


# -------------------------------------------------------------------------------
# Step 5.7: Flippin–Lodge Angle
# -------------------------------------------------------------------------------

def calculate_flippin_lodge(nuc_pos, c_pos, r1_pos, r2_pos) -> float:
    """
    Flippin–Lodge torsional pre-alignment angle in degrees.
    Ref: Lodge, E. P. & Heathcock, C. H. (1987) J. Am. Chem. Soc. 109(11):3353–3361.
    DOI: https://doi.org/10.1021/ja00245a028

    Measures the angle between the Nu–C vector projected onto the R1–C–R2 plane
    and the bisector of R1–C–R2.  Ideal value = 0° (eclipsed, minimises steric
    clash between the incoming nucleophile and the R substituents).

    AUXILIARY / REFERENCE metric only — NOT a tier gate. Reports the in-plane
    (lateral) component of the nucleophile approach, complementing the backside
    O–C–F attack angle. Requires TWO well-defined heavy spectator substituents
    (R1, R2) on the attacked carbon; for small substrates whose α-carbon carries
    only hydrogens (often absent in heavy-atom CIF output) the offset is undefined
    — callers must guard and report NaN rather than a spurious value.

    Returns 999.0 on numerical failure (degenerate geometry).
    """
    try:
        def _arr(p):
            if hasattr(p, "x"):
                return np.array([float(p.x), float(p.y), float(p.z)])
            return np.array(p, dtype=float)

        nuc = _arr(nuc_pos)
        c   = _arr(c_pos)
        r1  = _arr(r1_pos)
        r2  = _arr(r2_pos)

        vec_nuc = nuc - c
        vec_r1  = r1  - c
        vec_r2  = r2  - c

        # Unit-normalise with a 1e-6 zero-length guard. Degenerate geometry
        # (collinear R1–C–R2 or coincident atoms) gives a zero-length vector; raising
        # here routes such cases to the 999.0 fallback rather than emitting NaN.
        def _unit(v):
            n = np.linalg.norm(v)
            if n < 1e-6:
                raise ValueError("degenerate vector (zero length)")
            return v / n

        normal        = _unit(np.cross(vec_r1, vec_r2))
        vec_nuc_proj  = _unit(vec_nuc - np.dot(vec_nuc, normal) * normal)
        bisector      = _unit(_unit(vec_r1) + _unit(vec_r2))

        return float(
            np.degrees(np.arccos(np.clip(np.dot(vec_nuc_proj, bisector), -1.0, 1.0)))
        )
    except Exception:
        return 999.0


# =============================================================================
# SECTION 6: DATA PIPELINE UTILITIES
# =============================================================================

def standardise_dataframe_tiers(df, cfg):
    """
    Apply canonical tier names, fill missing values, and ensure categorical
    ordering based on Project Config.

    Parameters
    ----------
    df  : pd.DataFrame containing a 'degrader_tier' column.
    cfg : CFG dataclass instance.
    """
    import pandas as pd
    col = cfg.COL_TIER
    if col not in df.columns:
        return df

    # 1. Fill missing/None with Decoy
    df[col] = df[col].fillna(cfg.TIER_DECOY).replace("None", cfg.TIER_DECOY)

    # 2. Enforce Categorical Type with Config Order
    df[col] = pd.Categorical(df[col], categories=cfg.TIER_ORDER, ordered=True)
    return df


def get_alignment_grade(identity_pct, cfg) -> str:
    """
    Resolve a letter grade (A-I) for a given sequence identity percentage
    using bins defined in Project Config.
    """
    # Left-open, right-closed intervals (lo, hi] to match the pd.cut(right=True)
    # grading used for the CSV Alignment_Grade column, so a boundary value such as
    # 90.0% receives the same letter from both code paths.
    for lo, hi, letter in cfg.ALIGN_GRADE_DEFS:
        if lo < identity_pct <= hi:
            return letter
    return "I"  # Fallback for ultra-low identity (incl. exactly 0.0%)


# =============================================================================
# SECTION 7: LIGAND PHYSICO-CHEMICAL PROPERTIES (RDKit)
# =============================================================================

def compute_ligand_properties(smiles_file) -> dict:
    """
    Compute per-ligand carbon count (nC), fluorine count (nF) and molecular
    weight (g/mol) directly from a SMILES panel via RDKit.

    The SMILES file is the single source — nothing is hardcoded, so swapping the
    ligand set (or proteins) requires no code change. Each line is
    "<SMILES>\\t<name>" (whitespace-separated); blank lines, comment lines (#)
    and entries RDKit cannot parse are skipped silently.

    Returns
    -------
    dict
        {ligand_name: {"nC": int, "nF": int, "mw": float}}. Empty dict if the
        file is absent or RDKit is unavailable (callers handle the empty case).
    """
    props: dict = {}
    p = Path(smiles_file)
    if not p.exists():
        return props
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors
    except Exception:
        return props
    for line in p.read_text(encoding="utf-8", errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        smiles, name = parts[0], parts[1]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        n_c = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "C")
        n_f = sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "F")
        props[name] = {"nC": n_c, "nF": n_f, "mw": round(Descriptors.MolWt(mol), 1)}
    return props


# =============================================================================
# SECTION: SCHRÖDINGER JOB-SERVER SCRATCH LOCATION
# =============================================================================
def _fs_device(path: Path) -> "int | None":
    """st_dev of the nearest existing ancestor of `path` (identifies its mounted
    filesystem), or None if unreadable."""
    try:
        p = Path(path).resolve()
        while not p.exists():
            p = p.parent
        return os.stat(p).st_dev
    except Exception:
        return None


def _jobserver_addresses(jsc: str) -> "list[str]":
    """Registered local job-server addresses (e.g. ['localhost:40931']) parsed from
    ``jsc server-info`` — needed to reload the hosts file on a RUNNING server."""
    try:
        out = subprocess.run([jsc, "server-info"], capture_output=True, text=True, timeout=60)
        return _re.findall(r"^\s*(\S+:\d+)\s*$", out.stdout, _re.MULTILINE)
    except Exception:
        return []


def ensure_jobserver_on_working_disk(
    schrodinger: str,
    work_disk_anchor: "Path | str",
    jobserver_subdir: str = "_Schrodinger_JobServer",
    echo=None,
) -> bool:
    """Force ALL Schrödinger job scratch onto the working disk — live, without
    stopping the server or killing a running job.

    The huge per-subjob scratch a Prime MM-GBSA / QSite run stages (the multi-GB
    complexes copies, hundreds of GB on a 100k-frame trajectory) lands in the job
    server's ``tmpdir``, NOT its server directory. That ``tmpdir`` is the
    ``localhost`` entry of ``$SCHRODINGER/schrodinger.hosts``; when it is ``/tmp``
    (the default) the OS disk fills and the job dies with a silent
    ``copy_file_range: no space left on device`` (rc=1) after hours.

    The fix is to point that ``tmpdir`` at the working disk and apply it with
    ``jsc admin reload-hosts <server-address>`` — which a RUNNING server accepts
    for its NEXT jobs, so nothing is stopped and no live job is lost. (This is why
    the earlier stop-and-relocate approach was wrong: it needed an idle server that
    never comes, and it moved the small server dir, not the scratch tmpdir.)

    Returns True when the hosts ``tmpdir`` is on the working-disk filesystem
    (already, or set this call) so scratch-heavy jobs are safe to launch; False
    only if the hosts file cannot be read/written.

    Parameters
    ----------
    schrodinger      : path to the Schrödinger installation ($SCHRODINGER).
    work_disk_anchor : any path on the target (large) working disk; scratch goes to
                       ``<anchor>/<jobserver_subdir>``.
    jobserver_subdir : name of the working-disk scratch directory (from CFG).
    echo             : optional print-like callable for progress lines.
    """
    say = echo or (lambda m="": print(f"  {m}", flush=True))
    jsc = os.path.join(str(schrodinger), "jsc")
    hosts = Path(schrodinger) / "schrodinger.hosts"
    anchor = Path(work_disk_anchor)
    target = (anchor / jobserver_subdir).resolve()
    if not hosts.is_file():
        say(f"[job-server] hosts file not found at {hosts}; cannot redirect scratch — "
            f"Schrödinger scratch may land on /tmp (OS disk).")
        return False
    try:
        text = hosts.read_text()
    except Exception as e:
        say(f"[job-server] could not read {hosts} ({e}); scratch may stay on /tmp.")
        return False

    _m = _re.search(r"(?m)^\s*tmpdir:\s*(.*)$", text)
    _cur = _m.group(1).strip() if _m else ""
    anchor_dev = _fs_device(anchor)
    if _cur and anchor_dev is not None and _fs_device(Path(_cur)) == anchor_dev:
        say(f"[job-server] scratch OK — hosts tmpdir {_cur} is on the working disk.")
        return True

    say(f"[job-server] hosts tmpdir is '{_cur or '(unset → /tmp)'}' (OS disk); a large "
        f"Prime/QSite job would exhaust it. Redirecting → {target}")
    try:
        target.mkdir(parents=True, exist_ok=True)
        # Back up the hosts file next to the run's code backups, then set the
        # localhost tmpdir (its settings apply to every host entry). Replace the
        # first tmpdir: line (localhost is first in the file); insert one after the
        # localhost name if none exists.
        try:
            (hosts.parent / f"schrodinger.hosts.bak").write_text(text)
        except Exception:
            pass
        if _m:
            new_text = _re.sub(r"(?m)^(\s*tmpdir:\s*).*$",
                               lambda _mm: f"{_mm.group(1)}{target}", text, count=1)
        else:
            new_text = _re.sub(r"(?m)^(\s*name:\s*localhost\s*)$",
                               lambda _mm: f"{_mm.group(1)}\ntmpdir:      {target}", text, count=1)
        hosts.write_text(new_text)
    except Exception as e:
        say(f"[job-server] could not write hosts tmpdir ({e}); scratch may stay on /tmp.")
        return False

    # Apply on the running server(s) so the change takes effect for the next jobs
    # without a restart (no live job is killed).
    _addrs = _jobserver_addresses(jsc) if os.path.isfile(jsc) else []
    for _a in _addrs:
        try:
            subprocess.run([jsc, "admin", "reload-hosts", _a],
                           capture_output=True, text=True, timeout=60)
        except Exception:
            pass
    say(f"[job-server] ✔ scratch tmpdir → {target} (applied to {len(_addrs) or 'no'} "
        f"running server(s) via reload-hosts; no job stopped).")
    return True


