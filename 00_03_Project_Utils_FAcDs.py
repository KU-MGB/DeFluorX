"""
===============================================================================
FAcDs Pipeline  |  MODULE 00_03  |  Shared Utilities
===============================================================================
Canonical source for console styling, logging infrastructure, matplotlib
spine helpers, MIC vector arithmetic, and geometric angle/dihedral functions.
All downstream scripts import from here — never duplicate these definitions.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 June 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Module        : 00_03_Project_Utils_FAcDs.py
  Role          : Shared utility library; no executable entry point.
  Imported by   : 01_Merge_FAcDs.py, 02_Production_FAcDs.py, 03_Validation_Figures_FAcDs.py,
                  04_Phylogeny_FAcDs.py, 05_CIF-PDB_Preparation_FAcDs.py,
                  06_Top-N_Extraction_FAcDs.py,
                  08_MD_Thermodynamics_QMMM_Engine_FAcDs.py
  Reads         : (none — pure utility module)
  Writes        : (none — pure utility module)
───────────────────────────────────────────────────────────────────────────────
"""

from __future__ import annotations

import logging
import re as _re
import sys as _sys
import time
from pathlib import Path
from typing import Any

import numpy as np

_ANSI_ESCAPE_RE = _re.compile(r'\033\[[0-9;]*[mKABCDEFGHJKSTfhilmnprsu]')


# ===============================================================================
# SECTION 1: CONSOLE STYLING
# ===============================================================================

class ConsoleColours:
    """ANSI terminal colour codes for pipeline console output."""
    OKGREEN = '\033[92m'   # green   — success / pass
    WARNING = '\033[93m'   # yellow  — caution
    FAIL    = '\033[91m'   # red     — error / fail
    OKBLUE  = '\033[94m'   # blue    — information
    MAGENTA = '\033[95m'   # magenta — script banners
    BOLD    = '\033[1m'    # bold    — section headers
    ENDC    = '\033[0m'    # reset   — end all formatting


# Horizontal separators — choose the weight that matches visual hierarchy.
SEPARATOR_HEAVY = "═" * 80   # major section boundary  (═══)
SEPARATOR_LIGHT = "─" * 80   # step / subsection       (───)
SEPARATOR_DASH  = "-" * 80   # info line / minor break  (---)


def safe_name(s: str) -> str:
    """Sanitises strings for secure usage as filenames, thereby preventing path-injection vulnerabilities."""
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in s)[:200]


# ===============================================================================
# SECTION 2: LOGGING INFRASTRUCTURE
# ===============================================================================

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


# ===============================================================================
# SECTION 3: CONSOLE OUTPUT FUNCTIONS
# ===============================================================================
# All functions accept an optional `logger` parameter.  Pass the module-level
# logger from the calling script so output goes to both terminal and log file.
# When `logger=None`, output is terminal-only (useful for standalone testing).

def _strip_ansi(s: str) -> str:
    """Remove all ANSI/VT100 escape sequences from a string."""
    return _ANSI_ESCAPE_RE.sub('', s)


def console_title(msg: str, logger: logging.Logger | None = None) -> None:
    """Bold section header — printed and optionally written to log file."""
    print(f"\n{ConsoleColours.BOLD}{msg}{ConsoleColours.ENDC}", flush=True)
    if logger:
        logger.info(_strip_ansi(msg))


def console_info(msg: str, logger: logging.Logger | None = None) -> None:
    """Two-space-indented info line — printed and optionally written to log file."""
    # Ensure every line of a multi-line message is indented by two spaces.
    for line in str(msg).split('\n'):
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


# ===============================================================================
# SECTION 4: MATPLOTLIB UTILITIES
# ===============================================================================

def clean_spines(ax) -> None:
    """
    Academic-style axes: remove top/right spines, thin the remaining borders.

    Canonical replacement for any ``apply_clean_spines`` defined locally in
    individual scripts — import and call this function instead.
    """
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.tick_params(direction='out', length=4, width=0.8, labelsize=10)
    ax.grid(axis='x', linestyle='--', alpha=0.3)


# --- Central Ramachandran Plotting Helpers ------------------------------------

def _rama_get_atom_pos(res, name: str):
    """Find atom coordinates in a Gemmi residue structure."""
    atom = res.find_atom(name, '*')
    return atom.pos if atom else None


def compute_ramachandran_angles(st) -> list[tuple[str, int, float, float]]:
    """Extract (resname, resnum, phi, psi) for every residue that has both angles."""
    import math
    import gemmi
    angles = []
    chain = st[0][0]
    residues = [r for r in chain if r.entity_type != gemmi.EntityType.NonPolymer
                and r.entity_type != gemmi.EntityType.Water]
    for i, res in enumerate(residues):
        N  = _rama_get_atom_pos(res, 'N')
        CA = _rama_get_atom_pos(res, 'CA')
        C  = _rama_get_atom_pos(res, 'C')
        if not (N and CA and C):
            continue
        phi = psi = None
        if i > 0:
            C_prev = _rama_get_atom_pos(residues[i - 1], 'C')
            if C_prev:
                phi = math.degrees(gemmi.calculate_dihedral(C_prev, N, CA, C))
        if i < len(residues) - 1:
            N_next = _rama_get_atom_pos(residues[i + 1], 'N')
            if N_next:
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
    if favored: return 'Favored'
    allowed = (
        (-180 <= phi <=   0 and -100 <= psi <=  80) or
        (-180 <= phi <= -30 and   80 <= psi <= 180) or
        (-180 <= phi <= -30 and -180 <= psi <= -130) or
        (   0 <= phi <= 130 and  -50 <= psi <= 100)
    )
    return 'Allowed' if allowed else 'Outlier'


def _rama_stats(angles: list[tuple]) -> dict[str, Any]:
    """Calculate percentages of residues in favored, allowed, and outlier regions."""
    total = len(angles)
    counts = {'Favored': 0, 'Allowed': 0, 'Outlier': 0}
    for _, _, phi, psi in angles:
        counts[_rama_classify(phi, psi)] += 1
    pct = {k: (v / total * 100 if total else 0.0) for k, v in counts.items()}
    return {'total': total, 'counts': counts, 'pct': pct}


def _draw_rama_background(ax) -> None:
    """Draw the standard alpha/beta/L region backgrounds on a Ramachandran axes."""
    import matplotlib.pyplot as plt
    fav_c = '#dcedc8' # light green
    all_c = '#fff9c4' # light yellow
    
    alpha_fav  = plt.Polygon([(-165,-70),(-30,-70),(-30,50),(-165,50)], closed=True, fc=fav_c, ec='none', zorder=0)
    beta_fav1  = plt.Polygon([(-180,110),(-50,110),(-50,180),(-180,180)], closed=True, fc=fav_c, ec='none', zorder=0)
    beta_fav2  = plt.Polygon([(-180,-180),(-50,-180),(-50,-155),(-180,-155)], closed=True, fc=fav_c, ec='none', zorder=0)
    lhand_fav  = plt.Polygon([(30,-25),(90,-25),(90,80),(30,80)], closed=True, fc=fav_c, ec='none', zorder=0)
    
    allowed1   = plt.Polygon([(-180,-100),(0,-100),(0,80),(-180,80)], closed=True, fc=all_c, ec='none', zorder=0)
    allowed2   = plt.Polygon([(-180,80),(-30,80),(-30,180),(-180,180)], closed=True, fc=all_c, ec='none', zorder=0)
    allowed3   = plt.Polygon([(-180,-180),(-30,-180),(-30,-130),(-180,-130)], closed=True, fc=all_c, ec='none', zorder=0)
    allowed4   = plt.Polygon([(0,-50),(130,-50),(130,100),(0,100)], closed=True, fc=all_c, ec='none', zorder=0)
    
    for patch in [allowed1, allowed2, allowed3, allowed4]:
        ax.add_patch(patch)
    for patch in [alpha_fav, beta_fav1, beta_fav2, lhand_fav]:
        ax.add_patch(patch)
        
    ax.axhline(0, color='#bdbdbd', lw=1, zorder=1)
    ax.axvline(0, color='#bdbdbd', lw=1, zorder=1)


def save_ramachandran_comparison(angles_ref: list[tuple], angles_con: list[tuple],
                                 label_ref: str, label_con: str, out_path: Path | str,
                                 critical_res: dict[int, tuple[str, str]] = None,
                                 dpi: int = 300) -> None:
    """Save a side-by-side comparison Ramachandran PNG."""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines
    
    plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans']})
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, angles, label, colour in [
        (axes[0], angles_ref, label_ref, '#6a1b9a'),
        (axes[1], angles_con, label_con, '#1565c0'),
    ]:
        _draw_rama_background(ax)
        if angles:
            for marker_type, m_shape in [('General', 'o'), ('Glycine', '^'), ('Proline', 's')]:
                m_phis, m_psis, m_cols = [], [], []
                for resname, resnum, phi, psi in angles:
                    if (marker_type == 'Glycine' and resname == 'GLY') or \
                       (marker_type == 'Proline' and resname == 'PRO') or \
                       (marker_type == 'General' and resname not in ['GLY', 'PRO']):
                        m_phis.append(phi)
                        m_psis.append(psi)
                        classification = _rama_classify(phi, psi)
                        if classification == 'Favored':
                            m_cols.append('#2e7d32')
                        elif classification == 'Allowed':
                            m_cols.append('#f57f17')
                        else:
                            m_cols.append('#c62828')
                if m_phis:
                    ax.scatter(m_phis, m_psis, c=m_cols, marker=m_shape, s=20, 
                               alpha=0.85, linewidths=0.4, edgecolors='white', zorder=3)
            
            if critical_res:
                triad_styles = {
                    'Acid': {'marker': 'D', 'fc': '#d81b60', 'ec': '#880e4f', 's': 80},
                    'Base': {'marker': 'p', 'fc': '#1e88e5', 'ec': '#0d47a1', 's': 100},
                    'Nuc':  {'marker': '*', 'fc': '#00bcd4', 'ec': '#006064', 's': 180}
                }
                for resname, resnum, phi, psi in angles:
                    if resnum in critical_res and critical_res[resnum][0] == resname:
                        role = critical_res[resnum][1]
                        style = triad_styles.get(role, {'marker': 'X', 'fc': '#9c27b0', 'ec': '#4a148c', 's': 100})
                        ax.scatter([phi], [psi], c=style['fc'], marker=style['marker'], s=style['s'], 
                                   alpha=1.0, linewidths=0.8, edgecolors=style['ec'], zorder=5)
        
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
                va='top', ha='right', multialignment='left', family='monospace',
                bbox=dict(fc='#ffffff', alpha=0.10, ec='#bdbdbd', boxstyle='round,pad=0.4'))
                
        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_aspect('equal')
        ax.set_xlabel('φ (phi) °', fontsize=11, fontweight='500')
        ax.set_ylabel('ψ (psi) °', fontsize=11, fontweight='500')
        ax.set_title(label, fontsize=12, fontweight='bold', pad=10)
        ax.set_xticks(range(-180, 181, 60))
        ax.set_yticks(range(-180, 181, 60))
        ax.tick_params(labelsize=9)
        ax.grid(True, linestyle=':', alpha=0.6, color='#9e9e9e', zorder=1)

    favoured_p = mpatches.Patch(color='#2e7d32', label='Favoured')
    allowed_p  = mpatches.Patch(color='#f57f17', label='Allowed')
    outlier_p  = mpatches.Patch(color='#c62828', label='Outlier')
    
    gen_m = mlines.Line2D([], [], color='none', marker='o', markerfacecolor='gray', markeredgecolor='white', markersize=7, label='General')
    gly_m = mlines.Line2D([], [], marker='^', color='none', markerfacecolor='gray', markeredgecolor='white', markersize=7, label='Glycine')
    pro_m = mlines.Line2D([], [], marker='s', color='none', markerfacecolor='gray', markeredgecolor='white', markersize=7, label='Proline')
    crit_handles = []
    if critical_res:
        triad_styles = {
            'Acid': {'marker': 'D', 'fc': '#d81b60', 'ec': '#880e4f'},
            'Base': {'marker': 'p', 'fc': '#1e88e5', 'ec': '#0d47a1'},
            'Nuc':  {'marker': '*', 'fc': '#00bcd4', 'ec': '#006064'}
        }
        for resnum, (resname, role) in sorted(critical_res.items(), key=lambda x: x[1][1]):
            style = triad_styles.get(role, {'marker': 'X', 'fc': '#9c27b0', 'ec': '#4a148c'})
            handle = mlines.Line2D([], [], color='none', marker=style['marker'], 
                                   markerfacecolor=style['fc'], markeredgecolor=style['ec'], 
                                   markersize=9, label=f"{role}: {resname}{resnum}")
            crit_handles.append(handle)

    all_handles = [favoured_p, allowed_p, outlier_p, gen_m, gly_m, pro_m] + crit_handles
    
    fig.legend(handles=all_handles,
               loc='upper center', ncol=len(all_handles), fontsize=9, frameon=True,
               bbox_to_anchor=(0.5, -0.005), columnspacing=0.8, handletextpad=0.4)

    fig.suptitle('Ramachandran Comparison', fontsize=15, fontweight='bold', y=1.05)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


def save_ramachandran_plot(angles: list[tuple], title: str, out_path: Path | str, dpi: int = 300) -> None:
    """Save a single-structure Ramachandran plot PNG."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 6))
    _draw_rama_background(ax)
    if angles:
        phis = [a[2] for a in angles]
        psis = [a[3] for a in angles]
        colours = ['#1b5e20' if _rama_classify(p, s) == 'Favored'
                   else ('#f9a825' if _rama_classify(p, s) == 'Allowed'
                         else '#c62828')
                   for p, s in zip(phis, psis)]
        ax.scatter(phis, psis, c=colours, s=14, alpha=0.75, linewidths=0, zorder=3)
    stats = _rama_stats(angles)
    legend_txt = (f"Favored  {stats['pct']['Favored']:.1f}%  ({stats['counts']['Favored']})\n"
                  f"Allowed  {stats['pct']['Allowed']:.1f}%  ({stats['counts']['Allowed']})\n"
                  f"Outlier   {stats['pct']['Outlier']:.1f}%  ({stats['counts']['Outlier']})\n"
                  f"Total: {stats['total']} residues")
    ax.text(0.98, 0.98, legend_txt, transform=ax.transAxes, fontsize=8,
            va='top', ha='right', family='monospace',
            bbox=dict(fc='white', alpha=0.7, ec='#cccccc', boxstyle='round,pad=0.3'))
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel('φ (phi) °', fontsize=11)
    ax.set_ylabel('ψ (psi) °', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.tick_params(labelsize=9)
    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)


# ===============================================================================
# SECTION 5: GEOMETRIC & BIO-MATHEMATICAL UTILITIES
# ===============================================================================

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
    if hasattr(pos1, 'x') and hasattr(pos1, 'y') and hasattr(pos1, 'z'):
        pos1 = np.array([pos1.x, pos1.y, pos1.z], dtype=float)
    else:
        pos1 = np.asarray(pos1, dtype=float)

    if hasattr(pos2, 'x') and hasattr(pos2, 'y') and hasattr(pos2, 'z'):
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
            b3   = np.asarray(box, dtype=float)
            diag = np.diag(b3) if b3.ndim == 2 else (b3 if b3.size == 3 else np.array([b3.flat[0], b3.flat[4], b3.flat[8]]))
            diag = np.where(diag == 0, 1e-6, diag)
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
            b3   = np.asarray(box, dtype=float)
            diag = np.diag(b3) if b3.ndim == 2 else (b3 if b3.size == 3 else np.array([b3.flat[0], b3.flat[4], b3.flat[8]]))
            diag = np.where(diag == 0, 1e-6, diag)
            vecs = vecs - diag * np.round(vecs / diag)
    return np.linalg.norm(vecs, axis=2)


def mic_vecs_1d(pos_batch: np.ndarray, pos_ref: np.ndarray, box) -> np.ndarray:
    """
    PBC-corrected displacement vectors (N, 3): pos_batch[i] − pos_ref for all i.
    """
    vecs = pos_batch - pos_ref[None, :]
    if box is not None:
        try:
            b3    = _ensure_box_3x3(box)
            inv_b = np.linalg.inv(b3)
            frac  = vecs @ inv_b
            frac -= np.round(frac)
            vecs  = frac @ b3
        except (np.linalg.LinAlgError, AttributeError, ValueError):
            b3   = np.asarray(box, dtype=float)
            diag = np.diag(b3) if b3.ndim == 2 else (b3 if b3.size == 3 else np.array([b3.flat[0], b3.flat[4], b3.flat[8]]))
            diag = np.where(diag == 0, 1e-6, diag)
            vecs = vecs - diag * np.round(vecs / diag)
    return vecs


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
    b1 /= np.linalg.norm(b1)
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

    Returns the angle between the p4 vector and the normal of the p2–p1–p3 plane
    (degrees from planarity).  Values near 0° indicate a near-planar TS geometry.
    """
    v1 = get_mic_vector(p2, p1, box)
    v2 = get_mic_vector(p3, p1, box)
    v3 = get_mic_vector(p4, p1, box)
    n = np.cross(v1, v2)
    n_len = np.linalg.norm(n)
    if n_len < 1e-6:
        return 0.0
    n /= n_len
    cos_theta = np.dot(v3, n) / np.linalg.norm(v3)
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
    best_dist = float('inf')
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

def calculate_flippin_lodge(nuc_pos, c_pos, o_pos, r1_pos, r2_pos) -> float:
    """
    Flippin–Lodge torsional pre-alignment angle in degrees.

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
            if hasattr(p, 'x'):
                return np.array([float(p.x), float(p.y), float(p.z)])
            return np.array(p, dtype=float)

        nuc = _arr(nuc_pos)
        c   = _arr(c_pos)
        r1  = _arr(r1_pos)
        r2  = _arr(r2_pos)

        vec_nuc = nuc - c
        vec_r1  = r1  - c
        vec_r2  = r2  - c

        normal = np.cross(vec_r1, vec_r2)
        normal /= np.linalg.norm(normal)

        vec_nuc_proj  = vec_nuc - np.dot(vec_nuc, normal) * normal
        vec_nuc_proj /= np.linalg.norm(vec_nuc_proj)

        bisector  = vec_r1 / np.linalg.norm(vec_r1) + vec_r2 / np.linalg.norm(vec_r2)
        bisector /= np.linalg.norm(bisector)

        return float(
            np.degrees(np.arccos(np.clip(np.dot(vec_nuc_proj, bisector), -1.0, 1.0)))
        )
    except Exception:
        return 999.0


# ===============================================================================
# SECTION 6: DATA PIPELINE UTILITIES
# ===============================================================================

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
    df[col] = df[col].fillna(cfg.TIER_DECOY).replace('None', cfg.TIER_DECOY)

    # 2. Enforce Categorical Type with Config Order
    df[col] = pd.Categorical(df[col], categories=cfg.TIER_ORDER, ordered=True)
    return df


def get_alignment_grade(identity_pct, cfg) -> str:
    """
    Resolve a letter grade (A-I) for a given sequence identity percentage
    using bins defined in Project Config.
    """
    for lo, hi, letter in cfg.ALIGN_GRADE_DEFS:
        if lo <= identity_pct <= hi:
            return letter
    return 'I'  # Fallback for ultra-low identity


def get_grade_full_label(grade: str, cfg) -> str:
    """Return the expanded label (e.g. 'A (≥90%)') for a grade letter."""
    for lo, hi, letter in cfg.ALIGN_GRADE_DEFS:
        if letter == grade:
            if lo >= 90: return f"{letter} (≥{int(lo)}%)"
            if lo == 0:  return f"{letter} (<{int(hi)}%)"
            return f"{letter} ({int(lo)}–{int(hi)}%)"
    return f"{grade} (Unknown)"

