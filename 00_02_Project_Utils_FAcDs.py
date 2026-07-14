"""
===============================================================================
FAcDs Pipeline  |  MODULE 00_02  |  Shared Utilities
===============================================================================
Canonical source for console styling, logging infrastructure, matplotlib
spine helpers, MIC vector arithmetic, and geometric angle/dihedral functions.
All downstream scripts import from here — never duplicate these definitions.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 10 July 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Module        : 00_02_Project_Utils_FAcDs.py
  Role          : Shared utility library; no executable entry point.
  Imported by   : 01_Merge_FAcDs.py, 02_Production_FAcDs.py, 03_Validation_Figures_FAcDs.py,
                  04_Dendrogram_FAcDs.py, 05_TopN_and_PDB_Preparation_FAcDs.py,
                  05_TopN_and_PDB_Preparation_FAcDs.py, 06_SID_Prime-MMGBSA_FAcDs.py,
                  07_MD_QMMM_Defluorination_FAcDs.py
                  (also referenced by 00_01 for a ConsoleColours drift check)
  Reads         : (none — pure utility module)
  Writes        : (none — pure utility module)
───────────────────────────────────────────────────────────────────────────────

-------------------------------------------------------------------------------
The Critic's Corner: Known Limitations & Failure Points
-------------------------------------------------------------------------------
  1. Pure library: no executable entry point and no input validation of its own —
     callers must pass well-formed arrays/paths.
  2. Geometry helpers assume Cartesian coordinates in Ångström; the MIC routines
     expect a Schrödinger/Desmond frame.box (3×3 or flat-9 vectors), NOT a
     6-parameter crystallographic cell.
  3. Heavy optional dependencies (gemmi, RDKit, pandas) are imported lazily on
     demand; the relevant helper raises cleanly if the dependency is absent.
-------------------------------------------------------------------------------

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


class _ConsoleRuleFilter:
    """A stdout wrapper that collapses consecutive separator rules.

    The logs grow triple rules — ═══ / ─── / ═══ stacked with nothing between them — because a caller
    prints a rule and then invokes a helper that prints its own. Chasing every call site is endless
    and the next new print re-introduces it, so the rule is enforced where the text is actually
    emitted: a separator that immediately follows another separator, with only blank lines between,
    is dropped. The heavier rule wins, so a section boundary is never demoted to a subsection one.

    Nothing else is touched: any line that is not a rule passes through byte for byte.
    """

    _RULES = {SEPARATOR_HEAVY: 3, SEPARATOR_LIGHT: 2, SEPARATOR_DASH: 1}

    def __init__(self, stream):
        self._stream = stream
        self._pending_rule = None      # rule waiting to be emitted (weight, text)
        self._last_was_rule = False

    def __getattr__(self, name):
        return getattr(self._stream, name)

    def _emit(self, text):
        self._stream.write(text)

    def write(self, data):
        for _line in data.splitlines(keepends=True):
            _bare = _line.strip()
            _w = self._RULES.get(_bare)
            if _w is not None:
                # print() emits the text and its newline as SEPARATE writes, so a buffered rule must
                # carry its own newline or it is glued onto whatever line is emitted next.
                _norm = _bare + "\n"
                if self._last_was_rule:
                    # Already inside a rule run: keep the heaviest, print nothing yet.
                    if self._pending_rule is None or _w > self._pending_rule[0]:
                        self._pending_rule = (_w, _norm)
                    continue
                self._last_was_rule = True
                self._pending_rule = (_w, _norm)
                continue
            if not _bare:
                # Blank lines inside a rule run do not end it — they are what makes the stacks look
                # like separate rules when they are not.
                if self._last_was_rule:
                    continue
                self._emit(_line)
                continue
            if self._pending_rule is not None:
                self._emit(self._pending_rule[1])
                self._pending_rule = None
            self._last_was_rule = False
            self._emit(_line)

    def flush(self):
        if self._pending_rule is not None:
            self._emit(self._pending_rule[1])
            self._pending_rule = None
            self._last_was_rule = False
        self._stream.flush()


def install_console_rule_filter() -> None:
    """Collapse stacked separator rules for the rest of this process's output."""
    import sys as _sys
    if not isinstance(_sys.stdout, _ConsoleRuleFilter):
        _sys.stdout = _ConsoleRuleFilter(_sys.stdout)


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

    def log_file_only(self, text: str):
        """Record a line in the log file WITHOUT printing it to the console.

        For output the console is already showing in another form. A long job draws a progress bar
        rewritten in place (\\r); a milestone line sent through log() would print on top of it and
        break the bar into a ladder of half-finished lines. The file still gets the milestone, which
        is where it is wanted — a log full of carriage returns is unreadable, and a console full of
        milestone lines is a bar that does not work.
        """
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

def apply_figure_style(cfg) -> None:
    """The pipeline's one typography and canvas definition, applied to matplotlib's rcParams.

    Every step that draws a figure calls this instead of setting its own font family, point sizes
    and grid colour: a figure from Step 06 must be indistinguishable in style from one out of Step
    03, and a restyle must be one edit in CFG rather than a hunt through three plotting scripts.

    Weight is deliberately absent from the axis labels. Bolding every label emphasises nothing; the
    point sizes (VIS_FONT_AXIS_LABEL > VIS_FONT_TICK > VIS_FONT_LEGEND > VIS_FONT_ANNOT) already
    carry the hierarchy. Weight is spent only where a label must survive being read against a filled
    bar, and that is set at the call site, not here.
    """
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": list(cfg.VIS_FONT_FAMILY),
        "axes.labelsize": cfg.VIS_FONT_AXIS_LABEL,
        "axes.labelweight": "normal",
        "axes.titlesize": cfg.VIS_FONT_AXIS_LABEL + 2.0,
        "xtick.labelsize": cfg.VIS_FONT_TICK,
        "ytick.labelsize": cfg.VIS_FONT_TICK,
        "legend.fontsize": cfg.VIS_FONT_LEGEND,
        "legend.framealpha": cfg.VIS_LEGEND_FRAME_ALPHA,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": cfg.VIS_GRID_COLOUR,
        "grid.linewidth": cfg.VIS_GRID_LINEWIDTH,
        "savefig.dpi": cfg.VIS_FIGURE_DPI,
        "savefig.bbox": "tight",
    })


def write_json_atomic(path, payload: dict) -> None:
    """Write a JSON file so a reader never sees a half-written one.

    The write goes to a temporary file beside the target and is then renamed over it — rename is
    atomic within a filesystem, so a run killed mid-write leaves either the old file or the new
    one, never a truncated hybrid another step would parse as truth.
    """
    import json
    from pathlib import Path as _Path
    path = _Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    _tmp = path.with_suffix(path.suffix + ".tmp")
    _tmp.write_text(json.dumps(payload, indent=2, sort_keys=True))
    _tmp.replace(path)


def auto_label_colour(cfg, bg, threshold: float = 0.5) -> str:
    """The contrast colour for a label written ON a filled mark, from the fill's luminance.

    Returns near-black on a light fill and white on a dark one, judged by WCAG relative
    luminance. One coefficient set and one cutoff for the whole pipeline: a value written inside a
    bar must stay legible whatever colour that bar happens to take, and hardcoding white works
    until the first pale fill.
    """
    import matplotlib.colors as _mc
    try:
        r, g, b = _mc.to_rgb(bg)
    except Exception:
        return cfg.VIS_BAR_LABEL_COLOURS_ON_LIGHT[0]
    lum = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return (cfg.VIS_BAR_LABEL_COLOURS_ON_LIGHT[0] if lum > threshold
            else cfg.VIS_BAR_LABEL_COLOURS_ON_DARK[0])




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
    """Calculate percentages of residues in favoured, allowed, and outlier regions."""
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
    and the bisector of R1–C–R2.  Ideal value = 0° (the nucleophile approaches ALONG the
    bisector — equidistant from both substituents — which minimises steric clash with them.
    This is a BISECTING arrangement, not an eclipsed one: eclipsing a substituent would
    maximise the clash).

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
# SECTION 8: SCHRÖDINGER JOB-SERVER SCRATCH LOCATION
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
    """Force ALL Schrödinger job scratch onto the working disk so 06/07 run
    entirely on the USB with no manual jsc setup.

    A Prime MM-GBSA / QSite run stages hundreds of GB of per-subjob scratch (the
    multi-GB complexes copies × a 100k-frame trajectory) under the job server's
    LOCAL-SERVER-DIR. When that dir is ``/tmp`` (the default) the OS disk fills and
    the job dies with a silent ``copy_file_range: no space left on device`` (rc=1)
    after hours.

    Step 0 (authoritative): relocate the server's local-server-dir to the working
    disk with ``jsc local-server-dir --set`` — the only lever that actually moves
    it. That command needs the server stopped, so this is done via stop→set→start,
    but ONLY when the server is IDLE (no ``Running`` job in ``jsc list``) — a busy
    server is never stopped, so a live Desmond MD or another user's job is safe; in
    that case the function warns and leaves the server alone.

    Steps 1-3 (fallback): point the ``localhost`` hosts ``tmpdir`` at the working
    disk and ``jsc admin reload-hosts``. NOTE this does NOT move a RUNNING daemon's
    scratch (verified) — it only helps a freshly-started server — so it is a
    best-effort fallback for when step 0 could not run (server busy / no jsc).

    Returns True when scratch is on the working-disk filesystem (server-dir
    relocated, or hosts tmpdir already/now on it); False if neither could be set.

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
    anchor = Path(work_disk_anchor)
    target = (anchor / jobserver_subdir).resolve()
    anchor_dev = _fs_device(anchor)
    global_hosts = Path(schrodinger) / "schrodinger.hosts"       # shared /opt install
    user_hosts   = Path.home() / ".schrodinger" / "schrodinger.hosts"  # per-user override

    def _reload():
        """Apply the hosts change to the running server(s) — no restart, no killed
        job. reload-hosts inherits SCHRODINGER_HOSTS from os.environ, so it reads
        whichever hosts file we selected."""
        for _a in (_jobserver_addresses(jsc) if os.path.isfile(jsc) else []):
            try:
                subprocess.run([jsc, "admin", "reload-hosts", _a],
                               capture_output=True, text=True, timeout=60)
            except Exception:
                pass

    def _tmpdir_on_disk(hp: Path):
        """(current tmpdir string, is-it-on-the-working-disk?) for a hosts file."""
        try:
            t = hp.read_text()
            _m = _re.search(r"(?m)^\s*tmpdir:\s*(.*)$", t)
            cur = _m.group(1).strip() if _m else ""
            return cur, (bool(cur) and anchor_dev is not None
                         and _fs_device(Path(cur)) == anchor_dev)
        except Exception:
            return "", False

    def _write_tmpdir(hp: Path, template: str):
        """Set the localhost tmpdir to `target` in hosts file `hp`, seeding from
        `template` (the global localhost entry) or a minimal localhost entry."""
        base = template if _re.search(r"(?m)^\s*name:\s*localhost", template) else \
            "name:        localhost\ntmpdir:      \n"
        if _re.search(r"(?m)^\s*tmpdir:\s*.*$", base):
            new = _re.sub(r"(?m)^(\s*tmpdir:\s*).*$",
                          lambda mm: f"{mm.group(1)}{target}", base, count=1)
        else:
            new = _re.sub(r"(?m)^(\s*name:\s*localhost\s*)$",
                          lambda mm: f"{mm.group(1)}\ntmpdir:      {target}", base, count=1)
        hp.parent.mkdir(parents=True, exist_ok=True)
        hp.write_text(new)

    target.mkdir(parents=True, exist_ok=True)

    # 0. AUTHORITATIVE fix — relocate the job server's LOCAL-SERVER-DIR to the
    #    working disk. This is where Prime MM-GBSA / QSite subjobs actually stage
    #    their hundreds of GB; the hosts `tmpdir` + `reload-hosts` route (steps 1-3
    #    below) does NOT move it on an already-running daemon (verified: the running
    #    server keeps its launch-time dir on /tmp regardless of reload). The only
    #    lever is `jsc local-server-dir --set` — which requires the server stopped.
    #    Do it automatically, but ONLY when the server is IDLE: never stop a server
    #    that has a RUNNING job (it would kill e.g. a live Desmond MD).
    def _jsc(*args, timeout=120):
        try:
            return subprocess.run([jsc, *args], capture_output=True, text=True, timeout=timeout)
        except Exception:
            return None

    def _server_up() -> bool:
        """True only when the local job server can actually ACCEPT a submission.

        A `jobserverd` process in /proc is NOT proof: the daemon can be alive while its
        registration in the server dir is stale or gone, and Prime then dies at submit time
        with 'Local job submission requires a locally running job server'. Ask the client
        the same way Prime does — `jsc list` — and trust only a clean reply.
        """
        r = _jsc("list", timeout=60)
        if r is None:
            return False
        """
        Judge the MESSAGE, not the exit code: `jsc list` exits 1 merely because there are no
        jobs to list ("No active jobs were found."), which is a perfectly healthy server. Only
        the explicit 'server not running' / 'Error locating ... server' text means it is down.
        """
        return not _re.search(r"not running|Error locating|requires a locally running",
                              f"{r.stdout}\n{r.stderr}", _re.I)

    def _server_dir() -> str:
        r = _jsc("local-server-dir")
        return (r.stdout.strip().splitlines() or [""])[0].strip() if r else ""

    def _server_dir_on_disk() -> bool:
        d = _server_dir()
        return bool(d) and anchor_dev is not None and _fs_device(Path(d)) == anchor_dev

    def _running_jobs() -> int:
        r = _jsc("list")
        # Count ANY non-terminal job, not just "Running": a job in Submitted / Queued /
        # Incorporating / Waiting / Launched / Started would also be killed by a server restart.
        _active = r"\b(Running|Submitted|Queued|Incorporating|Waiting|Launched|Started|Active)\b"
        return sum(1 for ln in r.stdout.splitlines() if _re.search(_active, ln)) if r else 0

    def _jobserverd_tmp() -> str:
        # The daemon creates every job's WORKING directory under its own
        # $SCHRODINGER_TMPDIR / $TMPDIR ($TMPDIR/$USER/jobs). local-server-dir only moves the
        # jobdb, NOT the per-job working dirs — so this env, read from the live daemon, is what
        # decides whether Prime/Desmond subjobs land on the working disk or on /tmp.
        try:
            for _pid in os.listdir("/proc"):
                if not _pid.isdigit():
                    continue
                try:
                    if Path(f"/proc/{_pid}/comm").read_text().strip() != "jobserverd":
                        continue
                    _env = Path(f"/proc/{_pid}/environ").read_bytes().split(b"\0")
                except Exception:
                    continue
                _kv = dict(e.split(b"=", 1) for e in _env if b"=" in e)
                for _key in (b"SCHRODINGER_TMPDIR", b"TMPDIR"):
                    if _kv.get(_key):
                        return _kv[_key].decode("utf-8", "replace")
                return ""   # daemon running but neither var set → it defaults to /tmp
        except Exception:
            pass
        return ""

    def _daemon_tmp_on_disk() -> bool:
        d = _jobserverd_tmp()
        return bool(d) and anchor_dev is not None and _fs_device(Path(d)) == anchor_dev

    if os.path.isfile(jsc):
        """
        Liveness FIRST. Scratch being configured correctly is worthless if the server cannot
        accept a submission: every Prime subjob then fails rc=1 at hand-off, after the frames
        have already been read. Probe with the client, not with /proc.
        """
        _up = _server_up()
        if _up and _server_dir_on_disk() and _daemon_tmp_on_disk():
            say("[job-server] ✔ server is up, and local-server-dir + daemon TMPDIR are both on "
                "the working disk.")
            return True
        if not _up:
            say("[job-server] local job server is DOWN — starting it. Prime and QSite cannot "
                "submit any subjob without it.")
        # A stopped server has no running jobs, so the stop→set→start path below is safe.
        _busy = _running_jobs() if _up else 0
        if _busy == 0:
            # Force the working-disk scratch into this process's env so the RESTARTED daemon
            # inherits it and runs every subjob there instead of /tmp. This is the real fix:
            # local-server-dir alone leaves the daemon's TMPDIR at /tmp, which is where a
            # 21-GB-per-subjob Prime MM-GBSA silently overflows the OS disk.
            for _k in ("SCHRODINGER_TMPDIR", "TMPDIR"):
                _v = os.environ.get(_k, "")
                if not (_v and anchor_dev is not None and _fs_device(Path(_v)) == anchor_dev):
                    os.environ[_k] = str(target)
            _jsc("local-server-stop"); time.sleep(2)
            _jsc("local-server-dir", "--set", str(target)); time.sleep(1)
            _jsc("local-server-start"); time.sleep(2)   # inherits the USB TMPDIR set above
            if _server_up() and _server_dir_on_disk() and _daemon_tmp_on_disk():
                say(f"[job-server] ✔ server started, scratch → working disk (local-server-dir + "
                    f"daemon TMPDIR relocated via stop→set→start; server idle so no job lost).")
                return True
            say("[job-server] relocation did not fully verify; trying hosts-tmpdir fallback.")
        else:
            say(f"[job-server] ⚠ {_busy} job(s) RUNNING (incl. any Desmond) — NOT restarting the "
                f"server (would kill them). Daemon TMPDIR = '{_jobserverd_tmp() or '/tmp (default)'}'; "
                f"jobs keep using it until the server is idle and this runs again. Free the server to "
                f"move scratch onto the working disk ({target}).")

    # 1. Already on the working disk (user-local first, then global)? Then done.
    for hp in (user_hosts, global_hosts):
        if hp.is_file():
            cur, ok = _tmpdir_on_disk(hp)
            if ok:
                # Only pin SCHRODINGER_HOSTS for a NON-default (user-local) file;
                # pointing it at the default $SCHRODINGER/schrodinger.hosts makes
                # Schrödinger emit a noisy "custom hosts file … is being ignored"
                # warning (it reads that path by default anyway).
                if hp == user_hosts:
                    os.environ["SCHRODINGER_HOSTS"] = str(hp)
                else:
                    os.environ.pop("SCHRODINGER_HOSTS", None)
                _reload()
                say(f"[job-server] scratch OK — {hp} tmpdir {cur} is on the working disk.")
                return True

    # 2. PREFER a user-local hosts file (~/.schrodinger). This NEVER edits the shared
    #    /opt hosts, so on a multi-tenant cluster it cannot reroute other users' scratch;
    #    SCHRODINGER_HOSTS points job control at it. tmpdir → the USB working disk.
    _tmpl = ""
    try:
        _tmpl = global_hosts.read_text() if global_hosts.is_file() else ""
    except Exception:
        pass
    try:
        _write_tmpdir(user_hosts, _tmpl)
        os.environ["SCHRODINGER_HOSTS"] = str(user_hosts)
        _reload()
        _cur, _ok = _tmpdir_on_disk(user_hosts)
        if _ok:
            say(f"[job-server] ✔ scratch → {target} via user-local hosts {user_hosts} "
                f"(multi-tenant safe; global /opt hosts untouched).")
            return True
    except Exception as e:
        say(f"[job-server] user-local hosts setup failed ({e}); trying global fallback.")

    # 3. Fallback: edit the GLOBAL hosts, but ONLY when the user owns it (a single-user
    #    install) — never on a shared cluster where it is root-owned. This guarantees a
    #    600 GB Prime/QSite job cannot fill the OS disk on the user's own machine.
    try:
        if global_hosts.is_file() and os.access(global_hosts, os.W_OK):
            _txt = global_hosts.read_text()
            try:
                (global_hosts.parent / "schrodinger.hosts.bak").write_text(_txt)
            except Exception:
                pass
            _write_tmpdir(global_hosts, _txt)
            os.environ.pop("SCHRODINGER_HOSTS", None)  # use the (now-correct) global file
            _reload()
            say(f"[job-server] ✔ scratch → {target} via global hosts (fallback; you own "
                f"this single-user install).")
            return True
    except Exception as e:
        say(f"[job-server] global hosts write failed ({e}).")

    say(f"[job-server] ⚠ could NOT put scratch on the working disk — a large job may "
        f"fill the OS disk. Set 'tmpdir: {target}' in {user_hosts} manually and reload.")
    return False


