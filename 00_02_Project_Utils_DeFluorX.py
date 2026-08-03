"""
===============================================================================
DeFluorX Pipeline  |  MODULE 00_02  |  Shared Utilities
===============================================================================
Canonical source for console styling, logging infrastructure, matplotlib
spine helpers, MIC vector arithmetic, and geometric angle/dihedral functions.
All downstream scripts import from here - never duplicate these definitions.

Author : Shaban Ahmad (https://orcid.org/0000-0001-9832-2830)
Date   : 05 August 2026 <────────────────────────────────────────────────────────

── Dependency Map ─────────────────────────────────────────────────────────────
  Module        : 00_02_Project_Utils_DeFluorX.py
  Role          : Shared utility library; no executable entry point.
  Imported by   : 01_Merge_DeFluorX.py, 02_Production_DeFluorX.py, 03_Validation_Figures_DeFluorX.py,
                  04_Dendrogram_DeFluorX.py, 05_TopN_and_PDB_Preparation_DeFluorX.py,
                  06_Physics_Validation_DeFluorX.py,
                  07_MD_QMMM_Defluorination_DeFluorX.py
                  (also referenced by 00_03 for a ConsoleColours drift check)
  Reads         : (none - pure utility module)
  Writes        : (none - pure utility module)
  Upstream      : 00_01_Project_Config_DeFluorX.py (CFG is passed in by callers).
  Downstream    : 01_Merge through 07_MD_QMMM (every pipeline step imports these helpers).
───────────────────────────────────────────────────────────────────────────────

── The Critic's Corner: Known Limitations & Failure Points ──────────────────
  1. Pure library: no executable entry point and no input validation of its own -
     callers must pass well-formed arrays/paths.
  2. Geometry helpers assume Cartesian coordinates in Ångström; the MIC routines
     expect a Schrödinger/Desmond frame.box (3×3 or flat-9 vectors), NOT a
     6-parameter crystallographic cell.
  3. Heavy optional dependencies (gemmi, RDKit, pandas) are imported lazily on
     demand; the relevant helper raises cleanly if the dependency is absent.
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
    definitions and their primary literature live in 00_01_Project_Config_DeFluorX.py.
-------------------------------------------------------------------------------
"""

from __future__ import annotations

import logging
import re as _re
import time
from pathlib import Path
from typing import Any

import numpy as np
# --- consolidated top-level imports (optional/heavy + Schrodinger stay function-local) ---
from datetime import datetime as _dt
import json
import math
import os
import pandas as pd
import sys as _sys
import uuid
import matplotlib
matplotlib.use("Agg")  # headless backend before pyplot
# --- consolidated matplotlib imports ---
import matplotlib.colors as _mc
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

_ANSI_ESCAPE_RE = _re.compile(r"\033\[[0-9;]*[mKABCDEFGHJKSTfhilmnprsu]")


# =============================================================================
# SECTION 1: CONSOLE STYLING
# =============================================================================

class ConsoleColours:
    """ANSI terminal colour codes for pipeline console output."""
    OKGREEN = "\033[92m"   # green   - success / pass
    WARNING = "\033[93m"   # yellow  - caution
    FAIL    = "\033[91m"   # red     - error / fail
    OKBLUE  = "\033[94m"   # blue    - information
    MAGENTA = "\033[95m"   # magenta - script banners
    CYAN    = "\033[96m"   # cyan    - phase / heartbeat headers
    BOLD    = "\033[1m"    # bold    - section headers
    DIM     = "\033[2m"    # dim     - de-emphasised detail
    ENDC    = "\033[0m"    # reset   - end all formatting


# Horizontal separators - choose the weight that matches visual hierarchy.
SEPARATOR_HEAVY = "═" * 80   # major section boundary  (═══)
SEPARATOR_LIGHT = "─" * 80   # step / subsection       (───)
SEPARATOR_DASH  = "-" * 80   # info line / minor break  (---)


class _ConsoleRuleFilter:
    """A stdout wrapper that collapses consecutive separator rules.

    The logs grow triple rules - ═══ / ─── / ═══ stacked with nothing between them - because a caller
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
                # Blank lines inside a rule run do not end it - they are what makes the stacks look
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
    """Bold section header - printed and optionally written to log file."""
    print(f"\n{ConsoleColours.BOLD}{msg}{ConsoleColours.ENDC}", flush=True)
    if logger:
        logger.info(_strip_ansi(msg))


def console_info(msg: str, logger: logging.Logger | None = None) -> None:
    """Two-space-indented info line - printed and optionally written to log file."""
    # Ensure every line of a multi-line message is indented by two spaces.
    for line in str(msg).split("\n"):
        print(f"  {line}", flush=True)
    if logger:
        logger.info(_strip_ansi(msg))


def console_separator(
    logger: logging.Logger | None = None,
    heavy: bool = False,
) -> None:
    """Horizontal rule - printed and optionally written to log file.

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
    script_name : Filename of the running script (e.g. "02_Production_DeFluorX.py").
    subtitle    : One-line description shown below the script name.
    """
    _now = time.strftime("%Y-%m-%d %H:%M:%S")
    print(f"\n{SEPARATOR_HEAVY}", flush=True)
    print(
        f"  {ConsoleColours.MAGENTA}{ConsoleColours.BOLD}▶  {script_name}"
        f"{ConsoleColours.ENDC}  │  DeFluorX Pipeline",
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
    script_name : Filename displayed in the output line (e.g. ``"02_Production_DeFluorX.py"``).
    """
    _el = time.perf_counter() - t0
    _h, _rem = divmod(int(_el), 3600)
    _m, _s   = divmod(_rem, 60)
    _fmt = (f"{_h}h {_m:02d}m {_s:02d}s" if _h else
            f"{_m}m {_s:02d}s"            if _m else
            f"{_s}s")
    print(f"\n{SEPARATOR_HEAVY}", flush=True)
    print(
        f"  {ConsoleColours.OKGREEN}✔  {script_name}  -  Pipeline Phase Complete"
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
        self.path = Path(log_path)
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
        is where it is wanted - a log full of carriage returns is unreadable, and a console full of
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

# Image suffixes the SSOT format router (installed in apply_figure_style) rewrites at savefig time.
_DEFLX_IMG_SUFFIXES = {".png", ".jpg", ".jpeg", ".tif", ".tiff", ".pdf", ".eps", ".ps", ".svg"}


def deflx_fig_name(path) -> str:
    """Return a figure path with its suffix swapped to the ACTIVE SSOT format (CFG.VIS_FIGURE_FORMAT,
    as recorded on matplotlib by apply_figure_style) - so a LOG line names the file savefig actually
    wrote, not the source-literal '.png'. Non-image suffixes (e.g. .csv) are returned unchanged, so it
    is safe to apply even where a path might not be a figure."""
    try:
        p = Path(path)
        if p.suffix.lower() in _DEFLX_IMG_SUFFIXES:
            fmt = str(getattr(matplotlib, "_deflx_fig_fmt", "png")).lower().lstrip(".")
            return str(p.with_suffix("." + fmt))
    except Exception:
        pass
    return str(path)


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
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": list(cfg.VIS_FONT_FAMILY),
        "axes.labelsize": cfg.VIS_FONT_AXIS_LABEL,
        "axes.labelweight": "normal",
        "axes.titlesize": cfg.VIS_FONT_AXIS_LABEL + 2.0,
        "xtick.labelsize": cfg.VIS_FONT_TICK,
        "ytick.labelsize": cfg.VIS_FONT_TICK,
        "legend.fontsize": cfg.VIS_FONT_LEGEND,
        "legend.title_fontsize": cfg.VIS_FONT_LEGEND_TITLE,
        "legend.framealpha": cfg.VIS_LEGEND_FRAME_ALPHA,
        "legend.fancybox": True,
        "legend.labelspacing": cfg.VIS_LEGEND_LABELSPACING,
        "legend.handletextpad": cfg.VIS_LEGEND_HANDLETEXTPAD,
        "legend.columnspacing": cfg.VIS_LEGEND_COLUMNSPACING,
        "legend.borderpad": cfg.VIS_LEGEND_BORDERPAD,
        "legend.handlelength": cfg.VIS_LEGEND_HANDLELENGTH,
        "legend.borderaxespad": cfg.VIS_LEGEND_BORDERAXESPAD,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.color": cfg.VIS_GRID_COLOUR,
        "grid.linewidth": cfg.VIS_GRID_LINEWIDTH,
        "savefig.dpi": cfg.VIS_FIGURE_DPI,
        "savefig.bbox": "tight",
    })

    # ---- Central figure-format SSOT ----------------------------------------------------------------
    # Every matplotlib save (fig.savefig / plt.savefig) is routed to cfg.VIS_FIGURE_FORMAT regardless
    # of the ".png" a call site happens to name: the suffix and the write format are rewritten here, so
    # one CFG edit re-targets all plots (svg / pdf / tiff / png / jpg). A call that passes an explicit
    # format= (a PIL-composited or PyMOL raster tile that must stay png) is left untouched. The patch is
    # installed once at the class level (idempotent guard) and reads the format from matplotlib each
    # call, so a later apply_figure_style with a different CFG re-targets without re-patching. It
    # composes with Step-03's savefig redirect: that wrapper rewrites the folder, then delegates to this
    # one which rewrites the extension.
    matplotlib._deflx_fig_fmt = str(getattr(cfg, "VIS_FIGURE_FORMAT", "png")).lower().lstrip(".")
    if not getattr(matplotlib, "_deflx_savefig_patched", False):
        _SWAP = {"png", "jpg", "jpeg", "tif", "tiff", "pdf", "svg", "eps", "ps"}
        _orig_fig_savefig = plt.Figure.savefig
        _orig_plt_savefig = plt.savefig
        def _deflx_route(fname):
            try:
                _p = Path(fname)
                if _p.suffix.lower().lstrip(".") in _SWAP:
                    return str(_p.with_suffix("." + getattr(matplotlib, "_deflx_fig_fmt", "png")))
            except Exception:
                pass
            return fname
        def _deflx_fig_savefig(self, fname, *a, **k):
            if "format" not in k:
                fname = _deflx_route(fname)
            return _orig_fig_savefig(self, fname, *a, **k)
        def _deflx_plt_savefig(*a, **k):
            if "format" not in k:
                if a:
                    a = (_deflx_route(a[0]),) + a[1:]
                elif "fname" in k:
                    k["fname"] = _deflx_route(k["fname"])
            return _orig_plt_savefig(*a, **k)
        plt.Figure.savefig = _deflx_fig_savefig
        plt.savefig = _deflx_plt_savefig
        matplotlib._deflx_savefig_patched = True


def latest_by_mtime(paths):
    """The newest existing file by modification time, or None.

    Selection is by mtime, NOT by name: a name sort ranks on the leading number, so a file with a higher
    leading digit sorts last even when it is older - which would silently feed a stale file downstream.
    mtime is prefix-agnostic and always returns the file written last. Callers pass the SSOT glob first
    and a wildcard fallback second, e.g.
        latest_by_mtime(prod.glob(CFG.GLOB_RANKED_CSV)) or latest_by_mtime(prod.glob("*Ranked*.csv"))
    """
    _ps = [Path(p) for p in paths if Path(p).exists()]
    return max(_ps, key=lambda p: p.stat().st_mtime) if _ps else None


def write_json_atomic(path, payload: dict) -> None:
    """Write a JSON file so a reader never sees a half-written one.

    The write goes to a scratch file beside the target and is then renamed over it - rename is
    atomic within a filesystem, so a run killed mid-write leaves either the previous file or the new
    one, never a truncated hybrid another step would parse as truth.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Per-process/per-call unique temp name: a fixed "<file>.tmp" would let two workers writing the same
    # path clobber each other's half-written temp before the rename. pid + a short uuid make it collision-free.
    _tmp = path.with_name(f"{path.name}.{os.getpid()}.{uuid.uuid4().hex[:8]}.tmp")
    _tmp.write_text(json.dumps(payload, indent=2, sort_keys=True))
    _tmp.replace(path)


def atomic_write_csv(df, path, **to_csv_kwargs) -> None:
    """Write a DataFrame to CSV so a reader never sees a half-written file.

    Same guarantee as write_json_atomic: the frame goes to a scratch file beside the target and is
    renamed over it (atomic within a filesystem), so a run killed mid-write leaves either the previous file
    or the complete new one, never a truncated hybrid a downstream step would parse as truth.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    # Per-process/per-call unique temp name (see write_json_atomic) so concurrent writers to one path
    # cannot clobber each other's temp before the atomic rename.
    _tmp = path.with_name(f"{path.name}.{os.getpid()}.{uuid.uuid4().hex[:8]}.tmp")
    to_csv_kwargs.setdefault("index", False)
    df.to_csv(_tmp, **to_csv_kwargs)
    _tmp.replace(path)


def auto_label_colour(cfg, bg, threshold: float = 0.5) -> str:
    """The contrast colour for a label written ON a filled mark, from the fill's luminance.

    Returns near-black on a light fill and white on a dark one, judged by WCAG relative
    luminance. One coefficient set and one cutoff for the whole pipeline: a value written inside a
    bar must stay legible whatever colour that bar happens to take, and hardcoding white works
    until the first pale fill.
    """
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
    individual scripts - import and call this function instead.
    """
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(direction="out", length=4, width=0.8)
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


_RAMA_FALLBACK = {
    "favoured": "#2e7d32", "allowed": "#f57f17", "outlier": "#c62828",
    "region_favoured": "#dcedc8", "region_allowed": "#fff9c4",
    "axis": "#bdbdbd", "grid": "#9e9e9e", "box_edge": "#bdbdbd",
    "acid_fc": "#d81b60", "acid_ec": "#880e4f",
    "base_fc": "#1e88e5", "base_ec": "#0d47a1",
    "nuc_fc": "#00bcd4", "nuc_ec": "#006064",
    "default_fc": "#9c27b0", "default_ec": "#4a148c",
    "title": "#1f4e79",
}


def _rama_palette(cfg=None) -> dict:
    """Ramachandran colours from CFG.VIS_RAMA (the single source of truth); the built-in fallback
    (identical values) is used only when no cfg is supplied, so behaviour is unchanged either way."""
    _p = dict(_RAMA_FALLBACK)
    if cfg is not None:
        _p.update(getattr(cfg, "VIS_RAMA", {}) or {})
    return _p


def _draw_rama_background(ax, cfg=None) -> None:
    """Draw the standard alpha/beta/L region backgrounds on a Ramachandran axes."""
    _P = _rama_palette(cfg)
    fav_c = _P["region_favoured"]
    all_c = _P["region_allowed"]

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

    ax.axhline(0, color=_P["axis"], lw=1, zorder=1)
    ax.axvline(0, color=_P["axis"], lw=1, zorder=1)


def save_ramachandran_comparison(angles_ref: list[tuple], angles_con: list[tuple],
                                 label_ref: str, label_con: str, out_path: Path | str,
                                 critical_res: dict[int, tuple[str, str]] = None,
                                 dpi: int = 300, cfg=None) -> None:
    """Save a side-by-side comparison Ramachandran figure. Colours come from CFG.VIS_RAMA when cfg is
    supplied (SSOT), else the identical built-in fallback."""

    _P = _rama_palette(cfg)
    _fam = list(getattr(cfg, "VIS_FONT_FAMILY", ("Arial", "Helvetica", "DejaVu Sans")))
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": _fam})

    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    for ax, angles, label in [
        (axes[0], angles_ref, label_ref),
        (axes[1], angles_con, label_con),
    ]:
        _draw_rama_background(ax, cfg)
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
                            m_cols.append(_P["favoured"])
                        elif classification == "Allowed":
                            m_cols.append(_P["allowed"])
                        else:
                            m_cols.append(_P["outlier"])
                if m_phis:
                    ax.scatter(m_phis, m_psis, c=m_cols, marker=m_shape, s=20,
                               alpha=0.85, linewidths=0.4, edgecolors="white", zorder=3)

            if critical_res:
                triad_styles = {
                    "Acid": {"marker": "D", "fc": _P["acid_fc"], "ec": _P["acid_ec"], "s": 80},
                    "Base": {"marker": "p", "fc": _P["base_fc"], "ec": _P["base_ec"], "s": 100},
                    "Nuc":  {"marker": "*", "fc": _P["nuc_fc"], "ec": _P["nuc_ec"], "s": 180}
                }
                for resname, resnum, phi, psi in angles:
                    if resnum in critical_res and critical_res[resnum][0] == resname:
                        role = critical_res[resnum][1]
                        style = triad_styles.get(role, {"marker": "X", "fc": _P["default_fc"], "ec": _P["default_ec"], "s": 100})
                        ax.scatter([phi], [psi], c=style["fc"], marker=style["marker"], s=style["s"],
                                   alpha=1.0, linewidths=0.8, edgecolors=style["ec"], zorder=5)

        st = _rama_stats(angles)
        fav_pct = f"({st['pct']['Favored']:.1f}%)"
        allowed_pct = f"({st['pct']['Allowed']:.1f}%)"
        out_pct = f"({st['pct']['Outlier']:.1f}%)"
        tot_pct = "(100.0%)"

        _rows = [f"Favoured  {st['counts']['Favored']:<4} {fav_pct:<8}",
                 f"Allowed   {st['counts']['Allowed']:<4} {allowed_pct:<8}",
                 f"Outlier   {st['counts']['Outlier']:<4} {out_pct:<8}",
                 f"Total     {st['total']:<4} {tot_pct:<8}"]
        _fs = plt.rcParams["xtick.labelsize"]
        # Structure name as a coloured, centred title over the left-justified stats (monospace columns).
        _leg = ax.legend([mlines.Line2D([], [], color="none") for _ in _rows], _rows,
                         loc="upper right", title=label, handlelength=0, handletextpad=0,
                         labelspacing=0.3, borderpad=0.6, prop={"family": "monospace", "size": _fs},
                         framealpha=0.85, edgecolor=_P["box_edge"], facecolor="white")
        _t = _leg.get_title()
        _t.set_color(_P["title"]); _t.set_fontfamily("monospace"); _t.set_fontsize(_fs); _t.set_fontweight("bold")
        _leg.set_zorder(6)

        ax.set_xlim(-180, 180)
        ax.set_ylim(-180, 180)
        ax.set_aspect("equal")
        ax.set_xlabel("φ (phi) °", fontweight="500")
        ax.set_ylabel("ψ (psi) °", fontweight="500")
        ax.set_xticks(range(-180, 181, 60))
        ax.set_yticks(range(-180, 181, 60))
        ax.tick_params()
        ax.grid(True, linestyle=":", alpha=0.6, color=_P["grid"], zorder=1)

    favoured_p = mpatches.Patch(color=_P["favoured"], label="Favoured")
    allowed_p  = mpatches.Patch(color=_P["allowed"], label="Allowed")
    outlier_p  = mpatches.Patch(color=_P["outlier"], label="Outlier")

    gen_m = mlines.Line2D([], [], color="none", marker="o", markerfacecolor="gray", markeredgecolor="white", markersize=7, label="General")
    gly_m = mlines.Line2D([], [], marker="^", color="none", markerfacecolor="gray", markeredgecolor="white", markersize=7, label="Glycine")
    pro_m = mlines.Line2D([], [], marker="s", color="none", markerfacecolor="gray", markeredgecolor="white", markersize=7, label="Proline")
    crit_handles = []
    if critical_res:
        triad_styles = {
            "Acid": {"marker": "D", "fc": _P["acid_fc"], "ec": _P["acid_ec"]},
            "Base": {"marker": "p", "fc": _P["base_fc"], "ec": _P["base_ec"]},
            "Nuc":  {"marker": "*", "fc": _P["nuc_fc"], "ec": _P["nuc_ec"]}
        }
        for resnum, (resname, role) in sorted(critical_res.items(), key=lambda x: x[1][1]):
            style = triad_styles.get(role, {"marker": "X", "fc": _P["default_fc"], "ec": _P["default_ec"]})
            handle = mlines.Line2D([], [], color="none", marker=style["marker"],
                                   markerfacecolor=style["fc"], markeredgecolor=style["ec"],
                                   markersize=9, label=f"{role}: {resname}{resnum}")
            crit_handles.append(handle)

    all_handles = [favoured_p, allowed_p, outlier_p, gen_m, gly_m, pro_m] + crit_handles

    fig.legend(handles=all_handles,
               loc="upper center", ncol=len(all_handles), fontsize=plt.rcParams["xtick.labelsize"], frameon=True,
               bbox_to_anchor=(0.5, -0.005), columnspacing=0.8, handletextpad=0.4)

    fig.tight_layout()
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def save_ramachandran_plot(angles: list[tuple], title: str, out_path: Path | str,
                           dpi: int = 300, cfg=None) -> None:
    """Save a single-structure Ramachandran figure. Colours from CFG.VIS_RAMA when cfg is supplied."""
    _P = _rama_palette(cfg)
    fig, ax = plt.subplots(figsize=(6, 6))
    _draw_rama_background(ax, cfg)
    if angles:
        phis = [a[2] for a in angles]
        psis = [a[3] for a in angles]
        _rama_cls = [_rama_classify(p, s) for p, s in zip(phis, psis)]
        colours = [_P["favoured"] if c == "Favored"
                   else (_P["allowed"] if c == "Allowed" else _P["outlier"])
                   for c in _rama_cls]
        ax.scatter(phis, psis, c=colours, s=14, alpha=0.75, linewidths=0, zorder=3)
    stats = _rama_stats(angles)
    _lines = [
        f"Favored {stats['pct']['Favored']:5.1f}%  ({stats['counts']['Favored']})",
        f"Allowed {stats['pct']['Allowed']:5.1f}%  ({stats['counts']['Allowed']})",
        f"Outlier {stats['pct']['Outlier']:5.1f}%  ({stats['counts']['Outlier']})",
        f"Total:  {stats['total']} residues",
    ]
    _fs = plt.rcParams["legend.fontsize"]
    # Stats as a boxed legend so the structure-name title carries its own colour, centred over the
    # left-justified rows (monospace keeps the percent/count columns aligned).
    _leg = ax.legend([mlines.Line2D([], [], color="none") for _ in _lines], _lines,
                     loc="upper right", title=title, handlelength=0, handletextpad=0,
                     labelspacing=0.3, borderpad=0.6, prop={"family": "monospace", "size": _fs},
                     framealpha=0.7, edgecolor=_P["box_edge"], facecolor="white")
    _t = _leg.get_title()
    _t.set_color(_P["title"]); _t.set_fontfamily("monospace"); _t.set_fontsize(_fs); _t.set_fontweight("bold")
    _leg.set_zorder(6)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("φ (phi) °")
    ax.set_ylabel("ψ (psi) °")
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.tick_params()
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
    if b1_len < 1e-6:          # coincident central atoms (p2 ≡ p3) - dihedral undefined
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
    Out-of-plane improper dihedral at centre p1 - for detecting transition-state
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
    if v3_len < 1e-6:          # p4 coincides with the centre p1 - undefined out-of-plane angle
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

    AUXILIARY / REFERENCE metric only - NOT a tier gate. In the FAcD pipeline the
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
    bisector - equidistant from both substituents - which minimises steric clash with them.
    This is a BISECTING arrangement, not an eclipsed one: eclipsing a substituent would
    maximise the clash).

    AUXILIARY / REFERENCE metric only - NOT a tier gate. Reports the in-plane
    (lateral) component of the nucleophile approach, complementing the backside
    O–C–F attack angle. Requires TWO well-defined heavy spectator substituents
    (R1, R2) on the attacked carbon; for small substrates whose α-carbon carries
    only hydrogens (often absent in heavy-atom CIF output) the offset is undefined
    - callers must guard and report NaN rather than a spurious value.

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

    The SMILES file is the single source - nothing is hardcoded, so swapping the
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
