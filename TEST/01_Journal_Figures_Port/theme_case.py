"""Shared theme, one colour system and helpers for the 20 figures — Python port of R/theme_case.R.

Design rules enforced everywhere so the panels read as one balanced system:
  • NO figure titles (panels carry only a small lower-case letter tag).
  • ONE font size scheme: axis labels and tick labels share a size, legends one notch smaller.
  • ONE colour system: TIER_COLORS (7 tiers) / OKABE_ITO (other categories) / SIG (up-down-ns) for
    categorical; a single sequential map (SEQ) and one diverging map (DIV) for every heatmap-like panel.
"""
from __future__ import annotations
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

# ---- palettes (Okabe-Ito, verbatim from R/theme_case.R) ---------------------------------------
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442",
             "#0072B2", "#D55E00", "#CC79A7", "#999999"]
SIG = {"up": "#C0392B", "down": "#2166AC", "ns": "#CCCCCC"}

TIER_ORDER = ["Tier_1A", "Tier_1B", "Tier_2A", "Tier_2B", "Tier_3", "Tier_4", "Tier_5_Decoy"]
TIER_COLORS = {
    "Tier_1A": "#08519C", "Tier_1B": "#3182BD", "Tier_2A": "#009E73",
    "Tier_2B": "#E6C700", "Tier_3": "#E69F00", "Tier_4": "#D55E00", "Tier_5_Decoy": "#9AA0A6",
}
# One sequential and one diverging map used by EVERY heatmap/continuous panel, for colour balance.
SEQ = "cividis"
DIV = "RdBu_r"

# ---- one font scheme --------------------------------------------------------------------------
FS_LABEL = 9.0      # x and y axis labels — identical
FS_TICK = 8.0       # tick labels
FS_TICK_DENSE = 6.0 # dense categorical tick labels (many ligands/metrics)
FS_LEGEND = 7.5
FS_ANNOT = 6.5
FS_TAG = 12.0       # panel letter

FIG_DIR = Path(__file__).resolve().parent / "figures"
FIG_DIR.mkdir(exist_ok=True)


def apply_theme() -> None:
    """rcParams mirroring ggplot2 theme_case, minus titles, with one uniform font scheme."""
    mpl.rcParams.update({
        "figure.facecolor": "white", "axes.facecolor": "white", "savefig.facecolor": "white",
        "savefig.dpi": 300,
        "font.family": "DejaVu Sans", "font.size": FS_LABEL,
        "axes.edgecolor": "black", "axes.linewidth": 0.5,
        "axes.grid": True, "axes.axisbelow": True,
        "grid.color": "#EBEBEB", "grid.linewidth": 0.3,
        "axes.labelsize": FS_LABEL, "axes.labelcolor": "black",
        "xtick.labelsize": FS_TICK, "ytick.labelsize": FS_TICK,
        "xtick.color": "black", "ytick.color": "black",
        "xtick.labelcolor": "black", "ytick.labelcolor": "black",
        "xtick.major.width": 0.4, "ytick.major.width": 0.4,
        "axes.titlesize": 0.0,           # titles suppressed by convention
        "legend.frameon": False, "legend.fontsize": FS_LEGEND,
    })


def style_panel(ax) -> None:
    """Full black border box, minor grid off — the per-axes part of theme_case."""
    ax.grid(which="minor", visible=False)
    for s in ax.spines.values():
        s.set_visible(True); s.set_color("black"); s.set_linewidth(0.5)


def tag(container, letter: str) -> None:
    """Small bold lower-case panel letter, top-left — a panel identifier, not a title."""
    container.text(0.008, 0.985, letter, ha="left", va="top", fontsize=FS_TAG, fontweight="bold")


def cbar(fig, mappable, ax, label: str, **kw):
    """One consistent colourbar style: slim, labelled at FS_LABEL, ticks at FS_TICK."""
    cb = fig.colorbar(mappable, ax=ax, fraction=kw.pop("fraction", 0.046), pad=kw.pop("pad", 0.03), **kw)
    cb.set_label(label, fontsize=FS_LABEL)
    cb.ax.tick_params(labelsize=FS_TICK)
    cb.outline.set_linewidth(0.4)
    return cb


def save_case(fig, name: str) -> Path:
    """Export one panel at 300 dpi PNG on white (individual review copies)."""
    out = FIG_DIR / f"{name}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return out


def squarify(sizes, x, y, dx, dy):
    """Bruls squarified-treemap layout (vendored, no external package) → list of {x,y,dx,dy} rects."""
    sizes = list(sizes); total = float(sum(sizes))
    if total <= 0:
        return []
    sizes = [s * dx * dy / total for s in sizes]

    def worst(row, length):
        s = sum(row); mx = max(row); mn = min(row)
        return max((length ** 2) * mx / (s ** 2), (s ** 2) / ((length ** 2) * mn))

    def layout_row(row, x, y, dx, dy):
        cov = sum(row); rects = []
        if dx >= dy:
            w = cov / dy; yy = y
            for r in row:
                hh = r / w; rects.append({"x": x, "y": yy, "dx": w, "dy": hh}); yy += hh
            return rects, x + w, y, dx - w, dy
        h = cov / dx; xx = x
        for r in row:
            ww = r / h; rects.append({"x": xx, "y": y, "dx": ww, "dy": h}); xx += ww
        return rects, x, y + h, dx, dy - h

    out = []; row = []; rx, ry, rdx, rdy = x, y, dx, dy
    for s in sizes:
        length = min(rdx, rdy)
        if not row or worst(row + [s], length) <= worst(row, length):
            row.append(s)
        else:
            rects, rx, ry, rdx, rdy = layout_row(row, rx, ry, rdx, rdy); out += rects; row = [s]
    if row:
        rects, *_ = layout_row(row, rx, ry, rdx, rdy); out += rects
    return out
