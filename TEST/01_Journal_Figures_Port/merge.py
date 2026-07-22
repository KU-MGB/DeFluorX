"""Merge all 20 panels into one master figure — vector SVG, panel letters a–t, no titles.

The panels are drawn as native matplotlib sub-figures (not pasted rasters), so the SVG is fully
vector: infinite zoom, no pixel loss, and it opens in any vector editor. A 4-column × 5-row grid.
"""
from __future__ import annotations
import string
import matplotlib.pyplot as plt
import cases
from theme_case import tag, apply_theme

apply_theme()
NCOL, NROW = 4, 5
LETTERS = list(string.ascii_lowercase[:20])


def build():
    fig = plt.figure(figsize=(15, 18.5), facecolor="white")
    subfigs = fig.subfigures(NROW, NCOL, wspace=0.11, hspace=0.11).ravel()
    for i, fn in enumerate(cases.ALL):
        sf = subfigs[i]
        fn(sf)                                   # draw the case into its own sub-figure
        tag(sf, f"({LETTERS[i]})")               # panel letter, top-left
    return fig


def main():
    print("Building merged master figure (20 panels) …")
    fig = build()
    svg = "figures/00_MASTER_20_panels.svg"
    fig.savefig(svg, format="svg", bbox_inches="tight", facecolor="white")   # vector — no pixel loss
    print(f"  ✔ vector → {svg}")
    plt.close(fig)

if __name__ == "__main__":
    main()
