"""Render each case as a standalone 300-dpi PNG (review copies) — no titles, panel-consistent style."""
from __future__ import annotations
import time, traceback
import matplotlib.pyplot as plt
import cases
from theme_case import save_case

NAMES = [
    "01_individualized_error_dotplot", "02_grouped_error_dotplot", "03_multigroup_volcano",
    "04_manhattan_twas", "05_paired_boxplot", "06_raincloud", "07_swimmer",
    "08_importance_streams", "09_variance_bars_letters", "10_circos_chord", "11_module_network",
    "12_sankey_enrichment_bubble", "13_discrete_heatmap", "14_mantel_composite",
    "15_correlation_heatmap", "16_polar_heatmap", "17_multilevel_sankey", "18_treemap",
    "19_mosaic_sunburst", "20_split_violin",
]

def main():
    print("Rendering 20 individual panels …")
    ok = fail = 0
    for i, fn in enumerate(cases.ALL):
        t0 = time.time()
        try:
            fig = plt.figure(figsize=cases.SIZE.get(i, (7, 5)))
            fn(fig)
            save_case(fig, NAMES[i])
            print(f"  ok  {NAMES[i]}  ({time.time()-t0:.1f}s)"); ok += 1
        except Exception as exc:
            print(f"  FAIL {fn.__name__}: {exc}"); traceback.print_exc(); fail += 1
    print(f"Done: {ok} ok, {fail} failed → figures/")
    return fail

if __name__ == "__main__":
    raise SystemExit(main())
