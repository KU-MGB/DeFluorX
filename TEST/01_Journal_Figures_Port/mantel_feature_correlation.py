"""Figure 08 (feature-correlation) + Mantel couple-links — TEST proof.

Reproduces the project's real Figure 08 — family-grouped lower-triangle Spearman heatmap with
BH-FDR significance stars, family block borders and family labels — from the Figure-Enriched dataset,
then ADDS the linkET `geom_couple` Mantel overlay in the (otherwise empty) upper-right triangle:
mechanistic driver variables as nodes on the right, each joined to the feature block by a curved link
whose COLOUR encodes the Mantel permutation p and WIDTH the Mantel r. Nothing else about the figure
changes. No title.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.lines import Line2D
from scipy import stats as sps
from scipy.stats import false_discovery_control

from theme_case import apply_theme, FS_LABEL, FS_TICK, FS_ANNOT, FS_LEGEND, save_case

ENRICHED = (Path(__file__).resolve().parent.parent /
            "Boltz-2_Run_20260309T085406Z" / "3_Validation_Figures" /
            "01_Analysis_Data" / "03_Figure_Enriched_Dataset.csv")

RNG = np.random.default_rng(7)
N_SAMPLE, N_PERM = 260, 999

# The 11 Figure-08 features, in family order (matches the production figure).
FEAT = [("tier_numeric", "Tier")]
FAMILIES = [
    ("Outcome",               [("tier_numeric", "Tier")]),
    ("Catalytic geometry",    [("SN2_Attack_Angle", "SN2°"), ("soft_catalytic_score", "SoftCat"),
                               ("SN2_Trajectory_Deviation_A", "TrajDev")]),
    ("Affinity / seq / rank", [("Chemical_Affinity_Score", "ChemAff"), ("identity_pct", "SeqID%"),
                               ("Pareto_Rank", "Pareto")]),
    ("Binding & confidence",  [("Binding_Probability", "BindP"), ("Interaction_Density_Norm", "IntDen"),
                               ("Boltz_Model_Confidence", "Conf"), ("iptm", "ipTM")]),
]
# External mechanistic drivers for the Mantel overlay (NOT in the matrix).
DRIVERS = [("scissile_cf_bde", "C–F BDE"), ("sn2_backside_occlusion", "Backside occ."),
           ("halide_stabilisation_score", "Halide stab."), ("competence_score", "Competence"),
           ("carboxylate_clamp_integrity", "Cbx clamp"), ("mechanistic_score_effective", "Mech. score")]

P_COL = [("< 0.01", "#1B9E77"), ("0.01 – 0.05", "#E69F00"), ("≥ 0.05", "#B8BCC2")]
R_W = [("< 0.2", 1.0), ("0.2 – 0.4", 2.6), ("≥ 0.4", 4.6)]


def _z(a):
    s = a.std(0); s[s == 0] = 1.0
    return (a - a.mean(0)) / s


def _pdist(x):
    d = x[:, None, :] - x[None, :, :]
    return np.sqrt((d ** 2).sum(-1))


def _mantel(dd, dm_tri, iu):
    r = np.corrcoef(dd[iu], dm_tri)[0, 1]
    n = dd.shape[0]; ge = 1
    for _ in range(N_PERM):
        p = RNG.permutation(n)
        if abs(np.corrcoef(dd[np.ix_(p, p)][iu], dm_tri)[0, 1]) >= abs(r):
            ge += 1
    return float(r), ge / (N_PERM + 1)


def _p_colour(p):
    return P_COL[0][1] if p < 0.01 else (P_COL[1][1] if p < 0.05 else P_COL[2][1])


def _r_width(r):
    a = abs(r)
    return R_W[0][1] if a < 0.2 else (R_W[1][1] if a < 0.4 else R_W[2][1])


def _couple(p0, p1, curve=0.40):
    dx, dy = p1[0] - p0[0], p1[1] - p0[1]
    c0 = (p0[0] + dx * curve, p0[1] + dy * 0.05)
    c1 = (p1[0] - dx * curve, p1[1])
    return PathPatch(MPath([p0, c0, c1, p1],
                     [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]), fc="none")


def _stars(p):
    return "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))


def build(host):
    order = [(c, lab) for _, items in FAMILIES for (c, lab) in items]
    cols = [c for c, _ in order]
    labs = [lab for _, lab in order]
    drivers = [(c, lab) for c, lab in DRIVERS]

    df = pd.read_csv(ENRICHED, usecols=lambda c: c in set(cols + [c0 for c0, _ in drivers]),
                     low_memory=False)
    sub_all = df[cols].apply(pd.to_numeric, errors="coerce")
    N = len(cols)

    # Spearman ρ + per-pair p, BH-FDR over unique pairs (as in the production figure).
    corr = sub_all.corr(method="spearman", min_periods=20).values
    p_mat = np.full((N, N), np.nan)
    for i in range(N):
        for j in range(N):
            if i != j and (sub_all.iloc[:, i].notna() & sub_all.iloc[:, j].notna()).sum() >= 20:
                _, pv = sps.spearmanr(sub_all.iloc[:, i], sub_all.iloc[:, j], nan_policy="omit")
                p_mat[i, j] = pv
    flat, idx = [], []
    for i in range(N):
        for j in range(i + 1, N):
            if np.isfinite(p_mat[i, j]):
                flat.append(p_mat[i, j]); idx.append((i, j))
    padj = np.ones((N, N))
    adj = np.clip(false_discovery_control(np.array(flat), method="bh"), 0, 1)
    for k, (i, j) in enumerate(idx):
        padj[i, j] = padj[j, i] = adj[k]

    annot = np.empty((N, N), dtype=object)
    for i in range(N):
        for j in range(N):
            if i == j:
                annot[i, j] = "1.00"
            elif np.isnan(corr[i, j]):
                annot[i, j] = "—"
            else:
                s = _stars(padj[i, j])
                annot[i, j] = f"{corr[i, j]:.2f}\n{s}" if s else f"{corr[i, j]:.2f}"
    mask = np.triu(np.ones((N, N), bool), k=1)

    ax = host.subplots()
    sns.heatmap(pd.DataFrame(corr, index=labs, columns=labs), mask=mask,
                annot=annot, fmt="", cmap="coolwarm", vmin=-1, vmax=1, center=0,
                square=True, linewidths=0.5, annot_kws={"size": FS_ANNOT + 0.5},
                cbar_kws={"label": "Spearman ρ  (−1 = perfect negative, +1 = perfect positive)",
                          "shrink": 0.6, "pad": 0.02}, ax=ax)
    ax.set_xticklabels(labs, fontsize=FS_TICK, rotation=45, ha="right")
    ax.set_yticklabels(labs, fontsize=FS_TICK, rotation=0, va="center")
    for t in ax.texts:                       # star colour on dark cells
        ln = t.get_text().split("\n")
        if len(ln) == 2 and ln[1]:
            try:
                t.set_color("white" if abs(float(ln[0])) > 0.5 else "#20242A")
            except ValueError:
                pass
    ax.text(0.5, 1.02, "★ p<0.05   ★★ p<0.01   ★★★ p<0.001  (BH FDR, unique pairs)   |   "
            "Features grouped by family   |   pairwise-complete Spearman, n ≤ {:,}".format(len(df)),
            transform=ax.transAxes, ha="center", va="bottom", fontsize=FS_ANNOT, style="italic",
            color="#6b6f76", bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#cfd3d8",
                                       alpha=0.9, lw=0.7))

    # family block borders + bottom family labels
    starts, s = [], 0
    for fam, items in FAMILIES:
        e = s + len(items)
        ax.add_patch(Rectangle((s, s), e - s, e - s, fill=False, edgecolor="#1b1b1b",
                               lw=1.6, zorder=6))
        ax.text((s + e) / 2, N + 0.9, fam, ha="center", va="top", fontsize=FS_TICK,
                fontweight="bold", color="#2b2b2b")
        starts.append((s, e)); s = e

    # ---- Mantel overlay in the upper-right triangle -------------------------------
    use = df[cols + [c for c, _ in drivers]].apply(pd.to_numeric, errors="coerce").dropna()
    samp = use.iloc[RNG.choice(len(use), size=min(N_SAMPLE, len(use)), replace=False)]
    dm = _pdist(_z(samp[cols].values)); iu = np.triu_indices_from(dm, k=1); dm_tri = dm[iu]
    mant = [(lab, *_mantel(_pdist(_z(samp[[c]].values)), dm_tri, iu)) for c, lab in drivers]

    x_node = N + 1.9
    y_node = np.linspace(0.7, N * 0.99, len(drivers))
    anchor = np.linspace(0.6, N - 0.5, len(drivers))
    for (lab, r, p), yn, a in zip(mant, y_node, anchor):
        col = _p_colour(p)
        p0 = (a + 0.5, a + 0.5)                       # diagonal cell centre
        patch = _couple(p0, (x_node, yn))
        patch.set_edgecolor(col); patch.set_linewidth(_r_width(r))
        patch.set_alpha(0.9 if p < 0.05 else 0.7)
        patch.set_zorder(5 if p < 0.05 else 4)
        ax.add_patch(patch)
        ax.scatter([x_node], [yn], s=46, color="#33373D", zorder=7, clip_on=False)
        ax.text(x_node + 0.28, yn, lab, ha="left", va="center", fontsize=FS_TICK, zorder=7,
                clip_on=False)

    r_h = [Line2D([0], [0], color="#5b5f66", lw=w, label=l) for l, w in R_W]
    leg1 = ax.legend(handles=r_h, title="Mantel's r  (link width)", loc="upper left",
                     bbox_to_anchor=(0.40, 0.955), fontsize=FS_LEGEND, title_fontsize=FS_LEGEND,
                     frameon=False, handlelength=2.4)
    ax.add_artist(leg1)
    p_h = [Line2D([0], [0], color=c, lw=3.2, label=l) for l, c in P_COL]
    ax.legend(handles=p_h, title="Mantel's p  (link colour)", loc="upper left",
              bbox_to_anchor=(0.63, 0.955), fontsize=FS_LEGEND, title_fontsize=FS_LEGEND,
              frameon=False, handlelength=2.4)

    ax.set_xlim(0, x_node + 4.0)
    ax.set_ylim(N + 1.6, -0.4)
    ax.set_xlabel(""); ax.set_ylabel("")


def main():
    apply_theme()
    fig = plt.figure(figsize=(12.5, 11), facecolor="white")
    build(fig)
    out = save_case(fig, "case14_mantel_feature_correlation")
    print(f"  saved -> {out}")
    plt.close(fig)


if __name__ == "__main__":
    main()
