"""The 20 journal-figure recipes on real FAcDs data — no titles, one font scheme, one colour system.

Each `caseNN(host)` draws into `host` (a matplotlib Figure or SubFigure) and returns nothing, so the
same function serves a standalone PNG and a panel of the merged master figure (true vector in SVG).
Figure TYPES mirror GeneticistHere/ggplot2-20-journal-cases; DATA is this project's own.
"""
from __future__ import annotations
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch, Patch
from matplotlib.lines import Line2D
from matplotlib.path import Path as MPath
from matplotlib.colors import ListedColormap, BoundaryNorm
from scipy import stats as sps
from scipy.cluster.hierarchy import linkage, leaves_list

import dataio as D
from theme_case import (apply_theme, style_panel, cbar, OKABE_ITO, SIG, TIER_ORDER, TIER_COLORS,
                        SEQ, DIV, FS_LABEL, FS_TICK, FS_TICK_DENSE, FS_LEGEND, FS_ANNOT)

warnings.filterwarnings("ignore")
apply_theme()
RNG = np.random.default_rng(42)
_TL = [t.replace("Tier_", "").replace("_Decoy", "-D") for t in TIER_ORDER]      # short tier labels


def _tiers_present(df):
    return [t for t in TIER_ORDER if (df.degrader_tier == t).any()]


# ── (a) 01 · Individualized error dot-plot ───────────────────────────────────────────────────
def case01(host):
    ax = host.subplots()
    g = D.per_ligand("iptm").head(20); y = np.arange(len(g))[::-1]
    colors = [SIG["up"] if d else SIG["down"] for d in g["is_degrader"]]
    ax.errorbar(g["mean"], y, xerr=g["se"] * 1.96, fmt="none", ecolor="#666666", lw=1.0, capsize=2, zorder=1)
    ax.scatter(g["mean"], y, c=colors, s=42, edgecolor="black", lw=0.5, zorder=2)
    ax.set_yticks(y); ax.set_yticklabels(g.index, fontsize=FS_TICK_DENSE)
    ax.set_xlabel("Interface pTM  (mean ± 95% CI)"); ax.set_ylabel("Ligand")
    ax.legend(handles=[Line2D([0],[0],marker="o",ls="",mfc=SIG["up"],mec="k",label="degrader"),
                       Line2D([0],[0],marker="o",ls="",mfc=SIG["down"],mec="k",label="non-degrader")],
              loc="lower right")
    style_panel(ax)


# ── (b) 02 · Grouped error dot-plot ──────────────────────────────────────────────────────────
def case02(host):
    ax = host.subplots(); df = D.ranked(); metric = "mechanistic_score_effective"
    sub = df.dropna(subset=[metric])
    grp = sub.groupby(["degrader_tier", "is_degrader"], observed=True)[metric].agg(["mean", "std", "count"])
    grp["se"] = grp["std"] / np.sqrt(grp["count"].clip(lower=1))
    tiers = _tiers_present(df); x = np.arange(len(tiers))
    for deg in [True, False]:
        m = [grp.loc[(t, deg), "mean"] if (t, deg) in grp.index else np.nan for t in tiers]
        e = [grp.loc[(t, deg), "se"] * 1.96 if (t, deg) in grp.index else 0 for t in tiers]
        off = (-0.11 if deg else 0.11); c = SIG["up"] if deg else SIG["down"]
        ax.errorbar(x + off, m, yerr=e, fmt="o", ms=6, color=c, ecolor=c, lw=1.1, capsize=3,
                    mec="black", mew=0.5, label="degrader" if deg else "non-degrader")
    ax.set_xticks(x); ax.set_xticklabels([t.replace("Tier_", "").replace("_Decoy", "-D") for t in tiers])
    ax.set_xlabel("Degrader tier"); ax.set_ylabel("Mechanistic score  (mean ± 95% CI)")
    ax.legend(loc="upper right"); style_panel(ax)


# ── (c) 03 · Multi-group volcano ─────────────────────────────────────────────────────────────
def case03(host):
    ax = host.subplots(); df = D.ranked(); deg = df["is_degrader"] == True
    rows = []
    for c in D.NUMERIC_METRICS:
        if c not in df: continue
        a = pd.to_numeric(df.loc[deg, c], errors="coerce").dropna()
        b = pd.to_numeric(df.loc[~deg, c], errors="coerce").dropna()
        if len(a) < 30 or len(b) < 30: continue
        _, p = sps.mannwhitneyu(a, b, alternative="two-sided")
        pooled = np.sqrt(((a.std()**2)+(b.std()**2))/2) or 1e-9
        rows.append((c, (a.mean()-b.mean())/pooled, p))
    v = pd.DataFrame(rows, columns=["metric", "effect", "p"])
    v["nlp"] = -np.log10(v["p"].clip(lower=1e-300))
    v["col"] = np.where((v.p < 0.05) & (v.effect > 0.1), SIG["up"],
                np.where((v.p < 0.05) & (v.effect < -0.1), SIG["down"], SIG["ns"]))
    ax.scatter(v.effect, v.nlp, c=v.col, s=48, edgecolor="black", lw=0.5)
    ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color="#888888")
    ax.axvline(0, ls="-", lw=0.5, color="#888888")
    top = v.sort_values("nlp", ascending=False).head(6)
    ys = np.linspace(v.nlp.max()*0.97, v.nlp.max()*0.55, len(top))
    for (_, r), yy in zip(top.iterrows(), ys):
        ax.annotate(D.mlabel(r.metric), (r.effect, r.nlp), (r.effect, yy), fontsize=FS_ANNOT,
                    ha="left", arrowprops=dict(arrowstyle="-", lw=0.3, color="#999999"))
    ax.set_xlabel("Effect size  (Cohen's d)"); ax.set_ylabel(r"$-\log_{10}$ p  (Mann–Whitney)")
    style_panel(ax)


# ── (d) 04 · Manhattan / TWAS ────────────────────────────────────────────────────────────────
def case04(host):
    ax = host.subplots(); df = D.ranked().copy()
    df["nlp"] = -np.log10((df["Scientific_Rank"].rank(pct=True)).clip(lower=1e-6))
    df = df.sort_values(["ligand_short", "Scientific_Rank"])
    ligs = list(df["ligand_short"].unique()); xpos = 0; ticks = []
    for i, lg in enumerate(ligs):
        d = df[df.ligand_short == lg]; x = np.arange(xpos, xpos + len(d))
        ax.scatter(x, d["nlp"], s=3, color=OKABE_ITO[i % len(OKABE_ITO)], alpha=0.7, edgecolor="none")
        ticks.append(xpos + len(d)/2); xpos += len(d) + 40
    ax.axhline(-np.log10(0.05), ls="--", lw=0.7, color=SIG["up"])
    ax.set_xticks(ticks); ax.set_xticklabels(ligs, rotation=90, fontsize=FS_TICK_DENSE)
    ax.set_ylabel(r"$-\log_{10}$ p  (rank-based)"); ax.set_xlabel("Ligand locus")
    ax.margins(x=0.01); style_panel(ax)


# ── (e) 05 · Paired boxplot ──────────────────────────────────────────────────────────────────
def case05(host):
    ax = host.subplots(); df = D.ranked().dropna(subset=["iptm", "ptm"])
    s = df.sample(min(120, len(df)), random_state=1)
    bp = ax.boxplot([df["iptm"].values, df["ptm"].values], positions=[1, 2], widths=0.45,
                    patch_artist=True, showfliers=False, medianprops=dict(color="black"))
    for patch, c in zip(bp["boxes"], [SIG["down"], SIG["up"]]):
        patch.set_facecolor(c); patch.set_alpha(0.5); patch.set_edgecolor("black")
    for _, r in s.iterrows():
        ax.plot([1, 2], [r["iptm"], r["ptm"]], color="#BBBBBB", lw=0.4, alpha=0.5, zorder=0)
    jx = RNG.normal(0, 0.03, len(s))
    ax.scatter(1 + jx, s["iptm"], s=7, color=SIG["down"], edgecolor="none", alpha=0.6, zorder=3)
    ax.scatter(2 + jx, s["ptm"], s=7, color=SIG["up"], edgecolor="none", alpha=0.6, zorder=3)
    ax.set_xticks([1, 2]); ax.set_xticklabels(["iPTM", "pTM"]); ax.set_ylabel("Score"); style_panel(ax)


# ── (f) 06 · Raincloud ───────────────────────────────────────────────────────────────────────
def case06(host):
    ax = host.subplots(); df = D.ranked(); metric = "mechanistic_score_effective"
    tiers = _tiers_present(df)
    for i, t in enumerate(tiers):
        vals = pd.to_numeric(df.loc[df.degrader_tier == t, metric], errors="coerce").dropna().values
        if len(vals) < 5: continue
        c = TIER_COLORS[t]
        kde = sps.gaussian_kde(vals); ys = np.linspace(vals.min(), vals.max(), 120); dens = kde(ys)
        dens = dens / dens.max() * 0.35
        ax.fill_betweenx(ys, i, i + dens, color=c, alpha=0.55, lw=0)
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        ax.add_patch(Rectangle((i-0.06, q1), 0.12, q3-q1, facecolor="white", edgecolor="black", lw=0.7, zorder=3))
        ax.plot([i-0.06, i+0.06], [med, med], color="black", lw=1.2, zorder=4)
        samp = RNG.choice(vals, min(250, len(vals)), replace=False)
        ax.scatter(i-0.18 + RNG.normal(0, 0.04, len(samp)), samp, s=3, color=c, alpha=0.3, edgecolor="none")
    ax.set_xticks(range(len(tiers))); ax.set_xticklabels([t.replace("Tier_", "").replace("_Decoy", "-D") for t in tiers])
    ax.set_xlabel("Degrader tier"); ax.set_ylabel("Mechanistic score (effective)"); style_panel(ax)


# ── (g) 07 · Swimmer plot (top candidates per tier, Tier_1A → Tier_5_Decoy) ──────────────────
def case07(host):
    ax = host.subplots(); df = D.ranked(); tiers = _tiers_present(df)
    picks = [df[df.degrader_tier == t].sort_values("Scientific_Rank").head(3) for t in tiers]
    d = pd.concat(picks).copy()                                    # 3 per tier, ordered 1A → decoy
    d["len"] = pd.to_numeric(d["mechanistic_score_effective"], errors="coerce").fillna(0)
    d = d.iloc[::-1]; y = np.arange(len(d))                        # reversed so Tier_1A sits at the top
    tcol = [TIER_COLORS.get(str(t), "#999999") for t in d["degrader_tier"]]
    ax.barh(y, d["len"], color=tcol, edgecolor="black", lw=0.4, height=0.7)
    md = (d["MD_Selected"] == True).values
    ax.scatter(d["len"].values[md] + 0.012, y[md], marker=">", s=34, color="black", zorder=3, label="MD-selected")
    deg = (d["is_degrader"] == True).values
    ax.scatter(np.full(deg.sum(), -0.012), y[deg], marker="o", s=18, color=SIG["up"], zorder=3, label="degrader")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{str(t).replace('Tier_','').replace('_Decoy','-D')} · {lg}"
                        for t, lg in zip(d["degrader_tier"], d["ligand_short"])], fontsize=FS_ANNOT)
    ax.set_xlabel("Mechanistic score (effective)"); ax.set_ylabel("Tier · ligand  (top 3 per tier)")
    # tier is already named on every y-tick and shown by bar colour, so the legend only decodes the markers
    handles = [Line2D([0],[0],marker=">",ls="",mfc="black",mec="black",label="MD-selected"),
               Line2D([0],[0],marker="o",ls="",mfc=SIG["up"],mec=SIG["up"],label="degrader")]
    leg = ax.legend(handles=handles, loc="lower right", fontsize=FS_ANNOT, handlelength=1,
                    frameon=True, framealpha=0.85, edgecolor="none")
    leg.get_frame().set_facecolor("white")
    ax.margins(y=0.01); style_panel(ax)


# ── (h) 08 · Importance bars + abundance streams ─────────────────────────────────────────────
def case08(host):
    axL, axR = host.subplots(1, 2, gridspec_kw={"width_ratios": [1, 1.5], "wspace": 0.5})
    imp = D.pca_loadings()["PC1"].abs().sort_values(); df = D.ranked()
    axL.barh(range(len(imp)), imp.values, color=SIG["down"], edgecolor="black", lw=0.4)
    axL.set_yticks(range(len(imp))); axL.set_yticklabels([D.mlabel(c) for c in imp.index], fontsize=FS_ANNOT)
    axL.set_xlabel("|PC1 loading|"); style_panel(axL)
    order = D.per_ligand("iptm").index.tolist()
    ct = (df.groupby(["ligand_short", "degrader_tier"], observed=True).size().unstack(fill_value=0)
          .reindex(order).fillna(0))
    ct = ct[[t for t in TIER_ORDER if t in ct.columns]]; x = np.arange(len(ct)); base = -ct.sum(axis=1).values/2
    for t in ct.columns:
        axR.fill_between(x, base, base + ct[t].values, color=TIER_COLORS[t],
                         label=t.replace("Tier_", "").replace("_Decoy", "-D"), lw=0); base = base + ct[t].values
    axR.set_xticks(x); axR.set_xticklabels(ct.index, rotation=90, fontsize=FS_ANNOT)
    axR.set_ylabel("Candidates (stream)"); axR.set_xlabel("Ligand")
    axR.legend(ncol=4, fontsize=FS_ANNOT, loc="upper center", handlelength=1); axR.grid(False)
    for s in axR.spines.values(): s.set_visible(False)


# ── (i) 09 · Grouped variance bars with letters ──────────────────────────────────────────────
def case09(host):
    ax = host.subplots(); df = D.ranked(); metric = "iptm"; tiers = _tiers_present(df)
    groups = {t: pd.to_numeric(df.loc[df.degrader_tier == t, metric], errors="coerce").dropna() for t in tiers}
    means = [groups[t].mean() for t in tiers]; sds = [groups[t].std() for t in tiers]
    order = sorted(tiers, key=lambda t: -groups[t].mean()); letters = {}; cur = ord("a")
    for t in order:
        placed = False
        for L in sorted(set(letters.values())) if letters else []:
            members = [k for k, v in letters.items() if v == L]
            if all(sps.mannwhitneyu(groups[t], groups[m])[1] > 0.05 for m in members):
                letters[t] = L; placed = True; break
        if not placed:
            letters[t] = chr(cur); cur += 1
    x = np.arange(len(tiers))
    ax.bar(x, means, yerr=sds, color=[TIER_COLORS[t] for t in tiers], edgecolor="black", lw=0.5,
           capsize=4, error_kw=dict(lw=1))
    for xi, t, m, s in zip(x, tiers, means, sds):
        ax.text(xi, m + (s or 0) + 0.012, letters[t], ha="center", fontsize=10, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels([t.replace("Tier_", "").replace("_Decoy", "-D") for t in tiers])
    ax.set_xlabel("Degrader tier"); ax.set_ylabel("iPTM  (mean ± SD)"); ax.set_ylim(0, max(means)+max(sds)+0.08)
    style_panel(ax)


# ── (j) 10 · Circos / chord ──────────────────────────────────────────────────────────────────
def _bezier(ax, a0, a1, color, lw):
    p0 = np.array([np.cos(a0), np.sin(a0)]); p1 = np.array([np.cos(a1), np.sin(a1)])
    path = MPath([p0, p0*0.15, p1*0.15, p1], [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4])
    ax.add_patch(PathPatch(path, fc="none", ec=color, lw=lw, alpha=0.5))


def case10(host):
    ax = host.subplots(); ax.set_aspect("equal"); ax.axis("off")
    df = D.ranked()
    ct = df.groupby(["ligand_short", "degrader_tier"], observed=True).size().unstack(fill_value=0)
    ct = ct[[t for t in TIER_ORDER if t in ct.columns]]
    ligs = ct.index.tolist(); tiers = ct.columns.tolist(); nodes = ligs + tiers; N = len(nodes)
    ang = {n: 2*np.pi*i/N for i, n in enumerate(nodes)}
    for n, a in ang.items():
        is_t = n in TIER_COLORS; col = TIER_COLORS.get(n, "#9AA0A6")
        ax.plot([0.97*np.cos(a), 1.03*np.cos(a)], [0.97*np.sin(a), 1.03*np.sin(a)],
                color=col, lw=8 if is_t else 3.5, solid_capstyle="butt")
        deg = np.degrees(a) % 360; rot = np.degrees(a); rot = rot-180 if 90 < deg < 270 else rot
        ax.text(1.08*np.cos(a), 1.08*np.sin(a), n.replace("Tier_", "T").replace("_Decoy", "5-D"),
                ha="left" if deg <= 90 or deg >= 270 else "right", va="center",
                fontsize=FS_ANNOT, rotation=rot, rotation_mode="anchor",
                fontweight="bold" if is_t else "normal")
    mx = ct.values.max()
    for lg in ligs:
        for t in tiers:
            w = ct.loc[lg, t]
            if w > 0: _bezier(ax, ang[lg], ang[t], TIER_COLORS[t], 0.3 + 3.0*w/mx)
    ax.set_xlim(-1.35, 1.35); ax.set_ylim(-1.35, 1.35)


# ── (k) 11 · Module interaction network ──────────────────────────────────────────────────────
def case11(host):
    import networkx as nx
    ax = host.subplots(); ax.axis("off")
    m = D.metric_matrix().sample(min(5000, len(D.metric_matrix())), random_state=3); corr = m.corr().abs()
    G = nx.Graph(); cols = corr.columns
    for i in range(len(cols)):
        for j in range(i+1, len(cols)):
            if corr.iloc[i, j] >= 0.5: G.add_edge(cols[i], cols[j], w=corr.iloc[i, j])
    G.add_nodes_from(cols); pos = nx.spring_layout(G, seed=7, k=0.9)
    comp = list(nx.connected_components(G))
    cmap = {n: OKABE_ITO[i % len(OKABE_ITO)] for i, c in enumerate(comp) for n in c}
    for u, v, d in G.edges(data=True):
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#CCCCCC",
                lw=0.5 + 2.2*(d["w"]-0.5), alpha=0.6, zorder=1)
    deg = dict(G.degree())
    for n, p in pos.items():
        ax.scatter(*p, s=90 + 70*deg.get(n, 0), color=cmap.get(n, "#999999"), edgecolor="black", lw=0.6, zorder=2)
        ax.text(p[0], p[1], D.mlabel(n), fontsize=FS_ANNOT, ha="center", va="center", zorder=3)
    ax.margins(0.08)


# ── (l) 12 · Sankey with enrichment bubble ───────────────────────────────────────────────────
def _sankey_two(ax, left, right, M, lcol, rcol):
    Ln, Rn = M.shape; ltot = M.sum(1); rtot = M.sum(0); tot = M.sum(); gap = 0.02
    lpos = {}; y = 1.0
    for i in range(Ln):
        h = ltot[i]/tot*(1-gap*(Ln-1)); lpos[i] = (y-h, y); y -= h+gap
    rpos = {}; y = 1.0
    for j in range(Rn):
        h = rtot[j]/tot*(1-gap*(Rn-1)); rpos[j] = (y-h, y); y -= h+gap
    lc = {i: lpos[i][1] for i in range(Ln)}; rc = {j: rpos[j][1] for j in range(Rn)}
    for i in range(Ln):
        for j in range(Rn):
            f = M[i, j]
            if f <= 0: continue
            h = f/tot; ys = np.linspace(0, 1, 40); sm = 3*ys**2-2*ys**3; xs = 0.14+0.72*ys
            top = lc[i] + (rc[j]-lc[i])*sm; bot = (lc[i]-h) + ((rc[j]-h)-(lc[i]-h))*sm
            ax.fill_between(xs, bot, top, color=lcol[i], alpha=0.45, lw=0); lc[i] -= h; rc[j] -= h
    for i in range(Ln):
        ax.add_patch(Rectangle((0.10, lpos[i][0]), 0.03, lpos[i][1]-lpos[i][0], color=lcol[i]))
        ax.text(0.085, sum(lpos[i])/2, left[i], ha="right", va="center", fontsize=FS_ANNOT)
    for j in range(Rn):
        ax.add_patch(Rectangle((0.87, rpos[j][0]), 0.03, rpos[j][1]-rpos[j][0], color=rcol[j]))
        ax.text(0.915, sum(rpos[j])/2, right[j], ha="left", va="center", fontsize=FS_ANNOT)
    ax.set_xlim(0, 1); ax.set_ylim(-0.02, 1.02); ax.axis("off")


def case12(host):
    axS, axB = host.subplots(1, 2, gridspec_kw={"width_ratios": [2, 1], "wspace": 0.35})
    df = D.ranked(); order = D.per_ligand("iptm").index.tolist()[:12]; sub = df[df.ligand_short.isin(order)]
    ct = sub.groupby(["ligand_short", "degrader_tier"], observed=True).size().unstack(fill_value=0)
    ct = ct.reindex(order)[[t for t in TIER_ORDER if t in ct.columns]]
    lcol = [OKABE_ITO[i % len(OKABE_ITO)] for i in range(len(ct.index))]; rcol = [TIER_COLORS[t] for t in ct.columns]
    _sankey_two(axS, list(ct.index), [t.replace("Tier_", "T").replace("_Decoy", "5-D") for t in ct.columns],
                ct.values, lcol, rcol)
    enr = df.groupby("degrader_tier", observed=True).agg(n=("iptm", "size"), frac=("is_degrader", "mean")).reindex(
        _tiers_present(df))
    yb = np.arange(len(enr))[::-1]
    axB.scatter(enr["frac"], yb, s=enr["n"]/enr["n"].max()*700+30, c=[TIER_COLORS[t] for t in enr.index],
                edgecolor="black", lw=0.6, alpha=0.85)
    axB.set_yticks(yb); axB.set_yticklabels([t.replace("Tier_", "").replace("_Decoy", "-D") for t in enr.index])
    axB.set_xlabel("Degrader fraction"); axB.set_ylabel("Tier"); style_panel(axB)


# ── (m) 13 · Discrete (binned) heatmap ───────────────────────────────────────────────────────
def case13(host):
    ax = host.subplots(); df = D.ranked()
    metrics = ["iptm", "mean_plddt", "interaction_density", "mechanistic_score_effective",
               "competence_score", "halide_stabilisation_score", "pocket_occupancy",
               "num_interactions", "count_hydrogen_bond"]
    order = D.per_ligand("iptm").index.tolist()
    agg = df.groupby("ligand_short", observed=True)[metrics].mean().reindex(order)
    z = agg.apply(lambda c: (c-c.min())/(c.max()-c.min()+1e-9)); binned = np.clip((z*5).astype(int), 0, 4)
    cmap = ListedColormap(plt.get_cmap(SEQ)(np.linspace(0.1, 0.95, 5)))
    im = ax.imshow(binned.values, cmap=cmap, norm=BoundaryNorm([-.5, .5, 1.5, 2.5, 3.5, 4.5], 5), aspect="auto")
    ax.set_xticks(range(len(metrics))); ax.set_xticklabels([D.mlabel(m) for m in metrics], rotation=45, ha="right", fontsize=FS_TICK_DENSE)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order, fontsize=FS_TICK_DENSE)
    cbar(host, im, ax, "Quintile band", ticks=[0, 1, 2, 3, 4], shrink=0.6)


# ── (n) 14 · Mantel composite heatmap ────────────────────────────────────────────────────────
def case14(host):
    ax = host.subplots(); m = D.metric_matrix().sample(min(6000, len(D.metric_matrix())), random_state=5)
    corr = m.corr(); n = corr.shape[0]; pv = np.ones((n, n))
    for i in range(n):
        for j in range(i+1, n):
            _, p = sps.spearmanr(m.iloc[:, i], m.iloc[:, j], nan_policy="omit"); pv[i, j] = pv[j, i] = p
    up = np.triu(np.ones((n, n), bool), 1)
    A = np.where(up, corr.values, np.nan); B = np.where(~up, -np.log10(np.clip(pv, 1e-50, 1)), np.nan)
    np.fill_diagonal(B, np.nan)
    im1 = ax.imshow(A, cmap=DIV, vmin=-1, vmax=1); im2 = ax.imshow(B, cmap=SEQ)
    ax.set_xticks(range(n)); ax.set_xticklabels([D.mlabel(c) for c in corr.columns], rotation=90, fontsize=FS_ANNOT)
    ax.set_yticks(range(n)); ax.set_yticklabels([D.mlabel(c) for c in corr.columns], fontsize=FS_ANNOT)
    cbar(host, im1, ax, "Pearson r (upper)", fraction=0.045, pad=0.02)
    cbar(host, im2, ax, r"$-\log_{10}$p (lower)", fraction=0.045, pad=0.10)


# ── (o) 15 · Multi-group correlation heatmap ─────────────────────────────────────────────────
def case15(host):
    ax = host.subplots(); m = D.metric_matrix().sample(min(8000, len(D.metric_matrix())), random_state=9)
    corr = m.corr(); idx = leaves_list(linkage(1 - corr.abs(), method="average")); corr = corr.iloc[idx, idx]
    im = ax.imshow(corr.values, cmap=DIV, vmin=-1, vmax=1)
    ax.set_xticks(range(len(corr))); ax.set_xticklabels([D.mlabel(c) for c in corr.columns], rotation=90, fontsize=FS_ANNOT)
    ax.set_yticks(range(len(corr))); ax.set_yticklabels([D.mlabel(c) for c in corr.index], fontsize=FS_ANNOT)
    for i in range(len(corr)):
        for j in range(len(corr)):
            if abs(corr.values[i, j]) >= 0.6 and i != j:
                ax.text(j, i, f"{corr.values[i,j]:.1f}", ha="center", va="center", fontsize=4.5,
                        color="white" if abs(corr.values[i, j]) > 0.8 else "black")
    cbar(host, im, ax, "Pearson r")


# ── (p) 16 · Polar-coordinate heatmap ────────────────────────────────────────────────────────
def case16(host):
    ax = host.add_subplot(projection="polar"); df = D.ranked()
    metrics = ["iptm", "mean_plddt", "interaction_density", "competence_score",
               "mechanistic_score_effective", "halide_stabilisation_score", "pocket_occupancy", "num_interactions"]
    order = D.per_ligand("iptm").index.tolist()
    agg = df.groupby("ligand_short", observed=True)[metrics].mean().reindex(order)
    z = agg.apply(lambda c: (c-c.min())/(c.max()-c.min()+1e-9)); nL, nM = z.shape
    theta = np.linspace(0, 2*np.pi, nL, endpoint=False); dt = 2*np.pi/nL; cmap = plt.get_cmap(SEQ)
    for ri in range(nM):
        for li in range(nL):
            ax.bar(theta[li], 1, width=dt*0.98, bottom=ri, color=cmap(z.iloc[li, ri]),
                   edgecolor="white", lw=0.3, align="edge")
    ax.set_xticks(theta + dt/2); ax.set_xticklabels(order, fontsize=FS_ANNOT)
    ax.set_yticks(np.arange(nM) + 0.5); ax.set_yticklabels([D.mlabel(m) for m in metrics], fontsize=FS_ANNOT); ax.set_ylim(0, nM)
    sm = plt.cm.ScalarMappable(cmap=cmap); sm.set_array([]); cbar(host, sm, ax, "min–max scaled", pad=0.12)


# ── (q) 17 · Multi-level Sankey ──────────────────────────────────────────────────────────────
def _stage(ax, x0, x1, M, lcol, tot, gap=0.015):
    Ln, Rn = M.shape; ltot = M.sum(1); rtot = M.sum(0)
    lpos = {}; y = 1.0
    for i in range(Ln):
        h = ltot[i]/tot; lpos[i] = (y-h, y); y -= h+gap
    rpos = {}; y = 1.0
    for j in range(Rn):
        h = rtot[j]/tot; rpos[j] = (y-h, y); y -= h+gap
    lc = {i: lpos[i][1] for i in range(Ln)}; rc = {j: rpos[j][1] for j in range(Rn)}
    for i in range(Ln):
        for j in range(Rn):
            f = M[i, j]
            if f <= 0: continue
            h = f/tot; ys = np.linspace(0, 1, 30); sm = 3*ys**2-2*ys**3; xs = x0+(x1-x0)*ys
            ax.fill_between(xs, (lc[i]-h)+((rc[j]-h)-(lc[i]-h))*sm, lc[i]+(rc[j]-lc[i])*sm,
                            color=lcol[i], alpha=0.4, lw=0); lc[i] -= h; rc[j] -= h
    return lpos, rpos


def case17(host):
    ax = host.subplots(); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(-0.03, 1.03)
    df = D.ranked().copy(); df["class"] = np.where(df.is_degrader, "degrader", "non-degr.")
    df["md"] = np.where(df.MD_Selected == True, "MD-sel", "not-sel"); tiers = _tiers_present(df); tot = len(df)
    A = df.groupby(["class", "degrader_tier"], observed=True).size().unstack(fill_value=0).reindex(
        ["degrader", "non-degr."])[tiers].values.astype(float)
    Bm = df.groupby(["degrader_tier", "md"], observed=True).size().unstack(fill_value=0).reindex(
        tiers).fillna(0)[["MD-sel", "not-sel"]].values.astype(float)
    ccol = [SIG["up"], SIG["down"]]; tcol = [TIER_COLORS[t] for t in tiers]; mcol = ["#009E73", "#9AA0A6"]
    lp1, rp1 = _stage(ax, 0.12, 0.44, A, ccol, tot); _stage(ax, 0.56, 0.88, Bm, tcol, tot)
    for i, (l, c) in enumerate(zip(["degrader", "non-degr."], ccol)):
        ax.add_patch(Rectangle((0.095, lp1[i][0]), 0.02, lp1[i][1]-lp1[i][0], color=c))
        ax.text(0.09, sum(lp1[i])/2, l, ha="right", va="center", fontsize=FS_ANNOT)
    for j, t in enumerate(tiers):
        ax.add_patch(Rectangle((0.445, rp1[j][0]), 0.02, rp1[j][1]-rp1[j][0], color=tcol[j]))
        ax.text(0.5, (rp1[j][0]+rp1[j][1])/2, t.replace("Tier_", "T").replace("_Decoy", "5-D"),
                ha="center", va="center", fontsize=FS_ANNOT)
    y = 1.0
    for j, (l, c) in enumerate(zip(["MD-sel", "not-sel"], mcol)):
        h = Bm.sum(0)[j]/tot; ax.add_patch(Rectangle((0.885, y-h), 0.02, h, color=c))
        ax.text(0.91, y-h/2, l, ha="left", va="center", fontsize=FS_ANNOT); y -= h+0.015
    for xx, lab in [(0.12, "class"), (0.5, "tier"), (0.88, "MD")]:
        ax.text(xx, -0.02, lab, ha="center", fontsize=FS_TICK, color="#555555")


# ── (r) 18 · Treemap ─────────────────────────────────────────────────────────────────────────
def case18(host):
    from theme_case import squarify
    ax = host.subplots(); ax.axis("off"); df = D.ranked()
    cnt = df.groupby(["degrader_tier", "ligand_short"], observed=True).size().reset_index(name="n")
    cnt = cnt[cnt.n > 0].sort_values("n", ascending=False)
    for r, (_, row) in zip(squarify(cnt["n"].values, 0, 0, 100, 100), cnt.iterrows()):
        c = TIER_COLORS.get(str(row["degrader_tier"]), "#999999")
        ax.add_patch(Rectangle((r["x"], r["y"]), r["dx"], r["dy"], facecolor=c, edgecolor="white", lw=0.8))
        if r["dx"]*r["dy"] > 42:
            ax.text(r["x"]+r["dx"]/2, r["y"]+r["dy"]/2, f"{row['ligand_short']}\n{int(row['n'])}",
                    ha="center", va="center", fontsize=FS_ANNOT,
                    color="white" if c not in ("#E6C700",) else "black")
    ax.set_xlim(0, 100); ax.set_ylim(0, 100)
    ax.legend(handles=[Patch(color=TIER_COLORS[t], label=t.replace("Tier_", "T").replace("_Decoy", "5-D"))
                       for t in TIER_ORDER], ncol=7, fontsize=FS_ANNOT, loc="upper center",
              bbox_to_anchor=(0.5, -0.01), handlelength=1)


# ── (s) 19 · Mosaic and sunburst ─────────────────────────────────────────────────────────────
def case19(host):
    df = D.ranked(); tiers = _tiers_present(df)
    axM = host.add_subplot(1, 2, 1); axS = host.add_subplot(1, 2, 2, projection="polar")
    ct = df.groupby("degrader_tier", observed=True).size().reindex(tiers).fillna(0)
    fracdeg = df.groupby("degrader_tier", observed=True)["is_degrader"].mean().reindex(tiers).fillna(0)
    x = 0; total = ct.sum()
    for t in tiers:
        w = ct[t]/total
        axM.add_patch(Rectangle((x, 0), w, fracdeg[t], facecolor=SIG["up"], edgecolor="white", lw=0.8))
        axM.add_patch(Rectangle((x, fracdeg[t]), w, 1-fracdeg[t], facecolor=SIG["down"], edgecolor="white", lw=0.8))
        axM.text(x+w/2, -0.05, t.replace("Tier_", "T").replace("_Decoy", "5-D"), ha="center", fontsize=FS_ANNOT)
        x += w
    axM.set_xlim(0, 1); axM.set_ylim(0, 1); axM.axis("off")
    axM.legend(handles=[Patch(color=SIG["up"], label="degrader"), Patch(color=SIG["down"], label="non-degrader")],
               loc="upper center", bbox_to_anchor=(0.5, -0.04), ncol=1, fontsize=FS_ANNOT)
    axS.axis("off"); ct2 = df.groupby(["degrader_tier", "ligand_short"], observed=True).size()
    tt = ct2.groupby(level=0).sum().reindex(tiers).fillna(0); total = tt.sum(); a0 = 0.0
    for t in tiers:
        span = tt[t]/total*2*np.pi
        axS.bar(a0+span/2, 0.5, width=span, bottom=0.5, color=TIER_COLORS[t], edgecolor="white", lw=0.8)
        aa = a0
        for lg, n in ct2.loc[t].sort_values(ascending=False).items():
            s = n/total*2*np.pi
            axS.bar(aa+s/2, 0.45, width=s*0.98, bottom=1.0, color=TIER_COLORS[t], edgecolor="white", lw=0.3, alpha=0.7)
            aa += s
        a0 += span
    axS.set_ylim(0, 1.5)


# ── (t) 20 · Split violin ────────────────────────────────────────────────────────────────────
def case20(host):
    ax = host.subplots(); df = D.ranked(); metric = "iptm"; tiers = _tiers_present(df)
    for i, t in enumerate(tiers):
        for side, deg, col in [(-1, True, SIG["up"]), (1, False, SIG["down"])]:
            vals = pd.to_numeric(df.loc[(df.degrader_tier == t) & (df.is_degrader == deg), metric],
                                 errors="coerce").dropna().values
            if len(vals) < 5: continue
            kde = sps.gaussian_kde(vals); ys = np.linspace(vals.min(), vals.max(), 100)
            dens = kde(ys); dens = dens/dens.max()*0.4
            ax.fill_betweenx(ys, i, i + side*dens, color=col, alpha=0.6, lw=0.4, edgecolor="black")
            ax.scatter(i, np.median(vals), color="white", edgecolor="black", zorder=3, s=13)
    ax.set_xticks(range(len(tiers))); ax.set_xticklabels([t.replace("Tier_", "").replace("_Decoy", "-D") for t in tiers])
    ax.set_xlabel("Degrader tier"); ax.set_ylabel("iPTM")
    ax.legend(handles=[Patch(color=SIG["up"], label="degrader"), Patch(color=SIG["down"], label="non-degrader")],
              loc="lower right"); style_panel(ax)


ALL = [case01, case02, case03, case04, case05, case06, case07, case08, case09, case10,
       case11, case12, case13, case14, case15, case16, case17, case18, case19, case20]
# preferred standalone aspect (w, h) inches — merged mode ignores this
SIZE = {0: (6, 5.5), 1: (6.5, 5), 2: (6.5, 5.5), 3: (8, 4.5), 4: (5.5, 5.5), 5: (7.5, 5),
        6: (7.5, 6.5), 7: (10, 5), 8: (6.5, 5), 9: (7, 7), 10: (8, 7), 11: (11, 6),
        12: (7.5, 8), 13: (8.5, 8), 14: (8.5, 8), 15: (8, 8), 16: (10, 6), 17: (9, 6.5),
        18: (12, 6.5), 19: (8, 5)}
