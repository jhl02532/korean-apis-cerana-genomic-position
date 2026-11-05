"""Supplementary Fig. S1. Reference individuals on the ancestry index and silhouette calibration.
Run from the repository root: python figures/make_figS1.py"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt
import os; os.makedirs("figures/out", exist_ok=True)
from scipy.stats import spearmanr
C = "#1f77b4"
# ---- S3
COL = pd.read_csv("results/two_axes_colonies.tsv", sep="\t")
AN = pd.read_csv("results/two_axes_anchors.tsv", sep="\t"); AN["g"] = AN.group.str.split(":").str[0]
AN.loc[AN.group.isin(["China:Jilin", "China:Liaoning"]), "g"] = "NE China"
SIL = pd.read_csv("results/silhouette_fst_calibration.tsv", sep="\t")
GC = {"colonies": "#7f7f7f", "China": "#D55E00", "NE China": "#F4A582", "Japan": "#E6AB02"}
LAB = {"colonies": "This study\n(19 colonies)", "China": "Chinese reference,\noutside the northeast", "NE China": "Northeast Chinese\nreference", "Japan": "Japanese reference"}
rng = np.random.default_rng(1)
fig = plt.figure(figsize=(6.7, 3.2), dpi=300); gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1], wspace=0.42, left=0.215, right=0.985, top=0.9, bottom=0.17)
for k, (axis, rows, xlab, title) in enumerate([("KC", ["colonies", "Japan", "NE China", "China"], "Ancestry index (35,825 sites)\n← northeast Asian    Chinese mainland →", "Diagnostic sites")]):
    ax = fig.add_subplot(gs[k])
    for i, gname in enumerate(rows[::-1]):
        v = COL[axis].values if gname == "colonies" else AN[AN.g == gname][axis].dropna().values
        ax.scatter(v, i + rng.uniform(-0.12, 0.12, len(v)), s=22, color=GC[gname], edgecolor="white", lw=0.4, zorder=3)
        ax.text(1.01, i, f"n = {len(v)}", transform=ax.get_yaxis_transform(), fontsize=6.5, color="#555555", va="center")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([LAB[g] for g in rows[::-1]], fontsize=6.5)
    ax.set_xlabel(xlab, fontsize=7); ax.set_title(title, fontsize=7.5, loc="left"); ax.tick_params(axis="x", labelsize=6.5)
    ax.spines[["top", "right"]].set_visible(False)
ax = fig.add_subplot(gs[1])
ax.plot(SIL.F_ST, SIL.sil_mean, "o-", color="#0072B2", ms=2.5, lw=0.7, label="two discrete lineages\n(simulated, 1/18 split; mean)")
ax.fill_between(SIL.F_ST, SIL.sil_p05, SIL.sil_p95, color="#0072B2", alpha=0.5, lw=0, label="5–95% of simulations\n(narrower than the line)")
for y, l, c, xpos, ha, dy in [(0.25, "0.25 heuristic", "#333333", 0.004, "left", -0.006), (0.0979, "observed SNP 0.098", "#C1440E", 0.252, "right", 0.005),
                              (0.0749, "observed SV 0.075", "#C1440E", 0.252, "right", -0.006)]:
    ax.axhline(y, ls="--", color=c, lw=0.9); ax.text(xpos, y + dy, l, ha=ha, va="bottom" if dy > 0 else "top", fontsize=6.5, color=c)
ax.set_xlabel("Fst between simulated lineages", fontsize=7); ax.set_ylabel("Mean silhouette at $k$ = 2", fontsize=7)
ax.set_xlim(0, 0.26); ax.set_ylim(0, 0.37)
ax.legend(frameon=False, fontsize=6.5, loc="upper left", bbox_to_anchor=(0.0, 1.02)); ax.tick_params(labelsize=6.5); ax.spines[["top", "right"]].set_visible(False)
for a, l in zip(fig.axes, "AB"): a.text(-0.2 if l == "B" else -0.36, 1.04, l, transform=a.transAxes, fontsize=10, fontweight="bold")
fig.savefig("figures/out/fig_s1_calibration.png", dpi=300); print("ok")
