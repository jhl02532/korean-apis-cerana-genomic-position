"""Supplementary Fig. S3. Gene content: genes by number of colonies lacking them, and chemoreceptor families.
Run from the repository root: python figures/make_figS3.py"""
import pandas as pd
import os; os.makedirs("figures/out", exist_ok=True)
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams.update({'font.size': 6.5, 'axes.titlesize': 7.5, 'axes.labelsize': 7, 'xtick.labelsize': 6.5,
                     'ytick.labelsize': 6.5, 'legend.fontsize': 6.5})
META_GREY = "#7a7a7a"; MUT = "#555555"; BLUE = "#0072B2"; RES = "#C1440E"

def panel_letter(ax, letter, dx=-0.13):
    ax.text(dx, 1.03, letter, transform=ax.transAxes, fontsize=10, fontweight='bold', va='bottom', ha='right')

# ---- A 자료: 유전자별 결손 군락 수 (읽기 깊이 < 전장 평균 10% = 결손)
GP = pd.read_csv('results/rs_gene_presence_read.tsv', sep="\t").set_index('g'); assert GP.shape == (9884, 19)
n_abs = (GP == 0).sum(axis=1); per_col = (GP == 0).sum(axis=0)
BINS = [("0", n_abs == 0), ("1", n_abs == 1), ("2–6", n_abs.between(2, 6)), ("7–19", n_abs >= 7)]
cnt = [int(m.sum()) for _, m in BINS]; assert cnt == [9875, 2, 7, 0], cnt
assert per_col.median() == 1 and per_col.min() == 0 and per_col.max() == 5

FAMV = [("OR",0,93),("GR",0,6),("OBP",0,12),("IR",0,11),("CSP",0,5),("background",9,9757)]   # 1개 이상 군락에서 결손
names = [f for f,_,_ in FAMV]; pcts = [100*k/n for _,k,n in FAMV]
cols = [BLUE]*5 + [META_GREY]

fig = plt.figure(figsize=(6.7, 2.6))
gt = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.42, left=0.1, right=0.97, top=0.88, bottom=0.18)

# A
ax = fig.add_subplot(gt[0]); xb = np.arange(len(BINS))
ax.bar(xb, [max(c, 0) for c in cnt], color=[BLUE, RES, RES, RES], width=0.6, edgecolor="white", lw=0.5, zorder=3)
for i, c in enumerate(cnt):
    ax.text(i, max(c, 1) * 1.25, f"{c:,}", ha="center", va="bottom", fontsize=6.5)
ax.set_yscale("log"); ax.minorticks_off(); ax.set_ylim(0.6, 6e4); ax.set_yticks([1, 10, 100, 1000, 10000]); ax.set_yticklabels(["1", "10", "100", "1,000", "10,000"])
ax.set_xticks(xb); ax.set_xticklabels([b for b, _ in BINS])
ax.set_xlabel("Colonies in which the gene is absent (of 19)"); ax.set_ylabel("Coding genes (log scale)")
ax.text(0.97, 0.60, f"per colony: median {int(per_col.median())}\nabsent gene (range {per_col.min()}–{per_col.max()})",
        transform=ax.transAxes, ha="right", va="center", fontsize=6.5, color=MUT)
ax.spines[["top", "right"]].set_visible(False)
ax.set_title(f"{sum(cnt):,} chromosomal coding genes (read depth)", loc="left", fontsize=7)
panel_letter(ax, "A", dx=-0.2)

# B
ax = fig.add_subplot(gt[1])
ax.bar(names, pcts, color=cols, width=0.58, edgecolor="white", lw=0.5, zorder=3)
for i,(f,k,n) in enumerate(FAMV):
    v = 100*k/n
    if v == 0: ax.plot([i],[0.002], marker="o", ms=2.6, color=cols[i], zorder=4)
    ax.text(i, max(v,0.002)+0.004, f"{k}/{n:,}", ha="center", va="bottom", fontsize=6.5)
ax.set_ylabel("Genes absent in ≥1 colony (%)")
ax.set_ylim(0, 0.12)
ax.set_xticks(range(len(names))); ax.set_xticklabels(["OR","GR","OBP","IR","CSP","back-\nground"])
ax.text(0.03, 0.97, "no family gene absent\nfrom any colony", transform=ax.transAxes, ha="left", va="top", linespacing=1.25, color=MUT)
ax.spines[["top", "right"]].set_visible(False)
panel_letter(ax, "B", dx=-0.2)

fig.savefig("figures/out/fig_s3_gene_content.png", dpi=300)
print("bins", cnt, "per-colony median", per_col.median(), "range", per_col.min(), per_col.max())
