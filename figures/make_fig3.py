"""Fig. 3. Joint PCA with reference individuals (A) and structure among colonies from SVs (B) and SNPs (C).
Run from the repository root: python figures/make_fig3.py"""
import json, pandas as pd, numpy as np, matplotlib as mpl, matplotlib.pyplot as plt, matplotlib.lines as ml
import os; os.makedirs("figures/out", exist_ok=True)
from scipy.stats import spearmanr
T = pd.read_csv("results/sites_19.tsv", sep="\t"); assert len(T) == 19
NORM = mpl.colors.Normalize(vmin=0.49, vmax=0.81); CMAP = mpl.cm.viridis
FS, FT = 9.0, 8.5
CC, CJ = "#D55E00", "#E6AB02"
STY = {"S/C China": dict(marker="^", fc=CC, ec="k", lab="China, outside the\nnortheast"),
       "NE China": dict(marker="D", fc="#F4A582", ec="k", lab="Northeast China\n(Jilin, Liaoning)"),
       "Japan N/C": dict(marker="v", fc=CJ, ec="k", lab="Japan, Honshu–Shikoku"),
       "Japan S": dict(marker="v", fc="white", ec="#B8860B", lab="Japan, Kyushu")}
MS = {"^": 1.0, "v": 1.0, "D": 0.7}
fig = plt.figure(figsize=(6.7, 7.0), dpi=300)
# ---- A: 결합 PCA (유사 반수체)
axA = fig.add_axes([0.105, 0.605, 0.50, 0.355])
PCA = pd.read_csv("results/joint_pca/pca_pseudohap_oriented.tsv", sep="\t"); VE = json.load(open("results/joint_pca/summary.json"))["pseudohap"]["var_pct"]
RM = {l.split("\t")[0]: l.split("\t")[1:] for l in open("results/ref_meta.txt").read().splitlines() if "\t" in l}
def rgrp(r):
    if r.group == "China": return "NE China" if RM[r.sample][0] in ("China:Jilin", "China:Liaoning") else "S/C China"
    if r.group == "Japan": return "Japan S" if RM[r.sample][1] == "Southern" else "Japan N/C"
    return "colony"
PCA["rg"] = [rgrp(r) for r in PCA.itertuples()]
NN = {}
for g, s in STY.items():
    s_ = PCA[PCA.rg == g]; NN[g] = len(s_)
    axA.scatter(s_.PC1, s_.PC2, marker=s["marker"], s=34 * MS[s["marker"]], facecolor=s["fc"], edgecolor=s["ec"], lw=0.5, zorder=3)
assert NN["NE China"] == 4 and NN["Japan S"] == 6 and sum(NN.values()) == 36
cs = PCA[PCA.rg == "colony"].merge(T[["sample", "china_index_poly"]], on="sample"); assert len(cs) == 19
assert NORM.vmin <= cs.china_index_poly.min() and cs.china_index_poly.max() <= NORM.vmax
axA.scatter(cs.PC1, cs.PC2, c=cs.china_index_poly, cmap=CMAP, norm=NORM, s=32, edgecolor="#333333", lw=0.4, zorder=4)
axA.set_xlabel(f"PC1 ({VE[0]:.1f}%)\n← Japanese          Chinese mainland →", fontsize=FS); axA.set_ylabel(f"PC2 ({VE[1]:.1f}%)", fontsize=FS)
axA.spines[["top", "right"]].set_visible(False); axA.tick_params(labelsize=FT)
H = [ml.Line2D([], [], marker=s["marker"], ls="", mfc=s["fc"], mec=s["ec"], mew=0.6, ms=6 * (0.85 if s["marker"] == "D" else 1),
               label=f"{s['lab']} (n = {NN[g]})") for g, s in STY.items()]
H += [ml.Line2D([], [], marker="o", ls="", mfc=CMAP(0.6), mec="#333333", mew=0.4, ms=5.5,
                label="Colonies, this study (n = 19);\ncoloured by ancestry index\n(A–C; scale at lower right)")]
fig.legend(handles=H, loc="upper left", bbox_to_anchor=(0.635, 0.965), fontsize=FT, frameon=False, handletextpad=0.4, labelspacing=0.7)
# ---- B, C: 군락 구조
S = pd.read_csv("results/snp_vs_sv_structure.tsv", sep="\t"); assert len(S) == 19
for c in ("PC1", "snpPC1"):
    if spearmanr(S.china_index_poly, S[c]).statistic < 0: S[c] = -S[c]
assert NORM.vmin <= S.china_index_poly.min() and S.china_index_poly.max() <= NORM.vmax
NAMED = {"CWN2_AC1", "WJ_AC_3"}
AX = []
for k, (x, y, vx, vy, title, left) in enumerate([("PC1", "PC2", 10.72, 7.48, "SVs (12,219 loci; Jaccard PCoA)", 0.105),
                                                 ("snpPC1", "snpPC2", 11.49, 8.00, "SNPs (97,609 thinned sites; PCA)", 0.535)]):
    ax = fig.add_axes([left, 0.07, 0.31, 0.36]); AX.append(ax)
    ax.scatter(S[x], S[y], c=S.china_index_poly, cmap=CMAP, norm=NORM, s=44, edgecolor="#333333", lw=0.5, zorder=3)
    r = spearmanr(S.china_index_poly, S[x])
    ax.text(0.03, 0.98, f"axis 1 vs index: ρ = {abs(r.statistic):.3f}", transform=ax.transAxes, ha="left", va="top", fontsize=FT)
    for rr in S[S["sample"].isin(NAMED)].itertuples():
        px, py = getattr(rr, x), getattr(rr, y)
        ax.annotate(rr.sample, (px, py), xytext=((2, 6) if rr.sample == "WJ_AC_3" else (-5, 0)), textcoords="offset points", ha="right",
                    va="bottom" if rr.sample == "WJ_AC_3" else "center", fontsize=FT, color="#333333")
    ax.margins(x=0.08, y=0.10); y0, y1 = ax.get_ylim(); ax.set_ylim(y0, y1 + 0.10 * (y1 - y0))   # ρ 문구와 CWN2_AC1 이름 간격
    ax.set_xlabel(f"Axis 1 ({vx:.2f}%)", fontsize=FS); ax.set_ylabel(f"Axis 2 ({vy:.2f}%)", fontsize=FS); ax.set_title(title, fontsize=FS, loc="left")
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=FT)
cb = fig.colorbar(mpl.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=fig.add_axes([0.875, 0.07, 0.016, 0.36]))
cb.set_label("Ancestry index\n← northeast Asian    Chinese mainland →", fontsize=FS); cb.ax.tick_params(labelsize=FT)
for a, l in [(axA, "A"), (AX[0], "B"), (AX[1], "C")]:
    a.text(-0.17 if a is not axA else -0.115, 1.03, l, transform=a.transAxes, fontsize=12, fontweight="bold", va="bottom")
fig.savefig("figures/out/fig3_position_structure.png", dpi=300)
print("ok", NN, "VE", [round(v, 2) for v in VE[:2]], "SV rho", spearmanr(S.china_index_poly, S.PC1).statistic, "SNP rho", spearmanr(S.china_index_poly, S.snpPC1).statistic)
