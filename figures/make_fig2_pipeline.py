"""Fig. 2. Analysis pipeline.
Run from the repository root: python figures/make_fig2_pipeline.py"""
import matplotlib.pyplot as plt
import os; os.makedirs("figures/out", exist_ok=True)
from matplotlib.patches import FancyBboxPatch
W, H = 6.7, 3.7
fig, ax = plt.subplots(figsize=(W, H), dpi=300); ax.set_xlim(0, 100); ax.set_ylim(0, 64); ax.axis("off"); fig.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
C = {"data": ("#E8F0F8", "#0072B2"), "var": ("#EFF6EE", "#2E8B57"), "ref": ("#FDF0E3", "#D55E00"), "ana": ("#F6F6F6", "#666666")}
def box(x, y, w, h, title, tool, kind):
    fc, ec = C[kind]
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.3,rounding_size=1.4", fc=fc, ec=ec, lw=1.0))
    if not tool: ax.text(x + w / 2, y + h / 2, title, ha="center", va="center", fontsize=7.5, fontweight="bold", color="#1a1a1a"); return
    ax.text(x + w / 2, y + h * 0.70, title, ha="center", va="center", fontsize=7.5, fontweight="bold", color="#1a1a1a")
    ax.text(x + w / 2, y + h * 0.31, tool, ha="center", va="center", fontsize=6.6, color="#444444", style="italic", linespacing=1.25)
def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1), arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#888888", shrinkA=0, shrinkB=0, mutation_scale=7))
def stage(y, label): ax.text(0.2, y, label, rotation=90, ha="center", va="center", fontsize=7, color="#888888", fontweight="bold")
# 1. data generation
Y1, H1 = 50, 11; stage(Y1 + H1 / 2, "Data")
box(3, Y1, 29, H1, "19 colonies", "one worker each", "data")
box(36, Y1, 29, H1, "PacBio HiFi sequencing", "Revio, SMRT Link v13.3", "data")
box(69, Y1, 29, H1, "Read alignment", "AcerK_1.0; minimap2 v2.26,\nsamtools v1.19", "data")
arrow(32.6, Y1 + H1 / 2, 35.4, Y1 + H1 / 2); arrow(65.6, Y1 + H1 / 2, 68.4, Y1 + H1 / 2)
# 2. variants
Y2, H2 = 27, 12; stage(Y2 + H2 / 2, "Variants")
V = {"SV": (3, 21), "SNP": (27.5, 23.5), "DS": (54.5, 20.5), "GP": (78.5, 19.5)}
for k, (t, tool, kind) in {"SV": ("Structural variants", "Sniffles2 v2.8.1", "var"),
                           "SNP": ("SNPs", "bcftools v1.24;\n36 reference individuals", "var"),
                           "DS": ("Diagnostic sites", "4 Chinese assemblies;\nminimap2 v2.26", "ref"),
                           "GP": ("Gene presence", "mosdepth v0.3.14", "var")}.items():
    x, w = V[k]; box(x, Y2, w, H2, t, tool, kind); arrow(83.5, Y1 - 0.7, x + w / 2, Y2 + H2 + 0.7)
cx = {k: x + w / 2 for k, (x, w) in V.items()}
# 3. analyses (위치는 위 변이 상자 바로 아래)
Y3, H3 = 4, 9; stage(Y3 + H3 / 2, "Analyses")
A = {"PS": (V["SV"][0], V["SV"][1], "Population\nstructure"), "DL": (V["SNP"][0], V["SNP"][1], "Declared\nlabels"),
     "AI": (V["DS"][0], V["DS"][1], "Ancestry\nindex"), "GC": (V["GP"][0], V["GP"][1], "Gene\ncontent")}
for x, w, t in A.values(): box(x, Y3, w, H3, t, "", "ana")
top3, bot2 = Y3 + H3 + 0.7, Y2 - 0.7
arrow(cx["SV"], bot2, cx["SV"], top3)                      # SV → population structure
arrow(cx["SNP"] - 3, bot2, cx["SV"] + 4, top3)             # SNP → population structure
arrow(cx["SV"] + 4, bot2, cx["SNP"] - 3, top3)             # SV distance → declared labels
arrow(cx["SNP"], bot2, cx["SNP"], top3)                    # SNP distance → declared labels
arrow(cx["DS"], bot2, cx["DS"], top3)                      # diagnostic sites → ancestry index
arrow(V["DS"][0] - 0.7, Y3 + H3 / 2, V["SNP"][0] + V["SNP"][1] + 0.7, Y3 + H3 / 2)   # ancestry index → declared labels
arrow(cx["GP"], bot2, cx["GP"], top3)                      # gene presence → gene content
ax.text(50, Y3 - 3.0, "statistics in Python 3 (NumPy, SciPy, pandas)", ha="center", va="center", fontsize=6.6, color="#444444", style="italic")
fig.savefig("figures/out/fig2_pipeline.png", dpi=300, bbox_inches="tight", pad_inches=0.04); print("ok")
