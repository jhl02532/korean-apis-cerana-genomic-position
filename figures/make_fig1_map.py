"""Fig. 1. Sampling sites and reference individuals.
Run from the repository root: python figures/make_fig1_map.py"""
import pandas as pd, matplotlib.pyplot as plt, matplotlib.lines as ml
import os; os.makedirs("figures/out", exist_ok=True)
import cartopy.crs as ccrs, cartopy.feature as cf
PC = ccrs.PlateCarree()
T = pd.read_csv("results/sites_19.tsv", sep="\t")
REG = {"Cheorwon": "Cheorwon", "Gapyeong": "Gapyeong", "Icheon": "Icheon", "Yangpyeong": "Yangpyeong", "Wonju": "Wonju",
       "Incheon (INU)": "Incheon (INU)"}
T["reg"] = T.site.map(REG); assert T.reg.notna().all() and len(T) == 19
COL = {"Icheon": "#009E73", "Gapyeong": "#56B4E9", "Cheorwon": "#0072B2", "Yangpyeong": "#CC79A7", "Wonju": "#D55E00",
       "Incheon (INU)": "#7F7F7F"}
LL = {"Cheorwon": (127.31, 38.15), "Gapyeong": (127.51, 37.83), "Yangpyeong": (127.49, 37.49), "Wonju": (127.92, 37.34),
      "Icheon": (127.44, 37.27), "Incheon (INU)": (126.63, 37.37)}
CN = [(105.29,27.30),(106.6,26.6),(108.32,22.82),(111.57,24.40),(113.93,35.30),(102.22,31.90),(102.06,31.48),(114.17,22.32),(110.3,19.6),
      (110.80,19.54),(126.42,41.94),(129.51,42.89),(120.84,40.71),(124.92,42.10),(100.97,22.79),(98.86,25.85),(118.36,35.10),(118.05,36.81)]
JP = [(141.15,39.70),(140.33,38.25),(140.55,35.95),(139.64,35.44),(138.88,37.07),(137.82,35.51),(136.51,34.72),(136.11,34.63),(135.77,35.01),
      (132.87,33.79),(132.87,33.79),(133.05,35.47),(131.43,33.18),(130.51,33.32),(130.51,33.32),(130.19,32.46),(130.76,31.74),(130.44,31.38)]
ASM_CN = [(109.4,19.2),(101.7,32.9),(110.5,31.6),(101.43,36.04)]
def cn_grp(lo, la): return "NE China" if la > 40 and lo > 118 else "S/C China"
def jp_grp(lo, la): return "Japan S" if lo < 131.8 and la < 34 else "Japan N/C"
CC, CJ = "#D55E00", "#E6AB02"
STY = {"S/C China": dict(marker="^", fc=CC, ec="k", lab="China, outside the\nnortheast"),
       "NE China": dict(marker="D", fc="#F4A582", ec="k", lab="Northeast China\n(Jilin, Liaoning)"),
       "Japan N/C": dict(marker="v", fc=CJ, ec="k", lab="Japan, Honshu–Shikoku"),
       "Japan S": dict(marker="v", fc="white", ec="#B8860B", lab="Japan, Kyushu")}
MS = {"^": 1.0, "v": 1.0, "D": 0.7}
NN = {g: 0 for g in STY}
for lo, la in CN: NN[cn_grp(lo, la)] += 1
for lo, la in JP: NN[jp_grp(lo, la)] += 1
assert NN["NE China"] == 4 and NN["Japan S"] == 6 and sum(NN.values()) == 36
FT = 6.5
fig = plt.figure(figsize=(6.7, 2.75), dpi=300)
# ---- A
axA = fig.add_axes([0.065, 0.08, 0.60, 0.89], projection=PC); axA.set_extent((97, 146, 17, 46), crs=PC)
axA.add_feature(cf.LAND.with_scale("50m"), facecolor="#F4F4F4"); axA.add_feature(cf.COASTLINE.with_scale("50m"), lw=0.35, edgecolor="#888888")
axA.add_feature(cf.BORDERS.with_scale("50m"), lw=0.3, edgecolor="#BBBBBB")
for pts, fn in [(CN, cn_grp), (JP, jp_grp)]:
    for lo, la in pts:
        s = STY[fn(lo, la)]; axA.scatter(lo, la, marker=s["marker"], s=22 * MS[s["marker"]], facecolor=s["fc"], edgecolor=s["ec"], lw=0.45, zorder=4, transform=PC)
for lo, la in ASM_CN: axA.scatter(lo, la, marker="*", s=65, color="#FFD700", edgecolor="k", lw=0.4, zorder=6, transform=PC)
for r, (lo, la) in LL.items(): axA.scatter(lo, la, s=10, color="#333333", edgecolor="white", lw=0.25, zorder=5, transform=PC)
axA.add_patch(plt.Rectangle((125.6, 36.6), 3.7, 2.1, fill=False, lw=0.6, ec="k", transform=PC, zorder=7))
axA.text(126.0, 39.0, "B", fontsize=7, fontweight="bold", transform=PC, va="bottom")
gl = axA.gridlines(draw_labels=True, lw=0.2, color="#DDDDDD", xlocs=[100, 110, 120, 130, 140], ylocs=[20, 30, 40])
gl.top_labels = gl.right_labels = False; gl.xlabel_style = gl.ylabel_style = {"size": FT}
# ---- B: Korea
axI = fig.add_axes([0.725, 0.575, 0.265, 0.395], projection=PC); axI.set_extent((125.6, 129.3, 36.6, 38.7), crs=PC)
axI.add_feature(cf.LAND.with_scale("10m"), facecolor="#F4F4F4"); axI.add_feature(cf.COASTLINE.with_scale("10m"), lw=0.35, edgecolor="#888888")
axI.add_feature(cf.BORDERS.with_scale("10m"), lw=0.3, edgecolor="#BBBBBB")
gi = axI.gridlines(draw_labels=True, lw=0.2, color="#DDDDDD", xlocs=[126, 128], ylocs=[37, 38]); gi.top_labels = gi.right_labels = False
gi.xlabel_style = gi.ylabel_style = {"size": FT}
for r, (lo, la) in LL.items():
    n = int((T.reg == r).sum()); axI.scatter(lo, la, s=10 + 9 * n, color=COL[r], edgecolor="white", lw=0.3, zorder=5, transform=PC)
LAB = {"Cheorwon": (127.0, 38.48), "Gapyeong": (129.22, 38.02), "Yangpyeong": (129.22, 37.66), "Wonju": (128.3, 37.12),
       "Icheon": (127.6, 36.75), "Incheon (INU)": (125.68, 36.72)}
HA = {"Gapyeong": "right", "Yangpyeong": "right"}
for r, (lo, la) in LL.items():
    n = int((T.reg == r).sum()); tx, ty = LAB[r]
    axI.annotate(f"{r} ({n})", xy=(lo, la), xytext=(tx, ty), fontsize=FT, color=COL[r], ha=HA.get(r, "left"), xycoords=PC._as_mpl_transform(axI),
                 textcoords=PC._as_mpl_transform(axI), arrowprops=dict(arrowstyle="-", lw=0.4, color=COL[r]), bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.4))
for s in axI.spines.values(): s.set_linewidth(0.6)
# ---- legend (지도 기호만)
H = [ml.Line2D([], [], marker=s["marker"], ls="", mfc=s["fc"], mec=s["ec"], mew=0.6, ms=5 * (0.85 if s["marker"] == "D" else 1),
               label=f"{s['lab']} (n = {NN[g]})") for g, s in STY.items()]
H += [ml.Line2D([], [], marker="*", ls="", color="#FFD700", mec="k", mew=0.4, ms=8, label="Chinese genome assemblies\ndefining the diagnostic sites"),
      ml.Line2D([], [], marker="o", ls="", mfc="#333333", mec="white", ms=4.5, label="Colony sites, this study (n = 19);\ncoloured by region in B")]
fig.legend(handles=H, loc="upper left", bbox_to_anchor=(0.695, 0.51), fontsize=FT, frameon=False, handletextpad=0.4, labelspacing=0.45)
for a, l, dx in [(axA, "A", -0.075), (axI, "B", -0.15)]: a.text(dx, 1.0, l, transform=a.transAxes, fontsize=10, fontweight="bold", va="top")
fig.savefig("figures/out/fig1_sampling_map.png", dpi=300); print("ok", NN)
