"""Fig. 4. Colonies ranked by ancestry index with apiary and declared-strain tiles, and label permutation tests.
Run from the repository root: python figures/make_fig4.py"""
import pandas as pd, numpy as np, matplotlib.pyplot as plt, matplotlib.patches as mp, matplotlib.lines as ml
import os; os.makedirs("figures/out", exist_ok=True)
K = pd.read_csv("results/colony_labels_index.tsv", sep="\t").dropna(subset=["china_index_poly"])
assert len(K) == 19
K = K.sort_values("china_index_poly").reset_index(drop=True); x = K.china_index_poly.values
BC = {"Tojong (native)": "#0072B2", "Halla": "#E69F00", "Cherry × Tojong": "#CC79A7", "Cherry": "#C1440E"}
BL = {"Tojong (native)": r"$\it{tojong}$", "Halla": r"$\it{Hallabeol}$", "Cherry × Tojong": r"$\it{Cherry}$ × $\it{tojong}$", "Cherry": r"$\it{Cherry}$"}
SH = {"Tojong (native)": "T", "Halla": "H", "Cherry × Tojong": "C×T", "Cherry": "C"}
multi = [g for g, n in K.apiary_code.value_counts().items() if n > 1]
AP = dict(zip(sorted(multi), ["#E3E3E3", "#C4C4C4", "#A3A3A3", "#7E7E7E", "#5A5A5A", "#3C3C3C"]))   # 회색 명암; 단독 양봉장은 흰 타일
APTXT = {g: ("white" if i >= 3 else "black") for i, g in enumerate(sorted(multi))}

# 순열 귀무분포 (p = 0.17, 0.69; 본문 17%·69%)
rng = np.random.default_rng(20261014); NA, NT = 20000, 200000
a = K.apiary_code.values; sizes = [int((a == g).sum()) for g in multi]
obs_a = max(np.ptp(x[a == g]) for g in multi); widest = max(multi, key=lambda g: np.ptp(x[a == g]))
null_a = np.empty(NA)
for r in range(NA):
    p = rng.permutation(x); i = 0; m = 0
    for s in sizes: m = max(m, np.ptp(p[i:i + s])); i += s
    null_a[r] = m
t = (K.breed_en == "Tojong (native)").values; obs_t = np.ptp(x[t])
null_t = np.array([np.ptp(rng.choice(x, t.sum(), replace=False)) for _ in range(NT)])
P_a, P_t = (null_a >= obs_a).mean(), (null_t >= obs_t).mean()

fig = plt.figure(figsize=(6.7, 5.4), dpi=300)
G = fig.add_gridspec(2, 1, height_ratios=[3.9, 2.0], hspace=0.34, left=0.115, right=0.98, top=0.875, bottom=0.085)
gs = G[0].subgridspec(3, 1, height_ratios=[3.0, 0.42, 0.42], hspace=0.08)
gb = G[1].subgridspec(1, 2, wspace=0.22)
ax = fig.add_subplot(gs[0]); xs = np.arange(19)
ax.vlines(xs, 0.48, x, color="#CCCCCC", lw=0.8, zorder=1)
ax.scatter(xs, x, s=34, c=[BC[b] for b in K.breed_en], edgecolor="white", lw=0.5, zorder=3)
ax.set_xlim(-0.6, 18.6); ax.set_ylim(0.48, 0.83); ax.set_xticks([])
ax.set_ylabel("Ancestry index\n← northeast Asian    Chinese mainland →", fontsize=7); ax.tick_params(labelsize=6.5)
ax.spines[["top", "right", "bottom"]].set_visible(False)
ax.set_title("A  Colonies ranked by ancestry index", loc="left", fontsize=8, fontweight="bold")
for row, (lab, col) in enumerate([("Apiary", "apiary_code"), ("Declared\nstrain", "breed_en")]):
    tx = fig.add_subplot(gs[1 + row], sharex=ax)
    for i, v in enumerate(K[col]):
        fc = (AP.get(v, "white") if col == "apiary_code" else BC[v])
        tx.add_patch(mp.Rectangle((i - 0.48, 0.04), 0.96, 0.92, fc=fc, ec="#9A9A9A" if fc == "white" else "white", lw=0.5))
        txt = v if col == "apiary_code" else SH[v]
        tx.text(i, 0.5, txt, ha="center", va="center", fontsize=6 if col == "apiary_code" else 6.5,
                color=APTXT.get(v, "black") if col == "apiary_code" else "white", fontweight="bold" if col != "apiary_code" else None)
    tx.set_ylim(0, 1); tx.set_yticks([0.5]); tx.set_yticklabels([lab], fontsize=6.5); tx.tick_params(length=0)
    tx.set_xticks([]); [s.set_visible(False) for s in tx.spines.values()]
tx.set_xlabel("Colonies, ranked from lowest to highest ancestry index", fontsize=7, labelpad=3)

def hist(pos, null, obs, P, L, title, lab, lo):
    h = fig.add_subplot(gb[pos]); w = 0.0125; bins = obs + w * np.arange(-int((obs - lo) / w), int((0.31 - obs) / w) + 2); wt = np.full(len(null), 100 / len(null))
    h.hist(null[null < obs], bins=bins, weights=wt[null < obs], color="#BBBBBB", ec="white", lw=0.3)
    h.hist(null[null >= obs], bins=bins, weights=wt[null >= obs], color="#7A7A7A", ec="white", lw=0.3, label="at least as wide")
    h.axvline(obs, color="black", lw=1.1)
    h.set_xlim(lo, 0.32); h.set_xticks([0.10, 0.15, 0.20, 0.25, 0.30]); h.text(obs - 0.004, h.get_ylim()[1] * 0.97, f"{lab}\np = {P:.2f} (dark bars)", ha="right", va="top", fontsize=6.5, bbox=dict(fc="white", ec="none", pad=1.2))
    h.set_xlabel("Span of ancestry index", fontsize=7); h.set_ylabel("Permutations (%)", fontsize=7) if pos == 0 else None; h.tick_params(labelsize=6.5)
    h.spines[["top", "right"]].set_visible(False)
    h.set_title(f"{L}  {title}", loc="left", fontsize=8, fontweight="bold")
hist(0, null_a, obs_a, P_a, "B", "Widest apiary span vs random labels", f"observed ({widest})", 0.10)
hist(1, null_t, obs_t, P_t, "C", r"$\it{tojong}$ span vs random 14 colonies", r"observed (14 $\it{tojong}$)", 0.10)

H = [ml.Line2D([], [], marker="s", ls="", color=c, ms=6, label=f"{SH[b]}  {BL[b]}") for b, c in BC.items()]
fig.legend(handles=H, loc="upper center", ncol=4, fontsize=7, frameon=False, bbox_to_anchor=(0.55, 0.995), title="Declared strain (point and lower tile colour)", title_fontsize=7)
fig.savefig("figures/out/fig4_labels.png", dpi=300)
print("ok", widest, round(obs_a, 4), round(P_a, 3), round(obs_t, 4), round(P_t, 3))
