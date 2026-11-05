#!/usr/bin/env python3
# 신고 계통·양봉장 표지와 조상 지수: 범위, Mann–Whitney, 양봉장별 범위, 최소검출효과(MDE)
# 입력: $KH_RESULTS/two_axes_colonies.tsv (08_index_bootstrap_and_summary.py), metadata/colonies.tsv
# 표지 순열 검정(tojong 범위, 최대 양봉장 범위)은 figures/make_fig4.py 에 있다.
import json, os, numpy as np, pandas as pd, scipy.stats as st
R = os.environ["KH_RESULTS"]; HERE = os.path.dirname(os.path.abspath(__file__))
K = pd.read_csv(f"{R}/two_axes_colonies.tsv", sep="\t").set_index("colony").KC
M = pd.read_csv(os.path.join(HERE, "../metadata/colonies.tsv"), sep="\t").set_index("sample")
M["index"] = K.reindex(M.index); assert M["index"].notna().all() and len(M) == 19
lo, hi = M["index"].min(), M["index"].max(); span = hi - lo
t = M[M.declared_strain == "tojong (native)"]["index"]; o = M[M.declared_strain != "tojong (native)"]["index"]
h = M[M.declared_strain == "Hallabeol"]["index"]
res = {"range": [lo, hi, span],
       "tojong": {"n": len(t), "min": t.min(), "max": t.max(), "span_pct": round(100 * (t.max() - t.min()) / span, 1), "median": t.median()},
       "halla": [h.min(), h.max()], "others": o.to_dict(),
       "MW_tojong_vs_other": [st.mannwhitneyu(t, o, alternative="two-sided", method="exact").statistic, st.mannwhitneyu(t, o, alternative="two-sided", method="exact").pvalue, o.median()],
       "MW_tojong_vs_halla": [st.mannwhitneyu(t, h, alternative="two-sided", method="exact").statistic, st.mannwhitneyu(t, h, alternative="two-sided", method="exact").pvalue],
       "apiary_spans": {a: round(100 * (g.max() - g.min()) / span, 1) for a, g in M.groupby("apiary")["index"] if len(g) > 1}}
SD = float(M[M.apiary != "INU"]["index"].std(ddof=1)); res["SD_farm"] = SD


def mde(sizes, shift_idx, nsim=4000, seed=20261014):
    """정규(SD_farm) 모의, 한 집단 이동, 양측 Mann–Whitney α=0.05, 검정력 80% 가 되는 최소 이동 (이분 탐색)."""
    rng = np.random.default_rng(seed); a, b = 0.0, 0.6
    for _ in range(14):
        d = (a + b) / 2; hit = 0
        for _ in range(nsim):
            g = [rng.normal(0, SD, n) for n in sizes]; g[shift_idx] = g[shift_idx] + d
            hit += st.mannwhitneyu(*g, alternative="two-sided").pvalue < 0.05
        a, b = (a, d) if hit / nsim >= 0.8 else (d, b)
    return round(b, 3)


res["MDE_tojong_vs_halla_14v3"] = mde([len(t), len(h)], 1)
os.makedirs(f"{R}/summary_tables", exist_ok=True)
json.dump(res, open(f"{R}/summary_tables/label_tests.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
