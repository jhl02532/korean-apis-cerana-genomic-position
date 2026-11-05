#!/usr/bin/env python3
# 신고 계통 효과를 교란과 분리해 다시 검정한다.
#  A. 토종벌(14) vs 한라벌(3) — 양쪽 모두 반복이 있는 유일한 대비
#  C. 지역을 먼저 넣은 뒤 계통이 남는가 (한라벌 3군 중 2군이 같은 양평 사육장)
import numpy as np, pandas as pd, os, json
N = os.environ["KH_ROOT"]
OUT = os.environ["KH_RESULTS"] + "/label_tests"
NPERM, SEED = 9999, 20261013
rng = np.random.default_rng(SEED)
meta = pd.read_csv(os.environ["KH_RESULTS"] + "/summary_tables/snp_vs_sv_structure.tsv", sep="\t").set_index("sample")

def load(p):
    d = pd.read_csv(p, sep="\t", index_col=0); d.columns = [c.strip() for c in d.columns]; return d
DS = {"SNP": load(os.environ["KH_RESULTS"] + "/snpstruct/snp_dist.tsv"),
      "SV":  load(os.environ["KH_RESULTS"] + "/sv_jaccard_distance.tsv")}

def gower(D):
    A = -0.5 * D**2; n = A.shape[0]; J = np.eye(n) - np.ones((n, n)) / n
    return J @ A @ J
def hat(X): return X @ np.linalg.pinv(X.T @ X) @ X.T
def term_test(G, Xred, Zcols):
    Xf = np.hstack([Xred, Zcols])
    Hf, Hr = hat(Xf), hat(Xred)
    tot = np.trace(G); ssf = np.trace(Hf @ G @ Hf); ssr = np.trace(Hr @ G @ Hr)
    dft = np.linalg.matrix_rank(Xf) - np.linalg.matrix_rank(Xred)
    dfr = G.shape[0] - np.linalg.matrix_rank(Xf)
    F = ((ssf - ssr) / dft) / ((tot - ssf) / dfr) if dft > 0 and dfr > 0 else np.nan
    return F, (ssf - ssr) / tot, dft
def dumm(s):
    cats = sorted(pd.unique(s))
    return np.column_stack([(np.asarray(s) == c).astype(float) for c in cats[1:]]) if len(cats) > 1 else None

rows = []
def run(tag, dname, D, samp, Zser, cov_names):
    m = meta.loc[samp]
    G = gower(D.loc[samp, samp].to_numpy(float)); n = len(samp)
    one = np.ones((n, 1))
    cov = [one]
    if "index" in cov_names:
        v = m["china_index_poly"].to_numpy(float).reshape(-1, 1)
        cov.append((v - v.mean()) / v.std())
    if "region" in cov_names:
        rd = dumm(m["region"].astype(str))
        if rd is not None: cov.append(rd)
    Xred = np.hstack(cov)
    Z = dumm(Zser.loc[samp].astype(str))
    if Z is None: return
    F, R2, dft = term_test(G, Xred, Z)
    null = []
    for _ in range(NPERM):
        p = rng.permutation(n)
        null.append(term_test(G, Xred, Z[p])[0])
    pv = (np.sum(np.array(null) >= F) + 1) / (NPERM + 1)
    rows.append(dict(test=tag, data=dname, n=n, df=dft, covariates="+".join(cov_names) or "none",
                     R2=round(R2, 4), pseudo_F=round(F, 3), p=round(pv, 4)))

for dname, D in DS.items():
    samp_all = [s for s in D.index if s in meta.index]
    br = meta["breed_en"].astype(str)
    # A. 토종벌 vs 한라벌
    sA = [s for s in samp_all if br[s] in ("Tojong (native)", "Halla")]
    for cov in ([], ["index"], ["region"], ["index", "region"]):
        run("A. Tojong(14) vs Halla(3)", dname, D, sA, br, cov)
    # C. 전체 4수준 (단독군 2개 포함) — 지역 통제
    for cov in (["index"], ["region"], ["index", "region"]):
        run("C. declared strain, 4 levels", dname, D, samp_all, br, cov)
    # D. 단독군을 하나로 묶은 3수준
    mg = br.copy()
    mg[:] = np.where(br.isin(["Tojong (native)", "Halla"]), br, "other(singletons)")
    for cov in (["index"], ["index", "region"]):
        run("D. 3 levels (singletons merged)", dname, D, samp_all, mg, cov)

res = pd.DataFrame(rows)
res.to_csv(OUT + "/strain_effect_tests.tsv", sep="\t", index=False)
print(res.to_string(index=False))

# E. 한라벌 3군이 서로 가까운가 — 조상지수 맞춤 3개 조합 귀무
ex = []
for dname, D in DS.items():
    samp = [s for s in D.index if s in meta.index]
    m = meta.loc[samp]; Dm = D.loc[samp, samp].to_numpy(float)
    halla = [i for i, s in enumerate(samp) if m.loc[s, "breed_en"] == "Halla"]
    obs = np.mean([Dm[i, j] for i in halla for j in halla if i < j])
    idxv = m["china_index_poly"].to_numpy(float)
    tgt = np.sort(idxv[halla])
    null, tried = [], 0
    while len(null) < 20000 and tried < 400000:
        tried += 1
        pick = rng.choice(len(samp), 3, replace=False)
        if np.abs(np.sort(idxv[pick]) - tgt).mean() > 0.05: continue   # 조상지수 맞춤
        null.append(np.mean([Dm[i, j] for i in pick for j in pick if i < j]))
    null = np.array(null)
    ex.append(dict(data=dname, observed_mean_within_Halla=round(float(obs), 5),
                   matched_null_mean=round(float(null.mean()), 5), n_null=len(null),
                   ratio=round(float(obs / null.mean()), 3),
                   p_emp=round(float((np.sum(null <= obs) + 1) / (len(null) + 1)), 4)))
pd.DataFrame(ex).to_csv(OUT + "/halla_cohesion.tsv", sep="\t", index=False)
print(pd.DataFrame(ex).to_string(index=False))
print("DONE STRAIN")
