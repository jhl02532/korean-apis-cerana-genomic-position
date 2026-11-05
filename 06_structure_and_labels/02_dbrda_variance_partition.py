#!/usr/bin/env python3
# 거리기반 분산분석 (McArdle & Anderson 2001): 조상 구배를 제거한 뒤에도
# 신고 계통 / 지역 / 사육장이 유전체 거리를 설명하는가?
import numpy as np, pandas as pd, json, itertools, os
N = os.environ["KH_ROOT"]
OUT = os.environ["KH_RESULTS"] + "/label_tests"; os.makedirs(OUT, exist_ok=True)
NPERM, SEED = 9999, 20261013
rng = np.random.default_rng(SEED)

meta = pd.read_csv(os.environ["KH_RESULTS"] + "/summary_tables/snp_vs_sv_structure.tsv", sep="\t")
meta = meta.set_index("sample")

def load_dist(path, sep="\t"):
    d = pd.read_csv(path, sep=sep, index_col=0)
    d.columns = [c.strip() for c in d.columns]
    return d

dists = {}
dists["SNP (97,609 sites)"] = load_dist(os.environ["KH_RESULTS"] + "/snpstruct/snp_dist.tsv")
try:
    dists["SV (12,219 loci, Jaccard)"] = load_dist(os.environ["KH_RESULTS"] + "/sv_jaccard_distance.tsv")
except Exception as e:
    print("SV distance not loaded:", e)

def gower(D):
    A = -0.5 * D**2
    n = A.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    return J @ A @ J

def hat(X):
    return X @ np.linalg.pinv(X.T @ X) @ X.T

def dbrda(G, Xfull, Xred, labels_perm_idx=None):
    """sequential: extra variance explained by the terms in Xfull beyond Xred"""
    Hf, Hr = hat(Xfull), hat(Xred)
    n = G.shape[0]
    ss_total = np.trace(G)
    ss_full = np.trace(Hf @ G @ Hf)
    ss_red  = np.trace(Hr @ G @ Hr)
    ss_term = ss_full - ss_red
    df_term = np.linalg.matrix_rank(Xfull) - np.linalg.matrix_rank(Xred)
    df_res  = n - np.linalg.matrix_rank(Xfull)
    ss_res  = ss_total - ss_full
    F = (ss_term / df_term) / (ss_res / df_res) if df_res > 0 and df_term > 0 else np.nan
    return F, ss_term / ss_total, df_term, df_res

def dummies(series):
    cats = sorted(series.unique())
    return np.column_stack([(series.values == c).astype(float) for c in cats[1:]]) if len(cats) > 1 else None

rows = []
for dname, D in dists.items():
    samp = [s for s in D.index if s in meta.index]
    D2 = D.loc[samp, samp].to_numpy(dtype=float)
    m = meta.loc[samp]
    G = gower(D2)
    n = len(samp)
    one = np.ones((n, 1))
    idxv = m["china_index_poly"].to_numpy(dtype=float).reshape(-1, 1)
    idxv = (idxv - idxv.mean()) / idxv.std()
    factors = {"declared strain": m["breed_en"].astype(str),
               "region": m["region"].astype(str),
               "apiary": m["apiary_code"].astype(str)}
    # 조상 지수 단독
    F, R2, dft, dfr = dbrda(G, np.hstack([one, idxv]), one)
    perm = []
    for _ in range(NPERM):
        p = rng.permutation(n)
        perm.append(dbrda(G, np.hstack([one, idxv[p]]), one)[0])
    pv = (np.sum(np.array(perm) >= F) + 1) / (NPERM + 1)
    rows.append(dict(data=dname, term="ancestry index", model="alone",
                     df=dft, R2=round(R2, 4), pseudo_F=round(F, 3), p=round(pv, 4), n=n))
    for fname, fser in factors.items():
        Xd = dummies(fser)
        if Xd is None: continue
        # (a) 단독
        F, R2, dft, dfr = dbrda(G, np.hstack([one, Xd]), one)
        perm = []
        for _ in range(NPERM):
            p = rng.permutation(n)
            perm.append(dbrda(G, np.hstack([one, Xd[p]]), one)[0])
        pv = (np.sum(np.array(perm) >= F) + 1) / (NPERM + 1)
        rows.append(dict(data=dname, term=fname, model="alone", df=dft,
                         R2=round(R2, 4), pseudo_F=round(F, 3), p=round(pv, 4), n=n))
        # (b) 조상 지수를 먼저 넣은 뒤 추가 설명력
        Xred = np.hstack([one, idxv]); Xfull = np.hstack([one, idxv, Xd])
        F, R2, dft, dfr = dbrda(G, Xfull, Xred)
        perm = []
        for _ in range(NPERM):
            p = rng.permutation(n)
            perm.append(dbrda(G, np.hstack([one, idxv, Xd[p]]), Xred)[0])
        pv = (np.sum(np.array(perm) >= F) + 1) / (NPERM + 1)
        rows.append(dict(data=dname, term=fname, model="after ancestry index", df=dft,
                         R2=round(R2, 4), pseudo_F=round(F, 3), p=round(pv, 4), n=n))
res = pd.DataFrame(rows)
res.to_csv(OUT + "/variance_partition.tsv", sep="\t", index=False)
print(res.to_string(index=False))
json.dump(dict(n_perm=NPERM, seed=SEED), open(OUT + "/variance_partition.json", "w"))
print("DONE PERMANOVA")
