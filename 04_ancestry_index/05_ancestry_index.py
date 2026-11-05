#!/usr/bin/env python
import pandas as pd, numpy as np, json, os
D=os.environ["KH_RESULTS"] + "/china_ref"
samples=[s for s in open(f"{D}/geno_samples.txt").read().strip().split("\t") if s]
rows=[]
for ln in open(f"{D}/ad_matrix.tsv"):
    p=ln.rstrip("\n").split("\t")
    if len(p)!=4+len(samples): continue
    ad=[]
    for cell in p[4:]:
        try:
            a,b = cell.split(",")[:2]; ad.append((int(a),int(b)))
        except Exception: ad.append((0,0))
    rows.append((p[0],int(p[1]),ad))
REF=np.array([[x[0] for x in r[2]] for r in rows], float)
ALT=np.array([[x[1] for x in r[2]] for r in rows], float)
DP=REF+ALT
assert DP.shape[1]==19, DP.shape

MIN_DP=8
ok = DP>=MIN_DP
with np.errstate(invalid="ignore", divide="ignore"):
    vaf = np.where(ok, ALT/np.where(DP>0,DP,np.nan), np.nan)   # ALT = 중국형 대립유전자

site_ok = np.nanmean(ok,axis=1) >= 0.90
site_mean = np.nanmean(np.where(ok,vaf,np.nan), axis=1)
# 패널 내 다형 사이트만 — AcerK 사유 대립유전자(전원 1.0) 와 비정보 사이트 제거
poly = site_ok & (site_mean>0.05) & (site_mean<0.95)
sel = vaf[poly]
idx_all = int(site_ok.sum()); idx_poly = int(poly.sum())

frac_all  = np.nanmean(np.where(ok,vaf,np.nan)[site_ok], axis=0)
frac_poly = np.nanmean(sel, axis=0)
out=pd.DataFrame({"sample":samples,
                  "china_index_all":np.round(frac_all,5),
                  "china_index_poly":np.round(frac_poly,5),
                  "n_sites_poly":np.sum(~np.isnan(sel),axis=0)})
out.to_csv(f"{D}/china_free_ancestry.tsv", sep="\t", index=False)
print(json.dumps({"n_sites_total":len(rows),"n_sites_callable90":idx_all,
                  "n_sites_polymorphic":idx_poly,"min_dp":MIN_DP}, indent=1))
print(out.sort_values("china_index_poly").to_string(index=False))
