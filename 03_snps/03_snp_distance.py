#!/usr/bin/env python
# SNP 거리행렬 (표준화 유전형 전체 공간) — SV Jaccard 거리와 같은 층위에서 실루엣 비교용
import numpy as np, pandas as pd, os
O=os.environ["KH_RESULTS"] + "/snpstruct"
hdr=open(f"{O}/thin2000.raw").readline().split()
iids=[]; rows=[]
for ln in open(f"{O}/thin2000.raw").readlines()[1:]:
    p=ln.split(); iids.append(p[1])
    rows.append([np.nan if x=="NA" else float(x) for x in p[6:]])
G=np.array(rows); mu=np.nanmean(G,axis=0)
G=np.where(np.isnan(G), mu, G); p=mu/2.0
k=(p>0.01)&(p<0.99); G=G[:,k]; p=p[k]
Z=(G-2*p)/np.sqrt(2*p*(1-p))
n=len(iids)
D=np.zeros((n,n))
for i in range(n):
    for j in range(i+1,n):
        d=np.sqrt(np.mean((Z[i]-Z[j])**2))   # SNP 당 정규화 유클리드
        D[i,j]=D[j,i]=d
pd.DataFrame(D, index=iids, columns=iids).to_csv(f"{O}/snp_dist.tsv", sep="\t")
print(f"거리행렬 {n}×{n} · SNP {Z.shape[1]:,}")
print(f"거리 범위 {D[np.triu_indices(n,1)].min():.4f}–{D[np.triu_indices(n,1)].max():.4f}")
