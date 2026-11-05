#!/usr/bin/env python
# 읽기기반 SV 유전형 행렬 분석 — 한국 A. cerana 19군락, AcerK_1.0 기준
import pandas as pd, numpy as np, json, sys, os
THR=float(sys.argv[1]); OUT=sys.argv[2]
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, fcluster

BASE=os.environ["KH_RESULTS"]
samples=[s for s in open(f"{BASE}/gt_header.txt").read().strip().split("\t") if s]
cols=["chrom","pos","svtype","svlen","support","filt"]
gt=[f"{s}:GT" for s in samples]; dr=[f"{s}:DR" for s in samples]; dv=[f"{s}:DV" for s in samples]
names=cols+[c for trio in zip(gt,dr,dv) for c in trio]
df=pd.read_csv(f"{BASE}/gt_matrix2.tsv", sep="\t", header=None, names=names, low_memory=False)
assert len(samples)==19 and df.shape[1]==63, (len(samples), df.shape)

def silhouette(D, lab):
    s=[]
    for i in range(len(lab)):
        same=(lab==lab[i]); same[i]=False
        if same.sum()==0: s.append(0.0); continue
        a=D[i,same].mean()
        b=min(D[i,lab==k].mean() for k in np.unique(lab) if k!=lab[i])
        s.append((b-a)/max(a,b))
    return float(np.mean(s))

n0=len(df)
df=df[df.chrom.str.startswith("NC_0838")].copy()           # 염색체 배치만
df["svlen"]=pd.to_numeric(df.svlen, errors="coerce")
df=df[df.svlen.abs()>=50].copy()                            # 50 bp 이상
n_len=len(df)

DR=df[dr].apply(pd.to_numeric, errors="coerce").values
DV=df[dv].apply(pd.to_numeric, errors="coerce").values
DP=DR+DV

cov=pd.read_csv(f"{BASE}/qc/all_coverage.tsv", sep="\t", header=None, names=["sample","depth","cov"])
cov=cov.drop_duplicates("sample").set_index("sample").reindex(samples)
assert cov.depth.notna().all(), cov[cov.depth.isna()].index.tolist()

# 깊이 대비 상대 문턱 — 깊은 샘플에서 SV 가 더 나오는 인공물 차단
thr=np.maximum(5.0, 0.15*cov.depth.values)[None,:]
callable_m = DP >= thr
with np.errstate(invalid="ignore", divide="ignore"):
    vaf=np.where(DP>0, DV/DP, np.nan)
present=(vaf>=0.25)&callable_m

keep = callable_m.mean(axis=1) >= THR                      # THR(인자) 비율 군락에서 호출가능
P=present[keep]; C=callable_m[keep]; meta=df.loc[keep, cols].reset_index(drop=True)
n_keep=int(P.shape[0])

per_sample=pd.DataFrame({"sample":samples,"depth":cov.depth.values,"coverage_pct":cov["cov"].values,
                         "n_sv":P.sum(axis=0),"n_callable":C.sum(axis=0)})
per_sample.to_csv(f"{OUT}/sv_per_sample_qc.tsv", sep="\t", index=False)

af=P.sum(axis=1)/C.sum(axis=1)
info=(af>=0.10)&(af<=0.90)                                  # 패널 내 다형만
X=P[info].astype(float).T
n_info=int(info.sum())
assert n_info>500, n_info

D=squareform(pdist(X, metric="jaccard"))
pd.DataFrame(D, index=samples, columns=samples).to_csv(f"{OUT}/sv_jaccard_distance.tsv", sep="\t")

D2=D**2; n=len(samples); J=np.eye(n)-np.ones((n,n))/n
B=-0.5*J@D2@J; w,v=np.linalg.eigh(B); o=np.argsort(w)[::-1]; w,v=w[o],v[:,o]
pos=w>1e-9; coords=v[:,pos]*np.sqrt(w[pos]); varexp=w[pos]/w[pos].sum()

Z=linkage(pdist(X, metric="jaccard"), method="average")
res={}
for k in (2,3,4):
    lab=fcluster(Z,k,criterion="maxclust")
    res[str(k)]={"silhouette":round(silhouette(D,lab),4),"sizes":np.bincount(lab)[1:].tolist()}

pc1=coords[:,0]; srt=np.sort(pc1); gaps=np.diff(srt)
out=pd.DataFrame({"sample":samples,"PC1":coords[:,0].round(4),"PC2":coords[:,1].round(4),
                  "cluster_k2":fcluster(Z,2,criterion="maxclust"),
                  "n_sv":P.sum(axis=0),"depth":cov.depth.values})
out.to_csv(f"{OUT}/sv_pcoa_clusters.tsv", sep="\t", index=False)
meta.to_csv(f"{OUT}/sv_sites_kept.tsv", sep="\t", index=False)

# 깊이 교란 점검: SV 수가 깊이와 상관되면 문턱이 부족하다는 뜻
from scipy.stats import spearmanr, pearsonr
rho_d,p_d=spearmanr(per_sample.depth, per_sample.n_sv)
rho_p,p_p=spearmanr(per_sample.depth, pc1)

summary={"n_sv_raw":n0,"n_after_chrom_len":n_len,"n_callable90":n_keep,"n_informative":n_info,
 "varexp_pc1":round(float(varexp[0]),4),"varexp_pc2":round(float(varexp[1]),4),
 "cluster":res,"max_pc1_gap":round(float(gaps.max()),4),"pc1_range":round(float(np.ptp(pc1)),4),
 "max_gap_frac_of_range":round(float(gaps.max()/np.ptp(pc1)),4),
 "depth_vs_nsv_rho":round(float(rho_d),3),"depth_vs_nsv_p":round(float(p_d),4),
 "depth_vs_pc1_rho":round(float(rho_p),3),"depth_vs_pc1_p":round(float(p_p),4),
 "svtype_kept":{k:int(v) for k,v in meta.svtype.value_counts().items()}}
json.dump(summary, open(f"{OUT}/sv_summary.json","w"), indent=1)
print(json.dumps(summary, indent=1))
print(out.sort_values("PC1").to_string(index=False))

K=pd.read_csv(os.environ["KH_RESULTS"] + "/two_axes_colonies.tsv",sep="\t").set_index("colony").KC
r=spearmanr(out.PC1,out["sample"].map(K)); print("RHO_INDEX",abs(round(float(r.statistic),4)),r.pvalue,"n_info",n_info,"var",summary["varexp_pc1"],"sil",{k:v["silhouette"] for k,v in res.items()})
json.dump({"thr":THR,"n_informative":n_info,"rho_index_abs":abs(float(r.statistic)),"p":float(r.pvalue),"varexp_pc1":summary["varexp_pc1"],"silhouette":{k:v["silhouette"] for k,v in res.items()}},open(f"{OUT}/sensitivity.json","w"),indent=1)
