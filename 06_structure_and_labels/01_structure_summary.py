#!/usr/bin/env python3
# 19군락 구조 요약: 구조표, 실루엣, 상관, k=2 구성, 최대 간격 + 귀무, 실루엣 모의 보정
import numpy as np, pandas as pd, json, os
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
S=os.environ["KH_RESULTS"]; R=os.environ["KH_RESULTS"]
meta=pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../metadata/colonies.tsv"),sep="\t").set_index("sample")
BR={"tojong (native)":"Tojong (native)","Hallabeol":"Halla","Cherry":"Cherry","Cherry × tojong":"Cherry × Tojong"}   # 하위 스크립트의 내부 표지
old=pd.DataFrame({"region":meta.region,"breed_en":meta.declared_strain.map(BR),"apiary_code":meta.apiary}); assert old.notna().all().all()
os.makedirs(f"{R}/summary_tables",exist_ok=True)
snp=pd.read_csv(f"{R}/snpstruct/snp_pca.tsv",sep="\t").set_index("sample")
sv=pd.read_csv(f"{R}/sv_pcoa_clusters.tsv",sep="\t").set_index("sample")
ix=pd.read_csv(f"{R}/china_ref/china_free_ancestry.tsv",sep="\t").set_index("sample")
T=snp.join(sv[["PC1","PC2"]]).join(ix[["china_index_poly"]]).join(old); assert len(T)==19 and T.notna().all().all()
T.reset_index().to_csv(f"{R}/summary_tables/snp_vs_sv_structure.tsv",sep="\t",index=False)
idx=T.china_index_poly
def sil(D,lab):
    s=[]
    for i in range(len(lab)):
        same=(lab==lab[i]); same[i]=False
        if same.sum()==0: s.append(0.0); continue
        a=D[i,same].mean(); b=min(D[i,lab==k].mean() for k in np.unique(lab) if k!=lab[i]); s.append((b-a)/max(a,b))
    return float(np.mean(s))
out={}
Ds=pd.read_csv(f"{R}/snpstruct/snp_dist.tsv",sep="\t",index_col=0).loc[T.index,T.index]
Dv=pd.read_csv(f"{R}/sv_jaccard_distance.tsv",sep="\t",index_col=0).loc[T.index,T.index]
for nm,D in [("SNP",Ds.values),("SV",Dv.values)]:
    Z=linkage(squareform(D,checks=False),"average"); r={}
    for k in (2,3,4):
        lab=fcluster(Z,k,"maxclust"); r[k]={"sil":round(sil(D,lab),4),"sizes":sorted(np.bincount(lab)[1:].tolist())}
        if k==2: r["k2_minor"]=[s for s,l in zip(T.index,lab) if (lab==l).sum()==min(np.bincount(lab)[1:])]
    out[nm]=r
c=lambda a,b: (round(float(spearmanr(a,b).statistic),4), float(spearmanr(a,b).pvalue))
out["rho_SVpc1_index"]=c(T.PC1,idx); out["rho_SNPpc1_index"]=c(T.snpPC1,idx); out["rho_SNPpc1_SVpc1"]=c(T.snpPC1,T.PC1)
iu=np.triu_indices(19,1); out["rho_dist"]=c(Ds.values[iu],Dv.values[iu]); out["n_pairs"]=len(iu[0])
# 부호 통일한 PC1 에서 최대 간격
rng=np.random.default_rng(20261015)
def gaps(x):
    s=np.sort(x); g=np.diff(s); r=s[-1]-s[0]; o=np.argsort(g)[::-1]
    return g.max()/r, g[o[1]]/r, np.median(g)/r, int(np.argmax(g))+1
null_u=[];null_n=[]
for _ in range(100000):
    u=rng.random(19); n=rng.normal(size=19)
    null_u.append(np.diff(np.sort(u)).max()/np.ptp(u)); null_n.append(np.diff(np.sort(n)).max()/np.ptp(n))
null_u=np.array(null_u); null_n=np.array(null_n)
for nm,x in [("SV",T.PC1.values),("SNP",T.snpPC1.values)]:
    if spearmanr(x,idx).statistic<0: x=-x
    gmax,g2,gmed,pos=gaps(x); order=T.index[np.argsort(x)]
    out[f"gap_{nm}"]={"max":round(gmax,4),"second":round(g2,4),"median":round(gmed,4),"after_rank":pos,"above":list(order[pos:]),
                      "P_uniform":round(float((null_u>=gmax).mean()),3),"P_normal":round(float((null_n>=gmax).mean()),3)}
# 실루엣 모의 보정 (Balding–Nichols, 관측 SNP k=2 분할, 마커 수 = 솎은 SNP 수)
nmk=int(json.load(open(f"{R}/snpstruct/snp_pca_meta.json"))["n_snp"]); sz=out["SNP"][2]["sizes"]
cal=[]
for F in [0,0.01,0.02,0.05,0.10,0.15,0.25]:
    vals=[]
    for rep in range(12):
        p=rng.uniform(0.05,0.95,nmk); G=[]
        for n in sz:
            q=rng.beta(p*(1-F)/F,(1-p)*(1-F)/F) if F>0 else p
            G.append(rng.binomial(2,q,(n,nmk)))
        G=np.vstack(G).astype(float); mu=G.mean(0); pp=mu/2; k=(pp>0.01)&(pp<0.99)
        Zz=(G[:,k]-2*pp[k])/np.sqrt(2*pp[k]*(1-pp[k]))
        D=np.sqrt(((Zz[:,None,:]-Zz[None,:,:])**2).mean(2))
        lab=fcluster(linkage(squareform(D,checks=False),"average"),2,"maxclust"); vals.append(sil(D,lab))
    cal.append({"F_ST":F,"sil_mean":float(np.mean(vals)),"sil_p05":float(np.percentile(vals,5)),"sil_p95":float(np.percentile(vals,95))})
pd.DataFrame(cal).to_csv(f"{R}/silhouette_fst_calibration.tsv",sep="\t",index=False); out["sil_calibration"]=cal
json.dump(out,open(f"{R}/structure_summary.json","w"),indent=1,default=str); print(json.dumps(out,indent=1,default=str)[:4000])
