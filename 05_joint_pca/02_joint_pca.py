#!/usr/bin/env python3
# 결합 PCA: (a) 유사 반수체(무작위 읽기 1개) — HiFi/단읽기 깊이 차이 제거, (b) 대립유전자 비율 (점검용)
import numpy as np, pandas as pd, os, json
from scipy.stats import spearmanr
P=os.environ["KH_RESULTS"]; O=f"{P}/joint_pca_snp"
S=[os.path.basename(s).replace("_dedup.bam","").replace(".bam","") for s in open(f"{O}/samples.txt").read().strip().split("\t")]
n=len(S); R=[];A=[];K=[]
for ln in open(f"{O}/ad.tsv"):
    p=ln.rstrip("\n").split("\t")
    if len(p)!=4+n or "," in p[3]: continue
    r=[];a=[]
    for c in p[4:]:
        try: x,y=c.split(",")[:2]; r.append(int(x)); a.append(int(y))
        except Exception: r.append(0); a.append(0)
    R.append(r);A.append(a);K.append((p[0],int(p[1])))
R=np.array(R,float);A=np.array(A,float);DP=R+A
cn=set(open(f"{O}/cn.ids").read().split()); jp=set(open(f"{O}/jp.ids").read().split())
G=["China" if s in cn else "Japan" if s in jp else "colony" for s in S]
rng=np.random.default_rng(20261015)
ph=np.where(DP>=1,(rng.random(DP.shape)<np.where(DP>0,A/np.maximum(DP,1),0)).astype(float),np.nan)
vf=np.where(DP>=3,A/np.maximum(DP,1),np.nan)
miss=pd.DataFrame({"sample":S,"group":G,"missing_pseudohap":np.isnan(ph).mean(0).round(4),"mean_depth":DP.mean(0).round(2)})
miss.to_csv(f"{O}/per_sample_missing.tsv",sep="\t",index=False)
KC=pd.read_csv(f"{P}/two_axes_colonies.tsv",sep="\t").set_index("colony").KC
summ={"sites_genotyped":len(K)}
for lab,X in [("pseudohap",ph),("vaf",vf)]:
    ok=(~np.isnan(X)).mean(1)>=0.9; f=np.nanmean(X,1); keep=ok&(np.minimum(f,1-f)>=0.05)
    Y=X[keep]; mu=np.nanmean(Y,1,keepdims=True); Y=np.where(np.isnan(Y),mu,Y)-mu
    U,s,Vt=np.linalg.svd(Y.T,full_matrices=False); ev=s**2/np.sum(s**2)
    d=pd.DataFrame({"sample":S,"group":G,**{f"PC{i+1}":np.round(U[:,i]*s[i],4) for i in range(4)}})
    d.to_csv(f"{O}/pca_{lab}.tsv",sep="\t",index=False)
    c=d[d.group=="colony"]; rho=spearmanr(c.PC1,c["sample"].map(KC))
    summ[lab]={"sites":int(keep.sum()),"var_pct":[round(100*x,2) for x in ev[:4]],"rho_PC1_KC":round(float(rho.statistic),4),"p":float(rho.pvalue)}
    print(lab,summ[lab]); print(d.groupby("group")[["PC1","PC2"]].agg(["min","median","max"]).round(2))
json.dump(summ,open(f"{O}/summary.json","w"),indent=1); print(miss.groupby("group")[["missing_pseudohap","mean_depth"]].agg(["min","median","max"]))
# 그림용 방향 고정: 군락이 PC2 양(+)쪽, 중국 본토 개체가 PC1 양(+)쪽이 되도록 부호만 맞춘다 (Fig. 3A)
d=pd.read_csv(f"{O}/pca_pseudohap.tsv",sep="\t")
if d[d.group=="colony"].PC2.median()<0: d["PC2"]=-d["PC2"]
if d[d.group=="China"].PC1.median()<0: d["PC1"]=-d["PC1"]
d.to_csv(f"{O}/pca_pseudohap_oriented.tsv",sep="\t",index=False)
