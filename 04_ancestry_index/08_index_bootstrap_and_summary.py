#!/usr/bin/env python3
# 19군락 지수: 블록 부트스트랩 SE + 기준 개체(중국 18, 일본 18) 지수(같은 35,825 사이트) + 중국 기준 개체 중앙 깊이(9.5×)로 깊이 맞춤
import numpy as np, pandas as pd, json, os
S=os.environ["KH_RESULTS"]; R=os.environ["KH_RESULTS"]
def read(path,names):
    n=len(names); K=[];Rr=[];A=[]
    for ln in open(path):
        p=ln.rstrip("\n").split("\t")
        if len(p)!=4+n: continue
        r=[];a=[]
        for c in p[4:]:
            try: x,y=c.split(",")[:2]; r.append(int(x)); a.append(int(y))
            except Exception: r.append(0); a.append(0)
        K.append((p[0],int(p[1]))); Rr.append(r); A.append(a)
    return K,np.array(Rr,float),np.array(A,float)
def nm(f): return [os.path.basename(s).replace("_dedup.bam","").replace(".bam","") for s in open(f).read().strip().split("\t") if s]
cs=nm(f"{R}/china_ref/geno_samples.txt"); K,Rc,Ac=read(f"{R}/china_ref/ad_matrix.tsv",cs)
DP=Rc+Ac; ok=DP>=8; V=np.where(ok,Ac/np.where(DP>0,DP,1),np.nan)
site_ok=ok.mean(1)>=0.9; sm=np.nanmean(np.where(ok,V,np.nan),1); poly=site_ok&(sm>0.05)&(sm<0.95)
kk=[k for k,m in zip(K,poly) if m]; Vp=V[poly]; kc=np.nanmean(Vp,0)
blk=np.array([f"{c}:{p//1_000_000}" for c,p in kk]); ub,inv=np.unique(blk,return_inverse=True)
Sm=np.zeros((len(ub),Vp.shape[1])); Nn=np.zeros_like(Sm); np.add.at(Sm,inv,np.nan_to_num(Vp)); np.add.at(Nn,inv,(~np.isnan(Vp)).astype(float))
rng=np.random.default_rng(20261015); bs=[]
for _ in range(1000):
    b=rng.integers(0,len(ub),len(ub)); bs.append(Sm[b].sum(0)/Nn[b].sum(0))
se=np.std(bs,axis=0,ddof=1)
COL=pd.DataFrame({"colony":cs,"KC":np.round(kc,5),"KC_se":np.round(se,5)}); COL.to_csv(f"{R}/two_axes_colonies.tsv",sep="\t",index=False)
POLY=set(kk)
# 기준 개체
meta=pd.read_csv(f"{S}/china_ref/calib/calib_meta.tsv",sep="\t").set_index("sample").group.to_dict()
cal=nm(f"{S}/china_ref/calib/calib_samples.txt"); K2,R2,A2=read(f"{S}/china_ref/calib/calib_ad.tsv",cal)
jp=nm(f"{S}/japan_ref/jp_on_kc/samples.txt"); K3,R3,A3=read(f"{S}/japan_ref/jp_on_kc/ad.tsv",jp)
jpsel=set(pd.read_csv(f"{S}/japan_ref/japan_anchor_selected.tsv",sep="\t").Unified_ID)
def idx(K_,R_,A_,cols):
    sel=np.array([k in POLY for k in K_]); D=R_[sel]+A_[sel]; v=np.where(D>=8,A_[sel]/np.maximum(D,1),np.nan)
    return {c:float(np.nanmean(v[:,i])) for i,c in enumerate(cols)}, R_[sel], A_[sel]
ic,Rcs,Acs=idx(K2,R2,A2,cal); ij,_,_=idx(K3,R3,A3,jp)
rows=[{"sample":s,"group":meta[s],"KC":round(v,5)} for s,v in ic.items() if str(meta.get(s,"")).startswith("China")]
rows+=[{"sample":s,"group":"Japan","KC":round(v,5)} for s,v in ij.items() if s in jpsel]
AN=pd.DataFrame(rows); AN.to_csv(f"{R}/two_axes_anchors.tsv",sep="\t",index=False)
cn=AN[AN.group.str.startswith("China")]; 
# 깊이 맞춤: 중국 기준 개체 중앙 깊이(9.5×) (군락 + 중국 개체, 이항 솎기 20회)
cidx=[i for i,s in enumerate(cal) if str(meta.get(s,"")).startswith("China")]
tgt=float(np.median((Rcs[:,cidx]+Acs[:,cidx]).mean(0)))
rng=np.random.default_rng(20261015); out=[]
Rk=Rc[poly]; Ak=Ac[poly]
for rep in range(20):
    def thin(Rx,Ax):
        D=Rx+Ax; m=D.mean(0); f=np.minimum(1,tgt/m); a=rng.binomial(Ax.astype(int),f); r=rng.binomial(Rx.astype(int),f); Dt=a+r
        return np.nanmean(np.where(Dt>=3,a/np.maximum(Dt,1),np.nan),0)
    out.append((thin(Rk,Ak), thin(Rcs[:,cidx],Acs[:,cidx])))
colT=np.mean([o[0] for o in out],0); cnT=np.mean([o[1] for o in out],0)
summ={"n_callable90":int(site_ok.sum()),"n_poly":int(poly.sum()),"colony_range":[float(kc.min()),float(kc.max())],
 "china_median":float(cn.KC.median()),"china_n_above_max":int((cn.KC>kc.max()).sum()),
 "japan_median":float(AN[AN.group=="Japan"].KC.median()),"japan_range":[float(AN[AN.group=="Japan"].KC.min()),float(AN[AN.group=="Japan"].KC.max())],
 "NE":{s:round(v,4) for s,v in zip(cn["sample"],cn.KC) if meta[s] in ("China:Jilin","China:Liaoning")},
 "depthmatch_target":tgt,"depthmatch_china_above_colony_max":int((cnT>colT.max()).sum()),"se_range":[float(se.min()),float(se.max())]}
json.dump(summ,open(f"{R}/index_summary.json","w"),indent=1); print(json.dumps(summ,indent=1))
