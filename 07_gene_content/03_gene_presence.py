#!/usr/bin/env python
# 읽기 기반 유전자 보유 분류. 유전자 CDS 길이가중 평균 깊이 / 군락 전장(16 염색체) 평균 깊이
# 결손: 비율 < 0.10 (양쪽 대립유전자 모두 없음에 해당), 반감: 0.10–0.35 (참고용)
import pandas as pd, numpy as np, gzip, json, os
P=os.environ["KH_ROOT"]; O=os.environ["KH_RESULTS"]+f"/read_gene_presence"
bams=[l.strip() for l in open(f"{O}/bams.txt") if l.strip()]
S=[os.path.basename(b).replace(".bam","") for b in bams]; assert len(S)==19
ratio={}; gmean={}
for s in S:
    sm=pd.read_csv(f"{O}/md/{s}.mosdepth.summary.txt",sep="\t")
    ch=sm[sm.chrom.str.startswith("NC_0838") & ~sm.chrom.str.endswith("_region")]
    gm=(ch["bases"].sum()/ch["length"].sum()); gmean[s]=round(gm,2)
    r=pd.read_csv(f"{O}/md/{s}.regions.bed.gz",sep="\t",header=None,names=["c","s","e","g","d"])
    r["L"]=r.e-r.s; r["dl"]=r.d*r.L
    g=r.groupby("g")[["dl","L"]].sum(); ratio[s]=(g.dl/g.L)/gm
R=pd.DataFrame(ratio); n=R.shape[1]
pres=(R>=0.10).astype(int); half=((R>=0.10)&(R<0.35)).sum(axis=0)
cnt=pres.sum(axis=1)
cls=pd.cut(cnt,[-1,0,1,n-2,n-1,n],labels=["absent_all","private","dispensable","softcore","core"])
cc={str(k):int(v) for k,v in cls.value_counts().items()}
R.round(3).to_csv(f"{O}/gene_depth_ratio.tsv",sep="\t"); pres.to_csv(f"{O}/gene_presence_read.tsv",sep="\t")
chemo=pd.read_csv(os.environ["KH_REFSEQ_DIR"] + "/RS2026_chemo_table.tsv",sep="\t")
tab=pd.DataFrame({"gene":R.index,"n_present":cnt.values,"class":cls.values})
tab["family"]=tab.gene.map(dict(zip(chemo.gene_id.astype(str),chemo.family))).fillna("background")
fam=[]
for f in ["OR","GR","OBP","IR","CSP","background"]:
    s=tab[tab.family==f]
    fam.append({"family":f,"n_gene":len(s),"absent_in_>=1":int((s["class"]!="core").sum()),"dispensable+private":int(s["class"].isin(["dispensable","private"]).sum())})
tab.to_csv(f"{O}/gene_class_read.tsv",sep="\t",index=False)
# 결손 기준 민감도: 0.05–0.30 에서 모든 군락에 있는 유전자 비율
sens={str(t):round(100*int(((R>=t).sum(axis=1)==n).sum())/len(R),2) for t in (0.05,0.10,0.20,0.30)}
absent_per=(pres==0).sum(axis=0).sort_values(ascending=False)
res={"n_colony":n,"genome_mean_depth":gmean,"class_counts":cc,"pct_core":round(100*cc.get("core",0)/len(R),2),
     "n_noncore":int(len(R)-cc.get("core",0)),"absent_per_colony_top":{k:int(v) for k,v in absent_per.head(6).items()},
     "absent_per_colony_median":float(absent_per.median()),"half_depth_genes_per_colony_median":float(half.median()),
     "family":fam,"pct_core_by_threshold":sens}
json.dump(res,open(f"{O}/summary.json","w"),indent=1,default=str); print(json.dumps(res,indent=1,default=str))
