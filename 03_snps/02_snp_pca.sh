#!/bin/bash
source "$(dirname "$0")/../config.sh"
set -u
O=$KH_RESULTS/snpstruct
PL=plink2
$PL --bfile $O/thin2000 --allow-extra-chr --export A --out $O/thin2000 > $O/export.log 2>&1
echo "raw 내보냄: $(ls -la $O/thin2000.raw | awk '{print $5}') bytes"
python3 - <<'PY'
import numpy as np, pandas as pd, json, os
O=os.environ["KH_RESULTS"]+"/snpstruct"
# plink --export A: FID IID PAT MAT SEX PHENO 뒤에 SNP 당 0/1/2 (NA 가능)
hdr=open(f"{O}/thin2000.raw").readline().split()
nsnp=len(hdr)-6
iids=[]; rows=[]
for ln in open(f"{O}/thin2000.raw").readlines()[1:]:
    p=ln.split()
    iids.append(p[1])
    rows.append([np.nan if x=="NA" else float(x) for x in p[6:]])
G=np.array(rows)
print(f"유전형 행렬 {G.shape[0]} 개체 × {G.shape[1]:,} SNP · 결측 {np.isnan(G).mean()*100:.2f}%")
# 패널 자체 빈도로 평균대입 후 표준화 (plink2 가 암묵적으로 거부하는 것을 명시적으로 수행)
mu=np.nanmean(G,axis=0)
G=np.where(np.isnan(G), mu, G)
p=mu/2.0
keep=(p>0.01)&(p<0.99)
G=G[:,keep]; p=p[keep]
print(f"표준화 대상 {G.shape[1]:,} SNP")
Z=(G-2*p)/np.sqrt(2*p*(1-p))
U,S,Vt=np.linalg.svd(Z-Z.mean(0), full_matrices=False)
var=(S**2)/np.sum(S**2)
PCs=U*S
out=pd.DataFrame(PCs[:,:10], columns=[f"snpPC{i+1}" for i in range(10)])
out.insert(0,"sample",iids)
out.to_csv(f"{O}/snp_pca.tsv", sep="\t", index=False)
json.dump({"n_ind":int(G.shape[0]),"n_snp":int(G.shape[1]),
           "var_explained":[round(float(v),5) for v in var[:10]]},
          open(f"{O}/snp_pca_meta.json","w"), indent=1)
print("설명분산 PC1-5:", [f"{v*100:.2f}%" for v in var[:5]])
PY
echo PCA_DONE
