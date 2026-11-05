#!/usr/bin/env python
# 공개 중국 조립체 4종의 4/4 합의 SNV → 코호트 비의존 대륙 진단부위
import pandas as pd, sys, os
D=os.environ["KH_RESULTS"] + "/china_ref"
labs=["CN_Hainan","CN_Sichuan","CN_Hubei","CN_Qinghai"]
sets={}
for L in labs:
    rows=[]
    for ln in open(f"{D}/{L}.vcf"):
        if ln[0]=="#": continue
        p=ln.split("\t",5)
        if not p[0].startswith("NC_0838"): continue      # 염색체만
        if len(p[3])!=1 or len(p[4])!=1: continue        # SNV만
        if p[3] not in "ACGT" or p[4] not in "ACGT": continue
        rows.append((p[0],int(p[1]),p[3],p[4]))
    sets[L]=pd.DataFrame(rows, columns=["chrom","pos","ref","alt"])
    print(L, len(sets[L]), "SNV", file=sys.stderr)

m=sets[labs[0]]
for L in labs[1:]:
    m=m.merge(sets[L], on=["chrom","pos","ref","alt"], how="inner")
print("4/4 합의 SNV:", len(m), file=sys.stderr)
m=m.sort_values(["chrom","pos"])
m.to_csv(f"{D}/china_diag_sites.tsv", sep="\t", index=False)
m[["chrom","pos"]].to_csv(f"{D}/china_diag_sites.pos", sep="\t", index=False, header=False)
with open(f"{D}/china_diag_alleles.tsv","w") as f:
    for r in m.itertuples(): f.write(f"{r.chrom}\t{r.pos}\t{r.ref}\t{r.alt}\n")
print("저장 완료", flush=True)
