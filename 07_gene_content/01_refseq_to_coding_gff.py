#!/usr/bin/env python3
# RefSeq RS_2026_08 GFF → coding gff3 (gene/mRNA/exon/CDS), gene2name, 화학수용체 계열표
import gzip, re, sys, collections, pandas as pd, os
D=os.environ["KH_REFSEQ_DIR"]
SRC=f"{D}/GCF_029169275.2_AcerK_1.0_genomic.gff.gz"   # RefSeq RS_2026_08 GFF (NCBI)
def attrs(s): return dict(kv.split("=",1) for kv in s.split(";") if "=" in kv)
genes={}; mrna={}; feats=collections.defaultdict(list)
for ln in gzip.open(SRC,"rt"):
    if ln[0]=="#": continue
    p=ln.rstrip("\n").split("\t")
    if len(p)<9: continue
    a=attrs(p[8])
    if p[2]=="gene" and a.get("gene_biotype")=="protein_coding":
        gid="RS_"+re.search(r"GeneID:(\d+)",a["Dbxref"]).group(1)
        genes[a["ID"]]=(gid,p,a)
    elif p[2]=="mRNA" and a.get("Parent") in genes:
        mrna[a["ID"]]=(a["Parent"],p)
    elif p[2] in ("exon","CDS") and a.get("Parent") in mrna:
        feats[a["Parent"]].append(p)
# CDS 가 있는 mRNA 만 유지 
keep=collections.defaultdict(list)
for mid,(gref,p) in mrna.items():
    if any(f[2]=="CDS" for f in feats[mid]): keep[gref].append(mid)
out=open(f"{D}/RS2026.coding.gff3","w"); out.write("##gff-version 3\n#!annotation RefSeq GCF_029169275.2-RS_2026_08 (converted)\n")
g2n=[]
for gref,mids in keep.items():
    gid,p,a=genes[gref]
    out.write("\t".join(p[:8]+[f"ID={gid};Name={a.get('gene','')}"])+"\n")
    g2n.append((gid,a.get("gene",""),a.get("description","").replace("%2C",",")))
    for i,mid in enumerate(sorted(mids),1):
        tid=f"{gid}.t{i}"; mp=mrna[mid][1]
        out.write("\t".join(mp[:8]+[f"ID={tid};Parent={gid};old_id={mid}"])+"\n")
        for f in sorted(feats[mid],key=lambda f:(f[2],int(f[3]))):
            out.write("\t".join(f[:8]+[f"Parent={tid}"])+"\n")
out.close()
G=pd.DataFrame(g2n,columns=["gene_id","gene_symbol","description"])
G.to_csv(f"{D}/RS2026_gene2name.tsv",sep="\t",index=False)
# 화학수용체 계열: RefSeq description 기준 
PAT=[("OBP",r"odorant[- ]binding protein|general odorant-binding|pheromone-binding"),
     ("CSP",r"chemosensory protein|ejaculatory bulb-specific protein"),
     ("GR", r"gustatory receptor"),
     ("IR", r"ionotropic receptor"),
     ("OR", r"odorant receptor")]
def fam(d):
    for f,pt in PAT:
        if re.search(pt,d,re.I): return f
    return None
G["family"]=G.description.map(fam)
G[G.family.notna()].rename(columns={"gene_id":"gene_id"})[["gene_id","family","gene_symbol","description"]]\
  .to_csv(f"{D}/RS2026_chemo_table.tsv",sep="\t",index=False)
print("genes", len(G), "chromosomal", sum(1 for g in keep if genes[g][1][0].startswith("NC_0838")))
print(G.family.value_counts().to_dict())
