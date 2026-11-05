
import numpy as np, pandas as pd, json, itertools, os
B=os.environ["KH_RESULTS"] + "/"
OUT=B+"summary_tables/"
def load(p):
    rows=[l.rstrip("\n").split("\t") for l in open(p)]
    return rows[0][1:], np.array([[float(x) for x in r[1:]] for r in rows[1:]])
AP={"Ac_GP1":"GP-a","Ac_GP2":"GP-a","Ac_GP3":"GP-b","Ac_GP4":"GP-c","Ac_GP5":"GP-d","Ac":"INU",
 "CWN2_AC1":"CWN-a","CWN2_AC2":"CWN-a","CWN2_AC3":"CWN-a","CWN4_AC1":"CWN-b",
 "IC_AC_1":"IC-a","IC_AC_2":"IC-a","WJ_AC_1":"WJ-a","WJ_AC_2":"WJ-a","WJ_AC_3":"WJ-a","WJ_AC_4":"WJ-a",
 "WJ_AC_5":"WJ-a","YP_AC2":"YP-a","YP_AC_3":"YP-a"}
res=[]
for tag,p in [("SNP",B+"snpstruct/snp_dist.tsv"),("SV",B+"sv_jaccard_distance.tsv")]:
    ids,D=load(p); g=np.array([AP[s] for s in ids])
    # AMOVA 식 분산성분: 전체 SS = sum d^2/(2N); 군내 SS = 봉장별 sum d^2/(2 n_k)
    N=len(ids); SS_T=(D**2).sum()/(2*N)
    SS_W=0.0; dfW=0
    multi=[a for a in set(g) if (g==a).sum()>1]
    for a in set(g):
        m=g==a; n=m.sum()
        if n>1: SS_W += (D[np.ix_(m,m)]**2).sum()/(2*n); dfW += n-1
    SS_A=SS_T-SS_W; dfA=len(set(g))-1
    n_bar=(N-sum((g==a).sum()**2 for a in set(g))/N)/dfA
    MS_W=SS_W/dfW; MS_A=SS_A/dfA
    var_A=max((MS_A-MS_W)/n_bar,0.0); var_W=MS_W
    res.append({"data":tag,"within_apiary_pct":round(100*var_W/(var_A+var_W),1),
      "among_apiary_pct":round(100*var_A/(var_A+var_W),1),"MS_within":round(MS_W,4),"MS_among":round(MS_A,4),
      "n_apiaries":len(set(g)),"n_multi_colony":len(multi)})
    print(json.dumps(res[-1], ensure_ascii=False))
pd.DataFrame(res).to_csv(OUT+"apiary_variance_components.tsv", sep="\t", index=False)

print("DONE")
