"""Supplementary Fig. S2. Declared-strain tests from distance-based RDA.
Run from the repository root: python figures/make_figS2.py"""
import pandas as pd, matplotlib.pyplot as plt
import os; os.makedirs("figures/out", exist_ok=True)
D="results/"
ST=pd.read_csv(D+"strain_effect_tests.tsv",sep="\t")
rows=[("A. Tojong(14) vs Halla(3)","index",r"$\it{Hallabeol}$ vs $\it{tojong}$"+"\n(ancestry controlled)"),("A. Tojong(14) vs Halla(3)","index+region",r"$\it{Hallabeol}$ vs $\it{tojong}$"+"\n(+ region controlled)"),
      ("C. declared strain, 4 levels","index","All 4 declared strains\n(ancestry controlled)"),("C. declared strain, 4 levels","region","All 4 declared strains\n(region controlled)"),
      ("C. declared strain, 4 levels","index+region","All 4 declared strains\n(ancestry + region)")]
fig,b=plt.subplots(figsize=(4.2,3.0),dpi=300)
for i,(t,cv,l) in enumerate(rows):
    for dn,c,m in [("SNP","#1f77b4","o"),("SV","#7f7f7f","^")]:
        p=ST[(ST.test==t)&(ST.covariates==cv)&(ST.data==dn)].p.iloc[0]
        b.scatter(p,len(rows)-1-i+(0.12 if dn=="SNP" else -0.12),color=c,marker=m,s=22,label=f"{dn} distance" if i==0 else None)
b.axvline(0.05,ls="--",color="#B2182B",lw=0.9); b.set_ylim(-0.5,4.75); b.text(0.055,4.5,"p = 0.05",color="#B2182B",fontsize=6.5)
b.set_xscale("log"); b.set_xlim(0.005,1.5); b.set_xticks([0.01,0.05,0.1,0.5,1]); b.set_xticklabels(["0.01","0.05","0.1","0.5","1"]); b.minorticks_off(); b.set_yticks(range(len(rows))[::-1]); b.set_yticklabels([r[2] for r in rows],fontsize=6.5)
b.set_xlabel("Permutation p\n(distance-based variance partition)",fontsize=7); b.tick_params(axis="x",labelsize=6.5); b.legend(fontsize=6.5,frameon=False,loc="upper center",bbox_to_anchor=(0.4,-0.25),ncol=2,columnspacing=0.8,handletextpad=0.3)
b.spines[["top","right"]].set_visible(False)
fig.subplots_adjust(left=0.38,right=0.97,top=0.96,bottom=0.3); fig.savefig("figures/out/fig_s2_strain_partition.png",dpi=300); print("ok")
