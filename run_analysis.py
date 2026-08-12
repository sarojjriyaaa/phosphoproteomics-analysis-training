import pandas as pd, numpy as np
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
from pathlib import Path
ROOT=Path(__file__).parent; df=pd.read_csv(ROOT/"Input_Data/phosphosite_matrix_training.csv")
out=ROOT/"Statistics"; fig=ROOT/"Figures"; out.mkdir(exist_ok=True); fig.mkdir(exist_ok=True)
qc=df.groupby(["SampleID","Group"]).agg(Median=("Phospho.Abundance","median"),Missing_Percent=("Phospho.Abundance",lambda x:x.isna().mean()*100)).reset_index()
qc.to_csv(out/"sample_QC_metrics.csv",index=False)
rows=[]
for (p,s),x in df.groupby(["Protein","Site"]):
    c=x.loc[x.Group=="Control","Phospho.Abundance"]; d=x.loc[x.Group=="Disease","Phospho.Abundance"]
    rows.append([p,s,d.mean()-c.mean(),ttest_ind(d,c,equal_var=False).pvalue,x["Localization.Probability"].mean()])
r=pd.DataFrame(rows,columns=["Protein","Site","log2FC_Disease_vs_Control","p_value","Mean_Localization_Probability"])
r["rank"]=r.p_value.rank(method="first"); r["FDR_BH"]=(r.p_value*len(r)/r["rank"]).clip(upper=1)
r["Significant"]=(r["FDR_BH"]<.05)&(r["log2FC_Disease_vs_Control"].abs()>=.5)&(r["Mean_Localization_Probability"]>=.75)
r.to_csv(out/"differential_phosphosites.csv",index=False)
plt.figure(figsize=(7,5)); plt.scatter(r.log2FC_Disease_vs_Control,-np.log10(r.p_value))
plt.axvline(.5,ls="--"); plt.axvline(-.5,ls="--"); plt.axhline(-np.log10(.05),ls="--")
plt.xlabel("log2 fold change"); plt.ylabel("-log10(p-value)"); plt.tight_layout()
plt.savefig(fig/"phosphosite_volcano.png",dpi=160); plt.close()
print("Project 2 analysis complete.")
