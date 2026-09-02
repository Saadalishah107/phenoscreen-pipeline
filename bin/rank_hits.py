#!/usr/bin/env python3
"""Rank predicted virtual-screening candidates and attach chemical descriptors."""
import argparse, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--affinity_matrix",required=True); p.add_argument("--descriptors",required=True)
    p.add_argument("--top_n",type=int,default=50); p.add_argument("--out_table",required=True)
    p.add_argument("--out_plot",required=True); a=p.parse_args()
    df=pd.read_csv(a.affinity_matrix,sep="\t",index_col=0)
    if df.empty or df.shape[1] == 0:
        pd.DataFrame(columns=["rank","target_id","compound_id","predicted_affinity"]).to_csv(a.out_table,sep="\t",index=False)
        fig,ax=plt.subplots(figsize=(9,6)); ax.set_title("No predicted candidates"); fig.savefig(a.out_plot,dpi=150); plt.close(fig); return
    long=df.reset_index().melt(id_vars=df.index.name or "target_id",
                               var_name="compound_id",value_name="predicted_affinity").dropna()
    desc=pd.read_csv(a.descriptors,sep="\t")
    out=long.merge(desc,on="compound_id",how="left").sort_values("predicted_affinity").head(a.top_n).reset_index(drop=True)
    out.insert(0,"rank",range(1,len(out)+1)); out.to_csv(a.out_table,sep="\t",index=False)
    fig,ax=plt.subplots(figsize=(9,6))
    labels=out["compound_id"].astype(str)+" / "+out.iloc[:,1].astype(str)
    ax.barh(labels,out["predicted_affinity"]); ax.invert_yaxis()
    ax.set_xlabel("Model score (lower is stronger only for compatible models)")
    ax.set_title(f"Top {len(out)} predicted candidates")
    fig.tight_layout(); fig.savefig(a.out_plot,dpi=150); plt.close(fig)

if __name__=="__main__": main()
