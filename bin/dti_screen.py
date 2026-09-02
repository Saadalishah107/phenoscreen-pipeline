#!/usr/bin/env python3
"""Sequence/SMILES DTI screening. Mock mode is deterministic integration testing."""
import argparse, math
import numpy as np
import pandas as pd

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--targets",required=True); p.add_argument("--descriptors",required=True)
    p.add_argument("--mode",choices=["mock","deepurpose"],default="mock")
    p.add_argument("--model",default="Morgan_CNN_BindingDB")
    p.add_argument("--batch_size",type=int,default=64)
    p.add_argument("--max_compounds",type=int,default=0)
    p.add_argument("--out",required=True)
    a=p.parse_args()

    t=pd.read_csv(a.targets,sep="\t")
    d=pd.read_csv(a.descriptors,sep="\t")
    if a.max_compounds > 0:
        d=d.head(a.max_compounds).copy()
    if t.empty or d.empty:
        pd.DataFrame(columns=["target_id","compound_id","predicted_affinity"]).to_csv(a.out,sep="\t",index=False)
        return

    if a.mode=="mock":
        vals=[]
        for _,tr in t.iterrows():
            for _,c in d.iterrows():
                score = - (0.002*len(tr.protein_sequence) + 0.01*c.tpsa
                           - 0.05*abs(c.logp-2.0) - 0.001*c.mol_weight)
                vals.append((tr.rsid,c.compound_id,score))
        df=pd.DataFrame(vals,columns=["target_id","compound_id","predicted_affinity"])
        df.pivot(index="target_id",columns="compound_id",values="predicted_affinity").to_csv(a.out,sep="\t")
        return

    try:
        from DeepPurpose import utils, DTI as models
    except Exception as e:
        raise SystemExit(f"DeepPurpose import failed: {e}")

    model=models.model_pretrained(model=a.model)
    drugs=d["smiles"].astype(str).tolist()
    drug_ids=d["compound_id"].astype(str).tolist()
    drug_enc=a.model.split("_")[0]
    target_enc=a.model.split("_")[1]
    rows=[]
    for _,tr in t.iterrows():
        seq=str(tr.protein_sequence)
        for start in range(0,len(drugs),a.batch_size):
            stop=min(start+a.batch_size,len(drugs))
            xdrug=np.array(drugs[start:stop])
            xtarget=np.array([seq]*(stop-start))
            X=utils.data_process(
                X_drug=xdrug, X_target=xtarget, y=np.zeros(stop-start),
                drug_encoding=drug_enc, target_encoding=target_enc,
                split_method="no_split"
            )
            pred=model.predict(X)
            rows.extend(zip([tr.rsid]*(stop-start),drug_ids[start:stop],pred))
    pd.DataFrame(rows,columns=["target_id","compound_id","predicted_affinity"]).pivot(
        index="target_id",columns="compound_id",values="predicted_affinity").to_csv(a.out,sep="\t")

if __name__=="__main__":
    main()
