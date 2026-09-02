#!/usr/bin/env python3
"""Standardize COCONUT/SMILES and calculate reproducible RDKit descriptors."""
import argparse, csv, sys
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors, QED
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog("rdApp.*")

def load(path):
    if path.lower().endswith(".sdf"):
        return [m for m in Chem.SDMolSupplier(path) if m is not None]
    out=[]
    with open(path) as fh:
        for line in fh:
            if not line.strip(): continue
            parts=line.split()
            m=Chem.MolFromSmiles(parts[0])
            if m:
                m.SetProp("_Name", parts[1] if len(parts)>1 else f"compound_{len(out)}")
                out.append(m)
    return out

def std(m):
    try:
        m=rdMolStandardize.Cleanup(m)
        m=rdMolStandardize.FragmentParent(m)
        return rdMolStandardize.Uncharger().uncharge(m)
    except Exception:
        return m

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",required=True); p.add_argument("--radius",type=int,default=2)
    p.add_argument("--bits",type=int,default=2048); p.add_argument("--out_desc",required=True)
    p.add_argument("--out_plot",required=True); a=p.parse_args()
    mols=load(a.input)
    if not mols: raise SystemExit("No valid molecules.")
    rows=[]
    for i,m in enumerate(mols):
        m=std(m); cid=m.GetProp("_Name") if m.HasProp("_Name") else f"compound_{i}"
        rows.append({
            "compound_id":cid,"smiles":Chem.MolToSmiles(m),"mol_weight":Descriptors.MolWt(m),
            "logp":Crippen.MolLogP(m),"h_bond_donors":rdMolDescriptors.CalcNumHBD(m),
            "h_bond_acceptors":rdMolDescriptors.CalcNumHBA(m),"tpsa":rdMolDescriptors.CalcTPSA(m),
            "rotatable_bonds":rdMolDescriptors.CalcNumRotatableBonds(m),
            "aromatic_rings":rdMolDescriptors.CalcNumAromaticRings(m),
            "heavy_atoms":m.GetNumHeavyAtoms(),"ring_count":rdMolDescriptors.CalcNumRings(m),
            "formal_charge":Chem.GetFormalCharge(m),"qed":QED.qed(m),
            "fraction_csp3":rdMolDescriptors.CalcFractionCSP3(m)
        })
    df=pd.DataFrame(rows); df.to_csv(a.out_desc,sep="\t",index=False)
    fig,axs=plt.subplots(2,2,figsize=(10,8))
    for ax,col,title in zip(axs.flat,["mol_weight","logp","tpsa","qed"],
                            ["Molecular weight","LogP","TPSA","QED"]):
        ax.hist(df[col],bins=min(20,max(5,len(df)//2))); ax.set_title(title)
    fig.tight_layout(); fig.savefig(a.out_plot,dpi=150); plt.close(fig)

if __name__=="__main__": main()
