#!/usr/bin/env python3
"""Fetch AlphaFold DB structures with ESMFold as fallback; mock backend is for CI."""
import argparse, csv, time
from pathlib import Path
import requests
from Bio import SeqIO
from Bio.PDB import PDBParser

def parse():
    p=argparse.ArgumentParser()
    p.add_argument("--annotation", required=True)
    p.add_argument("--fasta", required=True)
    p.add_argument("--backend", choices=["alphafold","esmfold","mock"], default="alphafold")
    p.add_argument("--outdir", required=True)
    p.add_argument("--max_retries", type=int, default=3)
    p.add_argument("--retry_delay", type=int, default=3)
    return p.parse_args()

def mock_pdb(seq_id, seq, path):
    # Minimal CA trace: enough to archive a deterministic structure artifact.
    aa3={"A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE","G":"GLY","H":"HIS",
         "I":"ILE","K":"LYS","L":"LEU","M":"MET","N":"ASN","P":"PRO","Q":"GLN",
         "R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR"}
    lines=[]; serial=1
    for i, aa in enumerate(seq[:200], 1):
        res=aa3.get(aa,"GLY"); x=(i-1)*3.8
        for atom_name, dx, elem in [("N",-1.2,"N"),("CA",0.0,"C"),("C",1.2,"C")]:
            lines.append(f"ATOM  {serial:5d} {atom_name:>4} {res:>3} A{i:4d}    {x+dx:8.3f}{0:8.3f}{0:8.3f}  1.00 50.00           {elem:>2}")
            serial += 1
    lines.append("END")
    path.write_text("\n".join(lines)+"\n")

def fetch(url, out, retries, delay):
    for i in range(retries):
        try:
            r=requests.get(url, timeout=60)
            if r.ok:
                out.write_text(r.text); return True
        except requests.RequestException:
            pass
        time.sleep(delay*(i+1))
    return False

def main():
    a=parse(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    ann=list(csv.DictReader(open(a.annotation), delimiter="\t"))
    manifest=[]
    for r in ann:
        sid=r["rsid"]; seq=r["protein_sequence"]; uid=r.get("uniprot_id","")
        path=out/f"{sid}.pdb"; status="failed"; source=""; global_confidence=""
        if a.backend=="mock":
            mock_pdb(sid,seq,path); status="ok"; source="mock"
        elif a.backend=="alphafold" and uid:
            try:
                meta = requests.get(
                    f"https://alphafold.ebi.ac.uk/api/prediction/{uid}",
                    headers={"Accept":"application/json"}, timeout=45
                )
                if meta.ok:
                    entries = meta.json()
                    if entries:
                        entries = sorted(entries, key=lambda x: x.get("latestVersion", 0), reverse=True)
                        global_confidence = entries[0].get("globalMetricValue", "")
                        pdb_url = entries[0].get("pdbUrl")
                        if pdb_url and fetch(pdb_url, path, a.max_retries, a.retry_delay):
                            status="ok"; source="AlphaFold DB API"
            except requests.RequestException:
                pass
        if status!="ok" and a.backend in ("esmfold","alphafold"):
            try:
                r=requests.post("https://api.esmatlas.com/foldSequence/v1/pdb/",
                                 data=seq, headers={"Content-Type":"text/plain"}, timeout=180)
                if r.ok:
                    path.write_text(r.text); status="ok"; source="ESMFold API"
            except requests.RequestException:
                pass
        mean_plddt = ""
        if status == "ok":
            try:
                structure = PDBParser(QUIET=True).get_structure(sid, str(path))
                values = [atom.get_bfactor() for atom in structure.get_atoms() if atom.get_name()=="CA"]
                if values: mean_plddt = sum(values)/len(values)
            except Exception:
                pass
        manifest.append({"rsid":sid,"protein_id":r["protein_id"],"uniprot_id":uid,
                         "structure_file":str(path) if status=="ok" else "",
                         "status":status,"source":source,
                         "global_confidence":global_confidence,
                         "mean_ca_bfactor":mean_plddt})
    with open("structure_manifest.tsv","w",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=manifest[0].keys() if manifest else
                         ["rsid","protein_id","uniprot_id","structure_file","status","source",
                          "global_confidence","mean_ca_bfactor"],delimiter="\t")
        w.writeheader(); w.writerows(manifest)

if __name__=="__main__": main()
