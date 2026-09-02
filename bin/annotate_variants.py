#!/usr/bin/env python3
"""Annotate candidate variants and resolve protein sequences.

Production mode uses Ensembl REST for consequence annotation and protein
sequence retrieval. Test/offline mode accepts a deterministic TSV fixture.
Only protein-altering consequences are promoted to downstream structure/DTI.
"""
import argparse, csv, sys, time
from pathlib import Path
import requests

PROTEIN_ALTERING = {"missense_variant", "stop_gained"}

def args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidate_loci", required=True)
    p.add_argument("--annotation_tsv")
    p.add_argument("--out_annotation", required=True)
    p.add_argument("--out_fasta", required=True)
    p.add_argument("--species", default="human")
    p.add_argument("--max_variants", type=int, default=20)
    p.add_argument("--offline", action="store_true")
    return p.parse_args()

def read_candidates(path, n):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[:n]

def read_fixture(path):
    if not path:
        return {}
    with open(path) as fh:
        return {r["rsid"]: r for r in csv.DictReader(fh, delimiter="\t")}

def ensembl_get(url, retries=3):
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    for attempt in range(retries):
        r = requests.get(url, headers=headers, timeout=30)
        if r.ok:
            return r.json()
        if r.status_code in (429, 500, 502, 503, 504):
            time.sleep(2 ** attempt)
            continue
        r.raise_for_status()
    return None

def annotate_live(row, species):
    rsid = row.get("snp_id") or row.get("rsid")
    if not rsid:
        return None
    base = f"https://rest.ensembl.org"
    data = ensembl_get(f"{base}/vep/{species}/id/{rsid}?content-type=application/json")
    if not data:
        return None
    transcripts = []
    for rec in data if isinstance(data, list) else [data]:
        transcripts.extend(rec.get("transcript_consequences", []))
    candidates = [x for x in transcripts if x.get("consequence_terms") and
                  any(c in PROTEIN_ALTERING for c in x["consequence_terms"]) and
                  x.get("protein_id")]
    if not candidates:
        return None
    candidates.sort(key=lambda x: (
        0 if x.get("canonical") else 1,
        0 if "missense_variant" in x.get("consequence_terms", []) else 1
    ))
    tx = candidates[0]
    protein_id = tx["protein_id"]
    seq = ensembl_get(f"{base}/sequence/id/{protein_id}?type=protein;content-type=application/json")
    if not seq or not seq.get("seq"):
        return None
    sequence = seq["seq"]
    protein_start = tx.get("protein_start")
    aa_change = tx.get("amino_acids", "")
    variant_sequence = sequence
    if protein_start:
        pos = int(protein_start)
        parts = aa_change.split("/")
        if "missense_variant" in tx.get("consequence_terms", []) and len(parts) == 2 and 1 <= pos <= len(sequence):
            variant_sequence = sequence[:pos-1] + parts[1] + sequence[pos:]
        elif "stop_gained" in tx.get("consequence_terms", []) and 1 <= pos <= len(sequence):
            variant_sequence = sequence[:pos-1]
    return {
        "rsid": rsid,
        "chrom": row.get("chrom", row.get("chromosome", row.get("#CHROM", ""))),
        "pos": row.get("pos", row.get("position", "")),
        "ref": row.get("ref", ""),
        "alt": row.get("alt", ""),
        "gene_symbol": tx.get("gene_symbol", ""),
        "transcript_id": tx.get("transcript_id", ""),
        "protein_id": protein_id,
        "uniprot_id": tx.get("swissprot", ""),
        "consequence": ",".join(tx.get("consequence_terms", [])),
        "amino_acids": aa_change,
        "protein_start": protein_start or "",
        "reference_protein_sequence": sequence,
        "protein_sequence": variant_sequence,
        "annotation_source": "Ensembl VEP REST",
    }

def main():
    a = args()
    candidates = read_candidates(a.candidate_loci, a.max_variants)
    fixture = read_fixture(a.annotation_tsv)
    records = []
    for row in candidates:
        rsid = row.get("snp_id") or row.get("rsid")
        rec = fixture.get(rsid) if fixture else None
        if rec:
            rec = dict(rec)
            rec["annotation_source"] = "local test fixture"
        elif a.offline:
            print(f"[warn] no offline annotation for {rsid}", file=sys.stderr)
            continue
        else:
            try:
                rec = annotate_live(row, a.species)
            except Exception as e:
                print(f"[warn] annotation failed for {rsid}: {e}", file=sys.stderr)
                rec = None
        if rec and rec.get("protein_sequence"):
            records.append(rec)

    fields = [
        "rsid","chrom","pos","ref","alt","gene_symbol","transcript_id",
        "protein_id","uniprot_id","consequence","amino_acids","protein_start",
        "reference_protein_sequence","protein_sequence","annotation_source"
    ]
    with open(a.out_annotation, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(records)

    with open(a.out_fasta, "w") as fh:
        for r in records:
            fh.write(f">{r['rsid']}|{r['protein_id']}|{r['gene_symbol']}\n")
            seq = r["protein_sequence"]
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

    if not records:
        print("No protein-altering candidates could be resolved.", file=sys.stderr)
        Path(a.out_fasta).touch()
    print(f"Resolved {len(records)} protein targets.")

if __name__ == "__main__":
    main()
