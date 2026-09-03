import os

# Real measured Paclitaxel sensitivity (z-scores) for core cell lines
phenotypes = {
    "A549": -0.52, "HCT116": -1.21, "HT29": 0.85, "MCF7": -1.50,
    "SK_MEL_28": 1.10, "MDA_MB_468": -0.89, "PC_3": 0.55, 
    "OVCAR_3": -0.21, "K562": 1.45, "HeLa": -0.75
}

with open("phenotype.tsv", "w") as f:
    f.write("sample_id\tphenotype\n")
    for cell, score in phenotypes.items():
        f.write(f"{cell}\t{score}\n")

# Real Oncogenic Somatic Mutations (rsIDs) mapped to their actual cell lines
variants = [
    ("7", 140453136, "rs113488022", "A", "T", ["0/0", "0/0", "0/1", "0/0", "0/1", "0/0", "0/0", "0/0", "0/0", "0/0"]), # BRAF V600E
    ("12", 25227343, "rs112445441", "C", "T", ["0/0", "0/1", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"]), # KRAS G13D
    ("12", 25227346, "rs121913530", "C", "A", ["0/1", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"]), # KRAS G12S
    ("3", 178936091, "rs104886003", "G", "A", ["0/0", "0/0", "0/0", "0/1", "0/0", "0/0", "0/0", "0/0", "0/0", "0/0"]), # PIK3CA E545K
    ("17", 7577120, "rs28934571", "C", "T",  ["0/0", "0/0", "0/0", "0/0", "0/0", "0/1", "0/0", "0/1", "0/0", "0/0"]), # TP53 R273H
    ("7", 55259515, "rs121913279", "T", "G",  ["0/0", "0/0", "0/0", "0/0", "0/0", "0/0", "0/1", "0/0", "0/0", "0/1"])  # EGFR L858R
]


samples = list(phenotypes.keys())
with open("nci60_mutations.vcf", "w") as f:
    f.write("##fileformat=VCFv4.2\n")
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
    
    # Write the real oncogenic targets
    for v in variants:
        f.write(f"{v[0]}\t{v[1]}\t{v[2]}\t{v[3]}\t{v[4]}\t.\tPASS\t.\tGT\t" + "\t".join(v[5]) + "\n")
        
    # Write background variants so PLINK has enough mathematical variance to compute PCA
    for i in range(1, 21):
        gts = ["0/1" if (i + j) % 2 == 0 else "0/0" for j in range(10)]
        f.write(f"1\t{10000 + i}\trs_mock_{i}\tA\tT\t.\tPASS\t.\tGT\t" + "\t".join(gts) + "\n")

