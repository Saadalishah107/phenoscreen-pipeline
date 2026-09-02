#!/usr/bin/env python3
"""
Generate a tiny synthetic dataset (VCF, phenotype file, mock COCONUT SMILES
library) so the pipeline can be smoke-tested end-to-end without real data.
Run: python3 generate_test_data.py
"""
import random
from pathlib import Path

random.seed(42)
HERE = Path(__file__).parent

N_SAMPLES = 20
N_VARIANTS = 200

# ---- VCF ---------------------------------------------------------------
vcf_lines = [
    "##fileformat=VCFv4.2",
    '##INFO=<ID=AF,Number=1,Type=Float,Description="Allele Frequency">',
    "##contig=<ID=chr1>",
]
samples = [f"SAMPLE{i:03d}" for i in range(N_SAMPLES)]
header = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples)
vcf_lines.append(header)

bases = ["A", "C", "G", "T"]
for i in range(N_VARIANTS):
    pos = 1000 + i * 500
    ref, alt = random.sample(bases, 2)
    rsid = f"rs{100000+i}"
    af = round(random.uniform(0.05, 0.5), 3)
    genotypes = []
    for _ in samples:
        gt = random.choice(["0/0", "0/1", "1/1"])
        genotypes.append(gt)
    line = f"chr1\t{pos}\t{rsid}\t{ref}\t{alt}\t99\tPASS\tAF={af}\tGT\t" + "\t".join(genotypes)
    vcf_lines.append(line)

(HERE / "sample.vcf").write_text("\n".join(vcf_lines) + "\n")

# ---- Phenotype ------------------------------------------------------------
pheno_lines = ["sample_id\tphenotype"]
for s in samples:
    pheno_lines.append(f"{s}\t{round(random.gauss(50, 10), 2)}")
(HERE / "sample_phenotype.tsv").write_text("\n".join(pheno_lines) + "\n")

# ---- Mock COCONUT natural-product SMILES library -------------------------
smiles_pool = [
    ("caffeine_like", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"),
    ("flavonoid_like", "COc1cc(O)c2c(=O)cc(-c3ccccc3)oc2c1"),
    ("terpenoid_like", "CC1=CCC2CC1C2(C)C"),
    ("alkaloid_like", "CN1CCC23C4Oc5c3c(CC1C2C=CC4O)ccc5O"),
    ("phenolic_like", "Oc1ccc(O)cc1"),
]
with open(HERE / "sample_coconut.smi", "w") as fh:
    for i in range(50):
        name, smi = random.choice(smiles_pool)
        fh.write(f"{smi} COCONUT_{i:03d}_{name}\n")

print("Synthetic test data written to test_data/")
