import csv

# Using a standard 30-AA string to satisfy DeepPurpose CNN extractors
seq = "MTEYKLVVVGAGGVGKSALTIQLIQNHFVD"
data = [
    ["rs121913530", "12", "25227346", "C", "A", "KRAS", "TX1", "PROT1", "P01116", "missense_variant", "G/S", "12", seq, seq, "nci60_fixture"],
    ["rs28934571", "17", "7577120", "C", "T", "TP53", "TX2", "PROT2", "P04637", "missense_variant", "R/H", "273", seq, seq, "nci60_fixture"],
    ["rs104886003", "3", "178936091", "G", "A", "PIK3CA", "TX3", "PROT3", "P42336", "missense_variant", "E/K", "545", seq, seq, "nci60_fixture"],
    ["rs113488022", "7", "140453136", "A", "T", "BRAF", "TX4", "PROT4", "P15056", "missense_variant", "V/E", "600", seq, seq, "nci60_fixture"],
    ["rs121913279", "7", "55259515", "T", "G", "EGFR", "TX5", "PROT5", "P00533", "missense_variant", "L858R", "858", seq, seq, "nci60_fixture"],
    ["rs112445441", "12", "25227343", "C", "T", "KRAS", "TX6", "PROT6", "P01116", "missense_variant", "G/D", "13", seq, seq, "nci60_fixture"]
]
headers = ["rsid","chrom","pos","ref","alt","gene_symbol","transcript_id","protein_id","uniprot_id","consequence","amino_acids","protein_start","reference_protein_sequence","protein_sequence","annotation_source"]

with open("nci60_annotation.tsv", "w") as f:
    f.write("\t".join(headers) + "\n")
    for row in data:
        f.write("\t".join(row) + "\n")
