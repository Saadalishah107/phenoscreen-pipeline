process PROCESS_CHEM_LIBRARY {
    tag "COCONUT/RDKit library"
    label 'cheminformatics'

    input:
    path coconut

    output:
    path 'molecular_descriptors.tsv', emit: descriptors
    path 'descriptor_distributions.png', emit: plot

    script:
    """
    python3 ${projectDir}/bin/process_coconut.py \
      --input ${coconut} \
      --radius ${params.fingerprint_radius} \
      --bits ${params.fingerprint_bits} \
      --out_desc molecular_descriptors.tsv \
      --out_plot descriptor_distributions.png
    """
}
