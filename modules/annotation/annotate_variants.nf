process ANNOTATE_VARIANTS {
    tag "Variant annotation"
    label 'annotation'

    input:
    path candidate_loci
    path annotation_tsv

    output:
    path 'variant_annotation.tsv', emit: annotation
    path 'protein_sequences.fasta', emit: protein_fasta

    script:
    def offline_flag = params.offline_annotation ? "--offline" : ""
    """
    python3 ${projectDir}/bin/annotate_variants.py \
      --candidate_loci ${candidate_loci} \
      --annotation_tsv ${annotation_tsv} \
      --out_annotation variant_annotation.tsv \
      --out_fasta protein_sequences.fasta \
      --species ${params.species} \
      --max_variants ${params.max_annotation_variants} \
      ${offline_flag}
    """
}
