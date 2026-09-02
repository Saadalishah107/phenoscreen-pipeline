process PREDICT_STRUCTURE {
    tag "Structure: ${params.struct_backend}"
    label 'structure'

    input:
    path annotation
    path protein_fasta

    output:
    path 'structures', emit: structures
    path 'structure_manifest.tsv', emit: structure_manifest

    script:
    """
    python3 ${projectDir}/bin/fetch_structures.py \
      --annotation ${annotation} \
      --fasta ${protein_fasta} \
      --backend ${params.struct_backend} \
      --outdir structures \
      --max_retries ${params.max_retries} \
      --retry_delay ${params.retry_delay}
    """
}
