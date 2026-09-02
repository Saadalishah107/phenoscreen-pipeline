process DTI_SCREENING {
    tag "DTI screening: ${params.dti_mode}"
    label 'cheminformatics'
    

    input:
    path protein_fasta
    path descriptors

    output:
    path 'affinity_matrix.tsv', emit: affinity_matrix

    script:
    """
    python3 ${projectDir}/bin/dti_screen.py \
      --targets ${protein_fasta} \
      --descriptors ${descriptors} \
      --mode ${params.dti_mode} \
      --model ${params.dti_model} \
      --batch_size ${params.dti_batch_size} \
      --max_compounds ${params.max_compounds} \
      --out affinity_matrix.tsv
    """
}