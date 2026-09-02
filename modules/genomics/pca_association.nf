process PCA_ASSOCIATION {
    tag "PCA + association"
    label 'genomics'

    input:
    tuple path(pruned_bed), path(pruned_bim), path(pruned_fam)
    path phenotype

    output:
    path 'candidate_loci.tsv', emit: candidate_loci
    path 'association_results.tsv', emit: association_results
    path 'pca_eigenvectors.tsv', emit: pca_eigenvectors
    path 'pca_scree_plot.png', emit: pca_plot
    path 'pca_scatter_plot.png', emit: pca_scatter
    path 'association_manhattan.png', emit: association_plot
    path 'association_qq.png', emit: qq_plot

    script:
    """
    Rscript ${projectDir}/bin/pca_association.R \
      --bfile_prefix ${pruned_bed.baseName} \
      --phenotype ${phenotype} \
      --phenotype_type ${params.phenotype_type} \
      --n_pcs ${params.pca_components} \
      --top_n ${params.candidate_top_n} \
      --fdr ${params.fdr_threshold}
    """
}
