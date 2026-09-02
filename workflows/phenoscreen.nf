nextflow.enable.dsl=2

include { LD_PRUNING } from '../modules/genomics/ld_pruning'
include { PCA_ASSOCIATION } from '../modules/genomics/pca_association'
include { ANNOTATE_VARIANTS } from '../modules/annotation/annotate_variants'
include { PREDICT_STRUCTURE } from '../modules/structure/predict_structure'
include { PROCESS_CHEM_LIBRARY } from '../modules/cheminformatics/process_library'
include { DTI_SCREENING } from '../modules/cheminformatics/dti_screening'
include { RANK_HITS } from '../modules/cheminformatics/rank_hits'

workflow PHENOSCREEN {
    take:
    vcf
    phenotype
    coconut
    annotation_tsv

    main:
    LD_PRUNING(vcf)
    PCA_ASSOCIATION(
        LD_PRUNING.out.pruned_bed,
        phenotype
    )

    ANNOTATE_VARIANTS(
        PCA_ASSOCIATION.out.candidate_loci,
        annotation_tsv
    )

    PREDICT_STRUCTURE(
        ANNOTATE_VARIANTS.out.annotation,
        ANNOTATE_VARIANTS.out.protein_fasta
    )

    PROCESS_CHEM_LIBRARY(coconut)

    DTI_SCREENING(
        ANNOTATE_VARIANTS.out.annotation,
        PROCESS_CHEM_LIBRARY.out.descriptors
    )

    RANK_HITS(DTI_SCREENING.out.affinity_matrix, PROCESS_CHEM_LIBRARY.out.descriptors)

    emit:
    candidate_loci      = PCA_ASSOCIATION.out.candidate_loci
    association_results = PCA_ASSOCIATION.out.association_results
    pca_eigenvectors    = PCA_ASSOCIATION.out.pca_eigenvectors
    annotation          = ANNOTATE_VARIANTS.out.annotation
    protein_fasta       = ANNOTATE_VARIANTS.out.protein_fasta
    structure_manifest  = PREDICT_STRUCTURE.out.structure_manifest
    chem_descriptors    = PROCESS_CHEM_LIBRARY.out.descriptors
    affinity_matrix     = DTI_SCREENING.out.affinity_matrix
    top_hits            = RANK_HITS.out.top_hits
    pca_plot            = PCA_ASSOCIATION.out.pca_plot
    pca_scatter         = PCA_ASSOCIATION.out.pca_scatter
    association_plot    = PCA_ASSOCIATION.out.association_plot
    qq_plot             = PCA_ASSOCIATION.out.qq_plot
    chem_plot           = PROCESS_CHEM_LIBRARY.out.plot
    hit_plot            = RANK_HITS.out.plot
}
