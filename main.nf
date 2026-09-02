#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PHENOSCREEN } from './workflows/phenoscreen.nf'

def validate_params() {
    if (!params.vcf) {
        error 'Required parameter missing: --vcf'
    }

    if (!params.phenotype) {
        error 'Required parameter missing: --phenotype'
    }

    if (!params.coconut) {
        error 'Required parameter missing: --coconut'
    }

    if (!(params.phenotype_type in ['continuous', 'binary'])) {
        error "--phenotype_type must be 'continuous' or 'binary'"
    }

    if (!(params.struct_backend in ['alphafold', 'esmfold', 'mock'])) {
        error "--struct_backend must be 'alphafold', 'esmfold', or 'mock'"
    }

    if (!(params.dti_mode in ['deepurpose', 'mock'])) {
        error "--dti_mode must be 'deepurpose' or 'mock'"
    }
}

workflow {
    main:

    validate_params()

    vcf_ch = channel.fromPath(
        params.vcf,
        checkIfExists: true
    )

    phenotype_ch = channel.fromPath(
        params.phenotype,
        checkIfExists: true
    )

    coconut_ch = channel.fromPath(
        params.coconut,
        checkIfExists: true
    )

    annotation_ch = params.annotation_tsv
        ? channel.fromPath(params.annotation_tsv, checkIfExists: true)
        : channel.fromPath(
            "${projectDir}/test_data/empty_annotation.tsv",
            checkIfExists: true
        )

    PHENOSCREEN(
        vcf_ch,
        phenotype_ch,
        coconut_ch,
        annotation_ch
    )

    publish:
    candidate_loci       = PHENOSCREEN.out.candidate_loci
    association_results  = PHENOSCREEN.out.association_results
    pca_eigenvectors     = PHENOSCREEN.out.pca_eigenvectors
    pca_plot             = PHENOSCREEN.out.pca_plot
    pca_scatter          = PHENOSCREEN.out.pca_scatter
    association_plot     = PHENOSCREEN.out.association_plot
    qq_plot              = PHENOSCREEN.out.qq_plot
    annotation           = PHENOSCREEN.out.annotation
    protein_fasta        = PHENOSCREEN.out.protein_fasta
    structure_manifest   = PHENOSCREEN.out.structure_manifest
    chem_descriptors     = PHENOSCREEN.out.chem_descriptors
    chem_plot            = PHENOSCREEN.out.chem_plot
    affinity_matrix      = PHENOSCREEN.out.affinity_matrix
    top_hits             = PHENOSCREEN.out.top_hits
    hit_plot             = PHENOSCREEN.out.hit_plot
}

output {
    candidate_loci {
        path '01_genomics/candidates'
    }

    association_results {
        path '01_genomics/association'
    }

    pca_eigenvectors {
        path '01_genomics/pca'
    }

    pca_plot {
        path '01_genomics/pca'
    }

    pca_scatter {
        path '01_genomics/pca'
    }

    association_plot {
        path '01_genomics/association'
    }

    qq_plot {
        path '01_genomics/association'
    }

    annotation {
        path '02_annotation'
    }

    protein_fasta {
        path '02_annotation'
    }

    structure_manifest {
        path '03_structure'
    }

    chem_descriptors {
        path '04_cheminformatics'
    }

    chem_plot {
        path '04_cheminformatics'
    }

    affinity_matrix {
        path '05_dti'
    }

    top_hits {
        path '06_results'
    }

    hit_plot {
        path '06_results'
    }
}