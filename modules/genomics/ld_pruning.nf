process LD_PRUNING {
    tag "LD pruning: ${vcf.simpleName}"
    label 'genomics'

    input:
    path vcf

    output:
    tuple path('pruned.bed'), path('pruned.bim'), path('pruned.fam'), emit: pruned_bed
    path 'plink.prune.in', emit: prune_list
    path 'plink.prune.out', emit: prune_excluded
    path 'raw_variants.bim', emit: raw_bim
    path 'raw_variants.fam', emit: raw_fam

    script:
    def (window, step, r2) = params.indep_pairwise.trim().split(/\s+/)
    """
    plink --vcf ${vcf} --make-bed --out raw_variants --allow-extra-chr --double-id

    plink --bfile raw_variants --indep-pairwise ${window} ${step} ${r2} --out plink --allow-extra-chr

    plink --bfile raw_variants --extract plink.prune.in --make-bed --out pruned --allow-extra-chr
    """
}
