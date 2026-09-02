process RANK_HITS {
    tag "Rank predicted candidates"
    label 'cheminformatics'

    input:
    path affinity_matrix
    path descriptors

    output:
    path 'top_predicted_candidates.tsv', emit: top_hits
    path 'hit_summary.png', emit: plot

    script:
    """
    python3 ${projectDir}/bin/rank_hits.py \
      --affinity_matrix ${affinity_matrix} \
      --descriptors ${descriptors} \
      --top_n ${params.top_n} \
      --out_table top_predicted_candidates.tsv \
      --out_plot hit_summary.png
    """
}
