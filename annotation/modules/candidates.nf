// Subconjunto codificante / splicing y raro, listo para el modelo de ranking.
process SELECT_CANDIDATES {
    tag "${sample}"
    label 'python'
    publishDir "${params.outdir}/candidates", mode: 'copy'

    input:
    tuple val(sample), path(multianno)

    output:
    tuple val(sample), path("${sample}.candidates.tsv"), emit: tsv
    path "${sample}.candidates.log", emit: log

    script:
    def syn = params.keep_synonymous ? '--keep-synonymous' : ''
    """
    select_candidates.py ${multianno} \\
        -o ${sample}.candidates.tsv \\
        --gene-col '${params.gene_col}' \\
        --func-col '${params.func_col}' \\
        --exonicfunc-col '${params.exonicfunc_col}' \\
        --clinsig-col '${params.clinsig_col}' \\
        --func-keep '${params.func_keep}' \\
        --max-af ${params.max_af} \\
        --max-af-clinvar ${params.max_af_clinvar} \\
        ${syn} 2> ${sample}.candidates.log
    cat ${sample}.candidates.log >&2
    """

    stub:
    """
    touch ${sample}.candidates.tsv ${sample}.candidates.log
    """
}
