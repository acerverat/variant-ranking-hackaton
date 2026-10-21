#!/usr/bin/env nextflow
// Anotación de variantes para priorización clínica (hackatón).
// VCF (WGS, GRCh38) -> NORMALIZE -> ANNOVAR -> SELECT_CANDIDATES
// Basado en Pipelines_Inmegen/Annotation (INMEGEN, CC BY-NC-SA 4.0).

include { NORMALIZE         } from './modules/normalize'
include { ANNOVAR           } from './modules/annovar'
include { SELECT_CANDIDATES } from './modules/candidates'

workflow {
    log.info """
    Anotación de variantes
    ======================
    sample_info : ${params.sample_info}
    humandb     : ${params.humandb}
    protocol    : ${params.annovar_protocol}
    regiones    : ${params.regions ?: 'genoma completo'}
    outdir      : ${params.outdir}
    """.stripIndent()

    // Verificación temprana: cada base de datos 'f' debe existir en humandb
    def protocols  = params.annovar_protocol.tokenize(',')
    def operations = params.annovar_operation.tokenize(',')
    if (protocols.size() != operations.size()) {
        error "annovar_protocol y annovar_operation deben tener el mismo número de elementos"
    }
    if (!workflow.stubRun) {
        def missing = protocols.findAll { p ->
            !file("${params.humandb}/${params.buildver}_${p}.txt").exists()
        }
        if (missing) {
            error "Faltan bases de ANNOVAR en ${params.humandb}: ${missing.join(', ')}"
        }
    }

    vcfs_ch = channel.fromPath(params.sample_info)
        .splitCsv(sep: '\t', header: true)
        .map { row -> tuple(row.Sample, file(row.Path, checkIfExists: true), file(row.Index, checkIfExists: true)) }

    // la referencia es opcional: sin ella se separan multialélicos pero no se alinean indels
    def fasta = params.ref_fasta ?: "${projectDir}/assets/NO_REF"
    ref_fasta = tuple(file(fasta, checkIfExists: true), file("${fasta}.fai", checkIfExists: true))
    regions   = params.regions   ? file(params.regions,   checkIfExists: true) : file("${projectDir}/assets/NO_REGIONS")

    NORMALIZE(vcfs_ch, ref_fasta, regions)
    ANNOVAR(NORMALIZE.out.vcf)
    SELECT_CANDIDATES(ANNOVAR.out.txt)
}
