// Normaliza y filtra el VCF antes de anotar:
//  - restringe a cromosomas canónicos (params.contigs): la referencia de Ensembl usada
//    (Homo_sapiens.GRCh38.dna.primary_assembly.fa) no incluye los scaffolds sin ubicar ni
//    los decoys de hs38d1 del VCF original, y llama 'MT' al mitocondrial en vez de 'M', lo que
//    hace fallar a `bcftools norm` con "sequence not found" si no se filtran antes. Estos contigs
//    son ~2.3% de las variantes (scaffolds/decoys sin genes anotados + MT, que un llamador diploide
//    estándar no maneja de forma confiable por heteroplasmia; requeriría un pipeline mitocondrial aparte)
//  - separa sitios multialélicos (ANNOVAR y el ranking trabajan por alelo)
//  - alinea a la izquierda los indels si se da la referencia
//  - elimina alelos '*' (deleciones que abarcan el sitio)
//  - conserva FILTER=PASS y genotipos con calidad/profundidad mínimas
//  - opcionalmente restringe además a regiones (BED), p. ej. CDS ± 50 pb
process NORMALIZE {
    tag "${sample}"
    label 'bcftools'
    publishDir "${params.outdir}/normalized", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi)
    tuple path(ref_fasta), path(ref_fai)
    path regions

    output:
    tuple val(sample), path("${sample}.norm.vcf.gz"), path("${sample}.norm.vcf.gz.tbi"), emit: vcf
    path "${sample}.norm.stats.txt", emit: stats

    script:
    def ref_opt     = ref_fasta.name != 'NO_REF'     ? "-f ${ref_fasta} --check-ref w" : ''
    def regions_opt = regions.name   != 'NO_REGIONS' ? "-R ${regions}"                 : ''
    def pass_opt    = params.pass_only ? "-f PASS" : ''
    def contigs_opt = params.contigs ? "-t ${params.contigs}" : ''
    """
    bcftools view ${pass_opt} ${contigs_opt} ${regions_opt} -Ou ${vcf} \\
    | bcftools norm -m -any ${ref_opt} -Ou \\
    | bcftools view -e 'ALT="*" || GT="ref" || GT="mis" || FMT/GQ<${params.min_gq} || FMT/DP<${params.min_dp}' \\
        --threads ${task.cpus} -Oz -o ${sample}.norm.vcf.gz
    tabix -p vcf ${sample}.norm.vcf.gz

    bcftools stats ${sample}.norm.vcf.gz > ${sample}.norm.stats.txt
    """

    stub:
    """
    touch ${sample}.norm.vcf.gz ${sample}.norm.vcf.gz.tbi ${sample}.norm.stats.txt
    """
}
