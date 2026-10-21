// Anotación con ANNOVAR (table_annovar.pl).
// Adaptado del módulo de INMEGEN (Pipelines_Inmegen/modules/annotation/annovar.nf,
// licencia CC BY-NC-SA 4.0): protocolos parametrizables, -hgvs, umbral de splicing
// configurable, -nastring '.' y sin la edición del encabezado con sed.
//
// Dos pasos en vez de --vcfinput directo: table_annovar.pl con --vcfinput regenera además un
// VCF anotado, y esa regeneración interna tiene un bug conocido ("prefield not defined ... with
// field=N and prefield=M") que revienta con ciertos sitios intergénicos/ncRNA cuando hay muchas
// bases con columnas variables (nos pasó con el genoma completo, protocolo de 9 bases). Como solo
// se usa el .txt en SELECT_CANDIDATES, se convierte primero a avinput (convert2annovar.pl,
// conservando el VCF original como columnas Otherinfo) y se anota con --otherinfo en vez de
// --vcfinput, evitando por completo el paso de regeneración de VCF que falla.
process ANNOVAR {
    tag "${sample}"
    label 'annovar'
    publishDir "${params.outdir}/annovar", mode: 'copy'

    input:
    tuple val(sample), path(vcf), path(tbi)

    output:
    tuple val(sample), path("${sample}.${params.buildver}_multianno.txt"), emit: txt

    script:
    def protocols  = params.annovar_protocol.tokenize(',')
    def operations = params.annovar_operation.tokenize(',')
    // Opciones por protocolo: HGVS y umbral de splicing solo para anotación génica (g)
    def args = operations.collect { op ->
        op == 'g' ? "-hgvs -splicing_threshold ${params.splicing_threshold}" : ''
    }.join(',')
    """
    ${params.annovar_dir}/convert2annovar.pl -format vcf4 -includeinfo \\
        ${vcf} --outfile ${sample}.avinput

    ${params.annovar_dir}/table_annovar.pl ${sample}.avinput ${params.humandb} \\
        --buildver ${params.buildver} \\
        --out ${sample} \\
        --remove \\
        --protocol ${protocols.join(',')} \\
        --operation ${operations.join(',')} \\
        --argument '${args}' \\
        --nastring . \\
        --polish \\
        --otherinfo \\
        --thread ${task.cpus} \\
        --maxgenethread ${task.cpus}
    """

    stub:
    """
    touch ${sample}.${params.buildver}_multianno.txt
    """
}
