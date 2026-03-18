process XLOCI {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' : 
        'ghcr.io/alejandrogzi/xloci:latest' }"

    input:
    tuple val(meta),  path(genome)
    tuple val(meta1), path(regions)

    output:
    tuple val(meta), path("${meta.id}/*.fa"),       optional: true, emit: fasta
    tuple val(meta), path("${meta.id}/*.fa.gz"),    optional: true, emit: fasta_gz
    tuple val(meta), path("${meta.id}/tmp/*.fa"),   optional: true, emit: fasta_chunks
    tuple val(meta), path("${meta.id}/tmp/*.bed"),  optional: true, emit: bed_chunks
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    xloci \\
        $args \\
        -s $genome \\
        -r $regions \\
        -p ${prefix} \\
        -o ${prefix}
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xloci: \$(xloci --version | sed -e "s/xloci v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}
    touch ${prefix}/*.fa
    touch ${prefix}/*.fa.gz
    touch ${prefix}/tmp/*.fa
    touch ${prefix}/tmp/*.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        xloci: \$(xloci --version | sed -e "s/xloci v//g")
    END_VERSIONS
    """
}
