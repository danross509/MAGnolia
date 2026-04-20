process BAKTA_COLLECT_ANNOTATION_STATS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bakta:1.11.4--pyhdfd78af_0'
        : 'biocontainers/bakta:1.11.4--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(jsons)

    output:
    tuple val(meta), path("*.tsv")
    //path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
   bakta_collect-annotation-stats.py \\
        ${jsons} \\
        ${args} \\
        --prefix ${prefix} \\
        --output .
    """

}