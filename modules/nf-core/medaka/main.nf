process MEDAKA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:1.4.4--py38h130def0_0' :
        'biocontainers/medaka:1.4.4--py38h130def0_0' }"

    publishDir "${launchDir}/ASSEMBLY/${meta.id}/Medaka", mode: 'symlink'

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*_consensus.fa")     , emit: assembly
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Modification to address [E::fai_build3_core] Cannot index files compressed with gzip, please use bgzip
    // 2025-09-22, DANRoss
    def gunzip_assembly = assembly.getExtension() == "gz" ? "gunzip -c $assembly > ${prefix}_tmp.fa" : "cp $assembly ${prefix}_tmp.fa"

    """
    $gunzip_assembly

    medaka_consensus \\
        -t $task.cpus \\
        $args \\
        -i $reads \\
        -d ${prefix}_tmp.fa \\
        -o ./

    mv consensus.fasta ${prefix}_consensus.fa

    #gzip ${prefix}_consensus.fa
    rm ${prefix}_tmp.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}