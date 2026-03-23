#!/usr/bin/env nextflow

process BOWTIE2_BUILD_INDEX {
    tag "${meta.assembler}-${meta.id}"

    container "community.wave.seqera.io/library/bowtie2:2.5.4--d51920539234bea7"
    conda "bioconda::bowtie2=2.5.4"

    input:
        tuple val(meta), path(assembly)

    output:
        //tuple val(meta), path(assembly), path("${meta.id}*")        , emit: assembly_index
        tuple val(meta), path("${meta.id}*")                        , emit: index
        //path "versions.yml"                                         , emit: versions

    script:

    def args = task.ext.args ?: ''

    """
    bowtie2-build \
    --threads $task.cpus \
    $assembly \
    "${meta.id}_index" \
    $args

    """

}