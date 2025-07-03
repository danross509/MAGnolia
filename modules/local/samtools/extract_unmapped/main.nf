#!/usr/bin/env nextflow

process SAMTOOLS_EXTRACT_UNMAPPED {
    tag "$meta.id"
    label 'process_low'

    container "community.wave.seqera.io/library/samtools:1.21--0d76da7c3cf7751c"
    conda "bioconda::samtools=1.21"

    publishDir "${launchDir}/CLEAN_READS/nanopore", mode: 'symlink'

    input:

    tuple val(meta), path(alignment)


    output:
    tuple val(meta), path("*.fastq.gz")       , emit: filtered
    //path "versions.yml"                     , emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    samtools fastq -f 4 \
    -c 6 \
    -@ $task.cpus \
    $alignment > ${meta.id}.fastq 

    gzip ${meta.id}*.fastq



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version)
    END_VERSIONS
    """
}