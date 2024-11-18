#!/usr/bin/env nextflow

process build_index {

    container "community.wave.seqera.io/library/bowtie2:2.5.4--d51920539234bea7"
    conda "bioconda::bowtie2=2.5.4"

    //publishDir '', mode: 'symlink'

    input:
        path reference

    output:
        path "${host_DB}*"

    script:

    host_DB = "${reference.getBaseName(1)}"

    """
    bowtie2-build $reference $host_DB

    """
}