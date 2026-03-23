#!/usr/bin/env nextflow

process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29"
    conda "bioconda::fastqc=0.12.1"

    //publishDir "${params.resultsDir}/QC/${meta.id}/${step}", mode: 'symlink'

    input:
        tuple val(meta), path(reads_fastq)

    output:
        path "*{1,2}_fastqc.{html,zip}"

    script:
    def args = task.ext.args ?: ''

    """
    fastqc -q $reads_fastq \
    $args
    """
}