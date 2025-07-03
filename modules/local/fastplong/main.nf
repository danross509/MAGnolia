#!/usr/bin/env nextflow

process FASTPLONG {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/fastplong:0.2.2--03217e568967f163"
    conda "fastplong=0.2.2"

    publishDir "${launchDir}/QC/${meta.id}/fastplong", mode: 'symlink'

    input:
        tuple val(meta), path(reads_fastq)

    output:
        tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: reads, optional: true
        path "fastplong.html"
        path "fastplong.json"

    script:
    fq_trimmed = "${meta.id}_trimmed.fastq.gz"

        """
        fastplong  \
        -i ${reads_fastq} \
        -o $fq_trimmed
        """

}