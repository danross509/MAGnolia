#!/usr/bin/env nextflow

process RENAME_CLEAN_READS {
    tag "$meta.id"
    label 'process_low'

    container ""
    conda ""

    input:
        tuple val(meta), path(reads)
        val strand


    output:
        tuple val(meta), path("${meta.id}${strand}.fastq.gz")   , emit: fastq

    script:

    def filename = "${meta.id}${strand}"

    """
    ln -s ${reads[0]} ${filename}.fastq.gz
    """

}