#!/usr/bin/env nextflow

process CONCATENATE_CLEAN_READS {
    tag "$meta.id"
    label 'process_low'

    container ""
    conda ""

    //publishDir "${params.resultsDir}/CLEAN_READS/${meta.sequencer}", mode: 'symlink'

    input:
        tuple val(meta), path(reads)
        val strand


    output:
        tuple val(meta), path("${meta.id}${strand}.fastq.gz")   , emit: fastq

    script:

    def filename = "${meta.id}${strand}"
    def command = reads.size() > 1 ? "cat $reads > ${filename}.fastq.gz" : reads.size() == 1 ? "ln -s ${reads[0]} ${filename}.fastq.gz" : ""

    """
    $command
    """

}