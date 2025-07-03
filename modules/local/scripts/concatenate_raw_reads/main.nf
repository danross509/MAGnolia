#!/usr/bin/env nextflow

process CONCATENATE_RAW_READS {
    tag "$meta.id"
    label 'process_low'

    container ""
    conda ""

    publishDir "${launchDir}/${sample_prep_directory}", mode: 'symlink'

    input:
        tuple val(meta), path(reads)
        val sample_prep_directory
        val strand


    output:
        tuple val(meta), path("${meta.id}${strand}.fastq.gz")                                       , emit: fastq, optional: true
        tuple val(meta), val("${launchDir}/${sample_prep_directory}/${meta.id}${strand}.fastq.gz")  , emit: setup_reads_csv, optional: true

    script:

    def command = reads.size() > 1 ? "cat $reads > ${meta.id}${strand}.fastq.gz" : reads.size() == 1 ? "ln -s ${reads[0]} ${meta.id}${strand}.fastq.gz" : ""

    """
    $command
    """

}