#!/usr/bin/env nextflow

process CONCATENATE_RAW_READS {
    tag "$meta.id"

    container ""
    conda ""

    publishDir "${launchDir}/${sample_prep_directory}", mode: 'symlink'

    input:
        tuple val(meta), path(reads)
        val sample_prep_directory
        val strand


    output:
        tuple val(meta), path("${filename}.fastq.gz")                                       , emit: fastq, optional: true
        tuple val(meta), val("${launchDir}/${sample_prep_directory}/${filename}.fastq.gz")  , emit: setup_reads_csv, optional: true

    script:

    def filename = "${meta.id}${strand}"
    def command = reads.size() > 1 ? "cat $reads > ${filename}.fastq.gz" : reads.size() == 1 ? "ln -s ${reads[0]} ${filename}.fastq.gz" : ""

    """
    $command
    """

}