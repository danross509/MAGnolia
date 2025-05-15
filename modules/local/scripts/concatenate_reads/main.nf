#!/usr/bin/env nextflow

process CONCATENATE_READS {
    tag "$filename"

    container ""
    conda ""

    publishDir "${launchDir}/${save_directory}", mode: 'symlink'

    input:
        tuple val(meta), path(reads)
        val new_filename
        val save_directory


    output:
        tuple val(meta), path("${filename}.fastq.gz")                                   , emit: symlinks, optional: true
        tuple val(meta), val("${launchDir}/${save_directory}${filename}.fastq.gz")      , emit: setup_reads_csv, optional: true

    script:

    filename = "${new_filename}"
    if (filename.size() < 1) {
        filename = "${meta.id}"
    }

    """
    cat $reads > ${filename}.fastq.gz
    """ 

}