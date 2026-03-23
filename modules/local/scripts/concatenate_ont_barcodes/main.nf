#!/usr/bin/env nextflow

process CONCATENATE_ONT_BARCODES {
    tag "$meta.id"
    label 'process_low'

    container ""
    conda ""

    publishDir "${file_directory}", mode: 'symlink'

    input:
        tuple val(meta), path(reads)
        val file_directory


    output:
        tuple val(meta), path("${meta.id}*.fastq.gz")                                       , emit: fastq, optional: true
        //tuple val(meta), val("${file_directory}/${meta.id}.fastq.gz")  , emit: setup_reads_csv, optional: true

    script:
    def reads_size = reads.size()

    // Will crash if single fastq.gz in barcode folder AND filename already == ${meta.id}.fastq.gz
    """
    if [[ ! -f ${file_directory}/${meta.id}.fastq.gz ]]; then
        file_name=${meta.id}.fastq.gz
    else
        count=1
        while [[ -f ${file_directory}/${meta.id}_\$count.fastq.gz ]]; do
            let count++
        done

        file_name=${meta.id}_\$count.fastq.gz
    fi

    if [[ $reads_size > 1 ]]; then
        cat $reads > \$file_name
    elif [[ $reads_size == 1 ]]; then
        ln -s ${reads[0]} \$file_name
    else 
        ""
    fi
    """

}