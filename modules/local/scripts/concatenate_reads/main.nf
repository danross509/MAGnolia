#!/usr/bin/env nextflow

process concatenate_reads {

    container ""
    conda ""

    publishDir "${launchDir}/CLEAN_READS", mode: 'symlink'

    input:
        tuple val(meta), path(clean_reads)

    output:
        tuple val(meta), path("${filename}.fastq.gz")

    script:
 
    //allMeta = [:]
    //allMeta.id = "allReads"
    //allMeta.paired_end = meta.paired_end

    //if (meta.is_paired == "PE"){
    //} else if (meta.is_paired == "SE") {
    //    allMeta.paired_end = false
    //}
    //allMeta.sense = meta.sense
    
    if (meta.sense == "R1") {
        filename = "all_reads_1"
    } else if (meta.sense == "R2"){
        filename = "all_reads_2"
    }

    """
    cat $clean_reads > ${filename}.fastq.gz
    """ 

}