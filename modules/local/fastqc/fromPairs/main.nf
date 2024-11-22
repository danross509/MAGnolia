#!/usr/bin/env nextflow

/*
*/

process fastqc {

    container "community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29"
    conda "bioconda::fastqc=0.12.1"

    publishDir "${launchDir}/QC/${step}/${meta.id}", mode: 'symlink'

    input:
        tuple val(meta), path(reads_fastq)
        val step

    output:
        path "*{1,2}_fastqc.{html,zip}"
        //path "${sampleID}_{1,2}_fastqc.{html,zip}"

    script:
    """
    fastqc -q $reads_fastq
    """
}