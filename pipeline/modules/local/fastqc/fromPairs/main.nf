#!/usr/bin/env nextflow

/*
*/

process fastqc {

    container "community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29"
    conda "bioconda::fastqc=0.12.1"

    publishDir "${params.projectDir}/QC/${step}/${sampleID}", mode: 'symlink'

    input:
        tuple val(sampleID), path(reads_fastq)
        val step

    output:
        path "${sampleID}_{1,2}_fastqc.{html,zip}"

    """
    fastqc -q $reads_fastq
    """
}