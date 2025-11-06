#!/usr/bin/env nextflow

process FASTQC {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29"
    conda "bioconda::fastqc=0.12.1"

    publishDir "${params.resultsDir}/QC/${meta.id}/${step}", mode: 'symlink'

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