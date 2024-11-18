#!/usr/bin/env nextflow

process fastp {

    container "community.wave.seqera.io/library/fastp:0.23.4--f8cefc1e5f7a782e"
    conda "bioconda::fastp=0.23.4"

    publishDir "${launchDir}/QC/fastp", mode: 'symlink'

    input:
        tuple val(meta), path(reads_fastq)

    output:
        tuple val(meta), path("${meta.id}_{1,2,se}.fastq.gz")

    script:
    fq_1_paired = "${meta.id}_trimmed_1.fastq.gz"
    fq_2_paired = "${meta.id}_trimmed_2.fastq.gz"
    fq_unpaired = "${meta.id}_trimmed_se.fastq.gz"

    if (meta.paired_end) {

        """
        fastp  \
        -i ${reads_fastq[0]} \
        -I ${reads_fastq[1]} \
        -o $fq_1_paired \
        -O $fq_2_paired \
        --unpaired1 $fq_unpaired 
        """
    } else if (!meta.paired_end) {

        """
        fastp  \
        ${reads_fastq[0]} \
        $fq_unpaired 
        """
    }
}