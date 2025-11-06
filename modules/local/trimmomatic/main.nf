#!/usr/bin/env nextflow

process TRIMMOMATIC {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/trimmomatic:0.39--a688969e471089d7"
    conda "bioconda::trimmomatic=0.39"

    publishDir "${params.resultsDir}/QC/${meta.id}/trimmomatic", mode: 'symlink'

    input:
        tuple val(meta), path(reads_fastq)

    output:
        tuple val(meta), path("${meta.id}_trimmed_{1,2}.fastq.gz")     , emit: reads_PE, optional: true
        tuple val(meta), path("${meta.id}_trimmed_se_{1,2}.fastq.gz")  , emit: reads_SE, optional: true

    script:
    def fq_1_paired = "${meta.id}_trimmed_1.fastq.gz"
    def fq_1_unpaired = "${meta.id}_trimmed_se_1.fastq.gz"
    def fq_2_paired = "${meta.id}_trimmed_2.fastq.gz"
    def fq_2_unpaired = "${meta.id}_trimmed_se_2.fastq.gz"

    if (meta.paired_end) {

        """
        trimmomatic PE  \
        ${reads_fastq[0]} ${reads_fastq[1]} \
        $fq_1_paired \
        $fq_1_unpaired \
        $fq_2_paired \
        $fq_2_unpaired \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
    } else if (!meta.paired_end) {

        """
        trimmomatic SE  \
        ${reads_fastq[0]} \
        $fq_1_unpaired \
        ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
        """
    }
}