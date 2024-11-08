#!/usr/bin/env nextflow

process trimmomatic {

    container "community.wave.seqera.io/library/trimmomatic:0.39--a688969e471089d7"
    conda "bioconda::trimmomatic=0.39"

    //publishDir 'QC/trimmomatic', mode: 'symlink'

    input:
        tuple val(sampleID), path(reads_fastq)
        val pairing

    output:
        tuple val(sampleID), path("${sampleID}_p_{1,2}.fastq")

    script:
    fq_1_paired = "${sampleID}_p_1.fastq"
    fq_1_unpaired = "${sampleID}_s_1.fastq"
    fq_2_paired = "${sampleID}_p_2.fastq"
    fq_2_unpaired = "${sampleID}_s_2.fastq"

    """
    trimmomatic $pairing  \
    ${reads_fastq[0]} ${reads_fastq[1]} \
    $fq_1_paired \
    $fq_1_unpaired \
    $fq_2_paired \
    $fq_2_unpaired \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \
    LEADING:3 TRAILING:3 MINLEN:36
    """
}