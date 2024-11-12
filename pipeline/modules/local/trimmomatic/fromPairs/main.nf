#!/usr/bin/env nextflow

process trimmomatic {

    container "community.wave.seqera.io/library/trimmomatic:0.39--a688969e471089d7"
    conda "bioconda::trimmomatic=0.39"

    publishDir "${params.projectDir}/QC/trimmomatic", mode: 'symlink'

    input:
        tuple val(sampleID), path(reads_fastq)
        val pairing

    output:
        tuple val(sampleID), path("${sampleID}_{1,2}.p.fastq")

    script:
    fq_1_paired = "${sampleID}_1.p.fastq"
    fq_1_unpaired = "${sampleID}_1.s.fastq"
    fq_2_paired = "${sampleID}_2.p.fastq"
    fq_2_unpaired = "${sampleID}_2.s.fastq"

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