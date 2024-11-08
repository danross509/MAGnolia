#!/usr/bin/env nextflow

process filter_contaminants {

    container "community.wave.seqera.io/library/bowtie2:2.5.4--d51920539234bea7"
    conda "bioconda::bowtie2=2.5.4"

    publishDir "${params.projectDir}/CLEAN_READS/", mode: 'symlink'

    input:
        tuple val(sampleID), path(reads_trimmed)
        path index_db
        val sensitivity

    output:
        tuple val(sampleID), path("${sampleID}_{1,2}.fastq.gz")

    script:

    index = "${index_db[0].getBaseName(2)}"

    """
    bowtie2 -p 8 -x $index \
    -1 ${reads_trimmed[0]} \
    -2 ${reads_trimmed[1]} \
    --${sensitivity} \
    --un-conc-gz \
    $sampleID \
    > ${sampleID}_mapped_and_unmapped.sam

    mv "${sampleID}.1" "${sampleID}_1.fastq.gz"
    mv "${sampleID}.2" "${sampleID}_2.fastq.gz"
    """
}