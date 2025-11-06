#!/usr/bin/env nextflow

process BBDuk {
    label 'process_medium'

    container ""
    conda ""

    publishDir "${params.resultsDir}/QC/BBDuk", mode: 'symlink'

    input:
        tuple val(meta), path(reads_fastq)

    output:
        tuple val(meta), path("${meta.id}_{1,2}.{p,s}.fastq")

    script:
    fq_1_paired = "${meta.id}_1.p.fastq"
    fq_1_unpaired = "${meta.id}_1.s.fastq"
    fq_2_paired = "${meta.id}_2.p.fastq"
    fq_2_unpaired = "${meta.id}_2.s.fastq"

    if (meta.paired_end) {

        """
        trimmomatic PE  \
        ${reads_fastq[0]} ${reads_fastq[1]} \
        $fq_1_paired \
        $fq_1_unpaired \
        $fq_2_paired \
        $fq_2_unpaired \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True \
        LEADING:3 TRAILING:3 MINLEN:36
        """
    } else if (!meta.paired_end) {

        """
        trimmomatic SE  \
        ${reads_fastq[0]} \
        $fq_1_unpaired \
        ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 \
        LEADING:3 TRAILING:3 MINLEN:36
        """
    }
}