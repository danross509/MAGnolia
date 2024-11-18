#!/usr/bin/env nextflow

process filter_contaminants {

    container "community.wave.seqera.io/library/bowtie2:2.5.4--d51920539234bea7"
    conda "bioconda::bowtie2=2.5.4"

    publishDir "${launchDir}/CLEAN_READS/", mode: 'symlink'

    input:
        tuple val(meta), path(reads_trimmed)
        path index_db
        val sensitivity

    output:
        tuple val(meta), path("${meta.id}_{1,2}.fastq.gz")

    script:

    index = "${index_db[0].getBaseName(2)}"

    if (meta.paired_end) {

        """
        bowtie2 -p 8 -x $index \
        -1 ${reads_trimmed[0]} \
        -2 ${reads_trimmed[2]} \
        --${sensitivity} \
        --un-conc-gz \
        $meta.id \
        > ${meta.id}_mapped_and_unmapped.sam

        mv "${meta.id}.1" "${meta.id}_1.fastq.gz"
        mv "${meta.id}.2" "${meta.id}_2.fastq.gz"
        """
    } else if (!meta.paired_end) {

        """
        bowtie2 -p 8 -x $index \
        -U ${reads_trimmed[0]} \
        --${sensitivity} \
        --un-conc-gz \
        $meta.id \
        > ${meta.id}_mapped_and_unmapped.sam

        mv "${meta.id}.1" "${meta.id}_1.fastq.gz"
        """
    }
}