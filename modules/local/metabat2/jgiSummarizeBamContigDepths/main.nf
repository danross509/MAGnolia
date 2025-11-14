#!/usr/bin/env nextflow

process JGISUMMARIZEBAMCONTIGDEPTHS {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/metabat2:15c68d548f9e9b8f"
    conda "bioconda::metabat2=2.15"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.txt"), emit: depth

    script:
    """
    jgi_summarize_bam_contig_depths \
        --outputDepth ${meta.id}.txt \
        $bam

    #bgzip --threads $task.cpus ${meta.id}.txt
    """
}