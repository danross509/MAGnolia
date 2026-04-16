#!/usr/bin/env nextflow

process RENAME_SOLO_CONTIGS {
    tag "$meta.id"
    label 'process_single'

    container ""
    conda ""

    input:
        tuple val(meta), path(contigs)

    output:
        tuple val(meta), path("${meta.id}_contigs.fa"), emit: concatenated_fasta

    script:
    //def extension = contigs.getExtension() == "gz" ? "fa.gz" : "fa"

    """
    ln -s ${contigs[0]} ${meta.id}_contigs.fa
    """
}