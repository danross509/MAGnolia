#!/usr/bin/env nextflow

process RENAME_SOLO_TAXONOMY {
    tag "${meta.id}-${meta.assembler}"
    label 'process_single'

    container ""
    conda ""

    input:
        tuple val(meta), path(taxonomies)

    output:
        tuple val(meta), path("${meta.id}_contigTax.tsv"), emit: concatenated_tax
        //path "versions.yml", emit: versions

    script:

    """
    ln -s ${taxonomies[0]} ${meta.id}_contigTax.tsv
    """
}