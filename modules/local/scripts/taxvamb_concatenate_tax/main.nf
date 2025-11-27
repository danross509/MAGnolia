#!/usr/bin/env nextflow

process TAXVAMB_CONCATENATE_TAXONOMY {
    tag "$meta.id"
    label 'process_low'

    container ""
    conda ""

    input:
        tuple val(meta), path(taxonomies)

    output:
        tuple val(meta), path("${meta.id}_.tsv"), emit: concatenated_tax
        //path "versions.yml", emit: versions

    script:

    if ( meta.cobinning ) {
        """
        taxvamb_concatenate_tax.py \
        --input $taxonomies \
        --output ${meta.id}_k2Concatenated.tsv
        """
    } else {
        """
        ln -s ${taxonomies[0]} ${meta.id}_k2Concatenated.tsv
        """
    }



}