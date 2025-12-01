#!/usr/bin/env nextflow

process TAXVAMB_CONCATENATE_TAXONOMY {
    tag "$meta.id"
    label 'process_single'

    container ""
    conda ""

    input:
        tuple val(meta), path(taxonomies)

    output:
        tuple val(meta), path("${meta.id}_contigTax.tsv"), emit: concatenated_tax
        //path "versions.yml", emit: versions

    script:

    if ( meta.cobinning ) {
        """
        concatenate_contig_taxonomies.py \
        -i $taxonomies \
        -o ${meta.id}_contigTax.tsv
        """
    } else {
        """
        ln -s ${taxonomies[0]} ${meta.id}_contigTax.tsv
        """
    }



}