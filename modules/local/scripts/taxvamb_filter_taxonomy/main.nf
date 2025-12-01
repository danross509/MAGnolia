#!/usr/bin/env nextflow

process VAMB_FILTER_TAXONOMY {
    tag "$meta.id"
    label 'process_single'

    conda ""
    container ""

    input:
    tuple val(meta), path(assembly), path(abundance_tsv), path(bams), path(taxonomy)
    val minLength

    output:
    tuple val(meta), path(assembly), path(abundance_tsv), path(bams), path("${meta.id}_filteredTax.tsv")

    script:
    """
    taxvamb_filter_taxonomy.py -f $assembly -t $taxonomy -l $minLength -o ${meta.id}_filteredTax.tsv
    """
}