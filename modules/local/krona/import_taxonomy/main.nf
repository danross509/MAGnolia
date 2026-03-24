#!/usr/bin/env nextflow

process KRONA_K2_IMPORT_TAXONOMY {
    tag "$meta.id"
    label 'process_low'

    container ""
    conda "${moduleDir}/environment.yml"

    input:
        tuple val(meta), path(report)
        val file_type
        val placeholder

    output:
        tuple val(meta), path("*_krona.html")                        , emit: kronagram, optional: true


    script:
    def args = task.ext.args ?: ''

    """
    ktImportTaxonomy \
    -m 3 -t 5 \
    $report \
    $args \
    -o ${meta.id}_${file_type}_krona.html
    """
}