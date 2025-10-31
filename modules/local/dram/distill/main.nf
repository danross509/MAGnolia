#!/usr/bin/env nextflow

process DRAM_DISTILL {
    //tag "${meta.id}"
    label 'process_medium'

    container ""
    conda "${moduleDir}/../setup/environment.yml"

    publishDir "${launchDir}/BIN_ANNOTATION/${meta.id}/", mode: 'symlink'

    input:
        tuple val(meta), path(annotations)
        tuple val(meta), path(trnas)
        tuple val(meta), path(rrnas)

    output:
        tuple val(meta), path("DRAM-distill/genome_stats.tsv")          , emit: genome_stats
        tuple val(meta), path("DRAM-distill/metabolism_summary.xlsx")   , emit: metabolism, optional: true
        tuple val(meta), path("DRAM-distill/product.tsv")               , emit: product_tsv, optional: true
        tuple val(meta), path("DRAM-distill/product.html")              , emit: product_html, optional: true
        tuple val(meta), path("DRAM-distill/distill.log")              , emit: log
        //tuple val(meta), versions, emit: versions


    script:

    """
    DRAM.py distill \
        -i $annotations \
        -o DRAM-distill \
        --trna_path $trnas \
        --rrna_path $rrnas

    """
}