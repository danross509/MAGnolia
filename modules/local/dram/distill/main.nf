#!/usr/bin/env nextflow

process DRAM_DISTILL {
    //tag "${meta.id}"
    label 'process_medium'

    container ""
    conda "${moduleDir}/environment.yml"

    //publishDir "${launchDir}/KRAKEN2/${file_type}/${meta.id}", mode: 'symlink'

    input:


    output:


    script:

    """
    DRAM.py distill \\

    """
}