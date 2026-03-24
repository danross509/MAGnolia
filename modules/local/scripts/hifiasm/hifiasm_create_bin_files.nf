#!/usr/bin/env nextflow

process HIFIASM_CREATE_BIN_FILES {
    tag "$meta.id"
    label 'process_single'

    container ""
    conda ""

    input:
        tuple val(meta), path(all_contigs), path(circle_contigs), path(bins)

    output:
        tuple val(meta), path("circular_MAGs/*.fa"),    emit: circular_mags
        tuple val(meta), path("bins/*.fa"),             emit: bins

    script:
    def prefix = task.ext.prefix ?: "${meta.id}-${meta.assembler}-hmBin"

    """
    if [[ ! -d circular_MAGs ]]; then
        mkdir circular_MAGs
    fi

    hifiasm_split_circular_contigs.py $circle_contigs circular_MAGs $prefix

    if [[ ! -d bins ]]; then
        mkdir bins
    fi

    hifiasm_extract_bins.py --contigs $all_contigs --bins $bins --prefix $prefix
    """
}