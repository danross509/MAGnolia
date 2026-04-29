#!/usr/bin/env nextflow

process LRBINNER_CREATE_BIN_FILES {
    tag "${meta.id}-${meta.assembler}"
    label 'process_single'

    container ""
    conda ""

    input:
        tuple val(meta), path(contigs), path(binFile)

    output:
        tuple val(meta), path("bins/*.fa"),             emit: bins

    script:
    def prefix = task.ext.prefix ?: "${meta.id}-${meta.assembler}-LRBinner"

    """
    if [[ ! -d bins ]]; then
        mkdir bins
    fi

    lrbinner_extract_bins.py --contigs $contigs --bins $binFile --prefix $prefix
    """
}