#!/usr/bin/env nextflow

process JGI_DEPTHS4MAXBIN {
    tag "${meta.assembler}-${meta.id}"

    conda "conda-forge::python=3.11"
    container ""

    input:
    tuple val(meta), path(fasta), path(depth)

    output:
    tuple val(meta), path(fasta), val([]), path("*.abund"), emit: output

    script:
    """
    jgi_depths4maxbin.py ${depth}
    """
}