#!/usr/bin/env nextflow

process DRAM_ANNOTATE {
    tag "${meta.id}"
    label 'process_high'

    container ""
    conda "${moduleDir}/../setup/environment.yml"

    input:
        tuple val(meta), path(bins, stageAs: 'bins_input/*')

    output:
        tuple val(meta), path("DRAM-annotation/annotations.tsv")    , emit: annotations
        tuple val(meta), path("DRAM-annotation/trnas.tsv")          , emit: trnas, optional: true
        tuple val(meta), path("DRAM-annotation/rrnas.tsv")          , emit: rrnas, optional: true
        tuple val(meta), path("DRAM-annotation/genes.faa")          , emit: faa, optional: true
        tuple val(meta), path("DRAM-annotation/genes.fna")          , emit: fna, optional: true
        tuple val(meta), path("DRAM-annotation/genes.gff")          , emit: gff3, optional: true
        tuple val(meta), path("DRAM-annotation/scaffolds.fna")      , emit: scaffolds, optional: true
        tuple val(meta), path("DRAM-annotation/genbank")            , emit: genbank, optional: true
        tuple val(meta), path("DRAM-annotation/annotate.log")       , emit: log
        //tuple val(meta), versions, emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    python \$CONDA_PREFIX/bin/DRAM.py annotate \
        -i 'bins_input/*' \
        -o DRAM-annotation \
        $args
    """
}