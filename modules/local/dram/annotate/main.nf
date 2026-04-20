#!/usr/bin/env nextflow

process DRAM_ANNOTATE {
    tag "${meta.id}"
    label 'process_high'

    container ""
    conda "${moduleDir}/../setup/environment.yml"

    input:
        tuple val(meta), path(bins, stageAs: 'bins_input/*')
        val runName

    output:
        tuple val(meta), path("annotate/*annotations.tsv")    , emit: annotations
        tuple val(meta), path("annotate/*trnas.tsv")          , emit: trnas, optional: true
        tuple val(meta), path("annotate/*rrnas.tsv")          , emit: rrnas, optional: true
        tuple val(meta), path("annotate/*genes.faa")          , emit: faa, optional: true
        tuple val(meta), path("annotate/*genes.fna")          , emit: fna, optional: true
        tuple val(meta), path("annotate/*genes.gff")          , emit: gff3, optional: true
        tuple val(meta), path("annotate/*scaffolds.fna")      , emit: scaffolds, optional: true
        tuple val(meta), path("annotate/*genbank")            , emit: genbank, optional: true
        tuple val(meta), path("annotate/*annotate.log")       , emit: log
        //tuple val(meta), versions, emit: versions

    script:
    def args = task.ext.args ?: ''

    """
    python \$CONDA_PREFIX/bin/DRAM.py annotate \
        -i 'bins_input/*' \
        -o annotate \
        $args

    cd annotate
    for file in *; do
        if [ -f "\$file" ]; then
            mv \$file ${meta.id}-\${file}
        fi
    done
    """
}