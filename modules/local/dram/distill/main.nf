#!/usr/bin/env nextflow

process DRAM_DISTILL {
    tag "${meta.id}"
    label 'process_medium'

    container ""
    conda "${moduleDir}/../setup/environment.yml"

    input:
        tuple val(meta), path(annotations), path(trnas), path(rrnas)

    output:
        tuple val(meta), path("distill/*genome_stats.tsv")          , emit: genome_stats
        tuple val(meta), path("distill/*metabolism_summary.xlsx")   , emit: metabolism, optional: true
        tuple val(meta), path("distill/*product.tsv")               , emit: product_tsv, optional: true
        tuple val(meta), path("distill/*product.html")              , emit: product_html, optional: true
        tuple val(meta), path("distill/*distill.log")               , emit: log
        //tuple val(meta), versions, emit: versions


    script:
    def args = task.ext.args ?: ''

    """
    python \$CONDA_PREFIX/bin/DRAM.py distill \
        -i $annotations \
        -o distill \
        --trna_path $trnas \
        --rrna_path $rrnas \
        $args

    cd distill
    for file in *; do
        if [ -f "\$file" ]; then
            mv \$file ${meta.id}-\${file}
        fi
    done
    """
}