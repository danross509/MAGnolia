#!/usr/bin/env nextflow

process TAXCONVERTER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container ""

    input:
    tuple val(meta), path(krak2)

    output:
    tuple val(meta), path("*.tsv")          , emit: converted
    path "versions.yml"                     , emit: versions


    script:

    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''

    """
    if [[ ! -d \$CONDA_PREFIX/lib/python3.11/site-packages/data/ ]]; then
        mkdir -p \$CONDA_PREFIX/lib/python3.11/site-packages/data/
    fi 
    
    if [[ ! -f \$CONDA_PREFIX/lib/python3.11/site-packages/data/clades.tsv ]]; then
        cp ${moduleDir}/clades.tsv \$CONDA_PREFIX/lib/python3.11/site-packages/data/
    fi

    taxconverter_k2_pre-conversion.py $krak2 ${prefix}_pre-conversion.tsv

    taxconverter kraken2 \\
        -i ${prefix}_pre-conversion.tsv \\
        -o ${prefix}_k2Converted.tsv \\
        $args

    rm ${prefix}_pre-conversion.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vamb: \$(vamb --version | sed 's/Vamb //')
    END_VERSIONS
    """
}