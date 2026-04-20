#!/usr/bin/env nextflow

process BRACKEN_ABUNDANCE_ESTIMATION {
    tag "${meta.id}-${meta.assembler}"
    label 'process_medium'

    container ""
    conda "${moduleDir}/environment.yml"

    input:
        tuple val(meta), path(input_files)
        path database
        val level
        val threshold
        val read_len
        val file_type
        val built

    output:
        tuple val(meta), path("*.bracken")  , emit: output


    script:
    def name = "${meta.id}_${file_type}"
    def args = task.ext.args ?: ''

    """
    bracken \
    -d $database \
    -i $input_files \
    -o ${name}.bracken \
    -r $read_len \
    -l $level \
    -t $threshold \
    $args
    """
}
