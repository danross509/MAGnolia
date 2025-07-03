#!/usr/bin/env nextflow

process BRACKEN_ABUNDANCE_ESTIMATION {
    tag "${meta.id}"
    label 'process_medium'

    container ""
    conda "${moduleDir}/environment.yml"

    publishDir "${launchDir}/KRAKEN2/${meta.id}/${file_type}", mode: 'symlink'

    input:
        tuple val(meta), path(input_files)
        path database
        val level
        val threshold
        val read_len
        val file_type
        val built

    output:
        tuple val(meta), path("*.bracken")  , emit: reports


    script:
    def name = "${meta.id}_${file_type}"

    """
    bracken \
    -d $database \
    -i $input_files \
    -o ${name}.bracken \
    -r $read_len \
    -l $level \
    -t $threshold
    """
}
