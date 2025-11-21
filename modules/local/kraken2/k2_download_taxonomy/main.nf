#!/usr/bin/env nextflow

process K2_DOWNLOAD_TAXONOMY {
    tag ""
    label 'process_medium'

    container ""
    conda "${moduleDir}/../classify/environment.yml"

    //publishDir "${params.resultsDir}/KRAKEN2/${meta.id}/${file_type}", mode: 'symlink'

    input:
        val db_dir

    output:
        val "${db_dir}/kraken2_db"                        , emit: directory


    script:
    def name = "${db_dir}/kraken2_db"

    """
    k2 download-taxonomy \
    --db $name
    """
}