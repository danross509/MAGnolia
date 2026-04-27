#!/usr/bin/env nextflow

process EGGNOG_DB_DOWNLOAD {
    tag ""
    label 'process_medium'

    container ""
    conda "${moduleDir}/metagenome/environment.yml"

    output:
        path "eggnog_db" , emit: directory


    script:
    def db_dir = "./eggnog_db"
    def args = task.ext.args ?: ''

    """
    download_eggnog_data.py \
    --data_dir $db_dir \
    $args
    """
}