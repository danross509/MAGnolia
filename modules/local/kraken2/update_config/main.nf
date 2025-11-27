#!/usr/bin/env nextflow

process KRAKEN2_UPDATE_CONFIG {
    label 'process_single'

    container ""
    conda "${moduleDir}/../classify/environment.yml"

    input:
        val db_dir          // $PATH/to/downloaded/database/dir

    script:
    def config_file = params.database_config

    """
    installation=\$(which kraken2)

    kraken2_update_config.py -d $db_dir -c $config_file -e \$installation
    """

}