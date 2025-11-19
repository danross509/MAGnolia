#!/usr/bin/env nextflow

process KRAKEN2_UPDATE_CONFIG {
    label 'process_single'

    container ""
    conda ""

    input:
        val db_dir          // $PATH/to/downloaded/database/dir

    output:        

    script:
    def config_file = params.database_config

    """
    kraken2_update_config.py -d ${db_dir} -c $config_file
    """

}