#!/usr/bin/env nextflow

process CHECKM_UPDATE_CONFIG {
    label 'process_single'

    container ""
    conda ""

    input:
        val db_dir          // $PATH/to/downloaded/database/dir

    script:
    def config_file = params.database_config

    """
    checkm_update_config.py -d ${db_dir} -c $config_file
    """

}