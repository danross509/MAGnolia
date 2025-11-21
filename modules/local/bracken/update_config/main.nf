#!/usr/bin/env nextflow

process BRACKEN_UPDATE_CONFIG {
    label 'process_single'

    container ""
    conda ""

    input:
        val bracken_built          // $PATH/to/downloaded/database/dir

    script:
    def config_file = params.database_config

    """
    bracken_update_config.py -d $bracken_built -c $config_file
    """

}