#!/usr/bin/env nextflow

process DRAM_UPDATE_CONFIG {
    label 'process_single'

    container ""
    conda "${moduleDir}/../setup/environment.yml"

    input:
        val db_dir          // $PATH/to/database/dir
        path placeholder    // ensures DRAM setup is finished before generating dram_config.txt

    output:        

    script:
    def config_file = "${launchDir}/configs/database_download.config"

    """
    DRAM-setup.py export_config > ${db_dir}/dram_config.txt

    dram_update_config.py -d ${db_dir}/dram_config.txt -c $config_file
    """

}