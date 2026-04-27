#!/usr/bin/env nextflow

process EGGNOG_UPDATE_CONFIG {
    label 'process_single'

    container ""
    conda ""

    input:
        path db             // $downloaded/database/dir
        val download_dir    // $desired/db/directory

    output:
        val "${download_dir}"           , emit: directory
        path "config_updated.txt"       , emit: verified, optional: true

    script:
    def config_file = params.database_config

    """
    if [[ -d ${download_dir} ]]; then
        echo "ERROR: existing eggNOG database directory found at ${download_dir}"
        exit 1
    else
        target=\$(readlink -e $db)
        mv \$target ${download_dir}
    fi

    db_dir=${download_dir}

    eggnog_update_config.py -d \${db_dir} -c $config_file
    touch config_updated.txt
    """

}