#!/usr/bin/env nextflow

process BAKTA_UPDATE_CONFIG {
    label 'process_single'

    container ""
    conda ""

    input:
        path db             // $downloaded/database/dir
        val download_dir    // $desired/db/directory

    output:
        val "${download_dir}/${db}"     , emit: db
        path "config_updated.txt"       , emit: verified, optional: true

    script:
    def config_file = params.database_config

    """
    if [[ ! -d $download_dir ]]; then
        mkdir -p $download_dir
    fi

    if [[ -d ${download_dir}/${db} ]]; then
        echo "ERROR: existing Bakta database directory found at ${download_dir}/${db}"
        exit 1
    else
        target=\$(readlink -e $db)
        mv \$target ${download_dir}/
    fi

    db_dir=${download_dir}/${db}

    bakta_update_config.py -d \${db_dir} -c $config_file
    touch config_updated.txt
    """

}