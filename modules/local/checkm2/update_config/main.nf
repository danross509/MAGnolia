#!/usr/bin/env nextflow

process CHECKM2_UPDATE_CONFIG {
    label 'process_single'

    container ""
    conda ""

    input:
        tuple val(meta), path(db)       // $downloaded/database/file
        val download_dir                // $desired/db/directory

    output:
        tuple val(meta), path("${download_dir}/${db}")      , emit: db
        path "config_updated.txt"                           , emit: verified, optional: true

    script:
    def config_file = params.database_config

    """
    if [[ ! -d $download_dir ]]; then
        mkdir -p $download_dir
    fi

    if [[ -f ${download_dir}/${db} ]]; then
        echo "ERROR: existing Checkm2 database found at ${download_dir}/${db}"
        exit 1
    else
        target=\$(readlink -e $db)
        mv \$target ${download_dir}/
    fi

    db_dir=${download_dir}/${db}

    checkm2_update_config.py -d \${db_dir} -c $config_file
    touch config_updated.txt
    """

}