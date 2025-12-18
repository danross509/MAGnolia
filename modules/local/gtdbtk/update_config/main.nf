#!/usr/bin/env nextflow

process GTDB_UPDATE_CONFIG {
    label 'process_single'

    container ""
    conda ""

    input:
        path db             // $downloaded/database/dir
        val download_dir    // $desired/db/directory

    output:
        val "${download_dir}/gtdbtk_${db}"     , emit: db
        path "config_updated.txt"       , emit: verified, optional: true

    script:
    def config_file = params.database_config

    """
    if [[ ! -d $download_dir ]]; then
        mkdir -p $download_dir
    fi

    if [[ -d ${download_dir}/gtdbtk_${db} ]]; then
        echo "ERROR: existing gtdbtk database directory found at ${download_dir}/gtdbtk_${db}"
        exit 1
    else
        target=\$(readlink -e $db)
        mv \$target ${download_dir}/gtdbtk_${db}
    fi

    db_dir=${download_dir}/gtdbtk_${db}

    gtdb_update_config.py -d \${db_dir} -c $config_file
    touch config_updated_gtdb.txt
    """

}