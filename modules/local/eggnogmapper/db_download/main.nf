#!/usr/bin/env nextflow

process EGGNOG_DB_DOWNLOAD {
    tag ""
    label 'process_medium'

    container ""
    conda "${moduleDir}/../metagenome/environment.yml"

    output:
        path "eggnog_db" , emit: directory


    script:
    def db_dir = "./eggnog_db"
    def args = task.ext.args ?: ''

    """
    #Replace download_eggnog_data.py with the fixed version in bin
    if [[ -f \$CONDA_PREFIX/bin/download_eggnog_data.py ]]; then
        rm \$CONDA_PREFIX/bin/download_eggnog_data.py
    fi 

    cp ${moduleDir}/download_eggnog_data_fixed.py \$CONDA_PREFIX/bin/download_eggnog_data.py


    if [[ ! -d ./eggnog_db ]]; then
        mkdir ./eggnog_db
    fi

    download_eggnog_data.py \
    --data_dir $db_dir \
    -y \
    $args
    """
}