#!/usr/bin/env nextflow

process CREATE_TMP_CSV {
    label 'process_single'

    container ""
    conda ""

    //publishDir "${launchDir}", mode: 'move'

    //input:
        //tuple val(meta), val(reads), val(count)

    output:
        path "empty_file.txt"

    script:

    // Creates empty file in nextflow work directory to "collect", ensuring all samples are processed before removing tmp folder
    """
    touch empty_file.txt

    if [[ -d ${launchDir}/samples_tmp ]]; then
        echo "Warning : removing existing samples_tmp folder"
        rm -r ${launchDir}/samples_tmp
    fi

    mkdir ${launchDir}/samples_tmp

    if [[ -f ${launchDir}/samples.csv ]]; then
        echo "Warning : existing samples.csv file will be overwritten"
    fi

    cd ${launchDir}/samples_tmp
    create_samples_csv.py

    """

}