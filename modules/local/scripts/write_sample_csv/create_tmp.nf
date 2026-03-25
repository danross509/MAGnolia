#!/usr/bin/env nextflow

process CREATE_TMP_CSV {
    label 'process_single'

    container ""
    conda ""

    input:
        val launchDirectory

    output:
        path "empty_file.txt"

    script:

    // Creates empty file in nextflow work directory to "collect", ensuring all samples are processed before removing tmp folder
    """
    touch empty_file.txt

    if [[ -d ${launchDirectory}/samples_tmp ]]; then
        echo "Warning : removing existing samples_tmp folder"
        rm -r ${launchDirectory}/samples_tmp
    fi

    mkdir ${launchDirectory}/samples_tmp

    if [[ -f ${launchDirectory}/samples.csv ]]; then
        echo "Warning : existing samples.csv file will be overwritten"
    fi

    cd ${launchDirectory}/samples_tmp
    create_samples_csv.py

    """

}