#!/usr/bin/env nextflow

process REMOVE_TMP_CSV {
    label 'process_single'

    container ""
    conda "conda-forge::pandas=2.3.2"

    input:
        val launchDirectory
        val reads_count 

    output:
        path "samples_csv_overwrite.txt"    , emit: samples_csv_overwritten, optional: true
        path "error_no_file.txt"            , emit: error_no_file, optional: true
        path "error_file_empty.txt"         , emit: error_file_empty, optional: true
        
 
    script:

    """
    if [[ ! -f ${launchDirectory}/samples_tmp/samples_tmp.csv ]]; then
        touch error_no_file.txt
    elif [[ ! -s ${launchDirectory}/samples_tmp/samples_tmp.csv ]]; then
        touch error_file_empty.txt
        mv ${launchDirectory}/samples_tmp/samples_tmp.csv ${launchDirectory}/samples.csv
        rm -r ${launchDirectory}/samples_tmp/
    else
        touch samples_csv_overwrite.txt
        cd ${launchDirectory}/samples_tmp
        reorder_samples.py
        mv samples.csv ${launchDirectory}/samples.csv
        rm -r ${launchDirectory}/samples_tmp/
    fi
    
    """

}