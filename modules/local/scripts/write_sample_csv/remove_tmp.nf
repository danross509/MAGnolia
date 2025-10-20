#!/usr/bin/env nextflow

process REMOVE_TMP_CSV {
    label 'process_single'

    container ""
    conda "conda-forge::pandas=2.3.2"

    //publishDir "${launchDir}", mode: 'move'

    input:
        val reads_count 

    output:
        path "samples_csv_overwrite.txt"    , emit: samples_csv_overwritten, optional: true
        path "error_no_file.txt"            , emit: error_no_file, optional: true
        path "error_file_empty.txt"         , emit: error_file_empty, optional: true
        
 
    script:

    """
    if [[ ! -f ${launchDir}/samples_tmp/samples_tmp.csv ]]; then
        touch error_no_file.txt
    elif [[ ! -s ${launchDir}/samples_tmp/samples_tmp.csv ]]; then
        touch error_file_empty.txt
        mv ${launchDir}/samples_tmp/samples_tmp.csv ${launchDir}/samples.csv
        rm -r ${launchDir}/samples_tmp/
    else
        touch samples_csv_overwrite.txt
        cd ${launchDir}/samples_tmp
        reorder_samples.py
        mv samples.csv ${launchDir}/samples.csv
        rm -r ${launchDir}/samples_tmp/
    fi
    
    """

}