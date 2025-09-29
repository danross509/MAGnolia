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
    if [[ ! -f ${launchDir}/tmp/samples_tmp.csv ]]; then
        touch error_no_file.txt
    elif [[ ! -s ${launchDir}/tmp/samples_tmp.csv ]]; then
        touch error_file_empty.txt
        mv ${launchDir}/tmp/samples_tmp.csv ${launchDir}/samples.csv
        rm -r ${launchDir}/tmp/
    else
        touch samples_csv_overwrite.txt
        cd ${launchDir}/tmp
        reorder_samples.py
        mv samples.csv ${launchDir}/samples.csv
        rm -r ${launchDir}/tmp/
    fi
    
    """

}