#!/usr/bin/env nextflow

process WRITE_CONFIG {
    label 'process_single'

    container ""
    conda ""

    input:
        val projectDirectory
        val launchDirectory
        val short_reads_count
        val short_reads_paired
        val nanopore_barcodes_count
        val pacbio_reads_count
        val reads_corrected
        val use_gpu

    script:
    // Checks if config file already exists before writing

    def config_file = "${projectDirectory}/nextflow.config"
    def config_folder = "${projectDirectory}/configs"

    """
    if [[ -f ${launchDirectory}/nextflow.config ]]; then
        echo "Warning: existing nextflow.config will be overwritten"
        rm ${launchDirectory}/nextflow.config
    fi
    
    cd ${launchDirectory}
    write_config.py -s $short_reads_count -p $short_reads_paired -n $nanopore_barcodes_count -pb $pacbio_reads_count -c $reads_corrected -f $config_file -g $use_gpu

    if [[ -d ${launchDirectory}/configs ]]; then
        echo "Warning: configs folder already exists, skipping"
    else 
        cp -r $config_folder ./
        rm configs/databases.config 
    fi

    """

}