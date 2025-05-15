#!/usr/bin/env nextflow

process WRITE_CONFIG {
    tag ""

    container ""
    conda ""

    //publishDir "${launchDir}", mode: 'move'

    input:
        val short_reads_count
        val short_reads_paired
        val nanopore_barcodes_count
        val pacbio_reads_count
        val reads_corrected

    output:
        //path "empty_file.txt"

    script:

    /*def sampleID = meta.id ? "$meta.id": ""
    def sequencer = meta.sequencer ? "$meta.sequencer" : ""
    def paired_end = meta.paired_end ? true : false
    def corrected = meta.corrected ? true : false
    def bin_group = meta.bin_group ? "$meta.bin_group" : ""
    def assembly_group = meta.assembly_group ? "$meta.assembly_group" : ""
    def reads_1 = reads ? "${reads[0]}" : ""
    def reads_2 = reads.size() == 2 ? "${reads[1]}" : "NA"
    */

    // Checks if config file already exists before writing
    """
    #if [[ -f ${launchDir}/nextflow.config ]]; then
    #    echo "Error: nextflow.config already exists"
    #    exit 1
    #else
        cd ${launchDir}
        write_config.py -s $short_reads_count -p $short_reads_paired -n $nanopore_barcodes_count -pb $pacbio_reads_count -c $reads_corrected
    #fi
    """

}