#!/usr/bin/env nextflow

process WRITE_SAMPLES_CSV {
    tag "$meta.id"
    label 'process_single'

    container ""
    conda ""

    input:
        val launchDirectory
        tuple val(meta), val(reads)
        val placeholder
        

    output:
        path "empty_file.txt"

    script:

    def sampleID = meta.id ? "$meta.id": ""
    def sequencer = meta.sequencer ? "$meta.sequencer" : ""
    def paired_end = meta.paired_end ? true : false
    def corrected = meta.corrected ? true : false
    def assembly_group = meta.assembly_group ? "$meta.assembly_group" : ""
    def bin_group = meta.bin_group ? "$meta.bin_group" : ""
    def reads_1 = reads ? "${reads[0]}" : ""
    def reads_2 = reads.size() == 2 ? "${reads[1]}" : "NA"


    // Creates empty file in nextflow work directory to "collect", ensuring all samples are processed before removing tmp folder
    """
    touch empty_file.txt

    cd ${launchDirectory}/samples_tmp
    write_samples_csv.py $sampleID -s $sequencer -p $paired_end -c $corrected -a $assembly_group -b $bin_group -1 $reads_1 -2 $reads_2

    """

}