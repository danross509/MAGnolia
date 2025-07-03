#!/usr/bin/env nextflow

process WRITE_SAMPLES_CSV {
    tag "$meta.id"
    label 'process_single'

    container ""
    conda ""

    //publishDir "${launchDir}", mode: 'move'

    input:
        tuple val(meta), val(reads), val(count)

    output:
        path "empty_file.txt"

    script:

    def sampleID = meta.id ? "$meta.id": ""
    def sequencer = meta.sequencer ? "$meta.sequencer" : ""
    def paired_end = meta.paired_end ? true : false
    def corrected = meta.corrected ? true : false
    def bin_group = meta.bin_group ? "$meta.bin_group" : ""
    def assembly_group = meta.assembly_group ? "$meta.assembly_group" : ""
    def reads_1 = reads ? "${reads[0]}" : ""
    def reads_2 = reads.size() == 2 ? "${reads[1]}" : "NA"


    // Creates empty file in nextflow work directory to "collect", ensuring all samples are processed before removing tmp folder
    """
    touch empty_file.txt

    if [[ ! -d ${launchDir}/tmp ]]; then
        mkdir ${launchDir}/tmp
    fi

    if [[ -f ${launchDir}/samples.csv ]]; then
        echo "Warning : existing samples.csv file will be overwritten"
    fi

    cd ${launchDir}/tmp
    write_samples_csv.py $sampleID -s $sequencer -p $paired_end -c $corrected -b $bin_group -a $assembly_group -1 $reads_1 -2 $reads_2 -# $count

    """

}