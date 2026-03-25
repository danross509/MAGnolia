#!/usr/bin/env nextflow

process CONFIRM_UNIQUE_FILES {
    label 'process_single'

    container ""
    conda ""

    input:
        tuple val(meta), path(file)
        val file_directory
        val new_extension

    output:        
        tuple val(meta), path("output/*"), optional: true

    script:
    def original_extension = file.getExtension()

    """
    mkdir output

    if [[ ! -f ${file_directory}/${meta.id}.${new_extension} ]]; then
        ln -s \$(readlink -e $file) output/${meta.id}.${original_extension}
    #else
        #count=1
        #while [[ -f ${file_directory}/${meta.id}_\$count.${new_extension} ]]; do
        #    let count++
        #done
        #ln -s \$(readlink -e $file) output/${meta.id}_\$count.${original_extension}
        #if [[ $original_extension == "bam" && -f ${file_directory}/${meta.id}.bam.pbi ]]; then
        #    mv ${file_directory}/${meta.id}.bam.pbi ${file_directory}/${meta.id}_\$count.bam.pbi
        #fi
    fi
    """

}