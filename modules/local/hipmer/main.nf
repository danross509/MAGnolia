#!/usr/bin/env nextflow

process HIPMER {
    tag "$meta.id"

    container 
    conda 

    publishDir "${launchDir}/Assembly/${meta.id}/", mode: 'symlink'

    input:
        tuple val(meta), path(reads)
        val depth_file

    output:
        tuple val(meta), path("final_assembly.fasta"), emit: final_contigs
        tuple val(meta), path("final_assembly_depths.txt"), emit: depth

    script:

    threads = bigThreads

    if (meta.paired_end) {

        //reads_1 = reads[0]
        //reads_2 = reads[1]
        //{--presets meta} not working
        //{--tmp-dir tmp} not working
        //gzip files after? pigz?

        """
        mhm2.py \
        -1 ${reads[0]} -2 ${reads[1]} \
        --post-asm-abd $depth_file
        """
    } else if (!meta.paired_end) {

        """
        megahit \
        -r ${reads[0]} \
        --post-asm-abd $depth_file
        """
    }
}