#!/usr/bin/env nextflow

process GATB_MINIA {
    tag "$meta.id"
    label 'process_high'

    container "community.wave.seqera.io/library/minia:3.2.6--92bae1756baab1ef", "community.wave.seqera.io/library/gatb:1.4.2--cf37d08b7005497e"
    conda "bioconda::minia=3.2.6", "bioconda::gatb=1.4.2"

    publishDir "${launchDir}/Assembly/GATB/${meta.id}/", mode: 'symlink'

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
        ./minia \
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