#!/usr/bin/env nextflow

include { COVERM_GENOME as COVERM_GENOME_ORIGINAL } from '../../../modules/nf-core/coverm/genome/main.nf'
include { COVERM_GENOME as COVERM_GENOME_GROUPED } from '../../../modules/nf-core/coverm/genome/main.nf'

/*
Calculate the coverage of each read sample from each MAG

If samples were co-asssembled or co-binned, 
*/

workflow BIN_COVERAGE {
    
    take:
        original_reads
        bin_group_reads
        bins
    
    main:

    grouped_reads_input = bin_group_reads
        .map { meta, reads ->
            if ( !meta.coassembly && !meta.cobinning ) {
                []
            } else {
                def meta_new = meta + [id: "${meta.id}_grouped"]
                [ meta_new, reads ]
            }
        }

    original_reads.view()
    grouped_reads_input.view()
    bins.view()

    COVERM_GENOME_ORIGINAL ( 
        original_reads,             // single channel of fastqs
        bins,                       // single channel of bins
        [],                         // val bam_input
        [],                         // val interleaved
        'auto'
    )

    //coverage_file = 

    //emit:

    //coverage_file

}