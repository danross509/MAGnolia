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

    // Calculate coverage for each final bin on the original input reads
    original_reads_input = original_reads
        .map {meta, reads ->
            def meta_new = [:]
            meta_new.id = 'original_reads'
            meta_new.paired_end = meta.paired_end
            [ meta_new, reads ]
        }
        .groupTuple()
        .filter { it.size() > 0 }
        .map { meta, reads ->
            [ meta, reads.flatten().sort { it.name } ]
        }


    // If samples have been co-assembled or co-binned, calculate coverage for the concatenated reads
    grouped_reads_input = bin_group_reads
        .map { meta, reads ->
            if ( !meta.coassembly && !meta.cobinning ) {
                []                                          // Do not include reads which have been processed individually
            } else {
                def meta_new = [:]
                meta_new.id = 'grouped_reads'
                meta_new.paired_end = meta.paired_end
                [ meta_new, reads ]                         // Include concatenated read groups
            }
        }
        .groupTuple()
        .filter { it.size() > 0 }
        .map { meta, reads ->
            [ meta, reads.flatten().sort { it.name } ]
        }



    //original_reads_input.view()
    //grouped_reads_input.view()
    //bins.view()


    COVERM_GENOME_ORIGINAL ( 
        original_reads_input,             // single channel of fastqs
        bins,                       // single channel of bins
        [],                         // val bam_input
        [],                         // val interleaved
        'auto'
    )

    COVERM_GENOME_GROUPED ( 
        grouped_reads_input,             // single channel of fastqs
        bins,                       // single channel of bins
        [],                         // val bam_input
        [],                         // val interleaved
        'auto'
    )

    //coverage_file = 

    //emit:

    //coverage_file

}