#!/usr/bin/env nextflow

include { COVERM_GENOME as COVERM_GENOME_ORIGINAL } from '../../../modules/nf-core/coverm/genome/main.nf'
include { COVERM_GENOME as COVERM_GENOME_GROUPED } from '../../../modules/nf-core/coverm/genome/main.nf'

/*
Calculate the coverage of each read sample from each MAG

If samples were co-asssembled or co-binned, 
*/

workflow BIN_COVERAGE {
    
    take:
        bins
        original_reads
        bin_group_reads
    
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
        .filter { _meta, reads -> reads.size() > 0 }
        .map { meta, reads ->
            [ meta, reads.flatten().sort { f -> f.toString().split('/')[-1] } ]
        }

    // If samples have been co-assembled or co-binned, calculate coverage for the concatenated reads
    grouped_reads_input = bin_group_reads
        .filter { meta, _reads -> meta.coassembly || meta.cobinning }
        .map { meta, reads ->
            def meta_new = [:]
            meta_new.id = 'grouped_reads'
            meta_new.paired_end = meta.paired_end
            [ meta_new, reads ]                         // Include concatenated read groups
        }
        .groupTuple()
        .filter { _meta, reads -> reads.size() > 0 }
        .map { meta, reads ->
            [ meta, reads.flatten().sort { f -> f.toString().split('/')[-1] } ]
        }

    COVERM_GENOME_ORIGINAL ( 
        original_reads_input,       // single channel of fastqs
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