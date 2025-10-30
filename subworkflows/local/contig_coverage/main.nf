#!/usr/bin/env nextflow

include { COVERM_CONTIG } from '../../../modules/nf-core/coverm/contig/main.nf'

/*
Calculate the coverage of each read sample from each MAG

If samples were co-asssembled or co-binned, 
*/

workflow CONTIG_COVERAGE {
    
    take:
        coverage_input
    
    main:

    // Separate reads from contigs
    reads = coverage_input
        .map { meta, contigs, reads ->
            [ meta, reads.flatten() ]
        }

    contigs = coverage_input
        .map { meta, contigs, reads ->
            [ meta, contigs ]
        }

    reads.view()
    contigs.view()


    COVERM_CONTIG ( 
        reads,                      // single channel of fastqs
        contigs,                    // single channel of contigs
        [],                         // val bam_input
        []                          // val interleaved
    )

    //coverage_file = 

    //emit:

    //coverage_file

}