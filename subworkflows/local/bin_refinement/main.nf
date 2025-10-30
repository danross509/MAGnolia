#!/usr/bin/env nextflow

include { DASTOOL_FASTATOCONTIG2BIN } from '../../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_DASTOOL } from '../../../modules/nf-core/dastool/dastool/main.nf'
include { RENAME_REFINED_BINS } from '../../../modules/local/scripts/rename_refined_bins/main.nf'


workflow BIN_REFINEMENT {
    
    take:
        contigs
        initial_bins
    

    main:

    // Generate DASTool contig2bin tsv files
    DASTOOL_FASTATOCONTIG2BIN ( initial_bins, "fa" )


    // Remove original binner info, group by bin_group
    fastatocontig2bin_grouped = Channel.empty()
    fastatocontig2bin_grouped = fastatocontig2bin_grouped.mix ( DASTOOL_FASTATOCONTIG2BIN.out.fastatocontig2bin )
        .map { meta, fastatocontig2bin ->
            def meta_new = meta - meta.subMap('binner')
            [ meta_new, fastatocontig2bin ]
        }
        .groupTuple( by: 0 )


    // Join with bin_group contigs
    dastool_input = contigs.join ( fastatocontig2bin_grouped, by: 0 )

    // Run DASTool
    DASTOOL_DASTOOL ( dastool_input, [], [] )
    //DASTOOL_DASTOOL.out.bins.view()

    RENAME_REFINED_BINS ( DASTOOL_DASTOOL.out.bins )

    //RENAME_REFINED_BINS.out.refined_bins.view()
    //RENAME_REFINED_BINS.out.refined_unbins.view()

    refined_bins = RENAME_REFINED_BINS.out.refined_bins
        .map { meta, bins ->
            def meta_new = meta + [binner: 'DASTool', refined: true]
            [meta_new, bins]
        }

    refined_unbinned = RENAME_REFINED_BINS.out.refined_unbins
        .map { meta, bins ->
            def meta_new = meta + [binner: 'DASTool', refined: true]
            [meta_new, bins]
        }

    emit:

    refined_bins
    refined_unbinned

}