#!/usr/bin/env nextflow

include { DASTOOL_FASTATOCONTIG2BIN } from '../../../modules/nf-core/dastool/fastatocontig2bin/main.nf'
include { DASTOOL_DASTOOL } from '../../../modules/nf-core/dastool/dastool/main.nf'
include { HMBIN_VERIFY_F2C2B } from '../../../modules/local/scripts/hifiasm/hmbin_verify_f2c2b.nf'


workflow BIN_REFINEMENT {
    
    take:
        contigs
        initial_bins
    
    main:
    // Generate DASTool contig2bin tsv files
    DASTOOL_FASTATOCONTIG2BIN ( initial_bins, "fa" )

    // If hmBin is used, modify circular MAGs in fastatocontig2bin to match expected format
    f2c2b = DASTOOL_FASTATOCONTIG2BIN.out.fastatocontig2bin.branch { meta, _tsv ->
        hmbin: meta['binner'] == 'hmBin'
        remaining: meta['binner'] != 'hmBin'
    }
    
    hmbin_verified_f2c2b = channel.empty()

    if ( !params.skip_binning && !params.skip_hmbin ) {
        HMBIN_VERIFY_F2C2B ( f2c2b.hmbin )

        hmbin_verified_f2c2b = hmbin_verified_f2c2b.mix ( HMBIN_VERIFY_F2C2B.out )
    }

    // Remove original binner info, group by bin_group
    fastatocontig2bin_grouped = channel.empty()
    fastatocontig2bin_grouped = fastatocontig2bin_grouped.mix ( f2c2b.remaining, hmbin_verified_f2c2b )
        .map { meta, fastatocontig2bin ->
            def meta_new = meta - meta.subMap('binner')
            [ meta_new, fastatocontig2bin ]
        }
        .groupTuple( by: 0 )


    // Join with bin_group contigs
    dastool_input = contigs.join ( fastatocontig2bin_grouped, by: 0 )
 
    // Run DASTool
    DASTOOL_DASTOOL ( dastool_input, [], [] )

    refined_bins = DASTOOL_DASTOOL.out.bins
        .map { meta, bins ->
            def meta_new = meta + [binner: 'DASTool', refined: true]
            [meta_new, bins]
        }

    refined_unbinned = DASTOOL_DASTOOL.out.unbins
        .map { meta, bins ->
            def meta_new = meta + [binner: 'DASTool', refined: true]
            [meta_new, bins]
        }

    emit:

    refined_bins
    refined_unbinned

}