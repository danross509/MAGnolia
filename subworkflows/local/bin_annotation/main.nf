#!/usr/bin/env nextflow

include { DRAM_ANNOTATE } from '../../../modules/local/dram/annotate/main.nf'
include { DRAM_DISTILL } from '../../../modules/local/dram/distill/main.nf'

include { BAKTA_BAKTA } from '../../../modules/nf-core/bakta/bakta/main.nf'


workflow BIN_ANNOTATION {
    
    take:
        bins
        bakta_db
    
    main:

    if ( !params.skip_dram ) {
        DRAM_ANNOTATE ( 
            bins
        )

        DRAM_DISTILL (
            DRAM_ANNOTATE.out.annotations,
            DRAM_ANNOTATE.out.trnas,
            DRAM_ANNOTATE.out.rrnas
        )
    }

    if ( !params.skip_bakta ) {
        bakta_bins_input = bins.transpose()
            .map { meta, bin ->
                def sampleID = bin.getBaseName()
                def meta_new = meta + [id: "${sampleID}"]
                [ meta_new, bin ]
            }

        BAKTA_BAKTA (
            bakta_bins_input,
            bakta_db,
            [],[],[],[]
        )
    }

    //emit:

}