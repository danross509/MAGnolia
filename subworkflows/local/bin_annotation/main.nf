#!/usr/bin/env nextflow

include { DRAM_ANNOTATE } from '../../../modules/local/dram/annotate/main.nf'
include { DRAM_DISTILL } from '../../../modules/local/dram/distill/main.nf'


workflow BIN_ANNOTATION {
    
    take:
        bins
    
    main:

    DRAM_ANNOTATE ( 
        bins
    )

    DRAM_DISTILL (
        DRAM_ANNOTATE.out.annotations,
        DRAM_ANNOTATE.out.trnas,
        DRAM_ANNOTATE.out.rrnas
    )


    //emit:

}