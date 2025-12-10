#!/usr/bin/env nextflow

include { DRAM_ANNOTATE } from '../../../modules/local/dram/annotate/main.nf'
include { DRAM_DISTILL } from '../../../modules/local/dram/distill/main.nf'

include { BAKTA_BAKTA } from '../../../modules/nf-core/bakta/bakta/main.nf'
include { BAKTA_COLLECT_ANNOTATION_STATS } from '../../../modules/local/bakta/collect_annotation_stats/main.nf'


workflow BIN_ANNOTATION {
    
    take:
        bins
        bakta_db
    
    main:

    if ( !params.skip_dram ) {
        DRAM_ANNOTATE ( 
            bins
        )

        distill_input = DRAM_ANNOTATE.out.annotations
            .join ( DRAM_ANNOTATE.out.trnas, by: 0 )
            .join ( DRAM_ANNOTATE.out.rrnas, by: 0 )

        DRAM_DISTILL (
            distill_input
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

        bakta_summarize = bins.transpose()        
            .map { meta, bin ->
                def sampleID = bin.getBaseName()
                def meta_new = meta + [id: "${sampleID}"]
                [ meta_new, meta ]
            }
            .join ( BAKTA_BAKTA.out.json )
            .map { _meta_bin, meta_all, json ->
                [ meta_all, json ]
            }
            .groupTuple().view()
            

        BAKTA_COLLECT_ANNOTATION_STATS (
            bakta_summarize
        )
    }

    //emit:

}