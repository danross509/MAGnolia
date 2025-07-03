#!/usr/bin/env nextflow

include { KRAKEN2 } from '../../../modules/local/kraken2/classify/main.nf'
include { KRONA_K2_IMPORT_TAXONOMY } from '../../../modules/local/krona/import_taxonomy/main.nf'
include { BRACKEN_ABUNDANCE_ESTIMATION as BRACKEN_READS } from '../../../modules/local/bracken/abundance_estimation/main.nf'
include { BRACKEN_SUMMARIZE_ABUNDANCE as BRACKEN_READS_SUMMARY } from '../../../modules/local/bracken/summarize_abundance/main.nf'
//include { METAPHLAN } from '../../../modules/local/biobakery/metaphlan/main.nf'

workflow READ_CONTIG_TAXONOMY {
    take:
    reads
    kraken_db
    bracken_built
    stage

    main:

    if ( !params.skip_kracken2 ) {
        KRAKEN2 (
            reads,
            kraken_db,
            params.kraken2_useNames,
            params.kraken2_confidence,
            params.kraken2_quick,
            params.kraken2_preload_db,
            params.mdThreads,
            params.kraken2_minimizer,
            params.kraken2_zero_counts,
            params.kraken2_minimum_hits,
            stage
        )

        if ( !params.skip_krona ) {
            KRONA_K2_IMPORT_TAXONOMY (
                KRAKEN2.out.reports,
                stage
            )
        }

        if ( !params.skip_bracken ) {
            BRACKEN_READS (
                KRAKEN2.out.reports,
                kraken_db,
                params.bracken_classification_level,
                params.bracken_threshold,
                params.bracken_read_length,
                stage,
                bracken_built
            )

            summarize_bracken_input = BRACKEN_READS.out.reports
                .collect { meta, report ->
                    [ report ]
                }.view()

            BRACKEN_READS_SUMMARY (
                summarize_bracken_input,
                stage,
                params.bracken_classification_level,
                params.bracken_summary_NA
            )
        }

    }

    if ( !params.skip_metaphlan4 ) {

    }

    //emit:
    
    
}