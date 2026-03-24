#!/usr/bin/env nextflow

include { KRAKEN2 } from '../../../modules/local/kraken2/classify/main.nf'
include { KRONA_K2_UPDATE_TAXONOMY } from '../../../modules/local/krona/update_taxonomy/main.nf'
include { KRONA_K2_IMPORT_TAXONOMY } from '../../../modules/local/krona/import_taxonomy/main.nf'
include { BRACKEN_ABUNDANCE_ESTIMATION as BRACKEN } from '../../../modules/local/bracken/abundance_estimation/main.nf'
include { BRACKEN_SUMMARIZE_ABUNDANCE as BRACKEN_SUMMARY } from '../../../modules/local/bracken/summarize_abundance/main.nf'
include { COMBINE_BRACKEN_OUTPUTS } from '../../../modules/local/bracken/combine_bracken_outputs/main.nf'
include { TAXCONVERTER } from '../../../modules/local/taxconverter/main.nf'

workflow READ_CONTIG_TAXONOMY {
    take:
    reads
    kraken_db
    bracken_built
    stage

    main:
    tax_4_vamb = channel.empty()

    if ( !params.kraken2_memory_mapping && ( params.kraken2_hashk2d_size + 5.GB ) <= params.maxMem ) {
        kraken_mm = false
    } else {
        kraken_mm = true
    }

    KRAKEN2 (
        reads,
        kraken_db,
        params.kraken2_useNames,
        params.kraken2_confidence,
        params.kraken2_quick,
        kraken_mm,
        params.kraken2_minimizer,
        params.kraken2_zero_counts,
        params.kraken2_minimum_hits,
        stage
    )

    if ( !params.skip_krona ) {
        KRONA_K2_UPDATE_TAXONOMY ()

        KRONA_K2_IMPORT_TAXONOMY (
            KRAKEN2.out.reports,
            stage,
            KRONA_K2_UPDATE_TAXONOMY.out.placeholder
        )
    }

    if ( !params.skip_bracken ) {
        BRACKEN (
            KRAKEN2.out.reports,
            kraken_db,
            params.bracken_classification_level,
            params.bracken_threshold,
            params.bracken_read_length,
            stage,
            bracken_built
        )

        summarize_bracken_input = BRACKEN.out.output
            .collect { _meta, report ->
                [ report ]
            }

        BRACKEN_SUMMARY (
            summarize_bracken_input,
            stage,
            params.bracken_classification_level,
            params.bracken_summary_NA
        )

        COMBINE_BRACKEN_OUTPUTS (
            summarize_bracken_input,
            stage
        )
    }

    if ( stage == "contigs" ) {
        TAXCONVERTER ( KRAKEN2.out.output )

        tax_4_vamb = tax_4_vamb.mix ( TAXCONVERTER.out.converted )
    } 

    emit:
    tax_4_vamb
    
}