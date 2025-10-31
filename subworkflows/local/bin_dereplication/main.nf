#!/usr/bin/env nextflow

include { DREP_COMPARE } from '../../../modules/nf-core/drep/compare/main.nf'
include { DREP_DEREPLICATE } from '../../../modules/nf-core/drep/dereplicate/main.nf'
include { RENAME_DEREPLICATED_BINS } from '../../../modules/local/scripts/rename_dereplicated_bins/main.nf'


workflow BIN_DEREPLICATION {
    
    take:
        refined_bins
    
    main:

    dereplicate_input = refined_bins.transpose()
        .map { meta, bin ->
            def meta_new = [:]
            meta_new.id = 'dereplicated_bins'
            [ meta_new, bin ]
        }
        .groupTuple( by: 0 )

    DREP_COMPARE ( dereplicate_input )

    DREP_DEREPLICATE (
        dereplicate_input,
        DREP_COMPARE.out.directory
    )

    dereplicated_bins = RENAME_DEREPLICATED_BINS ( DREP_DEREPLICATE.out.fastas )

    // How to connect bins back to original meta?

    dereplicated_bins

    emit:

    dereplicated_bins

}