#!/usr/bin/env nextflow

include { DREP_COMPARE } from '../../../modules/nf-core/drep/compare/main.nf'
include { DREP_DEREPLICATE } from '../../../modules/nf-core/drep/dereplicate/main.nf'
include { RENAME_DEREPLICATED_BINS } from '../../../modules/local/scripts/rename_dereplicated_bins/main.nf'
include { MODIFY_CHECKM_FOR_DREP } from '../../../modules/local/scripts/modify_checkm_for_drep/main.nf'


workflow BIN_DEREPLICATION {
    
    take:
        refined_bins
        bin_evaluations
    
    main:

    dereplicate_input = refined_bins.transpose()
        .map { _meta, bin ->
            def meta_new = [:]
            meta_new.id = 'dereplicated_bins'
            [ meta_new, bin ]
        }
        .groupTuple( by: 0 )

    DREP_COMPARE ( dereplicate_input )

    // If input exists, modify it to drep expectations
    // Base it on checkm vs checkm2 usage
    // If empty, []

    bin_evaluation_input = channel.empty()

    if ( bin_evaluations ) {
        MODIFY_CHECKM_FOR_DREP ( bin_evaluations )

        bin_evaluation_input = bin_evaluation_input.mix ( MODIFY_CHECKM_FOR_DREP.out.modified_bin_quality )
    }

    DREP_DEREPLICATE (
        dereplicate_input,
        DREP_COMPARE.out.directory,
        bin_evaluation_input
    )

    dereplicated_bins = RENAME_DEREPLICATED_BINS ( DREP_DEREPLICATE.out.fastas )

    // How to connect bins back to original meta?

    emit:

    dereplicated_bins

}