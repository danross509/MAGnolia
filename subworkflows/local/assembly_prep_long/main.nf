#!/usr/bin/env nextflow

include { CONCATENATE_CLEAN_READS as CONCATENATE_LONG_READS } from '../../../modules/local/scripts/concatenate_clean_reads/main.nf'

workflow ASSEMBLY_PREP_LONG {
    take:
    corrected_ont_reads
    corrected_pacbio_reads

    main:
    clean_long_reads = Channel.empty()

    concatenated_long_reads = Channel.empty()
    original_clean_long_reads = Channel.empty()

    //concatenated_ont_reads = Channel.empty()
    //original_clean_ont_reads = Channel.empty()

    //concatenated_pacbio_reads = Channel.empty()
    //original_clean_pacbio_reads = Channel.empty()

    /*
    Since meta.sequencer is unique, this could also be done with the long reads mixed
    I don't know how it would affect downstream processes though
    */

    clean_long_reads = clean_long_reads.mix ( corrected_ont_reads, corrected_pacbio_reads )
        .map { meta, reads ->
            if ( meta.id == meta.assembly_group ) {
                [ meta, reads ]
            } else {
                def meta_new = meta + [ id: meta.assembly_group, bin_group: 'NA' ]
                [ meta_new, reads ]
            }
        }
        .groupTuple()
        .map { meta, reads ->
            if ( reads.size() >= 2 ) {
                def meta_new = meta + [ coassembly: true ]
                [ meta_new, reads ]
            } else {
                def meta_new = meta + [ coassembly: false ]
                [ meta_new, reads ]
            }
        }
    
    CONCATENATE_LONG_READS (
        clean_long_reads,
        ""
    )

    concatenated_long_reads = concatenated_long_reads.mix ( CONCATENATE_LONG_READS.out.fastq )

    original_clean_long_reads = original_clean_long_reads.mix ( corrected_ont_reads, corrected_pacbio_reads )
        .map { meta, reads ->
            if ( meta.id == meta.assembly_group ) {
                [ meta, reads ]
            } else {
                def meta_new = meta + [ id: meta.assembly_group, bin_group: 'NA' ]
                [ meta_new, reads ]
            }
        }
        .groupTuple()
        .map { meta, reads ->
            if ( reads.size() >= 2 ) {
                def meta_new = meta + [ coassembly: true ]
                [ meta_new, reads ]
            } else {
                def meta_new = meta + [ coassembly: false ]
                [ meta_new, reads ]
            }
        }


    emit:
    concatenated_long_reads
    original_clean_long_reads
}