#!/usr/bin/env nextflow

include { CONCATENATE_CLEAN_READS as CONCATENATE_SHORT_READS_R1 } from '../../../modules/local/scripts/concatenate_clean_reads/main.nf'
include { CONCATENATE_CLEAN_READS as CONCATENATE_SHORT_READS_R2 } from '../../../modules/local/scripts/concatenate_clean_reads/main.nf'
include { RENAME_CLEAN_READS as RENAME_SHORT_READS_R1 } from '../../../modules/local/scripts/rename_clean_reads/main.nf'
include { RENAME_CLEAN_READS as RENAME_SHORT_READS_R2 } from '../../../modules/local/scripts/rename_clean_reads/main.nf'

workflow ASSEMBLY_PREP_SHORT {
    take:
    corrected_reads

    main:
    clean_reads_1 = channel.empty()
    clean_reads_2 = channel.empty()

    concatenated_reads = channel.empty()
    original_clean_reads = channel.empty()

    clean_reads_1 = clean_reads_1.mix ( corrected_reads )
        .map { meta, reads ->
            if ( meta.id == meta.assembly_group ) {
                [ meta, reads[0] ]
            } else {
                def meta_new = meta + [ id: meta.assembly_group, bin_group: 'NA' ]
                [ meta_new, reads[0] ]
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
    
    clean_reads_1_coassembly = clean_reads_1
        .filter { _meta, reads -> 
            reads.size() > 1
        }

    clean_reads_1_solo = clean_reads_1
        .filter { _meta, reads ->
            reads.size() == 1
        }

    CONCATENATE_SHORT_READS_R1 (
        clean_reads_1_coassembly,
        "_1"
    )

    RENAME_SHORT_READS_R1 (
        clean_reads_1_solo,
        "_1"
    )

    if ( params.paired_short_reads ) {
        clean_reads_2 = clean_reads_2.mix ( corrected_reads )
            .map { meta, reads ->
                if ( reads.size() == 2 ) {
                    if ( meta.id == meta.assembly_group ) {
                        [ meta, reads[1] ]
                    } else {
                        def meta_new = meta + [ id: meta.assembly_group, bin_group: 'NA' ]
                        [ meta_new, reads[1] ]
                    }
                } else {
                    exit 1, "ERROR: Short reads unpaired : see parameter 'paired_short_reads = ${params.paired_short_reads} in nextflow.config"
                }
            }
            .groupTuple()
            .map { meta, reads ->  //This might cause errors with unpaired reads
                if ( reads.size() >= 2 ) {
                    def meta_new = meta + [ coassembly: true ]
                    [ meta_new, reads ]
                } else {
                    def meta_new = meta + [ coassembly: false ]
                    [ meta_new, reads ]
                }
            }

        clean_reads_2_coassembly = clean_reads_2
        .filter { _meta, reads -> 
            reads.size() > 1
        }

        clean_reads_2_solo = clean_reads_2
        .filter { _meta, reads ->
            reads.size() == 1
        }

        CONCATENATE_SHORT_READS_R2 (
            clean_reads_2_coassembly,
            "_2"
        )

        RENAME_SHORT_READS_R2 (
            clean_reads_2_solo,
            "_2"
        )

        concatenated_reads = concatenated_reads
            .mix ( CONCATENATE_SHORT_READS_R1.out.fastq.join ( CONCATENATE_SHORT_READS_R2.out.fastq ))
            .mix ( RENAME_SHORT_READS_R1.out.fastq.join ( RENAME_SHORT_READS_R2.out.fastq ))
            .map { meta, R1, R2 ->
                [ meta, [R1, R2] ]
            }

        original_clean_reads = original_clean_reads.mix ( corrected_reads )
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

    } else {
        concatenated_reads = concatenated_reads.mix ( CONCATENATE_SHORT_READS_R1.out.fastq, RENAME_SHORT_READS_R1.out.fastq )

        original_clean_reads = original_clean_reads.mix ( corrected_reads )
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
    }

    emit:
    concatenated_reads
    original_clean_reads
}