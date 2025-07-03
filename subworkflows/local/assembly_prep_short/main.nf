#!/usr/bin/env nextflow

include { CONCATENATE_CLEAN_READS as CONCATENATE_SHORT_READS_R1 } from '../../../modules/local/scripts/concatenate_clean_reads/main.nf'
include { CONCATENATE_CLEAN_READS as CONCATENATE_SHORT_READS_R2 } from '../../../modules/local/scripts/concatenate_clean_reads/main.nf'

workflow ASSEMBLY_PREP_SHORT {
    take:
    corrected_reads

    main:
    clean_reads_1 = Channel.empty()
    clean_reads_2 = Channel.empty()

    concatenated_reads = Channel.empty()
    original_clean_reads = Channel.empty()

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
    
    CONCATENATE_SHORT_READS_R1 (
        clean_reads_1,
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

        CONCATENATE_SHORT_READS_R2 (
            clean_reads_2,
            "_2"
        )

        concatenated_reads = concatenated_reads.mix ( CONCATENATE_SHORT_READS_R1.out.fastq.join ( CONCATENATE_SHORT_READS_R2.out.fastq ))
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
        concatenated_reads = concatenated_reads.mix ( CONCATENATE_SHORT_READS_R1.out.fastq ).view()
            //.map { meta, R1 ->
            //    [ meta, [R1] ]
            //}

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

    /*
    if ( params.assembly_mode == 'coassembly' ) { 
        clean_reads_1 = corrected_reads
            .map { meta, reads ->
                def meta_new = meta + [ id: 'allReads', bin_group: 'coassembled' ]
                [ meta_new, reads[0] ] 
            }
            .groupTuple()
        
        CONCATENATE_SHORT_READS_R1 (
            clean_reads_1,
            "CLEAN_READS/illumina"
        )

        if ( params.paired_short_reads ) {
            clean_reads_2 = corrected_reads
                .map { meta, reads ->
                    if ( reads.size() == 2 ) {
                        def meta_new = meta + [ id: 'allReads', bin_group: 'coassembled' ]
                        [ meta_new, reads[1] ]
                    } else { 
                        exit 1, "ERROR: Paired-read co-assembly contains invalid number of samples"
                    }
                }
                .groupTuple()
            
            CONCATENATE_SHORT_READS_R2 (
                clean_reads_2,
                "CLEAN_READS/illumina"
            )

            concatenated_reads = concatenated_reads.mix ( CONCATENATE_SHORT_READS_R1.out.symlinks.join ( CONCATENATE_SHORT_READS_R2.out.symlinks ))
                .map { meta, R1, R2 ->
                    [ meta, [R1, R2] ]
                }

        } else {

            concatenated_reads = concatenated_reads.mix ( CONCATENATE_SHORT_READS_R1.out.symlinks )

        }

        original_clean_reads = original_clean_reads.mix ( corrected_reads )
            .map { meta, reads ->
                def meta_new = meta + [ id: 'allReads', assembly_group: "all", bin_group: 'coassembled' ]
                [ meta_new, reads ]
            }
            .groupTuple()
        
    } else if ( params.assembly_mode == 'grouped' ) {
        println "To do"
    } else if ( params.assembly_mode == 'per_sample' ) {
        original_clean_reads = original_clean_reads.mix ( corrected_reads )
    } else {
        exit 1, "Assembly mode <${params.assembly_mode}> invalid, exiting"
    }
    */

    emit:
    concatenated_reads
    original_clean_reads
}