#!/usr/bin/env nextflow

include { FASTQC as FASTQC_OUT } from '../../../modules/local/fastqc/fromPairs/main.nf'
include { FASTQC as FASTQC_IN } from '../../../modules/local/fastqc/fromPairs/main.nf'
include { MULTIQC as MULTIQC_IN } from '../../../modules/local/multiqc/main.nf'
include { MULTIQC as MULTIQC_OUT } from '../../../modules/local/multiqc/main.nf'
include { trimmomatic } from '../../../modules/local/trimmomatic/fromPairs/main.nf'
include { FASTP } from '../../../modules/local/fastp/main.nf'
include { MINIMAP2_INDEX as MINIMAP2_PHIX_INDEX} from '../../../modules/local/minimap2/index/main.nf'
include { MINIMAP2_INDEX as MINIMAP2_HOST_INDEX} from '../../../modules/local/minimap2/index/main.nf'
include { MINIMAP2_FILTER_HOST as MINIMAP2_FILTER_PHIX} from '../../../modules/local/minimap2/filter_host/main.nf'
include { MINIMAP2_FILTER_HOST as MINIMAP2_FILTER_HOST} from '../../../modules/local/minimap2/filter_host/main.nf'
include { SAMTOOLS_EXTRACT_UNMAPPED as SAMTOOLS_EXTRACT_PHIX_UNMAPPED} from '../../../modules/local/samtools/extract_unmapped/main.nf'
include { SAMTOOLS_EXTRACT_UNMAPPED as SAMTOOLS_EXTRACT_HOST_UNMAPPED} from '../../../modules/local/samtools/extract_unmapped/main.nf'
include { BOWTIE2_BUILD_INDEX as BOWTIE2_PHIX_INDEX } from '../../../modules/local/bowtie2/build_index/main.nf'
include { BOWTIE2_BUILD_INDEX as BOWTIE2_HOST_INDEX } from '../../../modules/local/bowtie2/build_index/main.nf'
include { FILTER_CONTAMINANTS as BOWTIE2_FILTER_PHIX } from '../../../modules/local/bowtie2/filter_contaminants/main.nf'
include { FILTER_CONTAMINANTS as BOWTIE2_FILTER_HOST } from '../../../modules/local/bowtie2/filter_contaminants/main.nf'
include { CONCATENATE_READS as CONCATENATE_SHORT_READS_R1 } from '../../../modules/local/scripts/concatenate_reads/main.nf'
include { CONCATENATE_READS as CONCATENATE_SHORT_READS_R2 } from '../../../modules/local/scripts/concatenate_reads/main.nf'

workflow QC_SHORT {
    take:
    short_reads     // channel: [ val(meta), path([R1, R2]) ]
    phiX
    host_genome

    main:
    // Create pre-QC fastqc report for input reads;
    // Collect pre-QC fastqc output;
    // Summarize with multiqc
    if ( !params.skip_fastqc ) {
        FASTQC_IN ( short_reads, "Pre-FastQC" )
        all_fastqc_in = FASTQC_IN.out.collect() //.ifEmpty([])
        MULTIQC_IN ( all_fastqc_in, "Pre-FastQC" )
    }

    // Trim adapters and poor quality reads
    //trimmomatic(short_reads)
    trimmed_reads = Channel.empty()
    if ( !params.skip_fastp ){
        FASTP (
            short_reads
        )

        if ( params.paired_short_reads ) {
            trimmed_reads = trimmed_reads.mix ( FASTP.out.reads_PE )
                .map { meta, reads ->
                    def meta_new = meta + [ corrected: true ]
                    [ meta_new, reads ] 
                }
        } else {
            trimmed_reads = trimmed_reads.mix ( FASTP.out.reads_SE )
                .map { meta, reads ->
                    def meta_new = meta + [corrected: true]
                    [ meta_new, reads ] 
                }
        }
    } else {
        trimmed_reads = trimmed_reads.mix ( nanopore_reads )
            .map { meta, reads ->
                if ( params.short_reads_corrected ) {
                    def meta_new = meta + [ corrected: true ]
                    [ meta_new, reads ] 
                } else {
                    def meta_new = meta + [ corrected: false ]
                    [ meta_new, reads ]
                }
            }
    }
    
    // Remove phiX reads
    phiX_filtered_reads = Channel.empty()
    if ( params.remove_phiX ) {
        if ( params.short_read_filter == 'minimap2' ) {
            phiX_index_input = phiX
                .map { meta, genome ->
                    def short_read_preset = ''
                    if ( !params.minimap2_short_read_preset ) {
                        short_read_preset = "sr" // (Default short-read)
                    } else {
                        preset = params.minimap2_short_read_preset
                    }
                    [ meta, genome, short_read_preset ]
                }

            // Build index for host genome
            MINIMAP2_PHIX_INDEX (
                phiX_index_input
            )
    
            phiX_filter_input = trimmed_reads.combine ( MINIMAP2_PHIX_INDEX.out.index )
    
            //Align ONT reads to host genome;
            MINIMAP2_FILTER_PHIX (
                phiX_filter_input
            )

            /*SAMTOOLS_EXTRACT_HOST_UNMAPPED (
                MINIMAP2_HOST_ALIGNMENT.out.alignment_sam
            )*/

            phiX_filtered_reads = phiX_filtered_reads.mix ( MINIMAP2_FILTER_PHIX.out.R1.join ( MINIMAP2_FILTER_PHIX.out.R2, remainder: true ))
                .map { meta, R1, R2 ->
                    def meta_new = meta + [ filtered_phiX: true ]
                    if ( R2 != null ) {
                        [ meta_new, [R1, R2] ]
                    } else if (R2 == null) {
                        [ meta_new, [R1] ]
                    }
                }
        } else if ( params.short_read_filter == 'bowtie2' ) {
            BOWTIE2_PHIX_INDEX (
                phiX
            )

            phiX_filter_input = trimmed_reads.combine ( BOWTIE2_PHIX_INDEX.out )

            BOWTIE2_FILTER_PHIX (
                phiX_filter_input,
                params.bowtie2_sensitivity
            )

            phiX_filtered_reads = phiX_filtered_reads.mix ( BOWTIE2_FILTER_PHIX.out.R1.join ( BOWTIE2_FILTER_PHIX.out.R2, remainder: true ))
                .map { meta, R1, R2 ->
                    def meta_new = meta + [ filtered_phiX: true ]
                    if ( R2 != null ) {
                        [ meta_new, [R1, R2] ]
                    } else if (R2 == null) {
                        [ meta_new, [R1] ]
                    }
                }
        } else {
            exit 1, "ERROR: Short read filter option <${params.short_read_filter}> is invalid"
        }
    } else {
        phiX_filtered_reads = phiX_filtered_reads.mix ( trimmed_reads )
            .map { meta, reads ->
                if ( params.short_reads_corrected ) {
                    def meta_new = meta + [ filtered_phiX: true ]
                    [ meta_new, reads ] 
                } else {
                    def meta_new = meta + [ filtered_phiX: false ]
                    [ meta_new, reads ]
                }
            }
    }

    // Remove host reads
    host_filtered_reads = Channel.empty()
    if ( params.host_genome ) {
        if ( params.short_read_filter == 'minimap2' ) {
            host_index_input = host_genome
                .map { meta, genome ->
                    def short_read_preset = ''
                    if ( !params.minimap2_short_read_preset ) {
                        short_read_preset = "sr" // (Default short-read)
                    } else {
                        preset = params.minimap2_short_read_preset
                    }
                    [ meta, genome, short_read_preset ]
                }

            // Build index for host genome
            MINIMAP2_HOST_INDEX (
                host_index_input
            )
    
            host_filter_input = phiX_filtered_reads.combine ( MINIMAP2_HOST_INDEX.out.index )
    
            //Align ONT reads to host genome;
            MINIMAP2_FILTER_HOST (
                host_filter_input
            )

            /*SAMTOOLS_EXTRACT_HOST_UNMAPPED (
                MINIMAP2_HOST_ALIGNMENT.out.alignment_sam
            )*/

            host_filtered_reads = host_filtered_reads.mix ( MINIMAP2_FILTER_HOST.out.R1.join ( MINIMAP2_FILTER_HOST.out.R2, remainder: true ))
                .map { meta, R1, R2 ->
                    def meta_new = meta + [ filtered_host: true ]
                    if ( R2 != null ) {
                        [ meta_new, [R1, R2] ]
                    } else if (R2 == null) {
                        [ meta_new, [R1] ]
                    }
                }
        } else if ( params.short_read_filter == 'bowtie2' ) {
            BOWTIE2_HOST_INDEX (
                host_genome
            )

            host_filter_input = phiX_filtered_reads.combine ( BOWTIE2_HOST_INDEX.out )

            BOWTIE2_FILTER_HOST (
                host_filter_input,
                params.bowtie2_sensitivity
            )

            host_filtered_reads = host_filtered_reads.mix ( BOWTIE2_FILTER_HOST.out.R1.join ( BOWTIE2_FILTER_HOST.out.R2, remainder: true ))
                .map { meta, R1, R2 ->
                    def meta_new = meta + [ filtered_host: true ]
                    if ( R2 != null ) {
                        [ meta_new, [R1, R2] ]
                    } else if (R2 == null) {
                        [ meta_new, [R1] ]
                    }
                }
        } else {
            exit 1, "ERROR: Short read filter option <${params.short_read_filter}> is invalid"
        }
    } else {
        host_filtered_reads = host_filtered_reads.mix ( phiX_filtered_reads )
            .map { meta, reads ->
                if ( params.short_reads_corrected ) {
                    def meta_new = meta + [ filtered_host: true ]
                    [ meta_new, reads ] 
                } else {
                    def meta_new = meta + [ filtered_host: false ]
                    [ meta_new, reads ]
                }
            }
    }


    // Create post-QC fastqc report for input reads;
    // Collect post-QC fastqc output;
    // Summarize with multiqc
    if ( !params.skip_fastqc ) {
        FASTQC_OUT ( host_filtered_reads, "Post-FastQC" )
        all_fastqc_out = FASTQC_OUT.out.collect()
        MULTIQC_OUT ( all_fastqc_out, "Post-FastQC" )
    }

    concatenated_reads = Channel.empty()
    original_clean_reads = Channel.empty()

    if ( params.assembly_mode == 'coassembly' ) { 
        clean_reads_1 = host_filtered_reads
            .map { meta, reads ->
                def meta_new = meta + [ id: 'allReads' ]
                /*def meta_new = [:]
                meta_new.id = "allReads"
                meta_new.sequencer = 'Illumina'
                if ( reads.size() == 1 ) {
                    meta_new.paired_end = false
                } else if ( reads.size() == 2 ) {
                    meta_new.paired_end = true
                }
                meta_new.corrected = true
                meta_new.assembly_group = "all"*/
                [ meta_new, reads[0] ] 
            }
            .groupTuple()
        
        CONCATENATE_SHORT_READS_R1 (
            clean_reads_1,
            "all_reads_1",
            "CLEAN_READS/illumina"
        )

        if ( params.paired_short_reads ) {
            clean_reads_2 = host_filtered_reads
                .map { meta, reads ->
                    if ( reads.size() == 2 ) {
                        def meta_new = meta + [ id: 'allReads' ]
                        /*def meta_new = [:]
                        meta_new.id = "allReads"
                        meta_new.sequencer = 'Illumina'
                        meta_new.paired_end = true
                        meta_new.corrected = true
                        meta_new.assembly_group = "all"*/
                        [ meta_new, reads[1] ]
                    } else { 
                        exit 1, "ERROR: Paired-read co-assembly contains invalid number of samples"
                     }
                }
                .groupTuple()
            
            CONCATENATE_SHORT_READS_R2 (
                clean_reads_2,
                "all_reads_2",
                "CLEAN_READS/illumina"
            )

            concatenated_reads = concatenated_reads.mix ( CONCATENATE_SHORT_READS_R1.out.join(CONCATENATE_SHORT_READS_R2.out) )
                .map { meta, R1, R2 ->
                    def meta_new = meta + [ assembly_group: "all" ]
                    [ meta_new, [R1, R2] ]
                }

        } else {

            concatenated_reads = concatenated_reads.mix ( CONCATENATE_SHORT_READS_R1.out )
                .map { meta, R1 ->
                    def meta_new = meta + [ assembly_group: "all" ]
                    [ meta_new, R1 ]
                }
        }

        original_clean_reads = original_clean_reads.mix ( host_filtered_reads )
            .map { meta, reads ->
                def meta_new = meta + [ id: 'allReads', assembly_group: "all" ]
                /*def meta_new = [:]
                meta_new.id = "allReads"
                meta_new.sequencer = 'Illumina'
                meta_new.paired_end = meta.paired_end
                meta_new.corrected = true
                meta_new.assembly_group = "all"*/
                [ meta_new, reads ]
            }
            .groupTuple()
        
    } else if ( params.assembly_mode == 'grouped' ) {
        println "To do"
    } else if ( params.assembly_mode == 'per_sample' ) {
        original_clean_reads = original_clean_reads.mix ( host_filtered_reads )
            .map { meta, reads ->
                def meta_new = meta + [ assembly_group: "self" ]
                [ meta_new, reads ]
            }
    } else {
        exit 1, "Assembly mode <${params.assembly_mode}> invalid, exiting"
    }

    // Output view checks
    //println clean_reads_2 == null // always false
    //println clean_reads_2.size() == null // always false
    
    //filter_contaminants.out.R2.count().view()
    //println filter_contaminants.out.size() // always 2
    //clean_reads_1.view()
    /*concatenate_reads_1.out.view()
    if (paired_reads){
        clean_reads_2.view()
        concatenate_reads_2.out.view()
    }*/

    //filter_contaminants.out.view()
    //clean_reads_ch.view()
    //clean_reads_1.map{meta, reads -> println reads.size()}

    //final_reads.view()

    println "QC timestamp"

    emit:
    concatenated_reads
    original_clean_reads

}