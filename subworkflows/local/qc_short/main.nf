#!/usr/bin/env nextflow

include { FASTQC as FASTQC_OUT } from '../../../modules/local/fastqc/fromPairs/main.nf'
include { FASTQC as FASTQC_IN } from '../../../modules/local/fastqc/fromPairs/main.nf'
include { MULTIQC as MULTIQC_IN } from '../../../modules/local/multiqc/main.nf'
include { MULTIQC as MULTIQC_OUT } from '../../../modules/local/multiqc/main.nf'
include { TRIMMOMATIC } from '../../../modules/local/trimmomatic/main.nf'
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

workflow QC_SHORT {
    take:
    short_reads     // channel: [ val(meta), path([R1, R2]) ]
    phiX_index
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

    //short_reads.view()

    // Trim adapters and poor quality reads
    //trimmomatic(short_reads)
    trimmed_reads = Channel.empty()
    if ( !params.skip_short_read_trimming ) {
        if ( params.short_read_trimmer == 'fastp' ) {
            FASTP (
                short_reads,
                params.fastp_adapter_sequences,
                params.fastp_auto_adapter_detection,
                params.fastp_qualified_quality_phred,
                params.fastp_unqualified_percent_limit,
                params.fastp_disable_length_filtering,
                params.fastp_length_required,
                params.fastp_deduplication
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

        } else if ( params.short_read_trimmer == 'trimmomatic' ) {
            TRIMMOMATIC (
                short_reads
            )

            if ( params.paired_short_reads ) {
                trimmed_reads = trimmed_reads.mix ( TRIMMOMATIC.out.reads_PE )
                    .map { meta, reads ->
                        def meta_new = meta + [ corrected: true ]
                        [ meta_new, reads ] 
                    }
            } else {
                trimmed_reads = trimmed_reads.mix ( TRIMMOMATIC.out.reads_SE )
                    .map { meta, reads ->
                        def meta_new = meta + [corrected: true]
                        [ meta_new, reads ] 
                    }
            }
        } else {
            exit 1, "ERROR: Short read trimmer option <${params.short_read_trimmer}> is invalid"
        }
    } else {
        trimmed_reads = trimmed_reads.mix ( short_reads )
            .map { meta, reads ->
                if ( params.short_reads_corrected ) {
                    def meta_new = meta + [ corrected: true ]
                    [ meta_new, reads ] 
                } else {
                    [ meta, reads ]
                }
            }
    }
    
    //trimmed_reads.view()

    if ( params.remove_phiX && !params.host_genome ) {
        suffix_phix = "corrected"
        suffix_host = "NA"
    } else {
        suffix_phix = "filtered_phix"
        suffix_host = "corrected"
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
                phiX_filter_input,
                suffix_phix
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
            /*BOWTIE2_PHIX_INDEX (
                phiX
            )*/

            //phiX_filter_input = trimmed_reads.combine ( BOWTIE2_PHIX_INDEX.out.index )
            phiX_filter_input = trimmed_reads.combine ( phiX_index )

            //phiX_filter_input.view()

            BOWTIE2_FILTER_PHIX (
                phiX_filter_input,
                params.bowtie2_sensitivity,
                suffix_phix
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

    //phiX_filtered_reads.view()

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
                host_filter_input,
                suffix_host
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

            host_filter_input = phiX_filtered_reads.combine ( BOWTIE2_HOST_INDEX.out.index )

            BOWTIE2_FILTER_HOST (
                host_filter_input,
                params.bowtie2_sensitivity,
                suffix_host
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

    //host_filtered_reads.view()

    emit:
    host_filtered_reads

}