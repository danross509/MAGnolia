#!/usr/bin/env nextflow

include { TOULLIGQC as TOULLIGQC_IN } from '../../../modules/nf-core/toulligqc/main.nf'
include { TOULLIGQC as TOULLIGQC_OUT } from '../../../modules/nf-core/toulligqc/main.nf'
include { MULTIQC as MULTIQC_IN } from '../../../modules/local/multiqc/main.nf'
include { MULTIQC as MULTIQC_OUT } from '../../../modules/local/multiqc/main.nf'
include { FASTPLONG } from '../../../modules/local/fastplong/main.nf'
include { MINIMAP2_INDEX } from '../../../modules/local/minimap2/index/main.nf'
include { MINIMAP2_FILTER_HOST } from '../../../modules/local/minimap2/filter_host/main.nf'
include { SAMTOOLS_EXTRACT_UNMAPPED } from '../../../modules/local/samtools/extract_unmapped/main.nf'
include { CONCATENATE_READS as CONCATENATE_ONT_READS } from '../../../modules/local/scripts/concatenate_reads/main.nf'

workflow QC_NANOPORE {
    take:
    nanopore_reads
    host_genome

    main:
    // Create pre-QC ToulligQC report for input reads;
    // Collect pre-QC ToulligQC output;
    // Summarize with multiqc
    all_toulligqc_in = Channel.empty()
    if (!params.skip_toulligqc){
        TOULLIGQC_IN(
            nanopore_reads, 
            "Pre-ToulligQC")
        /*all_toulligqc_in = TOULLIGQC_IN.out.report_html
            .map { meta, report ->
                return [ report ]
            }
            .collect()
        MULTIQC_IN(all_toulligqc_in, "pre-toulligQC")*/
    }

    // Trim reads with FastpLong
    trimmed_nanopore_reads = Channel.empty()
    if ( !params.skip_fastplong ) {
        FASTPLONG (
            nanopore_reads
        )

        trimmed_nanopore_reads = trimmed_nanopore_reads.mix ( FASTPLONG.out.reads )
            .map { meta, reads ->
                def meta_new = meta + [ corrected: true ]
                [ meta_new, reads ] 
            }
    } else {
        trimmed_nanopore_reads = trimmed_nanopore_reads.mix ( nanopore_reads )
            .map { meta, reads ->
                if ( params.nanopore_hq || params.nanopore_reads_corrected ) {
                    def meta_new = meta + [ corrected: true ]
                    [ meta_new, reads ] 
                } else {
                    def meta_new = meta + [ corrected: false ]
                    [ meta_new, reads ]
                }
            }
    }

    // Remove host reads with Minimap2 + Samtools
    filtered_nanopore_reads = Channel.empty()
    if ( params.host_genome ) {
        minimap2_index_input = host_genome
            .map { meta, genome ->
                def nanopore_preset = ''
                if ( !params.minimap2_nanopore_preset ) {
                    if ( !params.skip_fastplong || params.nanopore_hq || params.nanopore_reads_corrected ) {
                        nanopore_preset = "lr:hq" // (<1% error)
                    } else {
                        nanopore_preset = "map-ont" // (<10% error)
                    }
                } else {
                    preset = params.minimap2_nanopore_preset
                }
                [ meta, genome, nanopore_preset ]
            }

        // Build index for host genome
        MINIMAP2_INDEX (
            minimap2_index_input
        )
    
        minimap2_filter_input = trimmed_nanopore_reads.combine ( MINIMAP2_INDEX.out.index )
    
        //Align ONT reads to host genome;
        MINIMAP2_FILTER_HOST (
            minimap2_filter_input
        )

        /*SAMTOOLS_EXTRACT_UNMAPPED (
            MINIMAP2_HOST_ALIGNMENT.out.alignment_sam
        )*/

        filtered_nanopore_reads = filtered_nanopore_reads.mix ( MINIMAP2_FILTER_HOST.out.filtered )
            .map { meta, reads ->
                def meta_new = meta + [filtered: true]
                [ meta_new, reads ] 
            }
    } else {
        filtered_nanopore_reads = filtered_nanopore_reads.mix ( trimmed_nanopore_reads )
            .map { meta, reads ->
                if ( params.nanopore_reads_corrected ) {
                    def meta_new = meta + [filtered: true]
                    [ meta_new, reads ] 
                } else {
                    def meta_new = meta + [filtered: false]
                    [ meta_new, reads ]
                }
            }
    }



    // Create post-QC toulligQC report for input reads;
    // Collect post-QC touliggQC output;
    // Summarize with multiqc
    all_toulligqc_out = Channel.empty()
    if (!params.skip_toulligqc){
        TOULLIGQC_OUT(
            filtered_nanopore_reads, 
            "Post-ToulligQC")
        /*all_toulligqc_out = TOULLIGQC_OUT.out.report_html
            .map { meta, report ->
                return [ report ]
            }
            .collect()
        MULTIQC_OUT(all_toulligqc_out, "post-toulligQC")*/
    }

    concatenated_reads = Channel.empty()
    original_clean_reads = Channel.empty()

    if (params.assembly_mode == 'coassembly') { 
        clean_nanopore_reads = filtered_nanopore_reads
            .map { meta, reads ->
                def meta_new = [:]
                //println reads.size()
                meta_new.id = "allReadsONT"
                meta_new.sequencer = "ONT"
                meta_new.corrected = meta.corrected
                meta_new.filtered = meta.filtered
                return [ meta_new, reads ]
            }
            .groupTuple()

        CONCATENATE_ONT_READS(
            clean_nanopore_reads,
            "all_reads_ont",
            "CLEAN_READS/nanopore"
        )

        concatenated_reads = concatenated_reads.mix ( CONCATENATE_ONT_READS.out )
            .map { meta, reads ->
                def meta_new = meta + [ assembly_group: "all" ]
                [ meta_new, reads ]
            }
        original_clean_reads = original_clean_reads.mix ( clean_nanopore_reads )
            .map { meta, reads ->
                def meta_new = meta + [ assembly_group: "all" ]
                [ meta_new, reads ]
            }

    } else if ( params.assembly_mode == 'grouped' ) {
        println "To do"
    } else if ( params.assembly_mode == 'per_sample' ) {
        original_clean_reads = filtered_nanopore_reads
            .map { meta, reads ->
                def meta_new = meta + [ assembly_group: "self" ]
                [ meta_new, reads ]
            }
            
    } else {
        println "Assembly mode <${params.assembly_mode}> invalid, exiting"
    }

    println "QC timestamp"

    original_clean_reads.view()

    emit:
    concatenated_reads
    original_clean_reads

}