#!/usr/bin/env nextflow

include { TOULLIGQC as TOULLIGQC_IN } from '../../../modules/nf-core/toulligqc/main.nf'
include { TOULLIGQC as TOULLIGQC_OUT } from '../../../modules/nf-core/toulligqc/main.nf'
include { MULTIQC as MULTIQC_IN } from '../../../modules/local/multiqc/main.nf'
include { MULTIQC as MULTIQC_OUT } from '../../../modules/local/multiqc/main.nf'
include { FASTPLONG } from '../../../modules/local/fastplong/main.nf'
include { MINIMAP2_INDEX } from '../../../modules/local/minimap2/index/main.nf'
include { MINIMAP2_FILTER_HOST } from '../../../modules/local/minimap2/filter_host/main.nf'

workflow QC_NANOPORE {
    take:
    nanopore_reads
    host_genome

    main:
    // Create pre-QC ToulligQC report for input reads;
    // Collect pre-QC ToulligQC output;
    // Summarize with multiqc
    //all_toulligqc_in = channel.empty()
    if ( !params.skip_toulligqc ) {
        TOULLIGQC_IN ( nanopore_reads )
        /*all_toulligqc_in = TOULLIGQC_IN.out.report_html
            .map { meta, report ->
                return [ report ]
            }
            .collect()
        MULTIQC_IN(all_toulligqc_in, "pre-toulligQC")*/
    }

    // Trim reads with FastpLong
    trimmed_nanopore_reads = channel.empty()
    if ( !params.skip_ont_fastplong ) {
        FASTPLONG (
            nanopore_reads,
            params.fplong_disable_quality,
            params.fplong_n_base_limit,
            params.fplong_n_per_limit,
            params.fplong_min_phred,
            params.fplong_unqualified_per,
            params.fplong_mean_qual,
            params.fplong_min_length,
            params.fplong_max_length
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
                    [ meta, reads ]
                }
            }
    }

    // Remove host reads with Minimap2 + Samtools
    filtered_nanopore_reads = channel.empty()
    if ( params.host_genome ) {
        minimap2_index_input = host_genome
            .map { meta, genome ->
                def nanopore_preset = ''
                if ( !params.minimap2_nanopore_preset ) {
                    if ( !params.skip_ont_fastplong || params.nanopore_hq || params.nanopore_reads_corrected ) {
                        nanopore_preset = "lr:hq" // (<1% error)
                    } else {
                        nanopore_preset = "map-ont" // (<10% error)
                    }
                } else {
                    nanopore_preset = params.minimap2_nanopore_preset
                }
                [ meta, genome, nanopore_preset ]
            }

        // Build index for host genome
        MINIMAP2_INDEX (
            minimap2_index_input
        )
    
        minimap2_filter_input = trimmed_nanopore_reads.combine ( MINIMAP2_INDEX.out.index )

        minimap2_filter_input.view()
    
        //Align ONT reads to host genome;
        MINIMAP2_FILTER_HOST (
            minimap2_filter_input,
            "corrected"
        )

        filtered_nanopore_reads = filtered_nanopore_reads.mix ( MINIMAP2_FILTER_HOST.out.filtered )
            .map { meta, reads ->
                def meta_new = meta + [filtered_host: true]
                [ meta_new, reads ] 
            }
    } else {
        filtered_nanopore_reads = filtered_nanopore_reads.mix ( trimmed_nanopore_reads )
            .map { meta, reads ->
                if ( params.nanopore_reads_corrected ) {
                    def meta_new = meta + [filtered_host: true]
                    [ meta_new, reads ] 
                } else {
                    def meta_new = meta + [filtered_host: false]
                    [ meta_new, reads ]
                }
            }
    }



    // Create post-QC toulligQC report for input reads;
    // Collect post-QC touliggQC output;
    // Summarize with multiqc
    //all_toulligqc_out = channel.empty()
    if ( !params.skip_toulligqc ) {
        TOULLIGQC_OUT ( filtered_nanopore_reads )
        /*all_toulligqc_out = TOULLIGQC_OUT.out.report_html
            .map { meta, report ->
                return [ report ]
            }
            .collect()
        MULTIQC_OUT(all_toulligqc_out, "post-toulligQC")*/
    }

    emit:
    filtered_nanopore_reads

}