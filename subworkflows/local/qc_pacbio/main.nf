#!/usr/bin/env nextflow

//include { TOULLIGQC as TOULLIGQC_IN } from '../../../modules/nf-core/toulligqc/main.nf'
//include { TOULLIGQC as TOULLIGQC_OUT } from '../../../modules/nf-core/toulligqc/main.nf'
//include { MULTIQC as MULTIQC_IN } from '../../../modules/local/multiqc/main.nf'
//include { MULTIQC as MULTIQC_OUT } from '../../../modules/local/multiqc/main.nf'
include { FASTPLONG } from '../../../modules/local/fastplong/main.nf'
include { MINIMAP2_INDEX } from '../../../modules/local/minimap2/index/main.nf'
include { MINIMAP2_FILTER_HOST } from '../../../modules/local/minimap2/filter_host/main.nf'

workflow QC_PACBIO {
    take:
    pacbio_reads
    host_genome

    main:
    // Trim reads with FastpLong
    trimmed_pacbio_reads = channel.empty()
    if ( !params.skip_pb_fastplong ) {
        FASTPLONG (
            pacbio_reads,
            params.fplong_disable_quality,
            params.fplong_n_base_limit,
            params.fplong_n_per_limit,
            params.fplong_min_phred,
            params.fplong_unqualified_per,
            params.fplong_mean_qual,
            params.fplong_min_length,
            params.fplong_max_length
        )

        trimmed_pacbio_reads = trimmed_pacbio_reads.mix ( FASTPLONG.out.reads )
            .map { meta, reads ->
                def meta_new = meta + [ corrected: true ]
                [ meta_new, reads ] 
            }
    } else {
        trimmed_pacbio_reads = trimmed_pacbio_reads.mix ( pacbio_reads )
            .map { meta, reads ->
                if ( params.pacbio_hifi || params.pacbio_reads_corrected ) {
                    def meta_new = meta + [ corrected: true ]
                    [ meta_new, reads ] 
                } else {
                    [ meta, reads ]
                }
            }
    }

    // Remove host reads with Minimap2 + Samtools
    filtered_pacbio_reads = channel.empty()
    if ( params.host_genome ) {
        minimap2_index_input = host_genome
            .map { meta, genome ->
                def pacbio_preset = ''
                if ( !params.minimap2_pacbio_preset ) {
                    if ( !params.skip_pb_fastplong || params.pacbio_hifi || params.pacbio_reads_corrected ) {
                        pacbio_preset = "map-hifi" // (<1% error)
                    } else {
                        pacbio_preset = "map-pb" // (<10% error)
                    }
                } else {
                    pacbio_preset = params.minimap2_pacbio_preset
                }
                [ meta, genome, pacbio_preset ]
            }

        // Build index for host genome
        MINIMAP2_INDEX (
            minimap2_index_input
        )
    
        minimap2_filter_input = trimmed_pacbio_reads.combine ( MINIMAP2_INDEX.out.index )

        minimap2_filter_input.view()
    
        //Align ONT reads to host genome;
        MINIMAP2_FILTER_HOST (
            minimap2_filter_input,
            "corrected"
        )

        filtered_pacbio_reads = filtered_pacbio_reads.mix ( MINIMAP2_FILTER_HOST.out.filtered )
            .map { meta, reads ->
                def meta_new = meta + [filtered_host: true]
                [ meta_new, reads ] 
            }
    } else {
        filtered_pacbio_reads = filtered_pacbio_reads.mix ( trimmed_pacbio_reads )
            .map { meta, reads ->
                if ( params.pacbio_reads_corrected ) {
                    def meta_new = meta + [filtered_host: true]
                    [ meta_new, reads ] 
                } else {
                    def meta_new = meta + [filtered_host: false]
                    [ meta_new, reads ]
                }
            }
    }

    emit:
    filtered_pacbio_reads

}