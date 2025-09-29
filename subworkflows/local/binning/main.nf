#!/usr/bin/env nextflow


include { SEMIBIN2 } from '../../../modules/local/semibin/semibin2/main.nf'
include { METABAT2 } from '../../../modules/local/metabat2/metabat2/main.nf'
//include { build_index } from '../../../modules/local/bowtie2/build_index/main.nf'
//include { bowtie2_align } from '../../../modules/local/bowtie2/bowtie2_align/main.nf'
include { JGISUMMARIZEBAMCONTIGDEPTHS } from '../../../modules/local/metabat2/jgiSummarizeBamContigDepths/main.nf'
include { CONVERT_DEPTHS } from '../../../modules/nf-core_mag/convert_depths/main.nf'
include { MAXBIN2 } from '../../../modules/local/maxbin2/main.nf'
include { ADJUST_MAXBIN2_EXT } from '../../../modules/nf-core_mag/adjust_maxbin2_ext'
//include {McDevol}

include { FASTA_BINNING_CONCOCT } from '../../nf-core/fasta_binning_concoct/main.nf'

include { SPLIT_FASTA } from '../../../modules/nf-core_mag/split_fasta/main.nf'
include { GUNZIP as GUNZIP_BINS } from '../../../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_UNBINS } from '../../../modules/nf-core/gunzip/main.nf'


workflow BINNING {
    
    take:
    assembly_alignments
    assembly_graphs
    

    main:

    jgi_input_ch = assembly_alignments
        .map {meta, contigs, bams, bais ->
            [meta, bams, bais]
        }

    jgi_input_ch.view()

    JGISUMMARIZEBAMCONTIGDEPTHS(jgi_input_ch,
                                params.lowThreads)

    ch_metabat_depths = JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [ meta_new, depths ]
        }
        

    // Main bins for decompressing for MAG_DEPTHS
    final_bins_for_gunzip = Channel.empty()

    // Final gzipped bins
    binning_results_gzipped_final = Channel.empty()

    // Create SemiBin2/SemiBin input
    if ( !params.skip_semibin2 ) {
        env = ['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'chicken_caecum', 'global']
        if ( env.contains(params.semibin_environment) ) {
            // NEEDS TO KNOW: {single_sample | coassembly OR grouped assembly | cobinning OR grouped binning}, {short read | long read} 
            semibin2_input_ch = assembly_alignments
                .map { meta, contigs, bams, bais ->
                    def meta_new = meta + [binner: 'SemiBin2']
                    [ meta_new, contigs, bams ]
                }

            SEMIBIN2 (
                semibin2_input_ch,
                params.semibin_environment,
                params.use_semibin1
            )

            //final_bins_for_gunzip = final_bins_for_gunzip.mix( SEMIBIN2.out.bins.transpose() )
            //binning_results_gzipped_final = binning_results_gzipped_final.mix( SEMIBIN2.out.bins )

        } else {
            exit 1, "ERROR: SemiBin2 environment <${params.semibin_environment}> is invalid"
        }
    }

        

    // Create metabat2 input
    metabat2_input_ch = assembly_alignments
        .map {meta, contigs, bams, bais ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [ meta_new, contigs, bams, bais ]
        }
        .join(ch_metabat_depths, by: 0)
        .map { meta, contigs, bams, bais, depths ->
            [ meta, contigs, depths ]
        }

    // Bin assembled contigs with MetaBAT2;
    if ( !params.skip_metabat2 ) {
        METABAT2(metabat2_input_ch,
            params.lgThreads)
        final_bins_for_gunzip = final_bins_for_gunzip.mix( METABAT2.out.bins.transpose() )
        binning_results_gzipped_final = binning_results_gzipped_final.mix( METABAT2.out.bins )
    }

    // nf_core MAG convert_depths for MaxBin2
    if ( !params.skip_maxbin2 ) {
        CONVERT_DEPTHS(metabat2_input_ch)
        maxbin2_input_ch = CONVERT_DEPTHS.out.output
            .map { meta, contigs, reads, abund ->
                def meta_new = meta + [binner: 'MaxBin2']
                [meta_new, contigs, reads, abund]
            }

    //maxbin2_input_ch.view()
        MAXBIN2(maxbin2_input_ch,
            params.lgThreads)

        ADJUST_MAXBIN2_EXT(MAXBIN2.out.bins)

        final_bins_for_gunzip = final_bins_for_gunzip.mix( ADJUST_MAXBIN2_EXT.out.renamed_bins.transpose() )
        binning_results_gzipped_final = binning_results_gzipped_final.mix( ADJUST_MAXBIN2_EXT.out.renamed_bins )
    }

    // Create CONCOCT input
    if ( !params.skip_concoct ){
        concoct_input_ch = assembly_alignments
            .map {meta, contigs, bams, bais ->
                def meta_new = meta + [binner: 'CONCOCT']
                [meta_new, contigs, bams, bais]
            }
            .multiMap {
                meta, contigs, bams, bais ->
                    contigs: [meta, contigs]
                    bams: [meta, bams, bais]
            }

        FASTA_BINNING_CONCOCT (concoct_input_ch.contigs,
                            concoct_input_ch.bams)

        final_bins_for_gunzip = final_bins_for_gunzip.mix( FASTA_BINNING_CONCOCT.out.bins.transpose() )
        binning_results_gzipped_final = binning_results_gzipped_final.mix( FASTA_BINNING_CONCOCT.out.bins )
    }

    // decide which unbinned fasta files to further filter, depending on which binners selected
    // NOTE: CONCOCT does not produce 'unbins' itself, therefore not included here.
    if ( !params.skip_metabat2 && params.skip_maxbin2 ) {
        ch_input_splitfasta = METABAT2.out.unbinned
    } else if ( params.skip_metabat2 && !params.skip_maxbin2 ) {
        ch_input_splitfasta = MAXBIN2.out.unbinned_fasta
    } else if ( params.skip_metabat2 && params.skip_maxbin2 ) {
        ch_input_splitfasta = Channel.empty()
    } else {
        ch_input_splitfasta = METABAT2.out.unbinned.mix(MAXBIN2.out.unbinned_fasta)
    }

    SPLIT_FASTA ( ch_input_splitfasta )
    // large unbinned contigs from SPLIT_FASTA for decompressing for MAG_DEPTHS,
    // first have to separate and re-group due to limitation of GUNZIP module
    ch_split_fasta_results_transposed = SPLIT_FASTA.out.unbinned.transpose()
    //ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)

    GUNZIP_BINS ( final_bins_for_gunzip )
    binning_results_gunzipped = GUNZIP_BINS.out.gunzip
        .groupTuple(by: 0)

    GUNZIP_UNBINS ( ch_split_fasta_results_transposed )
    ch_splitfasta_results_gunzipped = GUNZIP_UNBINS.out.gunzip
        .groupTuple(by: 0)

    //ch_versions = ch_versions.mix(GUNZIP_BINS.out.versions.first())
    //ch_versions = ch_versions.mix(GUNZIP_UNBINS.out.versions.first())

    binning_results_gunzipped.view()

    emit:
    bins                                         = binning_results_gunzipped
    bins_gz                                      = binning_results_gzipped_final
    unbinned                                     = ch_splitfasta_results_gunzipped
    unbinned_gz                                  = SPLIT_FASTA.out.unbinned
    metabat2depths                               = JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    //versions                                     = ch_versions
}