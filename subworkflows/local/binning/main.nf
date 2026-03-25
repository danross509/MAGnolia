#!/usr/bin/env nextflow

include { JGISUMMARIZEBAMCONTIGDEPTHS } from '../../../modules/local/metabat2/jgiSummarizeBamContigDepths/main.nf'
include { CONVERT_DEPTHS } from '../../../modules/nf-core_mag/convert_depths/main.nf'
include { VAMB_CONVERT_ABUNDANCE } from '../../../modules/local/vamb/convert_abundance/main.nf'
include { VAMB_FILTER_TAXONOMY } from '../../../modules/local/scripts/taxvamb_filter_taxonomy/main.nf'

include { SEMIBIN2 } from '../../../modules/local/semibin/semibin2/main.nf'
include { COMEBIN_RUNCOMEBIN as COMEBIN } from '../../../modules/nf-core/comebin/runcomebin/main.nf'
include { VAMB_BIN } from '../../../modules/nf-core/vamb/bin/main.nf'
include { METABAT2 } from '../../../modules/local/metabat2/metabat2/main.nf'
include { METABINNER } from '../../../modules/local/metabinner/main.nf'
include { MAXBIN2 } from '../../../modules/nf-core/maxbin2/main.nf'
include { FASTA_BINNING_CONCOCT } from '../../nf-core/fasta_binning_concoct/main.nf'
include { LRBINNER } from '../../../modules/local/lrbinner/main.nf'

workflow BINNING {
    
    take:
    assembly_alignments
    hifiasm_bins
    
    main:
    jgi_input_ch = assembly_alignments
        .map { meta, _reads, _contigs, bams, bais, _gfa, _tax ->
            [ meta, bams, bais ]
        }

    JGISUMMARIZEBAMCONTIGDEPTHS (
        jgi_input_ch
        )

    ch_metabat_depths = JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [ meta_new, depths ]
        }

    initial_bins = channel.empty()

    if ( !params.skip_hmbin ) {
        initial_bins = initial_bins.mix ( hifiasm_bins )
    }

    // Create SemiBin2/SemiBin input
    if ( !params.skip_semibin2 ) {
        env = ['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'chicken_caecum', 'global']
        if ( env.contains( params.semibin_environment )) {
            // NEEDS TO KNOW: {single_sample | coassembly OR grouped assembly | cobinning OR grouped binning}, {short read | long read} 
            semibin2_input_ch = assembly_alignments
                .map { meta, _reads, contigs, bams, _bais, _gfa, _tax ->
                    def meta_new = meta + [binner: 'SemiBin2']
                    [ meta_new, contigs, bams ]
                }

            SEMIBIN2 (
                semibin2_input_ch,
                params.semibin_environment,
                params.use_semibin1
            )

            initial_bins = initial_bins.mix ( SEMIBIN2.out.bins )

        } else {
            exit 1, "ERROR: SemiBin2 environment <${params.semibin_environment}> is invalid"
        }
    }

    // Create COMEbin input
    if ( !params.skip_comebin ) {
        comebin_input = assembly_alignments
            .map { meta, _reads, contigs, bams, _bais, _gfa, _tax ->
                def meta_new = meta + [ binner: 'COMEbin' ]
                [ meta_new, contigs, bams ]
            }

        COMEBIN ( 
            comebin_input,
            params.comebin_batch_size 
            )

        initial_bins = initial_bins.mix ( COMEBIN.out.bins )

    }

    // Create metabat2 input
    metabat2_input_ch = assembly_alignments
        .map { meta, _reads, contigs, bams, bais, _gfa, _tax ->
            def meta_new = meta + [ binner: 'MetaBAT2' ]
            [ meta_new, contigs, bams, bais ]
        }
        .join ( ch_metabat_depths, by: 0 )
        .map { meta, contigs, _bams, _bais, depths ->
            [ meta, contigs, depths ]
        }

    // Create Vamb / TaxVamb input
    if ( !params.skip_vamb ) {
        vamb_input = assembly_alignments
            .map { meta, _reads, contigs, _bams, _bais, _gfa, tax ->
                def meta_new = meta + [ binner: 'MetaBAT2' ]
                [ meta_new, contigs, tax ]
            }
            .join ( ch_metabat_depths, by: 0 )
            .map { meta, contigs, tax, depths ->
                def meta_new = meta + [binner: 'Vamb']
                if ( tax ) {
                    meta_new.binner = 'TaxVamb'
                }     
                [ meta_new, contigs, depths, [], tax ] // nf-core Vamb input includes bams and taxonomy
            }

        VAMB_CONVERT_ABUNDANCE (
            vamb_input
            )

        VAMB_FILTER_TAXONOMY (
            VAMB_CONVERT_ABUNDANCE.out,
            params.min_contig_length
        )

        VAMB_BIN ( 
            VAMB_FILTER_TAXONOMY.out,
            params.min_contig_length
            )

        initial_bins = initial_bins.mix ( VAMB_BIN.out.bins )

    }

    // Bin assembled contigs with MetaBAT2;
    if ( !params.skip_metabat2 ) {
        METABAT2 ( metabat2_input_ch )
        initial_bins = initial_bins.mix ( METABAT2.out.bins )
    }

    if ( !params.skip_metabinner ) {
        metabinner_input = metabat2_input_ch
            .map { meta, contigs, depths ->
                def meta_new = meta + [ binner: 'MetaBinner' ]
                [ meta_new, contigs, depths ]
            }

        METABINNER (
            metabinner_input,
            params.metabinner_dataset_scale,
            params.min_contig_length,
            params.metabinner_kmer_length
        )

        initial_bins = initial_bins.mix ( METABINNER.out.bins )
    }

    // nf_core MAG convert_depths for MaxBin2
    if ( !params.skip_maxbin2 ) {
        CONVERT_DEPTHS ( metabat2_input_ch )
        maxbin2_input_ch = CONVERT_DEPTHS.out.output
            .map { meta, contigs, reads, abund ->
                def meta_new = meta + [ binner: 'MaxBin2' ]
                [ meta_new, contigs, reads, abund ]
            }

        MAXBIN2 ( maxbin2_input_ch )

        initial_bins = initial_bins.mix ( MAXBIN2.out.binned_fastas )
    }

    // Create CONCOCT input
    if ( !params.skip_concoct ) {
        concoct_input_ch = assembly_alignments
            .map { meta, _reads, contigs, bams, bais, _gfa, _tax ->
                def meta_new = meta + [ binner: 'CONCOCT' ]
                [ meta_new, contigs, bams, bais ]
            }
            .multiMap {
                meta, contigs, bams, bais ->
                    contigs: [ meta, contigs ]
                    bams: [ meta, bams, bais ]
            }

        FASTA_BINNING_CONCOCT ( 
            concoct_input_ch.contigs, 
            concoct_input_ch.bams 
            )

        initial_bins = initial_bins.mix ( FASTA_BINNING_CONCOCT.out.bins )
    }

    // Create LRBinner input
    if ( !params.skip_lrbinner ) {
        lrbinner_input = assembly_alignments
            .map { meta, reads, contigs, _bams, _bais, _gfa, _tax ->
                def meta_new = meta + [ binner: 'LRBinner' ]
                [ meta_new, reads, contigs ]
            }

        LRBINNER ( 
            lrbinner_input
            )

        initial_bins = initial_bins.mix ( LRBINNER.out.bins )
    }

    emit:
    bins                                         = initial_bins
    metabat2depths                               = JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    //versions                                     = ch_versions
}