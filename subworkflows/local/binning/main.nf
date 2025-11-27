#!/usr/bin/env nextflow

include { JGISUMMARIZEBAMCONTIGDEPTHS } from '../../../modules/local/metabat2/jgiSummarizeBamContigDepths/main.nf'
include { CONVERT_DEPTHS } from '../../../modules/nf-core_mag/convert_depths/main.nf'
include { SEMIBIN2 } from '../../../modules/local/semibin/semibin2/main.nf'
include { METABAT2 } from '../../../modules/local/metabat2/metabat2/main.nf'
include { FASTA_BINNING_CONCOCT } from '../../nf-core/fasta_binning_concoct/main.nf'
include { METABINNER } from '../../../modules/local/metabinner/main.nf'
include { MAXBIN2 } from '../../../modules/local/maxbin2/main.nf'
include { ADJUST_MAXBIN2_EXT } from '../../../modules/nf-core_mag/adjust_maxbin2_ext'
include { COMEBIN_RUNCOMEBIN as COMEBIN } from '../../../modules/nf-core/comebin/runcomebin/main.nf'
include { MCDEVOL } from '../../../modules/local/mcdevol/main.nf'
include { LRBINNER } from '../../../modules/local/lrbinner/main.nf'
include { VAMB_CONVERT_ABUNDANCE } from '../../../modules/local/vamb/convert_abundance/main.nf'
include { VAMB_BIN } from '../../../modules/nf-core/vamb/bin/main.nf'
//include { REPBIN }

include { SPLIT_FASTA } from '../../../modules/nf-core_mag/split_fasta/main.nf'
include { GUNZIP as GUNZIP_BINS } from '../../../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_UNBINS } from '../../../modules/nf-core/gunzip/main.nf'


workflow BINNING {
    
    take:
    assembly_alignments
    hifiasm_bins
    

    main:

    jgi_input_ch = assembly_alignments
        .map { meta, _reads, _contigs, bams, bais, _gfa, _tax ->
            [ meta, bams, bais ]
        }

    //assembly_alignments.view()
    //jgi_input_ch.view()

    JGISUMMARIZEBAMCONTIGDEPTHS (
        jgi_input_ch
        )

    ch_metabat_depths = JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
        .map { meta, depths ->
            def meta_new = meta + [binner: 'MetaBAT2'] // Change this?
            [ meta_new, depths ]
        }

    // Main bins for decompressing for MAG_DEPTHS
    initial_bins = channel.empty()
    ch_input_splitfasta = channel.empty()

    // Final gzipped bins
    //binning_results_gzipped_final = Channel.empty()

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
            //binning_results_gzipped_final = binning_results_gzipped_final.mix( SEMIBIN2.out.bins )

        } else {
            exit 1, "ERROR: SemiBin2 environment <${params.semibin_environment}> is invalid"
        }
    }

    // Create COMEbin input
    if ( !params.skip_comebin ) {
        comebin_input = assembly_alignments
            .map { meta, _reads, contigs, bams, _bais, _gfa, _tax ->
                def meta_new = meta + [binner: 'COMEbin']
                [ meta_new, contigs, bams ]
            }

        COMEBIN ( 
            comebin_input,
            params.comebin_batch_size 
            )

        initial_bins = initial_bins.mix( COMEBIN.out.bins )
        //binning_results_gzipped_final = binning_results_gzipped_final.mix( COMEBIN.out.bins )

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

    // Create Vamb / TaxVamb input       ****** Adjust to include taxonomy, name binner TaxVamb if provided
    if ( !params.skip_vamb ) {
        vamb_input = assembly_alignments
            .map { meta, _reads, contigs, _bams, _bais, _gfa, tax ->
                def meta_new = meta + [ binner: 'MetaBAT2' ]
                [ meta_new, contigs, tax ]
            }
            .join ( ch_metabat_depths, by: 0 )
            .map { meta, contigs, tax, depths ->
                println(tax)
                def meta_new = meta + [binner: 'Vamb']
                if ( tax ) {
                    meta_new.binner = 'TaxVamb'
                }     
                [ meta_new, contigs, depths, [], tax ] // nf-core Vamb input includes bams and taxonomy
            }

        VAMB_CONVERT_ABUNDANCE (
            vamb_input
            )

        VAMB_BIN ( 
            VAMB_CONVERT_ABUNDANCE.out,
            params.min_contig_length
            )

        initial_bins = initial_bins.mix ( VAMB_BIN.out.bins )
        //binning_results_gzipped_final = binning_results_gzipped_final.mix( VAMB.out.bins )

    }

    // Create McDevol input
    /*if ( !params.skip_mcdevol ) {
        mcdevol_input = assembly_alignments
            .map { meta, reads, contigs, bams, bais, gfa, tax ->
                def meta_new = meta + [binner: 'McDevol']
                [ meta_new, contigs, bams ]
            }

        MCDEVOL ( 
            mcdevol_input
            )

        //initial_bins = initial_bins.mix( MCDEVOL.out.bins )
        //binning_results_gzipped_final = binning_results_gzipped_final.mix( MCDEVOL.out.bins )

    }*/

    // Bin assembled contigs with MetaBAT2;
    if ( !params.skip_metabat2 ) {
        METABAT2 ( metabat2_input_ch )
        initial_bins = initial_bins.mix ( METABAT2.out.bins )
        ch_input_splitfasta = ch_input_splitfasta.mix ( METABAT2.out.unbinned )
        //binning_results_gzipped_final = binning_results_gzipped_final.mix ( METABAT2.out.bins )
    }

    if ( !params.skip_metabinner ) {
        metabinner_input = metabat2_input_ch
            .map { meta, contigs, depths ->
                def meta_new = meta + [ binner: 'MetaBinner' ]
                [ meta_new, contigs, depths ]
            }

        //metabat2_input_ch.view()

        METABINNER (
            metabinner_input,
            params.metabinner_dataset_scale,
            params.min_contig_length,
            params.metabinner_kmer_length
        )

        //METABINNER.out.bins.view()

        initial_bins = initial_bins.mix ( METABINNER.out.bins )
        //binning_results_gzipped_final = binning_results_gzipped_final.mix ( METABINNER.out.bins )
    }

    // nf_core MAG convert_depths for MaxBin2
    if ( !params.skip_maxbin2 ) {
        CONVERT_DEPTHS ( metabat2_input_ch )
        maxbin2_input_ch = CONVERT_DEPTHS.out.output
            .map { meta, contigs, reads, abund ->
                def meta_new = meta + [ binner: 'MaxBin2' ]
                [ meta_new, contigs, reads, abund ]
            }

    //maxbin2_input_ch.view()
        MAXBIN2 ( maxbin2_input_ch )

        ADJUST_MAXBIN2_EXT ( MAXBIN2.out.bins )

        initial_bins = initial_bins.mix ( ADJUST_MAXBIN2_EXT.out.renamed_bins )
        ch_input_splitfasta = ch_input_splitfasta.mix ( MAXBIN2.out.unbinned_fasta )
        //binning_results_gzipped_final = binning_results_gzipped_final.mix ( ADJUST_MAXBIN2_EXT.out.renamed_bins )
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
        //binning_results_gzipped_final = binning_results_gzipped_final.mix( FASTA_BINNING_CONCOCT.out.bins )
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

        //initial_bins = initial_bins.mix ( LRBINNER.out.bins )
        //binning_results_gzipped_final = binning_results_gzipped_final.mix( FASTA_BINNING_CONCOCT.out.bins )
    }

    // decide which unbinned fasta files to further filter, depending on which binners selected
    // NOTE: CONCOCT does not produce 'unbins' itself, therefore not included here.
    /*if ( !params.skip_metabat2 && params.skip_maxbin2 ) {
        ch_input_splitfasta = METABAT2.out.unbinned
    } else if ( params.skip_metabat2 && !params.skip_maxbin2 ) {
        ch_input_splitfasta = MAXBIN2.out.unbinned_fasta
    } else if ( params.skip_metabat2 && params.skip_maxbin2 ) {
        ch_input_splitfasta = Channel.empty()
    } else {
        ch_input_splitfasta = METABAT2.out.unbinned.mix ( MAXBIN2.out.unbinned_fasta )
    }*/

//    SPLIT_FASTA ( ch_input_splitfasta )
    // large unbinned contigs from SPLIT_FASTA for decompressing for MAG_DEPTHS,
    // first have to separate and re-group due to limitation of GUNZIP module
//    ch_split_fasta_results_transposed = SPLIT_FASTA.out.unbinned.transpose()
    //ch_versions = ch_versions.mix(SPLIT_FASTA.out.versions)

    //GUNZIP_BINS ( initial_bins )
    initial_binning_results = initial_bins
        .groupTuple(by: 0)

//    GUNZIP_UNBINS ( ch_split_fasta_results_transposed )
//    ch_splitfasta_results_gunzipped = GUNZIP_UNBINS.out.gunzip
//        .groupTuple(by: 0)

    //ch_versions = ch_versions.mix(GUNZIP_BINS.out.versions.first())
    //ch_versions = ch_versions.mix(GUNZIP_UNBINS.out.versions.first())

    //binning_results_gunzipped.view()

    //MAG bins output:
        // Transposed bins
        // Add [filename: bin.name] ie. unique name of bin, therefore id_assembler_binner
        // Join with stats to remove small bins
        // Final map [meta, bin]

    emit:
    bins                                         = initial_bins
    //bins_gz                                      = binning_results_gzipped_final
//    unbinned                                     = ch_splitfasta_results_gunzipped
//    unbinned_gz                                  = SPLIT_FASTA.out.unbinned
    metabat2depths                               = JGISUMMARIZEBAMCONTIGDEPTHS.out.depth
    //versions                                     = ch_versions
}