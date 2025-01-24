#!/usr/bin/env nextflow



include { metabat2 } from '../../../modules/local/metabat2/metabat2/main.nf'
//include { build_index } from '../../../modules/local/bowtie2/build_index/main.nf'
//include { bowtie2_align } from '../../../modules/local/bowtie2/bowtie2_align/main.nf'
include { jgiSummarizeBamContigDepths } from '../../../modules/local/metabat2/jgiSummarizeBamContigDepths/main.nf'
include { convert_depths } from '../../../modules/local/scripts/convert_depths/main.nf'
include { maxbin2 } from '../../../modules/local/maxbin2/main.nf'
include { ADJUST_MAXBIN2_EXT } from '../../../modules/local/scripts/rename_maxbin2_bins'

include { FASTA_BINNING_CONCOCT } from '../../nf-core/fasta_binning_concoct/main.nf'

include { SPLIT_FASTA } from '../../../modules/local/scripts/split_fasta/main.nf'
include { GUNZIP as GUNZIP_BINS } from '../../../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_UNBINS } from '../../../modules/nf-core/gunzip/main.nf'


workflow BINNING {
    take:
    assemblies
    clean_reads
    bigThreads

    main:

    //clean_reads.view()
    //final_assembly.view()

    // Binning preparation with bowtie2
    //build_index(final_assembly)

    //align_input_ch = build_index.out.join(clean_reads, remainder: true)
    //align_input_ch = build_index.out.combine(clean_reads)
    //    .map { meta_assembly, assembly, index, meta_reads, reads ->
    //       tuple(meta_assembly, assembly, index, reads)
    //    }

    //align_input_ch.view()

    //bowtie2_align(align_input_ch,
    //            bigThreads)



    jgi_input_ch = assemblies
        .map {meta, assembly, bams, bais ->
            [meta, bams, bais]
        }

    jgiSummarizeBamContigDepths(jgi_input_ch,
                                bigThreads)

    ch_metabat_depths = jgiSummarizeBamContigDepths.out.depth
        .map { meta, depths ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [ meta_new, depths ]
        }
        

    // main bins for decompressing for MAG_DEPTHS
    final_bins_for_gunzip = Channel.empty()

    // final gzipped bins
    binning_results_gzipped_final = Channel.empty()


    // Create metabat2 input
    metabat2_input_ch = assemblies
        .map {meta, assembly, bams, bais ->
            def meta_new = meta + [binner: 'MetaBAT2']
            [ meta_new, assembly, bams, bais ]
        }
        .join(ch_metabat_depths, by: 0)
        .map { meta, assembly, bams, bais, depths ->
            [ meta, assembly, depths ]
        }


    // Bin assembled contigs with MetaBAT2;
    metabat2(metabat2_input_ch,
            bigThreads)
    final_bins_for_gunzip = final_bins_for_gunzip.mix( metabat2.out.bins.transpose() )
    binning_results_gzipped_final = binning_results_gzipped_final.mix( metabat2.out.bins )


    // nf_core MAG convert_depths for MaxBin2
    convert_depths(metabat2_input_ch)
    maxbin2_input_ch = convert_depths.out.output
        .map { meta, assembly, reads, abund ->
            meta.binner = 'MaxBin2'
            [meta, assembly, reads, abund]
        }

    //maxbin2_input_ch.view()

    maxbin2(maxbin2_input_ch,
            bigThreads)

    ADJUST_MAXBIN2_EXT(maxbin2.out.bins)

    final_bins_for_gunzip = final_bins_for_gunzip.mix( ADJUST_MAXBIN2_EXT.out.renamed_bins.transpose() )
    binning_results_gzipped_final = binning_results_gzipped_final.mix( ADJUST_MAXBIN2_EXT.out.renamed_bins )
    

    // Create CONCOCT input
    concoct_input_ch = assemblies
        .map {meta, assembly, bams, bais ->
            def meta_new = meta + [binner: 'CONCOCT']
            [meta_new, assembly, bams, bais]
        }
        .multiMap {
            meta, assembly, bams, bais ->
                assembly: [meta, assembly]
                bams: [meta, bams, bais]
        }

    //concoct_input_ch.assembly.view()
    //concoct_input_ch.bams.view()

    FASTA_BINNING_CONCOCT (concoct_input_ch.assembly,
                            concoct_input_ch.bams)

    final_bins_for_gunzip = final_bins_for_gunzip.mix( FASTA_BINNING_CONCOCT.out.bins.transpose() )
    binning_results_gzipped_final = binning_results_gzipped_final.mix( FASTA_BINNING_CONCOCT.out.bins )


    // decide which unbinned fasta files to further filter, depending on which binners selected
    // NOTE: CONCOCT does not produce 'unbins' itself, therefore not included here.
    /*if ( !params.skip_metabat2 && params.skip_maxbin2 ) {
        ch_input_splitfasta = metabat2.out.unbinned
    } else if ( params.skip_metabat2 && !params.skip_maxbin2 ) {
        ch_input_splitfasta = maxbin2.out.unbinned_fasta
    } else if ( params.skip_metabat2 && params.skip_maxbin2 ) {
        ch_input_splitfasta = Channel.empty()
    } else {*/
        ch_input_splitfasta = metabat2.out.unbinned.mix(maxbin2.out.unbinned_fasta)
    //}

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

    println "Binning timestamp"

    emit:
    bins                                         = binning_results_gunzipped
    bins_gz                                      = binning_results_gzipped_final
    unbinned                                     = ch_splitfasta_results_gunzipped
    unbinned_gz                                  = SPLIT_FASTA.out.unbinned
    metabat2depths                               = jgiSummarizeBamContigDepths.out.depth
    //versions                                     = ch_versions
}