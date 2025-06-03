#!/usr/bin/env nextflow

/*
 * Binning preparation with Bowtie2
 */

include { BOWTIE2_BUILD_INDEX } from '../../../modules/local/bowtie2/build_index/main.nf'
include { BOWTIE2_ASSEMBLY_ALIGNMENT } from '../../../modules/local/bowtie2/assembly_alignment/main.nf'
include { BOWTIE2_ASSEMBLY_MAPPED_SORTED } from '../../../modules/local/bowtie2/assembly_mapped_sorted/main.nf'

include { MINIMAP2_INDEX } from '../../../modules/local/minimap2/index/main.nf'
include { MINIMAP2_ASSEMBLY_ALIGNMENT } from '../../../modules/local/minimap2/assembly_alignment/main.nf'
include { MINIMAP2_ASSEMBLY_MAPPED_SORTED } from '../../../modules/local/minimap2/assembly_mapped_sorted/main.nf'

workflow BINNING_PREPARATION {
    take:
    assemblies                 // channel: [ val( meta ), path( contigs ), path([ [grouped], [reads] ]) ]


    main:
    // build bowtie2 index for all contigs
    build_index_input = assemblies
        .map { contigs_meta, contigs, grouped_reads_meta, grouped_reads ->
            def preset = ''
            if ( contigs_meta.sequencer == 'Illumina' ) {
                if ( !params.minimap2_short_read_preset ) {
                    preset = "sr"
                } else {
                    preset = params.minimap2_short_read_preset
                }
            } else if ( contigs_meta.sequencer == 'ONT' ) {
                if ( !params.minimap2_nanopore_preset ) {
                    if ( contigs_meta.corrected ) {
                        preset = "lr:hq" // (<1% error)
                    } else {
                        preset = "map-ont" // (<10% error)
                    }
                } else {
                    preset = params.minimap2_nanopore_preset
                }
            } else if ( contigs_meta.sequencer = "PacBio" ) {
                if ( !params.minimap2_pacbio_preset ) {
                    if ( contigs_meta.corrected ) {
                        preset = "map-hifi" // (<1% error)
                    } else {
                        preset = "map-pb" // (<10% error)
                    }
                } else {
                    preset = params.minimap2_pacbio_preset 
                }
            }
            [ contigs_meta, contigs, preset ]
        }

    MINIMAP2_INDEX (
        build_index_input
    )
    /*BOWTIE2_BUILD_INDEX (
        build_index_input
    )*/

    // FINAL OUTPUT SHOULD BE ONE CONTIG.FA, *_SORTED.BAM OF EACH SAMPLE IN CONTIG.FA
    //      ONE CONTIG + ONE BAM FOR "PER SAMPLE"
    //      ONE CONTIG + ALL BAMS FOR "COASSEMBLY" (OR FOR EACH "GROUPED ASSEMBLY")
    //      ONE CONCATENATED CONTIG + ALL BAMS FOR COBINNING (OR FOR EACH "GROUPED BINNING")
    //          GIVEN SAMPLES WERE INDIVIDUALLY ASSEMBLED

    split_reads = assemblies
        .map { contigs_meta, contigs, grouped_reads_meta, grouped_reads ->
            //println grouped_reads[0][0]
            [ grouped_reads ]
        }
        .flatten()
        .map { reads ->
            def meta = [:]
            def sampleID = reads.getBaseName(2).replaceAll(/_.$/, '')
            meta.id = sampleID
            [ meta, reads ]
        }
        .groupTuple()

    split_reads.view()

    mapping_input = build_index_input
        .map { contigs_meta, contigs, preset ->
            [ contigs_meta, contigs ]
        }
        .join( MINIMAP2_INDEX.out.index )
        .combine ( split_reads )
    //mapping_input = BOWTIE2_BUILD_INDEX.out.assembly_index.combine ( split_reads )

    mapping_input.view()

    ch_semibin2_mappings = Channel.empty()
    ch_grouped_mappings = Channel.empty()
    if ( !params.skip_semibin2 ) {
        /*BOWTIE2_ASSEMBLY_MAPPED_SORTED( 
            mapping_input
        )*/
        MINIMAP2_ASSEMBLY_MAPPED_SORTED( 
            mapping_input
        )

        //ch_semibin2_mappings = ch_semibin2_mappings.mix( BOWTIE2_ASSEMBLY_MAPPED_SORTED.out.mappings )
        ch_semibin2_mappings = ch_semibin2_mappings.mix( MINIMAP2_ASSEMBLY_MAPPED_SORTED.out.mappings )
            .groupTuple(by: 0)
            .map { meta, contigs, bams, bais -> 
                [ meta, contigs.sort()[0], bams, bais ] 
            }

        //bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_MAPPED_SORTED.out.log.map { contigs_meta, reads_meta, log -> [ log ] }
        //bowtie2_assembly_multiqc = MINIMAP2_ASSEMBLY_MAPPED_SORTED.out.log.map { contigs_meta, reads_meta, log -> [ log ] }

    } else {
        /*BOWTIE2_ASSEMBLY_ALIGNMENT( 
            mapping_input
        )*/
        MINIMAP2_ASSEMBLY_ALIGNMENT( 
            mapping_input
        )

        // group mappings for one assembly
        //ch_grouped_mappings = ch_grouped_mappings.mix ( BOWTIE2_ASSEMBLY_ALIGNMENT.out.mappings )
        ch_grouped_mappings = ch_grouped_mappings.mix ( MINIMAP2_ASSEMBLY_ALIGNMENT.out.mappings )
            .groupTuple(by: 0)
            .map { meta, contigs, bams, bais -> 
                [ meta, contigs.sort()[0], bams, bais ] 
            }     // multiple symlinks to the same assembly -> use first of sorted list

        //bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGNMENT.out.log.map { contigs_meta, reads_meta, log -> [ log ] }
        //bowtie2_assembly_multiqc = MINIMAP2_ASSEMBLY_ALIGNMENT.out.log.map { contigs_meta, reads_meta, log -> [ log ] }
    }




    emit:
    //bowtie2_assembly_multiqc
    //bowtie2_version = BOWTIE2_ASSEMBLY_ALIGNMENT.out.versions
    grouped_mappings = ch_grouped_mappings
    semibin2_mappings = ch_semibin2_mappings    
}