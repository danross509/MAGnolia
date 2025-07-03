#!/usr/bin/env nextflow

/*
 * Binning preparation with Bowtie2
 */

include { CONCATENATE_FASTA } from '../../../modules/local/semibin/concatenate_fasta/main.nf'

//include { BOWTIE2_BUILD_INDEX } from '../../../modules/local/bowtie2/build_index/main.nf'
//include { BOWTIE2_ASSEMBLY_ALIGNMENT } from '../../../modules/local/bowtie2/assembly_alignment/main.nf'
//include { BOWTIE2_ASSEMBLY_MAPPED_SORTED } from '../../../modules/local/bowtie2/assembly_mapped_sorted/main.nf'

include { MINIMAP2_INDEX } from '../../../modules/local/minimap2/index/main.nf'
include { MINIMAP2_ASSEMBLY_ALIGNMENT } from '../../../modules/local/minimap2/assembly_alignment/main.nf'
//include { MINIMAP2_ASSEMBLY_MAPPED_SORTED } from '../../../modules/local/minimap2/assembly_mapped_sorted/main.nf'

workflow BINNING_PREPARATION {
    take:
    assemblies                 // channel: [ val( meta ), path( contigs ), path([ [grouped], [reads] ]) ]


    main:

    bin_group_contigs = Channel.empty()
    bin_group_reads = Channel.empty()

    bin_group_contigs = bin_group_contigs.mix ( assemblies )
        .map { meta, contigs, reads ->
            if ( !meta.coassembly && meta.id != meta.bin_group ) {
                def meta_new = meta + [ id: meta.bin_group, assembly_group: 'self' ]
                [ meta_new, contigs ]
            } else { 
                println contigs.getExtension()
                [ meta, contigs ]
            }
        }
        .groupTuple()
        .map { meta, contigs ->
            if ( contigs.size() >= 2 ) {
                def meta_new = meta + [ cobinning: true ]
                [ meta_new, contigs ]
            } else {
                def meta_new = meta + [ cobinning: false ]
                [ meta_new, contigs ]
            }
        }

    bin_group_reads = bin_group_reads.mix ( assemblies )
        .map { meta, contigs, reads ->
            if ( !meta.coassembly && meta.id != meta.bin_group ) {
                def meta_new = meta + [ id: meta.bin_group, assembly_group: 'self' ]
                [ meta_new, reads ]
            } else { 
                [ meta, reads ]
            }
        }
        .groupTuple()
        .map { meta, reads ->
            if ( reads.size() >= 2 ) {
                def meta_new = meta + [ cobinning: true ]
                [ meta_new, reads ]
            } else {
                def meta_new = meta + [ cobinning: false ]
                [ meta_new, reads ]
            }
        }

    bin_group_contigs.view()
    bin_group_reads.view()
    //concatenated_assemblies = Channel.empty()

    CONCATENATE_FASTA (
        bin_group_contigs
    )

    // Do I need this step?
    //concatenated_assemblies = concatenated_assemblies.mix ( CONCATENATE_FASTA.out.concatenated_fasta.join ( bin_group_reads ))

    //concatenated_assemblies.view()

    // Build minimap2 index for each contigs
    build_index_input = CONCATENATE_FASTA.out.concatenated_fasta
        .map { meta, contigs ->
            def preset = ''
            if ( meta.sequencer == 'Illumina' ) {
                if ( !params.minimap2_short_read_preset ) {
                    preset = "sr"
                } else {
                    preset = params.minimap2_short_read_preset
                }
            } else if ( meta.sequencer == 'ONT' ) {
                if ( !params.minimap2_nanopore_preset ) {
                    if ( meta.corrected ) {
                        preset = "lr:hq" // (<1% error)
                    } else {
                        preset = "map-ont" // (<10% error)
                    }
                } else {
                    preset = params.minimap2_nanopore_preset
                }
            } else if ( meta.sequencer = "PacBio" ) {
                if ( !params.minimap2_pacbio_preset ) {
                    if ( meta.corrected ) {
                        preset = "map-hifi" // (<1% error)
                    } else {
                        preset = "map-pb" // (<10% error)
                    }
                } else {
                    preset = params.minimap2_pacbio_preset 
                }
            }
            [ meta, contigs, preset ]
        }

    MINIMAP2_INDEX (
        build_index_input
    )

    MINIMAP2_INDEX.out.index.view()
    /*BOWTIE2_BUILD_INDEX (
        build_index_input
    )*/

    // FINAL OUTPUT SHOULD BE ONE CONTIG.FA, *_SORTED.BAM OF EACH SAMPLE IN CONTIG.FA
    //      ONE CONTIG + ONE BAM FOR "PER SAMPLE"
    //      ONE CONTIG + ALL BAMS FOR "COASSEMBLY" (OR FOR EACH "GROUPED ASSEMBLY")
    //      ONE CONCATENATED CONTIG + ALL BAMS FOR COBINNING (OR FOR EACH "GROUPED BINNING")
    //          GIVEN SAMPLES WERE INDIVIDUALLY ASSEMBLED



    split_reads = bin_group_reads
        .view { meta, reads -> "Before 1st unwrap: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        .map { meta, reads ->
            def readsList = reads as List
            def unwrapped = readsList.size() == 1 && readsList[0] instanceof List ? readsList[0] : readsList
            [meta, unwrapped]
        }
        .view { meta, reads -> "After 1st unwrap: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        .transpose()
        .view { meta, reads -> "After transpose: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        .map { meta, reads ->
            if ( reads instanceof Path ) {
                [ meta, [ reads ]]
            } else {
                def readsList2 = reads as List
                def unwrapped2 = readsList2.size() == 1 && readsList2[0] instanceof List ? readsList2[0] : readsList2
                [ meta, unwrapped2 ]
            }
        }    
        .view { meta, reads -> "After 2nd unwrap: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        .map { meta, reads ->
            def sampleID = reads[0].getBaseName(2).replaceAll(/_corrected/, '').replaceAll(/_trimmed/, '').replaceAll(/_filtered/, '').replaceAll(/_.$/, '')
            //def sampleID = reads.getBaseName(2).replaceAll(/_.$/, '')
            [ meta, sampleID, reads ]
        }


    /*
        .map { meta, reads ->
        // Recursively unwrap single-element lists until we get to the actual read data
        def unwrap = { list ->
            while (list.size() == 1 && list[0] instanceof List) {
                list = list[0]
            }
            return list
        }
        
        [meta, unwrap(reads)]
    }
    ____________________
        .map { meta, reads ->
        // Recursively flatten until we have a list of pairs
        def flattenToPairs
        flattenToPairs = { list ->
            if (list.every { it instanceof List && it.size() == 2 && !(it[0] instanceof List) }) {
                return list  // We have pairs, stop flattening
            } else if (list.size() == 1 && list[0] instanceof List) {
                return flattenToPairs(list[0])  // Remove one level of nesting
            } else {
                return list
            }
        }
        
        [meta, flattenToPairs(reads)]
    }
    */

    split_reads.view()

    mapping_input = CONCATENATE_FASTA.out.concatenated_fasta
        .join( MINIMAP2_INDEX.out.index )
        .combine ( split_reads, by: 0 )
    //mapping_input = BOWTIE2_BUILD_INDEX.out.assembly_index.combine ( split_reads )

    mapping_input.view()

    ch_grouped_mappings = Channel.empty()

    //BOWTIE2_ASSEMBLY_ALIGNMENT( 
    //    mapping_input
    //)
    MINIMAP2_ASSEMBLY_ALIGNMENT( 
        mapping_input,
        params.remove_unmapped
    )

    // group mappings for one assembly
    //ch_grouped_mappings = ch_grouped_mappings.mix ( BOWTIE2_ASSEMBLY_ALIGNMENT.out.mappings )
    ch_grouped_mappings = ch_grouped_mappings.mix ( MINIMAP2_ASSEMBLY_ALIGNMENT.out.mappings )
        .groupTuple(by: 0)
        .map { meta, contigs, bams, bais -> 
            [ meta, contigs.sort()[0], bams, bais ] 
        }     
    
    // multiple symlinks to the same assembly -> use first of sorted list
    //bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGNMENT.out.log.map { contigs_meta, reads_meta, log -> [ log ] }
    //bowtie2_assembly_multiqc = MINIMAP2_ASSEMBLY_ALIGNMENT.out.log.map { contigs_meta, reads_meta, log -> [ log ] }


    emit:
    //bowtie2_assembly_multiqc
    //bowtie2_version = BOWTIE2_ASSEMBLY_ALIGNMENT.out.versions
    
    grouped_mappings = ch_grouped_mappings 

}