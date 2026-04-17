#!/usr/bin/env nextflow

/*
 * Binning preparation with Bowtie2
 */

include { RENAME_SOLO_CONTIGS } from '../../../modules/local/scripts/rename_solo_contigs/main.nf'

include { SEMIBIN_CONCATENATE_FASTA } from '../../../modules/local/semibin/concatenate_fasta/main.nf'
include { TAXVAMB_CONCATENATE_TAXONOMY } from '../../../modules/local/scripts/taxvamb_concatenate_tax/main.nf'

include { MINIMAP2_INDEX as MINIMAP2_ASSEMBLY_INDEX } from '../../../modules/local/minimap2/index/main.nf'
include { MINIMAP2_ASSEMBLY_ALIGNMENT } from '../../../modules/local/minimap2/assembly_alignment/main.nf'

workflow BINNING_PREPARATION {
    take:
    assemblies      // channel: [ val( meta ), path( contigs ), path (gfa), path([ [grouped], [reads] ]), path(contig taxonomy) ]


    main:

    bin_group_contigs = channel.empty()
    bin_group_gfa = channel.empty()
    bin_group_reads = channel.empty()
    bin_group_tax = channel.empty()
  
    // Group together contigs according to bin group
    // If bin group != meta.id, rename
    bin_group_contigs = bin_group_contigs.mix ( assemblies )
        .map { meta, contigs, _gfa, _reads, _tax ->
            if ( !meta.coassembly && meta.id != meta.bin_group ) {
                def meta_new = meta + [ id: meta.bin_group, assembly_group: 'self' ] // Don't change assembly group
                [ meta_new, contigs ]
            } else { 
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

    bin_group_gfa = bin_group_gfa.mix ( assemblies )
        .map { meta, _contigs, gfa, _reads, _tax ->
            if ( !meta.coassembly && meta.id != meta.bin_group ) {
                def meta_new = meta + [ id: meta.bin_group, assembly_group: 'self' ]
                [ meta_new, gfa ]
            } else { 
                [ meta, gfa ]
            }
        }
        .groupTuple()
        .map { meta, gfa ->
            if ( gfa.size() >= 2 ) {
                def meta_new = meta + [ cobinning: true ]
                [ meta_new, gfa ]
            } else {
                def meta_new = meta + [ cobinning: false ]
                [ meta_new, gfa ]
            }
        }

    bin_group_reads = bin_group_reads.mix ( assemblies )
        .map { meta, _contigs, _gfa, reads, _tax ->
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

    bin_group_tax = bin_group_tax.mix ( assemblies )
        .map { meta, _contigs, _gfa, _reads, tax ->
            if ( !meta.coassembly && meta.id != meta.bin_group ) {
                def meta_new = meta + [ id: meta.bin_group, assembly_group: 'self' ]
                [ meta_new, tax ]
            } else { 
                [ meta, tax ]
            }
        }
        .groupTuple()
        .map { meta, tax ->
            if ( tax.size() >= 2 ) {
                def meta_new = meta + [ cobinning: true ]
                [ meta_new, tax ]
            } else {
                def meta_new = meta + [ cobinning: false ]
                [ meta_new, tax ]
            }
        }

    contigs_cobinning = bin_group_contigs
        .filter { meta, _contigs ->
            meta.cobinning
        }

    contigs_solo = bin_group_contigs
        .filter { meta, _contigs ->
            !meta.cobinning
        }

    SEMIBIN_CONCATENATE_FASTA (
        contigs_cobinning
    )

    RENAME_SOLO_CONTIGS (
        contigs_solo
    )

    // Recombine contigs into single channel to pass through the next steps
    contigs_recombined = SEMIBIN_CONCATENATE_FASTA.out.concatenated_fasta.mix ( RENAME_SOLO_CONTIGS.out.concatenated_fasta )

    // Build minimap2 index for each contigs
    build_index_input = contigs_recombined
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
            } else if ( meta.sequencer == "PacBio" ) {
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

    MINIMAP2_ASSEMBLY_INDEX (
        build_index_input
    )

    // FINAL OUTPUT SHOULD BE ONE CONTIG.FA, *_SORTED.BAM OF EACH SAMPLE IN CONTIG.FA
    //      ONE CONTIG + ONE BAM FOR "PER SAMPLE"
    //      ONE CONTIG + ALL BAMS FOR "COASSEMBLY" (OR FOR EACH "GROUPED ASSEMBLY")
    //      ONE CONCATENATED CONTIG + ALL BAMS FOR COBINNING (OR FOR EACH "GROUPED BINNING"), GIVEN SAMPLES WERE INDIVIDUALLY ASSEMBLED

    split_reads = bin_group_reads
        .map { meta, reads ->
            def readsList = reads as List
            def unwrapped = readsList.size() == 1 && readsList[0] instanceof List ? readsList[0] : readsList
            [meta, unwrapped]
        }
        .transpose()
        .map { meta, reads ->
            if ( reads instanceof Path ) {
                [ meta, [ reads ]]
            } else {
                def readsList2 = reads as List
                def unwrapped2 = readsList2.size() == 1 && readsList2[0] instanceof List ? readsList2[0] : readsList2
                [ meta, unwrapped2 ]
            }
        }
        .map { meta, reads ->
            def sampleID = reads[0].getBaseName(2).replaceAll(/_corrected/, '').replaceAll(/_trimmed/, '').replaceAll(/_filtered/, '').replaceAll(/_.$/, '')
            [ meta, sampleID, reads ]
        }

    // Group mappings for one assembly
    grouped_reads = split_reads
        .map { meta, _sampleID, reads ->
            [ meta, reads ]
        }
        .groupTuple()

    grouped_gfa = bin_group_gfa
        .map { meta, gfa ->
            if ( gfa.size() == 1 && gfa[0] instanceof Path ) {
                [ meta, gfa ]
            } else {
                def gfaList = gfa as List
                def unwrapped = gfaList.size() == 1 && gfaList[0] instanceof List ? gfaList[0] : gfaList
                [ meta, unwrapped ]
            }
        }

    grouped_tax = bin_group_tax
        .map { meta, tax ->
            if ( tax.size() == 1 && tax[0] instanceof Path ) {
                [ meta, tax ]
            } else {
                def taxList = tax as List
                def unwrapped = taxList.size() == 1 && taxList[0] instanceof List ? taxList[0] : taxList
                [ meta, unwrapped ]
            }
        }

    mapping_input = contigs_recombined
        .join( MINIMAP2_ASSEMBLY_INDEX.out.index )
        .combine ( split_reads, by: 0 )

    MINIMAP2_ASSEMBLY_ALIGNMENT( 
        mapping_input,
        params.remove_unmapped
    )

    TAXVAMB_CONCATENATE_TAXONOMY (
        grouped_tax
    )
    
    ch_grouped_mappings = channel.empty()
    ch_grouped_mappings = ch_grouped_mappings.mix ( MINIMAP2_ASSEMBLY_ALIGNMENT.out.mappings )
        .groupTuple(by: 0)
        .map { meta, contigs, bams, bais -> 
            [ meta, contigs.sort()[0], bams, bais ] 
        }
        .join ( grouped_gfa, by: 0 )
        .combine ( grouped_reads, by: 0 )
        .map { meta, contigs, bams, bais, gfa, reads ->
            def readsList = reads as List
            def unwrapped = readsList.size() == 1 && readsList[0] instanceof List ? readsList[0] : readsList
            [ meta, unwrapped, contigs, bams, bais, gfa ]
        }
        .join ( TAXVAMB_CONCATENATE_TAXONOMY.out.concatenated_tax, by: 0 )

    emit:
    grouped_mappings = ch_grouped_mappings 


}