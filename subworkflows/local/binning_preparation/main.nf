#!/usr/bin/env nextflow

/*
 * Binning preparation with Bowtie2
 */

include { SEMIBIN_CONCATENATE_FASTA } from '../../../modules/local/semibin/concatenate_fasta/main.nf'
include { TAXVAMB_CONCATENATE_TAXONOMY } from '../../../modules/local/scripts/taxvamb_concatenate_tax/main.nf'

//include { BOWTIE2_BUILD_INDEX } from '../../../modules/local/bowtie2/build_index/main.nf'
//include { BOWTIE2_ASSEMBLY_ALIGNMENT } from '../../../modules/local/bowtie2/assembly_alignment/main.nf'
//include { BOWTIE2_ASSEMBLY_MAPPED_SORTED } from '../../../modules/local/bowtie2/assembly_mapped_sorted/main.nf'

include { MINIMAP2_INDEX } from '../../../modules/local/minimap2/index/main.nf'
include { MINIMAP2_ASSEMBLY_ALIGNMENT } from '../../../modules/local/minimap2/assembly_alignment/main.nf'
//include { MINIMAP2_ASSEMBLY_MAPPED_SORTED } from '../../../modules/local/minimap2/assembly_mapped_sorted/main.nf'

workflow BINNING_PREPARATION {
    take:
    assemblies                 // channel: [ val( meta ), path( contigs ), path (gfa), path([ [grouped], [reads] ]), path(contig taxonomy) ]


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

    SEMIBIN_CONCATENATE_FASTA (
        bin_group_contigs
    )

    // Build minimap2 index for each contigs
    build_index_input = SEMIBIN_CONCATENATE_FASTA.out.concatenated_fasta
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

    MINIMAP2_INDEX (
        build_index_input
    )

    // FINAL OUTPUT SHOULD BE ONE CONTIG.FA, *_SORTED.BAM OF EACH SAMPLE IN CONTIG.FA
    //      ONE CONTIG + ONE BAM FOR "PER SAMPLE"
    //      ONE CONTIG + ALL BAMS FOR "COASSEMBLY" (OR FOR EACH "GROUPED ASSEMBLY")
    //      ONE CONCATENATED CONTIG + ALL BAMS FOR COBINNING (OR FOR EACH "GROUPED BINNING")
    //          GIVEN SAMPLES WERE INDIVIDUALLY ASSEMBLED



    split_reads = bin_group_reads
        //.view { meta, reads -> "Before 1st unwrap: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        .map { meta, reads ->
            def readsList = reads as List
            def unwrapped = readsList.size() == 1 && readsList[0] instanceof List ? readsList[0] : readsList
            [meta, unwrapped]
        }
        //.view { meta, reads -> "After 1st unwrap: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        .transpose()
        //.view { meta, reads -> "After transpose: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        .map { meta, reads ->
            if ( reads instanceof Path ) {
                [ meta, [ reads ]]
            } else {
                def readsList2 = reads as List
                def unwrapped2 = readsList2.size() == 1 && readsList2[0] instanceof List ? readsList2[0] : readsList2
                [ meta, unwrapped2 ]
            }
        }    
        //.view { meta, reads -> "After 2nd unwrap: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        .map { meta, reads ->
            def sampleID = reads[0].getBaseName(2).replaceAll(/_corrected/, '').replaceAll(/_trimmed/, '').replaceAll(/_filtered/, '').replaceAll(/_.$/, '')
            //def sampleID = reads.getBaseName(2).replaceAll(/_.$/, '')
            [ meta, sampleID, reads ]
        }

    // group mappings for one assembly
    grouped_reads = split_reads
        .map { meta, _sampleID, reads ->
            [ meta, reads ]
        }
        .groupTuple()

    grouped_gfa = bin_group_gfa
        //.view { meta, gfa -> "Before 1st unwrap: meta=${meta}, gfa=${gfa}, gfa.class=${gfa.class}, gfa.size=${gfa.size()}, gfa[0].class=${gfa[0].class}" }
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
        //.view { meta, tax -> "Before 1st unwrap: meta=${meta}, tax=${tax}, tax.class=${tax.class}, tax.size=${tax.size()}, tax[0].class=${tax[0].class}" }
        .map { meta, tax ->
            if ( tax.size() == 1 && tax[0] instanceof Path ) {
                [ meta, tax ]
            } else {
                def taxList = tax as List
                def unwrapped = taxList.size() == 1 && taxList[0] instanceof List ? taxList[0] : taxList
                [ meta, unwrapped ]
            }
        }
        //.view { meta, tax -> "After 1st unwrap: meta=${meta}, tax=${tax}, tax.class=${tax.class}, tax.size=${tax.size()}, tax[0].class=${tax[0].class}" }
        //.view { meta, tax -> "After 1st unwrap: meta=${meta}, tax=${tax}, tax.class=${tax.class}, tax.size=${tax.size()}" }
        //.view()
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

    //split_reads.view()

    mapping_input = SEMIBIN_CONCATENATE_FASTA.out.concatenated_fasta
        .join( MINIMAP2_INDEX.out.index )
        .combine ( split_reads, by: 0 )
    //mapping_input = BOWTIE2_BUILD_INDEX.out.assembly_index.combine ( split_reads )

    //mapping_input.view()

    //BOWTIE2_ASSEMBLY_ALIGNMENT( 
    //    mapping_input
    //)
    MINIMAP2_ASSEMBLY_ALIGNMENT( 
        mapping_input,
        params.remove_unmapped
    )

    //MINIMAP2_ASSEMBLY_ALIGNMENT.out.mappings.view()

    TAXVAMB_CONCATENATE_TAXONOMY (
        grouped_tax
    )
    
    ch_grouped_mappings = channel.empty()
    //ch_grouped_mappings = ch_grouped_mappings.mix ( BOWTIE2_ASSEMBLY_ALIGNMENT.out.mappings )
    ch_grouped_mappings = ch_grouped_mappings.mix ( MINIMAP2_ASSEMBLY_ALIGNMENT.out.mappings )
        .groupTuple(by: 0)
        .map { meta, contigs, bams, bais -> 
            [ meta, contigs.sort()[0], bams, bais ] 
        }
        .join ( grouped_gfa, by: 0 )
        //.combine ( split_reads, by: 0 )
        .combine ( grouped_reads, by: 0 )
        //.view { meta, contigs, bams, bais, gfa, reads -> "Bin prep output before unwrap: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        //.map { meta, contigs, bams, bais, gfa, sampleID, reads ->
        .map { meta, contigs, bams, bais, gfa, reads ->
            def readsList = reads as List
            def unwrapped = readsList.size() == 1 && readsList[0] instanceof List ? readsList[0] : readsList
            [ meta, unwrapped, contigs, bams, bais, gfa ]
        }
        //.view { meta, reads, contigs, bams, bais, gfa -> "Bin prep output: meta=${meta}, reads=${reads}, reads.class=${reads.class}, reads.size=${reads.size()}, reads[0].class=${reads[0].class}" }
        .join ( TAXVAMB_CONCATENATE_TAXONOMY.out.concatenated_tax, by: 0 )

    // multiple symlinks to the same assembly -> use first of sorted list
    //bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGNMENT.out.log.map { contigs_meta, reads_meta, log -> [ log ] }
    //bowtie2_assembly_multiqc = MINIMAP2_ASSEMBLY_ALIGNMENT.out.log.map { contigs_meta, reads_meta, log -> [ log ] }


    emit:
    //bowtie2_assembly_multiqc
    //bowtie2_version = BOWTIE2_ASSEMBLY_ALIGNMENT.out.versions
    
    grouped_mappings = ch_grouped_mappings 


}