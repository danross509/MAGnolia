#!/usr/bin/env nextflow

include { RENAME_CLEAN_READS as RENAME_PACBIO_READS } from '../../../modules/local/scripts/rename_clean_reads/main.nf'
include { RENAME_CLEAN_READS as RENAME_ONT_READS } from '../../../modules/local/scripts/rename_clean_reads/main.nf'
include { SPADES as METASPADES } from '../../../modules/local/spades/main.nf'
include { QUAST_CONTIGS } from '../../../modules/local/quast/quast_contigs/main.nf'
include { CONTIG_COVERAGE } from '../../../subworkflows/local/contig_coverage/main.nf'

workflow ASSEMBLY_HYBRID {
    take:
    clean_reads
    original_clean_reads
    clean_long_reads

    main:

    assembly_out = channel.empty()
    assembly_graph_out = channel.empty()
    
    // Hybrid assembly with metaSPAdes
    // Establish use of --meta option : paired reads only
    if ( params.paired_short_reads && params.spades_use_meta ) {
        use_meta=true
    } else {
        use_meta=null
    }

    // Format input channel for hybrid assembly
    pacbio_reads = clean_long_reads
        .filter { meta, _reads -> meta.sequencer == 'PacBio' }
        .map { meta, reads ->
            [ meta.assembly_group, reads ]
        }
    
    nanopore_reads = clean_long_reads
        .filter { meta, _reads -> meta.sequencer == 'ONT' }
        .map { meta, reads ->
            [ meta.assembly_group, reads ]
        }

    spades_input_ch = clean_reads
        .map { meta, reads ->
            def meta_new = meta + [assembler: 'SPAdes']
            [ meta.assembly_group, meta_new, reads ]
        }
        .join ( pacbio_reads, by: 0, remainder: true )
        .join ( nanopore_reads, by: 0, remainder: true )
        .map { _id, meta, illumina, pacbio, nanopore ->
            def use_pacbio = pacbio ?: []
            def use_nanopore = nanopore ?: []
            [ meta, illumina, use_pacbio, use_nanopore ]
        }

    original_clean_reads = original_clean_reads
        .map { meta, reads ->
            def meta_new = meta + [assembler: 'SPAdes']
            [ meta_new, reads ]
        }

    // Run metaSPAdes
    METASPADES (
        spades_input_ch,
        use_meta,
        [],
        []
    )

    assembly_out = assembly_out.mix ( METASPADES.out.scaffolds )
    assembly_graph_out = assembly_graph_out.mix ( METASPADES.out.gfa )
        

    if ( !params.skip_quast_contigs ) {
        QUAST_CONTIGS ( 
            assembly_out,
            params.quast_assembly_min_contig,
            params.quast_assembly_rna_finding,
            params.quast_assembly_gene_finding,
            params.quast_assembly_conserved_gene_finding,
            params.quast_assembly_min_alignment,
            params.quast_assembly_min_identity,
            params.quast_assembly_mem_efficient,
            params.quast_assembly_space_efficient,
            params.quast_assembly_blast_db
        )
    }
    
    coverage_input = channel.empty()
    if ( !params.skip_contig_coverage ) {
        coverage_input = coverage_input.mix ( assembly_out.join ( original_clean_reads ))

        CONTIG_COVERAGE ( coverage_input )
    }

    emit:
    contigs = assembly_out
    assembly_graph = assembly_graph_out
    reads = original_clean_reads
}