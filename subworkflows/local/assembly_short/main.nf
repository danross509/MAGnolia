#!/usr/bin/env nextflow

include { MEGAHIT } from '../../../modules/local/megahit/main.nf'
include { SPADES as METASPADES } from '../../../modules/local/spades/main.nf'
include { GATB_MINIA_PIPELINE } from '../../../modules/local/gatb/minia_pipeline/main.nf'
include { QUAST_CONTIGS } from '../../../modules/local/quast/quast_contigs/main.nf'
include { CONTIG_COVERAGE } from '../../../subworkflows/local/contig_coverage/main.nf'

workflow ASSEMBLY_SHORT {
    take:
    clean_reads
    original_clean_reads

    main:

    assembly_out = channel.empty()
    assembly_graph_out = channel.empty()
    
    // Short read assembly with megahit
    if (params.assembler_short_reads == 'megahit') {
        megahit_input_ch = clean_reads
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'MEGAHIT']
                [ meta_new, reads ]
            }

        original_clean_reads = original_clean_reads
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'MEGAHIT']
                [ meta_new, reads ]
            }

        // Run MEGAHIT    
        MEGAHIT (
            megahit_input_ch,
            params.megahit_preset
        )

        assembly_out = assembly_out.mix ( MEGAHIT.out.final_contigs )
        assembly_graph_out = assembly_graph_out.mix ( MEGAHIT.out.assembly_graph )

    } else if ( params.assembler_short_reads == 'spades' ) {
    // Short read assembly with metaSPAdes
        // Establish use of --meta option : paired reads only
        if ( params.paired_short_reads && params.spades_use_meta ) {
            use_meta=true
        } else {
            use_meta=null
        }

        spades_input_ch = clean_reads
        .map { meta, reads ->
            def meta_new = meta + [assembler: 'SPAdes']
            [ meta_new, reads, [], [] ]
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
        
    } else if ( params.assembler_short_reads == 'gatb' ){
        gatb_input_ch = clean_reads
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'GATB']
                [ meta_new, reads ]
            }

        original_clean_reads = original_clean_reads
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'GATB']
                [ meta_new, reads ]
            }

        // Run GATB-minia-pipeline  
        GATB_MINIA_PIPELINE (
            gatb_input_ch
        )

        assembly_out = assembly_out.mix ( GATB_MINIA_PIPELINE.out.final_contigs )
        assembly_graph_out = assembly_graph_out.mix ( GATB_MINIA_PIPELINE.out.assembly_graph )
    }

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