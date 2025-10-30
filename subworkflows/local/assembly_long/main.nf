#!/usr/bin/env nextflow

include { HIFIASM_META } from '../../../modules/local/hifiasm_meta/main.nf'
include { HIFIASM_CONTIG_GFA2FA } from '../../../modules/local/scripts/hifiasm/hifiasm_contig_gfa2fa.nf'
include { HIFIASM_CREATE_BIN_FILES } from '../../../modules/local/scripts/hifiasm/hifiasm_create_bin_files.nf'
include { FLYE as METAFLYE } from '../../../modules/local/flye/main.nf'
include { METAMDBG } from '../../../modules/local/metamdbg/main.nf'
include { CONTIG_POLISHING } from '../../../subworkflows/local/contig_polishing/main.nf'
include { QUAST_CONTIGS } from '../../../modules/local/quast/quast_contigs/main.nf'

workflow ASSEMBLY_LONG {
    take:
    clean_reads_long
    original_clean_reads_long

    main:

    assembly_out = Channel.empty()
    assembly_graph_out = Channel.empty()
    unitigs = Channel.empty()
    reads_out = Channel.empty()

    // Long read assembly with Hifiasm_meta
    hifiasm_bins = Channel.empty()
    if (params.assembler_long_reads == 'hifiasm') {
        hifiasm_input_ch = clean_reads_long
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'Hifiasm_meta']
                [meta_new, reads]
            }

        reads_out = reads_out.mix ( original_clean_reads_long )
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'Hifiasm_meta']
                [meta_new, reads]
            }

        HIFIASM_META (
            hifiasm_input_ch,
            params.skip_hmbin,
            params.hifiasm_read_selection,
            params.hifiasm_force_rs,
            params.hifiasm_rs_threshold
        )

        HIFIASM_CONTIG_GFA2FA (
            HIFIASM_META.out.primary_contig_graph
        )

        HIFIASM_CREATE_BIN_FILES (
            HIFIASM_CONTIG_GFA2FA.out.final_contigs,
            HIFIASM_META.out.circle_contigs,
            HIFIASM_META.out.bins
        )

        assembly_out = assembly_out.mix ( HIFIASM_CONTIG_GFA2FA.out.final_contigs )
        assembly_graph_out = assembly_graph_out.mix ( HIFIASM_META.out.primary_contig_graph )
        unitigs = unitigs.mix ( HIFIASM_META.out.cleaned_unitig_graph )
        hifiasm_bins = hifiasm_bins.mix ( HIFIASM_CREATE_BIN_FILES.out.circular_mags, HIFIASM_CREATE_BIN_FILES.out.bins )
            .map { meta, bins ->
                def meta_new = meta + [cobinning: false, binner: 'hmBin']
                [ meta_new, bins ]
            }

        HIFIASM_CREATE_BIN_FILES.out.circular_mags.view()

    // Long read assembly with Flye
    } else if (params.assembler_long_reads == 'flye') {
        // Run Flye with "--meta" option?
        if (params.flye_use_meta) {
            use_meta=true
        } else {
            use_meta=null
        }

        flye_input_ch = clean_reads_long
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'Flye']
                // Determine flye read type
                def read_type = ''
                if ( !params.flye_read_type ) {
                    if ( meta.sequencer == "ONT") {
                        if ( meta.corrected ) {
                            read_type = "nano-corr" // (<3% error)
                        } else if ( params.nanopore_hq ) {
                            read_type = "nano-hq" // (<5% error)
                        } else {
                            read_type = "nano-raw" // (<20% error)
                        }
                    } else if ( meta.sequencer == "PacBio" ) {
                        if ( params.pacbio_hifi ) {
                            read_type = "pacbio-hifi" // (<1% error)
                        } else if ( meta.corrected ) {
                            read_type = "pacbio-corr" // (<3% error)
                        } else {
                            read_type = "pacbio-raw" // (<20% error)
                        }
                    } else {
                        exit 1, "ERROR: Meta.sequencer for long-read <${meta.id}> is invalid"
                    }
                } else {
                    if ( params.flye_read_type ) {
                        read_type = params.flye_read_type
                    }
                }
                [ meta_new, reads, read_type ]
            }

        reads_out = reads_out.mix ( original_clean_reads_long )
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'Flye']
                [meta_new, reads]
            }

        // Run Flye    
        METAFLYE (
            flye_input_ch,
            use_meta,
            params.flye_genome_size
        )

        assembly_out = assembly_out.mix ( METAFLYE.out.final_contigs )
        assembly_graph_out = assembly_graph_out.mix ( METAFLYE.out.assembly_graph )

    // Long read assembly with metaMDBG
    } else if (params.assembler_long_reads == 'mdbg') {
        mdbg_input_ch = clean_reads_long
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'metaMDBG']
                [meta_new, reads]
            }

        reads_out = reads_out.mix ( original_clean_reads_long )
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'metaMDBG']
                [meta_new, reads]
            }

        METAMDBG (
            mdbg_input_ch,
            params.mdbg_generate_assembly_graph,
            params.mdbg_assembly_graph_k,
            params.mdbg_assembly_graph_contigs,
            params.mdbg_assembly_graph_reads
        )

        assembly_out = assembly_out.mix ( METAMDBG.out.final_contigs )
        assembly_graph_out = assembly_graph_out.mix ( METAMDBG.out.assembly_graph )

    //} else if (params.assembler_long_reads == 'canu') {
    
        // Canu assembly module to be incorporated later

    } else {
        println "Unknown long read assembler, please see nextflow.config for instructions"
    }

    // Polish nanopore contigs

    final_contigs = Channel.empty()
    assembly_graphs = Channel.empty()
    final_reads = Channel.empty()

    if ( params.nanopore_reads && !params.skip_contig_polishing ) {
        // Isolate ONT assemblies
        ont_polishing_contigs = assembly_out
            .map { meta, contigs ->
                if ( meta.sequencer == 'ONT' ) {
                    return [ meta, contigs ]
                } 
            }

        ont_polishing_gfa = assembly_graph_out
            .map { meta, graphs ->
                if ( meta.sequencer == 'ONT' ) {
                    return [ meta, graphs ]
                } 
            }

        ont_polishing_reads = reads_out
            .map { meta, reads ->
                if ( meta.sequencer == 'ONT' ) {
                    return [ meta, reads ]
                } 
            }

        CONTIG_POLISHING (
            ont_polishing_contigs, 
            ont_polishing_gfa,
            ont_polishing_reads
        )

        // Isolate non-ONT assemblies, remix with polished ONT
        final_contigs = final_contigs.mix ( assembly_out )
            .map { meta, contigs ->
                if ( meta.sequencer != 'ONT' ) {
                    return [ meta, contigs ]
                } 
            }
            .mix ( CONTIG_POLISHING.out.final_contigs )

        assembly_graphs = assembly_graphs.mix ( assembly_graph_out )
            .map { meta, graphs ->
                if ( meta.sequencer != 'ONT' ) {
                    return [ meta, graphs ]
                } 
            }
            .mix ( CONTIG_POLISHING.out.corresponding_gfa )

        final_reads = final_reads.mix ( reads_out )
            .map { meta, reads ->
                if ( meta.sequencer != 'ONT' ) {
                    return [ meta, reads ]
                } 
            }
            .mix (CONTIG_POLISHING.out.corresponding_reads )

    } else {
        final_contigs = final_contigs.mix ( assembly_out )
        assembly_graphs = assembly_graphs.mix ( assembly_graph_out )
        final_reads = final_reads.mix ( reads_out )
    }

    if ( !params.skip_quast_contigs ) {
        QUAST_CONTIGS ( 
            final_contigs,
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
    

    coverage_input = Channel.empty()
    if ( !params.skip_contig_coverage ) {
        coverage_input = coverage_input.mix ( final_contigs.join ( final_reads ))

        CONTIG_COVERAGE ( coverage_input )
    }

    emit:
    final_contigs
    assembly_graphs
    unitigs
    final_reads
    hifiasm_bins

}