#!/usr/bin/env nextflow

include { FLYE as METAFLYE } from '../../../modules/local/flye/main.nf'

workflow ASSEMBLY_LONG {
    take:
    clean_reads_long

    main:

    assembly_out = Channel.empty()
    assembly_graph_out = Channel.empty()

    // Long read assembly with Flye
    if (params.assembler_long_reads == 'flye') {
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
                    if ( meta.sequencer = "ONT") {
                        if ( meta.corrected ) {
                            read_type = "nano-corr" // (<3% error)
                        } else if ( params.nanopore_hq ) {
                            read_type = "nano-hq" // (<5% error)
                        } else {
                            read_type = "nano-raw" // (<20% error)
                        }
                    } else if ( meta.sequencer = "PacBio" ) {
                        if ( params.nanopore_hq ) {
                            read_type = "pacbio-hifi" // (<1% error)
                        } else if ( meta.corrected ) {
                            read_type = "pacbio-corr" // (<3% error)
                        } else {
                            read_type = "pacbio-raw" // (<20% error)
                        }
                    } else {
                        exit 1, "ERROR: Meta.sequencer for long-read <${sampleID}> is invalid"
                    }
                } else {
                    if ( params.flye_read_type ) {
                        read_type = params.flye_read_type
                    }
                }
                

                [ meta_new, reads, read_type ]
            }

        // Run Flye    
        METAFLYE(
            flye_input_ch,
            use_meta,
            params.flye_genome_size,
            params.lgThreads
        )

        assembly_out = assembly_out.mix( METAFLYE.out.final_contigs )
        assembly_graph_out = assembly_graph_out.mix( METAFLYE.out.assembly_graph )

    } else if (params.assembler_long_reads == 'canu') {
    
        assembly_out = CANU.out.________
    } else if (params.assembler_long_reads == 'mdbg') {
        mdbg_input_ch = clean_reads_long
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'MDBG']
                [meta_new, reads]
            }
        
        METAMDBG(
            mdbg_input_ch,
        )
    } else {
        println "Unknown long read assembler, please see nextflow.config for instructions"
    }

    println "Long assembly timestamp"

    emit:
    assembly_out
    assembly_graph_out

}