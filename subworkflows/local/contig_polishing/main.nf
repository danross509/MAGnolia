#!/usr/bin/env nextflow

include { MEDAKA } from '../../../modules/nf-core/medaka/main.nf'


workflow CONTIG_POLISHING {
    take:
    contigs
    assembly_graphs
    reads

    main:

    final_contigs = Channel.empty()
    corresponding_gfa = Channel.empty()
    corresponding_reads = Channel.empty()

    // Contig polishing with medaka
    if (params.contig_polisher == 'medaka') {
        medaka_input = reads.join ( contigs )
            .map { meta, reads, contigs ->
                def meta_new = meta + [polisher: 'medaka']
                [ meta_new, reads, contigs ]
            }

        MEDAKA (
            medaka_input
        )

        final_contigs = final_contigs.mix ( MEDAKA.out.assembly )
        corresponding_gfa = corresponding_gfa.mix ( assembly_graphs )
            .map { meta, graph ->
                def meta_new = meta + [polisher: 'medaka']
                [ meta_new, graph ]
            }
        corresponding_reads = corresponding_reads.mix ( reads )
            .map { meta, reads ->
                def meta_new = meta + [polisher: 'medaka']
                [ meta_new, reads ]
            }

    // Contig polishing with Racon
    /*} else if (params.contig_polisher == 'racon') {
        racon_input = 
  
        RACON (
            racon_input
        )

        final_contigs = final_contigs.mix ( RACON.out )
        corresponding_reads = corresponding_reads.mix ( reads )
            .map {}*/

    // Contig polishing with NanoPolish
    /*} else if (params.contig_polisher == 'nanoPolish') {
        nanopolish_input = 
        
        NANOPOLISH (
            nanopolish_input
        )

        final_contigs = final_contigs.mix ( NANOPOLISH.out )
        corresponding_reads = corresponding_reads.mix ( reads )
            .map {}*/

    } else {
        println "Unknown contig polisher, please see nextflow.config for instructions"
    }

    emit:
    final_contigs
    corresponding_gfa
    corresponding_reads

}