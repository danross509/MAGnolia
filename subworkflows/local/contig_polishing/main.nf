#!/usr/bin/env nextflow

include { MEDAKA } from '../../../modules/nf-core/medaka/main.nf'


workflow CONTIG_POLISHING {
    take:
    input_contigs
    assembly_graphs
    input_reads

    main:

    final_contigs = channel.empty()
    corresponding_gfa = channel.empty()
    corresponding_reads = channel.empty()

    // Contig polishing with medaka
    medaka_input = input_reads.join ( input_contigs )
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
    corresponding_reads = corresponding_reads.mix ( input_reads )
        .map { meta, reads ->
            def meta_new = meta + [polisher: 'medaka']
            [ meta_new, reads ]
        }

    emit:
    final_contigs
    corresponding_gfa
    corresponding_reads

}