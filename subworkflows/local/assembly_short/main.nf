#!/usr/bin/env nextflow

include { MEGAHIT } from '../../../modules/local/megahit/main.nf'
include { SPADES as METASPADES } from '../../../modules/local/spades/main.nf'
//include { metaHipMer }
//include { minia }

workflow ASSEMBLY_SHORT {
    take:
    clean_reads
    original_clean_reads

    main:

    assembly_out = Channel.empty()
    assembly_graph_out = Channel.empty()
    
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

        //megahit_input_ch.view()

        // Run MEGAHIT    
        MEGAHIT (
            megahit_input_ch,
            params.megahit_preset,
            params.lgThreads,
            params.lgMem
        )

        assembly_out = assembly_out.mix ( MEGAHIT.out.final_contigs )

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
            params.lgThreads,
            params.lgMem,
            [],
            []
        )

        assembly_out = assembly_out.mix ( METASPADES.out.scaffolds )
        
    } else if (params.assembler_short_reads == 'hipmer') {
        hipmer_input_ch = clean_reads
            .map { meta, reads ->
                def meta_new = meta + [assembler: 'HipMer']
                [ meta_new, reads ]
            }
        
        HIPMER(
            hipmer_input_ch,
            params.hipmer_depths
        )
    } else if ( params.assembler_short_reads == 'gatb' ){
        
    }

    //"SPAdes" for tiara domain classification

    emit:
    contigs = assembly_out
    assembly_graph = assembly_graph_out
    reads = original_clean_reads
}