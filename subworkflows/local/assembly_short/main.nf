#!/usr/bin/env nextflow

include { MEGAHIT } from '../../../modules/local/megahit/main.nf'
include { SPADES as METASPADES } from '../../../modules/local/spades/main.nf'
//include { metaHipMer } from '../../../modules/local/hipmer/main.nf'
include { GATB_MINIA } from '../../../modules/local/gatb/minia/main.nf'

workflow ASSEMBLY_SHORT {
    take:
    clean_reads
    original_clean_reads

    main:

    assembly_out = Channel.empty()
    assembly_graph_out = Channel.empty()
    unitigs_out = Channel.empty()
    
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
            params.megahit_preset
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

        gatb_input_ch.view()

        // Run GATB-minia-pipeline  
        GATB_MINIA (
            gatb_input_ch
        )

        assembly_out = assembly_out.mix ( GATB_MINIA.out.final_contigs )
        assembly_graph_out = assembly_graph_out.mix ( GATB_MINIA.out.assembly_graph )
        unitigs_out = unitigs_out.mix ( GATB_MINIA.out.unitigs )
    }

    //"SPAdes" for tiara domain classification

    emit:
    contigs = assembly_out
    assembly_graph = assembly_graph_out
    unitigs = unitigs_out
    reads = original_clean_reads
}