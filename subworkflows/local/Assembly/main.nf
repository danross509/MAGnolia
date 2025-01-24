#!/usr/bin/env nextflow

include { megahit } from '../../../modules/local/megahit/main.nf'


workflow assembly {
    take:
    clean_reads
    mh_preset
    bigThreads
    bigMem

    main:
    // Assembly clean reads with megahit;
    clean_reads.view()

    megahit_input_ch = clean_reads
        .map { meta, reads ->
            def meta_new = meta + [assembler: 'MEGAHIT']
            [meta_new, reads]
        }
        

    megahit(megahit_input_ch,
            mh_preset,
            bigThreads,
            bigMem
            )

    assembly_out = megahit.out.final_contigs 


    //"SPAdes" for tiara domain classification

    println "Assembly timestamp"

    emit:
    assembly_out
}