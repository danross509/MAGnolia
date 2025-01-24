#!/usr/bin/env nextflow

include { megahit } from '../modules/local/fastqc/fromPairs/main.nf'


workflow assembly {
    take:
    clean_reads
    mh_preset

    main:
    // Assembly clean reads with megahit;
    megahit(clean_reads,
            mh_preset
            )

    emit:

}