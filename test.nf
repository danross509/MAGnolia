#!/usr/bin/env nextflow

/*
To run this pipeline during development:
From the MAG_Pipeline/pipeline_test/[short_only | nanopore_only | pacbio_only | hybrid] folder, run the following command:
    nextflow run ../../test.nf -resume
*/

workflow {

    reads = Channel.fromPath ( params.reads_file )
        .splitCsv (header: true)
        .map { meta ->
            def meta_new = [:]
            def reads = []

            meta_new.id = meta.sampleID
            meta_new.sequencer = meta.sequencer
            if ( meta.paired_end == 'true' ) {
                meta_new.paired_end = true
            } else {
                meta_new.paired_end = false
            }
            if ( meta.corrected == 'true' ) {
                meta_new.corrected = true
            } else {
                meta_new.corrected = false
            }
            meta_new.bin_group = meta.bin_group
            meta_new.assembly_group = meta.assembly_group

            if ( meta_new.paired_end ) {
                return [ meta_new, [ meta.reads_R1, meta.reads_R2 ]]
            } else {
                return [ meta_new, [ meta.reads_R1 ]]
            }
        }
    
    //reads.view()
    short_reads = Channel.empty()
    short_reads = short_reads.mix ( reads )
        .map { meta, reads ->
            if ( meta.sequencer == 'Illumina' ) {
                return [ meta, reads ]
            } 
        }

    nanopore_reads = Channel.empty()
    nanopore_reads = nanopore_reads.mix ( reads )
        .map { meta, reads ->
            if ( meta.sequencer == 'ONT' ) {
                return [ meta, reads ]
            } 
        }

    pacbio_reads = Channel.empty()
    pacbio_reads = pacbio_reads.mix ( reads )
        .map { meta, reads ->
            if ( meta.sequencer == 'PacBio' ) {
                return [ meta, reads ]
            } 
        }

    //short_reads.view()
    //nanopore_reads.view()

    if ( short_reads.count() == 0 ) {
        short_reads = false
    }

    if ( nanopore_reads.count() == 0 ) {
        nanopore_reads = false
    }

    if ( pacbio_reads.count() == 0 ) {
        pacbio_reads = false
    } 

    short_reads.count().view()
    short_reads.view()
    nanopore_reads.count().view()
    nanopore_reads.view()
    pacbio_reads.count().view()

}