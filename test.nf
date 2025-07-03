#!/usr/bin/env nextflow

/*
To run this pipeline during development:
From the MAG_Pipeline/pipeline_test/[short_only | nanopore_only | pacbio_only | hybrid] folder, run the following command:
    nextflow run ../../test.nf -resume
*/

include { QC_SHORT } from './subworkflows/local/qc_short/main.nf'

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

    short_reads = Channel.empty()
    short_reads = short_reads.mix ( reads )
        .map { meta, reads ->
            if ( meta.sequencer == 'Illumina' ) {
                return [ meta, reads ]
            } 
        }
    
    phiX = Channel.empty()
    if ( !params.skip_qc && params.short_reads && params.remove_phiX ) {
        phiX = Channel.fromPath ( params.phiX )
            .map { reference ->
                def meta = [:]
                meta.id = reference.getBaseName()
                return [ meta, reference ]
            }

        phiX_index = Channel.fromPath ( "${projectDir}/reference_genomes/phiX/*.bt2" )
            .map { index ->
                def meta =[:]
                meta.id = 'phiX174'
                return [ meta, index ]
            }
            .groupTuple()
    }

    phiX_index.view()

    host_genome = Channel.empty()
    // Create channel of host genome to remove
    if ( !params.skip_qc && params.host_genome ) {
        host_genome = Channel.fromPath ( params.host_genome )
            .map { reference ->
                def meta = [:]
                meta.id = reference.getBaseName()
                return [ meta, reference ]
            }
    }

    /*corrected_reads = Channel.empty()
    if ( !params.skip_qc && params.short_reads ) {
        QC_SHORT ( 
            short_reads,
            phiX,
            host_genome
        )

        corrected_reads = corrected_reads.mix ( QC_SHORT.out.host_filtered_reads )
    }*/

}