#!/usr/bin/env nextflow

/*
To run this pipeline during development:
From the MAG_Pipeline/pipeline_test/short_only/ folder, run the following command:
    nextflow run ../../setup.nf --illumina data/paired_RX/ [--corrected]

From the MAG_Pipeline/pipeline_test/nanopore_only/ folder, run the following command:
    nextflow run ../../setup.nf --nanopore data/nanopore/ [--corrected]

From the MAG_Pipeline/pipeline_test/pacbio_only/ folder, run the following command:
    nextflow run ../../setup.nf --pacbio $path [--corrected]

From the MAG_Pipeline/pipeline_test/hybrid/ folder, run the following command:
    nextflow run ../../setup.nf --illumina data/paired_RX/ [--nanopore data/nanopore/ --pacbio $path] [--corrected]
*/

include { CONCATENATE_RAW_READS as CONCATENATE_ONT_BARCODES } from './modules/local/scripts/concatenate_raw_reads/main.nf'
include { WRITE_SAMPLES_CSV } from './modules/local/scripts/write_sample_csv/main.nf'
include { REMOVE_TMP_CSV } from './modules/local/scripts/write_sample_csv/remove_tmp.nf'
include { REPLACE_SYMLINKS as REPLACE_CONCATENATED_SYMLINKS_ONT } from './modules/local/scripts/replace_symlinks/main.nf'
include { WRITE_CONFIG } from './modules/local/scripts/write_config/main.nf'


// Default parameters
params.illumina = false
params.nanopore = false
params.pacbio = false
params.corrected = false

workflow {

    /****************
        Data input
     ****************/

    write_csv_input = Channel.empty()

    //Create channel of raw short reads
    short_reads = Channel.empty()
    if ( params.illumina ) {
        def illumina_path = params.illumina.endsWith("/") ? params.illumina : params.illumina + "/"
        short_reads = Channel.fromFilePairs ( "${illumina_path}*{1,2}.fastq.gz", size: -1 )
            .map { sample, reads ->
                def sampleID = sample.replaceAll(/_R$/, '')
                return [ sampleID, reads ]
            }
            .map { sampleID, reads ->
                def meta = [:]
                // Set meta.id
                meta.id = sampleID
                meta.sequencer = 'Illumina'
                // Set meta.paired_end
                if ( reads.size() == 2 ) {
                    meta.paired_end = true
                } else if ( reads.size() == 1 ) {
                    meta.paired_end = false
                } else {
                    exit 1, "ERROR: Check short read input files -> ${sampleID} contains invalid number of samples"
                }
                // Set default correctedness
                if ( !params.corrected ) {
                    meta.corrected = false
                } else {
                    meta.corrected = true
                }
                    
                // Set bin group
                meta.bin_group = sampleID
                // Set assembly group
                meta.assembly_group = sampleID
                reads_path_1 = reads[0].toAbsolutePath().toString()
                if ( reads.size() == 2 ) {
                    reads_path_2 = reads[1].toAbsolutePath().toString()
                    return [ meta, [ reads_path_1, reads_path_2 ]]
                } else if ( reads.size() == 1 ) {
                    return [ meta, [ reads_path_1 ]]
                }
            }

        //short_reads.view()

        short_sample_count = short_reads.count()

        //println short_sample_count

        write_csv_input = write_csv_input.mix ( short_reads.combine ( short_sample_count ) )
    }
    
    //Create channel of concatenated raw ONT (Nanopore) reads
    nanopore_barcodes = Channel.empty()
    if ( params.nanopore ) {
        def nanopore_path = params.nanopore.endsWith("/") ? params.nanopore : params.nanopore + "/"
        nanopore_barcodes = Channel.fromPath ( "${nanopore_path}*/*" )
            .map { reads ->
                def meta = [:]
                // Set barcode folder name as meta.id
                def barcode = reads.getParent().getName()
                meta.id = barcode
                meta.sequencer = "ONT"
                // Set paired_end as false for long reads
                meta.paired_end = false
                // Set default correctedness
                if ( !params.corrected ) {
                    meta.corrected = false
                } else {
                    meta.corrected = true
                }
                // Set bin group
                meta.bin_group = barcode
                // Set assembly group
                meta.assembly_group = barcode
                return [ meta, reads ]
            }
            .groupTuple()

        CONCATENATE_ONT_BARCODES ( 
            nanopore_barcodes,
            nanopore_path,
            ""
            )

        // Replace the symlinks of the concatenated reads with the originals
        move_concatenated_input = CONCATENATE_ONT_BARCODES.out.count()
        REPLACE_CONCATENATED_SYMLINKS_ONT (
            move_concatenated_input,
            "${launchDir}/${nanopore_path}"
        )

        nanopore_reads = CONCATENATE_ONT_BARCODES.out.setup_reads_csv
            .map { meta, path ->
                return [meta, [ path ]]
            }

        nanopore_sample_count = nanopore_barcodes.count()

        write_csv_input = write_csv_input.mix ( nanopore_reads.combine ( nanopore_sample_count ) )         
    }

    //Create channel of concatenated raw PacBio reads
    pacbio_reads = Channel.empty()
    if ( params.pacbio ) {
        pacbio_reads = Channel.fromPath ( params.pacbio_reads )

        pacbio_sample_count = pacbio_reads.count()

        write_csv_input = write_csv_input.mix ( pacbio_reads.combine ( pacbio_sample_count ) )  
    }

    write_csv_input.view()

    // Write csv file of input files
    // This is a channel of input, adding each file in no particular order
    WRITE_SAMPLES_CSV {
        write_csv_input
    }

    // Counting the output of WRITE_SAMPLES_CSV ensures the csv is complete before moving
    remove_tmp_input = WRITE_SAMPLES_CSV.out.count()
    REMOVE_TMP_CSV { 
        remove_tmp_input
    }

    short_config_count = short_reads.count()
    short_config_paired = false
    if ( params.illumina ) {
        short_config_paired = short_reads
            .map { meta, reads ->
                if ( !meta.paired_end ) {
                    return "false"
                } else {
                    return "true"
                }
            }
            .unique()
    }
    nanopore_config_count = nanopore_barcodes.count()
    pacbio_config_count = pacbio_reads.count()
    config_corrected = params.corrected

    short_config_count.view()
    //short_config_paired.view()
    nanopore_config_count.view()
    pacbio_config_count.view()
    //println config_corrected

    WRITE_CONFIG (
        short_config_count,
        short_config_paired,
        nanopore_config_count,
        pacbio_config_count,
        config_corrected

    )
}

