#!/usr/bin/env nextflow

/*
To run this pipeline during development:
From the MAG_Pipeline/pipeline_test/short_only/ folder, run the following command:
    nextflow run ../../setup.nf --illumina data/paired_RX/ [--corrected]

From the MAG_Pipeline/pipeline_test/nanopore_only/ folder, run the following command:
    nextflow run ../../setup.nf --nanopore data/nanopore/ [--corrected]

From the MAG_Pipeline/pipeline_test/pacbio_only/ folder, run the following command:
    nextflow run ../../setup.nf --pacbio data/pacbio/ [--corrected]

From the MAG_Pipeline/pipeline_test/hybrid/ folder, run the following command:
    nextflow run ../../setup.nf --illumina data/paired_RX/ [--nanopore data/nanopore/ --pacbio $path] [--corrected]
*/

include { CONCATENATE_ONT_BARCODES } from './modules/local/scripts/concatenate_ont_barcodes/main.nf'
include { PBTK_BAM2FASTQ } from './modules/local/pbtk/bam2fastq/main.nf'
//include { CONFIRM_UNIQUE_FILES as CONFIRM_SAMPLE_ID_ONT } from './modules/local/scripts/confirm_unique_files/main.nf'
include { CONFIRM_UNIQUE_FILES as CONFIRM_SAMPLE_ID_PB } from './modules/local/scripts/confirm_unique_files/main.nf'
include { CREATE_TMP_CSV } from './modules/local/scripts/write_sample_csv/create_tmp.nf'
include { WRITE_SAMPLES_CSV } from './modules/local/scripts/write_sample_csv/main.nf'
include { REMOVE_TMP_CSV } from './modules/local/scripts/write_sample_csv/remove_tmp.nf'
include { REPLACE_SYMLINKS as REPLACE_CONCATENATED_SYMLINKS_ONT } from './modules/local/scripts/replace_symlinks/main.nf'
include { REPLACE_SYMLINKS as REPLACE_FASTQ_SYMLINKS_PB } from './modules/local/scripts/replace_symlinks/main.nf'
include { WRITE_CONFIG } from './modules/local/scripts/write_config/main.nf'


// Default parameters
params.illumina = false
params.nanopore = false
params.pacbio = false

params.corrected = false
params.coassembly = false
params.cobinning = false
params.use_gpu = false

workflow {

    /****************
        Data input
     ****************/

    write_csv_input = channel.empty()

    //Create channel of raw short reads
    short_reads = channel.empty()
    if ( params.illumina ) {
        def illumina_path = params.illumina.endsWith("/") ? params.illumina : params.illumina + "/"
        short_reads = channel.fromFilePairs ( "${illumina_path}*{1,2}.{fastq,fq}.gz", size: -1 )
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
                // Set assembly group
                if ( !params.coassembly ) {
                    meta.assembly_group = sampleID
                } else {
                    meta.assembly_group = "allReads"
                }
                // Set bin group
                if ( !params.cobinning ) {
                    meta.bin_group = sampleID
                } else {
                    meta.bin_group = "allReads"
                }
                

                def reads_path_1 = reads[0].toAbsolutePath().toString()
                if ( reads.size() == 2 ) {
                    def reads_path_2 = reads[1].toAbsolutePath().toString()
                    return [ meta, [ reads_path_1, reads_path_2 ]]
                } else if ( reads.size() == 1 ) {
                    return [ meta, [ reads_path_1 ]]
                }
            }

        //short_reads.view()

        //short_sample_count = short_reads.count()
        //println short_sample_count
        //write_csv_input = write_csv_input.mix ( short_reads.combine ( short_sample_count ) )
        
        write_csv_input = write_csv_input.mix ( short_reads )
    }
    
    //Create channel of concatenated raw ONT (Nanopore) reads
    nanopore_fastqs = channel.empty()
    nanopore_barcodes = channel.empty()
    nanopore_reads = channel.empty()
    if ( params.nanopore ) {
        def nanopore_path = params.nanopore.endsWith("/") ? params.nanopore : params.nanopore + "/"
        // Add option if presented as single fastq's?
        nanopore_fastqs = channel.fromPath ( "${nanopore_path}*.{fastq,fq}.gz" )
            .map { reads ->
                def meta = [:]
                def sampleID = reads.getBaseName(2)
                meta.id = sampleID
                meta.sequencer = "ONT"
                // Set paired_end as false for long reads
                meta.paired_end = false
                // Set default correctedness
                if ( !params.corrected ) {
                    meta.corrected = false
                } else {
                    meta.corrected = true
                }
                // Set assembly group
                if ( !params.coassembly ) {
                    meta.assembly_group = sampleID
                } else {
                    meta.assembly_group = "allReads"
                }
                // Set bin group
                if ( !params.cobinning ) {
                    meta.bin_group = sampleID
                } else {
                    meta.bin_group = "allReads"
                }
                
                def reads_path = reads.toAbsolutePath().toString()
                return [ meta, [ reads_path ]] 
            }

        nanopore_barcodes = channel.fromPath ( "${nanopore_path}*/*.{fastq,fq}.gz" )
            .map { reads ->
                def meta = [:]
                // Set barcode folder name as meta.id
                def barcode = reads.getParent().getName()
                meta.id = barcode

                return [ meta, reads ]
            }
            .groupTuple()

        CONCATENATE_ONT_BARCODES ( 
            nanopore_barcodes,
            "${launchDir}/${nanopore_path}"
            )

        // Replace the symlinks of the concatenated reads with the originals
        move_concatenated_input = CONCATENATE_ONT_BARCODES.out.count()
        REPLACE_CONCATENATED_SYMLINKS_ONT (
            move_concatenated_input,
            "${launchDir}/${nanopore_path}"
        )

        CONCATENATE_ONT_BARCODES.out.fastq.view()

        nanopore_concatenated_barcodes = ( CONCATENATE_ONT_BARCODES.out.fastq )
            .map { meta, reads ->
                def sampleID = reads.getBaseName(2)
                def meta_new = meta + [id: sampleID]
                meta_new.sequencer = "ONT"
                // Set paired_end as false for long reads
                meta_new.paired_end = false
                // Set default correctedness
                if ( !params.corrected ) {
                    meta_new.corrected = false
                } else {
                    meta_new.corrected = true
                }
                // Set assembly group
                if ( !params.coassembly ) {
                    meta_new.assembly_group = sampleID
                } else {
                    meta_new.assembly_group = "allReads"
                }
                // Set bin group
                if ( !params.cobinning ) {
                    meta_new.bin_group = sampleID
                } else {
                    meta_new.bin_group = "allReads"
                }
                
                def reads_path = "${launchDir}/${nanopore_path}${meta_new.id}.fastq.gz"
                return [ meta_new, [ reads_path ]] 
            }

        nanopore_reads = nanopore_reads.mix ( nanopore_fastqs, nanopore_concatenated_barcodes )

        //nanopore_sample_count = nanopore_barcodes.count()
        //write_csv_input = write_csv_input.mix ( nanopore_reads.combine ( nanopore_sample_count ) )

        write_csv_input = write_csv_input.mix ( nanopore_reads )         
    }

    //Create channel of concatenated raw PacBio reads
    pacbio_fastqs = channel.empty()
    pacbio_bams = channel.empty()
    pacbio_reads = channel.empty()
    if ( params.pacbio ) {
        def pacbio_path = params.pacbio.endsWith("/") ? params.pacbio : params.pacbio + "/"
        pacbio_fastqs = channel.fromPath ( "${pacbio_path}*.{fastq,fq}.gz" )
            .map { reads ->
                def meta = [:]
                def sampleID = reads.getBaseName(2).replaceAll(/.hifi_reads$/, '')
                meta.id = sampleID
                meta.sequencer = "PacBio"
                // Set paired_end as false for long reads
                meta.paired_end = false
                // Set default correctedness
                if ( !params.corrected ) {
                    meta.corrected = false
                } else {
                    meta.corrected = true
                }
                // Set assembly group
                if ( !params.coassembly ) {
                    meta.assembly_group = sampleID
                } else {
                    meta.assembly_group = "allReads"
                }
                // Set bin group
                if ( !params.cobinning ) {
                    meta.bin_group = sampleID
                } else {
                    meta.bin_group = "allReads"
                }
                
                def reads_path = reads.toAbsolutePath().toString()
                return [ meta, [ reads_path ]] 
            }

        pacbio_bams = channel.fromPath ( "${pacbio_path}*.bam" )
            .map { reads ->
                def meta = [:]
                def sampleID = reads.getBaseName(1).replaceAll(/.hifi_reads$/, '')
                meta.id = sampleID
                
                return [ meta, reads ] 
            }

        // To avoid replacing fastq.gz and bam with same sample names
        CONFIRM_SAMPLE_ID_PB (
            pacbio_bams,
            pacbio_path,
            "fastq.gz"
        )

        // Convert bams to fastq.gz. Generate bam.pbi if needed
        pbtk_bam2fastq_input = CONFIRM_SAMPLE_ID_PB.out
            .map { meta, reads ->
                def sampleID = reads.getBaseName(1)
                def meta_new = meta + [id: sampleID]
                meta_new.sequencer = "PacBio"
                // Set paired_end as false for long reads
                meta_new.paired_end = false
                // Set default correctedness
                if ( !params.corrected ) {
                    meta_new.corrected = false
                } else {
                    meta_new.corrected = true
                }
                // Set assembly group
                if ( !params.coassembly ) {
                    meta_new.assembly_group = sampleID
                } else {
                    meta_new.assembly_group = "allReads"
                }
                // Set bin group
                if ( !params.cobinning ) {
                    meta_new.bin_group = sampleID
                } else {
                    meta_new.bin_group = "allReads"
                }
                return [ meta_new, reads ]
            }

        PBTK_BAM2FASTQ (
            pbtk_bam2fastq_input,
            "${launchDir}/${pacbio_path}"
        )

        // Replace the symlinks of the generated fastq files with the originals
        move_fastq_pb_input = PBTK_BAM2FASTQ.out.count()
        REPLACE_FASTQ_SYMLINKS_PB (
            move_fastq_pb_input,
            "${launchDir}/${pacbio_path}"
        )

        pacbio_converted_fastqs = PBTK_BAM2FASTQ.out.setup_reads_csv
            .map { meta, path ->
                return [meta, [ path ]]
            }

        //pacbio_sample_count = pacbio_reads.count()
        //write_csv_input = write_csv_input.mix ( pacbio_reads.combine ( pacbio_sample_count ) )
        pacbio_reads = pacbio_reads.mix ( pacbio_fastqs, pacbio_converted_fastqs )
        write_csv_input = write_csv_input.mix ( pacbio_reads )  
    }

    write_csv_input.view()
    

    // Write csv file of input files
    CREATE_TMP_CSV ()

    // This is a channel of input, adding each file in no particular order
    WRITE_SAMPLES_CSV (
        write_csv_input,
        CREATE_TMP_CSV.out
    )

    // Counting the output of WRITE_SAMPLES_CSV ensures the csv is completed before moving
    remove_tmp_input = WRITE_SAMPLES_CSV.out.count()
    REMOVE_TMP_CSV (
        remove_tmp_input
    )

    short_config_count = short_reads.count()
    short_config_paired = false
    if ( params.illumina ) {
        short_config_paired = short_reads
            .map { meta, _reads ->
                if ( !meta.paired_end ) {
                    return "false"
                } else {
                    return "true"
                }
            }
            .unique()
    }
    nanopore_config_count = nanopore_reads.count()
    pacbio_config_count = pacbio_reads.count()
    config_corrected = params.corrected
    use_gpu = params.use_gpu

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
        config_corrected,
        use_gpu
    )

    workflow.onComplete = {
        // any workflow property can be used here
        println "<project> setup complete"
        println "Please "
        println "Command line: $workflow.commandLine"
    }

    workflow.onError = {
        println "Error: something went wrong, please see nextflow.log for details"
    }

}

