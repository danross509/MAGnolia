#!/usr/bin/env nextflow

/*
From the pipeline directory, use the command:
    nextflow run QC/main.nf -params-file QC/params.json -resume
*/

include {fastqc as fastqc_in} from '../modules/local/fastqc/fromPairs/main.nf'
include {fastqc as fastqc_out} from '../modules/local/fastqc/fromPairs/main.nf'
include {multiqc as multiqc_in} from '../modules/local/multiqc/main.nf'
include {multiqc as multiqc_out} from '../modules/local/multiqc/main.nf'
include {trimmomatic} from '../modules/local/trimmomatic/fromPairs/main.nf'
include {build_index} from '../modules/local/bowtie2/build_index/main.nf'
include {filter_contaminants} from '../modules/local/bowtie2/filter_contaminants/main.nf'

/*
process deduplicate {

    container 
    conda 

    //publishDir '${projectDir}/QC/clumpify', mode: 'symlink'

    input:
        path 

    output:
        path

}
*/

workflow {

    // Specify the output directory
    // projectDir = Channel.fromPath(params.projectDir, type: 'dir')

    // Create input channel from a text file listing input file paths
    short_reads = Channel.fromFilePairs(params.short_reads, checkIfExists:true)
 
    // Create pre-QC fastqc report for input reads;
    // Collect pre-QC fastqc output;
    // Summarize with multiqc
    fastqc_in(short_reads, "pre-QC")
    all_fastqc_in = fastqc_in.out.collect().ifEmpty([])
    multiqc_in(all_fastqc_in, "pre-QC")

    // Remove duplicate reads for processing
    //DEDUPLICATE()

    // Trim adapters and poor quality reads
    trimmomatic(short_reads, params.pairing)

    // Build index for contaminant genome(s);
    // Remove contaminant reads
    build_index(params.phiX)

    filter_contaminants(
        trimmomatic.out,
        build_index.out,
        params.bowtie2_sensitivity
        )

    // Create post-QC fastqc report for input reads;
    // Collect post-QC fastqc output;
    // Summarize with multiqc
    fastqc_out(filter_contaminants.out, "post-QC")
    all_fastqc_out = fastqc_out.out.collect().ifEmpty([])
    multiqc_out(all_fastqc_out, "post-QC")

}