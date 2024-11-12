#!/usr/bin/env nextflow

/*
I have created a separate folder to run tests from, located in the same parent directory as this git repository
I am testing it with the following data:
wget https://zenodo.org/record/3992790/files/test_reads.tar.gz
tar -xzf test_reads.tar.gz

I am developing it using conda environments. To run it you need to simply have a conda environment with nextflow installed
From the MAG_Pipeline folder, run the following command:
    nextflow run pipeline/QC/main.nf -params-file pipeline/QC/params.json -resume

You will need to ensure that the path to your data in parans.json is accurate
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
    // If paired-read filenames are identified by R1/R2, the R will be removed for this pipeline
    short_reads = Channel.fromFilePairs(params.short_reads, checkIfExists:true)
        .map { sample, reads ->
            def sampleID = sample.replaceAll(/_R$/, '')
            return [sampleID, reads]
        }
 
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