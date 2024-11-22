#!/usr/bin/env nextflow

include {fastqc as fastqc_in} from '../modules/local/fastqc/fromPairs/main.nf'
include {fastqc as fastqc_out} from '../modules/local/fastqc/fromPairs/main.nf'
include {multiqc as multiqc_in} from '../modules/local/multiqc/main.nf'
include {multiqc as multiqc_out} from '../modules/local/multiqc/main.nf'
include {trimmomatic} from '../modules/local/trimmomatic/fromPairs/main.nf'
//include {fastp} from '../modules/local/fastp/main.nf'
include {build_index} from '../modules/local/bowtie2/build_index/main.nf'
include {filter_contaminants} from '../modules/local/bowtie2/filter_contaminants/main.nf'
include {concatenate_reads as concatenate_reads_1} from '../modules/local/scripts/concatenate_reads/main.nf'
include {concatenate_reads as concatenate_reads_2} from '../modules/local/scripts/concatenate_reads/main.nf'

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

workflow QC {
    take:
    short_reads
    phiX
    bowtie2_sensitivity
    paired_reads
    coassembly

    main:
    // Create pre-QC fastqc report for input reads;
    // Collect pre-QC fastqc output;
    // Summarize with multiqc
    fastqc_in(short_reads, "pre-QC")
    all_fastqc_in = fastqc_in.out.collect() //.ifEmpty([])
    multiqc_in(all_fastqc_in, "pre-QC")

    // Remove duplicate reads for processing
    //DEDUPLICATE()

    // Trim adapters and poor quality reads
    trimmomatic(short_reads)
    //fastp(short_reads)

    // Build index for contaminant genome(s);
    // Remove contaminant reads
    build_index(phiX)

    filter_contaminants(
        trimmomatic.out,
        build_index.out,
        bowtie2_sensitivity
        )

    clean_reads_ch = filter_contaminants.out.R1.join(filter_contaminants.out.R2, remainder: true)
        .map { meta, R1, R2 ->
            if (R2 != null ) {
                tuple(meta, [R1, R2])
            } else if (R2 == null) {
                tuple(meta, [R1])
            }
        }

    // Create post-QC fastqc report for input reads;
    // Collect post-QC fastqc output;
    // Summarize with multiqc
    fastqc_out(clean_reads_ch, "post-QC")
    all_fastqc_out = fastqc_out.out.collect()
    multiqc_out(all_fastqc_out, "post-QC")

    if (coassembly) { 
        clean_reads_1 = clean_reads_ch
            .map{ meta, reads ->
                def cMeta = [:]
                //println reads.size()
                cMeta.id = "allReads"
                cMeta.sense = "R1"
                if (reads.size() == 1){
                    cMeta.paired_end = false
                } else if (reads.size() == 2){
                    cMeta.paired_end = true
                }
                tuple(cMeta, reads[0])
            }
            .groupTuple()

        concatenate_reads_1(
            clean_reads_1
        )
        
        if (paired_reads) {
            clean_reads_2 = clean_reads_ch
                .map{ meta, reads ->
                    if (reads.size() == 2){
                        def cMeta = [:]
                        cMeta.id = "allReads"
                        cMeta.sense = "R2"
                        cMeta.paired_end = true
                        tuple(cMeta, reads[1])
                    } else { 
                        exit 1, "ERROR: Paired-read co-assembly contains invalid number of samples"
                     }
                }
                .groupTuple()

            concatenate_reads_2(
                clean_reads_2
            )

            reads_qc_out = concatenate_reads_1.out.concat(concatenate_reads_2.out)
        } else {
            reads_qc_out = concatenate_reads_1.out
        }

    } else {
        reads_qc_out = clean_reads_ch
    }

    //println clean_reads_2 == null // always false
    //println clean_reads_2.size() == null // always false
    
    //filter_contaminants.out.R2.count().view()
    //println filter_contaminants.out.size() // always 2
    clean_reads_1.view()
    /*concatenate_reads_1.out.view()
    if (paired_reads){
        clean_reads_2.view()
        concatenate_reads_2.out.view()
    }*/

    //filter_contaminants.out.view()
    //clean_reads_ch.view()
    //clean_reads_1.map{meta, reads -> println reads.size()}

    reads_qc_out.view()

println "QC timestamp"

    emit:
    reads_qc_out

}