#!/usr/bin/env nextflow

process LRBINNER {
    tag "${meta.assembler}-${meta.id}"
    //label 'process_high'

    errorStrategy = { 
        def status = task.exitStatus
        
        // Log for debugging
        println "Task ${task.name} failed with exit status: ${status}"
        
        // Handle null/undefined status
        if (status == null || status < 0) {
            println "Exit status is null or negative, retrying..."
            return 'retry'
        }
        
        // Handle known retriable codes
        if (status in (130..145) || status == 104) {
            println "Exit status ${status} is retriable"
            return 'retry'
        }
        
        return 'retry'
    }

    conda "${moduleDir}/environment.yml"
    container ""

    publishDir "${launchDir}/BINNING/${meta.id}/${meta.assembler}-LRBinner", mode: 'symlink'

    input:
    tuple val(meta), path(reads), path(assembly)
    //tuple val(meta), path(assembly), path(bams)

    output:
    //tuple val(meta), path("${prefix}/comebin_res_bins/*.fa.gz"), emit: bins
    //path "versions.yml"                                        , emit: versions

    script:
    //def bam_to_sam = meta.paired_end ? "samtools view -h -f 0x2 -o sam_files/\${basename}.sam \$file" : "samtools view -h -o sam_files/\${basename}.sam \$file"
    def unzip_contigs = assembly.getExtension() == "gz" ? "gunzip -c $assembly > contigs_unzipped.fa" : "cp $assembly contigs_unzipped.fa"
    def unzip_reads = reads.getExtension() == "gz" ? "gunzip -c $reads > reads_unzipped.fastq" : "cp $reads reads_unzipped.fastq"
    def cuda = params.use_gpu ? "--cuda" : "" // Make this detect if gpu is available
    def resume = task.attempt > 1 ? "--resume" : ""

    //-i sam_files \\
    """
    $unzip_reads
    $unzip_contigs

    lrbinner.py contigs \\
        --contigs contigs_unzipped.fa \\
        --reads-path reads_unzipped.fastq \\
        -k 4 \\
        -t $task.cpus \\
        $cuda \\
        --output ./output \\
        $resume

    rm contigs_unzipped.fa
    rm reads_unzipped.fastq
    """
}