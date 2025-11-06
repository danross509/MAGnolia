#!/usr/bin/env nextflow

process MAXBIN2 {
    tag "${meta.assembler}-${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/maxbin2:d4f776ef82533cd3"
    conda "bioconda::maxbin2=2.2.7"

    //publishDir "${params.resultsDir}/Binning/maxbin2/", mode: 'symlink'

    input:
        tuple val(meta), path(assembly), path(reads), path(abund)


    output:
        tuple val(meta), path("*.fasta")     , optional:true, emit: bins
        tuple val(meta), path("*.summary")      , emit: summary
        tuple val(meta), path("*.abundance")    , emit: abundance   , optional: true
        tuple val(meta), path("*.log.gz")       , emit: log
        tuple val(meta), path("*.marker.gz")    , emit: marker_counts
        tuple val(meta), path("*.noclass")   , emit: unbinned_fasta
        tuple val(meta), path("*.tooshort.gz")  , emit: tooshort_fasta
        tuple val(meta), path("*_bin.tar.gz")   , emit: marker_bins , optional: true
        tuple val(meta), path("*_gene.tar.gz")  , emit: marker_genes, optional: true

    script:

    def prefix = task.ext.prefix ?: "${meta.id}_${meta.assembler}_MaxBin2"
    def associate_files = ""
    if ( abund instanceof List ) {
        associate_files = "-abund ${abund[0]}"
        for (i in 2..abund.size()) { associate_files += " -abund$i ${abund[i-1]}" }
    } else {
        associate_files = "-abund $abund"
    }

    //Look into maxbin2 options

    """
    run_MaxBin.pl \
    -contig $assembly \
    $associate_files \
    -thread ${task.cpus} \
    -out $prefix

    gzip *.tooshort *log *.marker
    #gzip *.fasta *.noclass *.tooshort *log *.marker 
    """
    
}