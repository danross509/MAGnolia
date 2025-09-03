#!/usr/bin/env nextflow

process KRAKEN2 {
    tag "${meta.id}"

    container ""
    conda "${moduleDir}/environment.yml"

    publishDir "${launchDir}/KRAKEN2/${file_type}/${meta.id}", mode: 'symlink'

    input:
        tuple val(meta), path(input_files)
        path database
        val names
        val confidence
        val quick
        val mapping
        val minimizer
        val zero_counts
        val min_hits
        val file_type

    output:
        tuple val(meta), path("*.kreport")  , emit: reports
        tuple val(meta), path("*.krak2")    , emit: output
        tuple val(meta), path("*.log")      , emit: log


    script:
    def name = "${meta.id}_${file_type}"
    def paired = input_files.size() == 2 ? "--paired": ""
    def use_quick = quick ? "--quick" : ""
    def use_names = names ? "--use-names" : ""
    def memory_mapping = mapping ? "--memory-mapping" : ""
    def minimizer_data = minimizer ? "--report-minimizer-data" : ""
    def report_zero_counts = zero_counts ? "--report-zero-counts" : ""

    """
    mkdir tmp
    cp \$(readlink $input_files) tmp/
    gunzip tmp/* 

    k2 classify \
    $use_names \
    --db $database \
    $paired \
    --threads $task.cpus \
    --report ${name}.kreport \
    $use_quick \
    --output ${name}.krak2 \
    $memory_mapping \
    $minimizer_data \
    $report_zero_counts \
    --minimum-hit-groups $min_hits \
    --log ${name}.log \
    tmp/*

    rm -r tmp
    """
}


