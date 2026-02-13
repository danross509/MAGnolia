#!/usr/bin/env nextflow

process FASTPLONG {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/fastplong:0.2.2--03217e568967f163"
    conda "fastplong=0.2.2"

    publishDir "${params.resultsDir}/QC/${meta.id}/fastplong", mode: 'symlink'

    input:
        tuple val(meta), path(reads_fastq)
        val disable_quality_filter
        val n_bases
        val n_percent
        val min_phred
        val unqualified_percent
        val mean_qual_limit
        val min_length
        val max_length

    output:
        tuple val(meta), path("${meta.id}_trimmed.fastq.gz"), emit: reads, optional: true
        path "*fastplong.html"
        path "*fastplong.json"

    script:
    def fq_trimmed = "${meta.id}_trimmed.fastq.gz"
    def quality_filter = disable_quality_filter ? "--disable_quality_filtering" : ""
    def n_base_limit = n_bases ? "--n_base_limit ${n_bases}" : ""
    def n_percent_limit = n_percent ? "--n_percent_limit ${n_percent}" : ""
    def qualified_quality_phred = min_phred ? "--qualified_quality_phred ${min_phred}" : ""
    def unqualified_percent_limit = unqualified_percent ? "--unqualified_percent_limit ${unqualified_percent}" : ""
    def mean_qual = mean_qual_limit ? "--mean_qual ${mean_qual_limit}" : ""
    def length_required = min_length ? "--length_required ${min_length}" : ""
    def length_limit = max_length ? "--length_limit ${max_length}" : ""

    """
    fastplong  \
    -i ${reads_fastq} \
    -o $fq_trimmed \
    $quality_filter \
    $n_base_limit \
    $n_percent_limit \
    $qualified_quality_phred \
    $unqualified_percent_limit \
    $mean_qual \
    $length_required \
    $length_limit \
    -h ${meta.id}_fastplong.html \
    -j ${meta.id}_fastplong.json
    """

}