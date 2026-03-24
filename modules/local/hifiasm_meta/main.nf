#!/usr/bin/env nextflow

process HIFIASM_META {
    tag "$meta.id"
    label 'process_high'

    container ""
    conda "${moduleDir}/environment.yml"

    input:
        tuple val(meta), path(reads)
        val skip_hmbin
        val enable_rs
        val force_rs
        val rs_threshold

    output:
        tuple val(meta), path("*.p_ctg.gfa"),   emit: primary_contig_graph
        tuple val(meta), path("*.a_ctg.gfa"),   emit: alternate_contig_graph
        tuple val(meta), path("*.r_utg*.gfa"),  emit: raw_unitig_graph
        tuple val(meta), path("*.p_utg*.gfa"),  emit: cleaned_unitig_graph
        tuple val(meta), path("*.rescue.fa"),   emit: circle_contigs, optional: true
        tuple val(meta), path("*.bins.tsv"),    emit: bins, optional: true 
        path "${meta.id}_assembly.log", emit: log
        //path "", emit: verison

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_assembly"
    def binning = skip_hmbin ? "--no-binning" : ""
    def read_selection = enable_rs ? "-S" : ""
    def force_selection = force_rs ? "--force-rs" : ""
    def rs_quantile = rs_threshold ?: ""

    def ont_ec = meta.sequencer == "ONT" ? "hifiasm --ont -o $prefix -e -t$task.cpus ${reads} &> ${meta.id}_ha.log" : ""
    def ont_assembly = meta.sequencer == "ONT" ? "--use-ha-bin" : ""

    """
    $ont_ec

    hifiasm_meta \
    -t $task.cpus \
    $ont_assembly \
    -o $prefix \
    $read_selection \
    $force_selection \
    $rs_quantile \
    $binning \
    $args \
    ${reads}
    &> ${meta.id}_assembly.log
    """
}