#!/usr/bin/env nextflow

process METAMDBG {
    tag "$meta.id"
    label 'process_high'

    container ""
    conda "${moduleDir}/environment.yml"

    publishDir "${launchDir}/Assembly/${meta.id}/metaMDBG", mode: 'symlink'

    input:
        tuple val(meta), path(reads)
        val generate_assembly_graph
        val assembly_graph_k
        val assembly_graph_contigs
        val assembly_graph_reads
        val threads

    output:
        tuple val(meta), path("output/*_assembly.fa.gz")       , emit: final_contigs
        tuple val(meta), path("output/*assemblyGraph*.gfa")    , emit: assembly_graph, optional: true
        //tuple val(meta), path("*_assembly_info.txt")    , emit: info, optional: true

    script:
    def reads_input = meta.sequencer == "ONT" ? "--in-ont ${reads}" : meta.sequencer == "PacBio" ? "--in-hifi ${reads}" : ""
    def gfa_k = assembly_graph_k ?: "21"
    def contigpath = assembly_graph_contigs ? "--contigpath" : ""
    def readpath = assembly_graph_reads ? "--readpath" : ""
    def assembly_graph = generate_assembly_graph ? "metaMDBG gfa --assembly-dir output --k ${gfa_k} ${contigpath} ${readpath} --threads ${threads}" : ""

    """
    metaMDBG asm \
    $reads_input \
    --threads $threads \
    --out-dir output

    $assembly_graph

    mv output/contigs.fasta.gz output/${meta.id}_assembly.fa.gz
    """
}