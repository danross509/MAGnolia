#!/usr/bin/env nextflow

process FLYE {
    tag "$meta.id"
    label 'process_high'

    container "community.wave.seqera.io/library/flye:2.9.5--d577924c8416ccd8"
    conda "bioconda::flye=2.9.5"

    publishDir "${launchDir}/Assembly/${meta.id}/metaFlye", mode: 'symlink'

    input:
        tuple val(meta), path(reads), val(read_type)
        val run_metaflye
        val genome_size
        val bigThreads

    output:
        tuple val(meta), path("*_assembly.fa"), emit: final_contigs
        tuple val(meta), path("*_assembly_graph.gfa"), path("*_assembly_graph.gv"), emit: assembly_graph
        tuple val(meta), path("*_assembly_info.txt"), emit: info

    script:
    def reads_input = "--${read_type} ${reads}"
    def metaflye = run_metaflye ? "--meta" : ""
    def genome = genome_size ?: ""

    """
    flye \
    $reads_input \
    --threads $bigThreads \
    $metaflye \
    $genome \
    --out-dir ./

    mv assembly.fasta ${meta.id}_assembly.fa
    mv assembly_graph.gfa ${meta.id}_assembly_graph.gfa
    mv assembly_graph.gv ${meta.id}_assembly_graph.gv
    mv assembly_info.txt ${meta.id}_assembly_info.txt
    """
}