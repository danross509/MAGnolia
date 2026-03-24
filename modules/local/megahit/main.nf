#!/usr/bin/env nextflow

process MEGAHIT {
    tag "$meta.id"
    label 'process_high'
    label 'process_gpu'

    container "community.wave.seqera.io/library/megahit:1.2.9--23234b8da1e27898"
    conda "bioconda::megahit=1.2.9"

    input:
        tuple val(meta), path(reads)
        val mh_preset

    output:
        tuple val(meta), path("megahit/*_assembly.fa"), emit: final_contigs
        tuple val(meta), val([]),  emit: assembly_graph                         // megahit does not produce a gfa, this is a placeholder 
        tuple val(meta), path("megahit/intermediate_contigs/k*.addi.fa")
        tuple val(meta), path("megahit/intermediate_contigs/k*.contigs.fa")
        tuple val(meta), path("megahit/intermediate_contigs/k*.final.contigs.fa")
        tuple val(meta), path("megahit/intermediate_contigs/k*.local.fa")
        tuple val(meta), path('megahit/log')
        path "megahit/options.json"

    script:
    
    def read_files = meta.paired_end ? "-1 ${reads[0]} -2 ${reads[1]}" : "-r ${reads[0]}"
    def use_gpu = params.use_gpu ? "--use-gpu" : ""

    def args = task.ext.args ?: ''

    """
    megahit \
    $read_files \
    -o megahit \
    -t $task.cpus \
    -m $task.memory \
    --presets $mh_preset \
    $use_gpu \
    --verbose \
    $args \
    --continue

    mv ./megahit/final.contigs.fa ./megahit/${meta.id}_assembly.fa
    """
}