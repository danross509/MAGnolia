#!/usr/bin/env nextflow

process METABAT2 {
    tag "${meta.id}-${meta.assembler}"
    label 'process_medium'

    container "community.wave.seqera.io/library/metabat2:15c68d548f9e9b8f"
    conda "bioconda::metabat2=2.15"

    input:
        tuple val(meta), path(assembly), path(depth)


    output:
        tuple val(meta), path("*[!lowDepth|tooShort|unbinned].fa")           , optional:true, emit: bins
        tuple val(meta), path("*.tooShort.fa")                               , optional:true, emit: tooshort
        tuple val(meta), path("*.lowDepth.fa")                               , optional:true, emit: lowdepth
        tuple val(meta), path("*.unbinned.fa")                               , optional:true, emit: unbinned
        tuple val(meta), path("*.tsv")                                       , optional:true, emit: membership


    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.assembler}_MetaBAT2"
    def args = task.ext.args ?: ''

    """
    metabat2 \
    -i $assembly \
    -a $depth \
    -t $task.cpus \
    $args \
    --saveCls \
    --unbinned \
    -o ${prefix}

    mv ${prefix} ${prefix}.tsv
    """
    
}