#!/usr/bin/env nextflow

process METABAT2 {
    tag "${meta.assembler}-${meta.id}"
    label 'process_medium'

    container "community.wave.seqera.io/library/metabat2:15c68d548f9e9b8f"
    conda "bioconda::metabat2=2.15"

    publishDir "${launchDir}/BINNING/${meta.id}/${meta.assembler}-metabat2", mode: 'symlink'

    input:
        tuple val(meta), path(assembly), path(depth)


    output:
        tuple val(meta), path("*[!lowDepth|tooShort|unbinned].fa")           , optional:true, emit: bins
        tuple val(meta), path("*.tooShort.fa")                               , optional:true, emit: tooshort
        tuple val(meta), path("*.lowDepth.fa")                               , optional:true, emit: lowdepth
        tuple val(meta), path("*.unbinned.fa")                               , optional:true, emit: unbinned
        tuple val(meta), path("*.tsv")                                       , optional:true, emit: membership


    script:

    depth_file = depth.getBaseName()
    def prefix = task.ext.prefix ?: "bin"
    //Look into metabat options

    """
    gzip -d -f $depth

    metabat2 \
    -i $assembly \
    -a $depth_file \
    -t $task.cpus \
    --saveCls \
    --unbinned \
    -o ${prefix}

    #gzip -cn ${prefix} > ${prefix}.tsv.gz
    #find . -name "*.fa" -type f | xargs -t -n 1 bgzip
    """
    
}