#!/usr/bin/env nextflow

process METAMDBG {
    tag "$meta.id"
    label 'process_high'
    // Use gpu for gfa?

    container ""
    conda "${moduleDir}/environment.yml"

    publishDir "${launchDir}/Assembly/${meta.id}/metaMDBG", mode: 'symlink'

    input:
        tuple val(meta), path(reads)
        val generate_assembly_graph
        val assembly_graph_k
        val assembly_graph_contigs
        val assembly_graph_reads

    output:
        tuple val(meta), path("output/*_assembly.fa.gz")            , emit: final_contigs
        tuple val(meta), path("output/*_assembly_graph.gfa.gz")     , emit: assembly_graph, optional: true
        //tuple val(meta), path("*_assembly_info.txt")    , emit: info, optional: true

    script:
    def reads_input = meta.sequencer == "ONT" ? "--in-ont ${reads}" : meta.sequencer == "PacBio" ? "--in-hifi ${reads}" : ""
    def contigpath = assembly_graph_contigs ? "--contigpath" : ""
    def readpath = assembly_graph_reads ? "--readpath" : ""
    def gfa_k = assembly_graph_k == "max" || assembly_graph_k == "min" || assembly_graph_k.toString().isInteger() ? assembly_graph_k : { error "Invalid assembly_graph_k value: ${assembly_graph_k}. Use 'max', 'min', or an integer." }


    """
    metaMDBG asm \
        $reads_input \
        --threads $task.cpus \
        --out-dir output

    if [[ $generate_assembly_graph ]]; then
        metaMDBG gfa \
            --assembly-dir \
            output \
            --k 0 \
            2> available_k.txt

        max_gfa_k=\$(tail -n 1 available_k.txt | sed 's/.*- \\([0-9]\\+\\).*/\\1/')
        min_gfa_k=\$(sed -n '3s/.*- \\([0-9]\\+\\).*/\\1/p' available_k.txt)

        echo "max gfa k_value = \$max_gfa_k"
        echo "min gfa k_value = \$min_gfa_k"

        if [[ $gfa_k == max ]]; then
            k_value=\$max_gfa_k
        elif [[ $gfa_k == min ]]; then
            k_value=\$min_gfa_k
        elif [[ $gfa_k -ge \$max_gfa_k ]]; then
            k_value=\$((\$max_gfa_k-1))
            echo "Specified metaMDBG k value <${gfa_k}> is greater than the maximum <\$max_gfa_k>, using <\$max_gfa_k>"
        elif [[ $gfa_k < \$min_gfa_k ]]; then
            k_value=\$min_gfa_k
            echo "Specified metaMDBG k value <${gfa_k}> is below the minimun <\$min_gfa_k>, using <\$min_gfa_k>"
        else
            k_value=$gfa_k
        fi

        metaMDBG gfa \
        --assembly-dir output \
        --k \$k_value \
        ${contigpath} \
        ${readpath} \
        --threads $task.cpus
    fi 

    mv output/contigs.fasta.gz output/${meta.id}_assembly.fa.gz
    mv output/assemblyGraph*bps.gfa output/${meta.id}_assembly_graph.gfa
    gzip output/${meta.id}_assembly_graph.gfa
    """
}