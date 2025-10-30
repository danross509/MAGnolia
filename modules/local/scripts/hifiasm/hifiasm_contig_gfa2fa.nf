#!/usr/bin/env nextflow

process HIFIASM_CONTIG_GFA2FA {
    tag "$meta.id"
    label 'process_single'

    container ""
    conda ""

    publishDir "${launchDir}/ASSEMBLY/${meta.id}/Hifiasm_meta", mode: 'symlink'

    input:
        tuple val(meta), path(contig_graph)

    output:
        tuple val(meta), path("*_assembly.fa"), emit: final_contigs

    script:
    def prefix = "${meta.id}_assembly"
    def command = contig_graph.getExtension() == "gz" ? "zcat $contig_graph | awk '/^S/{print \">\"\$2\"\\n\"\$3}' - >${prefix}.fa" : "cat $contig_graph | awk '/^S/{print \">\"\$2\"\\n\"\$3}' - >${prefix}.fa"

    """
    $command

    #gzip ${prefix}.fa
    """
}