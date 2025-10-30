#!/usr/bin/env nextflow

process HIFIASM_CREATE_BIN_FILES {
    tag "$meta.id"
    label 'process_single'

    container ""
    conda ""

    publishDir "${launchDir}/BINNING/${meta.id}/${meta.assembler}-hmBin", mode: 'symlink'

    input:
        tuple val(meta), path(all_contigs)
        tuple val(meta), path(circle_contigs)
        tuple val(meta), path(bins)

    output:
        tuple val(meta), path("circular_MAGs/*.fa"),    emit: circular_mags
        tuple val(meta), path("bins/*.fa"),             emit: bins

    script:
    def prefix = "${meta.id}_${meta.assembler}_hmBin"
    //def command = contig_graph.getExtension() == "gz" ? "zcat $contig_graph | awk '/^S/{print \">\"\$2\"\\n\"\$3}' - >${prefix}.fa" : "cat $contig_graph | awk '/^S/{print \">\"\$2\"\\n\"\$3}' - >${prefix}.fa"

    """
    if [[ ! -d circular_MAGs ]]; then
        mkdir circular_MAGs
    fi

    hifiasm_split_circular_contigs.py $circle_contigs circular_MAGs $prefix

    if [[ ! -d bins ]]; then
        mkdir bins
    fi

    hifiasm_extract_bins.py --contigs $all_contigs --bins $bins --prefix $prefix
    """
}