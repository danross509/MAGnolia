#!/usr/bin/env nextflow

process COVERM_CONTIGS {
    tag "${meta.assembler}-${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container ""

    //publishDir "${launchDir}/BINNING/${meta.id}/${meta.assembler}-McDevol", mode: 'symlink'

    input:
    tuple val(meta), path(reads), path(bams)

    output:
    //tuple val(meta), path("${prefix}/comebin_res_bins/*.fa.gz"), emit: bins
    //path "versions.yml"                                        , emit: versions

    script:
    def multi_split = meta.cobinning ? "--multi-split" : ""
    def bam_to_sam = meta.paired_end ? "samtools view -h -f 0x2 -o sam_files/\${basename}.sam \$file" : "samtools view -h -o sam_files/\${basename}.sam \$file"

    """
    coverm contig \
        --coupled fastqs/* \
        --
        -t 6 \
        -m mean relative_abundance covered_fraction \
        -o output_coverm.tsv

    rm -r sam_files
    """
}