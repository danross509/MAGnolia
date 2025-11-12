#!/usr/bin/env nextflow

process HMBIN_VERIFY_F2C2B {
    tag "$meta.id"
    label 'process_single'

    container ""
    conda ""

    //publishDir "${params.resultsDir}/ASSEMBLY/${meta.id}/Hifiasm_meta", mode: 'symlink'

    input:
        tuple val(meta), path(f2c2b)

    output:
        tuple val(meta), path("*_verified.tsv"), optional: true

    script:
    def basename = f2c2b.getBaseName()
    def output_file = "${basename}_verified.tsv"
    
    """
    verify_hmbin_f2c2b.py $f2c2b $output_file
    """
}