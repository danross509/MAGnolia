#!/usr/bin/env nextflow

process BRACKEN_SUMMARIZE_ABUNDANCE {
    tag "$file_type"
    label 'process_single'

    container ""
    conda "conda-forge::python=3.13 conda-forge::pandas=2.3.0 conda-forge::seaborn=0.13 conda-forge::matplotlib=3.10"

    input:
        val abundances
        val file_type
        val level
        val zero

    output:
        path("*_abundance_summary.csv")     , emit: summary


    
    script:

    """
    summarize_bracken_abundance.py -n $file_type -l $level -0 $zero -r $abundances
    """
}