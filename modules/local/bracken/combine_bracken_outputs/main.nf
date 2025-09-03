#!/usr/bin/env nextflow

process COMBINE_BRACKEN_OUTPUTS {
    tag "$file_type"
    label 'process_single'

    container ""
    conda "${moduleDir}/../abundance_estimation/environment.yml"

    publishDir "${launchDir}/KRAKEN2/${file_type}", mode: 'symlink'

    input:
        val abundances
        val file_type

    output:
        path("*.tsv")  , emit: combined_output


    script:
    def name = "${file_type}_combined_bracken_output.tsv"
    def files = "${abundances}".strip("[]").replaceAll(/,/, '')

    """
    combine_bracken_outputs.py --files $files --output $name
    """
}
