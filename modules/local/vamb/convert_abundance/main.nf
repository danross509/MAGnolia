#!/usr/bin/env nextflow

process VAMB_CONVERT_ABUNDANCE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::pandas"
    container ""

    input:
    tuple val(meta), path(assembly), path(abundance_tsv), path(bams, stageAs: "bams/*"), path(taxonomy)

    output:
    tuple val(meta), path(assembly), path("converted_abundance.tsv"), path(bams), path(taxonomy)

    script:
    def unzip_depths = abundance_tsv.getExtension() == "gz" ? "gunzip -c $abundance_tsv > abundance_unzipped.tsv" : "cp $abundance_tsv abundance_unzipped.tsv"

    """
    $unzip_depths

    vamb_convert_abundance.py -i abundance_unzipped.tsv -o converted_abundance.tsv
    """
}