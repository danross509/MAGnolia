#!/usr/bin/env nextflow

process MODIFY_CHECKM_FOR_DREP {
    tag "$meta.id"
    label 'process_single'

    container ""
    conda "conda-forge::pandas"

    input:
        tuple val(meta), path(bin_evaluation)

    output:
        tuple val(meta), path("${meta.id}_quality_drep.csv"), emit: modified_bin_quality
        //path "versions.yml", emit: versions

    script:

    """
    modify_checkm_for_drep.py \
    $bin_evaluation \
    ${meta.id}_quality_drep.csv
    """

}