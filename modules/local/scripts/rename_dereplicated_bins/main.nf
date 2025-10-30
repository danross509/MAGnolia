#!/usr/bin/env nextflow

process RENAME_DEREPLICATED_BINS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::sed=4.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("*_dRep.*.fa"),         optional:true, emit: dereplicated_bins

    script:
    """
    for bin in $bins; do
        if [[ -f \$bin ]]; then
            filename=\${bin##*/}                # sampleR2_MEGAHIT_SemiBin2_Refined.1.fa
            basename=\${filename%%.*}           # sampleR2_MEGAHIT_SemiBin2_Refined
            remainder=\${filename#*.}           #                                  .1.fa
            newbasename=\${basename%_Refined}   # sampleR2_MEGAHIT_SemiBin2
            mv \$bin ./\${newbasename}_dRep.\${remainder}
        fi
    done
    """
}