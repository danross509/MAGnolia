#!/usr/bin/env nextflow

process MINIMAP2_INDEX {
    tag "$meta.id"
    label 'process_low'

    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"
    conda "bioconda::minimap2=2.28"

    //publishDir "${params.resultsDir}/QC/${meta.id}/minimap2_alignment/", mode: 'symlink'

    input:

    tuple val(meta), path(fasta), val(preset_input)


    output:
    tuple val(meta), path("*.mmi"), val(preset_input)       , emit: index
    path "versions.yml"                                     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def preset = task.ext.preset_input ?: "map-ont" // Minimap2 default value

    """
    minimap2 -x $preset \
    -t $task.cpus \
    -d ${meta.id}-${preset}.mmi \
    $args \
    $fasta 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}