#!/usr/bin/env nextflow

process MINIMAP2_ASSEMBLY_ALIGNMENT {
    tag "${meta.id}_${sampleID}"
    label 'process_high'

    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"
    conda "bioconda::minimap2=2.28 bioconda::samtools=1.21"

    //publishDir "${launchDir}/CLEAN_READS/nanopore", mode: 'symlink'

    input:
    tuple val(meta), path(assembly), path(assembly_index), val(preset_input), val(sampleID), path(reads)
    val unmapped

    output:
    tuple val(meta), path(assembly), path("${meta.id}_${sampleID}.sorted.bam"), path("${meta.id}_${sampleID}.sorted.bam.bai"), emit: mappings
    //tuple val(meta), val(sampleID), path("*.minimap2.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def name = "${meta.id}_${sampleID}"
    def preset = task.ext.preset_input ?: "map-ont"
    def query = reads.size() == 1 ? "${reads[0]}" : "${reads[0]} ${reads[1]}"
    def remove_unmapped = unmapped ? "samtools view -@ $task.cpus -b -F 4 |" : ""

    """
    minimap2 --eqx -a \\
    -x $preset \\
    -t $task.cpus \\
    $args \\
    $assembly_index \\
    $query \\
    2> ${name}.minimap2.log | \
    $remove_unmapped \
    samtools sort -@ $task.cpus -o ${name}.sorted.bam
    samtools index ${name}.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version)
    END_VERSIONS
    """

}
