#!/usr/bin/env nextflow

process MINIMAP2_ASSEMBLY_MAPPED_SORTED {
    label 'process_MEDIUM'
    tag "${assembly_meta.id}_${reads_meta.id}"

    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"
    conda "bioconda::minimap2=2.28 bioconda::samtools=1.21"

    //publishDir "${launchDir}/CLEAN_READS/nanopore", mode: 'symlink'

    input:

    tuple val(assembly_meta), path(assembly), path(assembly_index), val(preset_input), val(reads_meta), path(reads)

    output:
    tuple val(assembly_meta), path(assembly), path("${assembly_meta.id}_${reads_meta.id}_mapped.sorted.bam"), path("${assembly_meta.id}_${reads_meta.id}_mapped.sorted.bam.bai"), emit: mappings
    //tuple val(assembly_meta), val(reads_meta), path("*.minimap2.log"), emit: log
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def name = "${assembly_meta.id}_${reads_meta.id}"
    def preset = task.ext.preset_input ?: "map-ont"
    def query = reads.size() == 1 ? "${reads[0]}" : "${reads[0]} ${reads[1]}"

    """
    minimap2 -a \\
    -x $preset \\
    -t $task.cpus \\
    $args \\
    $assembly_index \\
    $query \\
    2> ${name}.minimap2.log | \
    samtools view -@ $task.cpus -b -F 4 | \
    samtools sort -@ $task.cpus -o ${name}_mapped.sorted.bam
    samtools index ${name}_mapped.sorted.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version)
    END_VERSIONS
    """

}