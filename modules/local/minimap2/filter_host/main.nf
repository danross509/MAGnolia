#!/usr/bin/env nextflow

process MINIMAP2_FILTER_HOST {
    tag "$meta.id"
    label 'process_high'

    container "community.wave.seqera.io/library/minimap2:2.28--78db3d0b6e5cb797"
    conda "bioconda::minimap2=2.28 bioconda::samtools=1.21"

    //publishDir "${params.resultsDir}/QC/${meta.id}/minimap2_alignment/", mode: 'symlink'

    input:

    tuple val(meta), path(reads), val(contaminant_meta), path(contaminant_index), val(preset_input)
    val suffix

    output:
    tuple val(meta), path("${meta.id}_${suffix}.fastq.gz")          , emit: filtered, optional: true
    tuple val(meta), path("${meta.id}_${suffix}_1.fastq.gz")        , emit: R1, optional: true
    tuple val(meta), path("${meta.id}_${suffix}_2.fastq.gz")        , emit: R2, optional: true
    path "versions.yml"                                             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def name = "${meta.id}_${contaminant_meta.id}"
    def preset = task.ext.preset_input ?: "map-ont"
    def query = meta.paired_end ? "${reads[0]} ${reads[1]}" : "${reads[0]}"
    def flags = meta.paired_end ? "-f 12 -F 256" : "-f 4"
    def output = (meta.sequencer == 'ONT' || meta.sequencer == 'PacBio') ? "| gzip > ${meta.id}_${suffix}.fastq.gz" : (meta.sequencer == 'Illumina' && reads.size() == 1) ? "| gzip > ${meta.id}_${suffix}_1.fastq.gz" : "-1 >(gzip > ${meta.id}_${suffix}_1.fastq.gz) -2 >(gzip > ${meta.id}_${suffix}_2.fastq.gz) -0 /dev/null -s /dev/null -n"

    """
    minimap2 -a \\
    -x $preset \\
    -t $task.cpus \\
    $args \\
    $contaminant_index \\
    $query \\
    2> "${name}.minimap2.log" | \
    samtools view -@ $task.cpus -b $flags - | \
    samtools fastq -@ $task.cpus $output

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(samtools --version)
    END_VERSIONS
    """
}