#!/usr/bin/env nextflow

process PBTK_BAM2FASTQ {
    tag "$meta.id"

    container ""
    conda "bioconda::pbtk=3.5.0"

    publishDir "${launchDir}/${sample_prep_directory}", mode: 'symlink'

    input:
        tuple val(meta), path(bam)
        val sample_prep_directory


    output:
        tuple val(meta), path("${meta.id}.fastq.gz")                                       , emit: fastq, optional: true
        tuple val(meta), val("${launchDir}/${sample_prep_directory}/${meta.id}.fastq.gz")  , emit: setup_reads_csv, optional: true
        // versions

    script:

    """
    if [[ -f ${launchDir}/${sample_prep_directory}/${meta.id}.bam.pbi ]]; then
        ln -s ${launchDir}/${sample_prep_directory}/${meta.id}.bam.pbi ./${meta.id}.bam.pbi
    else
        pbindex $bam
    fi

    bam2fastq -o $meta.id $bam
    """
}