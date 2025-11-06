#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_single'

    container "community.wave.seqera.io/library/multiqc:1.25.1--dc1968330462e945"
    conda "bioconda::multiqc=1.25.1"

    publishDir "${params.resultsDir}/QC/MultiQC/${step}", mode: 'symlink'

    input:
        path all_fastqc
        val step

    output:
        path "multiqc_report.html"
        path "multiqc_data"

    """
    multiqc $all_fastqc
    """ 
}