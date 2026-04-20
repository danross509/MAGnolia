#!/usr/bin/env nextflow

process MULTIQC {
    label 'process_single'

    container "community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b"
    conda "bioconda::multiqc=1.33"

    input:
        path all_fastqc
        val runName

    output:
        path "multiqc_report*.html"
        path "multiqc_data*"

    script:
    def args = task.ext.args ?: ''

    """
    multiqc $all_fastqc \
    $args

    mv multiqc_report.html multiqc_report_${runName}.html
    mv multiqc_data multiqc_data_${runName}
    """ 
}