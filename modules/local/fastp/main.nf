#!/usr/bin/env nextflow

process FASTP {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/fastp:0.23.4--f8cefc1e5f7a782e"
    conda "bioconda::fastp=0.24.3"

    publishDir "${launchDir}/QC/${meta.id}/fastp", mode: 'symlink'

    input:
        tuple val(meta), path(reads_fastq)
        val adapter_sequences
        val auto_adapter_detection
        val qualified_quality_phred
        val unqualified_percent_limit
        val disable_length_filtering
        val length_required
        val deduplication

    output:
        tuple val(meta), path("${meta.id}_trimmed_{1,2}.fastq.gz"), emit: reads_PE, optional: true
        tuple val(meta), path("${meta.id}_trimmed_se.fastq.gz"), emit: reads_SE, optional: true
        path "*_fastp.html"
        path "*_fastp.json"

    script:
    fq_1_paired = "${meta.id}_trimmed_1.fastq.gz"
    fq_2_paired = "${meta.id}_trimmed_2.fastq.gz"
    fq_unpaired = "${meta.id}_trimmed_se.fastq.gz"

    def adapter_R1 = adapter_sequences ? "--adapter_sequence ${adapter_sequences[0]}" : ""
    def adapter_R2 = adapter_sequences && adapter_sequences.size() == 2 ? "--adapter_sequence_r2 ${adapter_sequences[1]}" : ""
    def detect_pe_adapters = auto_adapter_detection ? "--detect_adapter_for_pe" : ""
    def minimum_quality = "${qualified_quality_phred}"
    def unqualified_limit = "${unqualified_percent_limit}"
    def disable_length_limit = disable_length_filtering ? "--disable_length_filtering" : ""
    def minimum_length = length_required ? "--length_required ${length_required}" : ""
    def dedup = deduplication ? "--dedup" : ""

    if (meta.paired_end) {

        """
        fastp  \
        -i ${reads_fastq[0]} \
        -I ${reads_fastq[1]} \
        -o $fq_1_paired \
        -O $fq_2_paired \
        -h ${meta.id}_fastp.html \
        -j ${meta.id}_fastp.json \
        --unpaired1 $fq_unpaired \
        --unpaired2 $fq_unpaired \
        $adapter_R1 \
        $adapter_R2 \
        $detect_pe_adapters \
        --qualified_quality_phred $minimum_quality \
        --unqualified_percent_limit $unqualified_limit \
        $disable_length_limit \
        $minimum_length \
        $dedup
        """
    } else if (!meta.paired_end) {

        """
        fastp  \
        -i ${reads_fastq[0]} \
        -o $fq_unpaired \
        -h ${meta.id}_fastp.html \
        -j ${meta.id}_fastp.json \
        $adapter_R1 \
        --qualified_quality_phred $minimum_quality \
        --unqualified_percent_limit $unqualified_limit \
        $disable_length_limit \
        $minimum_length \
        $dedup
        """
    }
}