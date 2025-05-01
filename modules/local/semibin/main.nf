#!/usr/bin/env nextflow

process SEMIBIN {

    container "community.wave.seqera.io/library/pip_semibin:b6a41dbb4d1296c7"
    conda "${moduleDir}/environment.yml"

    "${launchDir}/Binning/${meta.id}/${meta.assembler}-${semibin_version}", mode: 'symlink'

    input:
        tuple val(assembly_meta), path(assembly), path(bams), path(environment)
        val use_semibin1

    output:
        tuple 
        //path "versions.yml", emit: versions

    script:

    def semibin_version = use_semibin1 ? "SemiBin1": "SemiBin2"
    def mode = 
    def command = 
    def sequencing_type = 
    def reads_bam = 
    def name = "${assembly_meta.id}_${reads_meta.id}"
    def args = task.ext.args ?: ''


    if (mode == "single_sample_binning") {
        """
        $semibin_version single_easy_bin \
        --environment human_gut \
        --sequencing-type long_read \ 
        -i S1.fa \
        -b S1.sorted.bam \
        -o output
        """
    } else if (mode == "coassembly_binning") { // Co-assembly
        """
        $semibin_version single_easy_bin \
        --sequencing-type long_read \ 
        -i contig.fa \
        -b S1.sorted.bam S2.sorted.bam S3.sorted.bam \
        -o co-assembly_output
        """
    } else if (mode == "multisample_binning") { // Co-binning
        """
        $semibin_version multi_easy_bin \
        --sequencing-type long_read \
        -i concatenated.fa \
        -b S1.sorted.bam S2.sorted.bam S3.sorted.bam S4.sorted.bam S5.sorted.bam \
        -o multi_output
        """
    }


}