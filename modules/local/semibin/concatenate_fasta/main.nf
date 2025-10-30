#!/usr/bin/env nextflow

process SEMIBIN_CONCATENATE_FASTA {
    tag "$meta.id"
    label 'process_low'

    container "community.wave.seqera.io/library/pip_semibin:b6a41dbb4d1296c7"
    conda "${moduleDir}/../semibin2/environment.yml"

    //publishDir "${launchDir}/BINNING/${meta.id}/${meta.assembler}-${semibin_version}", mode: 'symlink'

    input:
        tuple val(meta), path(contigs)

    output:
        tuple val(meta), path("${meta.id}_contigs.fa"), emit: concatenated_fasta
        //path "versions.yml", emit: versions

    script:

    def extension = !meta.cobinning && contigs.getExtension() == "gz" ? "fa.gz" : "fa"

    if ( meta.cobinning ) {
        """
        SemiBin2 concatenate_fasta \
        --input-fasta $contigs \
        -o ./

        mv ./concatenated.fa.gz ./${meta.id}_contigs.fa.gz
        gunzip ./${meta.id}_contigs.fa.gz
        """
    } else {
        """
        ln -s ${contigs[0]} ${meta.id}_contigs.${extension}
        """
    }



}