#!/usr/bin/env nextflow

process GATB_MINIA_PIPELINE {
    tag "$meta.id"
    label 'process_high'

    container "community.wave.seqera.io/library/minia:3.2.6--92bae1756baab1ef", "community.wave.seqera.io/library/gatb:1.4.2--cf37d08b7005497e"
    conda "${moduleDir}/environment.yml"

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*_assembly.fa"), emit: final_contigs
        tuple val(meta), path("*_assembly.gfa"), emit: assembly_graph
        tuple val(meta), path("*_unitigs.fa"), emit: unitigs, optional: true
        //tuple val(meta), path("final_assembly_depths.txt"), emit: depth

    script:
    def read_files = meta.paired_end ? "-1 ${reads[0]} -2 ${reads[1]}" : "-s ${reads[0]}"
    def args = task.ext.args ?: ''

    """
    gatb_minia_pipeline.py $read_files --no-scaffolding --no-error-correction $args --minia \$CONDA_PREFIX/bin/minia

    # Replace assembly_final.contigs.fa symlink with contigs file
    target=\$(readlink -e assembly_final.contigs.fa)
    prefix=\${target%%.*}
    unitigs=\${prefix}.unitigs.fa
    kmer_size=\${prefix#*_k}

    sed -i '' assembly_final.contigs.fa
    rm \$target

    gatb_bcalm_convertToGFA.py assembly_final.contigs.fa ${meta.id}_assembly.gfa \$kmer_size

    mv assembly_final.contigs.fa ${meta.id}_assembly.fa
    mv \${unitigs} ${meta.id}_unitigs.fa
    """
}