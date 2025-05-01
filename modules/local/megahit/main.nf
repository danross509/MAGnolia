#!/usr/bin/env nextflow

process MEGAHIT {
    tag "$meta.id"

    container "community.wave.seqera.io/library/megahit:1.2.9--23234b8da1e27898"
    conda "bioconda::megahit=1.2.9"

    publishDir "${launchDir}/Assembly/${meta.id}/", mode: 'symlink'

    input:
        tuple val(meta), path(reads)
        val mh_preset
        val bigThreads
        val bigMem

    output:
        tuple val(meta), path("megahit/*.contigs.fa"), emit: final_contigs
        tuple val(meta), path("megahit/intermediate_contigs/k*.addi.fa")
        tuple val(meta), path("megahit/intermediate_contigs/k*.contigs.fa")
        tuple val(meta), path("megahit/intermediate_contigs/k*.final.contigs.fa")
        tuple val(meta), path("megahit/intermediate_contigs/k*.local.fa")
        tuple val(meta), path('megahit/log')
        path "megahit/options.json"

    script:

    threads = bigThreads

    if (meta.paired_end) {

        //reads_1 = reads[0]
        //reads_2 = reads[1]
        //{--presets meta} not working
        //{--tmp-dir tmp} not working
        //gzip files after? pigz?

        """
        megahit \
        -1 ${reads[0]} -2 ${reads[1]} \
        -o megahit \
        -t $threads \
        -m ${bigMem}000000000 \
        --presets $mh_preset \
        --verbose \
        --continue

        mv ./megahit/final.contigs.fa ./megahit/${meta.id}_final.contigs.fa
        """
    } else if (!meta.paired_end) {

        """
        megahit \
        -r ${reads[0]} \
        -o megahit \
        -t $threads \
        -m ${bigMem}000000000 \
        --presets $mh_preset \
        --verbose \
        --continue

        mv ./megahit/final.contigs.fa ./megahit/${meta.id}_final.contigs.fa
        """
    }
}