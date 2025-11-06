#!/usr/bin/env nextflow

process KRAKEN2_DB_DOWNLOAD {
    tag ""
    label 'process_medium'

    container ""
    conda "${moduleDir}/../classify/environment.yml"

    //publishDir "${params.resultsDir}/KRAKEN2/${meta.id}/${file_type}", mode: 'symlink'

    input:
        val db_dir
        val build
        val kmer_len
        val max_db_size

    output:
        val "${db_dir}/kraken2_db"                        , emit: directory


    script:
    def name = "${db_dir}/kraken2_db"
    def kmer_length = kmer_len ?: ""
    def max_size = max_db_size ?: ""

    """
    k2 build \
    $build \
    --db $name \
    $kmer_length \
    $max_size
    """
}
