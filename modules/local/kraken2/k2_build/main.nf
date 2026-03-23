#!/usr/bin/env nextflow

process K2_BUILD {
    tag ""
    label 'process_high'

    container ""
    conda "${moduleDir}/../classify/environment.yml"

    input:
        val db_dir
        val build
        val kmer_len
        val max_db_size

    output:
        val db_dir , emit: directory


    script:
    def kmer_length = kmer_len ? "--kmer-len ${kmer_len}" : ""
    def max_size = max_db_size ? "--max-db-size ${max_db_size}" : ""
    def args = task.ext.args ?: ''

    """
    k2 build \
    $build \
    --db $db_dir \
    $kmer_length \
    $max_size \
    $args
    """
}