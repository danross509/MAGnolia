#!/usr/bin/env nextflow

process BRACKEN_BUILD {
    tag ""
    label 'process_high'

    container ""
    conda "bioconda::kraken2=2.17 bioconda::bracken=3.1"

    //publishDir "${params.resultsDir}/KRAKEN2/${meta.id}/${file_type}", mode: 'symlink'

    input:
        val database
        val kmer_len
        val read_len

    output:
        val true     , emit: bracken_built


    script:

    """
    if [[ -f \$CONDA_PREFIX/bin/bracken-build ]]; then
        rm \$CONDA_PREFIX/bin/bracken-build
    fi 

    cp ${moduleDir}/bracken-build_fixed \$CONDA_PREFIX/bin/bracken-build

    bracken-build \
    -d $database \
    -t $task.cpus \
    -k $kmer_len \
    -l $read_len 
    """
}
