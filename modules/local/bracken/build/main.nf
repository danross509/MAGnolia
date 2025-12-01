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
        val bracken_mm

    output:
        val true     , emit: bracken_built


    script:
    def replacement_py = bracken_mm ? "bracken-build_fixedMM" : "bracken-build_fixed"

    """
    if [[ -f \$CONDA_PREFIX/bin/bracken-build ]]; then
        rm \$CONDA_PREFIX/bin/bracken-build
    fi 

    cp ${moduleDir}/${replacement_py} \$CONDA_PREFIX/bin/bracken-build

    bracken-build \
    -d $database \
    -t $task.cpus \
    -k $kmer_len \
    -l $read_len

    """
}
