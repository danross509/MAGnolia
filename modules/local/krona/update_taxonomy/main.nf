#!/usr/bin/env nextflow

process KRONA_K2_UPDATE_TAXONOMY {
    label 'process_low'

    container ""
    conda "${moduleDir}/../import_taxonomy/environment.yml"

    input:

    output:
        path("empty_file.txt"), emit: placeholder


    script:

    """
    if [[ ! -f \$CONDA_PREFIX/opt/krona/taxonomy/taxonomy.tab ]]; then
        working_directory=\$PWD
        cd \$CONDA_PREFIX/opt/krona/taxonomy/
        wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
        ktUpdateTaxonomy.sh --only-build
        cd \$working_directory
    fi

    touch empty_file.txt
    """
}