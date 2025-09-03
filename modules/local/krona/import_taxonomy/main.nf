#!/usr/bin/env nextflow

process KRONA_K2_IMPORT_TAXONOMY {
    tag "$meta.id"
    label 'process_low'

    container ""
    conda "${moduleDir}/environment.yml"

    publishDir "${launchDir}/KRAKEN2/${file_type}/${meta.id}", mode: 'symlink'

    input:
        tuple val(meta), path(report)
        val file_type
        val placeholder

    output:
        tuple val(meta), path("*_krona.html")                        , emit: kronagram, optional: true


    script:

    """
    #if [[ -f \$CONDA_PREFIX/opt/krona/taxonomy/taxonomy.tab ]]; then
        ktImportTaxonomy \
        -m 3 -t 5 \
        $report \
        -o ${meta.id}_${file_type}_krona.html
    #fi
    """
/*    if [[ -f \$CONDA_PREFIX/opt/krona/taxonomy/placeholder ]]; then
        rm \$CONDA_PREFIX/opt/krona/taxonomy/placeholder
        working_directory=\$PWD
        cd \$CONDA_PREFIX/opt/krona/taxonomy/
        wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
        ktUpdateTaxonomy.sh --only-build
        cd \$working_directory
    fi

    count=0
    until [[ -f \$CONDA_PREFIX/opt/krona/taxonomy/taxonomy.tab ]]; do
        sleep 1
        if [[ \$count == 600 ]]; then
            echo "Skipping Krona visualization, ktUpdateTaxonomy.sh timed out (10min)"
            break
        fi
        ((count++))
    done*/
}