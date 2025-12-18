#!/usr/bin/env nextflow

process GTDB_DB_DOWNLOAD {

    conda "bioconda::gtdbtk=${params.gtdb_version}"
    container ""

    input:
    val db_url

    output:
    path "release*", emit: db
    path "versions.yml", emit: versions

    script:
    """
    wget $db_url
    tar -xvzf *.tar.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(echo \$(gtdbtk --version -v 2>&1) | sed "s/gtdbtk: version //; s/ Copyright.*//")
    END_VERSIONS
    """
}