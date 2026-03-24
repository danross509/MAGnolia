#!/usr/bin/env nextflow

process QUAST_CONTIGS {
    tag "${meta.id}-${meta.assembler}"

    conda "bioconda::quast=5.3.0"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        //'https://depot.galaxyproject.org/singularity/quast:5.0.2--py37pl526hb5aa323_2' :
        //'biocontainers/quast:5.0.2--py37pl526hb5aa323_2' }"

    input:
    tuple val(meta), path(assembly)
    val min_contig
    val rna_finding
    val gene_finding
    val conserved_gene_finding
    val min_alignment
    val min_identity
    val mem_efficient
    val space_efficient
    val blast_db

    output:
    path "${meta.id}/*"                                         , emit: qc
    path "${meta.id}/${meta.id}_${meta.assembler}_report.tsv"   , emit: report
    path "versions.yml"                                         , emit: versions

    script:
    def min_contig_length = !min_contig ? "" : "-m $min_contig"
    def rna = !rna_finding ? "" : "--rna-finding"
    def genes = !gene_finding ? "" : "--gene-finding"
    def conserved_genes = !conserved_gene_finding ? "" : "--conserved-genes-finding"
    def min_alignment_length = !min_alignment ? "" : "--min-alignment $min_alignment"
    def min_identity_per = !min_identity ? "" : "--min-identity $min_identity"
    def mem_eff = !mem_efficient ? "" : "--memory-efficient"
    def space_eff = !space_efficient ? "" : "--space-efficient"
    def use_blast_db = !blast_db ? "" : "--blast-db $blast_db"
    def args = task.ext.args ?: ''

    """
    metaquast.py \\
        -t "${task.cpus}" \\
        $min_contig_length \\
        $rna \\
        $genes \\
        $conserved_genes \\
        $min_alignment_length \\
        $min_identity_per \\
        $mem_eff \\
        $space_eff \\
        $use_blast_db \\
        --max-ref-number 0 \\
        -l "${meta.id}_${meta.assembler}" \\
        "${assembly}" \\
        $args \\
        -o "${meta.id}"

    for entry in ${meta.id}/*; do
        if [[ -f \$entry ]]; then
            basename=\${entry##*/}
            if [[ \$basename != ${meta.id}_${meta.assembler}_* ]]; then
                mv \$entry ${meta.id}/${meta.id}_${meta.assembler}_\${basename}
            fi
        elif [[ -d \$entry ]]; then
            directory=\${entry##*/}
            for file in \${entry}/*; do
                if [[ -f \$file ]]; then
                    basename=\${file##*/}
                    if [[ \$basename != ${meta.id}_${meta.assembler}_* ]]; then
                        mv \$file ${meta.id}/\${directory}/${meta.id}_${meta.assembler}_\${basename}
                    fi
                fi
            done
        fi
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        metaquast: \$(metaquast.py --version | sed "s/QUAST v//; s/ (MetaQUAST mode)//")
    END_VERSIONS
    """
}