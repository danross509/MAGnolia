process EGGNOGMAPPER_METAGENOME {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container ""

    input:
    tuple val(meta), path(fasta)
    path(eggnog_db)
    val(search_mode)
    val(gene_prediction)
    val(overlaps)
    val(overlap_tol)

    output:
    tuple val(meta), path("*.emapper.annotations")   , emit: annotations
    tuple val(meta), path("*.emapper.seed_orthologs"), emit: orthologs, optional: true
    tuple val(meta), path("*.emapper.hits")          , emit: hits     , optional: true
    //tuple val("${task.process}"), val('eggnog-mapper'), eval("emapper.py --version 2>&1 | grep -o 'emapper-[0-9]\\+\\.[0-9]\\+\\.[0-9]\\+' | sed 's/emapper-//'"), topic: versions, emit: versions_eggnogmapper

    script:
    def genepred = gene_prediction ? "--genepred ${gene_prediction}" : ''
    def args = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    emapper.py \\
        -m $search_mode \\
        --itype metagenome \\
        -i $fasta \\
        --data_dir ${eggnog_db} \\
        $genepred \\
        --allow_overlaps $overlaps \\
        --overlap_tol $overlap_tol \\
        --cpu ${task.cpus} \\
        -o ${prefix} \\
        $args
    """
}