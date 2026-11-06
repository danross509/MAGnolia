process MYLOASM {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/67/67f705c8fa26ae809d73e5ee2c872c3087fe8abaad12e133bdba8ba4f8cdc553/data'
        : 'community.wave.seqera.io/library/myloasm_findutils_gzip:154a360acd2dca07'}"

    input:
    tuple val(meta), path(reads), val(use_hifi)     // Added 'use_hifi' DRoss 10-20-26

    output:
    tuple val(meta), path("${prefix}"), emit: results
    tuple val(meta), path("${prefix}/${prefix}.assembly_primary.fa.gz"), emit: contigs
    tuple val(meta), path("${prefix}/${prefix}.final_contig_graph.gfa.gz"), emit: gfa
    tuple val(meta), path("${prefix}/alternate_assemblies/${prefix}.assembly_alternate.fa.gz"), emit: contigs_alt
    tuple val(meta), path("${prefix}/alternate_assemblies/${prefix}.duplicated_contigs.fa.gz"), emit: contigs_dup
    tuple val(meta), path("${prefix}/3-mapping/${prefix}.map_to_unitigs.paf.gz"), emit: mapping
    tuple val(meta), path("${prefix}/*.log"), emit: log
    tuple val("${task.process}"), val('myloasm'), eval("myloasm --version | sed 's/.* //'"), emit: versions_myloasm, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    output_folder = task.ext.output_folder ?: "myloasm-out" // Added DRoss 10-28-26
    def kmc = params.myloasm_use_kmc ? "--kmc" : ''         // Added DRoss 10-20-26
    def hifi = use_hifi ? "--hifi" : ''                     // Added DRoss 10-20-26
    """
    myloasm \\
        ${reads} \\
        -o ${output_folder} \\
        -t ${task.cpus} \\
        ${kmc} \\
        ${hifi} \\
        ${args}

    cd ${output_folder}
    mylotools report --output mylotools-report

    cd ../
    find "${output_folder}" \\( -name "*.fa" -o -name "*.gfa" -o -name "*.edges" -o -name "*.paf.gz" \\) -type f | while read file; do
        dirname=\$(dirname "\${file}")
        basename=\$(basename "\${file}")
        mv "\$file" "\${dirname}/${prefix}.\${basename}"
        gzip \${dirname}/${prefix}.\${basename}
    done

    mv ${output_folder}/mylotools-report ./
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}/alternate_assemblies
    mkdir -p ${prefix}/3-mapping
    echo "" | gzip > ${prefix}/${prefix}.assembly_primary.fa.gz
    echo "" | gzip > ${prefix}/${prefix}.final_contig_graph.gfa.gz
    echo "" | gzip > ${prefix}/alternate_assemblies/${prefix}.assembly_alternate.fa.gz
    echo "" | gzip > ${prefix}/alternate_assemblies/${prefix}.duplicated_contigs.fa.gz
    echo "" | gzip > ${prefix}/3-mapping/${prefix}.map_to_unitigs.paf.gz
    touch ${prefix}/myloasm_1.log
    """
}