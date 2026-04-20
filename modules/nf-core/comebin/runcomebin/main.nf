process COMEBIN_RUNCOMEBIN {
    tag "${meta.id}-${meta.assembler}"
    label 'process_high'
    label 'process_gpu'

    conda params.use_gpu ? "${moduleDir}/environment_gpu.yml" : "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/comebin:1.0.4--hdfd78af_0':
        'biocontainers/comebin:1.0.4--hdfd78af_0' }"

    input:
    tuple val(meta), path(assembly), path(bam, stageAs: "bam/*")
    val input_batch_size

    output:
    tuple val(meta), path("${prefix}/comebin_res_bins/*.fa")   , emit: bins
    tuple val(meta), path("${prefix}/*.comebin_res.tsv")         , emit: tsv
    tuple val(meta), path("${prefix}/*.comebin.log")             , emit: log
    tuple val(meta), path("${prefix}/*.embeddings.tsv")          , emit: embeddings
    tuple val(meta), path("${prefix}/*.covembeddings.tsv")       , emit: covembeddings
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args               = task.ext.args ?: ''
    prefix                 = task.ext.prefix ?: "${meta.id}"
    def clean_assembly     = assembly.toString() - ~/\.gz$/
    def decompress_contigs = assembly.toString() == clean_assembly ? "" : "gunzip -q -f $assembly"
    def cleanup            = decompress_contigs ? "rm ${clean_assembly}" : ""
    """
    ${decompress_contigs}

    nContigs=\$(seqkit seq -m 1000 ${clean_assembly} | grep -c "^>")

    if (( \$nContigs <= ${input_batch_size} )); then
        batch=\$((\${nContigs}/2))
    else 
        batch=${input_batch_size}
    fi

    echo "There are \$nContigs contigs"
    echo "Using -b \$batch"

    run_comebin.sh \\
        -t ${task.cpus} \\
        -a ${clean_assembly} \\
        -p bam/ \\
        -o . \\
        -b \$batch \\
        $args

    mv comebin_res ${prefix}
    
    for bin in ${prefix}/*; do
        if [[ -f \$bin ]]; then
            bin_name=\${bin##*/}
            mv \$bin ${prefix}/${meta.id}_${meta.assembler}_COMEbin.\${bin_name}
            if [[ "\$bin_name" == *.gz ]]; then
                gunzip ${prefix}/${meta.id}_${meta.assembler}_COMEbin.\${bin_name}
            fi
        fi
    done

    #find ${prefix}/comebin_res_bins/*.fa -exec gzip {} \\;

    ${cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: \$(run_comebin.sh | sed -n 2p | grep -o -E "[0-9]+(\\.[0-9]+)+")
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/comebin_res_bins

    echo "" | gzip > ${prefix}/comebin_res_bins/1.fa.gz
    echo "" | gzip > ${prefix}/comebin_res_bins/2.fa.gz

    touch ${prefix}/comebin_res.tsv
    touch ${prefix}/comebin.log
    touch ${prefix}/embeddings.tsv
    touch ${prefix}/covembeddings.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        comebin: \$(run_comebin.sh | sed -n 2p | grep -o -E "[0-9]+(\\.[0-9]+)+")
    END_VERSIONS
    """
}