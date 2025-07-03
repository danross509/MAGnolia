process SPADES {
    tag "$meta.id"
    label 'process_high'

    container "community.wave.seqera.io/library/spades:4.1.0--77799c52e1d1054a"
    conda 'bioconda::spades=4.1.0'

    publishDir "${launchDir}/Assembly/${meta.id}/", mode: 'symlink'

    input:
        tuple val(meta), path(illumina), path(pacbio), path(nanopore)
        val run_metaspades
        val bigThreads
        val bigMem
        path yml
        path hmm

    output:
        tuple val(meta), path('metaspades/*.scaffolds.fa.gz')    , optional:true, emit: scaffolds
        tuple val(meta), path('metaspades/*.assembly.fa.gz')      , optional:true, emit: contigs
        tuple val(meta), path('metaspades/*.transcripts.fa.gz')  , optional:true, emit: transcripts
        tuple val(meta), path('metaspades/*.gene_clusters.fa.gz'), optional:true, emit: gene_clusters
        tuple val(meta), path('metaspades/*.assembly_graph.gfa.gz')    , optional:true, emit: gfa
        tuple val(meta), path('metaspades/*.warnings.log')       , optional:true, emit: warnings
        tuple val(meta), path('metaspades/*.spades.log')         , emit: log
        //path  "versions.yml"                          , emit: versions

    //when:
    //task.ext.when == null || task.ext.when

    script:
    //def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def short_reads = illumina ? ( meta.single_end ? "-s $illumina" : "-1 ${illumina[0]} -2 ${illumina[1]}" ) : ""
    def pacbio_reads = pacbio ? "--pacbio $pacbio" : ""
    def nanopore_reads = nanopore ? "--nanopore $nanopore" : ""
    def custom_hmms = hmm ? "--custom-hmms $hmm" : ""
    def reads = yml ? "--dataset $yml" : "$short_reads $pacbio_reads $nanopore_reads"
    def metaspades = run_metaspades ? "--meta" : ""
    """
    spades.py \\
        --threads $bigThreads \\
        --memory $bigMem \\
        $custom_hmms \\
        $reads \\
        $metaspades \\
        -o metaspades
    mv metaspades/spades.log metaspades/${prefix}.spades.log

    if [ -f metaspades/scaffolds.fasta ]; then
        mv metaspades/scaffolds.fasta metaspades/${prefix}.scaffolds.fa
        gzip -n metaspades/${prefix}.scaffolds.fa
    fi
    if [ -f metaspades/contigs.fasta ]; then
        mv metaspades/contigs.fasta metaspades/${prefix}.assembly.fa
        gzip -n metaspades/${prefix}.assembly.fa
    fi
    if [ -f metaspades/transcripts.fasta ]; then
        mv metaspades/transcripts.fasta metaspades/${prefix}.transcripts.fa
        gzip -n metaspades/${prefix}.transcripts.fa
    fi
    if [ -f metaspades/assembly_graph_with_scaffolds.gfa ]; then
        mv metaspades/assembly_graph_with_scaffolds.gfa metaspades/${prefix}.assembly.gfa
        gzip -n metaspades/${prefix}.assembly.gfa
    fi

    if [ -f metaspades/gene_clusters.fasta ]; then
        mv metaspades/gene_clusters.fasta metaspades/${prefix}.gene_clusters.fa
        gzip -n metaspades/${prefix}.gene_clusters.fa
    fi

    if [ -f metaspades/warnings.log ]; then
        mv metaspades/warnings.log metaspades/${prefix}.warnings.log
    fi

    """
    /*
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p')
    END_VERSIONS
    */
    /*stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxmem = task.memory.toGiga()
    def illumina_reads = illumina ? ( meta.single_end ? "-s $illumina" : "-1 ${illumina[0]} -2 ${illumina[1]}" ) : ""
    def pacbio_reads = pacbio ? "--pacbio $pacbio" : ""
    def nanopore_reads = nanopore ? "--nanopore $nanopore" : ""
    def custom_hmms = hmm ? "--custom-hmms $hmm" : ""
    def reads = yml ? "--dataset $yml" : "$illumina_reads $pacbio_reads $nanopore_reads"
    """
    echo "" | gzip > ${prefix}.scaffolds.fa.gz
    echo "" | gzip > ${prefix}.contigs.fa.gz
    echo "" | gzip > ${prefix}.transcripts.fa.gz
    echo "" | gzip > ${prefix}.gene_clusters.fa.gz
    echo "" | gzip > ${prefix}.assembly.gfa.gz
    touch ${prefix}.spades.log
    touch ${prefix}.warnings.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py --version 2>&1 | sed -n 's/^.*SPAdes genome assembler v//p')
    END_VERSIONS
    """
    */
}