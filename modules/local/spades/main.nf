process SPADES {
    tag "$meta.id"
    label 'process_high'

    container "community.wave.seqera.io/library/spades:4.1.0--77799c52e1d1054a"
    conda 'bioconda::spades=4.1.0'

    input:
        tuple val(meta), path(illumina), path(pacbio), path(nanopore)
        val run_metaspades
        path yml
        path hmm

    output:
        tuple val(meta), path('metaspades/*_scaffolds.fa')          , optional:true, emit: scaffolds
        tuple val(meta), path('metaspades/*_assembly.fa')           , optional:true, emit: contigs
        tuple val(meta), path('metaspades/*_transcripts.fa')        , optional:true, emit: transcripts
        tuple val(meta), path('metaspades/*_gene_clusters.fa')      , optional:true, emit: gene_clusters
        tuple val(meta), path('metaspades/*_assembly_graph.gfa')    , optional:true, emit: gfa
        tuple val(meta), path('metaspades/*_warnings.log')          , optional:true, emit: warnings
        tuple val(meta), path('metaspades/*_spades.log')            , emit: log
        //path  "versions.yml"                                      , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def short_reads = illumina ? ( meta.single_end ? "-s $illumina" : "-1 ${illumina[0]} -2 ${illumina[1]}" ) : ""
    def pacbio_reads = pacbio ? "--pacbio $pacbio" : ""
    def nanopore_reads = nanopore ? "--nanopore $nanopore" : ""
    def custom_hmms = hmm ? "--custom-hmms $hmm" : ""
    def reads = yml ? "--dataset $yml" : "$short_reads $pacbio_reads $nanopore_reads"
    def metaspades = run_metaspades ? "--meta" : ""
    def mem = task.memory.toGiga()

    """
    spades.py \\
        --threads $task.cpus \\
        --memory $mem \\
        $custom_hmms \\
        $reads \\
        $metaspades \\
        $args \\
        -o metaspades
        
    mv metaspades/spades.log metaspades/${prefix}_spades.log

    if [ -f metaspades/scaffolds.fasta ]; then
        mv metaspades/scaffolds.fasta metaspades/${prefix}_scaffolds.fa
        #gzip -n metaspades/${prefix}_scaffolds.fa
    fi
    if [ -f metaspades/contigs.fasta ]; then
        mv metaspades/contigs.fasta metaspades/${prefix}_assembly.fa
        #gzip -n metaspades/${prefix}_assembly.fa
    fi
    if [ -f metaspades/transcripts.fasta ]; then
        mv metaspades/transcripts.fasta metaspades/${prefix}_transcripts.fa
        #gzip -n metaspades/${prefix}_transcripts.fa
    fi
    if [ -f metaspades/assembly_graph_with_scaffolds.gfa ]; then
        mv metaspades/assembly_graph_with_scaffolds.gfa metaspades/${prefix}_assembly_graph.gfa
        #gzip -n metaspades/${prefix}_assembly_graph.gfa
    fi

    if [ -f metaspades/gene_clusters.fasta ]; then
        mv metaspades/gene_clusters.fasta metaspades/${prefix}_gene_clusters.fa
        #gzip -n metaspades/${prefix}_gene_clusters.fa
    fi

    if [ -f metaspades/warnings.log ]; then
        mv metaspades/warnings.log metaspades/${prefix}_warnings.log
    fi

    """
}