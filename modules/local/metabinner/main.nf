#!/usr/bin/env nextflow

process METABINNER {
    tag "${meta.id}-${meta.assembler}"
    label 'process_medium'

    container ""
    conda "${moduleDir}/environment.yml"

    publishDir "${params.resultsDir}/BINNING/${meta.id}/${meta.assembler}-MetaBinner", mode: 'symlink'

    input:
        tuple val(meta), path(assembly), path(depths)
        val scale
        val min_contig_len
        val kmer_len


    output:
        tuple val(meta), path("bins/*.fa")  , optional:true, emit: bins

    script:
    def dataset_scale = scale ?: "large"
    def unzip_contigs = assembly.getExtension() == "gz" ? "gunzip -c $assembly > contigs_unzipped.fa" : "cp $assembly contigs_unzipped.fa"
    def unzip_depths = depths.getExtension() == "gz" ? "gunzip -c $depths > abundance_unzipped.tsv" : "cp $depths abundance_unzipped.tsv"

    """
    metabinner_path=\$(dirname \$(which run_metabinner.sh))

    $unzip_contigs
    $unzip_depths

    cat abundance_unzipped.tsv | awk '{if (\$2>${min_contig_len}) print \$0 }' | cut -f -1,4- > coverage_profile.tsv

    gen_kmer.py contigs_unzipped.fa $min_contig_len $kmer_len

    kmer_profile=\${PWD}/contigs_unzipped_kmer_${kmer_len}_f${min_contig_len}.csv

    Filter_tooshort.py contigs_unzipped.fa $min_contig_len

    run_metabinner.sh \\
        -a \${PWD}/contigs_unzipped_${min_contig_len}.fa \\
        -o \${PWD}/output \\
        -t $task.cpus \\
        -d \${PWD}/coverage_profile.tsv \\
        -k \$kmer_profile \\
        -p \$metabinner_path \\
        -s $dataset_scale

    rm contigs_unzipped.fa
    rm abundance_unzipped.tsv

    bin_path=\${PWD}/output/metabinner_res/ensemble_res/greedy_cont_weight_3_mincomp_50.0_maxcont_15.0_bins/ensemble_3logtrans/addrefined2and3comps/greedy_cont_weight_3_mincomp_50.0_maxcont_15.0_bins
    mkdir bins
    for file in \${bin_path}/*; do
        if [[ -f \$file ]]; then
            filename=\${file##*_}
            num=\${filename%.*}
            mv \$file ./bins/${meta.id}_${meta.assembler}_MetaBinner.\${num}.fa
        fi
    done

    #gzip bins/*
    
    """
    
}