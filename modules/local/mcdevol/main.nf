#!/usr/bin/env nextflow

process MCDEVOL {
    tag "${meta.assembler}-${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container ""

    publishDir "${launchDir}/BINNING/${meta.id}/${meta.assembler}-McDevol", mode: 'symlink'

    input:
    tuple val(meta), path(assembly), path(depths)
    //tuple val(meta), path(assembly), path(bams)

    output:
    //tuple val(meta), path("${prefix}/comebin_res_bins/*.fa.gz"), emit: bins
    //path "versions.yml"                                        , emit: versions

    script:
    def multi_split = meta.cobinning ? "--multi-split" : ""
    //def bam_to_sam = meta.paired_end ? "samtools view -h -f 0x2 -o sam_files/\${basename}.sam \$file" : "samtools view -h -o sam_files/\${basename}.sam \$file"
    def unzip_contigs = assembly.getExtension() == "gz" ? "gunzip -c $assembly > contigs_unzipped.fa" : "cp $assembly contigs_unzipped.fa"
    def unzip_depths = depths.getExtension() == "gz" ? "gunzip -c $depths > abundance_unzipped.tsv" : "cp $depths abundance_unzipped.tsv"

    //-i sam_files \\
    """
    #if [[ ! -d sam_files ]]; then
    #    mkdir sam_files
    #fi

    #for file in \$bams; do
    #    basename=\${file%.sorted.bam}
    #    \$bam_to_sam
    #    samtools view -c sam_files/\${basename}.sam
    #done

    $unzip_contigs
    $unzip_depths

    mcdevol.py \\
        -a abundance_unzipped.tsv \\
        -c contigs_unzipped.fa \\
        -o output \\
        --abundformat metabat \\
        -n $task.cpus \\
        $multi_split

    #rm -r sam_files
    rm contigs_unzipped.fa
    rm abundance_unzipped.tsv
    """
}