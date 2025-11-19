#!/usr/bin/env nextflow

process SEMIBIN2 {
    tag "${meta.assembler}-${meta.id}"
    label 'process_high'

    container "community.wave.seqera.io/library/pip_semibin:b6a41dbb4d1296c7"
    conda "${moduleDir}/environment.yml"

    publishDir "${params.resultsDir}/BINNING/${meta.id}/${meta.assembler}-semibin2", mode: 'symlink'

    input:
        tuple val(meta), path(assembly), path(bams)
        val envr
        val use_semibin1

    output:
        tuple val(meta), path("{,co-assembly_,multi_}output/{,output_}bins/*.fa*")                  , optional:true, emit: bins
        //tuple val(meta), path("{,co-assembly_}output/output_recluster_bins/*.fa.gz")              , optional:true, emit: recluster_bins
        tuple val(meta), path("{,co-assembly_,multi_}output/{samples/}*_data_cov.csv")              , optional:true, emit: coverage
        tuple val(meta), path("{,co-assembly_,multi_}output/{samples/*/}*contig_bins.tsv")          , optional:true, emit: contig_bin_assignment
        tuple val(meta), path("{,co-assembly_,multi_}output/{samples/*/}*data_split.csv")           , optional:true, emit: training_data_split
        tuple val(meta), path("{,co-assembly_,multi_}output/{samples/*/}*data.csv")                 , optional:true, emit: training_data
        tuple val(meta), path("{,co-assembly_,multi_}output/{samples/*/}*model.pt")                 , optional:true, emit: model
        tuple val(meta), path("{,co-assembly_,multi_}output/{samples/*/}*recluster_bins_info.tsv")  , optional:true, emit: bin_stats
        tuple val(meta), path("{,co-assembly_,multi_}output/*SemiBinRun.log")                       , optional:true, emit: log
        //path "versions.yml", emit: versions

    script:

    def semibin_version = use_semibin1 ? "SemiBin1": "SemiBin2"
    def mode = meta.cobinning ? "multi_easy_bin" : "single_easy_bin"
    def environment = meta.coassembly ? "" : meta.cobinning ? "" : "--environment $envr"
    def sequencing_type = meta.sequencer == "Illumina" ? "" : "--sequencing-type long_read"
    def output = meta.coassembly ? "co-assembly_output" : meta.cobinning ? "multi_output" : "output"
    def args = task.ext.args ?: ''

    """
    $semibin_version $mode \\
    $environment \\
    $sequencing_type \\
    -i $assembly \\
    -b $bams \\
    -o $output

    if [[ -d multi_output ]]; then 
        for bin in multi_output/bins/*; do
            if [[ -f \$bin ]]; then
                bin_name=\${bin##*_}
                mv \$bin multi_output/bins/${meta.id}_${meta.assembler}_SemiBin2.\${bin_name}
                if [[ "\$bin_name" == *.gz ]]; then
                    gunzip multi_output/bins/${meta.id}_${meta.assembler}_SemiBin2.\${bin_name}
                fi
            fi
        done
        for folder in multi_output/samples/*; do
            if [[ -d \$folder ]]; then
                sample_name=\${folder##*/}
                for file in \${folder}/*; do
                    if [[ -f \$file ]]; then
                        filename=\${file##*/}
                        mv \$file multi_output/samples/\${sample_name}/${meta.id}_\${sample_name}_\${filename}
                    fi
                done
            fi
        done
        mv multi_output/SemiBinRun.log multi_output/${meta.id}_SemiBinRun.log

    elif [[ -d co-assembly_output ]]; then 
        for bin in co-assembly_output/output_bins/*; do
            if [[ -f \$bin ]]; then
                bin_name=\${bin##*_}
                mv \$bin co-assembly_output/output_bins/${meta.id}_${meta.assembler}_SemiBin2.\${bin_name}
                if [[ "\$bin_name" == *.gz ]]; then
                    gunzip co-assembly_output/output_bins/${meta.id}_${meta.assembler}_SemiBin2.\${bin_name}
                fi
            fi
        done
        for file in co-assembly_output/*; do
            if [[ -f \$file ]]; then
                filename=\${file##*/}
                if [[ \$filename != ${meta.id}_* ]]; then
                    mv \$file co-assembly_output/${meta.id}_\${filename}
                fi
            fi
        done

    elif [[ -d output ]]; then
        for bin in output/output_bins/*; do
            if [[ -f \$bin ]]; then
                bin_name=\${bin##*_}
                mv \$bin output/output_bins/${meta.id}_${meta.assembler}_SemiBin2.\${bin_name}
                if [[ "\$bin_name" == *.gz ]]; then
                    gunzip output/output_bins/${meta.id}_${meta.assembler}_SemiBin2.\${bin_name}
                fi
            fi
        done
        for file in output/*; do
            if [[ -f \$file ]]; then
                filename=\${file##*/}
                if [[ \$filename != ${meta.id}_* ]]; then
                    mv \$file output/${meta.id}_\${filename}
                fi
            fi
        done
    fi
    """


}