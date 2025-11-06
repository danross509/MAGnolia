#!/usr/bin/env nextflow

process BOWTIE2_ASSEMBLY_MAPPED_SORTED {
    tag "${assembly_meta.id}_${reads_meta.id}"

    container "community.wave.seqera.io/library/bowtie2:2.5.4--d51920539234bea7"
    conda "bioconda::bowtie2=2.5.4 bioconda::samtools=1.21"

    //publishDir "${params.resultsDir}/CLEAN_READS/", mode: 'symlink'

    input:
        tuple val(assembly_meta), path(assembly), path(index), val(reads_meta), path(reads)

    output:
        tuple val(assembly_meta), path(assembly), path("${assembly_meta.id}_${reads_meta.id}_mapped.sorted.bam"), path("${assembly_meta.id}_${reads_meta.id}_mapped.sorted.bam.bai"), emit: mappings
        tuple val(assembly_meta), val(reads_meta), path("*.bowtie2.log"), emit: log
        //path "versions.yml", emit: versions

    script:

    def args = task.ext.args ?: ''
    def name = "${assembly_meta.id}_${reads_meta.id}"
    //def input = params.single_end ? "-U \"${reads}\"" :  "-1 \"${reads[0]}\" -2 \"${reads[1]}\""

    if (assembly_meta.paired_end) {
        """
        INDEX=`find -L ./ -name "*.rev.1.bt2l" -o -name "*.rev.1.bt2" | sed 's/.rev.1.bt2l//' | sed 's/.rev.1.bt2//'`
        bowtie2 \\
        -p "${task.cpus}" \\
        -x \$INDEX \\
        -1 ${reads[0]} \\
        -2 ${reads[1]} \\
        2> "${name}.bowtie2.log" | \
        samtools view -@ "${task.cpus}" -bS | \
        samtools view -@ "${task.cpus}" -b -F 4 | \
        samtools sort -@ "${task.cpus}" -o "${name}_mapped.sorted.bam"
        samtools index "${name}_mapped.sorted.bam"
        """
    } else if (!assembly_meta.paired_end) {
        """
        INDEX=`find -L ./ -name "*.rev.1.bt2l" -o -name "*.rev.1.bt2" | sed 's/.rev.1.bt2l//' | sed 's/.rev.1.bt2//'`
        bowtie2 \\
        -p "${task.cpus}" \\
        -x \$INDEX \\
        -U ${reads[0]} \\
        2> "${name}.bowtie2.log" | \
        samtools view -@ "${task.cpus}" -bS | \
        samtools view -@ "${task.cpus}" -b -F 4 | \
        samtools sort -@ "${task.cpus}" -o "${name}_mapped.sorted.bam"
        samtools index "${name}_mapped.sorted.bam"
        """
    }

        /*if [ "${name}" = "${assembly_meta.id}_${assembly_meta.id}" ] ; then
            mv "${name}.bowtie2.log" "${assembly_meta.id}.bowtie2.log"
        fi*/

}