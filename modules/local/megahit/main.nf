#!/usr/bin/env nextflow

process megahit {

    container "community.wave.seqera.io/library/megahit:1.2.9--23234b8da1e27898"
    conda "bioconda::megahit=1.2.9"

    publishDir "${launchDir}/Assembly/megahit", mode: 'symlink'

    input:
        tuple val(meta), path(reads_fastq)
        val mh_preset

    output:
        path("final.contigs.fa")

    script:

    if (meta.paired_end) {

        reads_1 = 
        reads_2 = 
        """
        megahit \
        -1 $reads_1 -2 $reads_2 \
        -o megahit \
        --tmp-dir tmp \
        --presets $mh_preset \
        --verbose \
        --continue
        """
    }/* else if (!meta.paired_end) {

        """
        megahit \
        -1 $reads_1 -2 $reads_2 \
        -o megahit \
        --tmp-dir tmp \
        -t $threads \
        -m ${mem}000000000 \
        --presets meta-large \
        --verbose \
        --continue
        """
    }*/
}