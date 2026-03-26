#!/usr/bin/env nextflow

process RENAME_REFINED_BINS {
    tag "${meta.id}-${meta.assembler}"
    label 'process_single'

    conda ""
    container ""

    input:
    tuple val(meta), path(bins)

    output:
    tuple val(meta), path("${meta.id}_${meta.assembler}*Refined.*.fa"),         optional:true, emit: refined_bins
    tuple val(meta), path("${meta.id}_${meta.assembler}_DASToolUnbinned.fa"),   optional:true, emit: refined_unbins

    script:
    """
    for bin in $bins; do
        if [[ "\${bin}" == unbinned.fa ]]; then
            mv unbinned.fa ${meta.id}_${meta.assembler}_DASToolUnbinned.fa
        elif [[ -f \$bin ]]; then
            filename=\${bin##*/}
            basename=\${filename%%.*}
            remainder=\${filename#*.}
            mv \$bin ./\${basename}_Refined.\${remainder}
        fi
    done
    """
}