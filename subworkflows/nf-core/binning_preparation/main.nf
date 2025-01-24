#!/usr/bin/env nextflow

/*
 * Binning preparation with Bowtie2
 */

include { BOWTIE2_BUILD_INDEX } from '../../../modules/local/bowtie2/build_index/main.nf'
include { BOWTIE2_ASSEMBLY_ALIGNMENT } from '../../../modules/local/bowtie2/assembly_alignment/main.nf'

workflow BINNING_PREPARATION {
    take:
    assemblies           // channel: [ val(meta), path(assembly) ]
    reads                // channel: [ val(meta), [ reads ] ]

    main:
    // build bowtie2 index for all assemblies
    BOWTIE2_BUILD_INDEX ( assemblies )

    // combine assemblies with sample reads for binning depending on specified mapping mode
    /*if (params.binning_map_mode == 'all'){
        // combine assemblies with reads of all samples
        ch_bowtie2_input = BOWTIE2_BUILD_INDEX.out.assembly_index
            .combine(reads)
    } else if (params.binning_map_mode == 'group'){
        // combine assemblies with reads of samples from same group
        ch_reads_bowtie2 = reads.map{ meta, reads -> [ meta.group, meta, reads ] }
        ch_bowtie2_input = BOWTIE2_BUILD_INDEX.out.assembly_index
            .map { meta, assembly, index -> [ meta.group, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { group, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }

    } else {
    */    // i.e. --binning_map_mode 'own'
        // combine assemblies (not co-assembled) with reads from own sample
        ch_reads_bowtie2 = reads.map{ meta, reads -> [ meta.id, meta, reads ] }
        ch_bowtie2_input = BOWTIE2_BUILD_INDEX.out.assembly_index
            .map { meta, assembly, index -> [ meta.id, meta, assembly, index ] }
            .combine(ch_reads_bowtie2, by: 0)
            .map { id, assembly_meta, assembly, index, reads_meta, reads -> [ assembly_meta, assembly, index, reads_meta, reads ] }

    //}

    BOWTIE2_ASSEMBLY_ALIGNMENT ( ch_bowtie2_input )
    // group mappings for one assembly
    ch_grouped_mappings = BOWTIE2_ASSEMBLY_ALIGNMENT.out.mappings
        .groupTuple(by: 0)
        .map { meta, assembly, bams, bais -> [ meta, assembly.sort()[0], bams, bais ] }     // multiple symlinks to the same assembly -> use first of sorted list

    emit:
    bowtie2_assembly_multiqc = BOWTIE2_ASSEMBLY_ALIGNMENT.out.log.map { assembly_meta, reads_meta, log -> [ log ] }
    //bowtie2_version          = BOWTIE2_ASSEMBLY_ALIGNMENT.out.versions
    grouped_mappings         = ch_grouped_mappings
}