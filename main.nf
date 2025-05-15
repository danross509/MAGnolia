#!/usr/bin/env nextflow

/*
To run this pipeline during development:
From the MAG_Pipeline/pipeline_test/ folder, run the following command:
    nextflow run ../main.nf -resume

DO NOT NAME {NANOPORE BARCODE FOLDERS} OR {PACBIO FILES} AS ***_#
THERE IS MAPPING CODE FOR SHORT READ SAMPLES TO FIRSTLY, ON IMPUT, CONVERT {_R1,_R2} TO {_1,_2}, AND LATER TO REMOVE THE {_1,_2} FOR SAMPLE ID EXTRACTION
USING {_##} NOTATION (ie. _01, _02) FOR LONG READ SAMPLEIDs IS UNAFFECTED



I am testing this pipeline with the following data:
    wget https://zenodo.org/record/3992790/files/test_reads.tar.gz
    tar -xzf test_reads.tar.gz
It is available in the folder pipeline_test/data/* in 3 formats:
    paired_X = Paired end reads with filename format 'sample1_{1,2}.fastq.gz'
    paired_RX = Paired end reads with filename format 'sampleR1_R{1,2}.fastq.gz'
    unpaired = Single endreads with filename format 'sample1_SE_{1,2}.fastq.gz'
To use a different dataset, change the folder name of the short_reads parameter in params.json
**For now, the <paired_reads> parameter must also be changed accordingly

Nextflow principally uses either Docker containers or conda environments to launch individual programs. 
I am developing it using conda environments. To run it you need to have a conda environment with nextflow installed
*/


include { QC_SHORT } from './subworkflows/local/qc_short/main.nf'
include { CONCATENATE_READS as CONCATENATE_ONT_BARCODES } from './modules/local/scripts/concatenate_reads/main.nf'
include { QC_NANOPORE } from './subworkflows/local/qc_nanopore/main.nf'
//include { QC_PACBIO} from './subworkflows/local/qc_pacbio/main.nf'
include { ASSEMBLY_SHORT } from './subworkflows/local/assembly_short/main.nf'
include { ASSEMBLY_LONG } from './subworkflows/local/assembly_long/main.nf'
include { BINNING_PREPARATION } from './subworkflows/nf-core/binning_preparation/main.nf'
include { BINNING } from './subworkflows/local/binning/main.nf'
include { domain_classification } from './subworkflows/nf-core/domain_classification/main.nf'
include { BINNING_REFINEMENT } from './subworkflows/nf-core/binning_refinement/main.nf'
include { DEPTHS } from './subworkflows/nf-core/depths/main.nf'
include { BIN_QC } from './subworkflows/nf-core/bin_qc/main.nf'

include { QUAST_BINS } from './modules/nf-core_mag/quast/quast_bins/main.nf'
include { QUAST_BINS_SUMMARY } from './modules/nf-core_mag/quast/quast_bins_summary/main.nf'

include { GTDBTK } from './subworkflows/nf-core_mag/gtdbtk/main.nf'
include { BIN_SUMMARY } from './modules/nf-core_mag/bin_summary/main.nf'

/*
The following modules are currently in development:

include { kraken2 } from './Kraken2/main.nf'
include { bracken } from './Bracken/main.nf'
include { quast } from './QUAST/main.nf'
include { busco } from './Busco/main.nf'
include { checkm } from './CheckM/main.nf' (include both 1 and 2)
include { taxonony } from './Taxonomy/main.nf' (include both ncbi, gtdb)
include { bakta } from './Bakta/main.nf'
include { annotation } from './Annotation/main.nf'
include { galah }
*/

//errorStrategy = { task.exitStatus in [12,143,137,104,134,139] ? 'retry' : 'finish' }

/*
for file in CLEAN_READS/*; do 
	mv $(readlink $file) CLEAN_READS/
done
*/

workflow {
    // Specify the output directory
    // projectDir = Channel.fromPath(params.projectDir, type: 'dir')

    /*
     * Import quality control parameters
     * If paired-read filenames are identified by R1/R2, the R will be removed for simplicity
     * PE and SE reads will be identified and handled as necessary for each software
     */

    ch_versions = Channel.empty()

     ////////////////////////////////////////////////////
    /* --  Create channel for reference databases  -- */
    ////////////////////////////////////////////////////

    /*if (params.host_genome) {
        host_fasta = params.genomes[params.host_genome].fasta ?: false
        ch_host_fasta = Channel.value(file("${host_fasta}"))
        host_bowtie2index = params.genomes[params.host_genome].bowtie2 ?: false
        ch_host_bowtie2index = Channel.value(file("${host_bowtie2index}/*"))
    }
    else if (params.host_fasta) {
        ch_host_fasta = Channel.value(file("${params.host_fasta}"))
    }
    else {
        ch_host_fasta = Channel.empty()
    }

    if (params.kraken2_db) {
        ch_kraken2_db_file = file(params.kraken2_db, checkIfExists: true)
    }
    else {
        ch_kraken2_db_file = []
    }

    if (params.cat_db) {
        ch_cat_db_file = Channel.value(file("${params.cat_db}"))
    }
    else {
        ch_cat_db_file = Channel.empty()
    }

    if (params.krona_db) {
        ch_krona_db_file = Channel.value(file("${params.krona_db}"))
    }
    else {
        ch_krona_db_file = Channel.empty()
    }

    if (!params.keep_phix) {
        ch_phix_db_file = Channel.value(file("${params.phix_reference}"))
    }*/

    gtdb = params.skip_binqc || params.skip_gtdbtk ? false : params.gtdb_db

    if ( gtdb ) {
        gtdb = file( "${gtdb}", checkIfExists: true )
        gtdb_mash = params.gtdb_mash ? file( "${params.gtdb_mash}", checkIfExists: true ) : []
    }
    else {
        gtdb = []
    }

    /****************
        Data input
     ****************/

    reads = Channel.fromPath ( params.reads_file )
        .splitCsv (header: true)
        .map { meta ->
            def meta_new = [:]
            def reads = []

            meta_new.id = meta.sampleID
            meta_new.sequencer = meta.sequencer
            if ( meta.paired_end == 'true' ) {
                meta_new.paired_end = true
            } else {
                meta_new.paired_end = false
            }
            if ( meta.corrected == 'true' ) {
                meta_new.corrected = true
            } else {
                meta_new.corrected = false
            }
            meta_new.bin_group = meta.bin_group
            meta_new.assembly_group = meta.assembly_group

            if ( meta_new.paired_end ) {
                return [ meta_new, [ meta.reads_R1, meta.reads_R2 ]]
            } else {
                return [ meta_new, [ meta.reads_R1 ]]
            }
        }
    
    //reads.view()
    short_reads = Channel.empty()
    short_reads = short_reads.mix ( reads )
        .map { meta, reads ->
            if ( meta.sequencer == 'Illumina' ) {
                return [ meta, reads ]
            } 
        }

    nanopore_reads = Channel.empty()
    nanopore_reads = nanopore_reads.mix ( reads )
        .map { meta, reads ->
            if ( meta.sequencer == 'ONT' ) {
                return [ meta, reads ]
            } 
        }

    pacbio_reads = Channel.empty()
    pacbio_reads = pacbio_reads.mix ( reads )
        .map { meta, reads ->
            if ( meta.sequencer == 'PacBio' ) {
                return [ meta, reads ]
            } 
        }

    // Set params for each read set
    
    
/*    //Create channel of raw short reads
    short_reads = Channel.empty()
    if ( params.short_reads ) {
        short_reads = Channel.fromFilePairs ( params.short_reads, size: -1 )
            .map { sample, reads ->
                def sampleID = sample.replaceAll(/_R$/, '')
                return [ sampleID, reads ]
            }
            .map { sampleID, reads ->
                def meta = [:]
                // Set meta.id
                meta.id = sampleID
                meta.sequencer = 'Illumina'
                // Set meta.paired_end
                if ( reads.size() == 2 ) {
                    meta.paired_end = true
                } else if ( reads.size() == 1 ) {
                    meta.paired_end = false
                } else {
                    exit 1, "ERROR: Check short read input files -> ${sampleID} contains invalid number of samples"
                }
                // Set bin-group
                meta.bin_group = "all"
                return [ meta, reads ]
            }
    }
    short_reads.view()
    
    
    //Create channel of raw ONT (Nanopore) reads
    nanopore_reads = Channel.empty()
    if ( params.nanopore_reads ) {
        nanopore_barcodes = Channel.fromPath ( "${params.nanopore_reads}/*" )
            .map { reads ->
                def meta = [:]
                // Set barcode folder name as meta.id
                def barcode = reads.getParent().getName()
                meta.id = barcode
                return [ meta, reads ]
            }
            .groupTuple()
            .map { meta, reads ->
                if ( !params.nanopore_reads_corrected ) {
                    def meta_new = meta + [ sequencer: 'ONT', corrected: false ]
                    return [ meta_new, reads ]
                } else {
                    def meta_new = meta + [ sequencer: 'ONT', corrected: true ]
                    return [ meta_new, reads ]
                }
            }
        //nanopore_barcodes.view()

        CONCATENATE_ONT_BARCODES ( 
            nanopore_barcodes,
            "",
            "RAW_ONT_BARCODES"
            )

        nanopore_reads = CONCATENATE_ONT_BARCODES.out

        nanopore_reads.view()
    }
    

    //Create channel of raw PacBio reads
    pacbio_reads = Channel.empty()
    if ( params.pacbio_reads ) {
        pacbio_reads = Channel.fromPath ( params.pacbio_reads )
    }
*/    



    /*********************
        Quality control
     *********************/
    // Create channel of phiX genome to remove
    phiX = Channel.empty()
    if ( !params.skip_qc && params.short_reads && params.remove_phiX ) {
        phiX = Channel.fromPath ( params.phiX )
            .map { reference ->
                def meta = [:]
                meta.id = reference.getBaseName()
                return [ meta, reference ]
            }
    }

    host_genome = Channel.empty()
    // Create channel of host genome to remove
    if ( !params.skip_qc && params.host_genome ) {
        host_genome = Channel.fromPath ( params.host_genome )
            .map { reference ->
                def meta = [:]
                meta.id = reference.getBaseName()
                return [ meta, reference ]
            }
    }

    // Short read quality control
    concatenated_reads = Channel.empty()        // Channel to concatenate fastq files for assembly ('per_sample' will be ungrouped sample fastq's)
    individual_clean_reads = Channel.empty()    // Channel to group original individual sample fastq's according to assembly file (for use in binning)
    if ( !params.skip_qc && params.short_reads ) {
        QC_SHORT ( 
            short_reads,
            phiX,
            host_genome
        )
        concatenated_reads = concatenated_reads.mix ( QC_SHORT.out.concatenated_reads )
        individual_clean_reads = individual_clean_reads.mix ( QC_SHORT.out.original_clean_reads )
    }
        
    concatenated_reads.view() 
    individual_clean_reads.view()  

    concatenated_long_reads = Channel.empty()
    individual_clean_long_reads = Channel.empty()
    // ONT quality control
    if ( !params.skip_qc && params.nanopore_reads ) {
        QC_NANOPORE ( 
            nanopore_reads,
            host_genome
        )
        concatenated_long_reads = concatenated_long_reads.mix ( QC_NANOPORE.out.concatenated_reads )
        individual_clean_long_reads = individual_clean_long_reads.mix ( QC_NANOPORE.out.original_clean_reads )
    }

    // PacBio quality control
    if ( !params.skip_qc && params.pacbio_reads ) {
        QC_PACBIO ( pacbio_reads
            // Host species
        )
        concatenated_long_reads = concatenated_long_reads.mix ( QC_PACBIO.out.concatenated_reads )
        individual_clean_long_reads = individual_clean_long_reads.mix ( QC_PACBIO.out.original_clean_reads )
    }

    concatenated_long_reads.view()

    /**************
        Assembly
     **************/

    
    final_contigs = Channel.empty()
    final_assembly_graphs = Channel.empty()

    if ( !params.skip_assembly ) {
        if ( params.short_reads && !params.nanopore_reads && !params.pacbio_reads ) {
            ASSEMBLY_SHORT ( concatenated_reads )
            final_contigs = final_contigs.mix ( ASSEMBLY_SHORT.out )
        } else if ( !params.short_reads && ( params.nanopore_reads || params.pacbio_reads )) {
            ASSEMBLY_LONG ( concatenated_long_reads )
            final_contigs = final_contigs.mix ( ASSEMBLY_LONG.out.assembly_out )
            final_assembly_graphs = final_assembly_graphs.mix ( ASSEMBLY_LONG.out.assembly_graph_out )
        } /*else if (params.short_reads && (params.nanopore_reads || params.pacbio_reads)) {
            ASSEMBLY_HYBRID( concatenated_reads )
            final_contigs = ASSEMBLY_HYBRID.out
        }*/
        
    }

    final_contigs.view()
    final_assembly_graphs.view()

    /*************
        Binning
     *************/

    binning_prep_input = Channel.empty()
    if ( !params.skip_binning ) {

        // 'Per_assembly' binning associates the read files used in each contig file, for all assembly modes
        // 'Grouped' binning concatenates the contig files of the selected samples, for 'per_sample' assembly only (for now)
        // 'Cobinning' concatenates all contig files and collects all individual read files, for all assembly modes

        binning_prep_reads = Channel.empty()
        // If short reads only
        if ( params.short_reads && !params.nanopore_reads && !params.pacbio_reads ) { 
            binning_prep_reads = binning_prep_reads.mix ( individual_clean_reads )
        // If long reads only    
        } else if ( !params.short_reads && ( params.nanopore_reads || params.pacbio_reads )) {
            binning_prep_reads = binning_prep_reads.mix ( individual_clean_long_reads )
        // If hybrid assembly
        } else if ( params.short_reads && ( params.nanopore_reads || params.pacbio_reads )) {
            binning_prep_reads = 'to do' // individual_clean_hybrid_reads
        }
            
        if ( params.binning_mode == 'per_assembly') {                
            binning_prep_reads = binning_prep_reads
                .map { meta, reads ->
                    [ meta.id, meta, reads ]
                }
            binning_prep_input = final_contigs
                .map { meta, contigs ->
                    [ meta.id, meta, contigs ]
                }
                .combine ( binning_prep_reads, by: 0 )
                .map { id, contigs_meta, contigs, grouped_reads_meta, grouped_reads ->
                    [ contigs_meta, contigs, grouped_reads_meta, grouped_reads ]
                }
        } else if ( params.binning_mode == 'grouped' ) {

        } else if ( params.binning_mode == 'cobinning' ) {

        } else {
                exit 1, "ERROR: Binniing mode <${params.binning_mode}> is invalid"
        }

        //binning_prep_input.view()                    

        BINNING_PREPARATION ( binning_prep_input )

        BINNING_PREPARATION.out.semibin2_mappings.view()
        /*BINNING (
            BINNING_PREPARATION.out.semibin2_mappings,
            BINNING_PREPARATION.out.grouped_mappings,
            final_assembly_graphs
        )*/
    }

    //BINNING.out.bins.view()

    /*
            ch_binning_results_bins = binning.out.bins.map { meta, bins ->
                def meta_new = meta + [domain: 'unclassified']
                [meta_new, bins]
            }
            ch_binning_results_unbins = binning.out.unbinned.map { meta, bins ->
                def meta_new = meta + [domain: 'unclassified']
                [meta_new, bins]
            }
    */
    ch_binning_results_bins = Channel.empty()
    ch_binning_results_unbins = Channel.empty()

/*    if ( !params.skip_binning && !params.skip_bin_domain_classification ) {

        domain_classification ( final_contigs,
                                BINNING.out.bins,
                                BINNING.out.unbinned )

        binning_results_bins = domain_classification.out.classified_bins.map { meta, bins ->
                def meta_new = meta + [ refinement: 'unrefined' ]
                [ meta_new, bins ]
        }
        binning_results_unbins = domain_classification.out.classified_unbins.map { meta, bins ->
                def meta_new = meta + [ refinement: 'unrefined' ]
                [ meta_new, bins ]
        }

    } else if ( !params.skip_binning && params.skip_bin_domain_classification) {

        ch_binning_results_bins = BINNING.out.bins.map { meta, bins ->
                def meta_new = meta + [ domain: 'unclassified', refinement: 'unrefined' ]
                [ meta_new, bins ]
        }
        ch_binning_results_unbins = BINNING.out.unbinned.map { meta, bins ->
                def meta_new = meta + [ domain: 'unclassified', refinement: 'unbinned_unrefined' ]
                [ meta_new, bins ]
        }

    }

      // Only run if 2+ binners are used*************************************************************************************************************
    if ( !params.skip_bin_refinement ) {
        ch_prokarya_bins_dastool = ch_binning_results_bins.filter { meta, bins ->
            meta.domain != "eukarya"
        }

        ch_eukarya_bins_dastool = ch_binning_results_bins.filter { meta, bins ->
            meta.domain == "eukarya"
        }

        //if (params.ancient_dna) {
        //    ch_contigs_for_binrefinement = ANCIENT_DNA_ASSEMBLY_VALIDATION.out.contigs_recalled
        //}
        //else {
            ch_contigs_for_binrefinement = BINNING_PREPARATION.out.grouped_mappings.map { meta, contigs, bam, bai -> [ meta, contigs ] } //}
        //}

        BINNING_REFINEMENT( ch_contigs_for_binrefinement, ch_prokarya_bins_dastool )
        // ch_refined_bins = ch_eukarya_bins_dastool
        //     .map{ meta, bins ->
        //             def meta_new = meta + [refinement: 'eukaryote_unrefined']
        //             [meta_new, bins]
        //         }.mix( BINNING_REFINEMENT.out.refined_bins)

        ch_refined_bins = BINNING_REFINEMENT.out.refined_bins
        ch_refined_unbins = BINNING_REFINEMENT.out.refined_unbins
        ch_versions = ch_versions.mix( BINNING_REFINEMENT.out.versions )

        if ( params.postbinning_input == 'raw_bins_only' ) {
            ch_input_for_postbinning_bins = ch_binning_results_bins
            ch_input_for_postbinning_bins_unbins = ch_binning_results_bins.mix( ch_binning_results_unbins )
        }
        else if ( params.postbinning_input == 'refined_bins_only' ) {
            ch_input_for_postbinning_bins = ch_refined_bins
            ch_input_for_postbinning_bins_unbins = ch_refined_bins.mix( ch_refined_unbins )
        }
        else if ( params.postbinning_input == 'both' ) {
            ch_all_bins = ch_binning_results_bins.mix( ch_refined_bins )
            ch_input_for_postbinning_bins = ch_all_bins
            ch_input_for_postbinning_bins_unbins = ch_all_bins.mix( ch_binning_results_unbins ).mix( ch_refined_unbins )
        }
    } else {
        ch_input_for_postbinning_bins = ch_binning_results_bins
        ch_input_for_postbinning_bins_unbins = ch_binning_results_bins.mix( ch_binning_results_unbins )
    }

    ch_input_for_postbinning = params.exclude_unbins_from_postbinning
        ? ch_input_for_postbinning_bins
        : ch_input_for_postbinning_bins_unbins

    if (!params.skip_binning){
        DEPTHS(ch_input_for_postbinning, BINNING.out.metabat2depths, clean_reads)
        ch_input_for_binsummary = DEPTHS.out.depths_summary
    }
*/    
    
    

    /*
    * Bin QC subworkflows: for checking bin completeness with either BUSCO, CHECKM, CHECKM2, and/or GUNC
    */

    if (!params.skip_binqc) {
        BIN_QC(ch_input_for_postbinning)

        ch_bin_qc_summary = BIN_QC.out.qc_summary
        ch_versions = ch_versions.mix(BIN_QC.out.versions)
    }
    
/*    ch_quast_bins_summary = Channel.empty()
    if (!params.skip_quast) {
        ch_input_for_quast_bins = ch_input_for_postbinning
            .groupTuple()
            .map { meta, bins ->
                def new_bins = bins.flatten()
                [meta, new_bins]
            }

        QUAST_BINS(ch_input_for_quast_bins)
        ch_versions = ch_versions.mix(QUAST_BINS.out.versions.first())
        ch_quast_bin_summary = QUAST_BINS.out.quast_bin_summaries.collectFile(keepHeader: true) { meta, summary ->
            ["${meta.id}.tsv", summary]
        }
        QUAST_BINS_SUMMARY(ch_quast_bin_summary.collect())
        ch_quast_bins_summary = QUAST_BINS_SUMMARY.out.summary
    
        }
*/
    // If CAT is not run, then the CAT global summary should be an empty channel
        /*if (params.cat_db_generate || params.cat_db) {
            ch_cat_global_summary = CAT_SUMMARY.out.combined
        }
        else {*/
            ch_cat_global_summary = Channel.empty()
        //}

    /*
     * GTDB-tk: taxonomic classifications using GTDB reference
     */
    /*
    if (!params.skip_gtdbtk) {

        ch_gtdbtk_summary = Channel.empty()
        if (gtdb) {

            ch_gtdb_bins = ch_input_for_postbinning.filter { meta, bins ->
                meta.domain != "eukarya"
            }

            GTDBTK(
                ch_gtdb_bins,
                ch_bin_qc_summary,
                gtdb,
                gtdb_mash
            )
            ch_versions = ch_versions.mix(GTDBTK.out.versions.first())
            ch_gtdbtk_summary = GTDBTK.out.summary
        }
    }
    else {
        ch_gtdbtk_summary = Channel.empty()
    }

    if ((!params.skip_binqc) || !params.skip_quast || !params.skip_gtdbtk) {
        BIN_SUMMARY(
            ch_input_for_binsummary,
            ch_bin_qc_summary.ifEmpty([]),
            ch_quast_bins_summary.ifEmpty([]),
            ch_gtdbtk_summary.ifEmpty([]),
            ch_cat_global_summary.ifEmpty([]),
            params.binqc_tool
        )
    }
    */
    /*
     * Prokka: Genome annotation
     */

    /*if (!params.skip_prokka) {
        ch_bins_for_prokka = ch_input_for_postbinning
            .transpose()
            .map { meta, bin ->
                def meta_new = meta + [id: bin.getBaseName()]
                [meta_new, bin]
            }
            .filter { meta, bin ->
                meta.domain != "eukarya"
            }

        PROKKA(
            ch_bins_for_prokka,
            [],
            []
        )
        ch_versions = ch_versions.mix(PROKKA.out.versions.first())
    }
    */
    
}

//  In development:
/*  Quality Control
        Currently the QC module produces a fastqc/multiqc report for the data before and after this step
        Trimmomatic is used to remove adapters and low quality base pairs
        Bowtie2 is used to remove contaminant/host reads, including PhiX as a default

        I need to incorporate steps for long reads
        I need to modify bowtie2 filter contaminants to also process singleton reads from paired samples
        I would like to include a deduplication step, from the following: 
            bbtools' clumpify (which does not handle large datasets well)
            fastp deduplication (though this program seems to have poor support from the developers)

    Assembly
        I am using megahit to test the assembly module / parameters in the pipeline, due to it's speed
        The final assembly module will have the following assembly algorithms to choose from:
            megahit (https://github.com/voutcn/megahit)
            metaSPAdes (https://github.com/ablab/spades)
            metaHipMer (https://sourceforge.net/projects/hipmer/)
            GATB (https://gatb.inria.fr/software/de-novo-genome-assembly/)
            metaFlye (https://github.com/mikolmogorov/Flye)
            Canu (https://github.com/marbl/canu)

        Kracken/Bracken will provide an initial overview of the taxa present in the samples
        (https://github.com/DerrickWood/kraken2/wiki)
        (https://github.com/jenniferlu717/Bracken)
        
    Binning
        I will include the following as binning options:
            Metabat2 (https://bitbucket.org/berkeleylab/metabat/src/master/)
            Maxbin2 (https://sourceforge.net/projects/maxbin2/)
            CONCOCT (https://github.com/BinPro/CONCOCT)
            AAMB/AVAMB/VAMB (https://github.com/RasmussenLab/vamb?tab=readme-ov-file)
            Semibin2 (https://github.com/BigDataBiology/SemiBin)
            COMEBin (https://github.com/ziyewang/COMEBin)
        Along with one (or both) of the following refinement tools:
            metaWRAP bin refinement (https://github.com/bxlab/metaWRAP)
            DASTools (https://github.com/cmks/DAS_Tool)

        CheckM/CheckM2 will validate the contamination and completeness of each bin
        (https://github.com/Ecogenomics/CheckM)
        (https://github.com/chklovski/CheckM2)

        Sourmash's gather alorithm will provide information on MAG distribution across samples
        and what proportion of the samples don't map back to the MAGs
        (https://github.com/sourmash-bio/sourmash)

    Taxonomy
        Users will have the option to classify bins using either the GTDB (https://github.com/Ecogenomics/GTDBTk)
        or NCBI (https://github.com/ncbi/blast_plus_docs) databases
    
    Annotation
        This is the step I have looked into the least.
            KEGG
            Eggnog
            DRAM
            CAZy
            Bakta/Prokka
*/

