#!/usr/bin/env nextflow

/*
To run this pipeline during development:
From the MAG_Pipeline/pipeline_test/[short_only | nanopore_only | pacbio_only | hybrid] folder, run the following command:
    nextflow run ../../main.nf -resume

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
include { QC_NANOPORE } from './subworkflows/local/qc_nanopore/main.nf'
include { QC_PACBIO} from './subworkflows/local/qc_pacbio/main.nf'
include { READ_CONTIG_TAXONOMY as READ_TAXONOMY } from './subworkflows/local/read_contig_taxonomy/main.nf'
include { ASSEMBLY_PREP_SHORT } from './subworkflows/local/assembly_prep_short/main.nf'
include { ASSEMBLY_PREP_LONG } from './subworkflows/local/assembly_prep_long/main.nf'
include { ASSEMBLY_SHORT } from './subworkflows/local/assembly_short/main.nf'
include { ASSEMBLY_LONG } from './subworkflows/local/assembly_long/main.nf'
include { READ_CONTIG_TAXONOMY as CONTIG_TAXONOMY } from './subworkflows/local/read_contig_taxonomy/main.nf'
include { BINNING_PREPARATION } from './subworkflows/local/binning_preparation/main.nf'
include { BINNING } from './subworkflows/local/binning/main.nf'
//include { domain_classification } from './subworkflows/nf-core/domain_classification/main.nf'
include { BIN_REFINEMENT } from './subworkflows/local/bin_refinement/main.nf'
include { BIN_DEREPLICATION } from './subworkflows/local/bin_dereplication/main.nf'
include { BIN_COVERAGE } from './subworkflows/local/bin_coverage/main.nf'
include { BIN_ANNOTATION } from './subworkflows/local/bin_annotation/main.nf'

//include { DEPTHS } from './subworkflows/nf-core/depths/main.nf'
include { BIN_QC } from './subworkflows/nf-core/bin_qc/main.nf'

include { QUAST_BINS } from './modules/local/quast/quast_bins/main.nf'
include { QUAST_BINS_SUMMARY } from './modules/local/quast/quast_bins_summary/main.nf'

include { GTDBTK } from './subworkflows/nf-core_mag/gtdbtk/main.nf'
include { BIN_SUMMARY } from './modules/nf-core_mag/bin_summary/main.nf'

include { KRAKEN2_DB_DOWNLOAD } from './modules/local/kraken2/db_download/main.nf'
include { BRACKEN_BUILD } from './modules/local/bracken/build/main.nf'
include { DRAM_SETUP as DRAM_IMPORT_CONFIG } from './modules/local/dram/setup/main.nf'
include { DRAM_SETUP as DRAM_PREPARE_DB } from './modules/local/dram/setup/main.nf'

/*
The following modules are currently in development:

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
    db_download_dir = file("${params.databaseDownloadDir}").toAbsolutePath().toString()
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
    }*/

    // Kraken2 database
    // If there is a kraken database given
    if ( params.kraken2_db ) {
        kraken2_db_dir = file ( params.kraken2_db, checkIfExists: true )
        // If that database will be used AND if bracken will be run
        if ( (!params.skip_read_taxonomy || !params.skip_contig_taxonomy) && !params.skip_kracken2 && !params.skip_bracken ) {
            // If there is no bracken build, build it
            if ( !params.bracken_build_exists ) {
                BRACKEN_BUILD (
                kraken2_db_dir.toAbsolutePath().toString(),
                params.bracken_kmer_len,
                params.bracken_read_length,
                params.lowThreads
                )

                // Set run_bracken to true
                run_bracken = BRACKEN_BUILD.out.bracken_built
            // Otherwise if bracken is already built
            } else run_bracken = true
        }
    // If no database is specified but Kraken2 will be used
    } else if ( (!params.skip_read_taxonomy || !params.skip_contig_taxonomy) && !params.skip_kracken2 ) {
        KRAKEN2_DB_DOWNLOAD (
            db_download_dir,
            params.kraken2_build,
            params.kraken2_kmer_len, 
            params.kraken2_max_db_size
        )
        
        kraken2_db_dir = KRAKEN2_DB_DOWNLOAD.out.directory

        // If Bracken will also be used
        if ( !params.skip_bracken ) {
            BRACKEN_BUILD (
                kraken2_db_dir,
                params.bracken_kmer_len,
                params.bracken_read_length,
                params.lowThreads
            )

            run_bracken = BRACKEN_BUILD.out.bracken_built
        }

    // If no Kraken2 database is given and Kraken2 will not be used
    } else {
        kraken2_db_dir = []
        run_bracken = []
    }
        
    // DRAM databases
    // If there is a DRAM config given
    if ( params.dram_config_loc ) {
        dram_config_loc = file ( params.dram_config_loc, checkIfExists: true )
        // AND If DRAM will be used
        if ( !params.skip_annotation && !params.skip_dram ) {
            DRAM_IMPORT_CONFIG (
                dram_config_loc.toAbsolutePath().toString(),
                db_download_dir,
                params.dram_kegg_loc,
                params.dram_skip_uniref
            )
        }

    // If no DRAM config is given but DRAM will be used
    } else if ( !params.skip_annotation && !params.skip_dram ) {
        if ( params.dram_kegg_loc ) {
            kegg_file = file ( params.dram_kegg_loc, checkIfExists: true ).toAbsolutePath().toString()
        } else {
            kegg_file = false
        }

        DRAM_PREPARE_DB (
            params.dram_config_loc,
            db_download_dir,
            kegg_file,
            params.dram_skip_uniref
        )

    // If no DRAM config is given and DRAM will not be used
    } else {

    }


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
            if ( meta.paired_end == 'True' ) {
                meta_new.paired_end = true
            } else {
                meta_new.paired_end = false
            }
            if ( meta.corrected == 'True' ) {
                meta_new.corrected = true
            } else {
                meta_new.corrected = false
            }
            meta_new.assembly_group = meta.assembly_group
            meta_new.bin_group = meta.bin_group

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

    /*********************
        Quality control
     *********************/
    // Create channel of phiX genome to remove
    //phiX = Channel.empty()
    phiX_index = Channel.empty()
    if ( !params.skip_qc && params.short_reads && params.remove_phiX ) {
        /*phiX = Channel.fromPath ( params.phiX )
            .map { reference ->
                def meta = [:]
                meta.id = reference.getBaseName()
                return [ meta, reference ]
            }
        */
        phiX_index = Channel.fromPath ( "${projectDir}/reference_genomes/phiX/*.bt2" )
            .map { index ->
                def meta =[:]
                meta.id = 'phiX174'
                return [ meta, index ]
            }
            .groupTuple()
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
    corrected_reads = Channel.empty()
    if ( !params.skip_qc && params.short_reads ) {
        QC_SHORT ( 
            short_reads,
            phiX_index,
            host_genome
        )

        corrected_reads = corrected_reads.mix ( QC_SHORT.out.host_filtered_reads )
    }

    //concatenated_long_reads = Channel.empty()
    //original_clean_long_reads = Channel.empty()
    // ONT quality control
    corrected_ont_reads = Channel.empty()
    if ( !params.skip_qc && params.nanopore_reads ) {
        QC_NANOPORE ( 
            nanopore_reads,
            host_genome
        )
        
        corrected_ont_reads = corrected_ont_reads.mix ( QC_NANOPORE.out.filtered_nanopore_reads )
    }

    // PacBio quality control
    corrected_pacbio_reads = Channel.empty()
    if ( !params.skip_qc && params.pacbio_reads ) {
        QC_PACBIO ( 
            pacbio_reads,
            host_genome
        )

        corrected_pacbio_reads = corrected_pacbio_reads.mix ( QC_PACBIO.out.filtered_pacbio_reads )
    }

    all_corrected_reads = corrected_reads.mix ( corrected_ont_reads, corrected_pacbio_reads )
    // if this is empty, exit

    //concatenated_long_reads.view()

    /*******************
        Read taxonomy
     *******************/

    read_taxonomy_input = Channel.empty()
    if ( !params.skip_read_taxonomy ) {
        READ_TAXONOMY (
            all_corrected_reads, 
            kraken2_db_dir,
            run_bracken,
            "reads"
        )
    }

    /**************
        Assembly
     **************/
    // **************    
    concatenated_reads = Channel.empty()        // Channel to concatenate fastq files for assembly ('per_sample' will be ungrouped sample fastq's)
    original_clean_reads = Channel.empty()    // Channel to group original individual sample fastq's according to assembly file (for use in binning)
    if ( !params.skip_assembly && params.short_reads ) {
        ASSEMBLY_PREP_SHORT (
            corrected_reads
        )
        
        concatenated_reads = concatenated_reads.mix ( ASSEMBLY_PREP_SHORT.out.concatenated_reads )
        original_clean_reads = original_clean_reads.mix ( ASSEMBLY_PREP_SHORT.out.original_clean_reads )
    }

    concatenated_long_reads = Channel.empty()       // Channel to concatenate fastq files for assembly ('per_sample' will be ungrouped sample fastq's)
    original_clean_long_reads = Channel.empty()     // Channel to group original individual sample fastq's according to assembly file (for use in binning)
    if ( !params.skip_assembly && ( params.nanopore_reads || params.pacbio_reads ) ) {
        ASSEMBLY_PREP_LONG (
            corrected_ont_reads,
            corrected_pacbio_reads
        )
        
        concatenated_long_reads = concatenated_long_reads.mix ( ASSEMBLY_PREP_LONG.out.concatenated_long_reads )
        original_clean_long_reads = original_clean_long_reads.mix ( ASSEMBLY_PREP_LONG.out.original_clean_long_reads )
    }
    // *****************
    
    final_contigs = Channel.empty()
    assembly_graphs = Channel.empty()
    reads_post_assembly = Channel.empty()
    hifiasm_bins = Channel.empty()

    if ( !params.skip_assembly ) {
        if ( params.short_reads && !params.nanopore_reads && !params.pacbio_reads ) {
            ASSEMBLY_SHORT ( concatenated_reads, original_clean_reads )
            final_contigs = final_contigs.mix ( ASSEMBLY_SHORT.out.contigs )
            assembly_graphs = assembly_graphs.mix ( ASSEMBLY_SHORT.out.assembly_graph )
            reads_post_assembly = reads_post_assembly.mix ( ASSEMBLY_SHORT.out.reads )
        } else if ( !params.short_reads && ( params.nanopore_reads || params.pacbio_reads )) {
            ASSEMBLY_LONG ( concatenated_long_reads, original_clean_long_reads )
            final_contigs = final_contigs.mix ( ASSEMBLY_LONG.out.final_contigs )
            assembly_graphs = assembly_graphs.mix ( ASSEMBLY_LONG.out.assembly_graphs )
            reads_post_assembly = reads_post_assembly.mix ( ASSEMBLY_LONG.out.final_reads )
            hifiasm_bins = hifiasm_bins.mix ( ASSEMBLY_LONG.out.hifiasm_bins )
        } /*else if (params.short_reads && (params.nanopore_reads || params.pacbio_reads)) {
            ASSEMBLY_HYBRID( concatenated_reads )
            final_contigs = ASSEMBLY_HYBRID.out
        }*/
        
    }

    //final_contigs.view()
    //assembly_graphs.view()
    //reads_post_assembly.view()

    /*********************
        Contig taxonomy
     *********************/

    if ( !params.skip_contig_taxonomy ) {
        CONTIG_TAXONOMY (
            final_contigs, 
            kraken2_db_dir,
            run_bracken,
            "contigs"
        )
    }

    /*************
        Binning
     *************/

    binning_prep_input = Channel.empty()
    initial_bins = Channel.empty()
    post_refinement_bins = Channel.empty()
    if ( !params.skip_binning ) {

        /*binning_prep_reads = Channel.empty()
        // If short reads only
        if ( params.short_reads && !params.nanopore_reads && !params.pacbio_reads ) { 
            binning_prep_reads = binning_prep_reads.mix ( original_clean_reads )
        // If long reads only    
        } else if ( !params.short_reads && ( params.nanopore_reads || params.pacbio_reads )) {
            binning_prep_reads = binning_prep_reads.mix ( original_clean_long_reads )
        // If hybrid assembly
        } else if ( params.short_reads && ( params.nanopore_reads || params.pacbio_reads )) {
            binning_prep_reads = 'to do' // original_clean_hybrid_reads
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
        }*/

        binning_prep_input = binning_prep_input.mix ( final_contigs )
            .join ( assembly_graphs )
            .join ( reads_post_assembly )

        //binning_prep_input.view()                    

        BINNING_PREPARATION ( binning_prep_input )

        //BINNING_PREPARATION.out.grouped_mappings.view()

        BINNING (
            BINNING_PREPARATION.out.grouped_mappings,
            hifiasm_bins
        )

        //BINNING.out.bins.view()

        refinement_contigs = BINNING_PREPARATION.out.grouped_mappings
            .map { meta, reads, contigs, bams, bais, gfa ->
                def meta_new = meta + [refined: false]
                [ meta_new, contigs ]
            }

        initial_bins = initial_bins.mix ( BINNING.out.bins )
            .map { meta, bins ->
                def meta_new = meta + [refined: false]
                [meta_new, bins]
            }   

        if ( !params.skip_bin_refinement ) {
            BIN_REFINEMENT (
                refinement_contigs,
                initial_bins
            )

            post_refinement_bins = post_refinement_bins.mix ( BIN_REFINEMENT.out.refined_bins )

            //BIN_REFINEMENT.out.refined_bins.view()

        } else {
            post_refinement_bins = post_refinement_bins.mix ( initial_bins )
        }
            
    }

    //*****************
    // DAStool and drep arguments
    

    /********************
        Bin evaluation
     ********************/

    // CheckM

    /*******************
        Dereplication
     *******************/

    final_bins = Channel.empty()

    if ( !params.skip_bin_dereplication ) {
        BIN_DEREPLICATION ( post_refinement_bins )

        final_bins = final_bins.mix ( BIN_DEREPLICATION.out.dereplicated_bins )

    } else {
        final_bins = final_bins.mix ( post_refinement_bins )
    }

    /************************
        Bin classification
     ************************/

    // gtdb

    /********************
        Bin annotation
     ********************/

    if ( !params.skip_annotation ) {
       BIN_ANNOTATION ( final_bins )
    }

    /******************
        Bin coverage
     ******************/
    
    
    bin_group_reads = BINNING_PREPARATION.out.grouped_mappings
        .map { meta, reads, contigs, bams, bais, gfa ->
            [ meta, reads ]
        }

    if ( !params.skip_bin_coverage ) {
        BIN_COVERAGE (
            all_corrected_reads,        // Calculate coverage for each original read sample...
            bin_group_reads,            // ... as well as concatenated reads used in each binning instance (if applicable)
            final_bins                  
        )
    }



    // Manage file storage


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

    if ( !params.skip_binqc ) {
        BIN_QC ( final_bins )

        ch_bin_qc_summary = BIN_QC.out.qc_summary
        ch_versions = ch_versions.mix ( BIN_QC.out.versions )
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

/*  CheckM/CheckM2 will validate the contamination and completeness of each bin
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

