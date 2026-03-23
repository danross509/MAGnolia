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
include { BIN_REFINEMENT } from './subworkflows/local/bin_refinement/main.nf'
include { BIN_DEREPLICATION } from './subworkflows/local/bin_dereplication/main.nf'

include { BIN_EVALUATION } from './subworkflows/local/bin_evaluation/main.nf'
include { BIN_CLASSIFICATION } from './subworkflows/local/bin_classification/main.nf'
include { BIN_COVERAGE } from './subworkflows/local/bin_coverage/main.nf'
include { BIN_ANNOTATION } from './subworkflows/local/bin_annotation/main.nf'

//include { DEPTHS } from './subworkflows/nf-core/depths/main.nf'
include { BIN_QC } from './subworkflows/nf-core/bin_qc/main.nf'

include { QUAST_BINS } from './modules/local/quast/quast_bins/main.nf'
include { QUAST_BINS_SUMMARY } from './modules/local/quast/quast_bins_summary/main.nf'

//include { BIN_SUMMARY } from './modules/nf-core_mag/bin_summary/main.nf'

//include { KRAKEN2_DB_DOWNLOAD } from './modules/local/kraken2/db_download/main.nf'
include { K2_DOWNLOAD_TAXONOMY } from './modules/local/kraken2/k2_download_taxonomy/main.nf'
include { K2_BUILD } from './modules/local/kraken2/k2_build/main.nf'
include { KRAKEN2_UPDATE_CONFIG } from './modules/local/kraken2/update_config/main.nf'
include { BRACKEN_BUILD } from './modules/local/bracken/build/main.nf'
include { BRACKEN_UPDATE_CONFIG } from './modules/local/bracken/update_config/main.nf'
include { BAKTA_BAKTADBDOWNLOAD } from './modules/nf-core/bakta/baktadbdownload/main.nf'
include { BAKTA_UPDATE_CONFIG } from './modules/local/bakta/update_config/main.nf'
include { DRAM_SETUP as DRAM_IMPORT_CONFIG } from './modules/local/dram/setup/main.nf'
include { DRAM_SETUP as DRAM_PREPARE_DB } from './modules/local/dram/setup/main.nf'
include { DRAM_UPDATE_CONFIG } from './modules/local/dram/update_config/main.nf'
include { CHECKM2_DATABASEDOWNLOAD } from './modules/nf-core/checkm2/databasedownload/main.nf'
include { CHECKM2_UPDATE_CONFIG } from './modules/local/checkm2/update_config/main.nf'
include { UNTAR as CHECKM_UNTAR } from './modules/nf-core/untar/main.nf'
include { CHECKM_UPDATE_CONFIG } from './modules/local/checkm/update_config/main.nf'
include { GTDB_DB_DOWNLOAD } from './modules/local/gtdbtk/db_download/main.nf'
include { GTDB_UPDATE_CONFIG } from './modules/local/gtdbtk/update_config/main.nf'

workflow {
    ch_versions = channel.empty()

    /********************
        Database setup
     ********************/
    db_download_dir = file("${params.databaseDownloadDir}").toAbsolutePath().toString()

    // Kraken2 database
    // If there is a kraken database given
    if ( params.kraken2_db ) {
        kraken2_db_dir = file ( params.kraken2_db, checkIfExists: true )
        if ( kraken2_db_dir ) {
            println ( "Kraken2 database found at ${kraken2_db_dir}" )
        } else {
            exit 1 ( "ERROR: Kraken2 database at ${params.kraken2_db} not found" )
        }
    // If no database is specified but Kraken2 will be used
    } else if (( !params.skip_read_taxonomy || !params.skip_contig_taxonomy ) && !params.skip_kracken2 ) {
        println ( "Downloading Kraken2 database at ${db_download_dir}/kraken2_db" )
        K2_DOWNLOAD_TAXONOMY (
            db_download_dir
        )

        K2_BUILD (
            K2_DOWNLOAD_TAXONOMY.out.directory,
            params.kraken2_build,
            params.kraken2_kmer_len, 
            params.kraken2_max_db_size
        )

        //kraken2_db_dir = KRAKEN2_DB_DOWNLOAD.out.directory
        kraken2_db_dir = K2_BUILD.out.directory

        KRAKEN2_UPDATE_CONFIG (
            kraken2_db_dir
        )

    // If no Kraken2 database is given and Kraken2 will not be used
    } else {
        kraken2_db_dir = []
    }

    // If Kraken2 will be used AND if bracken will be run
    if (( !params.skip_read_taxonomy || !params.skip_contig_taxonomy ) && !params.skip_kracken2 && !params.skip_bracken ) {
        // If there is no bracken build, build it
        println(kraken2_db_dir.toAbsolutePath().toString())
        if ( !params.bracken_build_exists ) {
            if ( !params.bracken_memory_map && params.maxMem >= 240.GB ) {
                bracken_mm = false
            } else {
                bracken_mm = true
            }

            BRACKEN_BUILD (
            kraken2_db_dir.toAbsolutePath().toString(),
            params.bracken_kmer_len,
            params.bracken_read_length,
            bracken_mm
            )

            BRACKEN_UPDATE_CONFIG (
                BRACKEN_BUILD.out.bracken_built
            )

            // Set run_bracken to true
            run_bracken = BRACKEN_BUILD.out.bracken_built

        // Otherwise if bracken is already built
        } else {
            run_bracken = true
        }
    } else {
        run_bracken = false
    }
    
    // Bin evaluation (CheckM/CheckM2) databases
    if ( !params.skip_bin_evaluation) {
        // If running CheckM2
        useCheckM2 = ['checkm2', 'both']
        if ( useCheckM2.contains ( params.checkm_version ) ) {
            // If a database is specified
            if ( params.checkm2_db ) {
                checkm2_db = file ( params.checkm2_db, checkIfExists: true )

                if ( checkm2_db ) {
                    println ( "CheckM2 database found at ${params.checkm2_db}" )
                    checkm2_db_dir = [[ id: 'checkm2_db' ], checkm2_db ]
                    
                } else {
                    exit 1 ( "ERROR: checkm2_db path ${params.checkm2_db} does not exist" )
                }

            // If not, download it
            } else {
                CHECKM2_DATABASEDOWNLOAD ( params.checkm2_db_version )
                ch_versions = ch_versions.mix ( CHECKM2_DATABASEDOWNLOAD.out.versions )

                // tuple val(meta), path("checkm2_db_v${db_version}.dmnd")

                CHECKM2_UPDATE_CONFIG (
                    CHECKM2_DATABASEDOWNLOAD.out.database,
                    params.checkm2_download_dir
                )

                checkm2_db_dir = CHECKM2_UPDATE_CONFIG.out.db
            }

            if ( params.checkm_version == 'checkm2') {
                checkm_db_dir = []
            }

        } 
        
        // If running CheckM
        useCheckM = ['checkm', 'both']
        if ( useCheckM.contains ( params.checkm_version ) ) {
            // If a database is specified
            if ( params.checkm_db ) {
                //checkm_db_dir = [[ id: 'checkm_db' ], file ( params.checkm_db, checkIfExists: true )]
                checkm_db_dir = [file ( params.checkm_db, checkIfExists: true )]

            // If not, download it
            } else {
                checkm_untar_input = [[ id: 'checkm_db' ], file ( params.checkm_download_url, checkIfExists: true )]

                CHECKM_UNTAR ( checkm_untar_input )
                CHECKM_UNTAR.out.untar.view()
                //checkm_db_dir = CHECKM_UNTAR.out.untar
                 //   .map { it[1] }
                ch_versions = ch_versions.mix ( CHECKM_UNTAR.out.versions )

                CHECKM_UPDATE_CONFIG (
                    CHECKM_UNTAR.out.untar,
                    params.checkm2_download_dir
                )

                checkm_db_dir = CHECKM_UPDATE_CONFIG.out.db
            }

            if ( params.checkm_version == 'checkm') {
                checkm2_db_dir = []
            }
        }

        if ( !useCheckM2.contains ( params.checkm_version ) && !useCheckM.contains ( params.checkm_version ) ) {
            // Error message and quit
        }

    } else {
        checkm2_db_dir = []
        checkm_db_dir = []
    }


    // Bin annotation databases
    if ( !params.skip_annotation ) {
        // If DRAM will be used
        if ( !params.skip_dram ) {
            // AND If there is a DRAM config given
            if ( params.dram_config_loc ) {
                dram_config_loc = file ( params.dram_config_loc, checkIfExists: true )
                if ( dram_config_loc ) { 
                    println("DRAM database found at ${dram_config_loc}")

                    DRAM_IMPORT_CONFIG (
                        dram_config_loc.toAbsolutePath().toString(),
                        db_download_dir,
                        params.dram_kegg_loc,
                        params.dram_skip_uniref
                    )
                } else {
                    exit 1 ( "ERROR: DRAM database at ${dram_config_loc} not found" )
                }
            }

            // If no DRAM config is given but DRAM will be used
            } else if ( !params.dram_config_loc ) {
                // If a kegg file is available
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

                DRAM_UPDATE_CONFIG (
                    db_download_dir,
                    DRAM_PREPARE_DB.out
                )

        } 
        // If Bakta will be used
        if ( !params.skip_bakta ) {
            // AND If there is a Bakta directory given
            if ( params.bakta_db ) {
                bakta_db_dir = file ( params.bakta_db, checkIfExists: true )

                if ( bakta_db_dir ) {
                    //bakta_db_dir = bakta_db_dir.toAbsolutePath().toString()
                    println ( "Bakta database found at ${params.bakta_db}" )
                } else {
                    exit 1 ( "ERROR: bakta_db path ${params.bakta_db} does not exist" )
                }

            // If no Bakta directory is given
            } else {
                println ( "Bakta database not given, downloading to ${db_download_dir}" )

                BAKTA_BAKTADBDOWNLOAD ()

                //BAKTA_BAKTADBDOWNLOAD.out.db.toAbsolutePath().toString()

                BAKTA_UPDATE_CONFIG (
                    BAKTA_BAKTADBDOWNLOAD.out.db,
                    params.bakta_download_dir
                )

                bakta_db_dir = BAKTA_UPDATE_CONFIG.out.db // "${db_download_dir}/bakta/db"
                
            }
        } else {
            bakta_db_dir = []
        }
    }

    // Bin classification databases
    if ( !params.skip_classification ) {
        // If GTDB will be used
        if ( !params.skip_gtdbtk ) {
            // AND if there is a GTDB database given
            if ( params.gtdb_db ) {
                gtdb_db_dir = file ( params.gtdb_db, checkIfExists: true )

                if ( gtdb_db_dir ) {
                    println ( "GTDB database found at ${params.gtdb_db}" )
                } else {
                    exit 1 ( "ERROR: gtdb_db path ${params.gtdb_db} does not exist" )
                }

            // If no GTDB database is given
            } else {
                println ( "GTDB database not given, downloading to ${db_download_dir}" )

                GTDB_DB_DOWNLOAD (
                    params.gtdb_download
                )

                GTDB_UPDATE_CONFIG (
                    GTDB_DB_DOWNLOAD.out.db,
                    db_download_dir
                )

                gtdb_db_dir = GTDB_UPDATE_CONFIG.out.db
            }
        } else {
            gtdb_db_dir = []
        }
    } else {
        gtdb_db_dir = []
    }

    /****************
        Data input
     ****************/

    reads_input = channel.fromPath ( params.reads_file )
        .splitCsv (header: true)
        .map { meta ->
            def meta_new = [:]

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
    
    //reads_input.view()

    short_reads = channel.empty()
    short_reads = short_reads.mix ( reads_input )
        .map { meta, reads ->
            if ( meta.sequencer == 'Illumina' ) {
                return [ meta, reads ]
            } 
        }

    nanopore_reads = channel.empty()
    nanopore_reads = nanopore_reads.mix ( reads_input )
        .map { meta, reads ->
            if ( meta.sequencer == 'ONT' ) {
                return [ meta, reads ]
            } 
        }

    pacbio_reads = channel.empty()
    pacbio_reads = pacbio_reads.mix ( reads_input )
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
    //phiX = channel.empty()
    phiX_index = channel.empty()
    if ( !params.skip_qc && params.short_reads && params.remove_phiX ) {
        /*phiX = channel.fromPath ( params.phiX )
            .map { reference ->
                def meta = [:]
                meta.id = reference.getBaseName()
                return [ meta, reference ]
            }
        */
        phiX_index = channel.fromPath ( "${projectDir}/reference_genomes/phiX/*.bt2" )
            .map { index ->
                def meta =[:]
                meta.id = 'phiX174'
                return [ meta, index ]
            }
            .groupTuple()
    }

    host_genome = channel.empty()
    // Create channel of host genome to remove
    if ( !params.skip_qc && params.host_genome ) {
        host_genome = channel.fromPath ( params.host_genome )
            .map { reference ->
                def meta = [:]
                meta.id = reference.getBaseName()
                return [ meta, reference ]
            }
    }

    // Short read quality control
    corrected_reads = channel.empty()
    if ( !params.skip_qc && params.short_reads ) {
        QC_SHORT ( 
            short_reads,
            phiX_index,
            host_genome
        )

        corrected_reads = corrected_reads.mix ( QC_SHORT.out.host_filtered_reads )
    }

    // ONT quality control
    corrected_ont_reads = channel.empty()
    if ( !params.skip_qc && params.nanopore_reads ) {
        QC_NANOPORE ( 
            nanopore_reads,
            host_genome
        )
        
        corrected_ont_reads = corrected_ont_reads.mix ( QC_NANOPORE.out.filtered_nanopore_reads )
    }

    // PacBio quality control
    corrected_pacbio_reads = channel.empty()
    if ( !params.skip_qc && params.pacbio_reads ) {
        QC_PACBIO ( 
            pacbio_reads,
            host_genome
        )

        corrected_pacbio_reads = corrected_pacbio_reads.mix ( QC_PACBIO.out.filtered_pacbio_reads )
    }

    all_corrected_reads = corrected_reads.mix ( corrected_ont_reads, corrected_pacbio_reads )
    // if this is empty, exit

    //all_corrected_reads.view()

    /*******************
        Read taxonomy
     *******************/

    //read_taxonomy_input = channel.empty()
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
    concatenated_reads = channel.empty()        // Channel to concatenate fastq files for assembly ('per_sample' will be ungrouped sample fastq's)
    original_clean_reads = channel.empty()    // Channel to group original individual sample fastq's according to assembly file (for use in binning)
    if ( !params.skip_assembly && params.short_reads ) {
        ASSEMBLY_PREP_SHORT (
            corrected_reads
        )
        
        concatenated_reads = concatenated_reads.mix ( ASSEMBLY_PREP_SHORT.out.concatenated_reads )
        original_clean_reads = original_clean_reads.mix ( ASSEMBLY_PREP_SHORT.out.original_clean_reads )
    }

    concatenated_long_reads = channel.empty()       // Channel to concatenate fastq files for assembly ('per_sample' will be ungrouped sample fastq's)
    original_clean_long_reads = channel.empty()     // Channel to group original individual sample fastq's according to assembly file (for use in binning)
    if ( !params.skip_assembly && ( params.nanopore_reads || params.pacbio_reads ) ) {
        ASSEMBLY_PREP_LONG (
            corrected_ont_reads,
            corrected_pacbio_reads
        )
        
        concatenated_long_reads = concatenated_long_reads.mix ( ASSEMBLY_PREP_LONG.out.concatenated_long_reads )
        original_clean_long_reads = original_clean_long_reads.mix ( ASSEMBLY_PREP_LONG.out.original_clean_long_reads )
    }
    // *****************
    
    final_contigs = channel.empty()
    assembly_graphs = channel.empty()
    reads_post_assembly = channel.empty()
    hifiasm_bins = channel.empty()

    if ( !params.skip_assembly ) {
        if  ( !params.skip_spadeshybrid ) {
            if ( params.short_reads && (params.nanopore_reads || params.pacbio_reads )) {
                /*ASSEMBLY_HYBRID (
                    concatenated_reads, 
                    original_clean_reads,
                    concatenated_long_reads, 
                    original_clean_long_reads
                )

                final_contigs = ASSEMBLY_HYBRID.out.contigs
                assembly_graphs = assembly_graphs.mix ( ASSEMBLY_HYBRID.out.assembly_graph )
                reads_post_assembly = reads_post_assembly.mix ( ASSEMBLY_HYBRID.out.reads )
                */
            } else {
                exit 1 ("ERROR: Cannot perform hybrid assembly without both short and long reads")
            }

        } else if ( params.short_reads ) {
            ASSEMBLY_SHORT (
                concatenated_reads, 
                original_clean_reads
            )

            final_contigs = final_contigs.mix ( ASSEMBLY_SHORT.out.contigs )
            assembly_graphs = assembly_graphs.mix ( ASSEMBLY_SHORT.out.assembly_graph )
            reads_post_assembly = reads_post_assembly.mix ( ASSEMBLY_SHORT.out.reads )

        } else if ( params.nanopore_reads || params.pacbio_reads ) {
            ASSEMBLY_LONG (
                concatenated_long_reads, 
                original_clean_long_reads
            )

            final_contigs = final_contigs.mix ( ASSEMBLY_LONG.out.final_contigs )
            assembly_graphs = assembly_graphs.mix ( ASSEMBLY_LONG.out.assembly_graphs )
            reads_post_assembly = reads_post_assembly.mix ( ASSEMBLY_LONG.out.final_reads )
            hifiasm_bins = hifiasm_bins.mix ( ASSEMBLY_LONG.out.hifiasm_bins )
        }
    }

    //final_contigs.view()
    //assembly_graphs.view()
    //reads_post_assembly.view()

    /*********************
        Contig taxonomy
     *********************/
    tax_4_vamb = channel.empty()
    if ( !params.skip_contig_taxonomy ) {
        CONTIG_TAXONOMY (
            final_contigs, 
            kraken2_db_dir,
            run_bracken,
            "contigs"
        )

        tax_4_vamb = tax_4_vamb.mix ( CONTIG_TAXONOMY.out.tax_4_vamb )
    } else {
        // If contig taxonomy not performed, create a channel with empty placeholders
        tax_4_vamb = tax_4_vamb.mix ( final_contigs )
            .map { meta, _contigs ->
                [ meta, [] ]
            }
    }

    /*************
        Binning
     *************/

    binning_prep_input = channel.empty()
    initial_bins = channel.empty()
    post_refinement_bins = channel.empty()
    if ( !params.skip_binning ) {
        binning_prep_input = binning_prep_input.mix ( final_contigs )
            .join ( assembly_graphs )
            .join ( reads_post_assembly )
            .join ( tax_4_vamb )

        //binning_prep_input.view()                    

        BINNING_PREPARATION ( binning_prep_input )

        //BINNING_PREPARATION.out.grouped_mappings.view()

        BINNING (
            BINNING_PREPARATION.out.grouped_mappings,
            hifiasm_bins
        )

        //BINNING.out.bins.view()

        refinement_contigs = BINNING_PREPARATION.out.grouped_mappings
            .map { meta, _reads, contigs, _bams, _bais, _gfa, _tax ->
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

            BIN_REFINEMENT.out.refined_bins.view()

        } else {
            post_refinement_bins = post_refinement_bins.mix ( initial_bins )
        }
            
    }

    //*****************
    // DAStool and drep arguments
    
    /********************
        Bin evaluation
     ********************/

    bin_evaluations = channel.empty()
    bins_post_evaluation = channel.empty()

    // CheckM
    if ( !params.skip_bin_evaluation) {
        BIN_EVALUATION (
            post_refinement_bins,
            checkm2_db_dir,
            checkm_db_dir
        )

        bins_post_evaluation = bins_post_evaluation.mix ( BIN_EVALUATION.out.bins_output )
        bin_evaluations = bin_evaluations.mix( BIN_EVALUATION.out.bin_summary )

    } else {

        bins_post_evaluation = bins_post_evaluation.mix ( post_refinement_bins )

    }

    //bin_evaluations.view()

    /*******************
        Dereplication
     *******************/

    final_bins = channel.empty()

    if ( !params.skip_bin_dereplication ) {
        BIN_DEREPLICATION ( 
            bins_post_evaluation,
            bin_evaluations
        )

        final_bins = final_bins.mix ( BIN_DEREPLICATION.out.dereplicated_bins )

    } else {
        final_bins = final_bins.mix ( bins_post_evaluation )
    }

    /************************
        Bin classification
     ************************/

    if ( !params.skip_classification ) {
       BIN_CLASSIFICATION ( 
            final_bins,
            gtdb_db_dir
        )
    }

    /********************
        Bin annotation
     ********************/

    if ( !params.skip_annotation ) {
        BIN_ANNOTATION ( 
            final_bins,
            bakta_db_dir
        )
    }

    /******************
        Bin coverage
     ******************/
    
    if ( !params.skip_binning ) {
        bin_group_reads = BINNING_PREPARATION.out.grouped_mappings
            .map { meta, reads, _contigs, _bams, _bais, _gfa, _tax ->
                [ meta, reads ]
            }

        if ( !params.skip_bin_coverage ) {
            BIN_COVERAGE (
                all_corrected_reads,        // Calculate coverage for each original read sample...
                bin_group_reads,            // ... as well as concatenated reads used in each binning instance (if applicable)
                final_bins                  
            )
        }
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
    
/*    ch_quast_bins_summary = channel.empty()
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
            ch_cat_global_summary = channel.empty()
        //}

    /*
     * GTDB-tk: taxonomic classifications using GTDB reference
     */
    
    /*if (!params.skip_gtdbtk) {

        ch_gtdbtk_summary = channel.empty()
        if (gtdb) {

            GTDBTK(
                final_bins,
                ch_bin_qc_summary,
                gtdb,
                gtdb_mash
            )
            ch_versions = ch_versions.mix(GTDBTK.out.versions.first())
            ch_gtdbtk_summary = GTDBTK.out.summary
        }
    }
    else {
        ch_gtdbtk_summary = channel.empty()
    }*/
    /*
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

    workflow.onComplete = {
        // any workflow property can be used here
        println "<project> complete"
        println "Command line: $workflow.commandLine"
    }

    workflow.onError = {
        println "Error: something went wrong, please see nextflow.log for details"
    }
    
}

