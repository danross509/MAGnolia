#!/usr/bin/env nextflow

/*
To run this pipeline during development:
From the MAG_Pipeline/pipeline_test/ folder, run the following command:
    nextflow run ../main.nf -params-file ../params.json -resume



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

include { QC } from './QC/main.nf'
/*
The following modules are currently in development:
include { assembly } from './Assembly/main.nf'
include { binning } from './Binning/main.nf'
include { taxonony } from './Taxonomy/main.nf'
include { annotation } from './Annotation/main.nf'
*/

workflow {
    // Specify the output directory
    // projectDir = Channel.fromPath(params.projectDir, type: 'dir')

    /*
     * Import quality control parameters
     * If paired-read filenames are identified by R1/R2, the R will be removed for simplicity
     * PE and SE reads will be identified and handled as necessary for each software
     */
    short_reads = Channel.fromFilePairs(params.short_reads, size: -1, checkIfExists:true)
        .map { sample, reads ->
            def sampleID = sample.replaceAll(/_R$/, '')
            return [sampleID, reads]
        }
        .map { sampleID, reads ->           //Inspired by nf-core components meta-map
            def meta = [:]
            // Set meta.id
            meta.id = sampleID
            // Set meta.paired_end
            if (reads.size() == 2) {
                meta.paired_end = true
            } else if (reads.size() == 1) {
                meta.paired_end = false
            } else {
                exit 1, "ERROR: Check input files-> ${sampleID} contains invalid number of samples"
            }
            // Set bin-group
            meta.bin_group = "all"
            return [ meta, reads ]
        }
    //short_reads.view()
    
    phiX = params.phiX
    bowtie2_sensitivity = params.bowtie2_sensitivity
    paired_reads = params.paired_reads
    coassembly = params.coassembly

    QC( short_reads,
        phiX,
        bowtie2_sensitivity,
        paired_reads,
        coassembly
        )

    //QC.out[1].collect().view()




    /*
     *  Collect clean reads from QC, assemble
     */
    clean_reads = QC.out
    mh_preset = params.mh_preset

    /*assembly(   clean_reads,
                mh_preset
                )*/
    
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

