# MAGnolia
 This pipeline is intended for the recovery of Metagenome-Assembled Genomes (MAGs) using short reads, long reads or a hybrid approach.

## Installation

MAGnolia requires a conda environment with Nextflow installed (written in v25.10.4):

```
conda create -n magnolia bioconda::nextflow

conda activate magnolia
```

MAGnolia is installed by cloning the github repository:
```
git clone https://github.com/danross509/MAGnolia.git
```
## Databases
When first installing MAGnolia on a new server, please ensure you open the database configuration file `./MAGnolia/configs/databases.config` and identify any relevant databases currently installed on the server, to avoid installing duplicates. This file will serve as the default database configuration for all analyses using this installation of MAGnnolia; however, alternate database configuration files may be specified in your nextflow.config file. 

MAGnolia will download any databases required for the analysis, based on the settings in the main configuration file, that are not properly specified in the database configuration. Upon downloading a database, MAGolia will update the databases.config file.

## Running MAGnolia
MAGnolia uses the launch directory as the default. We recommend creating a new directory for each dataset (see tutorial.md)

It is designed to be run using two commands, where `$PATH` is the path to the cloned repository, and `$path` identifies the path to the directory containing each input type:

``` 
nextflow run $PATH/MAGnolia/setup.nf [--illumina $path] [--nanopore $path] [--pacbio $path] [--corrected] [--coassembly] [--cobinning] [--use_gpu]

Parameters:
--illumina:     Specify the folder containing illumina (short read) data: *_{R1,R2}.fastq.gz or *_{1,2}.fastq.gz

--nanopore:     Specify the folder containing nanopore data. Barcode subfolders will be concatenated if a corresponding fastq.gz is not found: ${path}/barcode/*.fastq.gz -> ${path}/barcode.fastq.gz

--pacbio:       Specify folder containing PacBio data, either in BAM or fastq.gz format. BAM files will be converted if a corresponding fastq.gz file is not found: ${path}/*.hifi_reads.bam -> ${path}/*.fastq.gz

--corrected:    Specify if all input reads have been previously corrected. This can be edited later in samples.csv

--coassembly:   Set all samples to be coassembled

--cobinning:    Set all samples to be cobinned

--use_gpu:      Programs with GPU compatibility will access a GPU
```

This will generate 2 files in the launch directory:

-`samples.csv`: a list of all samples, with default parameters
```
sampleID    sequencer   paired_end  corrected   assembly_group  bin_group  reads_R1    reads_R2
```

-`nextflow.config`: a detailed configuration sheet (in json format), to be reviewed and modified if necessary before running the main script :

```
nextflow run $PATH/MAGnolia/main.nf [-resume]
```
Which will by default run Quality Control, Assembly, Binning, Classification, and Annotation on all samples unless otherwise modified in nextflow.config

### Configuration options
We recommend reviewing the following parameters found in your nextflow.config file, and modifying them where necessary:

    -maxMem                 # The maximum RAM available for MAGnolia
    -maxCores               # The maximum number of cores
    -use_gpu                # Is there a GPU available?
    -host_genome            # Ifentify a host genome you would like to filter out
    -min_completeness       # Minimum completeness of final bins, to be used by multiple programs
    -max_contamination      # Maximum contamination of final bins
    -drep_sa                # Secondary ANI (%) for dereplication

# Data Input
### Short reads

Each fastq.gz file is expected to represent an individual sample. Paired and unpaired reads are supported. Interleaved files have not yet been tested.  

### Nanopore reads

Only fastq.gz format is supported, ensure your reads are basecalled and demultiplexed. 
 
Setup will look for fastq.gz files contained within the designated `$path`, as well as subfolders containing the demultiplexed fastq.gz files of each barcode. The files in each barcode subfolder `${path}/barcode/*.fastq.gz` will be concatenated into a single file named `${path}/barcode.fastq.gz`, where `barcode` is the name of the barcode folder, to reduce user preprocessing. Barcode subfolders will be ignored if a corresponding `${path}/barcode.fastq.gz` is present. 

### PacBio reads

Both BAM and fastq.gz formats are supported. 
 
Setup will look for BAM and fastq.gz files contained within the designated `$path`. BAM files will be converted to fastq.gz using pbtk_bam2fastq: `${path}/*.hifi_reads.bam` -> `${path}/*.fastq.gz`. BAM files will be ignored if a corresponding `${path}/*.fastq.gz` is present. 

# Sample management
To assemble and bin each sample on its own, the `assembly_group` and `bin_group` parameters should match the corresponding unique `sampleID`. 

## Assembly

MAGnolia supports *per sample*, *grouped sample*, and all sample *co-assembly*.

This should be adjusted to suit your analysis in the generated `samples.csv` file.
 
### Per sample
This is the default assembly method unless otherwise specified, where the `assembly_group` parameter is set identical to `sampleID`. Each sample will be assembled individually. 

In order to perform either grouped sample or all sample `co-binning` in the following step, `assembly_group` must be identical to `sampleID`.

### Co-assembly
Setup will set the `assembly_group` parameter as "allReads" if `--coassembly` is specified. You can manually set `assembly_group` of all samples to a unique name in order to co-assemble them. DO NOT USE AN `assembly_group` NAME THAT IS IDENTICAL TO ANY `sampleID`.

This will concatenate all reads into a single (or paired) fastq file(s) before assembly. It can be computationally demanding depending on the assembler. It can help to identify rare taxa, but at the cost of contaminating common ones. 

Only *per assembly* binning is available for *co-assembled* samples.

### Grouped sample
This is the same as *co-assembly*, but with user-specified subgroups. Use unique `assembly_group` names to identify the samples you want assembled together. DO NOT USE `assembly_group` NAMES THAT ARE IDENTICAL TO ANY `sampleID`.

This method will pass multiple *co-assemblies* through MAGnolia in parallel.

Only *per assembly* binning is available for *co-assembled* samples.


## Binning

Likewise, MAGnolia supports *per assembly*, *grouped sample*, and all sample *co-binning*. 

Again, we recommend adjusting this to suit your analysis in the `samples.csv` file.

### Per assembly
This is the default binning method, and the only method available for any group of *co-assembled* samples. 

By default, `bin_group` for each sample will be set identical to `sampleID`, and will be changed to match the new `sampleID` of any *co-assembled* contig files (ie. identical to `assembly_group`). Again, for *per sample* assembly and binning, ensure that `assembly_group` and `bin_group` match the unique corresponding `sampleID`.

### Co-binning
This method is available only for *per sample* assembled samples. 

Setup will set the `bin_group` parameter as "allReads" if `--cobinning` is specified. You can also manually set `bin_group` for all samples to a unique name in order to *co-bin* them. DO NOT USE A `bin_group` NAME THAT IS IDENTICAL TO ANY `sampleID`.

This will concatenate all assembled contigs into a single fasta before binning. It can be a less computationally demanding alternative to *co-assembly* to help to identify rare taxa (again, at the cost of contaminating common ones). 

### Grouped sample
This method is also available only for *per sample* assembled samples. 

This is the same as *co-binning*, but with user-specified subgroups. Use unique `bin_group` names to identify the samples you want binned together. DO NOT USE `bin_group` NAMES THAT ARE IDENTICAL TO ANY `sampleID`.

This method will pass multiple *co-binning* groups through MAGnolia in parallel.

