# MAG_Pipeline
 This pipeline is intended for the recovery of Metagenome-AssembledGenomes (MAGs) using short reads, long reads or a hybrid approach.

 Once installed, it is designed to be run using two commands:

``` 
nextflow run ../setup.nf [--illumina $path] [--nanopore $path] [--pacbio $path] [--corrected]
```
Where `$path` identifies the directory containing, for each input type:

```
--illumina:     Specify folder containing illumina (short read) fastq.gz files. E.g. `${path}/*{1,2}.fastq.gz`

--nanopore:     Specify folder containing nanopore barcode folders. E.g. `${path}/{barcodes}/*.fastq.gz`

--pacbio:       Specify folder containing PacBio reads... TO DO

--corrected:    Specify if all input reads have been previously corrected, otherwise can be edited in samples.csv
```

This will generate 2 files in the launch directory:

-`samples.csv`: a list of all samples, with default parameters
```
sampleID    sequencer   paired_end  corrected   bin_group   assembly_group  reads_R1    reads_R2
```

-`nextflow.config`: a detailed configuration sheet (in json format), to be reviewed and modified if necessary before running the main script :

```
nextflow run ../main.nf [-resume]
```
Which will by default run Quality Control, Assembly, Binning, Classification, and Annotation on all samples unless otherwise modified in nextflow.config

# Data Input
## Short reads

Each fastq.gz file is expected to represent an individual sample. 

## Nanopore reads

For now, only fastq.gz format is supported. Please ensure your reads are demultiplexed and basecalled. 
 
Setup will look for fastq.gz files contained within barcode folders in the designated `$path`, then concatenate them into a single file named `${barcode}.fastq.gz`, where `${barcode}` is the name of the barcode folder. This is intended to reduce user preprocessing, using the raw Nanopore output format. If you have already have concatenated `${barcode}.fastq.gz` files, they should be placed each within their own folder (appropriately named). 

## PacBio reads

TO DO

# Sample management
To assemble and bin each sample on its own, the `assembly_group` and `bin_group` parameters should match the corresponding unique `sampleID`. 

## Assembly

____ supports *per sample*, *grouped sample*, and all sample *co-assembly*.

This should be adjusted to suit your analysis in the generated `samples.csv` file.
 
### Per sample
This is the default assembly method for 4+ samples, where the `assembly_group` parameter is set identical to `sampleID`. Each sample will be assembled individually. 

In order to perform either `grouped sample` or all sample `co-binning` in the following step, `assembly_group` MUST be identical to `sampleID`.

### Co-assembly
Setup will default to using a *co-assembly* method for 3 or fewer samples, with `assembly_group` set as "allReads". Regardless of sample number, you can manually set `assembly_group` of all samples to a unique name in order to co-assemble them. DO NOT USE AN `assembly_group` NAME THAT IS IDENTICAL TO ANY `sampleID`.

This will concatenate all reads into a single (or paired) fastq file(s) before assembly. It can be computationally demanding depending on the assembler. It can help to identify rare taxa, but at the cost of contaminating common ones. 

Only *per assembly* binning is available for *co-assembled* samples.

### Grouped sample
This method is essentially the same as *co-assembly*, but with user-specified subgroups. Use unique `assembly_group` names to identify the samples you want assembled together. DO NOT USE `assembly_group` NAMES THAT ARE IDENTICAL TO ANY `sampleID`.

This method is intended to be able to pass, essentially, multiple *co-assemblies* through _____ in parallel.

Only *per assembly* binning is available for *co-assembled* samples.


## Binning

Likewise, ____ supports *per assembly*, *grouped sample*, and all sample *co-binning*. 

Before binning, `sampleID` will be replaced with `bin_group` in order to concatenate the appropriate contigs together.

Again, we recommend adjusting this to suit your analysis in the `samples.csv` file.

### Per assembly
This is the default binning method, and the only method available for any group of *co-assembled* samples. 

By default, `bin_group` for each sample will be set identical to `sampleID`, and will be changed to match the new `sampleID` of any *co-assembly* contig files (ie. identical to `assembly_group`). Again, for *per sample* assembly and binning, ensure that `assembly_group` and `bin_group` match the unique corresponding `sampleID`.

### Co-binning
This method is available only for *per sample* assembled samples. 

Set `bin_group` for all samples to a unique name in order to *co-bin* them. DO NOT USE A `bin_group` NAME THAT IS IDENTICAL TO ANY `sampleID`.

This will concatenate all assembled contigs into a single fasta before binning. It can be a less computationally demanding alternative to *co-assembly* to help to identify rare taxa (again, at the cost of contaminating common ones). 

### Grouped sample
This method is also available only for *per sample* assembled samples. 

This method is essentially the same as *co-binning*, but with user-specified subgroups. Use unique `bin_group` names to identify the samples you want binned together. DO NOT USE `bin_group` NAMES THAT ARE IDENTICAL TO ANY `sampleID`.

This method is intended to be able to pass, essentially, multiple *co-binning* groups through _____ in parallel.

