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
-"samples.csv": a list of all samples, with default parameters
-"nextflow.config": a detailed configuration sheet (in json format), to be reviewed and modified if necessary before running the main script :

```
nextflow run ../main.nf -resume
```
Which will by default run Quality Control, Assembly, Binning, Classification, and Annotation on all samples unless otherwise modified in nextflow.config

 # Data Input
 ## Short reads

 Each fastq.gz file is expected to represent an individual sample. 

 ## Nanopore reads

 For now, only fastq.gz format is supported. Please ensure your reads are demultiplexed and basecalled. 
 
 setup.nf will look for fastq.gz files contained within barcode folders in the designated `$path`, then concatenate them into a single file named `${barcode}.fastq.gz`, where `${barcode}` is the name of the barcode folder. This is intended to reduce user preprocessing, using the raw Nanopore output format. If you have already have concatenated `${barcode}.fastq.gz` files, they should be placed each within their own folder (appropriately named). 

 ## PacBio reads

 TO DO

 

