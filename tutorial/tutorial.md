# MAGnolia tutorial

### Installation

You will need to create a conda environment with Nextflow installed (if you don't already have one), and activate it:

```
conda create -n magnolia bioconda::nextflow
conda activate magnolia
```

Install MAGnolia by cloning the github repository. In this tutorial, we will use the mock pathway /home/username/ as the parent directory for both MAGnolia and our analysis (replace this with the correct paths on your server):
```
cd /home/username/
git clone https://github.com/danross509/MAGnolia.git
```

### Database setup
Upon first installing MAGnolia on your server, please ensure you open the database configuration file `./MAGnolia/configs/databases.config` and identify any relevant databases currently installed. Any required databases (according to the configuration settings) will be downloaded when running MAGnolia's main workflow, and the databases.config file will be updated automatically. Please see the provided example `MAGnolia/tutorial/example_databases.config`

### Data download
Create an empty directory for the tutorial, and for the tutorial dataset:
```
mkdir -p magnolia_tutorial/data
cd magnolia_tutorial/data
```
We will use the sample metagenomic data from the metaHIT gut survey for the tutorial analysis
```
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011347/ERR011347_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011348/ERR011348_2.fastq.gz

wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_1.fastq.gz
wget ftp.sra.ebi.ac.uk/vol1/fastq/ERR011/ERR011349/ERR011349_2.fastq.gz
```

### MAGnolia setup
With this done, ensuring you have your Nextflow (`magnolia`) environment activated, return to the tutorial directory and begin running MAGnolia:
```
cd ..       #to return to /home/username/magnolia_tutorial
nextflow run ../MAGnolia/setup.nf --illumina data
```
MAGnolia can be run using either the relative or the absolute path `/home/username/MAGnolia/setup.nf`. See README.md for additional setup parameters. For this tutorial, we will run MAGnolia without a GPU.

This script will generate 2 primary files - `samples.csv` and `nextflow.config` - and the folder `configs` (containing additional configuration files which can be modified for advanced useage) into the launch directory.

`samples.csv` is a list of all samples, with default assembly and binning parameters. Please check it to confirm the samples were correctly identified, please see the provided example `MAGnolia/tutorial/example_samples.csv`.

`nextflow.config` is the main file you will use to select which steps are run and choose which programs are used, as well as modify any other major parameters if needed. Please see the provided example `MAGnolia/tutorial/example_nextflow.config`. 

For this tutorial, we will be using default settings**. However, ensure that the `maxMem` and `maxCores` parameters reflect what is available on your server - we recommend setting these at least ~3-5GB and 1-2 cores below the actual maximum values. 

**The default workflow requires the following databases: `Kraken2`, `Bracken` (if maxMem found to be <240.GB), `CheckM2`, `Bakta`, and `GTDB`. Bin evaluation metrics from `CheckM` or `CheckM2` are required for dereplication - if MAGnolia's bin_evaluation is skipped, dRep will install its own CheckM. For the purpose of brevity, if you do not have these databases previously installed and would like to avoid downloading them during the tutorial, please modify the relevant parameters indicated below:

```
Kraken2 (this also skips Bracken)
    23 skip_read_taxonomy = true
        AND
    28 skip_contig_taxonomy = true

Bracken (if not skipping Kraken2)
    95 skip_bracken = true

Bakta
    37 skip_annotation = true
        OR
    156 skip_bakta = true (if running DRAM or eggNOG instead)

GTDB
    36 skip_classification = true
        OR
    161 skip_gtdbtk = true

CheckM2
    34 skip_bin_evaluation = true
        Or if you prefer to use CheckM instead, specify:
    152 checkm_version = 'checkm' 

dRep (to avoid installing CheckM if skipping bin_evaluation)
    33 skip_bin_dereplication = true
```

### MAGnolia main workflow
Once the setup is complete, again ensuring you have your Nextflow (`magnolia`) environment activated, and you are in the `magnolia_tutorial` directory containing `samples.csv` and `nextflow.config`, run the main workflow:
```
nextflow run ../MAGnolia/main.nf [--resume]
```
The `--resume` parameter is not necessary, however in case there is an error or you need to restart MAGnolia, is will ensure you are not rerunning programs that have already finished.

The results will be published in the `magnolia_tutorial` directory.

### Results
For reference, we have provided the results of this demo analysis (excluding the Nextflow `work` directory) in the folder `MAGnolia/tutorial/output/`. Please note that several files have been truncated to accomodate Github's file size limit, as indicated in `MAGnolia/output/truncated_files.txt`.