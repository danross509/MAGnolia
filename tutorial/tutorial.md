# MAGnolia tutorial

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

Upon first installing MAGnolia on your server, please ensure you open the database configuration file `./MAGnolia/configs/databases.config` and identify any relevant databases currently installed. Any required databases (according to the configuration settings) will be downloaded, and the databases.config file will be updated automatically. Please see the provided example `MAGnolia/tutorial/example_databases.config`

Create an empty directory for the tutorial, and for the toy dataset:
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
With this done, ensuring you have the `magnolia` environment activated, return to the tutorial directory and begin running MAGnolia:
```
cd ..       #to return to /home/username/magnolia_tutorial
nextflow run ../MAGnolia/setup.nf --illumina data
```
MAGnolia can be run using either the relative or the absolute path `/home/username/MAGnolia/setup.nf`. See README.md for additional setup parameters. For this tutorial, we will run MAGnolia without a GPU.

This script will output 2 primary files...
```
/home/username/magnolia_tutorial/samples.csv
/home/username/magnolia_tutorial/nextflow.config
```
... and the folder `/home/username/magnolia_tutorial/configs`, containing additional configuration files which can be modified for advanced useage.

`samples.csv` is a list of all samples, with default assembly and binning parameters. Please open it to confirm the samples were correctly identified.
