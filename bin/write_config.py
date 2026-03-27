#!/usr/bin/env python3

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./write_config.py -s $short_reads_count -p $short_reads_paired -n $nanopore_barcodes_count -pb $pacbio_reads_count -d $databases_location -c $reads_corrected

import argparse

parser = argparse.ArgumentParser(
                    prog='write_config',
                    description='Write nextflow.config filebased on samples and metadata',
                    epilog='')

parser.add_argument('-s', '--short_count')      # Number of short read files
parser.add_argument('-p', '--short_paired')     # Are short reads paired?
parser.add_argument('-n', '--ont_count')        # Number of nanopore barcodes
parser.add_argument('-pb', '--pacbio_count')    # Number of PacBio reads
parser.add_argument('-c', '--corrected')        # Are reads corrected
parser.add_argument('-f', '--default_file')     # Default config file path
parser.add_argument('-g', '--use_gpu')          # Use an available GPU
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

INPUT_CONFIG = args.default_file
OUTPUT_CONFIG = "nextflow.config"

# Read default nextflow.config
with open(INPUT_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the required parameters:
# Based on read presence / absence
if int(args.short_count) > 0 :
    text = text.replace("short_reads = false", "short_reads = true")
    if args.use_gpu == 'true': 
        text = text.replace("skip_semibin2 = true", "skip_semibin2 = false")
    else:
        text = text.replace("skip_maxbin2 = true", "skip_maxbin2 = false")
if args.short_paired == 'true' :
    text = text.replace("paired_short_reads = false", "paired_short_reads = true")

if int(args.ont_count) > 0 :
    text = text.replace("nanopore_reads = false", "nanopore_reads = true")
    text = text.replace("skip_contig_polishing = true", "skip_contig_polishing = false")

if int(args.pacbio_count) > 0 :
    text = text.replace("pacbio_reads = false", "pacbio_reads = true")

if int(args.ont_count) > 0 or int(args.pacbio_count) > 0 :
    text = text.replace("skip_lrbinner = true", "skip_lrbinner = false")
    if args.use_gpu == 'true': 
        text = text.replace("skip_semibin2 = true", "skip_semibin2 = false")

# Read correctedness
if args.corrected == 'true':
    if int(args.short_count) > 0 :
        text = text.replace("short_reads_corrected = false", "short_reads_corrected = true")
    if int(args.ont_count) > 0 :
        text = text.replace("nanopore_reads_corrected = false", "nanopore_reads_corrected = true")
    if int(args.pacbio_count) > 0 :
        text = text.replace("pacbio_reads_corrected = false", "pacbio_reads_corrected = true")

# Hybrid assembly
if int(args.short_count) > 0 and (int(args.ont_count) > 0 or int(args.pacbio_count) > 0) :
    text = text.replace("skip_spadeshybrid = true", "skip_spadeshybrid = false")

if args.use_gpu == 'true': 
    text = text.replace("use_gpu = false", "use_gpu = true")

# Write the file out again
with open(OUTPUT_CONFIG, "w") as file:
  file.write(text)