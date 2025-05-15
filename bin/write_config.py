#!/usr/bin/env python

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./write_config.py -s $short_reads_count -p $short_reads_paired -n $nanopore_barcodes_count -pb $pacbio_reads_count -d $databases_location -c $reads_corrected

import argparse

#OUTPUT_JSON = "samples.json"
INPUT_CONFIG = "../default.config"
OUTPUT_CONFIG = "../nextflow.config"

parser = argparse.ArgumentParser(
                    prog='write_config',
                    description='Write nextflow.config filebased on samples and metadata',
                    epilog='')

parser.add_argument('-s', '--short_count')      # Number of short read files
parser.add_argument('-p', '--short_paired')     # Are short reads paired?
parser.add_argument('-n', '--ont_count')        # Number of nanopore barcodes
parser.add_argument('-pb', '--pacbio_count')    # Number of PacBio reads
parser.add_argument('-c', '--corrected')        # Are reads corrected
parser.add_argument('-v', '--verbose',
                    action='store_true')        # on/off flag

args = parser.parse_args()

# Read default nextflow.config
with open(INPUT_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the required parameters
if int(args.short_count) > 0 :
    text = text.replace("short_reads = false", "short_reads = true")
if args.short_paired == 'true' :
    text = text.replace("paired_short_reads = false", "paired_short_reads = true")
if int(args.ont_count) > 0 :
    text = text.replace("nanopore_reads = false", "nanopore_reads = true")
if int(args.pacbio_count) > 0 :
    text = text.replace("pacbio_reads = false", "pacbio_reads = true")

if args.corrected == 'true':
    if int(args.short_count) > 0 :
        text = text.replace("short_reads_corrected = false", "short_reads_corrected = true")
    if int(args.ont_count) > 0 :
        text = text.replace("nanopore_reads_corrected = false", "nanopore_reads_corrected = true")
    if int(args.pacbio_count) > 0 :
        text = text.replace("pacbio_reads_corrected = false", "pacbio_reads_corrected = true")

if 0 <= int(args.short_count) <= 3 and 0 <= int(args.ont_count) <= 3 and 0 <= int(args.pacbio_count) <= 3 :
    text = text.replace("assembly_mode = 'per_sample'", "assembly_mode = 'coassembly'")

if int(args.short_count) > 0 and (int(args.ont_count) > 0 or int(args.pacbio_count) > 0) :
    text = text.replace("skip_spadeshybrid = true", "skip_spadeshybrid = false")

# Write the file out again
with open(OUTPUT_CONFIG, "w") as file:
  file.write(text)