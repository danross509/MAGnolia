#!/usr/bin/env python3

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./write_config.py -s $short_reads_count -p $short_reads_paired -n $nanopore_barcodes_count -pb $pacbio_reads_count -d $databases_location -c $reads_corrected

import argparse
import re
import subprocess

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

# Find local max cores
def get_lscpu_value(key):
    result = subprocess.check_output(
        f"lscpu | grep -E '^{key}' | awk '{{print $NF}}'",
        shell=True
    ).strip().decode()
    return int(result)

cpus = get_lscpu_value("CPU\(s\):")
threads_per_core = get_lscpu_value("Thread\(s\) per core:")
cores_per_socket = get_lscpu_value("Core\(s\) per socket:")
sockets = get_lscpu_value("Socket\(s\):")
expected = threads_per_core * cores_per_socket * sockets

if cpus == expected:
    if cpus < 4:
        maxcores = cpus
    else:
        maxcores = cpus - 2
    print(f"Setting the expected maxCores value to {maxcores} cores")
else:
    physcores = cores_per_socket * sockets
    if physcores < 4:
        maxcores = physcores
    else:
        maxcores = physcores - 2
    print(f"We found a hybrid CPU architecture, therefore will set the expected maxCores value to {maxcores} physical cores")

# Find local max memory
foundmem = None
with open("/proc/meminfo") as f:
    for line in f:
        if line.startswith("MemTotal:"):
            kb = int(re.search(r'\d+', line).group())
            gb = kb / (1024 ** 2)
            foundmem = int(gb)

if foundmem is None:
    raise ValueError("maxmem could not be assigned: MemTotal not found in /proc/meminfo")
elif not isinstance(foundmem, int):
    raise TypeError(f"maxmem is not an integer: got {type(foundmem).__name__}")
else:
    if foundmem < 16:
        maxmem = int(foundmem * 0.9)
    else:
        maxmem = int(foundmem - 6)
    print(f"Setting the expected maxMem value to {maxmem}.GB")

# Read default nextflow.config
with open(INPUT_CONFIG, "r", encoding="utf8") as file:
    text = file.read()

# Replace the required parameters:
text = text.replace("maxMem = 124.GB", f"maxMem = {maxmem}.GB")
text = text.replace("maxCores = 16", f"maxCores = {maxcores}")

if maxmem < 240:
    text = text.replace("skip_bracken = false", "skip_bracken = true")

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