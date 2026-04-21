#!/usr/bin/env python3

## Originally written by David Ross for use within MAGnolia

# USAGE: ./write_samples_csv.py $sampleID -s $sequencer -p $paired_end -c $corrected -b $bin_group -a $assembly_group -1 $reads_1 -2 $reads_2

import argparse
import os

OUTPUT_CSV = "samples_tmp.csv"

parser = argparse.ArgumentParser(
                    prog='write_samples_csv',
                    description='Use input folders to identify samples and metadata',
                    epilog='')

parser.add_argument('sampleID')                     # positional argument
parser.add_argument('-s', '--sequencer')            # sequencer used
parser.add_argument('-p', '--paired_end')           # are reads paired?
parser.add_argument('-c', '--corrected')            # are reads corrected?
parser.add_argument('-a', '--assembly_group')       # assembly group ID
parser.add_argument('-b', '--bin_group')            # bin group ID
parser.add_argument('-1', '--reads_1')              # path to illumina reads R1, and to nanopore and pacbio reads
parser.add_argument('-2', '--reads_2')              # path to paired illumina reads R2
#parser.add_argument('-#', '--sample_count')         # number of samples of given sequencer
parser.add_argument('-v', '--verbose')    

args = parser.parse_args()

# Open sample file, write line 

with open(OUTPUT_CSV, "a", encoding="utf8") as sample_file:
    sample_file.write(f"{args.sampleID},{args.sequencer},{args.paired_end},{args.corrected},{args.assembly_group},{args.bin_group},{args.reads_1},{args.reads_2}\n")

