#!/usr/bin/env python

## Originally written by David Ross for use within __
## See git repository (https://github.com/) for full license text.

# USAGE: ./create_samples_csv.py $sampleID -s $sequencer -p $paired_end -c $corrected -b $bin_group -a $assembly_group -1 $reads_1 -2 $reads_2

OUTPUT_CSV = "samples_tmp.csv"

with open(OUTPUT_CSV, "w", encoding="utf8") as sample_file:
    sample_file.write(f"sampleID,sequencer,paired_end,corrected,assembly_group,bin_group,reads_R1,reads_R2\n")

