#!/usr/bin/env python3

# Originally written by David Ross for use within __.
# See git repository (https://github.com/) for full license text.

# USAGE: ./summarize_bracken_abundances.py -n $file_type -r $abundances 

import sys
import argparse
import os.path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



parser = argparse.ArgumentParser(
    prog='summarize_bracken_abundances',
    description='Generate a sample x taxa matrix of relative abundances generated from kraken2 and bracken results',
    epilog='')

parser.add_argument('-n', '--name')                             # file type [reads, contigs]. Must be identical to string used in abundance estimation
parser.add_argument('-l', '--level')                            # Taxonomy level used in Bracken abundance estimation
parser.add_argument('-0', '--na_rep')                          # Empty cell representation (eg. 0, NA, NaN)
parser.add_argument('-v', '--verbose')    
parser.add_argument('-r', '--reports', nargs='+', default="")   # list of sample abundance reports

args = parser.parse_args()

# Create output filename
OUTPUT_CSV = f"{args.name}_taxa_{args.level}_abundance_summary.csv"

#reports = args.reports.strip('[]').split(', ')

df = pd.DataFrame()
all_taxa = set()

sample_data = {}

for file in args.reports:
    path = file.strip('[,]')
    sampleID = path.split('/')[-1].replace(f'_{args.name}.bracken', "")
        
    try:
        sample_df = pd.read_csv(path, sep='\t')
            
        # Check if required column exists
        if 'fraction_total_reads' not in sample_df.columns:
            print(f"Warning: 'fraction_total_reads' not found in {path}")
            continue
    
        # Create dictionary of taxon -> fraction_total_reads for this sample
        sample_taxa = dict(zip(sample_df['name'], sample_df['fraction_total_reads']))
        sample_data[sampleID] = sample_taxa
            
        # Add taxa to the set of all taxa
        all_taxa.update(sample_taxa.keys())
            
    except Exception as e:
        print(f"Error processing {path}: {str(e)}")
        continue
    
if not sample_data:
    raise ValueError("No valid Bracken files were processed")
    
# Create the final dataframe
# Initialize with NaN values
df = pd.DataFrame(index=list(sample_data.keys()), 
                            columns=sorted(all_taxa), 
                            dtype=float)
    
# Fill in the values
for sampleID, taxa_dict in sample_data.items():
    for taxon, abundance in taxa_dict.items():
        df.loc[sampleID, taxon] = abundance
    
# Sort the dataframe by sample ID for consistency
df = df.sort_index()
    
df.to_csv(OUTPUT_CSV, index=True, index_label=args.name, na_rep=args.na_rep)