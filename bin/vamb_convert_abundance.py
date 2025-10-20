#!/usr/bin/env python3

import sys
import argparse
import pandas as pd

parser = argparse.ArgumentParser(
                    prog='vamb_convert_abundance',
                    description='Convert metabat2 jgisummarizebamcontigdepths output to be compatible with Vamb v5.0.4',
                    epilog='')

parser.add_argument('-i', '--input') #Input tsv
parser.add_argument('-o', '--output') #Input tsv
    
args = parser.parse_args()

INPUT_FILE = args.input
OUTPUT_FILE = args.output


df = pd.read_csv(INPUT_FILE, sep='\t')
    
columns = df.columns.tolist()
    
abundance_cols = [col for col in columns 
                   if col.endswith('.sorted.bam') and not col.endswith('-var')]

sample_names = [col.replace('.sorted.bam', '') for col in abundance_cols]
    
output_df = pd.DataFrame()
output_df['contigname'] = df['contigName']
    
# Add sample columns with cleaned names
for abundance_col, sample_name in zip(abundance_cols, sample_names):
    output_df[sample_name] = df[abundance_col]
    
# Write to output file
output_df.to_csv(OUTPUT_FILE, sep='\t', index=False)

