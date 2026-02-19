#!/usr/bin/env python3
"""

"""

import argparse
import sys
import pandas as pd


def standardize_and_convert(input_file, output_file):
    df = pd.read_csv(input_file, sep='\t')
        
    rename_map = {}
        
    if 'Bin Id' in df.columns:
        rename_map['Bin Id'] = 'genome'
    elif 'Name' in df.columns:
        rename_map['Name'] = 'genome'
    else:
        print("Warning: Neither 'Bin Id' nor 'Name' column found in input file", 
                file=sys.stderr)
        
    if 'Completeness' in df.columns:
        rename_map['Completeness'] = 'completeness'
    else:
        print("Warning: 'Completeness' column not found in input file", 
                file=sys.stderr)
        
    if 'Contamination' in df.columns:
        rename_map['Contamination'] = 'contamination'
    else:
        print("Warning: 'Contamination' column not found in input file", 
                file=sys.stderr)
        
    df.rename(columns=rename_map, inplace=True)

    # nf-core's dRep module inports the fastas into a folder "input_fastas"
    if 'genome' in df.columns:
        df['genome'] = df['genome'].astype(str) + '.fa'
        
    df.to_csv(output_file, index=False)


def main():
    parser = argparse.ArgumentParser(
        description='Convert checkm or checkm2 output for drep input',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s input.tsv output.csv
        """
    )
    
    parser.add_argument('input', help='Input TSV file')
    parser.add_argument('output', help='Output CSV file')
    
    args = parser.parse_args()
    
    standardize_and_convert(args.input, args.output)


if __name__ == '__main__':
    main()