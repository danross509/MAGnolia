#!/usr/bin/env python3
"""
Concatenate TSV files with modified contig names.
Extracts meta.id from filename pattern ${meta.id}_*Converted.tsv
and prefixes contigs as ${meta.id}_assembly:${contig}
"""

import sys
import re
import argparse
from pathlib import Path


def extract_meta_id(filename):
    match = re.match(r'(.+)_.*Converted\.tsv$', Path(filename).name)
    if match:
        return match.group(1)
    raise ValueError(f"Filename {filename} doesn't match pattern *_*Converted.tsv")


def process_files(input_files, output_file):
    out = open(output_file, 'w')
    
    try:
        # Write header once
        out.write("contigs\tpredictions\n")
        
        for filepath in input_files:
            meta_id = extract_meta_id(filepath)
            
            with open(filepath, 'r') as f:
                # Skip header line
                next(f)
                
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        contig, predictions = parts[0], parts[1]
                        modified_contig = f"{meta_id}_assembly:{contig}"
                        out.write(f"{modified_contig}\t{predictions}\n")
    
    finally:
        if output_file:
            out.close()


def main():
    parser = argparse.ArgumentParser(
        description='Concatenate TSV files with modified contig names',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
        Example usage:
        %(prog)s -i sample1_k2Converted.tsv sample2_k2Converted.tsv -o output.tsv
        %(prog)s -i *_*Converted.tsv -o combined.tsv
        '''
    )
    
    parser.add_argument(
        '-i', '--input',
        nargs='+',
        required=True,
        metavar='FILE',
        help='Input TSV files (pattern: *_*Converted.tsv)'
    )
    
    parser.add_argument(
        '-o', '--output',
        required=True,
        metavar='FILE',
        help='Output file (pattern: binGroup_contigTax.tsv)'
    )
    
    args = parser.parse_args()
    
    try:
        process_files(args.input, args.output)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()