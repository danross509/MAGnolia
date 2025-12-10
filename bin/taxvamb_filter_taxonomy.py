#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path


def get_fasta_contigs(fasta_path):
    contig_names = set()
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                current_contig = line[1:].split()[0]
                contig_names.add(current_contig)
            else:
                continue
    
    return contig_names

def update_taxonomy(taxonomy_path, contig_names, output_path):
    classified_contigs = set()
    added_contigs = set()
    
    with open(taxonomy_path, 'r') as infile, open(output_path, 'w') as outfile:
        # Process header
        header = infile.readline()
        outfile.write(header)
        
        # First process each classified contig
        for line in infile:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split('\t')

            # Add classified contig names to classified_contigs    
            contig_name = parts[0]
            classified_contigs.add(contig_name)
            
            # Write each line to output file
            if len(parts) < 2:
                outfile.write(line + '\t' + '\n')
            else:
                outfile.write(line + '\n')

        # Then parse all contigs, add those not classified to ouput file with blank classification
        for contig in contig_names:
            if contig not in classified_contigs:
                added_contigs.add(contig)

                # Write unclassified contigs to output file
                outfile.write(contig + '\t' + '\n')

    return added_contigs



def main():
    parser = argparse.ArgumentParser(
        description='Filter taxonomy and FASTA files based on contig lengths and taxonomy predictions'
    )
    parser.add_argument('-f', '--fasta', required=True,
                        help='Input assembly file')
    parser.add_argument('-t', '--taxonomy', required=True,
                        help='Input taxonomy file')
    parser.add_argument('-o', '--output', required=True,
                        help='Output filtered taxonomy file')
    
    args = parser.parse_args()
    
    # Validate input files exist
    if not Path(args.fasta).exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.taxonomy).exists():
        print(f"Error: Taxonomy file not found: {args.taxonomy}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Reading contigs from {args.fasta}...", file=sys.stderr)
    contig_names = get_fasta_contigs(args.fasta)
    
    print(f"Adding unclassified contig names to taxonomy file...", file=sys.stderr)
    added_contigs = update_taxonomy(args.taxonomy, contig_names, args.output)

    print(f"{len(added_contigs)} contigs from {args.fasta} without a taxonomy prediction were added to {args.taxonomy}", file=sys.stderr)

if __name__ == '__main__':
    main()