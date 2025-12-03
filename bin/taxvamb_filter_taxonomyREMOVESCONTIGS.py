#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path


def parse_fasta_lengths(fasta_path):
    contig_lengths = {}
    current_contig = None
    current_length = 0
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous contig if exists
                if current_contig is not None:
                    contig_lengths[current_contig] = current_length
                
                # Start new contig (remove '>' and take first word)
                current_contig = line[1:].split()[0]
                current_length = 0
            else:
                # Add to current sequence length
                current_length += len(line)
        
        # Save last contig
        if current_contig is not None:
            contig_lengths[current_contig] = current_length
    
    return contig_lengths

def filter_taxonomy(taxonomy_path, contig_lengths, min_length, output_path):
    kept = 0
    filtered = 0
    kept_contigs = set()
    
    with open(taxonomy_path, 'r') as infile, open(output_path, 'w') as outfile:
        # Process header
        header = infile.readline()
        outfile.write(header)
        
        # Process each contig
        for line in infile:
            line = line.strip()
            if not line:
                continue
                
            # Split to get contig name (first column)
            parts = line.split('\t')
            if len(parts) < 2:
                continue
                
            contig_name = parts[0]
            
            # Check if contig exists and meets length requirement
            if contig_name in contig_lengths:
                if contig_lengths[contig_name] >= min_length:
                    outfile.write(line + '\n')
                    kept += 1
                    kept_contigs.add(contig_name)
                else:
                    filtered += 1
            else:
                # Contig not found in assembly
                filtered += 1
                print(f"Warning: Contig '{contig_name}' not found in assembly", 
                      file=sys.stderr)
    
    return kept, filtered, kept_contigs

def filter_fasta(fasta_path, kept_contigs, output_path):
    kept = 0
    filtered = 0
    writing = False
    
    with open(fasta_path, 'r') as infile, open(output_path, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # Extract contig name
                contig_name = line[1:].strip().split()[0]
                
                # Check if this contig should be kept
                if contig_name in kept_contigs:
                    writing = True
                    kept += 1
                    outfile.write(line)
                else:
                    writing = False
                    filtered += 1
            elif writing:
                # Write sequence lines for contigs we're keeping
                outfile.write(line)
    
    return kept, filtered

def filter_abundance(abundance_path, kept_contigs, output_path):
    kept = 0
    filtered = 0
    
    with open(abundance_path, 'r') as infile, open(output_path, 'w') as outfile:
        # Process header
        header = infile.readline()
        outfile.write(header)
        
        # Process each contig
        for line in infile:
            line = line.strip()
            if not line:
                continue
                
            # Split to get contig name (first column)
            parts = line.split('\t')
            if len(parts) < 2:
                continue
                
            contig_name = parts[0]
            
            # Check if contig exists and meets length requirement
            if contig_name in kept_contigs:
                outfile.write(line + '\n')
                kept += 1
            else:
                filtered += 1
    
    return kept, filtered


def main():
    parser = argparse.ArgumentParser(
        description='Filter taxonomy and FASTA files based on contig lengths and taxonomy predictions'
    )
    parser.add_argument('-f', '--fasta', required=True,
                        help='Input assembly file')
    parser.add_argument('-t', '--taxonomy', required=True,
                        help='Input taxonomy file')
    parser.add_argument('-a', '--abundance', required=True,
                        help='Input abundance file')
    parser.add_argument('-l', '--min-length', type=int, required=True,
                        help='Minimum contig length to keep')
    parser.add_argument('-ot', '--outputTax', required=True,
                        help='Output filtered taxonomy file')
    parser.add_argument('-oa', '--outputAbundance', required=True,
                        help='Output filtered abundance file')
    parser.add_argument('-of', '--outputFasta', required=True,
                        help='Output filtered FASTA file')
    
    args = parser.parse_args()
    
    # Validate input files exist
    if not Path(args.fasta).exists():
        print(f"Error: FASTA file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)
    
    if not Path(args.taxonomy).exists():
        print(f"Error: Taxonomy file not found: {args.taxonomy}", file=sys.stderr)
        sys.exit(1)
    
    # Parse contig lengths from FASTA
    print(f"Reading contig lengths from {args.fasta}...", file=sys.stderr)
    contig_lengths = parse_fasta_lengths(args.fasta)
    print(f"Found {len(contig_lengths)} contigs in assembly", file=sys.stderr)
    
    # Filter taxonomy file
    print(f"Filtering taxonomy with minimum length {args.min_length}...", 
          file=sys.stderr)
    kept_tax, filtered_tax, kept_contigs = filter_taxonomy(
        args.taxonomy, contig_lengths, args.min_length, args.outputTax
    )
    
    # Filter contig file
    print(f"Filtering contigs to keep only contigs with taxonomy...", 
          file=sys.stderr)
    kept_fasta, filtered_fasta = filter_fasta(
        args.fasta, kept_contigs, args.outputFasta
    )

    # Filter abundance file
    print(f"Filtering abundances to keep only contigs with taxonomy...", 
          file=sys.stderr)
    kept_abundance, filtered_abundance = filter_abundance(
        args.abundance, kept_contigs, args.outputAbundance
    )

if __name__ == '__main__':
    main()