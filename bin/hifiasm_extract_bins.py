#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_fasta(fasta_file):
    sequences = {}
    current_header = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header:
                    sequences[current_header] = ''.join(current_seq)
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_header:
            sequences[current_header] = ''.join(current_seq)
    
    return sequences


def parse_bins(bins_file):
    bins = {}
    
    with open(bins_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split('\t')
            if len(parts) < 2:
                continue
            
            bin_name = parts[0]
            contigs = parts[1:]
            bins[bin_name] = contigs
    
    return bins


def write_bin_fasta(bin_name, contigs, sequences, bin_prefix, output_dir):
    output_file = Path(output_dir) / f"{bin_prefix}.{bin_name}.fa"
    found = 0
    missing = 0
    
    with open(output_file, 'w') as f:
        for contig in contigs:
            if contig in sequences:
                f.write(f">{contig}\n")
                # Write sequence in 80 character lines
                seq = sequences[contig]
                for i in range(0, len(seq), 80):
                    f.write(seq[i:i+80] + '\n')
                found += 1
            else:
                print(f"Warning: Contig '{contig}' not found in FASTA file (bin: {bin_name})")
                missing += 1
    
    return found, missing


def main():
    parser = argparse.ArgumentParser(
        description='Extract metagenomic bins from FASTA file based on TSV bin assignments'
    )
    parser.add_argument('-c', '--contigs',  help='Input FASTA file')
    parser.add_argument('-b', '--bins',     help='Input TSV file with bin assignments')
    parser.add_argument('-p', '--prefix',   help='Bin prefix')
    parser.add_argument('-o', '--output-dir', default='bins',
                        help='Output directory for bin FASTA files (default: bins)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Reading contig file: {args.contigs}")
    sequences = parse_fasta(args.contigs)
    print(f"Loaded {len(sequences)} sequences")
    
    print(f"\nReading bins file: {args.bins}")
    bins = parse_bins(args.bins)
    print(f"Loaded {len(bins)} bins")
    
    print(f"\nExtracting bins to: {output_dir}")
    print("-" * 60)
    
    total_found = 0
    total_missing = 0
    
    for bin_name, contigs in bins.items():
        found, missing = write_bin_fasta(bin_name, contigs, sequences, args.prefix, output_dir)
        total_found += found
        total_missing += missing
        print(f"{bin_name}: {found} contigs written, {missing} missing")
    
    print("-" * 60)
    print(f"\nSummary:")
    print(f"  Total bins created: {len(bins)}")
    print(f"  Total contigs extracted: {total_found}")
    print(f"  Total contigs missing: {total_missing}")
    print(f"  Output directory: {output_dir}")


if __name__ == '__main__':
    main()