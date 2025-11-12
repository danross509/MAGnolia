#!/usr/bin/env python3

import sys
import os
from pathlib import Path

def split_fasta(input_file, output_dir, bin_prefix):
    input_path = Path(input_file)
    output_dir = Path(output_dir)
    
    with open(input_path, 'r') as f:
        current_header = None
        current_seq = []
        file_count = 0
        
        for line in f:
            line = line.strip()
            
            if line.startswith('>'):
                # Save previous sequence if exists
                if current_header:
                    save_sequence(current_header, current_seq, output_dir, bin_prefix)
                    file_count += 1
                
                # Start new sequence
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                # Accumulate sequence lines
                if line:  # Skip empty lines
                    current_seq.append(line)
        
        # Save last sequence
        if current_header:
            save_sequence(current_header, current_seq, output_dir, bin_prefix)
            file_count += 1

def save_sequence(header, seq_lines, output_dir, bin_prefix):
    """Save a single sequence to its own file."""
    # Create safe filename from header (take first word, remove special chars)
    safe_name = header.split()[0].replace('/', '_').replace('\\', '_')
    #output_file = output_dir / f"{safe_name}.fa"

    # Handle duplicate names
    counter = 0
    output_file = output_dir / f"{bin_prefix}.binCR{counter}.fa"
    while output_file.exists():
        output_file = output_dir / f"{bin_prefix}.binCR{counter}.fa"
        counter += 1
    
    with open(output_file, 'w') as f:
        f.write(f">{header}\n")
        f.write('\n'.join(seq_lines) + '\n')

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python hifiasm_split_circular_contigs.py <input.fa> [output_directory]")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    bin_prefix = sys.argv[3]
    
    split_fasta(input_file, output_dir, bin_prefix)