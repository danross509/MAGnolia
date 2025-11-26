#!/usr/bin/env python3

import re
import sys

def extract_taxid(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if not line.strip():
                continue
            
            parts = line.strip().split('\t')
            
            if len(parts) >= 3:
                classification = parts[0]  # C, U
                id_num = parts[1]  # Contig ID
                species_info = parts[2]  # Species
                
                match = re.search(r'taxid (\d+)', species_info)
                
                if match:
                    taxid = match.group(1)
                    outfile.write(f"{classification}\t{id_num}\t{taxid}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: taxconverter_k2_pre-conversion.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        extract_taxid(input_file, output_file)
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)