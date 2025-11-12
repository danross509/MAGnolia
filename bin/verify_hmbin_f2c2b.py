#!/usr/bin/env python3

import sys
import re

def process_tsv(input_file, output_file):
    """
    Process TSV file, splitting rows with multiple contigs.
    
    Args:
        input_file: Path to input TSV file
        output_file: Path to output TSV file
    """
    processed_rows = []
    
    with open(input_file, 'r') as f:
        for line_num, line in enumerate(f, 1):
            line = line.rstrip('\n')
            columns = line.split('\t')
            
            # Verify we have exactly 2 columns
            if len(columns) != 2:
                print(f"Warning: Line {line_num} has {len(columns)} columns, expected 2")
                continue
            
            first_col = columns[0]
            second_col = columns[1]
            
            # Check if first column contains multiple space-separated items
            items = first_col.split()
            
            if len(items) == 1:
                # Standard row with single contig
                processed_rows.append(f"{first_col}\t{second_col}")
            else:
                # Row with multiple contigs (e.g., "rc.9 s2.ctg000252l- s2.ctg006726l+")
                # Skip the first item (rc.X) and process the rest
                for item in items[1:]:
                    # Remove strand information (+ or -)
                    contig = re.sub(r'[+-]$', '', item)
                    processed_rows.append(f"{contig}\t{second_col}")
                
                print(f"Line {line_num}: Split {len(items)-1} contigs from combined row")
    
    # Write processed rows to output file
    with open(output_file, 'w') as f:
        for row in processed_rows:
            f.write(row + '\n')
    
    print(f"\nProcessing complete!")
    print(f"Total output rows: {len(processed_rows)}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python verify_hmbin_f2c2b.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    try:
        process_tsv(input_file, output_file)
    except FileNotFoundError:
        print(f"Error: File '{input_file}' not found")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)