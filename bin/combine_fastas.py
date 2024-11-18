#!/usr/bin/env python3

from Bio import SeqIO
import os
import argparse
import pandas as pd
from pathlib import Path

def get_genome_name(filename):
    return os.path.splitext(os.path.basename(filename))[0]

def is_valid_sequence(seq):
    # Check if sequence exists and isn't empty after stripping whitespace
    return bool(str(seq).strip())

def process_fasta_file(input_file, genome_name):
    new_records = {}
    with open(input_file, "r") as fasta_file:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # Skip records with empty sequences
            if not is_valid_sequence(record.seq):
                print(f"Warning: Skipping empty sequence for {record.id} in {genome_name}")
                continue
                
            # Remove trailing '*' if present
            if record.seq.endswith('*'):
                record.seq = record.seq[:-1]
            record.id = f"{genome_name}_id_{record.id}"
            new_records[record.id] = record
    return new_records

def combine_fastas(input_files, output_file, metadata_file=None, split_dir=None):
    # Read metadata if provided
    metadata = None
    if metadata_file:
        metadata = pd.read_csv(metadata_file, sep='\t')
        # Create dictionary of mag_id to fermented_food
        food_map = dict(zip(metadata['mag_id'], metadata['fermented_food']))
        # Create directory for split files if needed
        if split_dir:
            os.makedirs(split_dir, exist_ok=True)
            # Initialize dictionaries for food-specific records
            food_records = {}

    all_records = {}
    for input_file in input_files:
        genome_name = get_genome_name(input_file)
        new_records = process_fasta_file(input_file, genome_name)
        all_records.update(new_records)

        # If metadata provided, add sequences to food-specific records
        if metadata is not None and split_dir:
            food_type = food_map.get(genome_name)
            if food_type:
                if food_type not in food_records:
                    food_records[food_type] = []
                food_records[food_type].extend(list(new_records.values()))
            else:
                print(f"Warning: No food type found for {genome_name}")

    # Final validation before writing
    valid_records = []
    for record_id, record in all_records.items():
        if is_valid_sequence(record.seq):
            valid_records.append(record)
        else:
            print(f"Warning: Removing record with empty sequence: {record_id}")

    # Write combined FASTA
    with open(output_file, "w") as outfile:
        SeqIO.write(valid_records, outfile, "fasta")
        print(f"Wrote {len(valid_records)} valid sequences to {output_file}")

    # Write food-specific FASTAs if metadata provided
    if metadata is not None and split_dir:
        for food_type, records in food_records.items():
            # Clean filename by replacing problematic characters
            clean_food_type = food_type.replace('/', '_').replace(' ', '_')
            food_file = os.path.join(split_dir, f"{clean_food_type}.split.fasta")
            with open(food_file, "w") as outfile:
                SeqIO.write(records, outfile, "fasta")
                print(f"Wrote {len(records)} sequences to {food_file}")

def parse_arguments():
    parser = argparse.ArgumentParser(description='Rename FASTA headers for input files and concatenate to single FASTA')
    parser.add_argument('input_files', nargs='+', help='Input FASTA files')
    parser.add_argument('--output_file', help='Output concatenated FASTA file')
    parser.add_argument('--metadata', help='Metadata TSV file with mag_id and fermented_food columns')
    parser.add_argument('--split_dir', help='Directory to output food-specific FASTA files')
    return parser.parse_args()

def main():
    args = parse_arguments()
    combine_fastas(args.input_files, args.output_file, 
                  metadata_file=args.metadata, 
                  split_dir=args.split_dir)

if __name__ == '__main__':
    main()