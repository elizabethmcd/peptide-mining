#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import sys

def process_all_clusters(cluster_files, fasta_files):
    """Process all cluster assignments and sequences"""
    
    print(f"Processing {len(fasta_files)} FASTA files and {len(cluster_files)} cluster files")
    
    # Read all sequences from all FASTA files
    sequences = {}
    for fasta_file in fasta_files:
        try:
            sequences.update({
                record.id: str(record.seq) 
                for record in SeqIO.parse(fasta_file, "fasta")
            })
        except Exception as e:
            print(f"Error reading FASTA file {fasta_file}: {e}", file=sys.stderr)
    
    print(f"Loaded {len(sequences)} sequences")
    
    # Process all cluster files
    all_clusters = []
    for cluster_file in cluster_files:
        try:
            # Get substrate from filename
            substrate = cluster_file.split('_cluster.tsv')[0]
            print(f"Processing clusters for {substrate}")
            
            # Read cluster assignments
            cluster_df = pd.read_csv(cluster_file, sep='\t', header=None, 
                                   names=['cluster_name', 'protein_id'])
            
            # Process clusters
            cluster_data = []
            for cluster_name, group in cluster_df.groupby('cluster_name'):
                rep_seq = sequences.get(group.iloc[0]['protein_id'])
                if rep_seq is None:
                    print(f"Warning: No sequence found for representative {group.iloc[0]['protein_id']}")
                    continue
                    
                cluster_data.append({
                    'cluster_name': cluster_name,
                    'representative_sequence': rep_seq,
                    'protein_count': len(group),
                    'substrate': substrate
                })
            
            if cluster_data:
                summary_df = pd.DataFrame(cluster_data)
                all_clusters.append(summary_df)
                print(f"Found {len(cluster_data)} clusters for {substrate}")
            
        except Exception as e:
            print(f"Error processing cluster file {cluster_file}: {e}", file=sys.stderr)
    
    if not all_clusters:
        print("No valid clusters found!", file=sys.stderr)
        return pd.DataFrame(columns=['cluster_name', 'representative_sequence', 
                                   'protein_count', 'substrate'])
    
    return pd.concat(all_clusters, ignore_index=True)

def main():
    parser = argparse.ArgumentParser(description="Summarize MMseqs2 clusters")
    parser.add_argument('--cluster_files', nargs='+', required=True,
                      help='Cluster assignment TSV files')
    parser.add_argument('--fasta_files', nargs='+', required=True,
                      help='Representative sequence FASTA files')
    parser.add_argument('--output', required=True,
                      help='Output summary TSV file')
    
    args = parser.parse_args()
    
    print("Starting cluster processing...")
    print(f"Cluster files: {args.cluster_files}")
    print(f"FASTA files: {args.fasta_files}")
    
    # Process all clusters
    summary_df = process_all_clusters(args.cluster_files, args.fasta_files)
    
    # Save to file
    summary_df.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote summary for {len(summary_df)} total clusters to {args.output}")

if __name__ == '__main__':
    main()