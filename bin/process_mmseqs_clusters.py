#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict
import glob

def process_all_clusters(cluster_files, fasta_files):
    """Process all cluster assignments and sequences"""
    
    # Read all sequences from all FASTA files
    sequences = {}
    for fasta_file in fasta_files:
        sequences.update({
            record.id: str(record.seq) 
            for record in SeqIO.parse(fasta_file, "fasta")
        })
    
    # Process all cluster files
    all_clusters = []
    for cluster_file in cluster_files:
        # Get substrate from filename
        substrate = cluster_file.split('_cluster.tsv')[0]
        
        # Read cluster assignments
        cluster_df = pd.read_csv(cluster_file, sep='\t', header=None, 
                               names=['cluster_name', 'protein_id'])
        
        # Process clusters
        cluster_summary = defaultdict(lambda: {'count': 0, 'representative_sequence': None})
        
        for _, row in cluster_df.iterrows():
            cluster_name, protein_id = row['cluster_name'], row['protein_id']
            if cluster_summary[cluster_name]['count'] == 0:
                cluster_summary[cluster_name]['representative_sequence'] = sequences.get(protein_id, None)
            cluster_summary[cluster_name]['count'] += 1
        
        # Convert to DataFrame
        summary_df = pd.DataFrame({
            'cluster_name': list(cluster_summary.keys()),
            'representative_sequence': [v['representative_sequence'] for v in cluster_summary.values()],
            'protein_count': [v['count'] for v in cluster_summary.values()],
            'substrate': substrate
        })
        
        all_clusters.append(summary_df)
    
    # Combine all summaries
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
    
    # Process all clusters
    summary_df = process_all_clusters(args.cluster_files, args.fasta_files)
    
    # Save to file
    summary_df.to_csv(args.output, sep='\t', index=False)
    print(f"Wrote summary for {len(summary_df)} total clusters to {args.output}")

if __name__ == '__main__':
    main()