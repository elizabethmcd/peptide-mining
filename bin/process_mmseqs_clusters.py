import argparse
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

def load_data(cluster_file, fasta_file, metadata_file):
    cluster_df = pd.read_csv(cluster_file, sep='\t', header=None, names=['cluster_name', 'protein_id'])
    metadata_df = pd.read_csv(metadata_file, sep='\t')
    sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    return cluster_df, metadata_df, sequences

def process_clusters(cluster_df, sequences):
    cluster_summary = defaultdict(lambda: {'count': 0, 'representative_sequence': None})
    for _, row in cluster_df.iterrows():
        cluster_name, protein_id = row['cluster_name'], row['protein_id']
        if cluster_summary[cluster_name]['count'] == 0:
            cluster_summary[cluster_name]['representative_sequence'] = sequences.get(protein_id, None)
        cluster_summary[cluster_name]['count'] += 1
    return cluster_summary

def create_summary_df(cluster_summary):
    return pd.DataFrame({
        'cluster_name': list(cluster_summary.keys()),
        'representative_sequence': [v['representative_sequence'] for v in cluster_summary.values()],
        'protein_count': [v['count'] for v in cluster_summary.values()]
    })

def merge_metadata(cluster_df, metadata_df):
    cluster_df['mag_id'] = cluster_df['protein_id'].str.split('_id_', expand=True)[0]
    return pd.merge(cluster_df, metadata_df, left_on='mag_id', right_on='mag_id', how='left')

def main(cluster_file, fasta_file, metadata_file, output_summary, output_merged):
    # Load data
    cluster_df, metadata_df, sequences = load_data(cluster_file, fasta_file, metadata_file)
    
    # Process clusters and create summary
    cluster_summary = process_clusters(cluster_df, sequences)
    summary_df = create_summary_df(cluster_summary)
    
    # Merge with metadata
    merged_df = merge_metadata(cluster_df, metadata_df)
    
    # Save outputs
    summary_df.to_csv(output_summary, sep='\t', index=False)
    merged_df.to_csv(output_merged, sep='\t', index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process clusters and metadata")
    parser.add_argument('cluster_file', help='Path to the clusters TSV file')
    parser.add_argument('fasta_file', help='Path to the FASTA file containing sequences')
    parser.add_argument('metadata_file', help='Path to the metadata TSV file')
    parser.add_argument('output_summary', help='Path to the output summary TSV file')
    parser.add_argument('output_merged', help='Path to the output merged TSV file')

    args = parser.parse_args()

    main(args.cluster_file, args.fasta_file, args.metadata_file, 
         args.output_summary, args.output_merged)