#! /usr/bin/env nextflow

// Description
// Mine peptides from bacterial genomes or proteomic experiments for bioactivity.

log.info """\

MINE PEPTIDES FROM BACTERIAL GENOMES AND PREDICT THEIR BIOACTIVITY.
NOTE THAT THIS WORKFLOW DOES NOT HANDLE AUTOMATIC DOWNLOADING OF DATABASES.
YOU MUST INPUT PRE-PREPPED DATABASES FOR PEPTIDE SEQUENCES AND ML MODELS.
=================================================================
input_genomes                   : $params.input_genomes
genome_metadata                 : $params.genome_metadata
peptides_fasta                  : $params.peptides_fasta
models_dir                      : $params.models_dir
models_list                     : $params.models_list
outdir                          : $params.outdir
threads                         : $params.threads
"""

// define channels
// genome_fastas tuple with genome name and fasta filepath
genome_fastas = Channel.fromPath("${params.input_genomes}/*.fa")
    .map { file -> 
        def baseName = file.getBaseName()
        return [baseName, file]
    }
// metadata TSV
genome_metadata = channel.fromPath(params.genome_metadata)
// directory with peptide models
peptide_models_dir = channel.fromPath(params.models_dir)
// TXT file containing list of which models to run
peptide_models_list = channel.fromPath(params.models_list)
    .splitText()
    .map { it.trim() }
// database of peptides to compare against
peptides_db_ch = channel.fromPath(params.peptides_fasta)
// sequence identities for clustering
seq_identities = [50, 60, 70, 80, 90, 100]

// workflow steps
workflow {
    // get small ORF predictions with smorfinder
    smorfinder(genome_fastas)
    smorf_proteins = smorfinder.out.faa_file

    // combine smorf proteins into a single FASTA and food-specific split FASTAs
    combine_smorf_proteins(smorf_proteins.collect(), genome_metadata)
    combined_smorf_proteins = combine_smorf_proteins.out.combined_smorf_proteins
    food_split_smorf_proteins = combine_smorf_proteins.out.split_smorf_proteins
        .flatten() // convert list of files to channel of individual files

    // cluster food-specific proteins at 50% identity to reduce redundancy for oversampled foods and complexity of overall dataset
    cluster_combinations_ch = food_split_smorf_proteins
        .combine(Channel.from(seq_identities))
    mmseqs_cluster(cluster_combinations_ch)
    cluster_files_by_seq_id = mmseqs_cluster.out.clusters_tsv
        .groupTuple()
    rep_seqs_by_seq_id = mmseqs_cluster.out.rep_seqs
        .groupTuple()
    cluster_summary_input = cluster_files_by_seq_id
        .join(rep_seqs_by_seq_id)

    // summarize clusters
    summarize_clusters(cluster_summary_input)
    
    // all-v-all sequence identity comparisons of clustered proteins
    rep_seqs_50id = mmseqs_cluster.out.rep_seqs
        .filter { seq_id, files -> seq_id == 50 }
        .map { seq_id, files -> files }
        .collect()
    mmseqs_all_v_all(rep_seqs_50id)

    // deepsig predictions on combined, non-redundant smorf proteins
    deepsig(combined_smorf_proteins)
    deepsig_results = deepsig.out.deepsig_tsv

    // peptides.py sequence characterization on combined, non-redundant smorf proteins
    characterize_peptides(combined_smorf_proteins)
    peptides_results = characterize_peptides.out.peptides_tsv

    // DIAMOND seq similarity to peptide database with known bioactivities
    make_diamond_db(peptides_db_ch)
    peptides_dmnd_db = make_diamond_db.out.peptides_diamond_db
    diamond_blastp(combined_smorf_proteins, peptides_dmnd_db)
    blastp_results = diamond_blastp.out.blastp_hits_tsv

    // autopeptideml predictions
    model_combos_ch = combined_smorf_proteins
        .combine(peptide_models_dir)
        .combine(peptide_models_list)
    autopeptideml_predictions(model_combos_ch)

    // merge peptide stats from peptides.py, deepsig, blastp results, and autopeptideml results along with the metadata into one TSV output
    autopeptideml_results = autopeptideml_predictions.out.autopeptideml_tsv.collect()
    merge_peptide_stats(peptides_results, deepsig_results, blastp_results, genome_metadata, autopeptideml_results)

}

process smorfinder {
    tag "${genome_name}_smorfinder"
    publishDir "${params.outdir}/smorfinder", mode: 'copy'

    errorStrategy 'ignore'
    // rarely some genomes will fail for no discernible reason, skip over these

    memory = '10 GB'
    cpus = 1

    container "public.ecr.aws/v7p5x0i6/elizabethmcd/smorfinder:v0.2"
    conda "envs/smorfinder.yml"

    input:
    tuple val(genome_name), path(fasta)

    output:
    path("*.gff"), emit: gff_file
    path("*.faa"), emit: faa_file
    path("*.ffn"), emit: ffn_file
    path("*.tsv"), emit: tsv_file

    script:
    """
    smorf single ${fasta} -o ${genome_name}
    ln -s ${genome_name}/${genome_name}.gff
    ln -s ${genome_name}/${genome_name}.faa
    ln -s ${genome_name}/${genome_name}.ffn
    ln -s ${genome_name}/${genome_name}.tsv
    """

}

process combine_smorf_proteins {
    tag "combine_smorf_proteins"
    publishDir "${params.outdir}/combined_smorf_proteins", mode: 'copy'

    memory = '10 GB'
    cpus = 1
    
    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"
    conda "envs/biopython.yml"

    input:
    path(smorf_proteins)
    path(metadata_tsv)

    output: 
    path("combined_smorf_proteins.fasta"), emit: combined_smorf_proteins
    path("*.split.fasta"), emit: split_smorf_proteins

    script:
    """
    python ${baseDir}/bin/combine_fastas.py ${smorf_proteins.join(' ')} \\
    --output_file combined_smorf_proteins.fasta \\
    --metadata ${metadata_tsv} \\
    --split_dir ./
    """
    
}

process mmseqs_cluster {
    tag "mmseqs_${seq_id}id_cluster"
    publishDir "${params.outdir}/mmseqs_cluster_output/${seq_id}id_cluster_results", mode: 'copy'

    memory = '10 GB'
    cpus = 8
    
    container "public.ecr.aws/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_2"
    conda "envs/mmseqs2.yml"

    input:
    tuple path(protein_fasta_file), val(seq_id)
    
    output:
    tuple val(seq_id), path("*_rep_seq.fasta"), emit: rep_seqs
    tuple val(seq_id), path("*.tsv"), emit: clusters_tsv

    script:
    def substrate = protein_fasta_file.simpleName
    """
    mmseqs easy-cluster ${protein_fasta_file} ${substrate} tmp --min-seq-id ${seq_id/100} -c 0.8 --threads ${task.cpus}
    """   
}

process summarize_clusters {
    tag "summarize_clusters"
    publishDir "${params.outdir}/main_results/cluster_summaries", mode: 'copy'
    
    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"
    conda "envs/biopython.yml"
    
    input:
    tuple val(seq_id), path(cluster_files), path(rep_seq_files)
    
    output:
    path("*.tsv"), emit: clusters_summary
    
    script:
    """
    python ${baseDir}/bin/process_mmseqs_clusters.py \\
        --cluster_files *_cluster.tsv \\
        --fasta_files *_rep_seq.fasta \\
        --output clusters_${seq_id}id_summary.tsv
    """
}

process mmseqs_all_v_all {
    tag "mmseqs_all_v_all"
    publishDir "${params.outdir}/mmseqs_all_v_all", mode: 'copy'

    memory = '20 GB'
    cpus = 12

    container "public.ecr.aws/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_2"
    conda "envs/mmseqs2.yml"

    input:
    path("*_rep_seq.fasta")

    output:
    path("combined_representative_sequences.fasta"), emit: combined_representative_sequences
    path("*.tsv"), emit: mmseqs_easy_search_tsv

    script:
    """
    cat *_rep_seq.fasta > combined_representative_sequences.fasta
    mmseqs easy-search combined_representative_sequences.fasta combined_representative_sequences.fasta mmseqs-all-v-all-results.tsv tmp --threads ${task.cpus} --exhaustive-search

    """
}

process deepsig {
    tag "deepsig_predictions"
    publishDir "${params.outdir}/deepsig", mode: 'copy'
    
    accelerator 1, type: 'nvidia-t4'
    cpus = 8
    
    container "public.ecr.aws/biocontainers/deepsig:1.2.5--pyhca03a8a_1"
    conda "envs/deepsig.yml"

    input: 
    path(faa_file)

    output: 
    path("*.tsv"), emit: deepsig_tsv

    script: 
    """
    deepsig -f ${faa_file} -o nonredundant_smorf_proteins_deepsig.tsv -k gramp -t ${task.cpus}
    """
}

process characterize_peptides {
    tag "characterize_peptides"
    publishDir "${params.outdir}/peptide_characterization", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "elizabethmcd/peptides"
    conda "envs/peptides.yml"

    input:
    path(faa_file)

    output: 
    path("*.tsv"), emit: peptides_tsv

    script:
    """
    python ${baseDir}/bin/characterize_peptides.py ${faa_file} nonredundant_smorf_proteins_peptide_characteristics.tsv
    """
}

process make_diamond_db {
    tag "make_diamond_db"

    memory = "5 GB"
    cpus = 1

    container "public.ecr.aws/biocontainers/diamond:2.1.7--h43eeafb_1"
    conda "envs/diamond.yml"

    input:
    path(peptides_fasta)

    output:
    path("*.dmnd"), emit: peptides_diamond_db

    // ignore warnings about DNA only because of short peptides
    script:
    """
    diamond makedb --in ${peptides_fasta} -d peptides_db.dmnd --ignore-warnings 
    """
}

process diamond_blastp {
    tag "diamond_blastp"
    publishDir "${params.outdir}/diamond_blastp", mode: 'copy'

    memory = "10 GB"

    container "public.ecr.aws/biocontainers/diamond:2.1.7--h43eeafb_1"
    conda "envs/diamond.yml"

    input:
    path(faa_file)
    path(peptides_diamond_db)

    output:
    path("*.tsv"), emit: blastp_hits_tsv

    script:
    """
    diamond blastp -d ${peptides_diamond_db} \\
     -q ${faa_file} \\
     -o nonredundant_smorf_proteins_blast_results.tsv \\
     --header simple \\
     --max-target-seqs 1 \\
    --outfmt 6 qseqid sseqid full_sseq pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore
    """
}

process autopeptideml_predictions {
    tag "${model_name}_autopeptideml"
    publishDir "${params.outdir}/autopeptideml", mode: 'copy'

    memory = "10 GB"
    cpus = 6

    container "elizabethmcd/autopeptideml:latest"
    conda "envs/autopeptideml"

    input:
    tuple path(peptides_fasta), path(model_dir), val(model_name)

    output:
    path("*.tsv"), emit: autopeptideml_tsv

    script:
    """
    python3 ${baseDir}/bin/run_autopeptideml.py \\
        --input_fasta ${peptides_fasta} \\
        --model_folder "${model_dir}/${model_name}/ensemble" \\
        --model_name ${model_name} \\
        --output_tsv "autopeptideml_${model_name}.tsv"
    """
}

process merge_peptide_stats {
    tag "merge_peptide_stats"
    publishDir "${params.outdir}/main_results/peptide_results", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "public.ecr.aws/csgenetics/tidyverse:latest"
    conda "envs/tidyverse.yml"

    input:
    path(peptides_info_tsv)
    path(deepsig_tsv)
    path(blastp_hits_tsv)
    path(genome_metadata)
    path("autopeptideml_*.tsv")

    output:
    path("all_peptide_results.tsv"), emit: all_peptide_results

    script:
    """
    Rscript ${baseDir}/bin/merge_peptide_stats.R ${peptides_info_tsv} ${deepsig_tsv} ${blastp_hits_tsv} ${genome_metadata} ./ all_peptide_results.tsv
    """
}