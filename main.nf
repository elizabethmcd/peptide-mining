#! /usr/bin/env nextflow

// Description
// Mine peptides from bacterial genomes or proteomic experiments for bioactivity.

log.info """\

MINE PEPTIDES FROM BACTERIAL GENOMES AND PREDICT THEIR BIOACTIVITY.
NOTE THAT THIS WORKFLOW DOES NOT HANDEL AUTOMATIC DOWNLOADING OF DATABASES.
YOU MUST PREPARE AND INPUT DATABASES FOR PEPTIDE SEQUENCES AND ML MODELS.
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

genome_metadata = channel.fromPath(params.genome_metadata)
peptide_models_dir = channel.fromPath(params.models_dir)
peptide_models_list = channel.fromPath(params.models_list)
    .splitText()
    .map { it.trim() }
peptides_db_ch = channel.fromPath(params.peptides_fasta)

// workflow steps
workflow {
    // make genome STB
    all_genome_fastas_ch = genome_fastas.map{ it[1] }.collect()
    make_genome_stb(all_genome_fastas_ch)
    genome_stb_tsv = make_genome_stb.out.stb_tsv
    
    // get small ORF predictions with smorfinder
    smorfinder(genome_fastas)
    smorf_proteins = smorfinder.out.faa_file

    // combine smorf proteins into a single FASTA
    combine_smorf_proteins(smorf_proteins.collect())
    combined_smorf_proteins = combine_smorf_proteins.out.combined_smorf_proteins

    // cluster smorf proteins 95% identity and get representative seqs
    mmseqs_95id_cluster(combined_smorf_proteins)
    nonredundant_smorfs = mmseqs_95id_cluster.out.nonredundant_seqs_fasta
    mmseqs_clusters = mmseqs_95id_cluster.out.cluster_summary_tsv

    // mmseqs cluster summaries and stats, merging with metadata
    summarize_mmseqs_clusters(mmseqs_clusters, nonredundant_smorfs, genome_metadata)

    // deepsig predictions on combined, non-redundant smorf proteins
    deepsig(nonredundant_smorfs)
    deepsig_results = deepsig.out.deepsig_tsv

    // peptides.py sequence characterization on combined, non-redundant smorf proteins
    characterize_peptides(nonredundant_smorfs)
    peptides_results = characterize_peptides.out.peptides_tsv

    // DIAMOND seq similarity to Peptipedia peptide sequences of interest
    make_diamond_db(peptides_db_ch)
    peptides_dmnd_db = make_diamond_db.out.peptides_diamond_db
    diamond_blastp(nonredundant_smorfs, peptides_dmnd_db)
    blastp_results = diamond_blastp.out.blastp_hits_tsv

    // merge peptide stats from peptides.py, deepsig, and blastp results
    merge_peptide_stats(peptides_results, deepsig_results, blastp_results, genome_metadata)

    // autopeptideml predictions
    autopeptideml_predictions(peptide_models_dir, peptide_models_list, nonredundant_smorfs)

}

process make_genome_stb {
    tag "make_genome_stb"
    publishDir "${params.outdir}/genomestb", mode: 'copy'

    memory = '10 GB'
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"
    conda "envs/biopython.yml"

    input:
    path(fasta_files)

    output:
    path("*.tsv"), emit: stb_tsv

    script:
    """
    python ${baseDir}/bin/generate_genome_stb.py ${fasta_files.join(' ')} -o genomes_stb.tsv
    """

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

    output: 
    path("*.fasta"), emit: combined_smorf_proteins

    script:
    """
    python ${baseDir}/bin/combine_fastas.py ${smorf_proteins.join(' ')} combined_smorf_proteins.fasta
    """
    
}

process mmseqs_95id_cluster {
    tag "mmseqs_95id_cluster"
    publishDir "${params.outdir}/mmseqs_95id_cluster", mode: 'copy'

    memory = '10 GB'
    cpus = 8
    
    container "public.ecr.aws/biocontainers/mmseqs2:15.6f452--pl5321h6a68c12_2"
    conda "envs/mmseqs2.yml"

    input:
    path(protein_fasta_file)
    
    output:
    path("*_rep_seq.fasta"), emit: nonredundant_seqs_fasta
    path("*.tsv"), emit: cluster_summary_tsv

    script:
    """
    mmseqs easy-cluster ${protein_fasta_file} nonredundant_smorf_proteins tmp --min-seq-id 0.95 --threads ${task.cpus}
    """   
}

process summarize_mmseqs_clusters {
    tag "summarize_mmseqs_clusters"
    publishDir "${params.outdir}/main_results/mmseqs_clusters", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "quay.io/biocontainers/mulled-v2-949aaaddebd054dc6bded102520daff6f0f93ce6:aa2a3707bfa0550fee316844baba7752eaab7802-0"
    conda "envs/biopython.yml"

    input:
    path(mmseqs_cluster_file)
    path(mmseqs_nonredundant_seqs)
    path(genome_metadata_tsv)

    output:
    path("mmseqs_summary.tsv"), emit: mmseqs_summary
    path("mmseqs_metadata.tsv"), emit: mmseqs_metadata
    path("mmseqs_substrate_counts.tsv"), emit: mmseqs_substrate_counts
    path("mmseqs_phylo_groups_counts.tsv"), emit: mmseqs_phylo_groups_counts

    script:
    """
    python ${baseDir}/bin/process_mmseqs_clusters.py ${mmseqs_cluster_file} ${mmseqs_nonredundant_seqs} ${genome_metadata_tsv} mmseqs_summary.tsv mmseqs_metadata.tsv mmseqs_substrate_counts.tsv mmseqs_phylo_groups_counts.tsv
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

process merge_peptide_stats {
    tag "merge_peptide_stats"
    publishDir "${params.outdir}/main_results/peptide_stats", mode: 'copy'

    memory = "10 GB"
    cpus = 1

    container "public.ecr.aws/csgenetics/tidyverse:latest"
    conda "envs/tidyverse.yml"

    input:
    path(peptides_info_tsv)
    path(deepsig_tsv)
    path(blastp_hits_tsv)
    path(genome_metadata)

    output:
    path("all_peptide_stats.tsv"), emit: all_peptide_stats

    script:
    """
    Rscript ${baseDir}/bin/merge_peptide_stats.R ${peptides_info_tsv} ${deepsig_tsv} ${blastp_hits_tsv} ${genome_metadata} all_peptide_stats.tsv
    """
}

process autopeptideml_predictions {
    tag "${model_name}_autopeptideml"
    publishDir "${params.outdir}/autopeptideml", mode: 'copy'

    memory = "20 GB"
    cpus = 6

    container "elizabethmcd/autopeptideml:latest"
    conda "envs/autopeptideml"

    input:
    path(model_dir)
    val(model_name)
    path(peptides_fasta)

    output:
    path("*.tsv"), emit: autopeptideml_tsv

    script:
    """
    python3 ${baseDir}/bin/run_autopeptideml.py \\
        --input_fasta ${peptides_fasta} \\
        --model_folder ${model_dir} \\
        --model_name ${model_name} \\
        --output_tsv "autopeptideml_${model_name}.tsv"
    """
}