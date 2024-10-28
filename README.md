# Mining bacterial genomes for bioactive peptides

This repository contains the `peptide-mining` Nextflow workflow for mining bacterial genomes (either metagenome-assembled genomes or isolate genomes) for bioactive peptides. 

The main input is a directory of bacterial genomes in FASTA format. The workflow predicts small ORFs (smORFs) and compares these to known databases of peptides, reports summary statistics, and predicts bioactivity based on existing ML classification models. Note that this workflow is highly configured for custom purposes, such as inputting a specific genome metadata TSV for joining metadata with main results files. You can see an example of the required metadata TSV file in `test_data/metadata/`. 

## Workflow Usage

The pipeline can be run with either conda or docker using the `-profile` flag. All input genomes in the input directory should end in `.fa`. 

Importantly, the workflow does not handle automatic downloading of databases, so these need to be prepared beforehand and input as parameters to the workflow. 

```
nextflow run main.nf \\
--input_genomes <INPUT_DIRECTORY> \\
--outdir <OUTPUT_DIRECTORY> \\
--genome_metadata <GENOME_METADATA_TSV> \\
--peptides_fasta <PEPTIDES_FILE_FOR_COMPARISON> \\
--models_dir <MODELS_DIRECTORY> \\
--models_list <MODELS_LIST_FILE> \\ 
--threads <THREADS> \\
-profile <docker|conda>
```