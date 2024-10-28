# Mining bacterial genomes for bioactive peptides

This repository contains the `peptide-mining` Nextflow workflow for mining bacterial genomes (either metagenome-assembled genomes or isolate genomes) for bioactive peptides. 

The main input is a directory of bacterial genomes in FASTA format. The workflow predicts small ORFs (smORFs) and compares these to known databases of peptides, reports summary statistics, and predicts bioactivity based on existing ML classification models. Note that this workflow is highly configured for custom purposes, such as inputting a specific genome metadata TSV for joining metadata with main results files. You can see an example of the required metadata TSV file in `test_data/metadata/`. 

## Workflow Usage

The pipeline can be run with either conda or docker using the `-profile` flag. All input genomes in the input directory should end in `.fa`. 

Importantly, the workflow does not handle automatic downloading of databases, so these need to be prepared beforehand and input as parameters to the workflow. See below for preparing databases and examples of publicly available databases for download.

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

## Databases

The main external databases that are provided to this workflow are sequences of known peptides and machine-learning classification models for predicting bioactivity. We curated sets of peptide sequences from [Peptipedia](https://app.peptipedia.cl/) and used bioactivity classification models from [AutoPeptideML], which are available on [Zenodo](https://zenodo.org/records/13363975). You can provide any input set of peptide sequences as long as they are in FASTA format. By providing a TXT file of the list of models you want to perform bioactivity searches, you can select certain models or all models available through AutoPeptideML, or additionally your own models that you made. 