#! /usr/local/bin/R

library(tidyverse)

######################################################################################
# Merge peptide stats with metadata
# Merge results from DeepSig, peptides.py, and DIAMOND Blastp results against database
######################################################################################

# command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# input files
peptides_info_tsv <- args[1]
deep_sig_tsv <- args[2]
diamond_tsv <- args[3]
genome_metadata <- args[4]
autopeptideml_dir <- args[5]
output_tsv <- args[6]

# read in files
peptides_info <- read_tsv(peptides_info_tsv)
deepsig_info <- read_tsv(deep_sig_tsv, col_names = c("peptide_id", 
                                                    "tool",
                                                    "deepsig_feature",
                                                    "deepsig_feature_start",
                                                    "deepsig_feature_end",
                                                    "deepsig_feature_score",
                                                    "tmp1",
                                                    "tmp2",
                                                    "deepsig_description"))  %>% 
                select(-tool, -tmp1, -tmp2)

diamond_blast_results <- read_tsv(diamond_tsv)  %>% 
    mutate(peptide_id = qseqid)  %>% 
    select(-qseqid)

genome_metadata <- read_tsv(genome_metadata)  %>% 
    select(mag_id, substrate, species, group)

autopeptideml_files <- list.files(path = autopeptideml_dir, pattern = "autopeptideml_.*\\.tsv$", full.names=TRUE)

autopeptideml_df <- map_dfr(autopeptideml_files, function(file) {
    read_tsv(file) %>%
        # Rename the third column to 'value' and get bioactivity name from column name
        rename(bioactivity = 3) %>%
        # Add column for bioactivity name (original column name)
        mutate(bioactivity_name = names(read_tsv(file))[3])
}) %>%
    # Pivot wider to get one column per bioactivity
    pivot_wider(
        id_cols = c(peptide_id, sequence),
        names_from = bioactivity_name,
        values_from = bioactivity
    )

# merge peptides info with deepsig results
# merge with diamond blastp results
# create mag_id column for merging with metadata
# merge with metadata
all_peptide_info_metadata <- left_join(deepsig_info, peptides_info)  %>% 
    left_join(diamond_blast_results)  %>% 
    mutate(mag_id = str_extract(peptide_id, "^.*?(?=_id_)"))  %>% 
    left_join(genome_metadata)  %>% 
    left_join(autopeptideml_df)  %>% 
    select(mag_id, peptide_id, sequence, substrate, species, group, everything())

# write to tsv
write_tsv(all_peptide_info_metadata, output_tsv)

