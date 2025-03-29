
# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
input_dir <- args[1]
output_dir <- args[2]
aa_validation <- as.logical(args[3])

# input_dir <- "/bigdata/cutlerlab/zxing001/GeneAssembly/ForIrina/Take2_P5P3/Input"
# output_dir <- "/bigdata/cutlerlab/zxing001/GeneAssembly/ForIrina/Take2_P5P3"
# aa_validation <- TRUE

# Load the required packages ----------------------------------------------

library(tidyverse)
library(Biostrings)

# Load the Sample Info. ---------------------------------------------------

df_input <- read.table(paste0(input_dir, "/SampleInfo.tsv"), header = TRUE)

# Demultiplexing - Input Files --------------------------------------------

## Index info.
df_index <- df_input %>% 
  dplyr::select(SampleID, Fwindex, FwPrimer, Rvindex, RvPrimer) %>% 
  arrange(SampleID)
write.table(df_index, paste0(output_dir, "/Analysis_Results/InputFiles/IndexCombination.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

## Reference info.
df_ref <- df_input %>% 
  dplyr::select(SampleID, Reference) %>% 
  mutate(ReferenceFile = paste0(Reference, ".fasta")) %>% 
  dplyr::select(SampleID, ReferenceFile)
write.table(df_ref, paste0(output_dir, "/Analysis_Results/InputFiles/references.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

## Reference fasta files
df_ref <- df_input %>% dplyr::select(Reference, ReferenceSequence) %>% distinct()
  
for (i in 1:nrow(df_ref)) {
  # Creating a DNAStringSet with example sequences
  dna_sequences <- DNAStringSet(df_ref$ReferenceSequence[i])
  # Assigning names to the sequences, which will be used as headers in the FASTA file
  names(dna_sequences) <- df_ref$Reference[i]
  # Writing the DNAStringSet to a FASTA file
  writeXStringSet(dna_sequences, filepath = paste0(output_dir, "/Analysis_Results/InputFiles/references/", df_ref$Reference[i], ".fasta"))
}

if (isTRUE(aa_validation)) {
  
  # ORF Searching - Input Files ---------------------------------------------
  
  df_orf <- df_input %>% 
    mutate(Search_length = CDS_length-3) %>% 
    dplyr::select(Reference, Search_length) %>% 
    distinct()
  
  write.table(df_orf, paste0(output_dir, "/Analysis_Results/InputFiles/references_info.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
  
}

# Config file for array job submission ------------------------------------

df_config <- df_input %>% 
  mutate(ArrayTaskID = 1:nrow(df_input)) %>% 
  dplyr::select(ArrayTaskID, SampleID, Reference)
write.table(df_config, paste0(output_dir, "/Analysis_Results/InputFiles/config_sample.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

# Set Info for pairwise alignment -----------------------------------------

df_set <- df_input %>% 
  mutate(Set = gsub("(.*)_[0-9]+", "\\1", SampleID)) %>% 
  dplyr::select(Set, Reference) %>% 
  distinct()
write.table(df_set, paste0(output_dir, "/Analysis_Results/InputFiles/primer_index_info.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
