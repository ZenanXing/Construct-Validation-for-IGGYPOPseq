
# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Assign arguments to variables
input_dir <- args[1]
output_dir <- args[2]
aa_validation <- as.logical(args[3])

# Load the required packages ----------------------------------------------

library(tidyverse)
library(openxlsx)
library(Biostrings)
library(seqinr)
library(stringr)
library(purrr)

# Export the summary table ------------------------------------------------

setwd(paste0(output_dir, "/Analysis_Results/OutputFiles"))

## Demultiplexing
df_dmtplx <- read.table("Demultiplexing_ReadsSummary.tsv", header = TRUE, sep = "\t") %>% 
  mutate(SampleID = gsub("(.*)_filtered.fastq", "\\1", FileName))
colnames(df_dmtplx) <- c("FileName", "Reads_sorted_after_Demultiplexing", "SampleID")

## Alignment
df_align <- read.table("Alignment.tsv", header = TRUE, sep = "\t") %>% 
  mutate(SampleID = gsub("./Alignment/(.*).sorted.bam", "\\1", FileName),
         ref_DNA_Length = End - Start + 1)

## Assembly
df_asmbly <- read.table("ntAlignment.tsv", header = TRUE, sep = "\t") %>% 
  mutate(DNA_Identity = gsub("(.*) \\(.*%\\)", "\\1", Identity_nt),
         DNA_Identity_pct = as.numeric(gsub(".* \\((.*)%\\)", "\\1", Identity_nt)),
         DNA_matched = as.numeric(gsub("([0-9]+)\\/([0-9]+)", "\\1", DNA_Identity)),
         AssembledDNA_length = as.numeric(gsub("([0-9]+)\\/([0-9]+)", "\\2", DNA_Identity))) %>% 
  mutate(ntCorrect = ifelse(DNA_matched==AssembledDNA_length, "Y", NA)) %>% 
  dplyr::select(SampleID, PolishingStep, DNA_Identity, DNA_Identity_pct, AssembledDNA_length, ntCorrect) %>% 
  pivot_wider(names_from = PolishingStep, values_from = 3:6) %>% 
  mutate(CompleteAssembled = ifelse(DNA_Identity_pct_medaka_1rd >= 95, "Y", NA)) %>% 
  arrange(SampleID)


## Variant Calling
# Function to determine the variant
snp_fun <- function(df) {
  df_snp <- df %>% filter(Type %in% c("SNP", "SNB", "MNP"))
  if (nrow(df_snp) != 0) {
    ret <- NULL
    for (i in 1:nrow(df_snp)) {
      if (df_snp$Type[i] == "MNP") {
        if (!str_detect(df_snp$Alt_allele[i], "N|\\*")) {ret <- c(ret, "SNP")} else {ret <- c(ret, NA)}
      } else {ret <- c(ret, "SNP") }
    }
    if ("SNP" %in% ret) {ret <- "SNP"} else {ret <- NA}
  } else {
    ret <- NA
  }
  return(ret)
}
indel_fun <- function(df) {
  df_indel <- df %>% filter(!Type %in% c("SNP", "SNB"))
  if (nrow(df_indel) != 0) {
    ret <- NULL
    for (i in 1:nrow(df_indel)) {
      if (df_indel$Type[i] == "MNP") {
        if (str_detect(df_indel$Alt_allele[i], "N|\\*")) {ret <- c(ret, "INDEL")} else {ret <- c(ret, NA)}
      } else {ret <- c(ret, "INDEL") }
    }
    if ("INDEL" %in% ret) {ret <- "INDEL"} else {ret <- NA}
  } else {
    ret <- NA
  }
  return(ret)
}
# Results
# racon polishing
# racon_1rd
df_vc_11 <- read.table("VariantCalling_racon_1rd.tsv", header = TRUE, sep = "\t", colClasses = "character") %>% 
  group_by(SampleID) %>% nest() %>% 
  mutate(SNP = map_chr(data, snp_fun),
         INDEL = map_chr(data, indel_fun)) %>% 
  unite(SNP, INDEL, col = "Variant_racon_1rd", sep = " & ", na.rm = TRUE)
# racon_2rd
df_vc_12 <- read.table("VariantCalling_racon_2rd.tsv", header = TRUE, sep = "\t", colClasses = "character") %>% 
  group_by(SampleID) %>% nest() %>% 
  mutate(SNP = map_chr(data, snp_fun),
         INDEL = map_chr(data, indel_fun)) %>% 
  unite(SNP, INDEL, col = "Variant_racon_2rd", sep = " & ", na.rm = TRUE)
# racon_3rd
df_vc_13 <- read.table("VariantCalling_racon_3rd.tsv", header = TRUE, sep = "\t", colClasses = "character") %>% 
  group_by(SampleID) %>% nest() %>% 
  mutate(SNP = map_chr(data, snp_fun),
         INDEL = map_chr(data, indel_fun)) %>% 
  unite(SNP, INDEL, col = "Variant_racon_3rd", sep = " & ", na.rm = TRUE)
# medaka polishing
# medaka_1rd
df_vc_21 <- read.table("VariantCalling_medaka_1rd.tsv", header = TRUE, sep = "\t", colClasses = "character") %>% 
  group_by(SampleID) %>% nest() %>% 
  mutate(SNP = map_chr(data, snp_fun),
         INDEL = map_chr(data, indel_fun)) %>% 
  tidyr::unite(SNP, INDEL, col = "Variant_medaka_1rd", sep = " & ", na.rm = TRUE)

if (isTRUE(aa_validation)) {
  df_vc <- df_vc_11[, c(1, 3)] %>% 
    full_join(df_vc_12[, c(1, 3)], by = "SampleID") %>% 
    full_join(df_vc_13[, c(1, 3)], by = "SampleID") %>% 
    full_join(df_vc_21[, c(1, 3)], by = "SampleID") %>% 
    pivot_longer(2:5, names_to = "PolishingStep", names_pattern = "Variant_(.*)", values_to = "Variant") %>% 
    tidyr::unite(SampleID, PolishingStep, col = "comb", sep = "--")
  ## Mutation Identification
  df_ref_orf_length <- read.table("Ref_ORF_length.tsv", header = TRUE, sep = "\t") %>% 
    mutate(ReferenceName = gsub("./VariantIdentification/ReferenceORFs/(.*)_orf.fasta", "\\1", FileName)) %>% 
    left_join(read.table(paste0(input_dir, "/SampleInfo.tsv"), header = TRUE, colClasses = "character") %>% dplyr::select(SampleID, ReferenceName), by = "ReferenceName") %>% 
    dplyr::select(SampleID, RefORFlength)
  df_orf_length <- read.table("ORF_length.tsv", header = TRUE, sep = "\t") %>% 
    mutate(SampleID = gsub("(.*_[0-9]+)_(.*)", "\\1", chr),
           PolishingStep = gsub("(.*_[0-9]+)_(.*)", "\\2", chr)) %>% 
    dplyr::select(SampleID, PolishingStep, ORFlength) %>% 
    left_join(df_ref_orf_length, by = "SampleID") %>% 
    tidyr::unite(SampleID, PolishingStep, col = "comb", sep = "--")
  
  df_mutid <- read.table("aaAlignment.tsv", header = TRUE, sep = "\t") %>% 
    tidyr::unite(SampleID, PolishingStep, col = "comb", sep = "--") %>% 
    left_join(df_orf_length, by = "comb") %>% 
    mutate(Pro_Identity = gsub("(.*) \\(.*%\\)", "\\1", Identity_aa),
           Pro_Identity_pct = as.numeric(gsub(".* \\((.*)%\\)", "\\1", Identity_aa)),
           Pro_matched = as.numeric(gsub("([0-9]+)\\/([0-9]+)", "\\1", Pro_Identity)),
           Pro_length = as.numeric(gsub("([0-9]+)\\/([0-9]+)", "\\2", Pro_Identity))) %>% 
    mutate(aaCorrect = ifelse(Pro_matched == Pro_length, "Y", NA),
           aaSameLength = ifelse(RefORFlength == ORFlength, "Y", NA)) %>% 
    dplyr::select(comb, Pro_Identity, Pro_Identity_pct, ORFlength, aaCorrect, aaSameLength) %>% 
    full_join(df_vc, by = "comb") %>% 
    separate(comb, into = c("SampleID", "PolishingStep"), sep = "--") %>% 
    mutate(missenseMutation = ifelse(is.na(aaCorrect)&aaSameLength=="Y"&Variant=="SNP", "Y", NA)) %>% 
    pivot_wider(names_from = PolishingStep, values_from = 3:9) %>% 
    arrange(SampleID) %>% 
    full_join(df_ref_orf_length, by = "SampleID")
  df_mutid_misss <- df_mutid %>% 
    filter(missenseMutation_medaka_1rd == "Y") %>% 
    left_join(read.table(paste0(input_dir, "/SampleInfo.tsv"), header = TRUE) %>% dplyr::select(SampleID, ReferenceName), by = "SampleID")
  if (nrow(df_mutid_misss) != 0) {
    mat_ret <- NULL
    for (i in 1:nrow(df_mutid_misss)) {
      ref <- readAAStringSet(paste0("../VariantIdentification/ReferenceORFs/", df_mutid_misss$ReferenceName[i], "_orf.fasta"))
      smp <- readAAStringSet(paste0("../VariantIdentification/AssembledORFs/selectedORF/", df_mutid_misss$SampleID[i], "/", df_mutid_misss$SampleID[i], "_medaka_1rd_orf.fasta"))
      alignment <- pairwiseAlignment(ref, smp, type="global")
      pos <- alignment@subject@mismatch@unlistData
      #rawToChar(charToRaw(as.character(alignment@pattern))[pos])
      df_mut_temp <- data.frame(Pos = pos) %>% 
        mutate(ref_seq = as.character(alignment@pattern),
               smp_seq = as.character(alignment@subject),
               ref = substr(ref_seq, Pos, Pos),
               smp = substr(smp_seq, Pos, Pos),
               mut = paste0(ref, Pos, smp))
      ret <- data.frame(
        SampleID = df_mutid_misss$SampleID[i],
        aaMutation_medaka_1rd = paste(df_mut_temp$mut, collapse = ",")
      )
      mat_ret <- rbind(mat_ret, ret)
    }
  } else {
    mat_ret <- data.frame(
      SampleID = NA,
      aaMutation_medaka_1rd = NA
    )
  }
  df_mutid <- df_mutid %>% left_join(mat_ret, by = "SampleID")
  
} else {
  df_vc <- df_vc_11[, c(1, 3)] %>% 
    full_join(df_vc_12[, c(1, 3)], by = "SampleID") %>% 
    full_join(df_vc_13[, c(1, 3)], by = "SampleID") %>% 
    full_join(df_vc_21[, c(1, 3)], by = "SampleID")
}


## Hyperlinks to the pairwise alignment
df_hyper <- read.table(paste0(input_dir, "/SampleInfo.tsv"), header = TRUE) %>% 
  dplyr::select(SampleID) %>% 
  mutate(Set = gsub("(.*)_[0-9]+", "\\1", SampleID)) %>% 
  mutate(Alignment_nt = paste0("./ntAlignment/", SampleID, "_nt_align_muscle_f.html"),
         Alignment_nt_all_steps = paste0("./ntAlignment/", SampleID, "_nt_align_muscle.html"),
         Alignment_nt_all_clones = paste0("./ntAlignment/", Set, "_nt_align_muscle_clones.html"))
if (isTRUE(aa_validation)) {
  df_hyper <- df_hyper %>% 
    mutate(Alignment_aa = paste0("./aaAlignment/", SampleID, "_aa_align_muscle_f.html"),
           Alignment_aa_all_steps = paste0("./aaAlignment/", SampleID, "_aa_align_muscle.html"),
           Alignment_aa_all_clones = paste0("./aaAlignment/", Set, "_aa_align_muscle_clones.html"))
}

## Consensus sequence
df_csss <- data.frame(
  FileName_medaka = list.files(path = "../Assembly/medaka/1rd", pattern = ".*_medaka_1rd_f.fasta$"),
  FileName_sam = list.files(path = "../Assembly/bedtools/sam/", pattern = ".*_sam_f.fasta$"),
  #SeqID = NA,
  ConsensusSequence = NA, 
  ConsensusSequence_sam = NA
)

for (i in 1:nrow(df_csss)) {
  fasta_seq <- try(read.fasta(file = paste0("../Assembly/medaka/1rd/", df_csss$FileName_medaka[i]), as.string = TRUE, forceDNAtolower = FALSE))
  #df_csss$SeqID[i] <- names(fasta_seq)
  if (!inherits(fasta_seq, "try-error")) {temp_seq <- as.character(fasta_seq)} else {temp_seq <- ""}
  if (length(temp_seq) != 0) {df_csss$ConsensusSequence[i] <- temp_seq }
  fasta_seq <- try(read.fasta(file = paste0("../Assembly/bedtools/sam/", df_csss$FileName_sam[i]), as.string = TRUE, forceDNAtolower = FALSE))
  if (!inherits(fasta_seq, "try-error")) {temp_seq <- as.character(fasta_seq)} else {temp_seq <- ""}
  if (length(temp_seq) != 0) {df_csss$ConsensusSequence_sam[i] <- temp_seq }
}

df_csss <- df_csss %>% 
  mutate(SampleID = gsub("(.*)_medaka_1rd_f.fasta", "\\1", FileName_medaka),
         AssembledDNALength = nchar(ConsensusSequence)) %>%
  dplyr::select(SampleID, AssembledDNALength, ConsensusSequence, ConsensusSequence_sam)

## Final Table - full info.
df_exp <- read.table(paste0(input_dir, "/SampleInfo.tsv"), header = TRUE) %>% 
  dplyr::select(primer_index, SampleID, ReferenceName, n_frags) %>% 
  left_join(df_dmtplx %>% dplyr::select(SampleID, Reads_sorted_after_Demultiplexing), by = "SampleID") %>% 
  left_join(df_align %>% dplyr::select(SampleID, Reads_Mapped, Coverage_Percentage, Mean_Depth, Mean_Base_Quality, Mean_Mapping_Quality, ref_DNA_Length), by = "SampleID") %>% 
  left_join(df_asmbly, by = "SampleID")

if (isTRUE(aa_validation)) {
  df_exp <- df_exp %>% left_join(df_mutid, by = "SampleID")
} else {
  df_exp <- df_exp %>% left_join(df_vc, by = "SampleID")
}

df_exp <- df_exp %>% 
  mutate(Attention_LowDepth = ifelse(Mean_Depth<15, "Y", NA)) %>% 
  #mutate(Attention_PossibleSNPs = ifelse(is.na(ntCorrect_sam), "Y", NA)) %>% 
  left_join(df_hyper, by = "SampleID") %>% 
  left_join(df_csss, by = "SampleID") %>% 
  left_join(read.table(paste0(input_dir, "/SampleInfo.tsv"), header = TRUE) %>% dplyr::select(SampleID, ReferenceSequence), by = "SampleID")

# CheckINDEL column
ck_indel <- function(df, var, ref, ass, seq, cmplt) {
  #r <- 10
  #df <- df_ck_indel[[2]][[r]]
  #var <- df_ck_indel$Variant_medaka_1rd[r]
  #ref <- df_ck_indel$ref_DNA_Length[r]
  #ass <- df_ck_indel$AssembledDNALength[r]
  #seq <- df_ck_indel$ConsensusSequence[r]
  
  seq <- gsub("\\*", "N", seq)
  df_nt <- NULL
  for (pat in c("A", "T", "G", "C")) {  
    df_temp <- matchPattern(paste0(rep(pat, 4), collapse = ""), DNAString(seq))@ranges %>% as.data.frame() %>% mutate(PolyNT = pat)
    df_nt <- rbind(df_nt, df_temp)
  }
  
  if ((nrow(df)==1) && (df$Position%in%c(1,ref)) && (ref<ass)) {
    return("extra_nt")
    } else {
      if ((nrow(df)==1) && (!is.null(df_nt)) && any(df_nt$PolyNT %in% c(df$Ref_allele, df$Alt_allele)) && any(df$Position %in% c(df_nt$start, df_nt$end))) {
        ref_allele <- unique(unlist(strsplit(df$Ref_allele, "")))
        alt_allele <- unique(unlist(strsplit(df$Alt_allele, "")))
        if ((length(ref_allele)==1) && (length(alt_allele)==1) && (ref_allele == alt_allele)) {return("poly_nt")} else {return("")}
      } else {
        if ((var=="INDEL"||var=="SNP & INDEL") && (ass==ref) && (!is.na(cmplt))) {return("possible_SNPs")} else {return("")}
      }
  }
}
df_ck_indel <- read.table("VariantCalling_medaka_1rd.tsv", header = TRUE, sep = "\t", colClasses = "character") %>% 
  group_by(SampleID) %>% nest() %>% 
  left_join(df_exp %>% dplyr::select(SampleID, Variant_medaka_1rd, ref_DNA_Length, AssembledDNALength, ConsensusSequence, CompleteAssembled), by = "SampleID") %>% 
  mutate(CheckINDEL = purrr::pmap_chr(list(data), ck_indel, 
                                      var = Variant_medaka_1rd, ref = ref_DNA_Length, ass = AssembledDNALength, 
                                      seq = ConsensusSequence, cmplt = CompleteAssembled))


# Attention_PossibleHetSites column
het_site <- function (ref, sam, var) {
  # x <- 1 # 195 311 33 74 8 49 1
  # ref <- df_het_site$ReferenceSequence[x]
  # sam <- df_het_site$ConsensusSequence_sam[x]
  # var <- df_het_site$Variant_medaka_1rd[x]
  # possible mutations
  if (grepl("[^ACGTN\\*]", sam)){
    
    # indel
    df_temp_indel <- data.frame(mutation = unlist(str_extract_all(sam, "[a-z]+"))) %>% 
      cbind(as.data.frame(str_locate_all(sam, "[a-z]+"))) %>% 
      mutate(type = NA)
    if (nrow(df_temp_indel) != 0) {
      for (i in 1:nrow(df_temp_indel)) {
        if (length(unique(unlist(strsplit(df_temp_indel$mutation[i], "")))) == 1 && df_temp_indel$mutation[i] != "*") {
          srdseq <- substr(sam, df_temp_indel$start[i]-4, df_temp_indel$end[i]+4)
          if(str_detect(sam, regex(paste0(rep(unique(unlist(strsplit(df_temp_indel$mutation[i], ""))), 4), collapse = ""), ignore_case = TRUE))) {
            df_temp_indel$type[i] <- "polynt"
          } else {
            if (!is.na(var) && str_detect(var, "INDEL")) {df_temp_indel$type[i] <- NA} else {df_temp_indel$type[i] <- "indel"}
          }
        } else {
          rep_length <- df_temp_indel$end[i]-df_temp_indel$start[i]+1
          srdseq <- substr(sam, df_temp_indel$start[i]-4*rep_length, df_temp_indel$end[i]+4*rep_length)
          if(str_detect(sam, regex(paste0(rep(df_temp_indel$mutation[i], 3), collapse = ""), ignore_case = TRUE))) {
            df_temp_indel$type[i] <- "rep"
          } else {
            if (!is.na(var) && str_detect(var, "INDEL")) {df_temp_indel$type[i] <- NA} else {df_temp_indel$type[i] <- "indel"}
          }
        }
      }
    }
    
    # snp
    df_temp_snp <- data.frame(mutation = unlist(str_extract_all(sam, "[RYWSKMBDHV]"))) %>% 
      cbind(as.data.frame(str_locate_all(sam, "[RYWSKMBDHV]"))) %>% 
      mutate(type = "snp")
    
    if (nrow(df_temp_snp) != 0) {
      df_temp_snp <- df_temp_snp %>% 
        mutate(REF = substr(ref, start, end)) %>% 
        mutate(mutation = paste0(REF, start, mutation)) %>% 
        dplyr::select(-REF)
    }
    
    df_temp <- rbind(df_temp_indel, df_temp_snp) %>% 
      filter(!type %in% c("polynt", "rep")) %>% 
      drop_na() %>% 
      mutate(res = paste0(type, "(", mutation, ")"))
    
    res <- paste0(unique(df_temp$res), collapse = " & ")

  } else {
    
    res <- ""
    
  }
  
  return(res)

}

df_het_site <- df_csss %>% 
  left_join(df_exp %>% dplyr::select(SampleID, Variant_medaka_1rd, ReferenceSequence), by = "SampleID") %>% 
  mutate(Attention_PossibleHetSites = pmap_chr(list(ReferenceSequence, ConsensusSequence_sam, Variant_medaka_1rd), het_site))

# AssemblyStatus column
df_exp <- df_exp %>% 
  left_join(df_ck_indel %>% dplyr::select(SampleID, CheckINDEL), by = "SampleID") %>% 
  left_join(df_het_site %>% dplyr::select(SampleID, Attention_PossibleHetSites), by = "SampleID") %>% 
  mutate(Assembly_Status = ifelse(is.na(CompleteAssembled), "unassembled", 
                                  ifelse(ntCorrect_medaka_1rd == "Y"|(Variant_medaka_1rd == "INDEL" & CheckINDEL%in%c("poly_nt", "extra_nt")), 
                                         "error_free", "snp_indel"))) %>% 
  mutate(Assembly_Status = replace_na(as.character(Assembly_Status), "snp_indel")) %>% 
  dplyr::select(-Set) %>%
  arrange(primer_index)

# Assuming df_exp is your data frame
alignment_cols <- grep("^Alignment", names(df_exp), value = TRUE)
# Change class of each selected column to "hyperlink"
for (col in alignment_cols) { class(df_exp[[col]]) <- "hyperlink" }


## Table - with essential info.

if (isTRUE(aa_validation)) {
  df_exp_all <- df_exp %>% 
    dplyr::select(primer_index, SampleID, ReferenceName, n_frags, ref_DNA_Length, 
                  Reads_sorted_after_Demultiplexing, Reads_Mapped, Mean_Depth, 
                  Attention_LowDepth, Attention_PossibleHetSites, CompleteAssembled, CheckINDEL, Assembly_Status, 
                  ntCorrect_medaka_1rd, DNA_Identity_pct_medaka_1rd, Variant_medaka_1rd, 
                  aaCorrect_medaka_1rd, Pro_Identity_pct_medaka_1rd, missenseMutation_medaka_1rd, aaMutation_medaka_1rd, 
                  Alignment_nt, Alignment_nt_all_steps, Alignment_nt_all_clones, 
                  Alignment_aa, Alignment_aa_all_steps, Alignment_aa_all_clones, 
                  ConsensusSequence, ReferenceSequence)
} else {
  df_exp_all <- df_exp %>% 
    dplyr::select(primer_index, SampleID, ReferenceName, n_frags, ref_DNA_Length, 
                  Reads_sorted_after_Demultiplexing, Reads_Mapped, Mean_Depth, 
                  Attention_LowDepth, Attention_PossibleHetSites, CompleteAssembled, CheckINDEL, Assembly_Status, 
                  ntCorrect_medaka_1rd, DNA_Identity_pct_medaka_1rd, Variant_medaka_1rd, 
                  Alignment_nt, Alignment_nt_all_steps, Alignment_nt_all_clones, 
                  ConsensusSequence, ReferenceSequence)
}


## Table - final list
# Function to find the best colony if AA_validation is false
best_colony <- function(df) {
  #df <- df_exp_f[[2]][[3]]
  # no mutation in sam and final results
  df_select <- df %>% filter(ntCorrect_medaka_1rd == "Y", Attention_PossibleHetSites == "")
  if(nrow(df_select)!=0) {
    df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
    x <- df_select$SampleID[1]
  } else {
    # no mutation in sam but false indels
    df_select <- df %>% filter(CheckINDEL %in% c("poly_nt", "extra_nt"), Attention_PossibleHetSites == "")
    if(nrow(df_select)!=0) {
      # select the one annotated as "poly_nt"
      if(any(df_select$CheckINDEL=="poly_nt")) {
        df_select <- df_select %>% filter(CheckINDEL == "poly_nt") %>% dplyr::arrange(desc(Mean_Depth))
        x <- df_select$SampleID[1]
      } else {
        # select the one annotated as "extra_nt"
        df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
        x <- df_select$SampleID[1]
      }
    } else {
      # no mutation in final and no indel in sam
      df_select <- df %>% filter(ntCorrect_medaka_1rd == "Y", !str_detect(Attention_PossibleHetSites, "indel"))
      if(nrow(df_select)!=0) {
        df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
        x <- df_select$SampleID[1]
      } else {
        # no mutation in final but indel in sam
        df_select <- df %>% filter(ntCorrect_medaka_1rd == "Y")
        if(nrow(df_select)!=0) {
          df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
          x <- df_select$SampleID[1]
        } else {
          # mutation with low depth
          df_select <- df %>% filter(df$DNA_Identity_pct_medaka_1rd >= 99.9, df$Mean_Depth <= 10) %>% dplyr::arrange(desc(DNA_Identity_pct_medaka_1rd))
          if(nrow(df_select)!=0){
            x <- df_select$SampleID[1]
          } else {
            x <- "NA"
          }
        }
      }
    }
  }
  return(x)
}

# Function to find the best colony if AA_validation is true
best_colony_aa <- function(df) {
  #df <- df_exp_f[[2]][[3]]
  # no mutation in sam and final results
  df_select <- df %>% filter(ntCorrect_medaka_1rd == "Y", Attention_PossibleHetSites == "")
  if(nrow(df_select)!=0) {
    df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
    x <- df_select$SampleID[1]
  } else {
    # no mutation in sam but false indels
    df_select <- df %>% filter(CheckINDEL %in% c("poly_nt", "extra_nt"), Attention_PossibleHetSites == "")
    if(nrow(df_select)!=0) {
      # select the one annotated as "poly_nt"
      if(any(df_select$CheckINDEL=="poly_nt")) {
        df_select <- df_select %>% filter(CheckINDEL == "poly_nt") %>% dplyr::arrange(desc(Mean_Depth))
        x <- df_select$SampleID[1]
      } else {
        # select the one annotated as "extra_nt"
        df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
        x <- df_select$SampleID[1]
      }
    } else {
      # silent mutation but no mutation in sam
      df_select <- df %>% filter(aaCorrect_medaka_1rd == "Y", Attention_PossibleHetSites == "")
      if(nrow(df_select)!=0) {
        df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
        x <- df_select$SampleID[1]
      } else {
        # no mutation in final and no indel in sam
        df_select <- df %>% filter(ntCorrect_medaka_1rd == "Y", !str_detect(Attention_PossibleHetSites, "indel"))
        if(nrow(df_select)!=0) {
          df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
          x <- df_select$SampleID[1]
        } else {
          # missense mutation but no mutation in sam
          df_select <- df %>% filter(missenseMutation_medaka_1rd == "Y", Attention_PossibleHetSites == "")
          if(nrow(df_select)!=0) {
            df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
            x <- df_select$SampleID[1]
          } else {
            # no mutation in final but indel in sam
            df_select <- df %>% filter(ntCorrect_medaka_1rd == "Y")
            if(nrow(df_select)!=0) {
              df_select <- df_select %>% dplyr::arrange(desc(Mean_Depth))
              x <- df_select$SampleID[1]
            } else {
              # mutation with low depth
              df_select <- df %>% filter(df$DNA_Identity_pct_medaka_1rd >= 99.9, df$Mean_Depth <= 10) %>% dplyr::arrange(desc(DNA_Identity_pct_medaka_1rd))
              if(nrow(df_select)!=0){
                x <- df_select$SampleID[1]
              } else {
                x <- "NA"
              }
            }
          }
        }
      }
    }
  }
  return(x)
}


if (isTRUE(aa_validation)) {
  df_exp_f <- df_exp %>% 
    group_by(primer_index) %>% 
    nest() %>% 
    mutate(SampleID = purrr::map_chr(data, best_colony_aa)) %>% 
    mutate(Status = ifelse(SampleID != "NA", "Successful", "Failed")) %>% dplyr::select(primer_index, Status, SampleID)
  df_success <- df_exp_f %>% filter(Status == "Successful") %>% left_join(df_exp_all %>% dplyr::select(-c(1, 11, 13:18)), by = "SampleID")
  df_failed <- df_exp_f %>% 
    filter(Status == "Failed") %>% 
    left_join(df_exp_all %>% dplyr::select(-c(2, 11, 13:18)) %>% group_by(primer_index) %>% dplyr::slice(1), by = "primer_index") %>% 
    ungroup() %>% 
    mutate(across(c(3, 7:14, 21), ~NA)) %>% 
    mutate(across(c(15:16, 18:19), ~"NA"))
} else {
  df_exp_f <- df_exp %>% 
    group_by(primer_index) %>% 
    nest() %>% 
    mutate(SampleID = purrr::map_chr(data, best_colony)) %>% 
    mutate(Status = ifelse(SampleID != "NA", "Successful", "Failed")) %>% dplyr::select(primer_index, Status, SampleID)
  df_success <- df_exp_f %>% filter(Status == "Successful") %>% left_join(df_exp_all %>% dplyr::select(-c(1, 11, 13:16)), by = "SampleID")
  df_failed <- df_exp_f %>% 
    filter(Status == "Failed") %>% 
    left_join(df_exp_all %>% dplyr::select(-c(2, 11, 13:16)) %>% group_by(primer_index) %>% dplyr::slice(1), by = "primer_index") %>% 
    ungroup() %>% 
    mutate(across(c(3, 7:12, 16), ~NA)) %>% 
    mutate(across(c(13:14), ~"NA"))
}
alignment_cols <- grep("^Alignment", names(df_failed), value = TRUE)
for (col in alignment_cols) { class(df_failed[[col]]) <- "hyperlink" }

if ("Failed" %in% df_exp_f$Status) {
  df_exp_f <- rbind(df_success, df_failed)
} else {
  df_exp_f <- df_success
}

colnames(df_exp_f)[3] <- "clone_to_keep"

# Export the results
openxlsx::write.xlsx(df_exp_f, "Summary.xlsx", sheetName = "Summary_by_Gene", na.string = "")
wb <- loadWorkbook("Summary.xlsx")
if ("Failed" %in% df_exp_f$Status) {
  df_nt_inc <- df_exp_f %>% filter(Status == "Failed")
  df_incorrect <- df_exp_all %>% filter(primer_index %in% df_nt_inc$primer_index)
  addWorksheet(wb, "Sets_incorrect")
  writeData(wb, "Sets_incorrect", df_incorrect, na.string = "")
}
addWorksheet(wb, "Summary_all")
writeData(wb, "Summary_all", df_exp_all, na.string = "")
addWorksheet(wb, "All_Info")
writeData(wb, "All_Info", df_exp, na.string = "")
saveWorkbook(wb, "Summary.xlsx", overwrite = TRUE)

