#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2     # This will be the number of CPUs per individual array job
#SBATCH --mem=20G     # This will be the memory per individual array job
#SBATCH --time=0-04:00:00     # 4 hrs
#SBATCH --job-name="analysis"

# Set variables
SCRIPT_DIR=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
AA_Validation=$4
Medaka_Mod=$5

# Specify the path to the config file
config=$OUTPUT_DIR/Analysis_Results/InputFiles/config_sample.txt

# SampleID
SampleID=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
# ReferenceName
ReferenceName=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

### Load the required modules
module load parallel # Run the code in parallel
module load minimap2 samtools bedtools # Alignment & Assembly
module load seqtk # Polishing
module load bcftools # Variant Calling
module load emboss # Variant Identification
module load muscle/3.8.31 # Pairwise Alignment

### Alignment
## Change the directory and create required directory:
cd $OUTPUT_DIR/Analysis_Results/Alignment
## Alignment suing minimap2
minimap2 -ax map-ont --secondary=no ../InputFiles/references/${ReferenceName}.fasta ../Demultiplexing/final/${SampleID}_filtered.fastq > ${SampleID}.sam
# Convert SAM to BAM
samtools view -q 60 -F 2048 -S -b ${SampleID}.sam -o ${SampleID}.bam
# Sort BAM files
samtools sort ${SampleID}.bam -o ${SampleID}.sorted.bam
# Index Sorted BAM files
samtools index ${SampleID}.sorted.bam

### Variant calling - before Polishing - using bcftools
cd $OUTPUT_DIR/Analysis_Results/VariantCalling
## Variant calling using bcftools and samtools
bcftools mpileup -Ou -f ../InputFiles/references/${ReferenceName}.fasta ../Alignment/${SampleID}.sorted.bam | bcftools call -cv -Ov > ./beforePolishing/${SampleID}_bP.vcf

### Assembly
cd $OUTPUT_DIR/Analysis_Results/Assembly
## Draft Assembly
# assembly using samtools
samtools consensus -a -A -f fasta --show-ins no --show-del yes -o ./bedtools/sam/${SampleID}_sam.fasta ../Alignment/${SampleID}.sorted.bam
# apply the variant using bedtools
bedtools maskfasta -fi ./bedtools/sam/${SampleID}_sam.fasta -fo ./bedtools/bed/${SampleID}_draft.fasta -bed ../VariantCalling/beforePolishing/${SampleID}_bP.vcf

### Polishing by Racon
module load racon
## racon polishing - 1rd
minimap2 -ax map-ont ./bedtools/bed/${SampleID}_draft.fasta ../Demultiplexing/final/${SampleID}_filtered.fastq > ./racon/sam/${SampleID}_racon_1rd.sam
racon -m 8 -x -6 -g -8 -w 500 ../Demultiplexing/final/${SampleID}_filtered.fastq ./racon/sam/${SampleID}_racon_1rd.sam ./bedtools/bed/${SampleID}_draft.fasta > ./racon/1rd/${SampleID}_racon_1rd.fasta

## racon polishing - 2rd
minimap2 -ax map-ont ./racon/1rd/${SampleID}_racon_1rd.fasta ../Demultiplexing/final/${SampleID}_filtered.fastq > ./racon/sam/${SampleID}_racon_2rd.sam
racon -m 8 -x -6 -g -8 -w 500 ../Demultiplexing/final/${SampleID}_filtered.fastq ./racon/sam/${SampleID}_racon_2rd.sam ./racon/1rd/${SampleID}_racon_1rd.fasta > ./racon/2rd/${SampleID}_racon_2rd.fasta

## racon polishing - 3rd
minimap2 -ax map-ont ./racon/2rd/${SampleID}_racon_2rd.fasta ../Demultiplexing/final/${SampleID}_filtered.fastq > ./racon/sam/${SampleID}_racon_3rd.sam
racon -m 8 -x -6 -g -8 -w 500 ../Demultiplexing/final/${SampleID}_filtered.fastq ./racon/sam/${SampleID}_racon_3rd.sam ./racon/2rd/${SampleID}_racon_2rd.fasta > ./racon/3rd/${SampleID}_racon_3rd.fasta

### Polishing with medaka
conda init
source ~/bigdata/.conda/envs/chopper_env/bin/chopper

module load medaka
## 1rd
medaka_consensus -f -x -i ../Demultiplexing/final/${SampleID}_filtered.fastq -d ./racon/3rd/${SampleID}_racon_3rd.fasta -o ./medaka/1rd/${SampleID} -m ${Medaka_Mod}; mv ./medaka/1rd/${SampleID}/consensus.fasta ./medaka/1rd/${SampleID}_medaka_1rd.fasta

## Further processing the resulting fasta files so we can process the fasta files in the multiple alignments
## Change the header of polished fasta files and saved it in a new fasta files
# sam:
awk -v new_header="${SampleID}_sam" 'BEGIN {OFS = FS} NR == 1 {sub(/^>[^ ]*/, ">" new_header)} {print}' ./bedtools/sam/${SampleID}_sam.fasta > ./bedtools/sam/${SampleID}_sam_f.fasta
# draft:
awk -v new_header="${SampleID}_draft" 'BEGIN {OFS = FS} NR == 1 {sub(/^>[^ ]*/, ">" new_header)} {print}' ./bedtools/bed/${SampleID}_draft.fasta > ./bedtools/bed/${SampleID}_draft_f.fasta
# 1rd of racon polishing:
awk -v new_header="${SampleID}_racon_1rd" 'BEGIN {OFS = FS} NR == 1 {sub(/^>[^ ]*/, ">" new_header)} {print}' ./racon/1rd/${SampleID}_racon_1rd.fasta > ./racon/1rd/${SampleID}_racon_1rd_f.fasta
# 2rd of racon polishing:
awk -v new_header="${SampleID}_racon_2rd" 'BEGIN {OFS = FS} NR == 1 {sub(/^>[^ ]*/, ">" new_header)} {print}' ./racon/2rd/${SampleID}_racon_2rd.fasta > ./racon/2rd/${SampleID}_racon_2rd_f.fasta
# 3rd of racon polishing:
awk -v new_header="${SampleID}_racon_3rd" 'BEGIN {OFS = FS} NR == 1 {sub(/^>[^ ]*/, ">" new_header)} {print}' ./racon/3rd/${SampleID}_racon_3rd.fasta > ./racon/3rd/${SampleID}_racon_3rd_f.fasta
# 1rd of medaka polishing:
awk -v new_header="${SampleID}_medaka_1rd" 'BEGIN {OFS = FS} NR == 1 {sub(/^>[^ ]*/, ">" new_header)} {print}' ./medaka/1rd/${SampleID}_medaka_1rd.fasta > ./medaka/1rd/${SampleID}_medaka_1rd_f.fasta


### Variant Calling 
cd $OUTPUT_DIR/Analysis_Results/VariantCalling
module load medaka
## after racon polishing
# 1st round
medaka tools consensus2vcf --mode NW --out_prefix ./afterPolishing/racon/1rd/${SampleID}_racon_1rd ../Assembly/racon/1rd/${SampleID}_racon_1rd.fasta ../InputFiles/references/${ReferenceName}.fasta
medaka tools classify_variants ./afterPolishing/racon/1rd/${SampleID}_racon_1rd.vcf
# 2nd round
medaka tools consensus2vcf --mode NW --out_prefix ./afterPolishing/racon/2rd/${SampleID}_racon_2rd ../Assembly/racon/2rd/${SampleID}_racon_2rd.fasta ../InputFiles/references/${ReferenceName}.fasta
medaka tools classify_variants ./afterPolishing/racon/2rd/${SampleID}_racon_2rd.vcf
# 3rd round
medaka tools consensus2vcf --mode NW --out_prefix ./afterPolishing/racon/3rd/${SampleID}_racon_3rd ../Assembly/racon/3rd/${SampleID}_racon_3rd.fasta ../InputFiles/references/${ReferenceName}.fasta
medaka tools classify_variants ./afterPolishing/racon/3rd/${SampleID}_racon_3rd.vcf
## after medaka polishing
# 1st round
medaka tools consensus2vcf --mode NW --out_prefix ./afterPolishing/medaka_1rd/${SampleID}_medaka_1rd ../Assembly/medaka/1rd/${SampleID}_medaka_1rd.fasta ../InputFiles/references/${ReferenceName}.fasta
medaka tools classify_variants ./afterPolishing/medaka_1rd/${SampleID}_medaka_1rd.vcf


### Visualization of the pairwise alignment - nt
cd $OUTPUT_DIR/Analysis_Results/PairwiseAlignment
## Align all the seqs by muscle
# Combine all the fasta files in one file
cat ../InputFiles/references/${ReferenceName}.fasta ../Assembly/bedtools/sam/${SampleID}_sam_f.fasta ../Assembly/bedtools/bed/${SampleID}_draft_f.fasta ../Assembly/racon/1rd/${SampleID}_racon_1rd_f.fasta ../Assembly/racon/2rd/${SampleID}_racon_2rd_f.fasta ../Assembly/racon/3rd/${SampleID}_racon_3rd_f.fasta ../Assembly/medaka/1rd/${SampleID}_medaka_1rd_f.fasta > ./combined_fasta/${SampleID}_comb.fasta
# Alignments of all the assemblies using Muscle3
muscle -in ./combined_fasta/${SampleID}_comb.fasta -out ../OutputFiles/ntAlignment/${SampleID}_nt_align_muscle.html -html
## Align sam&medaka_1rd to Reference by muscle
# Combine all the fasta files in one file
cat ../InputFiles/references/${ReferenceName}.fasta ../Assembly/bedtools/sam/${SampleID}_sam_f.fasta ../Assembly/medaka/1rd/${SampleID}_medaka_1rd_f.fasta > ./combined_fasta/${SampleID}_comb.fasta
# Alignments of all the assemblies using Muscle3
muscle -in ./combined_fasta/${SampleID}_comb.fasta -out ../OutputFiles/ntAlignment/${SampleID}_nt_align_muscle_f.html -html

## align seqs with references individually by emboss
# Combine all the fasta files in one file
cat ../Assembly/bedtools/sam/${SampleID}_sam_f.fasta ../Assembly/bedtools/bed/${SampleID}_draft_f.fasta ../Assembly/racon/1rd/${SampleID}_racon_1rd_f.fasta ../Assembly/racon/2rd/${SampleID}_racon_2rd_f.fasta ../Assembly/racon/3rd/${SampleID}_racon_3rd_f.fasta ../Assembly/medaka/1rd/${SampleID}_medaka_1rd_f.fasta > ./combined_fasta/${SampleID}_comb.fasta
sed "s/\*/N/g" ./combined_fasta/${SampleID}_comb.fasta > ./combined_fasta/${SampleID}_comb_cleaned.fasta
# Alignments of all the assemblies using emboss
needle -asequence ../InputFiles/references/${ReferenceName}.fasta -bsequence ./combined_fasta/${SampleID}_comb_cleaned.fasta -gapopen 10 -gapextend 0.5 -outfile ../OutputFiles/ntAlignment/${SampleID}_nt_align_emboss.txt


if $AA_Validation; then
    ### Variant Identification - pairwise alignment - aa
    cd $OUTPUT_DIR/Analysis_Results/VariantIdentification
    
    ## Get all ORFs for the polished consensus sequences
    # Create folders for each sample
    mkdir -p $OUTPUT_DIR/Analysis_Results/VariantIdentification/AssembledORFs/{allORFs/${SampleID},selectedORF/${SampleID}}
    
    # sam
    sed "s/N//g; s/\*//g" ../Assembly/bedtools/sam/${SampleID}_sam_f.fasta | getorf -sequence stdin -outseq ./AssembledORFs/allORFs/${SampleID}/${SampleID}_sam_all_orf.fasta -minsize 3 -find 1 -reverse N 
    # draft
    sed "s/N//g; s/\*//g" ../Assembly/bedtools/bed/${SampleID}_draft_f.fasta | getorf -sequence stdin -outseq ./AssembledORFs/allORFs/${SampleID}/${SampleID}_draft_all_orf.fasta -minsize 3 -find 1 -reverse N
    # racon - 1rd
    sed "s/N//g; s/\*//g" ../Assembly/racon/1rd/${SampleID}_racon_1rd_f.fasta | getorf -sequence stdin -outseq ./AssembledORFs/allORFs/${SampleID}/${SampleID}_racon_1rd_all_orf.fasta -minsize 3 -find 1 -reverse N
    # racon - 2rd
    sed "s/N//g; s/\*//g" ../Assembly/racon/2rd/${SampleID}_racon_2rd_f.fasta | getorf -sequence stdin -outseq ./AssembledORFs/allORFs/${SampleID}/${SampleID}_racon_2rd_all_orf.fasta -minsize 3 -find 1 -reverse N
    # racon - 3rd
    sed "s/N//g; s/\*//g" ../Assembly/racon/3rd/${SampleID}_racon_3rd_f.fasta | getorf -sequence stdin -outseq ./AssembledORFs/allORFs/${SampleID}/${SampleID}_racon_3rd_all_orf.fasta -minsize 3 -find 1 -reverse N
    # medaka - 1rd
    sed "s/N//g; s/\*//g" ../Assembly/medaka/1rd/${SampleID}_medaka_1rd_f.fasta | getorf -sequence stdin -outseq ./AssembledORFs/allORFs/${SampleID}/${SampleID}_medaka_1rd_all_orf.fasta -minsize 3 -find 1 -reverse N

    ## Find the ORFs with the smallest start site:
    find ./AssembledORFs/allORFs/${SampleID}/ -name '*.fasta' | parallel --jobs 4 'awk '\''/^>/ {header = $0; match($0, /\[([0-9]+)/, m); start = m[1]; if (min_start == "" || start < min_start) {min_start = start; min_header = header; min_seq = "";}} /^[^>]/ {if (start == min_start) {min_seq = min_seq $0;}} END {output_file = "./AssembledORFs/selectedORF/" substr(FILENAME, length("./AssembledORFs/allORFs/") + 1); sub("_all_orf.fasta", "_orf.fasta", output_file); min_header = gensub(/_[0-9]+( \[.*\])/, "\\1", "1", min_header); print min_header "\n" min_seq > output_file;}'\'' {}'
    
    ## Protein sequences alignment to see if it is a missense mutation
    
    ## align all the seqs by muscle
    # Combine all the fasta files in one file
    cat ReferenceORFs/${ReferenceName}_orf.fasta AssembledORFs/selectedORF/${SampleID}/*.fasta > AssembledORFs/combinedORFs/${SampleID}_comb.fasta
    # Alignments of all the assemblies using Muscle3
    muscle -in ./AssembledORFs/combinedORFs/${SampleID}_comb.fasta -out ../OutputFiles/aaAlignment/${SampleID}_aa_align_muscle.html -html
    ## align sam&medaka_1rd to Reference by muscle
    # Combine all the fasta files in one file
    cat ReferenceORFs/${ReferenceName}_orf.fasta AssembledORFs/selectedORF/${SampleID}/${SampleID}_sam_orf.fasta  AssembledORFs/selectedORF/${SampleID}/${SampleID}_medaka_1rd_orf.fasta > AssembledORFs/combinedORFs/${SampleID}_comb.fasta
    # Alignments of all the assemblies using Muscle3
    muscle -in ./AssembledORFs/combinedORFs/${SampleID}_comb.fasta -out ../OutputFiles/aaAlignment/${SampleID}_aa_align_muscle_f.html -html
    
    ## align seqs with references individually by emboss
    # Combine all the fasta files in one file
    cat AssembledORFs/selectedORF/${SampleID}/*.fasta > AssembledORFs/combinedORFs/${SampleID}_comb.fasta
    # Alignments of all the assemblies using emboss
    needle -asequence ReferenceORFs/${ReferenceName}_orf.fasta -bsequence ./AssembledORFs/combinedORFs/${SampleID}_comb.fasta -gapopen 10 -gapextend 0.5 -outfile ../OutputFiles/aaAlignment/${SampleID}_aa_align_emboss.txt
fi
