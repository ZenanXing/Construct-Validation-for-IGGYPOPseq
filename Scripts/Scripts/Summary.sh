#!/bin/bash

## Set variables
SCRIPT_DIR=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
AA_Validation=$4

### Change the directory
cd $OUTPUT_DIR/Analysis_Results

### Demultiplexing -  count the number of reads after demultiplexing
# Create or overwrite the TSV file with the header
echo -e "FileName\tReads" > OutputFiles/Demultiplexing_ReadsSummary.tsv
# Use GNU Parallel to process each FASTQ file
parallel 'echo {/}"\t"$(($(wc -l < {}) / 4))' ::: ./Demultiplexing/final/*_filtered.fastq | sort >> OutputFiles/Demultiplexing_ReadsSummary.tsv

apt-get install pandoc
### Alignment
## Export alignment results using samtools:
echo -e "FileName\tReference\tStart\tEnd\tReads_Mapped\tCovered_Bases\tCoverage_Percentage\tMean_Depth\tMean_Base_Quality\tMean_Mapping_Quality" > OutputFiles/Alignment.tsv # Print the header (adjust according to the samtools coverage output)
parallel "samtools coverage {} | awk -v OFS='\t' -v file={} 'NR>1 {print file, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9}' >> OutputFiles/Alignment.tsv" ::: ./Alignment/*.sorted.bam


### Assembly
## Drafts
# delete the fasta files that is 0 bt
#find ./bedtools/sam -type f -name "*.fasta" -size 0 -delete
#find ./bedtools/bed -type f -name "*.fasta" -size 0 -delete
# export results of the draft assembly
echo -e "FileName\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts" > OutputFiles/Assembled_bP.tsv
parallel "seqtk comp {} | awk -v OFS='\t' -v file={} '{print file, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13}' >> OutputFiles/Assembled_bP.tsv" ::: ./Assembly/bedtools/bed/*.fasta
sort -t$'\t' -k1 OutputFiles/Assembled_bP.tsv

## Polishing by Racon
# racon polishing - 1rd
# delete the fasta files that is 0 bt
# find ./racon/1rd -type f -name "*.fasta" -size 0 -delete
# export results of the assembly after the first round of polishing
echo -e "FileName\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts" > OutputFiles/Racon_1rdP.tsv
parallel "seqtk comp {} | awk -v OFS='\t' -v file={} '{print file, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13}' >> OutputFiles/Racon_1rdP.tsv" ::: ./Assembly/racon/1rd/*.fasta
sort -t$'\t' -k1 OutputFiles/Racon_1rdP.tsv
# racon polishing - 2rd
# delete the fasta files that is 0 bt
# find ./racon/2rd -type f -name "*.fasta" -size 0 -delete
# export results of the assembly after the second round of polishing
echo -e "FileName\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts" > OutputFiles/Racon_2rdP.tsv
parallel "seqtk comp {} | awk -v OFS='\t' -v file={} '{print file, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13}' >> OutputFiles/Racon_2rdP.tsv" ::: ./Assembly/racon/2rd/*.fasta
sort -t$'\t' -k1 ../OutputFiles/Racon_2rdP.tsv

## racon polishing - 3rd
# delete the fasta files that is 0 bt
# find ./racon/3rd -type f -name "*.fasta" -size 0 -delete
# export results of the assembly after the third round of polishing
echo -e "FileName\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts" > OutputFiles/Racon_3rdP.tsv
parallel "seqtk comp {} | awk -v OFS='\t' -v file={} '{print file, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13}' >> OutputFiles/Racon_3rdP.tsv" ::: ./Assembly/racon/3rd/*.fasta
sort -t$'\t' -k1 OutputFiles/Racon_3rdP.tsv

## Polishing with medaka
# medaka polishing - 1rd
# delete the fasta files that is 0 bt
# find ./medaka/1rd -type f -name "*.fasta" -size 0 -delete
# export results of the assembly after the first round of polishing
echo -e "FileName\tchr\tlength\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts" > OutputFiles/Medaka_1rdP.tsv
parallel "seqtk comp {} | awk -v OFS='\t' -v file={} '{print file, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$9, \$10, \$11, \$12, \$13}' >> OutputFiles/Medaka_1rdP.tsv" ::: ./Assembly/medaka/1rd/*.fasta
sort -t$'\t' -k1 OutputFiles/Medaka_1rdP.tsv


### Variant Calling
## after racon polishing
# 1st round
echo -e "SampleID\tReference\tType\tPosition\tRef_allele\tAlt_allele" > OutputFiles/VariantCalling_racon_1rd.tsv
parallel 'bname=$(basename {} _racon_1rd.all.vcf); bcftools query -f "%CHROM\t%TYPE\t%POS\t%REF\t%ALT\n" {} | awk -v bname="$bname" "{print bname \"\t\" \$0}"' ::: ./VariantCalling/afterPolishing/racon/1rd/*.all.vcf >> OutputFiles/VariantCalling_racon_1rd.tsv
# 2nd round
echo -e "SampleID\tReference\tType\tPosition\tRef_allele\tAlt_allele" > OutputFiles/VariantCalling_racon_2rd.tsv
parallel 'bname=$(basename {} _racon_2rd.all.vcf); bcftools query -f "%CHROM\t%TYPE\t%POS\t%REF\t%ALT\n" {} | awk -v bname="$bname" "{print bname \"\t\" \$0}"' ::: ./VariantCalling/afterPolishing/racon/2rd/*.all.vcf >> OutputFiles/VariantCalling_racon_2rd.tsv
# 3rd round
echo -e "SampleID\tReference\tType\tPosition\tRef_allele\tAlt_allele" > OutputFiles/VariantCalling_racon_3rd.tsv
parallel 'bname=$(basename {} _racon_3rd.all.vcf); bcftools query -f "%CHROM\t%TYPE\t%POS\t%REF\t%ALT\n" {} | awk -v bname="$bname" "{print bname \"\t\" \$0}"' ::: ./VariantCalling/afterPolishing/racon/3rd/*.all.vcf >> OutputFiles/VariantCalling_racon_3rd.tsv
## after medaka polishing
# 1st round
echo -e "SampleID\tReference\tType\tPosition\tRef_allele\tAlt_allele" > OutputFiles/VariantCalling_medaka_1rd.tsv
parallel 'bname=$(basename {} _medaka_1rd.all.vcf); bcftools query -f "%CHROM\t%TYPE\t%POS\t%REF\t%ALT\n" {} | awk -v bname="$bname" "{print bname \"\t\" \$0}"' ::: ./VariantCalling/afterPolishing/medaka_1rd/*.all.vcf >> OutputFiles/VariantCalling_medaka_1rd.tsv


### Pairwise alignment - nt
# Initialize the TSV file with headers
echo -e "SampleID\tPolishingStep\tIdentity_nt" > OutputFiles/ntAlignment.tsv
parallel -j 100 --no-notice "
    process_file() {
        file=\$1
        echo \"Processing file: \$file\"  # Print the current file being processed
        # Extract SampleName - only remove the known suffix
        sample_name=\$(basename \"\$file\" | sed -E 's/_nt_align_emboss.txt\$//')
        # Loop through polishing steps
        for step in sam draft racon_1rd racon_2rd racon_3rd medaka_1rd; do
            polishing_step=\$(grep -oP \"\\# 2\\:.*_\\K\${step}\" \"\$file\")
            identity=\$(awk \"/\${step}/{flag=1} flag\" \"\$file\" | grep -m 1 -oP '# Identity:\\s+\\K.*')
            # Append the extracted information to the TSV file
            if [ -n \"\$sample_name\" ] && [ -n \"\$polishing_step\" ] && [ -n \"\$identity\" ]; then
                # Use echo without -e to prevent unwanted characters
                echo \"\$sample_name\t\$polishing_step\t\$identity\" >> OutputFiles/ntAlignment.tsv
            else
                echo \"Warning: Some data not found for file: \$file step: \$step\"  # Log missing data
            fi
        done
    }
    process_file {}
" ::: OutputFiles/ntAlignment/*_emboss.txt

### Visualization of the pairwise alignment - nt - samplewise
cd $OUTPUT_DIR/Analysis_Results/PairwiseAlignment
mkdir clones
## align all the seqs by muscle
# Combine all the fasta files in one file
parallel --colsep '\t' 'cat ../InputFiles/references/{2}.fasta ../Assembly/medaka/1rd/{1}_*_medaka_1rd_f.fasta > ./clones/{1}_comb.fasta' :::: <(tail -n +2 ../InputFiles/primer_index_info.tsv)
# Alignments of all the assemblies using Muscle3
parallel 'bname=$(basename {} _comb.fasta); muscle -in {} -out ../OutputFiles/ntAlignment/${bname}_nt_align_muscle_clones.html -html' ::: ./clones/*.fasta
cd ..

### Variant Identification - pairwise alignment - aa
if $AA_Validation; then
    # Initialize the TSV file with headers
    echo -e "SampleID\tPolishingStep\tIdentity_aa" > OutputFiles/aaAlignment.tsv
    parallel -j 100 --no-notice "
    process_file() {
        file=\$1
        echo \"Processing file: \$file\"  # Print the current file being processed
        # Extract SampleName - only remove the known suffix
        sample_name=\$(basename \"\$file\" | sed -E 's/_aa_align_emboss.txt\$//')
        # Loop through polishing steps
        for step in sam draft racon_1rd racon_2rd racon_3rd medaka_1rd; do
            polishing_step=\$(grep -oP \"\\# 2\\:.*_\\K\${step}\" \"\$file\")
            identity=\$(awk \"/\${step}/{flag=1} flag\" \"\$file\" | grep -m 1 -oP '# Identity:\\s+\\K.*')
            # Append the extracted information to the TSV file
            if [ -n \"\$sample_name\" ] && [ -n \"\$polishing_step\" ] && [ -n \"\$identity\" ]; then
                # Use echo without -e to prevent unwanted characters
                echo \"\$sample_name\t\$polishing_step\t\$identity\" >> OutputFiles/aaAlignment.tsv
            else
                echo \"Warning: Some data not found for file: \$file step: \$step\"  # Log missing data
            fi
        done
    }
    process_file {}
    " ::: OutputFiles/aaAlignment/*_emboss.txt
    
    # Export the length of ORFs
    # Reference
    echo -e "FileName\tchr\tRefORFlength" > OutputFiles/Ref_ORF_length.tsv
    parallel "seqtk comp {} | awk -v OFS='\t' -v file={} '{print file, \$1, \$2}' >> OutputFiles/Ref_ORF_length.tsv" ::: ./VariantIdentification/ReferenceORFs/*.fasta
    sort -t$'\t' -k1 OutputFiles/Ref_ORF_length.tsv
    # Assembly
    echo -e "FileName\tchr\tORFlength" > OutputFiles/ORF_length.tsv
    parallel "seqtk comp {} | awk -v OFS='\t' -v file={} '{print file, \$1, \$2}' >> OutputFiles/ORF_length.tsv" ::: ./VariantIdentification/AssembledORFs/selectedORF/*/*.fasta
    sort -t$'\t' -k1 OutputFiles/ORF_length.tsv
    
    # Visualization of the pairwise alignment - aa - samplewise
    cd $OUTPUT_DIR/Analysis_Results/VariantIdentification
    mkdir -p AssembledORFs/clones
    ## align all the seqs by muscle
    # Combine all the fasta files in one file
    parallel --colsep '\t' 'bname=$(basename {2} .fasta); cat ReferenceORFs/${bname}_orf.fasta AssembledORFs/selectedORF/{1}_*/{1}_*_medaka_1rd_orf.fasta > AssembledORFs/clones/{1}_comb.fasta' :::: <(tail -n +2 ../InputFiles/primer_index_info.tsv)
    # Alignments of all the assemblies using Muscle3
    parallel 'bname=$(basename {} _comb.fasta); muscle -in {} -out ../OutputFiles/aaAlignment/${bname}_aa_align_muscle_clones.html -html' ::: ./AssembledORFs/clones/*.fasta
    cd ..
 
fi

## Export the results
chmod 775 $SCRIPT_DIR/Summary.R
Rscript $SCRIPT_DIR/Summary.R "$INPUT_DIR" "$OUTPUT_DIR" "$AA_Validation"
# Report.html
chmod 775 $SCRIPT_DIR/Report.Rmd
Rscript -e "rmarkdown::render('$SCRIPT_DIR/Report.Rmd', output_dir = '$OUTPUT_DIR/Analysis_Results/OutputFiles', params = list(input_dir = '$INPUT_DIR', output_dir = '$OUTPUT_DIR', aa_validation = '$AA_VALIDATION'))"
