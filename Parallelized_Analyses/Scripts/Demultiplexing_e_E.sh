#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2     # This will be the number of CPUs per individual array job
#SBATCH --mem=20G     # This will be the memory per individual array job
#SBATCH --time=0-06:00:00     # 4 hrs
#SBATCH --job-name="demtplx_para_test"

## Set variables
SCRIPT_DIR=$1
INPUT_DIR=$2
OUTPUT_DIR=$3
q_Chopper=$4
# Specify the path to the config file
config=$OUTPUT_DIR/Analysis_Results/Demultiplexing/config.txt
# Extract the e for the current $SLURM_ARRAY_TASK_ID
i=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
# Extract the E for the current $SLURM_ARRAY_TASK_ID
j=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)

# Load modules
module load minimap2 samtools bedtools parallel # Alignment & Assembly

## Create required folders
mkdir -p $OUTPUT_DIR/Analysis_Results/Demultiplexing/$SLURM_ARRAY_TASK_ID/{rawReads,filteredReads,tempMapping}
mkdir -p $OUTPUT_DIR/Analysis_Results/Demultiplexing/$SLURM_ARRAY_TASK_ID/tempMapping/{sam,bam,sorted_bam}

## Demultiplexing
cd $OUTPUT_DIR/Analysis_Results/Demultiplexing/$SLURM_ARRAY_TASK_ID/rawReads
chmod 775 $SCRIPT_DIR/minibar.py
$SCRIPT_DIR/minibar.py $OUTPUT_DIR/Analysis_Results/InputFiles/IndexCombination.txt $INPUT_DIR/passed_all.fastq -e $i -E $j -T -F -P ""

## Filtering with chopper
# Create the chopper_env using conda
conda init
source ~/bigdata/.conda/envs/chopper_env/bin/chopper
conda activate chopper_env
parallel 'chopper -q '$q_Chopper' -i {} > ../filteredReads/{.}_filtered.fastq' ::: ./*.fastq
cd ..
## Export the results for demultiplexing -  count the number of reads after demultiplexing
# Create or overwrite the TSV file with the header
echo -e "FileName\tReads" > temp_readsum.tsv
# Use GNU Parallel to process each FASTQ file
parallel 'echo -e {/}"\t"$(cat {} | wc -l | awk "{{print int(\$1/4)}}")' ::: ./filteredReads/*_filtered.fastq | sort >> temp_readsum.tsv

### Alignment
## Change the directory and create required directory:
cd tempMapping
## Alignment suing minimap2
parallel --colsep '\t' 'output=$(basename {1} _filtered.fastq).sam; minimap2 -ax map-ont -N 0 --secondary=no {2} {1} > "./sam/$output"' :::: <(awk -F'\t' 'NR > 1 {print "../filteredReads/" $1 "_filtered.fastq\t../../../InputFiles/references/" $2}' $OUTPUT_DIR/Analysis_Results/InputFiles/references.tsv)
# Convert SAM to BAM
parallel 'bname=$(basename {} .sam); samtools view -q 60 -F 2048 -S -b {} -o ./bam/${bname}.bam' ::: ./sam/*.sam
# Sort BAM files
parallel 'bname=$(basename {} .bam); samtools sort {} -o ./sorted_bam/${bname}.sorted.bam' ::: ./bam/*.bam
# Index Sorted BAM files
# parallel 'samtools index {}' ::: ./sorted_bam/*.sorted.bam
cd ..
## Export alignment results using samtools:
echo -e "FileName\tReads_Mapped\tCoverage_Percentage\tMean_Depth" > temp_align.tsv # Print the header (adjust according to the samtools coverage output)
parallel "samtools coverage {} | awk -v OFS='\t' -v file={} 'NR>1 {print file, \$4, \$6, \$7}' >> temp_align.tsv" ::: ./tempMapping/sorted_bam/*.sorted.bam

### Extract the results and stored it in corresponding variables
Multi=$(awk -F'\t' '$1 == "Multiple_Matches_filtered.fastq" {print $2}' temp_readsum.tsv)
Unknown=$(awk -F'\t' '$1 == "unk_filtered.fastq" {print $2}' temp_readsum.tsv)
TotalDmtplx=$(awk -F'\t' 'NR > 1 {sum += $2} END {print sum}' temp_readsum.tsv)
Sorted=$((TotalDmtplx - Multi - Unknown))
Sorted_pct=$(echo "scale=4; $Sorted / $TotalDmtplx" | bc)
TotalMapped=$(awk -F'\t' 'NR > 1 {sum += $2} END {print sum}' temp_align.tsv)
Mapped_pct=$(echo "scale=4; $TotalMapped / $TotalDmtplx" | bc)
TotalSamples=$(awk -F '\t' 'NR>1 && $2 != 0 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort |wc -l)
TotalCompleteAssembledSamples=$(awk -F '\t' 'NR>1 && $2>2 && $3>95 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort |wc -l)
TotalSets=$(awk -F '\t' 'NR>1 && $2 != 0 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort | uniq |wc -l)
TotalCompleteAssembledSets=$(awk -F '\t' 'NR>1 && $2>2 && $3>95 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort | uniq |wc -l)
Sets_MeanDepth_5=$(awk -F '\t' 'NR>1 && $4>=5 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort | uniq |wc -l)
Samples_MeanDepth_5=$(awk -F '\t' 'NR>1 && $4>=5 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort |wc -l)
Sets_MeanDepth_10=$(awk -F '\t' 'NR>1 && $4>=10 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort | uniq |wc -l)
Samples_MeanDepth_10=$(awk -F '\t' 'NR>1 && $4>=10 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort |wc -l)
Sets_MeanDepth_15=$(awk -F '\t' 'NR>1 && $4>=15 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort | uniq |wc -l)
Samples_MeanDepth_15=$(awk -F '\t' 'NR>1 && $4>=15 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort |wc -l)
Sets_MeanDepth_20=$(awk -F '\t' 'NR>1 && $4>=20 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort | uniq |wc -l)
Samples_MeanDepth_20=$(awk -F '\t' 'NR>1 && $4>=20 {print $1}' temp_align.tsv | sed 's#.*/\(set[0-9]*\)_[0-9]*\.sorted\.bam#\1#' | sort |wc -l)

### Export the value of all variables to the .tsv file
echo -e "$i\t$j\t$Multi\t$Unknown\t$Sorted\t$TotalDmtplx\t$Sorted_pct\t$TotalMapped\t$Mapped_pct\t$TotalSamples\t$TotalCompleteAssembledSamples\t$TotalSets\t$TotalCompleteAssembledSets\t$Sets_MeanDepth_5\t$Samples_MeanDepth_5\t$Sets_MeanDepth_10\t$Samples_MeanDepth_10\t$Sets_MeanDepth_15\t$Samples_MeanDepth_15\t$Sets_MeanDepth_20\t$Samples_MeanDepth_20" >> $OUTPUT_DIR/Analysis_Results/Demultiplexing/e_n_E_Combination.tsv

rm -rf $OUTPUT_DIR/Analysis_Results/Demultiplexing/$SLURM_ARRAY_TASK_ID


