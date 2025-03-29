#!/bin/bash

### Define options
# Define the default values
SCRIPT_DIR=""
INPUT_DIR=""
OUTPUT_DIR=""
OPTIMIZE_eE=false
e_Pro=8
E_Pro=11
AA_Validation=false
q_Chopper=13
Medaka_Mod="r1041_e82_400bps_hac_v4.3.0"

# Function to display help
usage() {
  echo "Usage:"
  echo "./ConstructValidation.sh --script <DIR> -i <DIR> -o <DIR> [options]"
  echo "Options:"
  echo "  --script                      <DIR>   Scripts directory"
  echo "  -i, --input                   <DIR>   Input directory"
  echo "  -o, --output                  <DIR>   Output directory"
  echo "  -q_chopper                    <INT>   Sets a minimum average Phred quality score for reads filtered by Chopper (default:13)"
  echo "  -e                            <INT>   Barcode edit distance value (default:8, only if optimizing_eE is false)"
  echo "  -E                            <INT>   Primer edit distance value (default:11, only if optimizing_eE is false)"
  echo "  --medaka_model                <STR>   Base-calling model used for Medaka polishing (default:r1041_e82_400bps_hac_v4.3.0)"
  echo "  --optimizing_eE                       Use the built-in function to optimize demultiplexing parameters (e & E) for maximum mapped reads"
  echo "  --amino_acid_seq_validation           Translate the first ORF to its amino acid sequence to detect missense/silent mutations"
  echo "  -h, --help                            Display this help"
  exit 1
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --script) SCRIPT_DIR="$2"; shift ;;
        -i|--input) INPUT_DIR="$2"; shift ;;
        -o|--output) OUTPUT_DIR="$2"; shift ;;
        -q_chopper) q_Chopper="$2"; shift ;;
        --medaka_model) Medaka_Mod="$2"; shift ;;
        -e) 
            if [ "$OPTIMIZE_eE" = true ]; then
                echo "Error: Cannot define parameter e when optimizing_eE is true."
                exit 1
            fi
            e_Pro="$2"
            shift ;;
        -E)
            if [ "$OPTIMIZE_eE" = true ]; then
                echo "Error: Cannot define parameter E when optimizing_eE is true."
                exit 1
            fi
            E_Pro="$2"
            shift ;;
        --optimizing_eE) OPTIMIZE_eE=true ;;
        --amino_acid_seq_validation) AA_Validation=true ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if mandatory parameters are set
if [[ -z "$INPUT_DIR" || -z "$OUTPUT_DIR" ]]; then
    echo "Error: Input and output directories must be specified."
    usage
fi

### Making directory structure
mkdir -p $OUTPUT_DIR/Analysis_Results/{Demultiplexing,VariantCalling/{beforePolishing,afterPolishing},Alignment,Assembly,PairwiseAlignment,InputFiles/references,OutputFiles/ntAlignment}
mkdir -p $OUTPUT_DIR/Analysis_Results/Assembly/{bedtools/{sam,bed},racon/{1rd,2rd,3rd,sam},medaka/1rd}
mkdir -p $OUTPUT_DIR/Analysis_Results/VariantCalling/afterPolishing/{racon/{1rd,2rd,3rd},medaka_1rd}
mkdir -p $OUTPUT_DIR/Analysis_Results/PairwiseAlignment/combined_fasta


### Generate InputFiles
chmod 775 $SCRIPT_DIR/InputFiles.R
Rscript $SCRIPT_DIR/InputFiles.R "$INPUT_DIR" "$OUTPUT_DIR" "$AA_Validation"

if $AA_Validation; then
    mkdir -p $OUTPUT_DIR/Analysis_Results/{VariantIdentification/{ReferenceORFs,AssembledORFs/{allORFs,selectedORF,combinedORFs}},OutputFiles/aaAlignment}
    ### Variant Identification - pairwise alignment - aa - Get ORF for reference
    parallel --jobs 100 --colsep '\t' "getorf -sequence ${OUTPUT_DIR}/Analysis_Results/InputFiles/references/{1}.fasta -outseq ${OUTPUT_DIR}/Analysis_Results/VariantIdentification/ReferenceORFs/{1}_orf.fasta -minsize {2} -find 1 -reverse N" :::: <(tail -n +2 "$OUTPUT_DIR/Analysis_Results/InputFiles/references_info.tsv")
fi

### Set the working directory
cd $OUTPUT_DIR/Analysis_Results

### Demultiplexing

cd $OUTPUT_DIR/Analysis_Results/Demultiplexing
chmod 775 $SCRIPT_DIR/minibar.py

if $OPTIMIZE_eE; then
    ## Get the edit distance between the most similar barcodes and primers (smallest edit distance).
    # edit distance of closest primers
    edP=$($SCRIPT_DIR/minibar.py $OUTPUT_DIR/Analysis_Results/InputFiles/IndexCombination.txt -info primer | grep -o 'edit distance [0-9]*' | awk '{print $3}')
    # length of the primers
    maxP=$($SCRIPT_DIR/minibar.py $OUTPUT_DIR/Analysis_Results/InputFiles/IndexCombination.txt -info primer | grep -oP 'with lengths of \K[0-9]+')
    # edit distance of closest barcodes
    edB=$($SCRIPT_DIR/minibar.py $OUTPUT_DIR/Analysis_Results/InputFiles/IndexCombination.txt -info both | grep -o 'edit distance of [0-9]*' | awk '{print $4}')
    
    ## Apply all possible combinations of e & E
    echo -e "e\tE\tMulti\tUnknown\tSorted\tTotalDmtplx\tSorted_pct\tTotalMapped\tMapped_pct\tTotalSamples\tTotalCompleteAssembledSamples\tTotalSets\tTotalCompleteAssembledSets\tSets_MeanDepth_5\tSamples_MeanDepth_5\tSets_MeanDepth_10\tSamples_MeanDepth_10\tSets_MeanDepth_15\tSamples_MeanDepth_15\tSets_MeanDepth_20\tSamples_MeanDepth_20" > e_n_E_Combination.tsv
    ARRAY_TASK_ID=0
    for i in $(seq 1 $edB); do
      for j in $(seq 1 $maxP); do
        ARRAY_TASK_ID=$((ARRAY_TASK_ID+1))
        chmod 775 $SCRIPT_DIR/Demultiplexing_e_E.sh
        $SCRIPT_DIR/Demultiplexing_e_E.sh "$SCRIPT_DIR" "$INPUT_DIR" "$OUTPUT_DIR" "$q_Chopper" "$ARRAY_TASK_ID" "$i" "$j"
      done
    done

    ## Extract the best combination of e & E
    # Step 1: Sort the TSV file by TotalMapped and then TotalSets (descending) and get the top row
    top_row=$(sort -t$'\t' -k19,19nr -k8,8nr e_n_E_Combination.tsv | head -n 1)
    # Step 2: Extract the e and E columns from the top row
    e_slt=$(echo "$top_row" | cut -f1)
    E_slt=$(echo "$top_row" | cut -f2)

else
    # If not optimizing, use predefined e and E values
    e_slt="$e_Pro"
    E_slt="$E_Pro"
fi


## Sort the reads using minibar
mkdir final
cd final
$SCRIPT_DIR/minibar.py ../../InputFiles/IndexCombination.txt $INPUT_DIR/passed_all.fastq -e $e_slt -E $E_slt -T -F -P ""

## Filtering with chopper
#conda init
#conda activate chopper_env
parallel 'chopper -q '$q_Chopper' -i {} > {.}_filtered.fastq' ::: *.fastq
#conda deactivate

### Run analysis for individual samples

cd $OUTPUT_DIR/Analysis_Results
config=$OUTPUT_DIR/Analysis_Results/InputFiles/config_sample.txt
maxArrayid=$(awk 'NR > 1 { if ($1 > max) max = $1 } END { print max }' "$OUTPUT_DIR/Analysis_Results/InputFiles/config_sample.txt")
chmod 775 $SCRIPT_DIR/Construct_Validation_per_Sample.sh
for ArrayTaskID in $(seq 1 $maxArrayid); do 
  # SampleID
  SampleID=$(awk -v ArrayTaskID=$ArrayTaskID '$1==ArrayTaskID {print $2}' $config)
  # ReferenceName
  Reference=$(awk -v ArrayTaskID=$ArrayTaskID '$1==ArrayTaskID {print $3}' $config)
  # Run the script
  $SCRIPT_DIR/Construct_Validation_per_Sample.sh "$SCRIPT_DIR" "$INPUT_DIR" "$OUTPUT_DIR" "$AA_Validation" "$Medaka_Mod" "$SampleID" "$Reference"
done

## Run the Summary.sh
cd $OUTPUT_DIR/Analysis_Results
chmod 775 $SCRIPT_DIR/Summary.sh
$SCRIPT_DIR/Summary.sh "$SCRIPT_DIR" "$INPUT_DIR" "$OUTPUT_DIR" "$AA_Validation"

