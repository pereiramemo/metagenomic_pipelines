#!/bin/bash

###############################################################################
### 1. Load env
###############################################################################

set -o pipefail 
source "${HOME}/workspace/repositories/metagenomic_pipelines/preprocess/\
preprocess_pipeline_conf.bash"

###############################################################################
# 2. Define help
###############################################################################

show_usage(){
  cat <<EOF
Usage: ./quality_check_plots.bash <options>
--help                          print this help
--input_dir CHAR                directory with input fastq files
--output_dir CHAR               directory to output plots
--r1_pattern CHAR               pattern of R1 fastq files
--r2_pattern CHAR               pattern of R2 fastq files
--nslots NUM                    number of threads used (default 12)
--overwrite t|f                 overwrite previous output (default f)
EOF
}

###############################################################################
# 3. Parse input parameters
###############################################################################

while :; do
  case "${1}" in
    --help) # Call a "show_help" function to display a synopsis, then exit.
    show_usage
    exit 1;
    ;;
#############
  --input_dir)
  if [[ -n "${2}" ]]; then
    INPUT_DIR="${2}"
    shift
  fi
  ;;
  --input_dir=?*)
  INPUT_DIR="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --input_dir=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --output_dir)
  if [[ -n "${2}" ]]; then
    OUTPUT_DIR="${2}"
    shift
  fi
  ;;
  --output_dir=?*)
  OUTPUT_DIR="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --output_dir=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --nslots)
  if [[ -n "${2}" ]]; then
    NSLOTS="${2}"
    shift
  fi
  ;;
  --nslots=?*)
  NSLOTS="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --nslots=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --overwrite)
  if [[ -n "${2}" ]]; then
    OVERWRITE="${2}"
    shift
  fi
  ;;
  --overwrite=?*)
  OVERWRITE="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --overwrite=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;; 
#############
  --r1_pattern)
  if [[ -n "${2}" ]]; then
    R1_PATTERN="${2}"
    shift
  fi
  ;;
  --r1_pattern=?*)
  R1_PATTERN="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --r1_pattern=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;; 
#############
  --r2_pattern)
  if [[ -n "${2}" ]]; then
    R2_PATTERN="${2}"
    shift
  fi
  ;;
  --r2_pattern=?*)
  R2_PATTERN="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --r2_pattern=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;   
############ End of all options.
  --)       
  shift
  break
  ;;
  -?*)
  printf 'WARN: Unknown option (ignored): %s\n' "$1" >&2
  ;;
  *) # Default case: If no more options, then break out of the loop.
  break
  esac
  shift
done  

###############################################################################
# 4. Define defaults
###############################################################################

if [[ -z "${NSLOTS}" ]]; then
  NSLOTS="12"
fi

if [[ -z "${OVERWRITE}" ]]; then
  OVERWRITE="f"
fi

if [[ -z "${R1_PATTERN}" ]]; then
  R1_PATTERN="R1_001.fastq.gz"
fi  

if [[ -z "${R2_PATTERN}" ]]; then
  R2_PATTERN="R2_001.fastq.gz"
fi  

###############################################################################
# 5. Check input dir
###############################################################################

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "no input dir"
  exit 1
fi

###############################################################################
# 6. Check output 
###############################################################################

if [[ -d "${OUTPUT_DIR}" && "${OVERWRITE}" == f  ]]; then
  echo "output directory ${OUTPUT_DIR} already exists"
  echo "use \"--overwrite t\" to overwrite"
  exit 0
fi  
  
if [[ -d "${OUTPUT_DIR}" && "${OVERWRITE}" == t  ]]; then  
  rm -r "${OUTPUT_DIR}"/*.png
fi  
  
if [[ ! -d "${OUTPUT_DIR}" ]]; then  
  mkdir "${OUTPUT_DIR}"
fi    

###############################################################################
# 7. Run rscript
###############################################################################

Rscript --vanilla \
"${quality_check_plots}" \
"${INPUT_DIR}" \
"${OUTPUT_DIR}" \
"${NSLOTS}" \
"${R1_PATTERN}" \
"${R2_PATTERN}"

if [[ $? != "0" ]]; then
  echo "quality_check_plots.R failed"
  rm -r "${OUTPUT_DIR}"/*.png
  exit 1
fi

