#!/bin/bash

###############################################################################
### 1. Load env
###############################################################################

set -o pipefail 
# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/conf.sh"

###############################################################################
# 2. Define help
###############################################################################

show_usage(){
  cat <<EOF
Usage: ./quality_check_plots_runner.sh <options>
--help                          print this help
--input_dir CHAR                directory with input fastq files
--output_dir CHAR               directory to output plots
--r1_pattern CHAR               pattern of R1 fastq files (default: R1_001.fastq.gz)
--r2_pattern CHAR               pattern of R2 fastq files (default: R2_001.fastq.gz)
--nslots NUM                    number of threads used (default: 12)
--overwrite TRUE|FALSE          overwrite previous output (default: FALSE)

Note: All validation and directory creation is handled by the R script.
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
  OVERWRITE="FALSE"
fi

if [[ -z "${R1_PATTERN}" ]]; then
  R1_PATTERN="R1_001.fastq.gz"
fi  

if [[ -z "${R2_PATTERN}" ]]; then
  R2_PATTERN="R2_001.fastq.gz"
fi  

###############################################################################
# 5. Run rscript
###############################################################################

Rscript --vanilla \
"${quality_check_plots}" \
--input_dir "${INPUT_DIR}" \
--output_dir "${OUTPUT_DIR}" \
--nslots "${NSLOTS}" \
--r1_pattern "${R1_PATTERN}" \
--r2_pattern "${R2_PATTERN}" \
--overwrite "${OVERWRITE}"

if [[ $? != "0" ]]; then
  echo "quality_check_plots.R failed"
  exit 1
fi

