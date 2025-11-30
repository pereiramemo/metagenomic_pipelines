# Load conda environment
# Make sure the metagenomic_pipeline environment is activated before running the pipeline
# Run: mamba activate metagenomic_pipeline

# dirs
# Get the directory where this script is located
PREPROCESS_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# tools (using conda environment)
bbduk="bbduk.sh"
bbmerge="bbmerge.sh"
seqtk="seqtk"
pear="pear"
bzip2="bzip2"
gunzip="gunzip"
fq2fa="${PREPROCESS_DIR}/resources/fq2fa.sh"
plots="${PREPROCESS_DIR}/resources/plots.R"
quality_check_plots="${PREPROCESS_DIR}/resources/quality_check_plots.R"

# files
# BBMap adapters file location in conda environment
ADAPTERS="${CONDA_PREFIX}/opt/bbmap*/resources/adapters.fa"

# Check that all required tools are installed
function check_dependencies {
  local missing_tools=()
  
  # Check command-line tools
  for tool in bbduk.sh bbmerge.sh seqtk pear bzip2 gunzip; do
    if ! command -v ${tool} &> /dev/null; then
      missing_tools+=("${tool}")
    fi
  done
  
  # Check R packages
  if command -v Rscript &> /dev/null; then
    if ! Rscript -e "library(tidyverse)" &> /dev/null; then
      missing_tools+=("r-tidyverse")
    fi
  else
    missing_tools+=("R")
  fi
  
  # Check local resources
  if [[ ! -f "${fq2fa}" ]]; then
    missing_tools+=("fq2fa.sh")
  fi
  if [[ ! -f "${plots}" ]]; then
    missing_tools+=("plots.R")
  fi
  if [[ ! -f "${quality_check_plots}" ]]; then
    missing_tools+=("quality_check_plots.R")
  fi
  
  # Report missing tools
  if [[ ${#missing_tools[@]} -gt 0 ]]; then
    echo "ERROR: The following required tools are not installed or not in PATH:"
    printf '  - %s\n' "${missing_tools[@]}"
    echo ""
    echo "Please activate the metagenomic_pipeline conda environment:"
    echo "  mamba activate metagenomic_pipeline"
    return 1
  fi
  
  return 0
}

# Run dependency check
check_dependencies || exit 1

# functions
function count_fastq {

  N=$(wc -l  "${1}" | cut -f1 -d" ")
  N=$( echo "${N}" / 4 | bc -l )
  echo "${N}"

}

function count_fasta {

  N=$(egrep -c ">" "${1}")
  echo "${N}"

}

function compute_mean_length {

  N=$(infoseq "${1}" | \
      awk 'NR > 1 {tot = $6 + tot; n++} END {print tot/n}')
  echo "${N}"

}

