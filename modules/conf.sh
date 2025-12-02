###############################################################################
# conf.sh
# Configuration file for metagenomic pipeline
###############################################################################


###############################################################################
# Load conda environment
###############################################################################

# Check if metagenomic_pipeline environment is activated, if not, activate it
if [[ "${CONDA_DEFAULT_ENV}" != "metagenomic_pipeline" ]]; then
  echo "Activating metagenomic_pipeline conda environment..."
  
  # Try to find conda/mamba
  if command -v mamba &> /dev/null; then
    eval "$(mamba shell.bash hook)"
    mamba activate metagenomic_pipeline
  elif command -v conda &> /dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate metagenomic_pipeline
  else
    echo "ERROR: Neither mamba nor conda found in PATH"
    echo "Please install mamba or conda and ensure it's in your PATH"
    exit 1
  fi
  
  # Verify activation was successful
  if [[ "${CONDA_DEFAULT_ENV}" != "metagenomic_pipeline" ]]; then
    echo "ERROR: Failed to activate metagenomic_pipeline environment"
    echo "Please ensure the environment exists. Create it with:"
    echo "  mamba env create -f environment.yml"
    exit 1
  fi
  
  echo "Environment activated successfully"
fi

###############################################################################
# dirs
###############################################################################

# Get the directory where this script is located
MODULES_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# tools (using conda environment)
bbduk="bbduk.sh"
bbmerge="bbmerge.sh"
seqtk="seqtk"
pear="pear"
bzip2="bzip2"
gunzip="gunzip"
rscript="Rscript"
pigz="pigz"
fq2fa="${MODULES_DIR}/resources/fq2fa.sh"
plots="${MODULES_DIR}/resources/plots.R"
quality_check="${MODULES_DIR}/0-quality_check.R"

###############################################################################
# files
###############################################################################

# BBMap adapters file location in conda environment
ADAPTERS=$(find "${CONDA_PREFIX}/opt/bbmap"* -name "adapters.fa")
if [[ ! -f "${ADAPTERS}" ]]; then
  echo "ERROR: BBMap adapters file not found at expected location: ${ADAPTERS}"
  echo "Please ensure BBMap is installed in the metagenomic_pipeline conda environment."
fi

###############################################################################
# functions
###############################################################################

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
    # Check tidyverse
    if ! Rscript -e "library(tidyverse)" &> /dev/null; then
      missing_tools+=("r-tidyverse")
    fi
    # Check ShortRead (Bioconductor)
    if ! Rscript -e "library(ShortRead)" &> /dev/null; then
      missing_tools+=("bioconductor-shortread")
    fi
    # Check doParallel
    if ! Rscript -e "library(doParallel)" &> /dev/null; then
      missing_tools+=("r-doparallel")
    fi
    # Check dada2 (Bioconductor)
    if ! Rscript -e "library(dada2)" &> /dev/null; then
      missing_tools+=("bioconductor-dada2")
    fi
    # Check optparse
    if ! Rscript -e "library(optparse)" &> /dev/null; then
      missing_tools+=("r-optparse")
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
  if [[ ! -f "${quality_check}" ]]; then
    missing_tools+=("quality_check.R")
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

