#!/usr/bin/env bash
###############################################################################
### mg-clust_module-1.bash
### De novo assembly + read mapping + duplicate removal
###############################################################################

set -euo pipefail

###############################################################################
### 1. Set env
###############################################################################

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Load configuration (must define paths to tools: megahit, bwa, samtools, picard, etc.)
source "${SCRIPT_DIR}/conf.sh"

###############################################################################
### 1.1 Helper functions
###############################################################################

RED='\033[0;31m'
YELLOW='\033[1;33m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

log()        { echo -e "[INFO] $*"; }
log_warn()   { echo -e "${YELLOW}[WARN]${NC} $*" >&2; }
log_error()  { echo -e "${RED}[ERROR]${NC} $*" >&2; }

# Check that a command exists in PATH
check_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    log_error "Required command '$1' not found in PATH."
    exit 1
  fi
}

# Check that a file exists
check_file() {
  local f="$1"
  local label="${2:-file}"
  if [[ ! -f "$f" ]]; then
    log_error "$label '$f' does not exist or is not a regular file."
    exit 1
  fi
}

###############################################################################
### 2. Define help
###############################################################################

show_usage() {
  cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Required:
  --reads1 CHAR              Input R1 metagenome data (fastq/fa)
  --reads2 CHAR              Input R2 metagenome data (fastq/fa)
  --sample_name CHAR         Sample name used to name the files

Optional:
  --assem_dir CHAR           Directory with previously computed assemblies
                             (format: dirname/SAMPLE_NAME/SAMPLE_NAME.contigs.fa)
  --assem_preset CHAR        MEGAHIT preset to generate assembly
                             (default: meta-sensitive)
  --nslots NUM               Number of threads used (default: 12)
  --min_contig_length NUM    Minimum length of contigs to keep (default: 250)
  --output_dir CHAR          Output directory (default: mg-clust_output-1)
  --overwrite t|f            Overwrite previous folder if present (default: f)
  --help                     Print this help and exit

Examples:
  $(basename "$0") \\
    --reads1 sample_R1.fastq.gz --reads2 sample_R2.fastq.gz \\
    --sample_name Sample1 --nslots 16 --output_dir Sample1_map

EOF
}

###############################################################################
### 3. Parse input parameters (using getopt)
###############################################################################

# Defaults (can be overridden by CLI)
ASSEM_DIR=""
ASSEM_PRESET="meta-sensitive"
NSLOTS=12
MIN_CONTIG_LENGTH=250
OUTPUT_DIR="mg-clust_output-1"
OVERWRITE="f"
R1=""
R2=""
SAMPLE_NAME=""

# Ensure getopt is available
check_cmd getopt

PARSED_ARGS=$(
  getopt -o '' \
    --long help,assem_dir:,assem_preset:,nslots:,min_contig_length:,output_dir:,overwrite:,reads1:,reads2:,sample_name: \
    -n "$(basename "$0")" -- "$@" \
) || {
  log_error "Failed to parse arguments."
  show_usage
  exit 1
}

eval set -- "${PARSED_ARGS}"

while true; do
  case "$1" in
    --help)
      show_usage
      exit 0
      ;;
    --assem_dir)
      ASSEM_DIR="$2"; shift 2 ;;
    --assem_preset)
      ASSEM_PRESET="$2"; shift 2 ;;
    --nslots)
      NSLOTS="$2"; shift 2 ;;
    --min_contig_length)
      MIN_CONTIG_LENGTH="$2"; shift 2 ;;
    --output_dir)
      OUTPUT_DIR="$2"; shift 2 ;;
    --overwrite)
      OVERWRITE="$2"; shift 2 ;;
    --reads1)
      R1="$2"; shift 2 ;;
    --reads2)
      R2="$2"; shift 2 ;;
    --sample_name)
      SAMPLE_NAME="$2"; shift 2 ;;
    --)
      shift
      break
      ;;
    *)
      log_error "Unexpected argument: $1"
      show_usage
      exit 1
      ;;
  esac
done

###############################################################################
### 3.1 Validate parameters
###############################################################################

# Mandatory parameters
for v in R1 R2 SAMPLE_NAME; do
  if [[ -z "${!v}" ]]; then
    log_error "Missing mandatory parameter: --$(echo "$v" | tr '[:upper:]' '[:lower:]')"
    show_usage
    exit 1
  fi
done

# Simple validation for overwrite flag
if [[ "${OVERWRITE}" != "t" && "${OVERWRITE}" != "f" ]]; then
  log_error "--overwrite must be 't' or 'f' (got '${OVERWRITE}')"
  exit 1
fi

# Numeric validations
if ! [[ "${NSLOTS}" =~ ^[0-9]+$ ]]; then
  log_error "--nslots must be an integer (got '${NSLOTS}')"
  exit 1
fi

if ! [[ "${MIN_CONTIG_LENGTH}" =~ ^[0-9]+$ ]]; then
  log_error "--min_contig_length must be an integer (got '${MIN_CONTIG_LENGTH}')"
  exit 1
fi

###############################################################################
### 3.2 Check tools and conf variables
###############################################################################

# Ensure paths/commands from conf.sh are defined
: "${megahit:?megahit path not set in conf.sh}"
: "${bwa:?bwa path not set in conf.sh}"
: "${samtools:?samtools path not set in conf.sh}"
: "${picard:?picard jar not set in conf.sh}"

check_cmd "${megahit}"
check_cmd "${bwa}"
check_cmd "${samtools}"
check_cmd java

check_file "${picard}" "Picard jar"

###############################################################################
### 4. Check mandatory files
###############################################################################

check_file "${R1}" "read1"
check_file "${R2}" "read2"

###############################################################################
### 5. Create output folder
###############################################################################

if [[ -d "${OUTPUT_DIR}" ]]; then
  if [[ "${OVERWRITE}" == "t" ]]; then
    log_warn "Output directory '${OUTPUT_DIR}' exists and will be overwritten."
    rm -rf "${OUTPUT_DIR}" || {
      log_error "Failed to remove existing output directory '${OUTPUT_DIR}'."
      exit 1
    }
  else
    log_warn "Output directory '${OUTPUT_DIR}' already exists; use --overwrite t to overwrite."
    exit 0
  fi
fi

# Create base output directory if assembly is not run by megahit (when ASSEM_DIR is given)
if [[ -n "${ASSEM_DIR}" ]]; then
  mkdir -p "${OUTPUT_DIR}" || {
    log_error "Failed to create output directory '${OUTPUT_DIR}'."
    exit 1
  }
fi

###############################################################################
### 6. De novo assembly (if ASSEM_DIR not provided)
###############################################################################

ASSEMBLY_FILE=""

if [[ -z "${ASSEM_DIR}" ]]; then
  log "Running MEGAHIT assembly..."

  "${megahit}" \
    --num-cpu-threads "${NSLOTS}" \
    -1 "${R1}" \
    -2 "${R2}" \
    --presets "${ASSEM_PRESET}" \
    --min-contig-len "${MIN_CONTIG_LENGTH}" \
    --out-prefix "${SAMPLE_NAME}" \
    --out-dir "${OUTPUT_DIR}"

  log "MEGAHIT completed."

  ASSEMBLY_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}.contigs.fa"
else
  ASSEMBLY_FILE="${ASSEM_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.contigs.fa"
fi

# Check assembly file exists
check_file "${ASSEMBLY_FILE}" "assembly file"

###############################################################################
### 7. Map short reads
###############################################################################

log "Checking number of assembled contigs in '${ASSEMBLY_FILE}'..."
ASSEMBLY_FILE_NSEQ=$(grep -c ">" "${ASSEMBLY_FILE}" || echo 0)

if (( ASSEMBLY_FILE_NSEQ < 5 )); then
  log_warn "Not enough assembled sequences to continue (${ASSEMBLY_FILE_NSEQ} < 5). Exiting gracefully."
  # Distinguish from success but not a hard error
  exit 2
fi

log "Found ${ASSEMBLY_FILE_NSEQ} contigs. Proceeding with mapping."

# Index assembly
log "Indexing assembly with BWA..."
"${bwa}" index "${ASSEMBLY_FILE}" || {
  log_error "bwa index failed."
  exit 1
}

# Map reads
log "Mapping reads with BWA-MEM..."
"${bwa}" mem -M -t "${NSLOTS}" "${ASSEMBLY_FILE}" "${R1}" "${R2}" > \
  "${OUTPUT_DIR}/${SAMPLE_NAME}.sam" || {
  log_error "bwa mem failed."
  exit 1
}

# Convert to BAM and filter high-quality primary alignments
log "Converting SAM to BAM and filtering (q>=10, primary alignments)..."
"${samtools}" view \
  -@ "${NSLOTS}" \
  -q 10 \
  -F 260 \
  -b "${OUTPUT_DIR}/${SAMPLE_NAME}.sam" > "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" || {
  log_error "samtools view (SAM→BAM) failed."
  exit 1
}

# Sort BAM
log "Sorting BAM..."
"${samtools}" sort \
  -@ "${NSLOTS}" \
  "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" > "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" || {
  log_error "samtools sort failed."
  exit 1
}

# Index sorted BAM
log "Indexing sorted BAM..."
"${samtools}" index -@ "${NSLOTS}" "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" || {
  log_error "samtools index on sorted BAM failed."
  exit 1
}

# Remove duplicates with Picard
log "Marking and removing duplicates with Picard..."
mkdir -p "${OUTPUT_DIR}/tmp" || {
  log_error "Failed to create temporary directory '${OUTPUT_DIR}/tmp'."
  exit 1
}

java -jar "${picard}" MarkDuplicates \
  INPUT="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" \
  OUTPUT="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.bam" \
  METRICS_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.metrics.txt" \
  REMOVE_DUPLICATES=TRUE \
  ASSUME_SORTED=TRUE \
  MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=900 \
  TMP_DIR="${OUTPUT_DIR}/tmp" || {
  log_error "Picard MarkDuplicates failed."
  exit 1
}

# Index final markdup BAM
log "Indexing final duplicate-removed BAM..."
"${samtools}" index -@ "${NSLOTS}" "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.bam" || {
  log_error "samtools index on markdup BAM failed."
  exit 1
}

###############################################################################
### 8. Clean
###############################################################################

log "Removing intermediate mapping files..."
rm -f "${OUTPUT_DIR}/${SAMPLE_NAME}.sam" \
      "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" || {
  log_error "Failed to remove intermediate SAM/BAM files."
  exit 1
}

rm -rf "${OUTPUT_DIR}/tmp" || {
  log_error "Failed to remove temporary directory '${OUTPUT_DIR}/tmp'."
  exit 1
}

if [[ -z "${ASSEM_DIR}" ]]; then
  # Remove MEGAHIT intermediate contigs but keep main assembly
  if [[ -d "${OUTPUT_DIR}/intermediate_contigs" ]]; then
    log "Removing MEGAHIT intermediate contigs..."
    rm -rf "${OUTPUT_DIR}/intermediate_contigs" || {
      log_warn "Failed to remove 'intermediate_contigs' directory (non-fatal)."
    }
  fi
fi

###############################################################################
### 9. Exit
###############################################################################

log "${GREEN}mg-clust_module-1.bash exited successfully.${NC}"
exit 0
