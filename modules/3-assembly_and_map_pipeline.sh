#!/usr/bin/env bash
###############################################################################
### 3-assembly_and_map_pipeline.sh
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

# Find contigs file in assembly directory
# Searches for contigs files with various extensions in multiple locations
# Supports both compressed (.gz) and uncompressed files
find_contigs_file() {
  local assem_dir="$1"
  local sample_name="$2"
  local extensions=("contigs.fa" "contigs.fasta" "contigs.fna" "fa" "fasta" "fna" \
                    "contigs.fa.gz" "contigs.fasta.gz" "contigs.fna.gz" "fa.gz" "fasta.gz" "fna.gz")
  local search_paths=(
    "${assem_dir}/${sample_name}"
    "${assem_dir}"
  )

  # First, look for files with sample name prefix
  for dir in "${search_paths[@]}"; do
    for ext in "${extensions[@]}"; do
      local candidate="${dir}/${sample_name}.${ext}"
      if [[ -f "${candidate}" ]]; then
        echo "${candidate}"
        return 0
      fi
    done
  done

  # If not found, look for any .contigs.{fa,fasta,fna}[.gz] file in the directories
  # Enable nullglob so unmatched patterns expand to nothing
  shopt -s nullglob
  for dir in "${search_paths[@]}"; do
    for ext in fa fasta fna; do
      # Check uncompressed files
      for candidate in "${dir}"/*.contigs."${ext}"; do
        if [[ -f "${candidate}" ]]; then
          log_warn "Found contigs file '${candidate}' but it doesn't match expected pattern '${sample_name}.contigs.*'"
          log_warn "Using it anyway. Consider renaming or using --contigs to specify explicitly."
          shopt -u nullglob  # Disable nullglob before returning
          echo "${candidate}"
          return 0
        fi
      done
      # Check compressed files
      for candidate in "${dir}"/*.contigs."${ext}".gz; do
        if [[ -f "${candidate}" ]]; then
          log_warn "Found contigs file '${candidate}' but it doesn't match expected pattern '${sample_name}.contigs.*'"
          log_warn "Using it anyway. Consider renaming or using --contigs to specify explicitly."
          shopt -u nullglob  # Disable nullglob before returning
          echo "${candidate}"
          return 0
        fi
      done
    done
  done
  shopt -u nullglob  # Disable nullglob before returning

  # Not found
  return 1
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
  --contigs CHAR             Path to pre-assembled contigs file (FASTA format)
                             Supports both compressed (.gz) and uncompressed files
                             Takes precedence over --assem_dir
  --assem_dir CHAR           Directory with previously computed assemblies
                             Will search for: SAMPLE_NAME.contigs.{fa,fasta,fna}[.gz]
                             in ASSEM_DIR/ or ASSEM_DIR/SAMPLE_NAME/
                             Supports both compressed (.gz) and uncompressed files
  --assem_preset CHAR        MEGAHIT preset to generate assembly
                             (default: meta-sensitive)
  --nslots NUM               Number of threads used (default: 12)
  --min_contig_length NUM    Minimum length of contigs to keep (default: 250)
  --output_dir CHAR          Output directory (default: mg-clust_output-1)
  --overwrite t|f            Overwrite previous folder if present (default: f)
  --remove_duplicates t|f    Remove PCR duplicates with Picard (default: f)
  --help                     Print this help and exit

Examples:
  # Run de novo assembly with MEGAHIT:
  $(basename "$0") \\
    --reads1 sample_R1.fastq.gz --reads2 sample_R2.fastq.gz \\
    --sample_name Sample1 --nslots 16 --output_dir Sample1_map

  # Use pre-assembled contigs (direct path):
  $(basename "$0") \\
    --reads1 sample_R1.fastq.gz --reads2 sample_R2.fastq.gz \\
    --sample_name Sample1 --contigs /path/to/Sample1.contigs.fa \\
    --output_dir Sample1_map

  # Use compressed pre-assembled contigs:
  $(basename "$0") \\
    --reads1 sample_R1.fastq.gz --reads2 sample_R2.fastq.gz \\
    --sample_name Sample1 --contigs /path/to/Sample1.contigs.fa.gz \\
    --output_dir Sample1_map

  # Use pre-assembled contigs (search in directory):
  $(basename "$0") \\
    --reads1 sample_R1.fastq.gz --reads2 sample_R2.fastq.gz \\
    --sample_name Sample1 --assem_dir /path/to/assemblies \\
    --output_dir Sample1_map

EOF
}

###############################################################################
### 3. Parse input parameters (using getopt)
###############################################################################

# Defaults (can be overridden by CLI)
CONTIGS=""
ASSEM_DIR=""
ASSEM_PRESET="meta-sensitive"
NSLOTS=12
MIN_CONTIG_LENGTH=250
OUTPUT_DIR="mg-clust_output-1"
OVERWRITE="f"
REMOVE_DUPLICATES="f"
R1=""
R2=""
SAMPLE_NAME=""

# Ensure getopt is available
check_cmd getopt

ARGS=$(
  getopt -o '' \
    --long help,contigs:,assem_dir:,assem_preset:,nslots:,min_contig_length:,output_dir:,overwrite:,remove_duplicates:,reads1:,reads2:,sample_name: \
    -n "$(basename "$0")" -- "$@" \
) || {
  log_error "Failed to parse arguments."
  show_usage
  exit 1
}

eval set -- "${ARGS}"

while true; do
  case "$1" in
    --help)
      show_usage
      exit 0
      ;;
    --contigs)
      CONTIGS="$2"; shift 2 ;;
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
    --remove_duplicates)
      REMOVE_DUPLICATES="$2"; shift 2 ;;
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
### 4 Validate parameters
###############################################################################

# Mandatory parameters
for v in R1 R2 SAMPLE_NAME; do
  if [[ -z "${!v}" ]]; then
    log_error "Missing mandatory parameter: --$(echo "$v" | tr '[:upper:]' '[:lower:]')"
    show_usage
    exit 1
  fi
done

# Simple validation for boolean flags
if [[ "${OVERWRITE}" != "t" && "${OVERWRITE}" != "f" ]]; then
  log_error "--overwrite must be 't' or 'f' (got '${OVERWRITE}')"
  exit 1
fi

if [[ "${REMOVE_DUPLICATES}" != "t" && "${REMOVE_DUPLICATES}" != "f" ]]; then
  log_error "--remove_duplicates must be 't' or 'f' (got '${REMOVE_DUPLICATES}')"
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
### 5. Check tools and conf variables
###############################################################################

check_dependencies

# Check picard - can be either a command or jar file
if ! command -v "${picard}" >/dev/null 2>&1 && [[ ! -f "${picard}" ]]; then
    log_error "Picard not found in PATH or as jar file: ${picard}"
    exit 1
fi

###############################################################################
### 6. Check mandatory files
###############################################################################

check_file "${R1}" "read1"
check_file "${R2}" "read2"

###############################################################################
### 7. Create output folder
###############################################################################

if [[ -d "${OUTPUT_DIR}" ]]; then
  if [[ "${OVERWRITE}" == "t" ]]; then
    log_warn "Output directory '${OUTPUT_DIR}' exists and will be overwritten."
    rm -rf "${OUTPUT_DIR}" || {
      log_error "Failed to remove existing output directory '${OUTPUT_DIR}'."
      exit 1
    }
  else
    log_error "Output directory '${OUTPUT_DIR}' already exists; use --overwrite t to overwrite."
    exit 1
  fi
fi

# Create base output directory if assembly is not run by megahit (when CONTIGS or ASSEM_DIR is given)
if [[ -n "${CONTIGS}" || -n "${ASSEM_DIR}" ]]; then
  mkdir -p "${OUTPUT_DIR}" || {
    log_error "Failed to create output directory '${OUTPUT_DIR}'."
    exit 1
  }
fi

# Set up trap to clean up tmp directory on exit (success or failure)
trap 'rm -rf "${OUTPUT_DIR}/tmp" 2>/dev/null' EXIT

###############################################################################
### 8. De novo assembly or identify external contigs
###############################################################################

ASSEMBLY_FILE=""
ASSEMBLY_FILE_DECOMPRESSED=""  # Track if we created a decompressed temp file

# Priority 1: Use explicitly provided contigs file
if [[ -n "${CONTIGS}" ]]; then
  log "Using externally provided contigs file: ${CONTIGS}"
  ASSEMBLY_FILE="${CONTIGS}"
  check_file "${ASSEMBLY_FILE}" "contigs file"

# Priority 2: Search for contigs in provided assembly directory
elif [[ -n "${ASSEM_DIR}" ]]; then
  log "Searching for contigs file in assembly directory: ${ASSEM_DIR}"

  # Validate assembly directory exists
  if [[ ! -d "${ASSEM_DIR}" ]]; then
    log_error "Assembly directory '${ASSEM_DIR}' does not exist or is not a directory."
    exit 1
  fi

  # Try to find contigs file
  ASSEMBLY_FILE=$(find_contigs_file "${ASSEM_DIR}" "${SAMPLE_NAME}")

  if [[ -z "${ASSEMBLY_FILE}" ]]; then
    log_error "Could not find contigs file for sample '${SAMPLE_NAME}' in '${ASSEM_DIR}'."
    log_error "Searched in:"
    log_error "  - ${ASSEM_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}.{contigs.fa,contigs.fasta,contigs.fna,fa,fasta,fna}[.gz]"
    log_error "  - ${ASSEM_DIR}/${SAMPLE_NAME}.{contigs.fa,contigs.fasta,contigs.fna,fa,fasta,fna}[.gz]"
    log_error "  - ${ASSEM_DIR}/${SAMPLE_NAME}/*.contigs.{fa,fasta,fna}[.gz]"
    log_error "  - ${ASSEM_DIR}/*.contigs.{fa,fasta,fna}[.gz]"
    log_error ""
    log_error "Use --contigs to specify the exact path to your contigs file."
    exit 1
  fi

  log "Found contigs file: ${ASSEMBLY_FILE}"

# Priority 3: Run de novo assembly with MEGAHIT
else
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
  check_file "${ASSEMBLY_FILE}" "assembly file"
fi

# Handle compressed assembly files
if [[ "${ASSEMBLY_FILE}" == *.gz ]]; then
  log "Detected gzipped assembly file. Decompressing..."

  # Create temporary decompressed file in output directory
  ASSEMBLY_FILE_DECOMPRESSED="${OUTPUT_DIR}/$(basename "${ASSEMBLY_FILE}" .gz)"

  "${gunzip}" -c "${ASSEMBLY_FILE}" > "${ASSEMBLY_FILE_DECOMPRESSED}" || {
    log_error "Failed to decompress assembly file."
    exit 1
  }

  log "Decompressed to: ${ASSEMBLY_FILE_DECOMPRESSED}"

  # Update ASSEMBLY_FILE to point to decompressed version
  ASSEMBLY_FILE="${ASSEMBLY_FILE_DECOMPRESSED}"
  check_file "${ASSEMBLY_FILE}" "decompressed assembly file"
fi

###############################################################################
### 9. Map short reads
###############################################################################

log "Checking number of assembled contigs in '${ASSEMBLY_FILE}'..."
ASSEMBLY_FILE_NSEQ=$(grep -c ">" "${ASSEMBLY_FILE}" || echo 0)

if (( ASSEMBLY_FILE_NSEQ < 5 )); then 
  log_warn "Not enough assembled sequences to continue (${ASSEMBLY_FILE_NSEQ} < 5). Exiting gracefully."
  # Distinguish from success but not a hard error
  exit 0
fi

log "Found ${ASSEMBLY_FILE_NSEQ} contigs. Proceeding with mapping."

# Index assembly
log "Indexing assembly with BWA..."
"${bwa}" index "${ASSEMBLY_FILE}" || {
  log_error "bwa index failed."
  exit 1
}

# Map reads and convert directly to BAM (avoid huge SAM files)
log "Mapping reads with BWA-MEM and converting to BAM (q>=10, primary alignments)..."
"${bwa}" mem -M -t "${NSLOTS}" "${ASSEMBLY_FILE}" "${R1}" "${R2}" | \
  "${samtools}" view -@ "${NSLOTS}" -q 10 -F 260 -b > "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" || {
  log_error "bwa mem or samtools view failed."
  exit 1
}

# Sort BAM
log "Sorting BAM..."
"${samtools}" sort \
  -@ "${NSLOTS}" \
  -o "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" \
  "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" || {
  log_error "samtools sort failed."
  exit 1
}

# Index sorted BAM
log "Indexing sorted BAM..."
"${samtools}" index -@ "${NSLOTS}" "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" || {
  log_error "samtools index on sorted BAM failed."
  exit 1
}

# Remove duplicates with Picard (optional)
if [[ "${REMOVE_DUPLICATES}" == "t" ]]; then
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

  FINAL_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted_markdup.bam"
else
  log "Skipping duplicate removal (--remove_duplicates=f)"
  FINAL_BAM="${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam"
fi

###############################################################################
### 10. Clean
###############################################################################

log "Removing intermediate mapping files..."
# Always remove the initial unsorted BAM
rm -f "${OUTPUT_DIR}/${SAMPLE_NAME}.bam" || {
  log_error "Failed to remove intermediate BAM file."
  exit 1
}

# If duplicates were removed, also remove the sorted BAM (it's intermediate)
if [[ "${REMOVE_DUPLICATES}" == "t" ]]; then
  rm -f "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam" \
        "${OUTPUT_DIR}/${SAMPLE_NAME}_sorted.bam.bai" || {
    log_warn "Failed to remove intermediate sorted BAM files (non-fatal)."
  }
fi

# Remove BWA index files (always cleanup since they were created by this script)
log "Removing BWA index files..."
rm -f "${ASSEMBLY_FILE}".{amb,ann,bwt,pac,sa} || {
  log_warn "Failed to remove some BWA index files (non-fatal)."
}

# Remove temporary decompressed assembly file if created
if [[ -n "${ASSEMBLY_FILE_DECOMPRESSED}" && -f "${ASSEMBLY_FILE_DECOMPRESSED}" ]]; then
  log "Removing temporary decompressed assembly file..."
  rm -f "${ASSEMBLY_FILE_DECOMPRESSED}" || {
    log_warn "Failed to remove temporary decompressed assembly file (non-fatal)."
  }
fi

# Note: tmp directory is cleaned up automatically by trap on exit

if [[ -z "${CONTIGS}" && -z "${ASSEM_DIR}" ]]; then
  # Remove MEGAHIT intermediate contigs but keep main assembly
  if [[ -d "${OUTPUT_DIR}/intermediate_contigs" ]]; then
    log "Removing MEGAHIT intermediate contigs..."
    rm -rf "${OUTPUT_DIR}/intermediate_contigs" || {
      log_warn "Failed to remove 'intermediate_contigs' directory (non-fatal)."
    }
  fi
fi

###############################################################################
### 11. Exit
###############################################################################

log "${GREEN}assembly_and_map_pipeline.sh exited successfully.${NC}"
log "Final BAM file: ${FINAL_BAM}"
exit 0
