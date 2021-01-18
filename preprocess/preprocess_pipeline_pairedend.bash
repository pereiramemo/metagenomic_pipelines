###############################################################################
### 1. Set env
###############################################################################

set -o pipefail 
source "${HOME}/workspace/repositories/metagenomic_pipelines/preprocess/\
preprocess_pipeline_conf.bash"

###############################################################################
# 2. Define help
###############################################################################

show_usage(){
  cat <<EOF
Usage: ./preprocess_pipeline.bash <options>
--help                          print this help
--clean t|f                     remove all intermediate files
--min_qual NUM                  minimum quality score to trim reads (default 20)
--nslots NUM                    number of threads used (default 12)
--output_dir CHAR               directory to output generated data (i.e., preprocessed data, plots, tables)
--overwrite t|f                 overwrite previous directory (default f)
--reads CHAR                    input R1 reads
--reads2 CHAR                   input R2 reads
--sample_name CHAR              sample name (default metagenomex)
--subsample t|f                 subsample metagenome to 10K to test execution (default f)
--trim_adapters t|f             check for adapters and trim (default f)
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
  --clean)
  if [[ -n "${2}" ]]; then
    CLEAN="${2}"
    shift
  fi
  ;;
  --clean=?*)
  CLEAN="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --clean=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;
#############
  --min_qual)
  if [[ -n "${2}" ]]; then
    MIN_QUAL="${2}"
    shift
  fi
  ;;
  --min_qual=?*)
  MIN_QUAL="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --min_qual=) # Handle the empty case
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
  --reads)
  if [[ -n "${2}" ]]; then
    R1="${2}"
    shift
  fi
  ;;
  --reads=?*)
  R1="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --reads=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;  
#############
  --reads2)
  if [[ -n "${2}" ]]; then
    R2="${2}"
    shift
  fi
  ;;
  --reads2=?*)
  R2="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --reads2=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;    
#############
  --sample_name)
  if [[ -n "${2}" ]]; then
    SAMPLE_NAME="${2}"
    shift
  fi
  ;;
  --sample_name=?*)
  SAMPLE_NAME="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --sample_name=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;  
#############
  --subsample)
  if [[ -n "${2}" ]]; then
    SUBSAMPLE="${2}"
    shift
  fi
  ;;
  --subsample=?*)
  SUBSAMPLE="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --subsample=) # Handle the empty case
  printf 'Using default environment.\n' >&2
  ;;    
#############
  --trim_adapters)
  if [[ -n "${2}" ]]; then
    TRIM_ADAPTERS="${2}"
    shift
  fi
  ;;
  --trim_adapters=?*)
  TRIM_ADAPTERS="${1#*=}" # Delete everything up to "=" and assign the remainder.
  ;;
  --trim_adapters=) # Handle the empty case
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
### 4. Check mandatory parameters
###############################################################################

if [[ -z "${OUTPUT_DIR}" ]]; then
  echo "missing output dir"
  exit 1
fi

if [[ ! -a "${R1}" || ! -a "${R2}" ]]; then
  echo "missing input reads"
  exit 1
fi

###############################################################################
### 5. Define defaults
###############################################################################

if [[ -z "${CLEAN}" ]]; then
  CLEAN="f"
fi

if [[ -z "${MIN_QUAL}" ]]; then
  MIN_QUAL=20
fi  

if [[ -z "${NSLOTS}" ]]; then
  NSLOTS=12
fi  

if [[ -z "${OVERWRITE}" ]]; then
  OVERWRITE="f"
fi

if [[ -z "${SAMPLE_NAME}" ]]; then
  SAMPLE_NAME="metagenomex"
fi  

if [[ -z "${SUBSAMPLE}" ]]; then
  SUBSAMPLE="f"
fi

if [[ -z "${TRIM_ADAPTERS}" ]]; then
  TRIM_ADAPTERS="f"
fi

###############################################################################
### 6. Create output dir
###############################################################################

### test parameters
# DATA_DIR="/home/epereira/workspace/dev/indicators_contaminants_2018/data/"
# PREPRO_DIR="/home/epereira/workspace/dev/indicators_contaminants_2018/data/metagenomics/prepro"
# R1="${DATA_DIR}/metagenomics/AdapterClipped/Sample_1/1_R1_clipped.fastq.bz2"
# R2="${DATA_DIR}/metagenomics/AdapterClipped/Sample_1/1_R2_clipped.fastq.bz2"
# NSLOTS=12
# OUTPUT_DIR="${PREPRO_DIR}/sample_1_paired_end"
# SAMPLE_NAME="sample_1"
# OVERWRITE="f"
#####################

if [[ "${OVERWRITE}" == "t" ]]; then

  if [[ -d "${OUTPUT_DIR}" ]]; then
    rm -r "${OUTPUT_DIR}"
  fi
   mkdir -p "${OUTPUT_DIR}"
  
fi

if [[ "${OVERWRITE}" == "f" ]]; then

  if [[ -d "${OUTPUT_DIR}" ]]; then
    echo "${OUTPUT_DIR} already exists"
    exit 1
  else
    mkdir -p "${OUTPUT_DIR}"
  fi
  
fi  

###############################################################################
### 7. Check for compressed data
###############################################################################

UNCOMPRESSED="f"
TEST_FILE_BZIP2=$(file "${R1}"| egrep "bzip2")
TEST_FILE_GZIP=$(file "${R1}"| egrep "gzip")

if [[ -n "${TEST_FILE_BZIP2}" ]]; then

  echo "uncompressing ..."
  R1_UNCOMPRESSED="${OUTPUT_DIR}/${SAMPLE_NAME}_R1-00.fastq"
  R2_UNCOMPRESSED="${OUTPUT_DIR}/${SAMPLE_NAME}_R2-00.fastq"
  
  "${bzip2}" \
  --decompress \
  --keep \
  --stdout "${R1}" > "${R1_UNCOMPRESSED}"
  
  if [[ $? -ne 0 ]]; then
    echo "bzip2 R1 failed"
    exit 1
  fi 
  
  "${bzip2}" \
  --decompress \
  --keep \
  --stdout "${R2}" > "${R2_UNCOMPRESSED}"

  if [[ $? -ne 0 ]]; then
    echo "bzip2 R2 failed"
    exit 1
  fi 
  
  R1="${R1_UNCOMPRESSED}"
  R2="${R2_UNCOMPRESSED}"
  UNCOMPRESSED="t"
  
fi

if [[ -n "${TEST_FILE_GZIP}" ]]; then

  echo "uncompressing ..."
  R1_UNCOMPRESSED="${OUTPUT_DIR}/${SAMPLE_NAME}_R1-00.fastq"
  R2_UNCOMPRESSED="${OUTPUT_DIR}/${SAMPLE_NAME}_R2-00.fastq"
  
  "${gunzip}" \
  --keep \
  --stdout "${R1}" > "${R1_UNCOMPRESSED}"
  
  if [[ $? -ne 0 ]]; then
    echo "gunzip R1 failed"
    exit 1
  fi 
  
  "${gunzip}" \
  --keep \
  --stdout "${R2}" > "${R2_UNCOMPRESSED}"

  if [[ $? -ne 0 ]]; then
    echo "gunzip R2 failed"
    exit 1
  fi 
  
  R1="${R1_UNCOMPRESSED}"
  R2="${R2_UNCOMPRESSED}"
  UNCOMPRESSED="t"
  
fi

if [[ "${UNCOMPRESSED}" == "f" ]]; then

  ln -s "${R1}" "${OUTPUT_DIR}/${SAMPLE_NAME}_R1-00.fastq"
  if [[ $? -ne 0 ]]; then
    echo "ln -s R1 failed"
    exit 1
  fi 
  
  ln -s "${R2}" "${OUTPUT_DIR}/${SAMPLE_NAME}_R2-00.fastq"
  if [[ $? -ne 0 ]]; then
    echo "ln -s R2 failed"
    exit 1
  fi
  
  R1="${OUTPUT_DIR}/${SAMPLE_NAME}_R1-00.fastq"
  R2="${OUTPUT_DIR}/${SAMPLE_NAME}_R2-00.fastq"
  
fi  
  
##############################################################################
### 8. Subsample data
###############################################################################

if [[ "${SUBSAMPLE}" == "t" ]]; then

  echo "subsampling ..."
  
  R1_REDU="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_redu-00.fastq"
  R2_REDU="${OUTPUT_DIR}/${SAMPLE_NAME}_R2_redu-00.fastq"
   
  "${seqtk}" sample -s123 "${R1}" 10000 > "${R1_REDU}" 
  if [[ $? -ne 0 ]]; then
    echo "seqtk subsampling R1 failed"
    exit 1
  fi
  
  "${seqtk}" sample -s123 "${R2}" 10000 > "${R2_REDU}"
  if [[ $? -ne 0 ]]; then
    echo "seqtk subsampling R2 failed"
    exit 1
  fi
  
  R1="${R1_REDU}"
  R2="${R2_REDU}"
  
fi  
  
###############################################################################
### 9. Adapter trimming
###############################################################################

if [[ "${TRIM_ADAPTERS}" == "t" ]]; then

  echo "trimming adapters ..."

  R1_AT="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_at-01.fastq"
  R2_AT="${OUTPUT_DIR}/${SAMPLE_NAME}_R2_at-01.fastq"

  "${bbduk}" \
  in="${R1}" \
  in2="${R2}" \
  out="${R1_AT}" \
  out2="${R2_AT}" \
  threads="${NSLOTS}" \
  ktrim=r \
  k=23 \
  mink=11 \
  hdist=1 \
  tpe \
  tbo \
  ref="${ADAPTERS}"

  if [[ $? -ne 0 ]]; then
    echo "bbduk quality trimming failed"
    exit 1
  fi  

  R1_INPUT="${R1}"
  R2_INPUT="${R2}"
  R1="${R1_AT}"
  R2="${R2_AT}"
  
fi

###############################################################################
### 10. Quality check paired-end reads
###############################################################################

R1_QC="${OUTPUT_DIR}/${SAMPLE_NAME}_R1_qc-02.fastq"
R2_QC="${OUTPUT_DIR}/${SAMPLE_NAME}_R2_qc-02.fastq"

"${bbduk}" \
in="${R1}" \
in2="${R2}" \
out="${R1_QC}" \
out2="${R2_QC}" \
minlength=50 \
threads="${NSLOTS}" \
qtrim=rl \
trimq="${MIN_QUAL}"

if [[ $? -ne 0 ]]; then
  echo "bbduk quality trimming paired-end reads failed"
  exit 1
fi  

###############################################################################
### 11. Compute stats
###############################################################################

echo "computing statistics ..."

if [[ -a "${R1_AT}" ]]; then

  FILES=$(ls "${R1_INPUT}" "${R2_INPUT}" "${R1_AT}" "${R2_AT}" "${R1_QC}" "${R2_QC}")

else

  FILES=$(ls "${R1_INPUT}" "${R2_INPUT}" "${R1_QC}" "${R2_QC}")

fi

echo "${FILES}" | \
while read FILE; do

  if [[ -a "${FILE}" ]]; then
  
    N=$(count_fastq "${FILE}")
  
    if [[ $? -ne 0 ]]; then
      echo "count_fastq failed for ${FILE}"
      exit 1
    fi
  
    if [[ -s "${FILE}" ]]; then
    
      L=$(compute_mean_length "${FILE}")
      
      if [[ $? -ne 0 ]]; then
        echo "compute_mean_length failed for ${FILE}"
        exit 1
      fi
      
    else 
    
      L=0
      
    fi
    
    echo -e "${SAMPLE_NAME}\t$(basename ${FILE})\tnum_seq\t${N}" >> \
    "${OUTPUT_DIR}/stats_tmp.tsv"
    echo -e "${SAMPLE_NAME}\t$(basename ${FILE})\tmean_length\t${L}" >> \
    "${OUTPUT_DIR}/stats_tmp.tsv"
  
  fi 
  
done

awk -v "OFS=\t" '{
    n = gensub(".*-([0-9]+).fastq","\\1","g",$2); 
    print $0,n}' \
    "${OUTPUT_DIR}/stats_tmp.tsv" > "${OUTPUT_DIR}/stats.tsv"

if [[ $? -ne 0 ]]; then
  echo "awk formatting stats file failed"
  exit 1
fi

rm "${OUTPUT_DIR}/stats_tmp.tsv"

###############################################################################
### 12. plot stats
###############################################################################

echo "creating plot ..."

Rscript --vanilla \
"${PREPROCESS_DIR}/plots.R" \
"${OUTPUT_DIR}/stats.tsv" \
"${OUTPUT_DIR}/stats_plots.png"

if [[ $? -ne 0 ]]; then
  echo "R plotting failed"
  exit 1
fi

###############################################################################
### 17. Clean
###############################################################################

if [[ "${CLEAN}" == "t" ]]; then

  if [[ "${UNCOMPRESSED}" == "t" ]]; then
 
    rm "${R1_UNCOMPRESSED}" "${R2_UNCOMPRESSED}"
    
    if [[ $? -ne 0 ]]; then
      echo "rm R1 and R2 files failed"
      exit 1
    fi
  fi  
    
    
  if [[ "${SUBSAMPLE}" == "t" ]]; then
   
    rm "${R1_REDU}" "${R2_REDU}"
     
    if [[ $? -ne 0 ]]; then
      echo "rm R1 and R2 redu files failed"
      exit 1
    fi
  fi

  
  if [[ "${TRIM_ADAPTERS}" == "t" ]]; then
  
    rm "${R1_AT}" "${R2_AT}" 
    
    if [[ $? -ne 0 ]]; then
      echo "rm R1_AT and R2_AT files failed"
      exit 1
    fi
  fi
fi
