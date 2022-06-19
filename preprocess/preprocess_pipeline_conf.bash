# dirs
BIN="/home/bioinf/bin"
LOCAL="/home/epereira/.local/"
PREPROCESS_DIR="${HOME}/workspace/repositories/tools/metagenomic_pipelines/preprocess"

# tools
bbmap_version="38.79"
bbduk="${BIN}/bbmap/bbmap-${bbmap_version}/bbduk.sh"
bbmerge="${BIN}/bbmap/bbmap-${bbmap_version}/bbmerge.sh"
seqtk="${BIN}/seqtk/seqtk"
pear_version="0.9.8"
pear="${BIN}/pear/pear-${pear_version}/bin/pear"
cutadapt="${LOCAL}/bin/cutadapt"
vsearch="/usr/bin/vsearch"
bzip2="/bin/bzip2"
gunzip="/bin/gunzip"
fq2fa="${PREPROCESS_DIR}/resources/fq2fa.sh"
plots="${PREPROCESS_DIR}/resources/plots.R"
quality_check_plots="${PREPROCESS_DIR}/resources/quality_check_plots.R"
pigz="/usr/bin/pigz"

# files
ADAPTERS="${BIN}/bbmap/bbmap-${bbmap_version}/resources/adapters.fa"

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

