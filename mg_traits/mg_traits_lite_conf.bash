# dirs # general
BIOINF_DIR="/home/bioinf"
BIN="${BIOINF_DIR}/bin"
RESOURCES_DIR="${BIOINF_DIR}/resources"
MG_TRAITS_DIR="${HOME}/workspace/repositories/metagenomic_pipelines/mg_traits"
MG_TRAITS_RESOURCES="${MG_TRAITS_DIR}/resources_mg_traits"
LOCAL="/home/epereira/.local/"

# files
TRANSFACT_ACC="${MG_TRAITS_RESOURCES}/TF.txt"
RESFAM_HMM="${RESOURCES_DIR}/resfam/Resfams.hmm"
PFAM_MODEL_DIR="${RESOURCES_DIR}/pfam/uproc/model"
PFAM_DB="${RESOURCES_DIR}/pfam/uproc/pfam28_db"
BGC_MODEL_DIR="${RESOURCES_DIR}/bgc/model"
BGC_DB="${RESOURCES_DIR}/bgc/bgc13062014"

# tools
# bbduk
bbmap_version="38.79"
bbduk="${BIN}/bbmap/bbmap-${bbmap_version}/bbduk.sh"
bbmerge="${BIN}/bbmap/bbmap-${bbmap_version}/bbmerge.sh"
filterbyname="${BIN}/bbmap/bbmap-${bbmap_version}/filterbyname.sh"

# pear
pear_version="0.9.8"
pear="${BIN}/pear/pear-${pear_version}/bin/pear"

# fraggenescan
fraggenescanplusplus="${BIN}/FragGeneScan/FragGeneScanPlusPlus-master/FGSpp"
TRAIN="${BIN}/FragGeneScan/FragGeneScanPlusPlus-master/train"
fraggenescan_version="1.31"
fraggenescan="${BIN}/FragGeneScan/FragGeneScan-${fraggenescan_version}/run_FragGeneScan.pl"

#uproc
uproc_version="1.2.0"
uproc_prot="${BIN}/uproc/uproc-${uproc_version}/uproc-prot"
MODEL_DIR="${RESOURCES_DIR}/pfam/uproc/model"
DB="${RESOURCES_DIR}/pfam/uproc/pfam28_db"

# emboss
infoseq="/usr/bin/infoseq"
compseq="/usr/bin/compseq"
cusp="/usr/bin/cusp"

# ags and acn
ags="${BIN}/ags_n_acn/ags.sh"
acn="${BIN}/ags_n_acn/acn.sh"

# vsearch
vsearch_version="2.15.1"
vsearch="${BIN}/vsearch/vsearch-${vsearch_version}/bin/vsearch"

# hmmer
hmmer_version="3.3"
hmmsearch="${BIN}/hmmer/hmmer-${hmmer_version}/bin/hmmsearch"

# R
r_interpreter="/usr/bin/R"
r_script="/usr/bin/Rscript"
taxa_annot="${MG_TRAITS_RESOURCES}/taxa_annot_DADA2.R"

# other
seqtk="${BIN}/seqtk/seqtk"
cutadapt="${LOCAL}/bin/cutadapt"
seq_num_and_length_counter="${SCRIPTS}/seq_num_and_length_counter.bash"
bzip2="/bin/bzip2"
gunzip="/bin/gunzip"
fq2fa="/home/epereira/bin/fq2fa.sh"

