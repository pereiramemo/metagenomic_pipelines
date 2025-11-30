# Metagenomic pipelines
This repository contains the code to preprocess and quality check metagenomic data.

# Repository structure

```
.
├── LICENSE
├── README.md
├── environment.yml                          # Conda environment specification
└── preprocess/                              # Preprocessing pipeline
    ├── preprocess_pipeline.sh               # Main preprocessing pipeline script
    ├── preprocess_pipeline_conf.sh          # Configuration file with tool paths and functions
    ├── quality_check_plots_runner.sh        # Script to generate quality check plots
    └── resources/                           # Additional resources and scripts
        ├── fq2fa.sh                         # FASTQ to FASTA conversion script
        ├── plots.R                          # R script for generating statistics plots
        └── quality_check_plots.R            # R script for quality assessment plots
```

**Dependencies**:  
[bzip2](http://www.bzip.org)  
[gzip](https://www.gzip.org)  
[seqtk](https://github.com/lh3/seqtk)  
[BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide)  
[PEAR](https://cme.h-its.org/exelixis/web/software/pear)  
[R](https://www.r-project.org)  
[tidyverse](https://www.tidyverse.org) R package  
[ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) R/Bioconductor package  
[doParallel](https://cran.r-project.org/web/packages/doParallel/index.html) R package  
[dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html) R/Bioconductor package  

**Installation with mamba**:  
All dependencies can be installed using mamba (or conda). 

First, check if mamba is installed:
```bash
command -v mamba
```

If mamba is not installed, you can install it via:
- **Miniforge** (recommended): [https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge)
- **Mambaforge**: [https://github.com/conda-forge/miniforge#mambaforge](https://github.com/conda-forge/miniforge#mambaforge)
- **Or install mamba into an existing conda installation**:
  ```bash
  conda install -n base -c conda-forge mamba
  ```

Once mamba is installed, create a new environment using the provided `environment.yml` file:
```bash
mamba env create -f environment.yml
```

Then activate the environment:
```bash
mamba activate metagenomic_pipeline
```

# **How to use**

The folder **preprocess** contains a pipeline programmed in BASH, that can be executed from the command line to preprocess raw Illumina paired-end reads obtained from metagenomic samples (named [preprocess_pipeline.sh](https://github.com/pereiramemo/metagenomic_pipelines/blob/main/preprocess/preprocess_pipeline.sh)). The main tasks consist of checking for the presence of adapters, merging the paired-end reads, and quality trimming the merged and unmerged reads. It additionally computes and plots the number of reads and mean read length of the intermediate files to trace the preprocessing tasks and detect potential irregularities. The output consists of a fasta file (i.e., ```*workable.fasta```), ready to use in downstream analyses, and table, and a plot of the number of sequences and mean read length of the intermediate files (i.e., ```stats.tsv``` and ```stats_plots.png```).  
Optionally, the pipeline can be used to produce quality checked paired-end reads.


To see the help run ```./preprocess_pipeline.sh --help```

```
Usage: ./preprocess_pipeline.sh <options>
--help                          print this help
--clean t|f                     remove all intermediate files (default f)
--compress t|f                  output data as .gz files (default f)
--merger CHAR                   tool to merge paired-end reads, one of "pear" of "bbmerge" (default "pear")
--min_length NUM                minimum length of (PE or merged) reads (default 75)
--min_overlap NUM               minimum overlap to merge paired-end reads with pear
--min_qual NUM                  minimum quality score to trim reads (default 20)
--nslots NUM                    number of threads used (default 12)
--output_dir CHAR               directory to output generated data (i.e., preprocessed data, plots, tables)
--output_pe t|f                 output quality checked paired-end reads as fastq files (default f)
--output_merged t|f             output quality checked merged reads as fasta file (i.e., workable.fasta) (default t)
--overwrite t|f                 overwrite previous directory (default f)
--pvalue NUM                    p value used to run pear. See pear help for valid p values (default: 0.01)
--plot t|f                      create statistics barplot (default f)
--reads CHAR                    input R1 reads
--reads2 CHAR                   input R2 reads
--sample_name CHAR              sample name (default metagenomex)
--subsample t|f                 subsample metagenome to 10K to test execution (default f)
--trim_adapters t|f             check for adapters and trim (default f)
```

