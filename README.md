# Metagenomic pipelines
This repository contains scripts for quality checking, preprocessing, assembling, and mapping metagenomic data.

# Repository structure

```
.
├── LICENSE
├── README.md
├── environment.yml                          # Conda environment specification
└── modules/                                 # Pipeline modules
    ├── 1.1-quality_check_fastp.sh           # Quality check using fastp
    ├── 1.2-quality_check.R                  # Quality check plots using R
    ├── 2-preprocess_pipeline.sh             # Main preprocessing pipeline script
    ├── 3-assembly_and_map_pipeline.sh       # De novo assembly and read mapping
    ├── conf.sh                              # Configuration file with tool paths
    └── resources/                           # Additional resources and scripts
        ├── fq2fa.sh                         # FASTQ to FASTA conversion script
        └── plots.R                          # R script for generating statistics plots
```

# **Installation instructions**

**Clone the repository**:  
```bash
git clone https://github.com/pereiramemo/metagenomic_pipelines.git
cd metagenomic_pipelines
```

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

## 1.1-quality_check_fastp.sh

[1.1-quality_check_fastp.sh](modules/1.1-quality_check_fastp.sh): This script performs comprehensive quality assessment of raw Illumina paired-end reads using fastp. It generates HTML and JSON reports for quality control.

**Main features:**
- Quality assessment without modifying input files
- HTML and JSON report generation
- Configurable quality and length thresholds
- Batch processing of multiple samples
- PhiX contamination detection

To see the help run ```./modules/1.1-quality_check_fastp.sh --help```

**Key options:**
- `--input_dir=CHAR`: Directory containing input FASTQ files (required)
- `--output_dir=CHAR`: Directory to output QC reports and plots (required)
- `--pattern_r1=CHAR`: Pattern for R1 FASTQ files [default=_R1_001.fastq.gz]
- `--pattern_r2=CHAR`: Pattern for R2 FASTQ files [default=_R2_001.fastq.gz]
- `--nslots=NUM`: Number of threads to use [default=12]
- `--min_length=NUM`: Minimum read length filter [default=50]
- `--qualified_quality_phred=NUM`: Minimum quality value [default=20]
- `--overwrite=t|f`: Overwrite previous directory [default=f]

## 1.2-quality_check.R

[1.2-quality_check.R](modules/1.2-quality_check.R): This R script performs comprehensive quality assessment of raw Illumina paired-end reads. It generates multiple plots to evaluate read quality and identify potential issues.

**Main analyses:**
- Calculation of mean quality scores for R1 and R2 reads
- Plotting quality scores versus read counts
- Generating histograms of read count distributions
- Detection and quantification of PhiX contamination

**Output files:**
- `r1_mean_q_vs_nseq.png`: Quality score vs read count for R1
- `r2_mean_q_vs_nseq.png`: Quality score vs read count for R2
- `samples_hist.png`: Histogram of read counts per sample
- `samples_hist_log.png`: Log-transformed histogram of read counts
- `samples_perc_phix_barplot.png`: PhiX contamination levels per sample

To see the help run ```Rscript modules/1.2-quality_check.R --help```

**Key options:**
- `--input_dir=CHARACTER`: Input directory with FASTQ files (required)
- `--output_dir=CHARACTER`: Output directory for plots (required)
- `--nslots=INTEGER`: Number of threads to use [default=12]
- `--r1_pattern=CHARACTER`: Pattern for R1 FASTQ files [default=R1_001.fastq.gz]
- `--r2_pattern=CHARACTER`: Pattern for R2 FASTQ files [default=R2_001.fastq.gz]
- `--overwrite=LOGICAL`: Overwrite previous output [default=FALSE]

## 2-preprocess_pipeline.sh

[2-preprocess_pipeline.sh](modules/2-preprocess_pipeline.sh): This bash pipeline preprocesses raw Illumina paired-end reads from metagenomic samples.

**Main tasks:**
- Checking for the presence of adapters and trimming them
- Merging paired-end reads using PEAR or BBMerge
- Quality trimming of merged and unmerged reads
- Computing and plotting read statistics (number and mean length) for intermediate files

**Output files:**
- `*workable.fasta`: FASTA file ready for downstream analyses (merged reads)
- `stats.tsv`: Table with read statistics
- `stats_plots.png`: Plot showing number of sequences and mean read length

Optionally, the pipeline can output quality-checked paired-end reads as FASTQ files.

To see the help run ```./modules/2-preprocess_pipeline.sh --help```

**Required options:**
- `--reads FILE`: R1 file
- `--reads2 FILE`: R2 file
- `--output_dir DIR`: Output directory

**Optional parameters:**
- `--clean t|f`: Remove intermediates [default=f]
- `--compress t|f`: Compress outputs with pigz [default=f]
- `--merger STR`: pear|bbmerge [default=pear]
- `--min_length NUM`: Minimum read length after trimming [default=75]
- `--min_overlap NUM`: Minimum PE overlap for PEAR [default=10]
- `--min_qual NUM`: Quality trim threshold [default=20]
- `--nslots NUM`: Threads [default=12]
- `--output_pe t|f`: Output QC'ed paired-end reads [default=f]
- `--output_merged t|f`: Output merged QC'ed reads [default=t]
- `--overwrite t|f`: Replace existing output dir [default=f]
- `--pvalue NUM`: p-value for PEAR [default=0.01]
- `--plot t|f`: Produce QC plots [default=f]
- `--sample_name STR`: Name prefix [default=metagenomex]
- `--seed NUM`: Random seed for subsampling [default=123]
- `--subsample t|f`: Subsample to 10k reads [default=f]
- `--trim_adapters t|f`: Remove adapters [default=f]

## 3-assembly_and_map_pipeline.sh

[3-assembly_and_map_pipeline.sh](modules/3-assembly_and_map_pipeline.sh): This pipeline performs de novo assembly of metagenomic reads and maps reads back to the assembled contigs.

**Main tasks:**
- De novo assembly using MEGAHIT
- Read mapping using BWA
- Duplicate removal using Picard (optional)
- Contig length filtering

To see the help run ```./modules/3-assembly_and_map_pipeline.sh --help```

**Required options:**
- `--reads1 CHAR`: Input R1 metagenome data (fastq/fa)
- `--reads2 CHAR`: Input R2 metagenome data (fastq/fa)
- `--sample_name CHAR`: Sample name used to name the files

**Optional parameters:**
- `--assem_dir CHAR`: Directory with previously computed assemblies (format: dirname/SAMPLE_NAME/SAMPLE_NAME.contigs.fa)
- `--assem_preset CHAR`: MEGAHIT preset to generate assembly [default=meta-sensitive]
- `--nslots NUM`: Number of threads used [default=12]
- `--min_contig_length NUM`: Minimum length of contigs to keep [default=250]
- `--output_dir CHAR`: Output directory [default=mg-clust_output-1]
- `--overwrite t|f`: Overwrite previous folder if present [default=f]
- `--remove_duplicates t|f`: Remove PCR duplicates with Picard [default=f]

# **Dependencies**

All dependencies are specified in the [environment.yml](environment.yml) file and can be installed via mamba/conda.

**Core tools:**
- [bzip2](http://www.bzip.org) - File compression
- [gzip](https://www.gzip.org) - File compression
- [pigz](https://zlib.net/pigz/) - Parallel gzip compression
- [seqtk](https://github.com/lh3/seqtk) - Sequence processing toolkit
- [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) - Suite including BBDuk and BBMerge for adapter trimming and read merging
- [PEAR](https://cme.h-its.org/exelixis/web/software/pear) - Paired-end read merger
- [fastp](https://github.com/OpenGene/fastp) - Fast all-in-one preprocessing tool
- [MEGAHIT](https://github.com/voutcn/megahit) - De novo assembler for metagenomes
- [BWA](https://github.com/lh3/bwa) - Burrows-Wheeler Aligner for read mapping
- [SAMtools](http://www.htslib.org/) - SAM/BAM file manipulation
- [Picard](https://broadinstitute.github.io/picard/) - Java tools for manipulating sequencing data
- [EMBOSS](http://emboss.sourceforge.net/) - Sequence analysis tools (provides infoseq)

**R and R packages:**
- [R](https://www.r-project.org) - Statistical computing environment
- [tidyverse](https://www.tidyverse.org) - Data manipulation and visualization
- [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) - FASTQ file handling (Bioconductor)
- [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html) - Parallel processing
- [dada2](https://bioconductor.org/packages/release/bioc/html/dada2.html) - Sequence quality profiling (Bioconductor)
- [optparse](https://cran.r-project.org/web/packages/optparse/index.html) - Command-line argument parsing  

# **License**

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

Copyright (C) 2025 Emiliano Pereira

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  



