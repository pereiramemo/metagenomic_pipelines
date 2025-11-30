# Metagenomic pipelines
This repository contains the scripts used to preprocess and quality-check metagenomic data: `preprocess_pipeline.sh` and `quality_check.R`, written in Bash and R, respectively.

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

## preprocess_pipeline.sh

[preprocess_pipeline.sh](https://github.com/pereiramemo/metagenomic_pipelines/blob/main/preprocess/preprocess_pipeline.sh): This bash pipeline preprocesses raw Illumina paired-end reads from metagenomic samples. The main tasks include:
- Checking for the presence of adapters and trimming them
- Merging paired-end reads using PEAR or BBMerge
- Quality trimming of merged and unmerged reads
- Computing and plotting read statistics (number and mean length) for intermediate files

The output consists of:
- ```*workable.fasta```: FASTA file ready for downstream analyses (merged reads)
- ```stats.tsv```: Table with read statistics
- ```stats_plots.png```: Plot showing number of sequences and mean read length

Optionally, the pipeline can output quality-checked paired-end reads as FASTQ files.

To see the help run ```./preprocess_pipeline.sh --help```

```
Usage: ./preprocess_pipeline.sh <options>

Options:
	--help
		print this help

	--clean=t|f
		remove all intermediate files [default=f]

	--compress=t|f
		output data as .gz files [default=f]

	--merger=CHAR
		tool to merge paired-end reads, one of "pear" or "bbmerge" [default=pear]

	--min_length=NUM
		minimum length of (PE or merged) reads [default=75]

	--min_overlap=NUM
		minimum overlap to merge paired-end reads with pear

	--min_qual=NUM
		minimum quality score to trim reads [default=20]

	--nslots=NUM
		number of threads used [default=12]

	--output_dir=CHAR
		directory to output generated data (preprocessed data, plots, tables)

	--output_pe=t|f
		output quality checked paired-end reads as fastq files [default=f]

	--output_merged=t|f
		output quality checked merged reads as fasta file (workable.fasta) [default=t]

	--overwrite=t|f
		overwrite previous directory [default=f]

	--pvalue=NUM
		p value used to run pear [default=0.01]

	--plot=t|f
		create statistics barplot [default=f]

	--reads=CHAR
		input R1 reads

	--reads2=CHAR
		input R2 reads

	--sample_name=CHAR
		sample name [default=metagenomex]

	--subsample=t|f
		subsample metagenome to 10K to test execution [default=f]

	--trim_adapters=t|f
		check for adapters and trim [default=f]
```

## quality_check.R

[quality_check.R](https://github.com/pereiramemo/metagenomic_pipelines/blob/main/preprocess/quality_check.R): This R script performs comprehensive quality assessment of raw Illumina paired-end reads. It generates multiple plots to evaluate read quality and identify potential issues. The main analyses include:
- Calculation of mean quality scores for R1 and R2 reads
- Plotting quality scores versus read counts
- Generating histograms of read count distributions
- Detection and quantification of PhiX contamination

The output consists of five PNG files:
- ```r1_mean_q_vs_nseq.png```: Quality score vs read count for R1
- ```r2_mean_q_vs_nseq.png```: Quality score vs read count for R2
- ```samples_hist.png```: Histogram of read counts per sample
- ```samples_hist_log.png```: Log-transformed histogram of read counts
- ```samples_perc_phix_barplot.png```: PhiX contamination levels per sample

To see the help run ```Rscript quality_check_plots.R --help```

```
Usage: quality_check_plots.R [options]

Options:
	--input_dir=CHARACTER
		Input directory with FASTQ files

	--output_dir=CHARACTER
		Output directory for plots

	--nslots=INTEGER
		Number of threads to use [default=12]

	--r1_pattern=CHARACTER
		Pattern for R1 FASTQ files [default=R1_001.fastq.gz]

	--r2_pattern=CHARACTER
		Pattern for R2 FASTQ files [default=R2_001.fastq.gz]

	--overwrite=LOGICAL
		Overwrite previous output [default=FALSE]

	-h, --help
		Show this help message and exit
```

# **Dependencies**:  
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

# **License**

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

Copyright (C) 2025 Emiliano Pereira

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.  



