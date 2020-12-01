# Metagenomic pipelines
This repository contains the code to preprocess and analyze metagenomic data, organized in two different pipelines. 
The first one, *preprocess*, contains a pipeline that can be executed from the command line to preprocess raw Illumina paired-end reads obtained from metagenomic samples. The main tasks consist of checking for the presence of adapters, merging the paired-end reads, and quality trimming the merged and unmerged reads.
It additionally computes and plots the number of reads and mean read length of the intermediate files to trace the preprocessing tasks and detect potential irregularities.

The output consists of a workable.fasta file, ready to use in downstream analyses, and table, and a plot of the number of sequences and mean read length of the intermediate files. 
