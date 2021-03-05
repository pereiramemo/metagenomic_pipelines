# Metagenomic pipelines
This repository contains the code to preprocess and analyze metagenomic data, organized in two different pipelines: **preprocess_pipeline.bash**, and **mg_traits_lite.bash**. 

1. The folder **preprocess** contains a pipeline programmed in BASH, that can be executed from the command line to preprocess raw Illumina paired-end reads obtained from metagenomic samples (named [preprocess_pipeline.bash](https://github.com/pereiramemo/metagenomic_pipelines/blob/main/preprocess/preprocess_pipeline.bash)). The main tasks consist of checking for the presence of adapters, merging the paired-end reads, and quality trimming the merged and unmerged reads. It additionally computes and plots the number of reads and mean read length of the intermediate files to trace the preprocessing tasks and detect potential irregularities. The output consists of a fasta file (i.e., ```*workable.fasta```), ready to use in downstream analyses, and table, and a plot of the number of sequences and mean read length of the intermediate files (i.e., ```stats.tsv``` and ```stats_plots.png```).  
Optionally, the pipeline can be used to produce quality checked paired-end reads.

    **Dependencies**:  
[bzip2](http://www.bzip.org)  
[gzip](https://www.gzip.org)  
[seqtk](https://github.com/lh3/seqtk)  
[BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide)  
[PEAR](https://cme.h-its.org/exelixis/web/software/pear)  
[R](https://www.r-project.org)  
[tidyverse](https://www.tidyverse.org) R package  

    Once these tools are installed, the [preprocess_pipeline_conf.bash](https://github.com/pereiramemo/metagenomic_pipelines/blob/main/preprocess/preprocess_pipeline_conf.bash) has to be edited to set the correct paths.

    To see the help run ```./preprocess_pipeline.bash --help```

```
Usage: ./preprocess_pipeline.bash <options>
--help                          print this help  
--clean t|f                     remove all intermediate files  
--merger CHAR                   tool to merge paired-end reads, one of "pear" of "bbmerge" (default "pear")  
--min_qual NUM                  minimum quality score to trim reads (default 20)  
--min_overlap NUM               minimum overlap to merge paired-end reads with pear  
--nslots NUM                    number of threads used (default 12)  
--output_dir CHAR               directory to output generated data (i.e., preprocessed data, plots, tables)  
--output_pe t|f                 output quality checked paired-end reads as fastq files (default f)  
--output_merged t|f             output quality checked merged reads as fasta files (i.e., workable.fasta) (default t)  
--overwrite t|f                 overwrite previous directory (default f)
--pvalue NUM                    p value used to run pear. See pear help for valid p values (default: 0.01)  
--reads CHAR                    input R1 reads  
--reads2 CHAR                   input R2 reads  
--sample_name CHAR              sample name (default metagenomex)  
--subsample t|f                 subsample metagenome to 10K to test execution (default f)  
--trim_adapters t|f             check for adapters and trim (default f)  
```

2. The folder **mg_traits** contains the pipeline [mg_traits_lite](https://github.com/pereiramemo/metagenomic_pipelines/blob/main/mg_traits/mg_traits_lite.bash). This tool is a command line application programmed in BASH, AWK, and R, dedicated to the computation of 25 (and counting) metagenomic traits, ranging from GC variance and amino acid composition to functional diversity and average genome size. It takes as an input a preprocessed (unassembled) metagenomic sample and outputs the computed metagenomic traits organized in different tables and grouped in separate folders according to the type of data source
(see Fig. 1). 

    **Dependencies**:  
[BBDuk](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide)  
[PEAR](https://cme.h-its.org/exelixis/web/software/pear)  
[FragGeneScan](https://omics.informatics.indiana.edu/FragGeneScan/)  
[UProC](http://uproc.gobics.de/)  
[EMBOSS](http://emboss.sourceforge.net/)  
[VSEARCH](https://github.com/torognes/vsearch)  
[AGS and ACN tools](https://github.com/pereiramemo/AGS-and-ACN-tools)  
[HMMER](http://hmmer.org)  
[R](https://www.r-project.org)  
[tidyverse](https://www.tidyverse.org) R package  


    Similarly as mentioned above, once these tools are installed the configuration file [mg_traits_lite_conf.bash](https://github.com/pereiramemo/metagenomic_pipelines/blob/main/mg_traits/mg_traits_lite_conf.bash) has to be edited to set the correct paths to define the variables.

    To see the help run ```./mg_traits_lite.bash --help```  


```
Usage: ./mg_traits_lite.bash <options>
--help                          print this help
--clean t|f                     remove all intermediate files
--confidence NUM                confidence value to run rdp bayes classifier (from 0 to 100; default 50)
--evalue_acn NUM                evalue to filter reads for for AGS computaton (default 1e-15)
--evalue_div NUM                evalue to filter reads for diversity estimation (default 1e-15)
--input_file CHAR               input workable fasta file
--nslots NUM                    number of threads used (default 12)
--max_length NUM                maximum read length used to trim reads (from the 3' end) for AGS computaton (default 180)
--min_length NUM                minimum read length used to estimate taxonomic diversity (default 100)
--output_dir CHAR               directory to output generated data (i.e., preprocessed data, plots, tables)
--overwrite t|f                 overwrite previous directory (default f)
--ref_db CHAR                   refernce database to run NBC (default silva_nr99_v138_train_set.fa.gz) 
--sample_name CHAR              sample name (default metagenomex)
--train_file_name CHAR          train file name to run FragGeneScan, see FragGeneScan help for options (default illumina_5)
```



![Figure 1](mg_traits/figures/Mg-Traits2.png)