###############################################################################
### 1. Set env
###############################################################################

library(tidyverse)
library(ShortRead)
library(doParallel)
library(dada2)

args <- commandArgs(trailingOnly=TRUE)

###############################################################################
### 2. Load data
###############################################################################

# INPUT_DIR <- "/home/epereira/workspace/dev/samo/metagenomics/data_redu/"
# OUTPUT_DIR <- "/home/epereira/workspace/dev/samo/metagenomics/results/quality_analysis/"
# NSLOTS <- 40
# registerDoParallel(cores = NSLOTS)
# PATTERN_R1 <- "_R1_clipped_redu.fastq"
# PATTERN_R2 <- "_R2_clipped_redu.fastq"

INPUT_DIR <- args[1]
OUTPUT_DIR <- args[2]
NSLOTS <- args[3] %>% as.numeric
registerDoParallel(cores = NSLOTS)
PATTERN_R1 <- args[4]
PATTERN_R2 <- args[5]

rawR1 <- sort(list.files(INPUT_DIR, pattern = PATTERN_R1, full.names = T))
rawR2 <- sort(list.files(INPUT_DIR, pattern = PATTERN_R2, full.names = T))
SAMPLE_NAMES <- basename(rawR1) %>%
                sub(pattern = "_.*", replacement = "")

###############################################################################
### 3. R1 count
###############################################################################

count_seqs <- function(p) {
  n <- readFastq(p) %>% sread %>% as.character() %>% length()
  sample <- basename(p) %>% sub(x = ., pattern = "_R1.*", replacement = "")
  output <- data.frame(sample = sample, nseq = n)
  return(output)
}  
  
seq_counts_df <- foreach(i = rawR1, .combine=rbind) %dopar% {
  count_seqs(i)
}

###############################################################################
### 4. R1 quality vs nseq
###############################################################################

x_r1 <- qa(dirPath = INPUT_DIR, pattern = PATTERN_R1, sample = T, n = 5000)
qa_df <- x_r1[["perCycle"]][["quality"]]

qa_means <- qa_df %>%
            group_by(lane) %>%
            summarize(mean_q = sum(Score*Count)/sum(Count))

qa_means$lane <- qa_means$lane %>% sub(x = ., pattern = "_R1.*", replacement = "")
qa_means2counts <- left_join(x = qa_means, y = seq_counts_df, by = c("lane" = "sample"))

text_size <- 2
p_r1 <- ggplot(data = qa_means2counts, aes(x = mean_q, y = nseq)) +
        geom_point() +
        scale_y_log10() +
        ylab("Read counts (log)") +
        xlab("Mean quality score (R1)") +
        geom_text(aes(label=as.character(lane)), hjust=0.5, vjust=-1, size = text_size)


file_p_r1 <- paste(OUTPUT_DIR,"r1_mean_q_vs_nseq.png", sep = "/")
ggsave(p_r1, filename = file_p_r1, 
       device = "png", width = 5, height = 4, dpi = 300)
  
###############################################################################
### 4. R2 quality vs nseq
###############################################################################

x_r2 <- qa(dirPath = INPUT_DIR, pattern = PATTERN_R2, sample = T, n = 5000)
qa_df <- x_r2[["perCycle"]][["quality"]]
qa_means <- qa_df %>%
            group_by(lane) %>%
            summarize(mean_q = sum(Score*Count)/sum(Count))

qa_means$lane <- qa_means$lane %>% sub(x = ., pattern = "_R2.*", replacement = "")
qa_means2counts <- left_join(x = qa_means, y = seq_counts_df, by = c("lane" = "sample"))

text_size <- 2
p_r2 <- ggplot(data = qa_means2counts, aes(x = mean_q, y = nseq)) +
        geom_point() +
        scale_y_log10() +
        ylab("Read counts (log)") +
        xlab("Mean quality score (R2)") +
        geom_text(aes(label=as.character(lane)), hjust=0.5, vjust=-1, size = text_size)

file_p_r2 <- paste(OUTPUT_DIR,"r2_mean_q_vs_nseq.png", sep = "/")
ggsave(p_r2, filename = file_p_r2, 
       device = "png", width = 5, height = 4, dpi = 300)

###############################################################################
### 4. Plot nseq hist
###############################################################################

samples_hist_p <- ggplot(data = qa_means2counts, aes(nseq)) +
                  geom_histogram(bins = 40) +
                  ylab("Num of samples")

samples_hist_log_p <- ggplot(data = qa_means2counts, aes(log(nseq))) +
                      geom_histogram(bins = 40) +
                      ylab("Num of samples")

file_p_hist <- paste(OUTPUT_DIR,"samples_hist.png", sep = "/")
ggsave(samples_hist_p, filename = file_p_hist, 
       device = "png", width = 5, height = 4, dpi = 300)

file_p_hist_log <- paste(OUTPUT_DIR,"samples_hist_log.png", sep = "/")
ggsave(samples_hist_log_p, filename = file_p_hist_log, 
       device = "png", width = 5, height = 4, dpi = 300)


###############################################################################
### 5. Estimate % of Phix sequences
###############################################################################

count_phix_seqs <- function(p) {

  fastq <- readFastq(p) %>% sread %>% as.character()
  fastq_nphix <- isPhiX(seqs = fastq, wordSize = 16, minMatches = 2)
  output <- data.frame(file = p, n = length(fastq), nphix = sum(fastq_nphix))
  return(output)
}

phix_counts_df <- foreach(i = rawR1, .combine=rbind) %dopar% {
  count_phix_seqs(i)
}

###############################################################################
### Plot % of Phix sequences
###############################################################################

text_size <- 6

X <- data.frame(sample = phix_counts_df$file %>% sub(x = ., pattern = ".*/", replacement = ""),
                perc = 100*phix_counts_df$nphix/phix_counts_df$n,
                nphix = phix_counts_df$nphix + 1)

# nphix_barplot <- ggplot(data = X, aes(x = sample, y = nphix)) +
#                  geom_bar(stat = "identity") +
#                  scale_y_log10() +
#                  theme_bw() +
#                  theme(
#                    axis.text.x = element_text(size = text_size, angle = 45, hjust = 1)
#                  )
                
i <- X$sample %>% grep(pattern = "Undetermined")
if (length(i) > 0) {
  X <- X[-i,]
}

j <- X$perc > 0.005
X$color <- "gray50" 
X$color[j] <- "indianred"

perc_phix_barplot <- ggplot(data = X, aes(x = sample, y = perc, fill =  color)) +
                      geom_bar(stat = "identity") +
                      theme_bw() +
                      scale_fill_manual(values = c("gray50","indianred"), 
                                        labels = c("Phix <= 0.005%","Phix > 0.005%"),
                                        name = "") +
                      theme(
                        axis.text.x = element_text(size = text_size, angle = 45, hjust = 1),
                        legend.position = "top"
                      )

file_perc_phix_barplot <- paste(OUTPUT_DIR,"samples_perc_phix_barplot.png", sep = "/")
ggsave(perc_phix_barplot, filename = file_perc_phix_barplot, 
       device = "png", width = 10, height = 4, dpi = 300)

