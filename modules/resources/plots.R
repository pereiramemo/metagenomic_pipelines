###############################################################################
# 1. Set env
###############################################################################

library("tidyverse")
args = commandArgs(trailingOnly=TRUE)

###############################################################################
# 2. Load data
###############################################################################

TBL_FILE <- args[1]
OUTPUT_FILE <- args[2]

TBL <- read_tsv(TBL_FILE, col_names = F)


# TBL <- read_tsv("~/workspace/dev/indicators_contaminants_2018/data/metagenomics/prepro/sample_redu_1/stats.tsv",
#                 col_names = F)
# OUTPUT <- "~/workspace/dev/indicators_contaminants_2018/data/metagenomics/prepro/sample_redu_1/stats_plot.png"

colnames(TBL) <- c("sample","file","stat","value","order")

###############################################################################
# 3. Order
###############################################################################

TBL$order <- as.numeric(TBL$order)

###############################################################################
# 4. Plot
###############################################################################

lebelsx <- c("mean_length" = "Mean length", "num_seq" = "Number of sequences")

p <- ggplot(TBL, aes(x = reorder(file, order), y = value)) +
     facet_wrap(~stat, nrow = 2, ncol = 1, scales = "free_y", labeller = as_labeller(lebelsx)) +
     geom_bar(stat = "identity", alpha = 0.8) +
     scale_y_continuous(labels = scales::scientific) +
     theme_bw() +
     theme(
       axis.title.x = element_blank(),
       axis.title.y = element_blank(),
       axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
       plot.margin = margin(10,10,10,20),
       strip.background = element_rect(fill="white", size=1, color="darkblue"),
       strip.text = element_text(size = 10)
     )

###############################################################################
# 5. Save plot
###############################################################################

ggsave(p, filename = OUTPUT_FILE, device = "png", dpi = 500, width = 12, height = 8)

