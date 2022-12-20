#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}



library(ggplot2)
r <- read.table(args[1], header = T, sep = "\t", quote = "")
terms <- read.table(args[2], header = F, sep = "\t")

r <- r[r$name %in% terms$V1,]

pdf(args[3], width = 20, height = as.numeric(args[4]))
ggplot(data=r, aes(x=log2(`Binom_Fold_Enrichment`), y=name,
                   color=-log10(Binom_Adjp_BH), size=`Hyper_Observed_Gene_Hits`)) +
  geom_point(show.legend = TRUE) + ylab("") + xlab("log2(Fold change)") +
  scale_color_gradient(low = "#fb6a4a", high = "#a50f15") +
  scale_size_continuous(range = c(5, 20)) +
  theme(legend.text = element_text(size=15), legend.title = element_text(size=15), plot.title = element_text(size=20),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 15)) + ggtitle(args[1]) +
  geom_vline(xintercept=log2(2))
dev.off()


