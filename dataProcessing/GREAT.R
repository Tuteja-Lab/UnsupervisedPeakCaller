#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}


library(rGREAT)

bed <- read.table(args[1], header = F)
bed <- bed[,1:4]
colnames(bed) <- c("chr", "start", "end", "name")
job <- submitGreatJob(bed, species = args[3])
tbl <- getEnrichmentTables(job)
bp <- tbl[[2]]
write.table(bp, args[4], row.names = F, sep = "\t", quote = F)

bp <- bp[bp$Binom_Fold_Enrichment >= 2 &
           bp$Binom_Adjp_BH <= 0.05 &
           bp$Hyper_Observed_Gene_Hits >= 5,]
write.table(bp, args[2], row.names = F, sep = "\t", quote = F)

