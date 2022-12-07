#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}



library(dplyr)

final <- data.frame(V1=NA, V2=NA, V3=NA, V4=NA)
chr <- read.table(args[1], header = F)

for (i in chr$V1) {
if (file.size(paste0(args[2], "/chr", i, "/", i, "regions2.txt")) > 0) {
a <- read.table(paste0(args[2], "/chr", i, "/", i, "regions2.txt"), header = F)
a1 <- sub("_.*", "", a$V4)
b <- read.table(paste0(args[2], "/chr", i, "/", "temp", i, ".txt"), header = F)
b$V4 <- paste0(b$V1, b$V4)
b <- b[b$V4 %in% a1,]
b <- b[,1:4]
b <- distinct(b)
#b$V4 <- paste0(b$V1, b$V4)
final <- rbind(final, b)
}
}
final <- final[!is.na(final$V1),]
write.table(final, paste0(args[2], "/bigInputs.txt"), quote = F, sep = "\t", row.names = F)

