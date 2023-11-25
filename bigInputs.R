#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least four arguments: if not, return an error
if (length(args) < 4) {
	stop("At least five argument must be supplied.", call.=FALSE)
}
if (length(args) == 4) {
	ref.prefix = ""
} else {
	ref.prefix = args[5]
}

dir=args[1]
chr.file=args[2]
prefix=args[3]
out.file=args[4]

suppressPackageStartupMessages(library(dplyr))

final <- data.frame(V1=NA, V2=NA, V3=NA, V4=NA)
chr <- read.table(paste0(dir, "/", chr.file), header = F)

for (i in chr$V1) {
	# there are 1000 bp regions outside blacklist regions
	if (file.size(paste0(dir, "/", ref.prefix, i, "/", prefix, "_", i, "-regions4.txt")) > 0) {
		a <- read.table(paste0(dir, "/", ref.prefix, i, "/", prefix, "_", i, "-regions4.txt"), header = F)
		a1 <- sub("_.*", "", a$V4)
		# these are original regions passing threshold
		b <- read.table(paste0(dir, "/", ref.prefix, i, "/", prefix, "_", i, "-regions2.txt"), header = F)
		b$V4 <- paste0(b$V1, b$V4)
		b <- b[b$V4 %in% a1,]
		b <- b[,1:4]
		b <- distinct(b)
		#b$V4 <- paste0(b$V1, b$V4)
		final <- rbind(final, b)
	}
}
final <- final[!is.na(final$V1),]
write.table(final, out.file, quote = F, sep = "\t", row.names = F, col.names = F)

