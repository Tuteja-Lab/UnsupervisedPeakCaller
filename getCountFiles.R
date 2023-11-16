#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

#full new
process3_1 <- function(sub, s) {
  sub <- dplyr::distinct(sub[,1:3])
  sub$mid <- round((sub$end - sub$start)/2) + sub$start
  sub$start <- sub$mid - inputLength/2
  sub$end <- sub$mid + inputLength/2
  sub <- sub[,1:3]
  sub$name <- s
  colnames(sub) <- c("chr", "start", "end", "name")

  return(sub)
}

getSegments <- function(s) {
  sub <- file[file$name == s,]
 
  if (sub$end[1] - sub$start[1] < inputLength) {
  	sub <- process3_1(sub, s)
  } else {
  	sub1 <- sub[which(sub$counts >= quantile(sub$counts, 0.95)),]
  	sub1c <- sub1
  	sub1 <- sub1[,c("chr.c", "start.c", "end.c")]
  	colnames(sub1) <- c("chr", "start", "end")
  
  	if (nrow(sub1) > 1) {
    		sub1.sort   <- bedr.sort.region(sub1, check.chr = FALSE)
    		sub1.merge <- bedr.merge.region(sub1.sort, distance = inputLength, verbose = T, check.chr = FALSE)
    		sub1.merge$mid <- round((sub1.merge$end - sub1.merge$start)/2) + sub1.merge$start
    		sub1.merge$start <- sub1.merge$mid - inputLength/2
    		sub1.merge$end <- sub1.merge$mid + inputLength/2
    
		for (i in 1:nrow(sub1.merge)) {
		      sub1.merge$name[i] <- paste0(s, "_", i)
		}
  	sub <- sub1.merge[,c("chr", "start", "end", "name")]
  	} else {
    		sub$mid <- round((sub$end - sub$start)/2) + sub$start
    		sub$start <- sub$mid - inputLength/2
    		sub$end <- sub$mid + inputLength/2
    		sub <- sub[,1:4]
    		colnames(sub) <- c("chr", "start", "end", "name")
  	}
   }	
 
  return(sub)
}



suppressPackageStartupMessages(library(bedr))
suppressPackageStartupMessages(library(doParallel))

print(args[1])

if (file.size(args[1]) > 0) {
file <- read.table(args[1], header = F)
colnames(file) <- c("chr", "start", "end", "name", "chr.c", "start.c", "end.c", "counts")

if (length(grep("chr", file$chr[1])) == 0) {
  print("Adding `chr` to chromosome name")
  file$name <- paste0(file$chr, file$name)
  file$chr <- paste0("chr", file$chr)
  file$chr.c <- paste0("chr", file$chr.c)
}

inputLength <- as.numeric(args[2])

no_cores <- 10
cl <- makeCluster(no_cores, type = "FORK")
registerDoParallel(cl)
new <- data.frame(matrix(ncol = 4))
colnames(new) <- c("chr", "start", "end", "name")

result <- foreach(i=unique(file$name)) %dopar% {getSegments(i)}
#save(result, file=args[3])

for (i in 1:length(result)) {
        new <- rbind(new, result[[i]])
}
new <- new[!is.na(new$chr),]
write.table(new, file=args[3], quote = F, sep = "\t", row.names = F)
} else {
print("Input has 0 line.")
}

