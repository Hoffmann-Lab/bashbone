#! /usr/bin/env Rscript
# (c) Konstantin Riege
library("DEXSeq")
library("ggplot2")
library("gplots")
library("BiocParallel")

# samples <- "/path/to/a1.c","/path/to/a2.c","/path/to/a3.c","/path/to/b1.c","/path/to/b2.c","/path/to/b3.c"
# conditions <- "a","a","a","b","b","b"
# labels <- "a1","a2","a3","b1","b2","b3"
# libTypes <- "paired-end", "paired-end", "paired-end","single-end", "single-end", "single-end
# dexseqannotation <- /path/to/dexseqFlattened.gtf
# threads <- 48
# out <- /path/to/outdir
# profilescores <- /path/to/file (geneid score) - optional
args <- commandArgs(TRUE)

threads <- as.numeric(args[1])
out <- args[2]
samples <- unlist(strsplit(args[3],","))
labels <- unlist(strsplit(args[4],","))
conditions <- unlist(strsplit(args[5],","))
libTypes <- unlist(strsplit(args[6],","))
dexseqannotation <- args[7]
profilescores <- args[8] #optional

BPPARAM <- MulticoreParam(workers = threads)

sampleTable <- data.frame(row.names = labels, condition = conditions, libType = libTypes)
ddxDEXSeq <- DEXSeqDataSetFromHTSeq(countfiles = samples, sampleData = sampleTable, design = ~ sample + exon + condition:exon, flattenedfile = dexseqannotation)
file <- paste(out, "ddxDEXSeq.Rdata", sep = "/")
# DEXSeq does:
# ddx <- estimateSizeFactors(ddxDEXSeq)
# ddx <- estimateDispersions(ddx, BPPARAM = BPPARAM)
# plotDispEsts(ddx)
# ddx <- testForDEU(ddx, BPPARAM = BPPARAM)
# ddx <- estimateExonFoldChanges(ddx, BPPARAM = BPPARAM)
# ddxr <- DEXSeqResults(ddx)
ddxr <- DEXSeq(ddxDEXSeq, BPPARAM = BPPARAM)
file <- paste(out, "ddxr.Rdata", sep = "/")
save(ddxr, file = file)

ddxrNoNA <- ddxr[ !is.na(ddxr$padj) , ]
ddxrNoNA <- ddxrNoNA[order(ddxrNoNA$padj , -abs(ddxrNoNA[[10]])) , ] #padj and descending abs log2fc
ddxrNoNA05 <- ddxrNoNA[ ddxrNoNA$padj <= 0.05 , ]
ddxrNoNA01 <- ddxrNoNA[ ddxrNoNA$padj <= 0.01 , ]

csv <- paste(out, "full.csv", sep = "/")
o <- as.data.frame(ddxr)
o$transcripts <- ''
# unless removed, output file format gets corrupted
# "ENSG00000000001:E001","ENSG00000000001",...,c("ENST00000000001" "ENST00000000002")
# gets
# "ENSG00000000001:E001","ENSG00000000001",...,c("ENST00000000001" 
# "ENST00000000002")
write.csv(o, file = csv)

csv <- paste(out, "filtered_padjNA.csv", sep = "/")
o <- as.data.frame(ddxrNoNA)
o$transcripts <- ''
write.csv(o, file = csv)

csv <- paste(out, "filtered_padj05.csv", sep = "/")
o <- as.data.frame(ddxrNoNA05)
o$transcripts <- ''
write.csv(o, file = csv)

csv <- paste(out, "filtered_padj01.csv", sep = "/")
o <- as.data.frame(ddxrNoNA01)
o$transcripts <- ''
write.csv(o, file = csv)

pdf(paste(out, "ma_plot.pdf", sep = "/"))
plotMA(ddxr)
graphics.off()

IDsNoNA <- row.names(ddxrNoNA05)
tmp <- unlist(strsplit(IDsNoNA,":"))
IDsNoNA <- tmp[seq(1,length(tmp),2)]
uniques <- numeric()
for (i in IDsNoNA) {
	if(length(strsplit(i, "\\+")[[1]])  == 1){
		uniques <- unique(c(uniques, i))
		if(length(uniques)==50){
			break
		}
	}
}

o <- file.path(out, "top50logFC")
dir.create(o, showWarnings = FALSE)
for (i in uniques) {
	postscript(file.path(o, paste(i, ".ps", sep = "")))
	plotDEXSeq(ddxr, i, legend = TRUE, lwd = 1)
	graphics.off()
    err <- tryCatch({
        postscript(file.path(o, paste(i, "_transcripts.ps", sep = "")))
        plotDEXSeq(ddxr, i, legend = TRUE, displayTranscripts=TRUE, lwd = 1)
    }, error=function(e){
        print(paste(e, "warning: cannot create tranctipt plot for ", i, sep="" ))
    })
	graphics.off()
}

if (!is.na(profilescores)){
	o <- file.path(out, "top50profilescores")
	dir.create(o, showWarnings = FALSE)
	IDs <- read.csv(profilescores, header = FALSE, sep = "\t", dec = ".", stringsAsFactors = FALSE, strip.white = T)[[1]]
	uniques <- numeric()
	for (i in IDs) {
		basemeans <- ddxr[ ddxr$groupID==i , ]$exonBaseMean
		if(length(strsplit(i, "\\+")[[1]]) == 1 && sum(basemeans)/length(basemeans)>=10 && i %in% IDsNoNA) {
			uniques <- unique(c(uniques, i))
			if(length(uniques)==50){
				break
			}
		}
	}

	for (i in uniques) {
		postscript(file.path(o, paste(i, ".ps", sep = "")))
		plotDEXSeq(ddxr, i, legend = TRUE, lwd = 1)
		graphics.off()
	    err <- tryCatch({
	        postscript(file.path(o, paste(i, "_transcripts.ps", sep = "")))
	        plotDEXSeq(ddxr, i, legend = TRUE, displayTranscripts=TRUE, lwd = 1)
	    }, error=function(e){
	        print(paste(e, "warning: cannot create tranctipt plot for ", i, sep="" ))
	    })
		graphics.off()
	}
}