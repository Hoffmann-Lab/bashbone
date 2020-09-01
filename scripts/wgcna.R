#! /usr/bin/env Rscript
# (c) Arne Sahm, Konstantin Riege
suppressMessages(library('WGCNA'))
suppressMessages(library('DGCA'))
suppressMessages(library('ggplot2'))

args = commandArgs(TRUE)
# threads = as.numeric(args[1])
# may cause Error in sendMaster(try(lapply(X = S, FUN = FUN, ...), silent = TRUE))
# The error seems to be from sendMaster which is used by the multicore package to send results from child processes back to the master process. From this StackOverflow thread, it seems that the results of the child processes are too large for mclapply to send back to the master.
# Dev uses foreach, which, if you use the parallel backend (via doParallel), calls mclapply.
# https://github.com/aloysius-lim/bigrf/issues/15
memory = as.numeric(args[1])
datatype = args[2] # 'TPM' (will be log transformed) or deseq 'VSC'
filter = as.logical(args[3])
counts_raw = read.table(args[4], header=T, sep='\t', stringsAsFactors=F, row.names=1)
outdir = args[5]

#counts_raw[is.na(counts_raw)] = 0
#counts_filtered=counts_raw[which(rowSums(counts_raw)!=0),]
#q = c()
#for (i in 1:ncol(counts_filtered)){
#  q = c(q, unname(quantile(counts_filtered[,i][counts_filtered[,i]>0],c(0.05))) )
#}
#counts=log2(counts_filtered+1)

if (filter) {
    counts = filterGenes(counts_raw, filterTypes='central', filterCentralType='median', filterCentralPercentile = 0.3)
} else {
    counts = counts_raw
}
if (datatype == "TPM") counts = log2(counts+1)
ggplot(counts[1], aes(x=counts[,1])) + 
  geom_histogram(binwidth = 0.05, aes(fill=..count..)) + 
  theme_bw() + 
  theme(legend.box = "horizontal", plot.title = element_text(hjust = 0.5)) + 
  ggtitle(paste('log2',datatype,"Histogram",sep=" ")) + 
  xlab(datatype) + 
  ylab('#Genes')
suppressMessages(ggsave(file.path(outdir, "wgcna.histogram.pdf")))

# enableWGCNAThreads(threads)

sftThreshold = pickSoftThreshold(
  t(counts),
  dataIsExpr = TRUE,
  weights = NULL,
  RsquaredCut = 0.85,
  powerVector = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  removeFirst = FALSE,
  nBreaks = 10,
  blockSize = NULL,
  corFnc = bicor,
  corOptions = list(use = 'p'),
  networkType = "signed",
  moreNetworkConcepts = FALSE,
  gcInterval = NULL,
  verbose = 0,
  indent = 0
)

pdf(file.path(outdir, "wgcna.scale.pdf"))
plot(sftThreshold$fitIndices[,1], -sign(sftThreshold$fitIndices[,3])*sftThreshold$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'))
text(sftThreshold$fitIndices[,1], -sign(sftThreshold$fitIndices[,3])*sftThreshold$fitIndices[,2],
     labels=c(seq(1, 10, by = 1), seq(12, 30, by = 2)),cex=1,col='red'); abline(h=0.90,col='red')
graphics.off()

x = -sign(sftThreshold$fitIndices[,3])*sftThreshold$fitIndices[,2]
power = which(x>=0.9)[1]
power = ifelse(!is.na(power),power,ifelse(ncol(counts)<20,18,ifelse(ncol(counts)<31,16,ifelse(ncol(counts)<41,14,12)))) 

#~10k blocksize returns different results but is 10 times faster (blocksize -> parallel blocks is memory dependent)
wgcnar = blockwiseModules(
  t(counts),
#  nThreads = threads,
  power = power,
  maxPOutliers = 0.1,
  corType = "bicor",
  TOMType = "signed",
  networkType = "signed",
  replaceMissingAdjacencies = T,
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 3,
  saveTOMs = F,
  maxBlockSize = blockSize(50000, rectangularBlocks = TRUE, maxMemoryAllocation = 2^31*memory)
)

pdf(file.path(outdir, "wgcna.dendrogram.pdf"))
plotDendroAndColors(main="Gene dendrogram and modules",wgcnar$dendrograms[[1]], labels2colors(wgcnar$unmergedColors)[wgcnar$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
graphics.off()

save(wgcnar, file=file.path(outdir, "wgcnar.Rdata"))

write.table(as.data.frame(wgcnar$colors, row.names=row.names(counts)), file=file.path(outdir, "wgcna.cluster.tsv"), quote=FALSE, sep='\t', col.names = F)
write.table(as.data.frame(wgcnar$unmergedColors, row.names=row.names(counts)), file=file.path(outdir, "wgcna.modules.tsv"), quote=FALSE, sep='\t', col.names = F)

sink(file.path(outdir, "wgcna.cluster2modules"))
for (i in 0:max(wgcnar$colors)){
  cat(i,unique(sort(wgcnar$unmergedColors[which(wgcnar$colors %in% i)])),"\n")
}
sink()
