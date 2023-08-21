#! /usr/bin/env Rscript
# (c) Konstantin Riege

args = commandArgs(TRUE)

if(length(args)<6){
  cat("WGCNA from normalized feature counts\n")
  cat("\n")
  cat("usage parameter: <i:memory-mb> <b:log-transform> <b:filter-percentile> <f:matrix> <f:outdir> [<f:experiments>]\n")
  cat('example: 16000 TPM FALSE "/path/to/matrix.tsv" "/path/to/outdir"\n')
  cat("\n")
  cat("matrix: tab separated with header and feature ids/label.\n")
  cat("id       sample1 sample2 sample3 ..\n")
  cat("feature1 value1  value2  value3  ..\n")
  cat("..\n")
  cat("\n")
  cat("experiments (optional): comma separated WITH header. 4 or more columns. column 2 ignored. cluster-condition relationship merged by replicate (column 4).\n")
  cat("sample,path,condition,replicate[,..]\n")
  cat("sample1,foo,condition1,N1[,..]\n")
  cat("sample2,foo,condition1,N2[,..]\n")
  cat("sample3,foo,condition2,N1[,..]\n")
  cat("..\n")
  quit("no",1)
}

options(warn=-1)

suppressMessages({
  library(WGCNA)
  library(DGCA)
  library(ggplot2)
  library(dendextend)
})

# threads = as.numeric(args[1])
# may cause Error in sendMaster(try(lapply(X = S, FUN = FUN, ...), silent = TRUE))
# The error seems to be from sendMaster which is used by the multicore package to send results from child processes back to the master process. From this StackOverflow thread, it seems that the results of the child processes are too large for mclapply to send back to the master.
# Dev uses foreach, which, if you use the parallel backend (via doParallel), calls mclapply.
# https://github.com/aloysius-lim/bigrf/issues/15
memory = as.numeric(args[1])
dolog = as.logical(args[2]) # e.g. TPM will be log transformed, whereas VSC or RLD from deseq are log-scaled
filter = as.logical(args[3])
counts_raw = read.table(args[4], header=T, sep='\t', row.names=1, stringsAsFactors=F, check.names=F, quote="")
outdir = args[5]
experiments = data.frame()
if(length(args)>5){
  experiments = read.table(args[6], header=T, sep=',', stringsAsFactors=F, check.names=F, quote="")
}
#counts_raw[is.na(counts_raw)] = 0
#counts_filtered=counts_raw[which(rowSums(counts_raw)!=0),]
#q = c()
#for (i in 1:ncol(counts_filtered)){
#  q = c(q, unname(quantile(counts_filtered[,i][counts_filtered[,i]>0],c(0.05))) )
#}
#counts=log2(counts_filtered+1)

if (filter) {
  counts = filterGenes(counts_raw, filterTypes='central', filterCentralType='median', filterCentralPercentile = 0.3)
  # alternative filter by variance
  # vars = apply(df, 1, var)
  # counts = counts_raw[vars > (max(vars) - min(vars))/3 ,]
} else {
  counts = counts_raw
}
if (dolog) counts = log2(counts+1)
# ggplot(counts[1], aes(x=counts[,1])) +
#   geom_histogram(binwidth = 0.05, aes(fill=..count..)) +
#   theme_bw() +
#   theme(legend.box = "horizontal", plot.title = element_text(hjust = 0.5)) +
#   ggtitle(paste('log2',datatype,"Histogram",sep=" ")) +
#   xlab(datatype) +
#   ylab('#Genes')
# suppressMessages(ggsave(file.path(outdir, "wgcna.histogram.pdf")))


# enableWGCNAThreads(threads)
sftThreshold = pickSoftThreshold(
  t(counts),
  dataIsExpr = TRUE,
  weights = NULL,
  #RsquaredCut = 0.85,
  RsquaredCut = 0.9,
  # powerVector = c(seq(1, 10, by = 1), seq(12, 30, by = 2)),
  powerVector = c(1:30),
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

# power = sftThreshold$powerEstimate # better re-compute to handle NA
x = -sign(sftThreshold$fitIndices[,3])*sftThreshold$fitIndices[,2]
# power = which(x>0.9)[1] # does not work for non-continuous powerVector. instead use the power column from sftThreshold object
power = sftThreshold$fitIndices[,1][which(x>0.9)[1]]
# NA handler: fallback according to manual
power = ifelse(!is.na(power),power,ifelse(ncol(counts)<20,18,ifelse(ncol(counts)<31,16,ifelse(ncol(counts)<41,14,12))))
cat(paste0("scale free topology index reaches 0.9 at power: ",power,"\n"))

pdf(file.path(outdir, "wgcna.scale.pdf"))
plot(sftThreshold$fitIndices[,1], -sign(sftThreshold$fitIndices[,3])*sftThreshold$fitIndices[,2],
  xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit, signed R^2',
  type='n', main = paste('Scale independence'))
text(sftThreshold$fitIndices[,1], -sign(sftThreshold$fitIndices[,3])*sftThreshold$fitIndices[,2],
  # labels=c(seq(1, 10, by = 1), seq(12, 30, by = 2)),cex=1,col='red')
  labels=c(1:30),cex=1,col='red')
#abline(h=0.85,col='red')
abline(h=0.90,col='red')
graphics.off()


wgcnar = blockwiseModules(
  t(counts),
#  nThreads = threads,
  randomSeed = 12345,
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
  saveTOMs = T,
  saveTOMFileBase = file.path(outdir, "wgcnar.TOM"),
  maxBlockSize = blockSize(50000, rectangularBlocks = TRUE, maxMemoryAllocation = 2^31*memory)
)
#~10k blocksize returns different results but is 10 times faster (blocksize -> parallel blocks is memory dependent)
save(wgcnar, file=file.path(outdir, "wgcnar.Rdata"))


###########################
# save unmerged colors and merged colors as modules and clusters map
sink(file.path(outdir, "wgcna.cluster2modules"))
for (i in 0:max(wgcnar$colors)){
  cat(i,unique(sort(wgcnar$unmergedColors[which(wgcnar$colors %in% i)])),"\n")
}
sink()

###########################
# re-calculate eigengenes (first principal component within a module) with correct color codes (from merged modules aka cluster) and order by similarity
MEs = orderMEs(moduleEigengenes(t(counts), colors=labels2colors(wgcnar$colors))$eigengenes)
colnames(MEs) = unlist(lapply(colnames(MEs),function(x) sub("ME", "", x)))
# convert colors to cluster ids
MEcolnames = match(colnames(MEs),labels2colors(0:(ncol(MEs)-1)))-1

# order by ids ???
# MEs = MEs[,order(MEcolnames)]
# MEcolnames = 0:(length(MEcolnames)-1)

MEcolnames = sprintf("Cluster.%03d",MEcolnames)

# measure correlations between each gene and each ME (intramodular connectivity values kME)
MEcor = signedKME(t(counts), MEs)
colnames(MEcor) = unlist(lapply(MEcolnames,paste0,".cor"))
# get p-values for correlations
dfp = corPvalueStudent(as.matrix(MEcor),ncol(counts))
colnames(dfp) = unlist(lapply(MEcolnames,paste0,".pvalue"))
dfq = apply(dfp, 2, p.adjust, method="BH")
colnames(dfq) = unlist(lapply(MEcolnames,paste0,".padj"))
MEcor = do.call(cbind,lapply(1:ncol(MEs), function(x) cbind( subset(MEcor,select=x),subset(dfp,select=x),subset(dfq,select=x) )))
MEcor = cbind(as.data.frame(wgcnar$colors),MEcor)
colnames(MEcor)[1] = "cluster"
MEcor$cluster = sprintf("Cluster.%03d",MEcor$cluster)

writeLines(paste(c("id",colnames(MEcor)),collapse="\t") ,con = file.path(outdir, "wgcna.cluster.full.tsv"))
writeLines(paste(c("id",colnames(MEcor)),collapse="\t") ,con = file.path(outdir, "wgcna.cluster.tsv"))
writeLines(paste(c("id",colnames(MEcor)),collapse="\t") ,con = file.path(outdir, "wgcna.cluster.top.tsv"))
for (i in MEcolnames){
  # order genes of a cluster by correlation values or padj
  df = subset(MEcor[MEcor$cluster==i,],select=c(paste0(i,".cor"),paste0(i,".padj")))
  df = df[order(abs(df[,1]),decreasing = T),]
  write.table(MEcor[match( rownames(df), rownames(MEcor) ) ,], file=file.path(outdir, "wgcna.cluster.full.tsv"), quote=FALSE, sep='\t', col.names = F, append=T)
  df = df[df[,2]<=0.05,]
  write.table(MEcor[match( rownames(df), rownames(MEcor) ) ,], file=file.path(outdir, "wgcna.cluster.tsv"), quote=FALSE, sep='\t', col.names = F, append=T)
  df = head(df,n=50)
  write.table(MEcor[match( rownames(df), rownames(MEcor) ) ,], file=file.path(outdir, "wgcna.cluster.top.tsv"), quote=FALSE, sep='\t', col.names = F, append=T)
}

###########################
# plot ME hclust tree and pca
distPC1 = 1-abs(cor(MEs, use="p"))
distPC1 = ifelse(is.na(distPC1), 0, distPC1)
treePC1 = as.dendrogram(hclust(as.dist(distPC1), method="a"))
labels_colors(treePC1) = labels(treePC1)
pdf(file.path(outdir, "wgcna.cluster.tree.pdf"))
plot(treePC1, xlab="",ylab="",main="",sub="")
graphics.off()
mds = cmdscale(as.dist(distPC1),2)
pdf(file.path(outdir, "wgcna.cluster.PCA.pdf"))
plot(mds, col = names(MEs),cex=2, pch=19, xlab="PC1", ylab="PC2", main="Module eigengenes PCA")
graphics.off()

###########################
# compute module/cluster-traits relationships and plot as heatmap
if(ncol(experiments)<4){
  # traits info can be categorial or quantitative
  #         condition1 condition2 signal
  # sample1     1          0        10.1
  # sample2     1          0        15.5
  # sample3     0          1        40.3
  # sample4     0          1        73.0
  traits = matrix(0,nrow=ncol(counts),ncol=ncol(counts))
  diag(traits) = 1
  traits = as.data.frame(traits)
  colnames(traits) = colnames(counts)
  rownames(traits) = colnames(counts)

  Tcor = cor(MEs, traits, use= "p")
  dfp = corPvalueStudent(Tcor, ncol(counts))
  dfq = apply(dfp, 2, p.adjust, method="BH")
  rownames(dfq) = rownames(dfp)

  txt = paste(signif(Tcor, 2), "\n(",signif(dfq, 1), ")", sep= "")
  dim(txt) = dim(Tcor)
  pdf(file.path(outdir,"wgcna.cluster.heatmap.pdf"))
  labeledHeatmap(
    Matrix= Tcor,
    xLabels= colnames(traits),
    #yLabels= paste0("ME",colnames(MEs)),
    #ySymbols= sprintf("Cluster.%03d",MEcolnames),
    #colorLabels= T,
    yLabels= MEcolnames,
    colors= blueWhiteRed(50),
    textMatrix= txt,
    setStdMargins= FALSE,
    cex.text= 0.5,
    zlim= c(-1,1),
    main= paste("Cluster-condition correlation with FDR")
  )
  graphics.off()

} else {

  traits = matrix(0,nrow=nrow(experiments),ncol=nrow(experiments))
  diag(traits) = 1
  traits = as.data.frame(traits)
  colnames(traits) = paste(experiments$condition,experiments$replicate,sep=".")
  rownames(traits) = paste(experiments$condition,experiments$replicate,sep=".")

  Tcor = cor(MEs, traits, use= "p")
  dfp = corPvalueStudent(Tcor, ncol(counts))
  dfq = apply(dfp, 2, p.adjust, method="BH")
  rownames(dfq) = rownames(dfp)

  txt = paste(signif(Tcor, 2), "\n(",signif(dfq, 1), ")", sep= "")
  dim(txt) = dim(Tcor)

  pdf(file.path(outdir,"wgcna.cluster.heatmap.pdf"),width=max(7,as.integer(ncol(Tcor)/3)),height=max(7,as.integer(nrow(Tcor)/2)))
  par(mar= c(5, 8, 3, 3))
  labeledHeatmap(
    Matrix= Tcor,
    xLabels= colnames(traits),
    #yLabels= paste0("ME",colnames(MEs)),
    #ySymbols= sprintf("Cluster.%03d",MEcolnames),
    #colorLabels= T,
    yLabels= MEcolnames,
    colors= blueWhiteRed(50),
    textMatrix= txt,
    setStdMargins= FALSE,
    cex.text= 0.5,
    zlim= c(-1,1),
    main= paste("Cluster-condition correlation with FDR")
  )
  graphics.off()

  traits = unlist(lapply(unique(experiments$condition),function (x) as.integer(experiments$condition==x)))
  traits = as.data.frame(matrix(traits,nrow=nrow(experiments)))
  colnames(traits) = unique(experiments$condition)
  rownames(traits) = paste(experiments$condition,experiments$replicate,sep=".")

  Tcor = cor(MEs, traits, use= "p")
  dfp = corPvalueStudent(Tcor, ncol(counts))
  dfq = apply(dfp, 2, p.adjust, method="BH")
  rownames(dfq) = rownames(dfp)

  txt = paste(signif(Tcor, 2), "\n(",signif(dfq, 1), ")", sep= "")
  dim(txt) = dim(Tcor)
  pdf(file.path(outdir,"wgcna.cluster.mean.heatmap.pdf"),width=max(7,as.integer(ncol(Tcor)/3)),height=max(7,as.integer(nrow(Tcor)/2)))
  par(mar= c(5, 8, 3, 3))
  labeledHeatmap(
    Matrix= Tcor,
    xLabels= colnames(traits),
    #yLabels= paste0("ME",colnames(MEs)),
    #ySymbols= sprintf("Cluster.%03d",MEcolnames),
    #colorLabels= T,
    yLabels= MEcolnames,
    colors= blueWhiteRed(50),
    textMatrix= txt,
    setStdMargins= FALSE,
    cex.text= 0.5,
    zlim= c(-1,1),
    main= paste("Cluster-condition correlation with FDR")
  )
  graphics.off()
}

###########################
# plot gene tree (prior merge dendrograms in case of multiple blocks due to limited memory for large datasets)
if(length(wgcnar$blockGenes)==1){
  dendro = wgcnar$dendrograms[[1]]
  colors = labels2colors(wgcnar$colors)
} else {
  # requires ulimit -s $(ulimit -Hs)
  dendro = as.hclust(Reduce(function(x, y) merge(x, y, height = 1.0001), lapply(wgcnar$dendrograms, function(x) as.dendrogram(x))))
  colors = labels2colors(wgcnar$colors)[unlist(wgcnar$blockGenes)]
}
pdf(file.path(outdir,"wgcna.dendrogram.pdf"))
plotDendroAndColors(main="Dendrogram",dendro, colors,"Cluster colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
graphics.off()

###########################
# plot topological overlap matrix as heatmap (prior subsample clusters and re-compute tom in case of multiple blocks due to limited memory for large datasets)
if(nrow(counts)>1000){
  n=1000/nrow(counts)
  set.seed(12345)
  select = unlist(lapply(unique(wgcnar$colors), function(x) (sample(which(wgcnar$colors==x), size = max(1,sum(wgcnar$colors==x)*n)) )))
} else {
  select = 1:nrow(counts)
}
if(length(wgcnar$blockGenes)==1){
  load(wgcnar$TOMFiles[[1]])
  dissTOM = 1-as.matrix(TOM)[select,select]
} else {
  dissTOM = 1-TOMsimilarityFromExpr(
    t(counts[select,]), power = power,
    corType = "bicor", TOMType = "signed",
    maxPOutliers = 0.1, networkType = "signed",
    replaceMissingAdjacencies = T
  )
}
diag(dissTOM) = NA
dendro = hclust(as.dist(dissTOM), method="average")
pdf(file.path(outdir,"wgcna.tom.heatmap.pdf"))
TOMplot(1-dissTOM^4, dendro, labels2colors(wgcnar$colors)[select], main = paste0("TOM heatmap plot, n = ",length(select)))
graphics.off()


###########################
# adjacency=WGCNA::bicor(t(counts_raw)) with rownames=colanmes == rownames(counts_raw), call adjacency["name1/id1","name2/id2"] to check if pos or neg correlated [-1,+1]
# TOM = bicor ^ power * tomtype formula. use TOM to check for correlation of two disjoint gene sets e.g. pathways
# sample a normal distribution from mean values of n rounds and calculate p-values for both sets
# load(/path/to/tom/); adjacency_matrix=as.matrix(TOM) set rownames=colanmes == rownames(counts_raw)
bipartiteConnectivity<-function(adjacency_matrix,a,b,rounds=100000){
  c=a
  a=setdiff(a,b)
  b=setdiff(b,c)
  x=sapply(1:rounds,function(x){
    sample_a=sample(nrow(adjacency_matrix),length(a))
    sample_b=sample(setdiff(1:nrow(adjacency_matrix),sample_a),length(b))
    sum(adjacency_matrix[sample_a,sample_b])/(length(a)*length(b))
  })
  stat=sum(adjacency_matrix[a,b])/(length(a)*length(b))
  stat2=abs(stat-median(x))
  x2=abs(x-median(x))
  return(c(stat=stat,pvalue_one_sided=length(which(x>=stat))/length(x),mean=mean(x),median=median(x),stat2=stat2,pvalue_two_sided=length(which(x2>=stat2))/length(x)))
}