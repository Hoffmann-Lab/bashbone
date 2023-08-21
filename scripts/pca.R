#! /usr/bin/env Rscript
# (c) Konstantin Riege

args = commandArgs(TRUE)

if(length(args)<3){
  cat("pca 1vs2, 1vs3, 2vs3 from up to top 10000 most variable features\n")
  cat("\n")
  cat("usage parameter: <f:experiments> <f:matrix> <f:outdir> [<i:ignore_column> ..]\n")
  cat('example: "/path/to/experiments.tsv" "/path/to/matrix.tsv" "/path/to/outdir" 1\n')
  cat("\n")
  cat("matrix: tab separated with header. missing values can be '.' or 'NA'\n")
  cat("sample1 sample2 sample3 ..\n")
  cat("value1  value2  value3  ..\n")
  cat("..\n")
  cat("\n")
  cat("experiments: tab separated WITHOUT header. column 3 ignored. color by condition (column 2). dot shape by replicate (column 4) or factor (optional).\n")
  cat("sample1 condition1 foo N1 factorA\n")
  cat("sample2 condition1 foo N2 factorB\n")
  cat("sample3 condition2 foo N1 factorA\n")
  cat("..\n")
  quit("no",1)
}

cat("about to run pca\n")

options(warn=-1)

suppressMessages({
  library(stats)
  library(ggplot2)
  library(ggrepel)
})

# cmp file with rownames == colnames of args[2] i.e. sample condition foo replicate [factor..]
experiments = read.table(args[1], sep="\t", quote="", stringsAsFactors=F, header=F, check.names=F)
# matrix of normed values
df = read.table(args[2], sep="\t", quote="", stringsAsFactors=F, header=T, na.strings=c(".","na","NA","Na","nA"), check.names = F)
outdir = args[3]

if(ncol(experiments)==4){
  colnames(experiments) = c("sample","condition","foo","replicate")
  experiments$factor=experiments$replicate
} else {
  colnames(experiments) = c("sample","condition","foo","replicate","factor")
}

if(length(args)>3){
  # df = df[, ! colnames(df) %in% args[4:length(args)]]
  df = df[,-1*as.integer(args[4:length(args)])]
}
df = df[rowSums(is.na(df))==0,]

vars = apply(df, 1, var)
df = head(df[rev(order(vars)),],n=10000)
pca = prcomp(t(df), scale = F)

loadings = pca$rotation
for (i in 1:3){
  ids = rownames(loadings[order(abs(loadings[,i]), decreasing = TRUE),])
  ids = head(ids,n=max(1,length(ids)*0.05)) # variables that drive variation in PC1
  sink(file.path(outdir,paste("pca_pc",i,".variables",sep="")))
  lapply(ids, cat, "\n")
  sink()
}

percentVar = round(100*pca$sdev^2/sum(pca$sdev^2),1)
data = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], sample = rownames(pca$x))
experiments$sample = paste(experiments$condition, experiments$replicate, sep='.')
data = merge(data,experiments,by = "sample")

ggplot(data, aes(PC1, PC2, color = condition, group = condition, shape = factor, label=sample)) +
  ggtitle("PCA plot - PC1 vs PC2") +
  scale_shape_manual(values = c(1:length(unique(data$factor)) )) +
  # coord_fixed() +
  theme_bw() +
  theme(aspect.ratio=1, legend.box = "horizontal", legend.title=element_blank()) +
  geom_point(size = 3) +
  # geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1], "% variance")) +
  ylab(paste0("PC2: ",percentVar[2], "% variance"))
suppressMessages(ggsave(file.path(outdir,"pca_12.pdf")))

ggplot(data, aes(PC1, PC3, color = condition, group = condition, shape = factor, label=sample)) +
  ggtitle("PCA plot - PC1 vs PC3") +
  scale_shape_manual(values = c(1:length(unique(data$factor)) )) +
  # coord_fixed() +
  theme_bw() +
  theme(aspect.ratio=1, legend.box = "horizontal", legend.title=element_blank()) +
  geom_point(size = 3) +
  # geom_text_repel() +
  xlab(paste0("PC1: ",percentVar[1], "% variance")) +
  ylab(paste0("PC3: ",percentVar[3], "% variance"))
suppressMessages(ggsave(file.path(outdir,"pca_13.pdf")))

ggplot(data, aes(PC2, PC3, color = condition, group = condition, shape = factor, label=sample)) +
  ggtitle("PCA plot - PC2 vs PC3") +
  scale_shape_manual(values = c(1:length(unique(data$factor)) )) +
  # coord_fixed() +
  theme_bw() +
  theme(aspect.ratio=1, legend.box = "horizontal", legend.title=element_blank()) +
  geom_point(size = 3) +
  # geom_text_repel() +
  xlab(paste0("PC2: ",percentVar[2], "% variance")) +
  ylab(paste0("PC3: ",percentVar[3], "% variance"))
suppressMessages(ggsave(file.path(outdir,"pca_23.pdf")))
