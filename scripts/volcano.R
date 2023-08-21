#! /usr/bin/env Rscript
# (c) Konstantin Riege

args = commandArgs(TRUE)

if(length(args)<2){
  cat("volcano plot from annotated DESeq2 results table\n")
  cat("\n")
  cat("usage parameter: <f:matrix> <f:outfile>\n")
  cat('example: "/path/to/matrix.tsv" "/path/to/plot.pdf"\n')
  cat("\n")
  cat("matrix: tab separated with header containing at least the columns: id log2FoldChange padj name\n")
  cat("id baseMean log2FoldChange lfcSE stat pvalue padj name\n")
  cat("feature1 value1 value2 value3  ..\n")
  cat("..\n")
  quit("no",1)
}

options(warn=-1)

suppressMessages({
  library(ggplot2)
  library(ggrepel)
})

ddsr = args[1]
outfile = args[2]

df = read.table(ddsr, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="")

df$Regulation = rep("NA",nrow(df))
df$Regulation[!is.na(df$log2FoldChange) & df$log2FoldChange<0 & !is.na(df$padj) & df$padj<=0.05] = "Down"
df$Regulation[!is.na(df$log2FoldChange) & df$log2FoldChange>0 & !is.na(df$padj) & df$padj<=0.05] = "Up"
df$label = rep(NA,nrow(df))

mycolors=c()
if (sum(df$Regulation=="Down")>0){
  mycolors=c(mycolors,"blue")
}
if (sum(df$Regulation=="NA")>0){
  mycolors=c(mycolors,"darkgrey")
}
if (sum(df$Regulation=="Up")>0){
  mycolors=c(mycolors,"red")
}

df = df[rev(order(abs(df$log2FoldChange),na.last = F)),]
topidx = head(which(! is.na(df$padj) & df$padj<0.05),n=10)
if (sum(colnames(df)=="name") == 1) {
  df[topidx,]$label = df[topidx,]$name
} else {
  df[topidx,]$label = df[topidx,]$id
}

df = df[order(abs(df$padj),na.last = T),]
topidx = head(which(! is.na(df$padj) & df$padj<0.05),n=10)
if (sum(colnames(df)=="name") == 1) {
  df[topidx,]$label = df[topidx,]$name
} else {
  df[topidx,]$label = df[topidx,]$id
}

myxlim=max(abs(df$log2FoldChange)+0.5)
ggplot(df, aes(x=log2FoldChange, y=-log10(pvalue), col=Regulation, label=label)) +
  geom_point(alpha=0.65) + 
  theme_classic() +
  scale_color_manual(values=mycolors) +
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="longdash") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype="longdash") +
  ggtitle("Volcano plot") +
  xlab("log2 FoldChange") +
  ylab("-log10(p-value)") +
  xlim(-myxlim, myxlim) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(min.segment.length = 0, box.padding = 0.5, fill="white", force=1, max.overlaps =Inf, na.rm=T, colour = "black", size=3)
suppressMessages(ggsave(file.path(dirname(ddsr),outfile)))
