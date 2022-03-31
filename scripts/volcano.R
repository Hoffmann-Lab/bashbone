#! /usr/bin/env Rscript
# (c) Konstantin Riege
options(warn=-1)

suppressMessages({
  library(ggplot2)
  library(ggrepel)
})

args = commandArgs(TRUE)
ddsr = args[1]
outfile = args[2]
  
df = read.table(ddsr, header=T, quote="", stringsAsFactors = F, check.names = F, sep = "\t")

df$Regulation = rep("NA",nrow(df))
df$Regulation[!is.na(df$log2FoldChange) & df$log2FoldChange<0 & !is.na(df$padj) & df$padj<=0.05] = "Down"
df$Regulation[!is.na(df$log2FoldChange) & df$log2FoldChange>0 & !is.na(df$padj) & df$padj<=0.05] = "Up"
df$label = rep(NA,nrow(df))

df = df[rev(order(abs(df$log2FoldChange),na.last = F)),]
topidx = head(which(! is.na(df$padj) & df$padj<0.05),n=10)
if (sum(colnames(df)=="geneName") == 1) {
  df[topidx,]$label = df[topidx,]$geneName
} else {
  df[topidx,]$label = df[topidx,]$id
}

df = df[order(abs(df$padj),na.last = T),]
topidx = head(which(! is.na(df$padj) & df$padj<0.05),n=10)
if (sum(colnames(df)=="geneName") == 1) {
  df[topidx,]$label = df[topidx,]$geneName
} else {
  df[topidx,]$label = df[topidx,]$id
}

ggplot(df, aes(x=log2FoldChange, y=-log10(pvalue), col=Regulation, label=label)) +
  geom_point(alpha=0.65) + 
  theme_classic() +
  scale_color_manual(values=c("blue", "darkgrey", "red")) +
  geom_vline(xintercept=c(-0.5, 0.5), col="black", linetype="longdash") +
  geom_hline(yintercept=-log10(0.05), col="black", linetype="longdash") +
  ggtitle("Volcano plot") +
  xlab("log2 FoldChange") +
  ylab("-log10(p-value)") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(min.segment.length = 0, box.padding = 0.5, fill="white", force=1, max.overlaps =Inf, na.rm=T, colour = "black", size=3)
suppressMessages(ggsave(file.path(dirname(ddsr),outfile)))
