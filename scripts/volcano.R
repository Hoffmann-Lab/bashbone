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
if (is.null(df$name)){
  df$name=df$id
}
df$name[is.na(df$name)] = df$id[is.na(df$name)]
df$name[df$name==""] = df$id[df$name==""]

df$Regulation = rep("n.s.",nrow(df))
df$Regulation[!is.na(df$log2FoldChange) & df$log2FoldChange<0 & !is.na(df$padj) & df$padj<=0.05] = "Down"
df$Regulation[!is.na(df$log2FoldChange) & df$log2FoldChange>0 & !is.na(df$padj) & df$padj<=0.05] = "Up"
df$label = rep(NA,nrow(df))

e = df[1,]
e$id = "_fake_down_"
e$name = e$id
e$log2FoldChange = 0
e$padj = 1
e$pvalue = 10
e$Regulation="Down"
df = rbind(e,df)
e$id = "_fake_up_"
e$name = e$id
e$log2FoldChange = 0
e$padj = 1
e$pvalue = 10
e$Regulation="Up"
df = rbind(e,df)

mycolors=c()
if (sum(df$Regulation=="Down")>0){
  mycolors=c(mycolors,"blue")
}
if (sum(df$Regulation=="n.s.")>0){
  mycolors=c(mycolors,"darkgrey")
}
if (sum(df$Regulation=="Up")>0){
  mycolors=c(mycolors,"red")
}

nups=10
ndowns=10
ups=sum(!is.na(df$padj) & df$padj<0.05 & df$log2FoldChange > 0)
downs=sum(!is.na(df$padj) & df$padj<0.05 & df$log2FoldChange < 0)
if (ups < downs){
  if (downs > 0 && ups/downs*100 < 20) nups=5
} else {
  if (ups > 0 && downs/ups*100 < 20) ndowns=5
}

df = df[order(df$log2FoldChange),]
topidx = head(which(! is.na(df$padj) & df$padj<0.05),n=ndowns)
df[topidx,]$label = df[topidx,]$name
topidx = tail(which(! is.na(df$padj) & df$padj<0.05),n=nups)
df[topidx,]$label = df[topidx,]$name
df = df[order(df$padj,na.last = T),]
topidx = head(which(! is.na(df$padj) & df$padj<0.05),n=10)
df[topidx,]$label = df[topidx,]$name

df = df[! is.na(df$padj),]
myxlim=max(abs(df$log2FoldChange))+0.5
myylim=max(0.5+-log10(df$pvalue))
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
  annotate("text", x=myxlim, y=myylim, color="red", label= paste0(ups,paste0(rep(" ",10-nchar(ups)),collapse="") )) +
  annotate("text", x=-myxlim, y=myylim, color="blue", label= paste0(paste0(rep(" ",10-nchar(downs)),collapse=""),downs)) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  geom_label_repel(min.segment.length = 0, box.padding = 0.5, fill="white", force=1, max.overlaps =Inf, na.rm=T, colour = "black", size=3) +
  scale_y_continuous(limits = c(0,myylim),breaks = scales::pretty_breaks(n = round(myylim/2))) +
  scale_x_continuous(limits = c(-myxlim,myxlim),breaks = scales::pretty_breaks(n = round(myxlim/2)*2))
suppressMessages(ggsave(file.path(dirname(ddsr),outfile)))
