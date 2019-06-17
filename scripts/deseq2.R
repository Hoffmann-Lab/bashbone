#! /usr/bin/env Rscript
# (c) Konstantin Riege

suppressMessages(library("DESeq2"));
suppressMessages(library("BiocParallel"));
suppressMessages(library("ggplot2"));
suppressMessages(library("gplots"));
suppressMessages(library("RColorBrewer"));
args <- commandArgs(TRUE);
threads <- as.numeric(args[1]);
incsv <- args[2]; # sample,countfile,condition,replicate[,pairs]
outdir <- args[3];
args <- args[4:length(args)];
ctr <- args[rep(c(T,F),length(args)/2)]; # wt t1 wt t2 t1 t2
treat <- args[rep(c(F,T),length(args)/2)];
BPPARAM <- MulticoreParam(workers = threads);
df <- read.table(incsv, header=T, sep=",", stringsAsFactors=F);

if(is.null(df$pairs) || length(Reduce(intersect, split(df$pairs,df$condition))) == 0){
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=df, directory="", design= ~ condition);
} else {
	ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=df, directory="", design= ~ pairs + condition);
};

dds <- DESeq(ddsHTSeq, parallel = TRUE, BPPARAM = BPPARAM);
save(dds, file = file.path(outdir,"dds.Rdata"));

pdf(file.path(outdir,"dispersion.pdf"));
plotDispEsts(dds);
graphics.off();

rld <- rlog(dds, blind=FALSE);
save(rld, file = file.path(outdir,"rld.Rdata"));

data <- plotPCA(rld, intgroup = c("condition", "replicate"), returnData = T);
percentVar <- round(100 * attr(data, "percentVar"));
pdf(file.path(outdir,"pca.pdf"));
ggplot(data, aes(PC1, PC2, color = condition, group = condition, shape = replicate)) +
	ggtitle(paste("PC1 vs PC2: ", length(rownames(rld)), " genes")) +
	scale_shape_manual(values=1:nrow(df)) +
	coord_fixed() + 
	theme_bw() +
	theme(legend.box = "horizontal") +
	geom_point(size = 3) +
	stat_ellipse() +
	xlab(paste("PC1:",percentVar[1],"% variance",sep=" ")) +
	ylab(paste("PC2:",percentVar[2],"% variance",sep=" "));
graphics.off();

vsd <- varianceStabilizingTransformation(dds, blind=FALSE);
save(vsd, file = file.path(outdir,"vsd.Rdata"));

for (i in 1:length(ctr)){
	odir <- file.path(outdir,paste(ctr[i],"-vs-",treat[i],sep=""));
	dir.create(odir, recursive = T);

	ddsr <- results(dds, contrast=c("condition",treat[i],ctr[i]), parallel = TRUE, BPPARAM = BPPARAM);

	ddsr <- ddsr[order(ddsr$padj) , ];
	write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F, 
		file=file.path(odir,"deseq.full.tsv"), quote=F, sep="\t"
	);

	ddsr <- ddsr[!is.na(ddsr$log2FoldChange) , ];
	ddsr <- ddsr[!is.na(ddsr$padj) , ];
	ddsr <- ddsr[ddsr$baseMean > 0 , ];
	ddsr <- ddsr[rev(order(abs(ddsr$log2FoldChange))) , ];
	write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F, 
		file=file.path(odir,"deseq.noNA.tsv"), quote=F, sep="\t"
	);

	ddsr <- ddsr[ddsr$padj <= 0.05 , ];
	write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F, 
		file=file.path(odir,"deseq.tsv"), quote=F, sep="\t"
	);
	pdf(file.path(odir,"ma_plot.pdf"));
	plotMA(ddsr);
	graphics.off();

	vsdr <- vsd[,vsd$condition %in% c(ctr[i],treat[i])];

	vsc <- as.data.frame(assay(vsdr));
	write.table(data.frame(id=rownames(vsc),vsc), row.names = F, 
		file=file.path(odir,"experiments.vsc"), quote=F, sep="\t"
	);

	colnamesvsc <- colnames(vsc)
	colnames(vsc) <- as.vector(vsdr$condition);
	meanvsc <- t(apply(vsc, 1, function(x) tapply(x, colnames(vsc), mean)));
	meanvsc <- meanvsc[,c(ctr[i],treat[i])];
	colnames(vsc) <- colnamesvsc
	write.table(data.frame(id=rownames(meanvsc),meanvsc), row.names = F, 
		file = file.path(odir,"experiments.mean.vsc"), quote=F, sep="\t"
	);

	color <- colorRampPalette(brewer.pal(9, "GnBu"))(100);
	topids <- rownames(ddsr)[1:min(50,nrow(ddsr))];
	vsc <- vsc[rownames(vsc) %in% topids , ];
	write.table(data.frame(id=rownames(vsc),vsc), row.names = F, 
		file=file.path(odir,"heatmap.vsc"), quote=F, sep="\t"
	);
	postscript(file.path(odir,"heatmap.vsc.ps"));
	heatmap.2(as.matrix(vsc), col = color, Rowv = T, Colv = F, scale = "none",
		dendrogram = "row", trace = "none", margin = c(8, 8),
		cexRow = 0.6, cexCol = 1, key.title = NA, key.ylab = NA,
		key.xlab = "vsc", main = "Top 50 FC"
	);
	graphics.off();
	meanvsc <- meanvsc[rownames(meanvsc) %in% topids , ];
	write.table(data.frame(id=rownames(meanvsc),meanvsc), row.names = F, 
		file = file.path(odir,"heatmap.mean.vsc"), quote=F, sep="\t"
	);
	postscript(file.path(odir,"heatmap.mean.vsc.ps"));
	heatmap.2(as.matrix(meanvsc), col = color, Rowv = T, Colv = F, scale = "none",
		dendrogram = "row", trace = "none", margin = c(8, 8),
		cexRow = 0.6, cexCol = 1, key.title = NA, key.ylab = NA,
		key.xlab = "vsc", main = "Top 50 FC"
	);
	graphics.off();
};

quit()

pdf("pca_top500.pdf")
topN <- 500
Pvars <- rowVars(assay(rld))
select <- order(Pvars, decreasing = TRUE)[seq_len(min(topN, length(Pvars)))]
PCA <- prcomp(t(assay(rld)[select, ]), scale = F)
percentVar <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
data <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2], PC3 = PCA$x[,3], PC4 = PCA$x[,4], sampleNO = colData(rld)$type, condition = colData(rld)$condition)
rownames(data) <- data$sampleNO
ggplot(data, aes(PC1, PC2, color = condition, group = condition, shape = labels)) +
	ggtitle(paste("PC1 vs PC2: top ", topN, " variable genes")) +
	scale_shape_manual(values=1:length(labels)) +
	coord_fixed() + 
	theme_bw() +
	theme(legend.box = "horizontal") +
    geom_point(size = 3) +
    stat_ellipse() +
    xlab(paste0("PC1: ",percentVar[1], "% variance")) +
    ylab(paste0("PC2: ",percentVar[2], "% variance"))
graphics.off()
