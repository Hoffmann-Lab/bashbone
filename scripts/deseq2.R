#! /usr/bin/env Rscript
# (c) Konstantin Riege
options(warn=-1)

cat("about to run deseq2 and pca\n")

suppressMessages(library("DESeq2"))
suppressMessages(library("BiocParallel"))
suppressMessages(library("ggplot2"))
suppressMessages(library("gplots"))
suppressMessages(library("RColorBrewer"))

args = commandArgs(TRUE)
threads = as.numeric(args[1])
incsv = args[2] # sample,countfile,condition,replicate[,factor1,factor2,..]
outdir = args[3]
args = args[4:length(args)]
ctr = args[rep(c(T,F),length(args)/2)] # wt t1 wt t2 t1 t2
treat = args[rep(c(F,T),length(args)/2)]

dir.create(outdir, recursive = T, showWarnings = F)
BPPARAM = MulticoreParam(workers = threads)

df = read.table(incsv, header=T, sep=",", stringsAsFactors=F)
colnames(df)[1:4] = c("sample","countfile","condition","replicate")

# e.g. ~ factor1 + factor2 + condition
get_design = function(df){
	factors = c()
	if(length(colnames(df)) > 4){
		for (f in colnames(df)[5:length(colnames(df))]){
			v = df[,f]
			# levels(xyxy) -> 2 > 1 && xyxy @ aabb -> a:xy , b:xy -> {x,y} > 0 i.e. no linear combination
			# whereas levels(xxxx) -> 1 and xyzz @ aabb -> a:xy , b:z -> {} == 0 i.e. deseq error
			if (length(levels(as.factor(v))) > 1 && length(Reduce(intersect, split(v,df$condition))) > 0){
				factors = c(factors,f)
			}
		}
	}
	return(paste("~",paste(c(factors,"condition"),collapse=" + ")))
}

design = get_design(df)
cat(paste("computing pca with deseq normalized read counts based on design formula: ",design,"\n",sep=""))

ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = df, directory = "", design = as.formula(design))
dds = DESeq(ddsHTSeq, parallel = TRUE, BPPARAM = BPPARAM)
save(dds, file = file.path(outdir,"dds.Rdata"))

pdf(file.path(outdir,"dispersion.pdf"))
plotDispEsts(dds)
graphics.off()

###### pca

log = DESeqTransform(SummarizedExperiment(log2(counts(dds, normalized=T) + 1), colData=colData(dds)))
save(log, file = file.path(outdir,"log.Rdata"))
vsd = varianceStabilizingTransformation(dds, blind=FALSE)
save(vsd, file = file.path(outdir,"vsd.Rdata"))
rld = rlog(dds, blind=FALSE);
save(rld, file = file.path(outdir,"rld.Rdata"))

for (method in c("log","vsd","rld")){
	# deseq method
	# data = plotPCA(get(normed), intgroup = c("condition", "replicate"), returnData = T);
	# percentVar = round(100 * attr(data, "percentVar"));

	normed = assay(get(method))
	vars = rowVars(normed)
	n = 100
	#n = length(vars)
	topidx = order(vars, decreasing = TRUE)[1:min(n,length(vars))]
	pca = prcomp(t(normed[topidx, ]), scale = F)
	percentVar = round(100*pca$sdev^2/sum(pca$sdev^2),1)
	data = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], PC4 = pca$x[,4], replicate = df$replicate, condition = df$condition)
	write.table(data.frame(id=rownames(data),data), row.names = F,
		file=file.path(outdir,paste("pca_12_",method,".tsv",sep="")), quote=F, sep="\t"
	)

	suppressMessages({
		ggplot(data, aes(PC1, PC2, color = condition, group = condition, shape = replicate)) +
			ggtitle("PCA plot - PC1 vs PC2") +
			scale_shape_manual(values = c(1:length(unique(data$replicate)) )) +
			coord_fixed() +
			theme_bw() +
			theme(legend.box = "horizontal") +
			geom_point(size = 3) +
			stat_ellipse() +
			xlab(paste0("PC1: ",percentVar[1], "% variance")) +
			ylab(paste0("PC2: ",percentVar[2], "% variance"))
		ggsave(file.path(outdir,paste("pca_12_",method,".pdf",sep="")))

		ggplot(data, aes(PC1, PC3, color = condition, group = condition, shape = replicate)) +
			ggtitle("PCA plot - PC1 vs PC3") +
			scale_shape_manual(values = c(1:length(unique(data$replicate)) )) +
			coord_fixed() +
			theme_bw() +
			theme(legend.box = "horizontal") +
			geom_point(size = 3) +
			stat_ellipse() +
			xlab(paste0("PC1: ",percentVar[1], "% variance")) +
			ylab(paste0("PC3: ",percentVar[3], "% variance"))
		ggsave(file.path(outdir,paste("pca_13_",method,".pdf",sep="")))

		ggplot(data, aes(PC2, PC3, color = condition, group = condition, shape = replicate)) +
			ggtitle("PCA plot - PC2 vs PC3") +
			scale_shape_manual(values = c(1:length(unique(data$replicate)) )) +
			coord_fixed() +
			theme_bw() +
			theme(legend.box = "horizontal") +
			geom_point(size = 3) +
			stat_ellipse() +
			xlab(paste0("PC2: ",percentVar[2], "% variance")) +
			ylab(paste0("PC3: ",percentVar[3], "% variance"))
		ggsave(file.path(outdir,paste("pca_23_",method,".pdf",sep="")))
	})
}

###### deseq

get_table = function(dds){
	ddsr = results(dds, contrast=c("condition",treat[i],ctr[i]), parallel = TRUE, BPPARAM = BPPARAM)

	ddsr = ddsr[order(ddsr$padj) , ]
	write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F,
		file=file.path(odir,"deseq.full.tsv"), quote=F, sep="\t"
	);

	ddsr = ddsr[!is.na(ddsr$log2FoldChange) , ]
	ddsr = ddsr[!is.na(ddsr$padj) , ]
	ddsr = ddsr[ddsr$baseMean > 0 , ]
	ddsr = ddsr[rev(order(abs(ddsr$log2FoldChange))) , ];
	write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F,
		file=file.path(odir,"deseq.noNA.tsv"), quote=F, sep="\t"
	)

	ddsr = ddsr[ddsr$padj <= 0.05 , ]
	write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F,
		file=file.path(odir,"deseq.tsv"), quote=F, sep="\t"
	)

	if(nrow(ddsr)>0){
		pdf(file.path(odir,"ma_plot.pdf"))
		plotMA(ddsr)
		graphics.off()
	}

	rldr = rld[,rld$condition %in% c(ctr[i],treat[i])]
	data = plotPCA(rldr, intgroup = c("condition", "replicate"), returnData = T)
	write.table(data.frame(id=rownames(data),data), row.names = F,
		file=file.path(odir,"pca.tsv"), quote=F, sep="\t"
	)
	percentVar = round(100 * attr(data, "percentVar"))

	suppressMessages({
		ggplot(data, aes(PC1, PC2, color = condition, group = condition, shape = replicate)) +
			ggtitle(paste("PC1 vs PC2: ", length(rownames(rldr)), " genes")) +
			scale_shape_manual(values = c(1:length(unique(data$replicate)) )) +
			coord_fixed() +
			theme_bw() +
			theme(legend.box = "horizontal") +
			geom_point(size = 3) +
			stat_ellipse() +
			xlab(paste("PC1:",percentVar[1],"% variance",sep=" ")) +
			ylab(paste("PC2:",percentVar[2],"% variance",sep=" "))
		ggsave(file.path(odir,"pca.pdf"))
	})

	vsdr = vsd[,vsd$condition %in% c(ctr[i],treat[i])];

	vsc = as.data.frame(assay(vsdr))
	write.table(data.frame(id=rownames(vsc),vsc), row.names = F,
		file=file.path(odir,"experiments.vsc"), quote=F, sep="\t"
	);
	zscores = log(vsc+1)
	zscores = zscores-rowMeans(zscores)
	zscores = zscores/apply(zscores,1,sd)
	zscores[is.na(zscores)] = 0
	write.table(data.frame(id=rownames(zscores),zscores), row.names = F,
		file=file.path(odir,"experiments.vsc.zscores"), quote=F, sep="\t"
	)

	colnamesvsc = colnames(vsc)
	colnames(vsc) = as.vector(vsdr$condition)
	meanvsc = t(apply(vsc, 1, function(x) tapply(x, colnames(vsc), mean)))
	meanvsc = meanvsc[,c(ctr[i],treat[i])]
	colnames(vsc) = colnamesvsc
	write.table(data.frame(id=rownames(meanvsc),meanvsc), row.names = F,
		file = file.path(odir,"experiments.mean.vsc"), quote=F, sep="\t"
	)
	meanzscores = log(meanvsc+1)
	meanzscores = meanzscores-rowMeans(meanzscores)
	meanzscores = meanzscores/apply(meanzscores,1,sd)
	meanzscores[is.na(meanzscores)] = 0
	write.table(data.frame(id=rownames(meanzscores),meanzscores), row.names = F,
		file=file.path(odir,"experiments.mean.vsc.zscores"), quote=F, sep="\t"
	)

	if(nrow(ddsr)>1){
		color = colorRampPalette(brewer.pal(9, "GnBu"))(100)
		topids = rownames(ddsr)[1:min(50,nrow(ddsr))]

		vsc = vsc[rownames(vsc) %in% topids , ]
		write.table(data.frame(id=rownames(vsc),vsc), row.names = F,
			file=file.path(odir,"heatmap.vsc"), quote=F, sep="\t"
		)
		postscript(file.path(odir,"heatmap.vsc.ps"))
		heatmap.2(as.matrix(vsc), col = color, Rowv = T, Colv = F, scale = "none",
			dendrogram = "row", trace = "none", margin = c(10, 8),
			cexRow = 0.6, cexCol = 0.8, key.title = NA, key.ylab = NA,
			key.xlab = "vsc", main = paste("VSC Heatmap of top ",length(topids)," genes by fold-change")
		)
		graphics.off()

		zscores = zscores[rownames(zscores) %in% topids , ]
		write.table(data.frame(id=rownames(zscores),zscores), row.names = F,
			file=file.path(odir,"heatmap.vsc.zscores"), quote=F, sep="\t"
		)
		postscript(file.path(odir,"heatmap.vsc.zscores.ps"))
		heatmap.2(as.matrix(zscores), col = color, Rowv = T, Colv = F, scale = "none",
			dendrogram = "row", trace = "none", margin = c(10, 8),
			cexRow = 0.6, cexCol = 0.8, key.title = NA, key.ylab = NA,
			key.xlab = "zscore", main = paste("Z-score heatmap of top ",length(topids)," genes by fold-change")
		)
		graphics.off()

		meanvsc = meanvsc[rownames(meanvsc) %in% topids , ]
		write.table(data.frame(id=rownames(meanvsc),meanvsc), row.names = F,
			file = file.path(odir,"heatmap.mean.vsc"), quote=F, sep="\t"
		)
		postscript(file.path(odir,"heatmap.mean.vsc.ps"))
		heatmap.2(as.matrix(meanvsc), col = color, Rowv = T, Colv = F, scale = "none",
			dendrogram = "row", trace = "none", margin = c(10, 8),
			cexRow = 0.6, cexCol = 0.8, key.title = NA, key.ylab = NA,
			key.xlab = "vsc", main = paste("VSC heatmap of top ",length(topids)," genes by fold-change")
		)
		graphics.off()

		meanzscores = meanzscores[rownames(meanzscores) %in% topids , ]
		write.table(data.frame(id=rownames(meanzscores),meanzscores), row.names = F,
			file = file.path(odir,"heatmap.mean.vsc.zscores"), quote=F, sep="\t"
		);
		postscript(file.path(odir,"heatmap.mean.vsc.zscores.ps"))
		heatmap.2(as.matrix(meanzscores), col = color, Rowv = T, Colv = F, scale = "none",
			dendrogram = "row", trace = "none", margin = c(10, 8),
			cexRow = 0.6, cexCol = 0.8, key.title = NA, key.ylab = NA,
			key.xlab = "zscore", main = paste("Z-score heatmap of top ",length(topids)," genes by fold-change")
		)
		graphics.off()
	}

	cat(paste("number of significantly differentially expressed genes for ",ctr[i]," vs ",treat[i]," :",nrow(ddsr),"\n",sep=""))
}

for (i in 1:length(ctr)){
	odir = file.path(outdir,paste(ctr[i],"-vs-",treat[i],sep=""))
	dir.create(odir, recursive = T, showWarnings = F)

	thisdf = df[df$condition %in% c(ctr[i],treat[i]),]
	thisdesign = get_design(thisdf)

	if (design != thisdesign){
		cat(paste("warning: executing ",ctr[i]," vs ",treat[i]," independently from other samples with design formula: ",thisdesign,"\n",sep=""))
		ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = thisdf, directory = "", design = as.formula(design))
		thisdds = DESeq(ddsHTSeq, parallel = TRUE, BPPARAM = BPPARAM)
		save(thisdds, file = file.path(odir,"dds.Rdata"))
		get_table(thisdds)
	} else {
		cat(paste("executing ",ctr[i]," vs ",treat[i]," with design formula: ",design,"\n",sep=""))
		get_table(dds)
	}
}
