#! /usr/bin/env Rscript
# (c) Konstantin Riege
options(warn=-1)

cat("about to run pca\n")

suppressMessages({
	library("DESeq2")
	library("ggplot2")
	library("RColorBrewer")
})

args = commandArgs(TRUE)
incsv = args[1] # sample,countfile,condition,replicate[,factor1,factor2,..] - condition can be NA or empty
outdir = args[2]

dir.create(outdir, recursive = T, showWarnings = F)
setEPS(width=8, height=8, onefile=T)

experiments = read.table(incsv, header=T, sep=",", stringsAsFactors=F, quote="")
colnames(experiments)[1:4] = c("sample","countfile","condition","replicate")

# method 1
# raw = matrix(..) # of read counts - equal to assay(DESeqDataSetFromHTSeqCount(..))
# varianceStabilizingTransformation(raw) # works from matrix
# in order to get simple pca of counts devided by size factor, estimateSizeFactors requires deseq object
# dds <- DESeqDataSetFromMatrix(countData = raw, colData = data.frame(row.names = colnames(raw)), design = ~1)
# ...
# method 2
dds = DESeqDataSetFromHTSeqCount(sampleTable = experiments, directory = "", design = ~1)
dds = estimateSizeFactors(dds)

log = DESeqTransform(SummarizedExperiment(log2(counts(dds, normalized=T) + 1), colData=colData(dds))) # normalized=T devides by library size factors
save(log, file = file.path(outdir,"log.Rdata"))
vsd = varianceStabilizingTransformation(dds, blind=FALSE)
save(vsd, file = file.path(outdir,"vsd.Rdata"))
rld = rlog(dds, blind=FALSE)
save(rld, file = file.path(outdir,"rld.Rdata"))

for (method in c("log","vsd","rld")){
	# deseq method
	# data = plotPCA(get(normed), intgroup = c("condition", "replicate"), returnData = T)
	# percentVar = round(100 * attr(data, "percentVar"))

	normed = assay(get(method))
	vars = rowVars(normed)
	n = 100
	#n = length(vars)
	topidx = order(vars, decreasing = TRUE)[1:min(n,length(vars))]
	pca = prcomp(t(normed[topidx, ]), scale = F)
	percentVar = round(100*pca$sdev^2/sum(pca$sdev^2),1)
	data = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], PC4 = pca$x[,4], replicate = experiments$replicate, condition = experiments$condition)
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
			xlab(paste0("PC1: ",percentVar[1], "% variance")) +
			ylab(paste0("PC2: ",percentVar[2], "% variance"))
		suppressMessages(ggsave(file.path(outdir,paste("pca_12_",method,".pdf",sep=""))))
		# stat_ellipse() +

		ggplot(data, aes(PC1, PC3, color = condition, group = condition, shape = replicate)) +
			ggtitle("PCA plot - PC1 vs PC3") +
			scale_shape_manual(values = c(1:length(unique(data$replicate)) )) +
			coord_fixed() +
			theme_bw() +
			theme(legend.box = "horizontal") +
			geom_point(size = 3) +
			xlab(paste0("PC1: ",percentVar[1], "% variance")) +
			ylab(paste0("PC3: ",percentVar[3], "% variance"))
		suppressMessages(ggsave(file.path(outdir,paste("pca_13_",method,".pdf",sep=""))))
		# stat_ellipse() +

		ggplot(data, aes(PC2, PC3, color = condition, group = condition, shape = replicate)) +
			ggtitle("PCA plot - PC2 vs PC3") +
			scale_shape_manual(values = c(1:length(unique(data$replicate)) )) +
			coord_fixed() +
			theme_bw() +
			theme(legend.box = "horizontal") +
			geom_point(size = 3) +
			xlab(paste0("PC2: ",percentVar[2], "% variance")) +
			ylab(paste0("PC3: ",percentVar[3], "% variance"))
		suppressMessages(ggsave(file.path(outdir,paste("pca_23_",method,".pdf",sep=""))))
		# stat_ellipse() +
	})
}