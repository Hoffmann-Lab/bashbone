#! /usr/bin/env Rscript
# (c) Konstantin Riege
options(warn=-1)

cat("about to run pca\n")

suppressMessages({
	library("DESeq2")
	library("BiocParallel")
	library("ggplot2")
	library("gplots")
	library("pheatmap")
	library("RColorBrewer")
	library("dplyr")
})

args = commandArgs(TRUE)
threads = as.numeric(args[1])
incsv = args[2] # sample,countfile,condition,replicate[,factor1,factor2,..]
outdir = args[3]

dir.create(outdir, recursive = T, showWarnings = F)
BPPARAM = MulticoreParam(workers = threads)
setEPS(width=8, height=8, onefile=T)

experiments = read.table(incsv, header=T, sep=",", stringsAsFactors=F)
colnames(experiments)[1:4] = c("sample","countfile","condition","replicate")

# create design formula from factors under exclusion of potential linear combinations
# e.g. ~ factor1 + factor2 + condition
get_design = function(experiments){
	factors = c()
	if(length(colnames(experiments)) > 4){
		for (f in colnames(experiments)[5:length(colnames(experiments))]){
			v = experiments[,f]
			# levels(xyxy) -> 2 > 1 && xyxy @ aabb -> a:xy , b:xy -> {x,y} > 0 i.e. no linear combination
			# whereas levels(xxxx) -> 1 and xyzz @ aabb -> a:xy , b:z -> {} == 0 i.e. deseq error
			if (length(levels(as.factor(v))) > 1 && length(Reduce(intersect, split(v,experiments$condition))) > 0){
				factors = c(factors,f)
			}
		}
	}
	return(paste("~",paste(c(factors,"condition"),collapse=" + ")))
}

design = get_design(experiments)
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable = experiments, directory = "", design = as.formula(design))
dds = DESeq(ddsHTSeq, parallel = TRUE, BPPARAM = BPPARAM)
dds = dds[ rowSums(counts(dds)) > 1, ]

log = DESeqTransform(SummarizedExperiment(log2(counts(dds, normalized=T) + 1), colData=colData(dds)))
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
		ggsave(file.path(outdir,paste("pca_12_",method,".pdf",sep="")))
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
		ggsave(file.path(outdir,paste("pca_13_",method,".pdf",sep="")))
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
		ggsave(file.path(outdir,paste("pca_23_",method,".pdf",sep="")))
		# stat_ellipse() +
	})
}