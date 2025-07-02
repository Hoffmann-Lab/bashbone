#! /usr/bin/env Rscript
# (c) Konstantin Riege
args = commandArgs(TRUE)

if(length(args)<2){
	cat("pca 1vs2, 1vs3, 2vs3 from raw feature counts, normalized by DESeq2 vst, rlog or estimated library size\n")
	cat("\n")
	cat("usage parameter: <f:experiments> <f:outdir> [<b:isnormalized>]\n")
	cat('example: "/path/to/experiments.csv" "/path/to/outdir"\n')
	cat("\n")
	cat("experiments: comma separated WITH header. 4 or more columns. color by condition (column 2). dot shape by replicate (column 4).\n")
	cat("sample,path,condition,replicate[,..]\n")
	cat("sample1,/path/to/countfile1,condition1,N1[,..]\n")
	cat("sample2,/path/to/countfile2,condition1,N2[,..]\n")
	cat("sample3,/path/to/countfile3,condition2,N1[,..]\n")
	cat("..\n")
	quit("no",1)
}

cat("about to run pca\n")

options(warn=-1)
options(scipen = 999)

suppressMessages({
	library("DESeq2")
	library("ggpubr")
	library("RColorBrewer")
	# library("ggrepel")
	library("dplyr")
	library("gtools")
})

incsv = args[1]
outdir = args[2]
if(length(args)==3) isnormalized=args[3] else isnormalized=FALSE

experiments = read.table(incsv, header=T, sep=",", stringsAsFactors=F, check.names=F, quote="")
colnames(experiments)[1:4] = c("sample","countfile","condition","replicate")
experiments$condition=factor(experiments$condition,levels = unique(experiments$condition))
experiments$replicate=factor(experiments$replicate,levels = unique(mixedsort(experiments$replicate))) # use version sort via mixedsort() from gtools
if (!is.null(experiments$factor1)) {
	experiments$factor1 = factor(experiments$factor1,levels = unique(experiments$factor1))
}

dds = DESeqDataSetFromHTSeqCount(sampleTable = experiments, directory = "", design = ~1)
dds = estimateSizeFactors(dds)

dir.create(outdir, recursive = T, showWarnings = F)

if(isnormalized) {
	raw = DESeqTransform(SummarizedExperiment(counts(dds, normalized=F), colData=colData(dds))) # counts is memberfunction of DESeqDataSet
	save(raw, file = file.path(outdir,"raw.RData"))
	methods=c("raw")
} else {
	log = DESeqTransform(SummarizedExperiment(log2(counts(dds, normalized=T) + 1), colData=colData(dds))) # normalized=T devides by library size factors
	save(log, file = file.path(outdir,"log.RData"))
	vsd = varianceStabilizingTransformation(dds, blind=TRUE)
	save(vsd, file = file.path(outdir,"vsd.RData"))
	rld = rlog(dds, blind=FALSE)
	save(rld, file = file.path(outdir,"rld.RData"))
	#methods=c("log","vsd","rld")
	methods=c("vsd","rld")
}

pca_plot = function(df,method,n,PCx="PC1",PCy="PC2",outdir,shaped=TRUE,pal="Set1"){
	dir.create(outdir, recursive = T, showWarnings = F)
	maxcol = length(palette(pal))
	ncols = length(unique(df$condition))

	PCxVar=df[1,colnames(df)==paste0(PCx,"v")]
	PCyVar=df[1,colnames(df)==paste0(PCy,"v")]

	p = ggplot(df, aes(x=!!sym(PCx), y=!!sym(PCy), color = condition, shape = if(shaped) !!sym("replicate") else NULL)) +
		ggtitle(paste0("PCA of top ",n," variable features")) +
		theme_minimal() +
		theme(aspect.ratio=1, legend.box="horizontal", legend.title=element_blank()) +
		geom_point(size = 1.5) +
		scale_color_manual(values = colorRampPalette(brewer.pal(min(ncols,maxcol), pal))(ncols)) +
		scale_fill_manual(values = colorRampPalette(brewer.pal(min(ncols,maxcol), pal))(ncols)) +
		scale_shape_manual(values = c(1:length(unique(df$replicate)))) +
		xlab(paste0(PCx,": ",PCxVar, "% variance")) +
		ylab(paste0(PCy,": ",PCyVar, "% variance")) +
		geom_hline(yintercept=0, col="black", linetype="dashed", size=0.3) +
		geom_vline(xintercept=0, col="black", linetype="dashed", size=0.3)
	xl = max(abs(df[[PCx]]))
	yl = max(abs(df[[PCy]]))
	p + xlim(-xl,xl) +
		ylim(-yl,yl)
	ggsave(file.path(outdir,paste0("pca_",PCx,PCy,"_",method,"_top",format(n, scientific=F),".pdf")), width = 8, height = 6)

	if(!is.null(df$facet)){
		p + facet_wrap(~facet , scales = "free") +
			xlim(-xl,xl) +
			ylim(-yl,yl)
		ggsave(file.path(outdir,paste0("pca_",PCx,PCy,"_",method,"_top",format(n, scientific=F),".faceted.pdf")), width = 4*nlevels(df$facet), height = 4)
	}

	# + stat_ellipse(geom = "polygon", alpha=0.1, level = 0.85, aes(fill=condition))
	# => slightly different from ellipse::ellipse. does not work for less then 4 data-points, sometimes even less than 3
	# + ggforce::geom_mark_ellipse(aes(fill=condition)) +
	# => works for > 1 points, but completely differently computed i.e. rather geometric enclosing shape based than covariance based
	# for full control also regarding axis limits, better calculate manually
	ellipses = df %>% group_by(condition) %>% group_map(~ {
		df = data.frame(x=.x[[PCx]],y=.x[[PCy]])
		as.data.frame(ellipse::ellipse(cov(df), centre = colMeans(df), level = 0.85)) %>% mutate(condition = .y$condition)
	}) %>% bind_rows()
	xl = max(abs(ellipses$x))
	yl = max(abs(ellipses$y))
	p + geom_polygon(data = ellipses, aes(x, y, color=condition, fill=condition), alpha = 0.1, inherit.aes = F) +
		xlim(-xl,xl) +
		ylim(-yl,yl)
	ggsave(file.path(outdir,paste0("pca_",PCx,PCy,"_",method,"_top",format(n, scientific=F),".ellipses.pdf")), width = 8, height = 6)

	if(!is.null(df$facet)){
		ellipses = df %>% group_by(condition,facet) %>% group_map(~ {
			df = data.frame(x=.x[[PCx]],y=.x[[PCy]])
			as.data.frame(ellipse::ellipse(cov(df), centre = colMeans(df), level = 0.85)) %>% mutate(condition = .y$condition, facet = .y$facet)
		}) %>% bind_rows()
		p + geom_polygon(data = ellipses, aes(x, y, color=condition, fill=condition), alpha = 0.1, inherit.aes = F) +
			facet_wrap(~facet , scales = "free") +
			xlim(-xl,xl) +
			ylim(-yl,yl)
		ggsave(file.path(outdir,paste0("pca_",PCx,PCy,"_",method,"_top",format(n, scientific=F),".faceted.ellipses.pdf")), width = 4*nlevels(df$facet), height = 4)
	}
}

pca_list_plot = function(dfs,method,n,PCx="PC1",PCy="PC2",outdir,shaped=TRUE,pal="Set1"){
	dir.create(outdir, recursive = T, showWarnings = F)
	maxcol = length(palette(pal))

	plots = lapply(dfs, function(df) {
		ncols = length(unique(df$condition))
		PCxVar=df[1,colnames(df)==paste0(PCx,"v")]
		PCyVar=df[1,colnames(df)==paste0(PCy,"v")]

		p = ggplot(df, aes(x=!!sym(PCx), y=!!sym(PCy), color = condition, shape = if(shaped) !!sym("replicate") else NULL)) +
			ggtitle(df$factor[1]) +
			theme_minimal() +
			theme(aspect.ratio=1, legend.box = "horizontal", legend.title=element_blank(), plot.title = element_text(hjust = 0.5, face = "plain", size=11)) +
			geom_point(size = 1.5) +
			scale_color_manual(values = colorRampPalette(brewer.pal(min(ncols,maxcol), pal))(ncols)) +
			scale_fill_manual(values = colorRampPalette(brewer.pal(min(ncols,maxcol), pal))(ncols)) +
			xlab(paste0(PCx,": ",PCxVar, "% variance")) +
			ylab(paste0(PCy,": ",PCyVar, "% variance")) +
			geom_hline(yintercept=0, col="black", linetype="dashed", size=0.3) +
			geom_vline(xintercept=0, col="black", linetype="dashed", size=0.3)
		xl = max(abs(df[[PCx]]))
		yl = max(abs(df[[PCy]]))
		p + xlim(-xl,xl) +
			ylim(-yl,yl)
	})
	annotate_figure(
		ggarrange(plotlist = plots, nrow = 1, common.legend = T, legend = "right", align = "hv"),
		top = text_grob(paste0("PCA of top ",n," variable genes"), hjust = 1.46)
	)
	ggsave(file.path(outdir,paste0("pca_",PCx,PCy,"_",method,"_top",format(n, scientific=F),".splitted.pdf")), width = 1+4*length(dfs), height = 4)

	plots = lapply(dfs, function(df) {
		ncols = length(unique(df$condition))
		PCxVar=df[1,colnames(df)==paste0(PCx,"v")]
		PCyVar=df[1,colnames(df)==paste0(PCy,"v")]

		p = ggplot(df, aes(x=!!sym(PCx), y=!!sym(PCy), color = condition, shape = if(shaped) !!sym("replicate") else NULL)) +
			ggtitle(df$factor[1]) +
			theme_minimal() +
			theme(aspect.ratio=1, legend.box = "horizontal", legend.title=element_blank(), plot.title = element_text(hjust = 0.5, face = "plain", size=11)) +
			geom_point(size = 1.5) +
			scale_color_manual(values = colorRampPalette(brewer.pal(min(ncols,maxcol), pal))(ncols)) +
			scale_fill_manual(values = colorRampPalette(brewer.pal(min(ncols,maxcol), pal))(ncols)) +
			xlab(paste0(PCx,": ",PCxVar, "% variance")) +
			ylab(paste0(PCy,": ",PCyVar, "% variance")) +
			geom_hline(yintercept=0, col="black", linetype="dashed", size=0.3) +
			geom_vline(xintercept=0, col="black", linetype="dashed", size=0.3)
		ellipses = df %>% group_by(condition) %>% group_map(~ {
			df = data.frame(x=.x[[PCx]],y=.x[[PCy]])
			as.data.frame(ellipse::ellipse(cov(df), centre = colMeans(df), level = 0.85)) %>% mutate(condition = .y$condition)
		}) %>% bind_rows()
		xl = max(abs(ellipses$x))
		yl = max(abs(ellipses$y))
		p + geom_polygon(data = ellipses, aes(x, y, color=condition, fill=condition), alpha = 0.1, inherit.aes = F) +
			xlim(-xl,xl) +
			ylim(-yl,yl)
	})
	annotate_figure(
		ggarrange(plotlist = plots, nrow = 1, common.legend = T, legend = "right", align = "hv"),
		top = text_grob(paste0("PCA of top ",n," variable genes"), hjust = 1.46)
	)
	ggsave(file.path(outdir,paste0("pca_",PCx,PCy,"_",method,"_top",format(n, scientific=F),".splitted.ellipses.pdf")), width = 1+4*length(dfs), height = 4)
}

for (method in methods){
	data = as.data.frame(assay(get(method)))
	write.table(data.frame(id=rownames(data),data,check.names=F), row.names = F,
		file=file.path(outdir,paste0("experiments.",method)), quote=F, sep="\t"
	)

	normed = assay(get(method))
	vars = order(rowVars(normed), decreasing = TRUE)

	for (n in c(500,1000,2000,5000,10000,20000,50000,100000)){
		if(n>length(vars)) break
		topidx = vars[1:n]

		pca = prcomp(t(normed[topidx, ]), scale = F)
		loadings = pca$rotation
		for (i in 1:3){
			ids = rownames(loadings[order(abs(loadings[,i]), decreasing = TRUE),])
			ids = head(ids,n=max(1,length(ids)*0.05)) # 5% top variables that drive variation in PC
			sink(file.path(outdir,paste0("pca_PC",i,"_",method,"_top",format(n, scientific=F),".variables")))
			lapply(ids, cat, "\n")
			sink()
		}

		df = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], condition = experiments$condition, replicate = experiments$replicate)
		write.table(data.frame(id=rownames(df),df), row.names = F,
			file=file.path(outdir,paste0("pca_",method,"_top",format(n, scientific=F),".tsv")), quote=F, sep="\t"
		)

		if (!is.null(experiments$factor1)) {
			df$facet = experiments$factor1
		}
		percentVar = round(100*pca$sdev^2/sum(pca$sdev^2),1)
		df$PC1v=percentVar[1]
		df$PC2v=percentVar[2]
		df$PC3v=percentVar[3]

		suppressMessages({
			pca_plot(df, method, n, "PC1", "PC2", file.path(outdir,"pca_plots_shaped"))
			pca_plot(df, method, n, "PC1", "PC3", file.path(outdir,"pca_plots_shaped"))
			pca_plot(df, method, n, "PC2", "PC3", file.path(outdir,"pca_plots_shaped"))
			pca_plot(df, method, n, "PC1", "PC2", file.path(outdir,"pca_plots"), FALSE)
			pca_plot(df, method, n, "PC1", "PC3", file.path(outdir,"pca_plots"), FALSE)
			pca_plot(df, method, n, "PC2", "PC3", file.path(outdir,"pca_plots"), FALSE)
		})

		if (!is.null(experiments$factor1)) {
			dfs = lapply(levels(experiments$factor1), function(f) {
				exp = experiments[experiments$factor1 == f,]
				nor = normed[,colnames(normed) %in% exp$sample]
				v = order(rowVars(nor), decreasing = TRUE)
				topidx = v[1:n]
				pca = prcomp(t(nor[topidx, ]), scale = F)
				percentVar = round(100*pca$sdev^2/sum(pca$sdev^2),1)
				data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], condition = exp$condition, replicate = exp$replicate,
					PC1v=percentVar[1], PC2v=percentVar[2], PC3v=percentVar[3], factor=f
				)
			})

			suppressMessages({
				pca_list_plot(dfs, method, n, "PC1", "PC2", file.path(outdir,"pca_plots_shaped"))
				pca_list_plot(dfs, method, n, "PC1", "PC3", file.path(outdir,"pca_plots_shaped"))
				pca_list_plot(dfs, method, n, "PC2", "PC3", file.path(outdir,"pca_plots_shaped"))
				pca_list_plot(dfs, method, n, "PC1", "PC2", file.path(outdir,"pca_plots"), FALSE)
				pca_list_plot(dfs, method, n, "PC1", "PC3", file.path(outdir,"pca_plots"), FALSE)
				pca_list_plot(dfs, method, n, "PC2", "PC3", file.path(outdir,"pca_plots"), FALSE)
			})
		}
	}
}
