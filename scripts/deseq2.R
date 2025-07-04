#! /usr/bin/env Rscript
# (c) Konstantin Riege

args = commandArgs(TRUE)

if(length(args)<5){
	cat("DESeq2 from raw feature counts plus extra effects, interaction terms, locally/globally clustered vst/z-score heatmaps by condition, pca from vst or rlog or estimated library size\n")
	cat("\n")
	cat("usage parameter: <i:threads> <f:experiments> <f:outdir> <s:condition1> <s:condition2> [<s:condition1> <s:condition3> ..]\n")
	cat('example: 16 "/path/to/experiments.csv" "/path/to/outdir" "ctr" "treat"\n')
	cat("\n")
	cat("experiments: comma separated with header. 4 or more columns. color by condition (column 2). dot shape by replicate (column 4). design formula and interaction terms by factor columns\n")
	cat("sample,countfile,condition,replicate[,factor1..]\n")
	cat("sample1,/path/to/countfile1,condition1,N1[,factorA..]\n")
	cat("sample2,/path/to/countfile2,condition1,N2[,factorB..]\n")
	cat("sample3,/path/to/countfile3,condition2,N1[,factorA..]\n")
	cat("..\n")
	quit("no",1)
}

cat("about to run deseq2 and pca\n")

options(warn=-1)

suppressMessages({
	library("DESeq2")
	library("BiocParallel")
	library("ggpubr")
	library("gplots")
	library("pheatmap")
	library("RColorBrewer")
	library("dplyr")
	library("gtools")
})

threads = as.numeric(args[1])
incsv = args[2] # sample,countfile,condition,replicate[,factor1,factor2,..]
outdir = args[3]
args = args[4:length(args)]
ctr = args[rep(c(T,F),length(args)/2)] # wt t1 wt t2 t1 t2
treat = args[rep(c(F,T),length(args)/2)]

dir.create(outdir, recursive = T, showWarnings = F)
BPPARAM = MulticoreParam(workers = threads)
setEPS(width=8, height=8, onefile=T)


##### deseq

experiments = read.table(incsv, header=T, sep=",", stringsAsFactors=F, check.names=F, quote="")
colnames(experiments)[1:4] = c("sample","countfile","condition","replicate")
experiments$condition=factor(experiments$condition,levels = unique(experiments$condition))
experiments$replicate=factor(experiments$replicate,levels = unique(mixedsort(experiments$replicate))) # use version sort via mixedsort() from gtools
if (!is.null(experiments$factor1)) {
	experiments$factor1 = factor(experiments$factor1,levels = unique(experiments$factor1))
}

# create design formula from factors under exclusion of potential linear combinations
# e.g. ~ factor1 + factor2 + condition
# note: deseq always needs more samples than coefficents
# e.g. with four samples, you can only fit three coefficients and have a residual degree of freedom for estimating variance
# y <- rnorm(4)
# > dat <- data.frame(a=factor(c(0,0,1,1)),b=factor(c(0,1,1,0)),c=factor(c(0,1,1,1)))
# > summary(lm(y ~ 1 + a + b, data=dat))$sigma
# [1] 0.0001074134
# > summary(lm(y ~ 1 + a + b + c, data=dat))$sigma
# [1] NaN
# The last line has no estimate of variance because the X matrix has four columns. The fitted values equal the observed:
# > all(y == lm(y ~ 1 + a + b + c, data=dat)$fitted)
# [1] TRUE
# 13:19
get_design = function(experiments,interactionterms=F){
	factors = c()
	if(length(colnames(experiments)) > 4){
		for (f in colnames(experiments)[5:length(colnames(experiments))]){
			v = experiments[,f]
			# levels(xyxy) -> 2 > 1 && xyxy @ aabb -> a:xy , b:xy -> {x,y} > 0 i.e. no linear combination
			# whereas levels(xxxx) -> 1 and xyzz @ aabb -> a:xy , b:z -> {} == 0 i.e. deseq error

			# turnoed off for now, should be handled by user
			#if (length(levels(as.factor(v))) > 1 && length(Reduce(intersect, split(v,experiments$condition))) > 0){
				factors = c(factors,f)
			#}
		}
	}
	if(interactionterms){
		return(paste("~",paste(c(factors,"condition"),collapse=" * ")))
	} else {
		return(paste("~",paste(c(factors,"condition"),collapse=" + ")))
	}
}

design = get_design(experiments)
cat(paste("using design formula: ",design,"\n",sep=""))

suppressMessages({
	dds = DESeqDataSetFromHTSeqCount(sampleTable = experiments, directory = "", design = as.formula(design))
	dds = tryCatch(
		{
			DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM, fitType="parametric")
		},
		error = function(e){
			# in case of too less genes/data points for parametric overdispersion fitting, use simple mean fitting
			DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM, fitType="mean")
		}
	)
	# for new deseq versions, that use local as fallback without throwing an error if parametric cannot be applied
	fit = attr(dispersionFunction(dds), "fitType")
	if(fit != "parametric" && fit != "mean"){
		dds = DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM, fitType="mean")
	}
})
save(dds, file = file.path(outdir,"dds.RData"))
# load(file = file.path(outdir,"dds.RData"))

pdf(file.path(outdir,"dispersion.pdf"))
plotDispEsts(dds)
graphics.off()

###### pca functions

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

######## pca

# log = DESeqTransform(SummarizedExperiment(log2(counts(dds, normalized=T) + 1), colData=colData(dds)))
# save(log, file = file.path(outdir,"log.RData"))
vsd = varianceStabilizingTransformation(dds, blind=FALSE)
save(vsd, file = file.path(outdir,"vsd.RData"))
rld = rlog(dds, blind=FALSE)
save(rld, file = file.path(outdir,"rld.RData"))
# load(file = file.path(outdir,"vsd.RData"))
# load(file = file.path(outdir,"rld.RData"))

# for (method in c("log","vsd","rld")){
for (method in c("vsd","rld")){
	# deseq method
	# data = plotPCA(get(normed), intgroup = c("condition", "replicate"), returnData = T)
	# percentVar = round(100 * attr(data, "percentVar"))

	normed = assay(get(method))
	vars = order(rowVars(normed), decreasing = TRUE)

	for (n in c(500,1000,2000,5000,10000,20000,50000,100000)){
		if(n>length(vars)) break
		topidx = vars[1:n]

		pca = prcomp(t(normed[topidx, ]), scale = F)
		loadings = pca$rotation
		for (pc in 1:3){
			ids = rownames(loadings[order(abs(loadings[,pc]), decreasing = TRUE),])
			ids = head(ids,n=max(1,length(ids)*0.05))
			sink(file.path(outdir,paste0("pca_PC",pc,"_",method,"_top",format(n, scientific=F),".variables")))
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

##### heatmap functions


gplots_heatmap = function(){
	if(min(df)<0) {
		# use red (+) and blue (-)
		color = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
		# define breaks from min to 0 and 0 to max according to number colors/2
		breaks = c(seq(min(df), 0, length.out=51), seq(max(df)/100, max(df), length.out=50))
		keylabel = "Z-score "
	} else {
		# use simple gradient
		color = colorRampPalette(brewer.pal(9, "GnBu"))(100)
		# define breaks equal to heatmap default
		breaks = seq(min(df), max(df), length.out=101)
		keylabel = "VSC "
	}

	# get number conditions colors and repeat them according to colannotation
	# colannotationcolor = colorRampPalette(brewer.pal(3, "Blues"))(length(unique(colannotation)))
	colannotationcolor = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(colannotation)))
	j=0
	colannotationcolor.full=c()
	for (g in unique(colannotation)){
		j=full_join+1
		colannotationcolor.full=c(colannotationcolor.full,rep(colannotationcolor[j],length(colannotation[colannotation %in% g])))
	}

	m = as.matrix(df)
	rownames(m) = rownames(df)

	postscript(paste0(path,".localclust.ps"))
	par(mar=c(0,0,0,0))
	# par(mar=c(0,0,4.2,0)) if topleft, to put legend below key
	heatmap.2(m, col = color, Rowv = T, Colv = coldendro, dendrogram = "both",
		breaks = breaks, colsep = colsep, sepcolor="black", sepwidth=0.01,
		trace = "none", margin = c(8, 8), lhei=c(1,6), lwid=c(1,3),
		symkey = F, key.title = NA, key.ylab = NA, key.par = list(cex=0.6), density.info="none", key.xlab = keylabel,
		cexRow = 0.8, cexCol = 0.9, srtCol=45,
		ColSideColors = colannotationcolor.full
	)
	legend("topright", inset = c(0.01,0), legend=unique(colannotation), col=colannotationcolor, lty = 1, lwd = 8, cex=0.72, bty="n", title="Group")
	# coords <- locator(1)      + mouse click , then use values e.g x=0 and y=1 (topleft) instead of keywords
	graphics.off()
}


get_heatmap = function(input,path){
	# perform inner-group column clustering

	cluster=TRUE
	colclustlist = list()
	for (condition in c(ctr[i],treat[i])) {
		m = as.matrix(input[ , experiments$sample[experiments$condition == condition]])
		rownames(m) = input$id

		# dist(df, method = "euclidean")
		# hclust(df, method = "complete")
		if(ncol(m)>1){
			colclust = hclust(dist(t(m)))
			m = m[ , colclust$order]
		} else {
			cluster=FALSE
			break
		}

		if(length(colclustlist)==0){
			df = data.frame(id=rownames(m),m,check.names=F)
			colclustlist = list(colclust)
			colsep = ncol(m)
			colannotation = rep(condition,ncol(m))
		} else {
			df = full_join(df, data.frame(id=rownames(m),m,check.names=F), by = 'id')
			colclustlist = c(colclustlist, list(colclust))
			colsep = c(colsep,tail(colsep,1)+ncol(m))
			colannotation = c(colannotation,rep(condition,ncol(m)))
		}
	}

	if(! cluster){
		get_heatmap_mean(input,paste0(path,".localclust"),FALSE)
		get_heatmap_mean(input,paste0(path,".globalclust"),TRUE)
		return()
	}

	# remove last separator (for gplots) and id-column
	colsep = colsep[1:length(colsep)-1]
	df = df[,2:length(df)]
	rownames(df) = input$id

	# merge dendrograms
	# multiforcation at +10% height, i.e. merge( as.dendrogram(colclustlist[[1]]) , as.dendrogram(colclustlist[[2]]), as.dendrogram(colclustlist[[3]]) )
	# can be handled by gplots, but not converted to hclust object, due to missing binary list structure
	# coldendro = do.call(merge, lapply(colclustlist, function(x) as.dendrogram(x)))
	# binary
	# coldendro = Reduce(merge, lapply(colclustlist, function(x) as.dendrogram(x)))
	# binary at same heights
	maxheight = max(unlist(lapply(colclustlist, function(x) max(x$height))))
	maxheight = maxheight + maxheight/10
	coldendro = Reduce(function(x, y) merge(x, y, height = maxheight), lapply(colclustlist, function(x) as.dendrogram(x)))
	# convert to hclust object for pheatmap
	colclust = as.hclust(coldendro)


	if(min(df)<0) {
		# use red (+) and blue (-)
		color = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
		# define breaks from min to 0 and 0 to max according to number colors/2
		breaks = c(seq(min(df), 0, length.out=51), seq(max(df)/100, max(df), length.out=50))
		# define breaks around 0 and labels
		legendbreaks = unique(c(seq(min(df), 0, length.out = round(-min(df)+1)) , seq(0, max(df), length.out = round(max(df)+1))))
		legendlabels = round(legendbreaks)
		legendlabels[legendlabels==0] = "Z-score "
	} else {
		# use simple gradient
		color = colorRampPalette(brewer.pal(9, "GnBu"))(100)
		# define breaks equal to heatmap default
		breaks = seq(min(df), max(df), length.out=101)
		# define breaks and labels
		legendbreaks = seq(min(df), max(df), length.out=7)
		legendlabels = round(legendbreaks)
		legendlabels[1] = "VSC "
	}

	# get number conditions colors and name them by conditions
	# colannotationcolor = colorRampPalette(brewer.pal(3, "Blues"))(length(unique(colannotation)))
	colannotationcolor = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(colannotation)))
	names(colannotationcolor) = unique(colannotation)
	# create dataframe of column annotations with rownames
	colannotations = data.frame(row.names = colnames(df), Group = colannotation)
	# create list of annotation colors according to dataframes of column and or row annotations
	annotationcolors = list(Group = colannotationcolor)

	p = pheatmap(df, color = color,	cluster_rows = T, cluster_cols = colclust, cutree_cols = length(colclustlist),
			annotation_col = colannotations, annotation_colors = annotationcolors, annotation_names_col = F,
			border_color = NA,  fontsize_row = 7, fontsize_col = 7, angle_col = 45,
			breaks = breaks, legend_breaks = legendbreaks, legend_labels = legendlabels,
			main = paste0("Heatmap of ",nrow(df)," most differentially expressed genes\n")
	)
	# filename = "/../../png|pdf|tiff|bmp"
	# workaround to store figure as different file type - use ps to later replace row and column names + ps2pdf
	postscript(paste0(path,".localclust.ps"))
	grid::grid.newpage()
	grid::grid.draw(p$gtable)
	graphics.off()


	# default like column clustered heatmap
	p = pheatmap(df, color = color, cluster_rows = T, cluster_cols = T,
		annotation_col = colannotations, annotation_colors = annotationcolors, annotation_names_col = F,
		border_color = NA,  fontsize_row = 7, fontsize_col = 7, angle_col = 45,
		breaks = breaks, legend_breaks = legendbreaks, legend_labels = legendlabels,
		main =paste0("Heatmap of ",nrow(df)," most differentially expressed genes\n")
	)
	postscript(paste0(path,".globalclust.ps"))
	grid::grid.newpage()
	grid::grid.draw(p$gtable)
	graphics.off()
}


get_heatmap_mean = function(input,path,cluster=F){
	# remove id-column
	df = input[,2:length(input)]
	rownames(df) = input$id

	if(min(df)<0) {
		# use red (+) and blue (-)
		color = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
		# define breaks from min to 0 and 0 to max according to number colors/2
		breaks = c(seq(min(df), 0, length.out=51), seq(max(df)/100, max(df), length.out=50))
		# define breaks around 0 and labels
		legendbreaks = unique(c(seq(min(df), 0, length.out = round(-min(df)+1)) , seq(0, max(df), length.out = round(max(df)+1))))
		legendlabels = round(legendbreaks)
		legendlabels[legendlabels==0] = "Z-score "
	} else {
		# use simple gradient
		color = colorRampPalette(brewer.pal(9, "GnBu"))(100)
		# define breaks equal to heatmap default
		breaks = seq(min(df), max(df), length.out=101)
		# define breaks and labels
		legendbreaks = seq(min(df), max(df), length.out=7)
		legendlabels = round(legendbreaks)
		legendlabels[1] = "VSC "
	}

	colannotation = experiments$condition[experiments$sample %in% colnames(df)]
	if(length(colannotation)==0){
		colannotation = colnames(df)
	}
	# get number conditions colors and name them by conditions
	# colannotationcolor = colorRampPalette(brewer.pal(3, "Blues"))(length(unique(colannotation)))
	colannotationcolor = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(colannotation)))
	names(colannotationcolor) = unique(colannotation)
	# create dataframe of column annotations with rownames
	colannotations = data.frame(row.names = colnames(df), Group = colannotation)
	# create list of annotation colors according to dataframes of column and or row annotations
	annotationcolors = list(Group = colannotationcolor)

	p = pheatmap(df, color = color, cluster_rows = T, cluster_cols = cluster,
			border_color = NA,  fontsize_row = 7, fontsize_col = 7, angle_col = 45,
			annotation_col = colannotations, annotation_colors = annotationcolors, annotation_names_col = F,
			breaks = breaks, legend_breaks = legendbreaks, legend_labels = legendlabels,
			main = paste0("Heatmap of ",nrow(df)," most differentially expressed genes\n")
	)
	# filename = "/../../png|pdf|tiff|bmp"
	# workaround to store figure as different file type - use ps to later replace row and column names + ps2pdf
	postscript(paste0(path,".ps"))
	grid::grid.newpage()
	grid::grid.draw(p$gtable)
	graphics.off()
}


###### fetch pairwise deseq results and plot heatmaps


get_table = function(dds){
	ddsrf = results(dds, contrast=c("condition",treat[i],ctr[i]), parallel = TRUE, BPPARAM = BPPARAM)
	ddsr = ddsrf

	ddsr = ddsr[order(abs(ddsr$log2FoldChange),decreasing=T) , ]
	write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F,
		file=file.path(odir,"deseq.full.tsv"), quote=F, sep="\t"
	)

	# ddsr = results(dds, contrast=c("condition",treat[i],ctr[i]), alpha = 0.05, parallel = TRUE, BPPARAM = BPPARAM) # alpha (default: 0.1) should be set to FDR cutoff
	# shrinkage can be used for data visualization and ranking of RNA-Seq data (remove high LFCs from lowly expressed genes with high variability among samples and estimate moderate FCs more close to reality) if DESeq() was executed with betaPrior=FALSE, which is the default since v1.16
	# 'apeglm' and 'ashr' outperform the original 'normal' shrinkage estimator. but ‘apeglm’ requires use of ‘coef’ i.e. a name or index from resultsNames(dds)
	# for ashr, if res is provided, then coef and contrast are ignored
	# ashr may be a little less aggressive than apeglm
	# use shrunken data for GSEA
	suppressMessages({
		ddsrshrunkf <- lfcShrink(dds, contrast=c("condition",treat[i],ctr[i]), res=ddsr, type="ashr", parallel = TRUE, BPPARAM = BPPARAM)
		ddsrshrunk = ddsrshrunkf
	})

	ddsrshrunk = ddsrshrunk[order(abs(ddsrshrunk$log2FoldChange),decreasing=T) , ]
	write.table(data.frame(id=rownames(ddsrshrunk),ddsrshrunk), row.names = F,
		file=file.path(odir,"deseq.full.fcshrunk.tsv"), quote=F, sep="\t"
	)

	ddsrshrunk = ddsrshrunk[!is.na(ddsrshrunk$log2FoldChange) , ]
	ddsrshrunk = ddsrshrunk[!is.na(ddsrshrunk$padj) , ]
	ddsrshrunk = ddsrshrunk[ddsrshrunk$baseMean > 0 , ]
	ddsrshrunk = ddsrshrunk[ddsrshrunk$padj <= 0.05 , ]
	write.table(data.frame(id=rownames(ddsrshrunk),ddsrshrunk), row.names = F,
		file=file.path(odir,"deseq.fcshrunk.tsv"), quote=F, sep="\t"
	)

	ddsr = ddsr[!is.na(ddsr$log2FoldChange) , ]
	ddsr = ddsr[!is.na(ddsr$padj) , ]
	ddsr = ddsr[ddsr$baseMean > 0 , ]
	ddsr = ddsr[ddsr$padj <= 0.05 , ]
	write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F,
		file=file.path(odir,"deseq.tsv"), quote=F, sep="\t"
	)

	if(nrow(ddsr)>0){
		pdf(file.path(odir,"ma_plot.pdf"))
		plotMA(ddsrf)
		graphics.off()
		# pdf(file.path(odir,"ma_plot.pdf"))
		# plotMA(ddsr)
		# graphics.off()
		pdf(file.path(odir,"ma_plot.fcshrunk.pdf"))
		plotMA(ddsrshrunkf)
		graphics.off()
	}

	#for (method in c("log","vsd","rld")){
	for (method in c("vsd","rld")){

		# old: use topmost variable features across all samples
		# normed = get(method)
		# normed = normed[,normed$condition %in% c(ctr[i],treat[i])]
		# topids = rownames(normed)[order(rowVars(assay(normed)), decreasing=T)]
		# inloop:
		# pca = prcomp(t(assay(normed[rownames(normed) %in% head(topids,n=n) , ])), scale = F)
		# percentVar = round(100*pca$sdev^2/sum(pca$sdev^2),1)

		# new: use topmost variable features of samples of current contrast
		normed = get(method)
		normed = assay(normed[,normed$condition %in% c(ctr[i],treat[i])])
		vars = order(rowVars(normed), decreasing=T)
		exp = experiments[experiments$condition %in% c(ctr[i],treat[i]), ]

		for (n in c(500,1000,2000,5000,10000,20000,50000,100000)){
			if(n>length(vars)) break
			topidx = vars[1:n]

			pca = prcomp(t(normed[topidx, ]), scale = F)
			loadings = pca$rotation
			for (pc in 1:3){
				ids = rownames(loadings[order(abs(loadings[,pc]), decreasing = TRUE),])
				ids = head(ids,n=max(1,length(ids)*0.05)) # variables that drive variation in PC
				sink(file.path(odir,paste0("pca_PC",pc,"_",method,"_top",format(n, scientific=F),".variables")))
				lapply(ids, cat, "\n")
				sink()
			}

			df = data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], condition = exp$condition, replicate = exp$replicate)
			write.table(data.frame(id=rownames(df),df), row.names = F,
				file=file.path(odir,paste0("pca_",method,"_top",format(n, scientific=F),".tsv")), quote=F, sep="\t"
			)

			if (!is.null(exp$factor1)) {
				df$facet = exp$factor1
			}
			percentVar = round(100*pca$sdev^2/sum(pca$sdev^2),1)
			df$PC1v=percentVar[1]
			df$PC2v=percentVar[2]
			df$PC3v=percentVar[3]

			suppressMessages({
				pca_plot(df, method, n, "PC1", "PC2", file.path(odir,"pca_plots_shaped"))
				pca_plot(df, method, n, "PC1", "PC3", file.path(odir,"pca_plots_shaped"))
				pca_plot(df, method, n, "PC2", "PC3", file.path(odir,"pca_plots_shaped"))
				pca_plot(df, method, n, "PC1", "PC2", file.path(odir,"pca_plots"), FALSE)
				pca_plot(df, method, n, "PC1", "PC3", file.path(odir,"pca_plots"), FALSE)
				pca_plot(df, method, n, "PC2", "PC3", file.path(odir,"pca_plots"), FALSE)
			})

			if (!is.null(exp$factor1)) {
				dfs = lapply(levels(exp$factor1), function(f) {
					e = exp[exp$factor1 == f,]
					nor = normed[,colnames(normed) %in% e$sample]
					v = order(rowVars(nor), decreasing = TRUE)
					topidx = v[1:n]
					pca = prcomp(t(nor[topidx, ]), scale = F)
					percentVar = round(100*pca$sdev^2/sum(pca$sdev^2),1)
					data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], PC3 = pca$x[,3], condition = e$condition, replicate = e$replicate,
						PC1v=percentVar[1], PC2v=percentVar[2], PC3v=percentVar[3], factor=f
					)
				})

				suppressMessages({
					pca_list_plot(dfs, method, n, "PC1", "PC2", file.path(odir,"pca_plots_shaped"))
					pca_list_plot(dfs, method, n, "PC1", "PC3", file.path(odir,"pca_plots_shaped"))
					pca_list_plot(dfs, method, n, "PC2", "PC3", file.path(odir,"pca_plots_shaped"))
					pca_list_plot(dfs, method, n, "PC1", "PC2", file.path(odir,"pca_plots"), FALSE)
					pca_list_plot(dfs, method, n, "PC1", "PC3", file.path(odir,"pca_plots"), FALSE)
					pca_list_plot(dfs, method, n, "PC2", "PC3", file.path(odir,"pca_plots"), FALSE)
				})
			}
		}
	}

	vsdr = vsd[,vsd$condition %in% c(ctr[i],treat[i])]

	vsc = as.data.frame(assay(vsdr))
	write.table(data.frame(id=rownames(vsc),vsc,check.names=F), row.names = F,
		file=file.path(odir,"experiments.vsc"), quote=F, sep="\t"
	)

	zscores = vsc #log(vsc+1)
	zscores = zscores-rowMeans(zscores)
	zscores = zscores/apply(zscores,1,sd)
	zscores[is.na(zscores)] = 0
	write.table(data.frame(id=rownames(zscores),zscores,check.names=F), row.names = F,
		file=file.path(odir,"experiments.vsc.zscores"), quote=F, sep="\t"
	)

	colnamesvsc = colnames(vsc)
	colnames(vsc) = as.vector(vsdr$condition)
	meanvsc = t(apply(vsc, 1, function(x) tapply(x, colnames(vsc), mean)))
	meanvsc = meanvsc[,c(ctr[i],treat[i])]
	colnames(vsc) = colnamesvsc
	write.table(data.frame(id=rownames(meanvsc),meanvsc,check.names=F), row.names = F,
		file = file.path(odir,"experiments.mean.vsc"), quote=F, sep="\t"
	)
	meanzscores = meanvsc #log(meanvsc+1)
	meanzscores = meanzscores-rowMeans(meanzscores)
	meanzscores = meanzscores/apply(meanzscores,1,sd)
	meanzscores[is.na(meanzscores)] = 0
	write.table(data.frame(id=rownames(meanzscores),meanzscores,check.names=F), row.names = F,
		file=file.path(odir,"experiments.mean.vsc.zscores"), quote=F, sep="\t"
	)

	if(nrow(ddsr)>1){
		color = colorRampPalette(brewer.pal(9, "GnBu"))(100)

		topids = head(rownames(ddsrshrunk),n=50)

		vsc = vsc[rownames(vsc) %in% topids , ]
		write.table(data.frame(id=rownames(vsc),vsc,check.names=F), row.names = F,
			file=file.path(odir,"heatmap.vsc"), quote=F, sep="\t"
		)
		get_heatmap(data.frame(id=rownames(vsc),vsc,check.names=F),file.path(odir,"heatmap.vsc"))

		zscores = zscores[rownames(zscores) %in% topids , ]
		write.table(data.frame(id=rownames(zscores),zscores,check.names=F), row.names = F,
			file=file.path(odir,"heatmap.vsc.zscores"), quote=F, sep="\t"
		)
		get_heatmap(data.frame(id=rownames(zscores),zscores,check.names=F),file.path(odir,"heatmap.vsc.zscores"))

		meanvsc = meanvsc[rownames(meanvsc) %in% topids , ]
		write.table(data.frame(id=rownames(meanvsc),meanvsc,check.names=F), row.names = F,
			file = file.path(odir,"heatmap.mean.vsc"), quote=F, sep="\t"
		)
		get_heatmap_mean(data.frame(id=rownames(meanvsc),meanvsc,check.names=F),file.path(odir,"heatmap.mean.vsc"))

		meanzscores = meanzscores[rownames(meanzscores) %in% topids , ]
		write.table(data.frame(id=rownames(meanzscores),meanzscores,check.names=F), row.names = F,
			file = file.path(odir,"heatmap.mean.vsc.zscores"), quote=F, sep="\t"
		)
		get_heatmap_mean(data.frame(id=rownames(meanzscores),meanzscores,check.names=F),file.path(odir,"heatmap.mean.vsc.zscores"))
	}

	cat(paste("number of significantly differentially expressed genes for ",ctr[i]," vs ",treat[i],": ",nrow(ddsr),"\n",sep=""))
}

######

makename = function(name){
	name = make.names(paste0("ADAPTER",name)) # replaces \W characters by '.' and prepend X if starting with \W
	return(sub("ADAPTER","",name))
}

get_interactionterm = function(dds,experiments,outdir,ctr,treat){

	design_interactionterms = get_design(experiments,T)
	# each condition needs all factors - twice to be able to run Ca_Fx vs Cb_Fx with replicates
	# use try catch to ensure this
	suppressMessages({
		try = tryCatch(
			{
				design(dds) = as.formula(design_interactionterms)
				DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM)
			},
			error = function(e){
				cat("WARNING: ",file=stderr())
				message(e)
				return(NA)
			}
		)
	})
	if(is.na(try)){
		message(". calculating effects (interaction terms) of secondary factors not possible.")
		return()
	} else {
		cat("calculating effects (interaction terms) of secondary factors\n")
	}

	for (f in 5:ncol(experiments)){
		fac = colnames(experiments)[f]
		terms = unique(experiments[,f])

		if(! is.factor(dds[[fac]])){
			next # in case removed from design formula because not present in this experiments
		}

		for (i in 1:(length(terms)-1)){
			ctr_term = as.character(terms[i]) # need to remove factor for releveling
			# avoid alphabetically sorted levels. ensure first factor/condition listed is the one to be used as main/reference
			dds$condition = relevel(dds$condition, ref = ctr)
			dds[[fac]] = relevel(dds[[fac]], ref = ctr_term)

			cat(paste("calculating effects (interaction terms) of secondary factors with reference ",ctr_term," on ",ctr," vs ",treat," with design formula: ",design_interactionterms,"\n",sep=""))
            suppressMessages({
                dds = DESeq(dds, parallel = TRUE, BPPARAM = BPPARAM)
			})
			# as suggested upon relevel: run nbinomWaldTest or DESeq() i.e. estimate size factors <- automatically re-used, estimate dispersion <- no need to re-do it, nbinomWaldTest/nbinomLRT
			# LRT tests multiple terms at once e.g. multiple levels of a factor or all interactions between two variables
			# dds <- nbinomLRT(dds, reduced = ???)

			for (j in (i+1):length(terms)){
				term = terms[j]

				odir = file.path(outdir,"interactionterms",paste0(ctr_term,"-vs-",term))
				dir.create(odir, recursive = T, showWarnings = F)
				save(dds, file = file.path(odir,"dds.RData"))

				#cat(paste("calculating effect (interaction term) of ",ctr_term," vs ",term," on ",ctr," vs ",treat," with design formula: ",design_interactionterms,"\n",sep=""))
				ddsr = results(dds, name=paste0(fac,makename(term),".condition",makename(treat)), parallel = TRUE, BPPARAM = BPPARAM)

                ddsr = ddsr[order(abs(ddsr$log2FoldChange),decreasing=T) , ]
				write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F,
					file=file.path(odir,"deseq.full.tsv"), quote=F, sep="\t"
				)

				suppressMessages({
					ddsrshrunkf <- lfcShrink(dds, coef=paste0(fac,makename(term),".condition",makename(treat)), res=ddsr, type="ashr", parallel = TRUE, BPPARAM = BPPARAM)
					ddsrshrunk = ddsrshrunkf
				})

				ddsrshrunk = ddsrshrunk[order(abs(ddsrshrunk$log2FoldChange),decreasing=T) , ]
				write.table(data.frame(id=rownames(ddsrshrunk),ddsrshrunk), row.names = F,
					file=file.path(odir,"deseq.full.fcshrunk.tsv"), quote=F, sep="\t"
				)

				ddsrshrunk = ddsrshrunk[!is.na(ddsrshrunk$log2FoldChange) , ]
				ddsrshrunk = ddsrshrunk[!is.na(ddsrshrunk$padj) , ]
				ddsrshrunk = ddsrshrunk[ddsrshrunk$baseMean > 0 , ]
				ddsrshrunk = ddsrshrunk[ddsrshrunk$padj <= 0.05 , ]
				write.table(data.frame(id=rownames(ddsrshrunk),ddsrshrunk), row.names = F,
					file=file.path(odir,"deseq.fcshrunk.tsv"), quote=F, sep="\t"
				)

				ddsr = ddsr[!is.na(ddsr$log2FoldChange) , ]
				ddsr = ddsr[!is.na(ddsr$padj) , ]
				ddsr = ddsr[ddsr$baseMean > 0 , ]
				ddsr = ddsr[ddsr$padj <= 0.05 , ]
				write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F,
					file=file.path(odir,"deseq.tsv"), quote=F, sep="\t"
				)

				cat(paste("number of significantly affected genes by interaction term ",ctr_term," vs ",term," on ",ctr," vs ",treat,": ",nrow(ddsr),"\n",sep=""))
			}
		}
	}
}


##### iterate over pairs of conditions, recompute diff expr. according to experiment factors, get deseq results + vsc + z-scores + heatmaps

toplus=list()
for (i in 1:length(ctr)){
	if(length(toplus[[ctr[i]]])==0){
		toplus[[ctr[i]]]=treat[i]
	} else {
		toplus[[ctr[i]]]=unique(c(toplus[[ctr[i]]],treat[i]))
	}

	odir = file.path(outdir,paste0(ctr[i],"-vs-",treat[i]))
	dir.create(odir, recursive = T, showWarnings = F)

	thisexperiments = experiments[experiments$condition %in% c(ctr[i],treat[i]),]
	# thisdesign = get_design(thisexperiments)
	thisdds = dds

	cat(paste("calculating differential expression of ",ctr[i]," vs ",treat[i]," with design formula: ",design,"\n",sep=""))
	# cat(paste("calculating differential expression of ",ctr[i]," vs ",treat[i]," with design formula: ",thisdesign,"\n",sep=""))
	# if (design != thisdesign){ # use in case of matrix not full rank error
	# 	design(thisdds) = as.formula(thisdesign)
	# 	thisdds = DESeq(thisdds, parallel = TRUE, BPPARAM = BPPARAM)
	# 	save(thisdds, file = file.path(odir,"dds.RData"))
	# }
	get_table(thisdds)
	if(ncol(experiments) > 4){
		get_interactionterm(thisdds,thisexperiments,odir,ctr[i],treat[i])
	}
}

for (ctr in names(toplus)){
	treats=toplus[[ctr]]
	dds$condition = relevel(dds$condition,ref=ctr)
	suppressMessages({
		dds = nbinomWaldTest(dds)
	})

	if(length(treats)>1){
		for (i in 1:(length(treats))){
			for (j in 1:length(treats)){
				if(i==j){
					next
				}
				odir = file.path(outdir,paste0(ctr,"-vs-",treats[i]),"extra_effects",paste0(ctr,"-vs-",treats[j]))
				dir.create(odir, recursive = T, showWarnings = F)
				save(dds, file = file.path(odir,"dds.RData"))

				cat(paste("calculating extra effects of ",ctr," vs ",treats[j]," on top of ",ctr," vs ",treats[i],"\n",sep=""))
				ddsr = results(dds, contrast=list(paste0("condition_",makename(treats[i]),"_vs_",makename(ctr)),paste0("condition_",makename(treats[j]),"_vs_",makename(ctr))), parallel = TRUE, BPPARAM = BPPARAM)
				ddsr$log2FoldChange = -1*ddsr$log2FoldChange

				ddsr = ddsr[order(abs(ddsr$log2FoldChange),decreasing=T) , ]
				write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F,
					file=file.path(odir,"deseq.full.tsv"), quote=F, sep="\t"
				)

				suppressMessages({
					ddsrshrunkf <- lfcShrink(dds, contrast=list(paste0("condition_",makename(treats[i]),"_vs_",makename(ctr)),paste0("condition_",makename(treats[j]),"_vs_",makename(ctr))), res=ddsr, type="ashr", parallel = TRUE, BPPARAM = BPPARAM)
					ddsrshrunk = ddsrshrunkf
				})

				ddsrshrunk = ddsrshrunk[order(abs(ddsrshrunk$log2FoldChange),decreasing=T) , ]
				write.table(data.frame(id=rownames(ddsrshrunk),ddsrshrunk), row.names = F,
					file=file.path(odir,"deseq.full.fcshrunk.tsv"), quote=F, sep="\t"
				)

				ddsrshrunk = ddsrshrunk[!is.na(ddsrshrunk$log2FoldChange) , ]
				ddsrshrunk = ddsrshrunk[!is.na(ddsrshrunk$padj) , ]
				ddsrshrunk = ddsrshrunk[ddsrshrunk$baseMean > 0 , ]
				ddsrshrunk = ddsrshrunk[ddsrshrunk$padj <= 0.05 , ]
				write.table(data.frame(id=rownames(ddsrshrunk),ddsrshrunk), row.names = F,
					file=file.path(odir,"deseq.fcshrunk.tsv"), quote=F, sep="\t"
				)

				ddsr = ddsr[!is.na(ddsr$log2FoldChange) , ]
				ddsr = ddsr[!is.na(ddsr$padj) , ]
				ddsr = ddsr[ddsr$baseMean > 0 , ]
				ddsr = ddsr[rev(order(abs(ddsr$log2FoldChange))) , ]
				ddsr = ddsr[ddsr$padj <= 0.05 , ]
				write.table(data.frame(id=rownames(ddsr),ddsr), row.names = F,
					file=file.path(odir,"deseq.tsv"), quote=F, sep="\t"
				)

				#cat(paste("number of significantly extra effects of ",ctr," vs ",treats[j]," on top of ",ctr," vs ",treats[i],": ",nrow(ddsr),"\n",sep=""))
			}
		}
	}
}
