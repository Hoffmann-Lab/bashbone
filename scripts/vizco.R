#! /usr/bin/env Rscript
# (c) Konstantin Riege
args = commandArgs(TRUE)

if(length(args)<4){
	cat("plot feature heatmap and TPM/VSC trajectories per wgcna cluster\n")
	cat("\n")
	cat("usage parameter: <[Cluster|Module]:type> <[VSC|TPM|Z-score|s]:data-type> <f:matrix> <f:outbase> [<f:matrix>]\n")
	cat('example: Cluster Z-score "/path/to/matrix.tsv" "/path/to/outbase" "/path/to/corr_p_padj_matrix.tsv"\n')
	cat("\n")
	cat("matrix: tab separated with header and feature ids/label.\n")
	cat("id       sample1 sample2 sample3 ..\n")
	cat("feature1 value1  value2  value3  ..\n")
	cat("..\n")
	cat("\n")
	cat("(optional) corr_p_padj_matrix.tsv: tab separated with header and feature ids/label.\n")
	cat("id       cluster      Cluster001.cor [..] Cluster002.cor ..\n")
	cat("feature1 Cluster.001  value1         [..] value2  ..\n")
	cat("feature2 Cluster.002  value1         [..] value2  ..\n")
	cat("..\n")
	quit("no",1)
}

options(warn=-1)

suppressMessages({
	library(pheatmap)
	library(RColorBrewer)
	library(ggplot2)
	library(scales)
	library(reshape2)
	library(data.table)
})

type <- args[1] #e.g. Cluster or Module/Modules
datatype <- args[2] # 'TPM' (will be log transformed) or deseq 'VSC' else Z-Score
intsv <- args[3]
outbase <- args[4]

MEcor = data.frame()
if(length(args)==5){
	MEcor = read.table(args[5], header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="")
	colnames(MEcor)[2] = "type"
}

df = read.table(intsv, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="")
colnames(df)[1]="id"
colnames(df)[ncol(df)]="type"

ar <- c()
ac <- c()
gr <- c()
m <- matrix()
for (t in unique(df$type)){
	dfg <- df[df$type %in% t , ]
	l <- nrow(dfg)
	if(l<2) next
	tmp <- as.matrix(dfg[ , 2:(ncol(dfg)-1)])
	rownames(tmp) <- paste0(dfg$type,".",dfg$id)
	c <- hclust(dist(tmp))
	tmp <- tmp[c$order , ]
	ar <- c(ar, rownames(tmp))
	ac <- c(ac, rep(t, l))
	gr <- c(gr, gr[length(gr)]+l)
	if(length(m[1,])>1){
		m <- rbind(m,tmp)
	} else {
		m <- tmp
	}
}
ar <- data.frame( row.names = ar, dummy = ac)
if (nrow(ar) == 0) quit() #necessary if script is called with emty input
colnames(ar)[1] <- type
#if(length(unique(ac))==1) ar <- NA

labeldatatype <- datatype
if (datatype == "TPM") {
	m <- log2(m+1)
	labeldatatype <- paste("log2",datatype,sep=" ")
}

write.table(data.frame(id=rownames(m),m,type=ac,check.names=F), quote=FALSE, row.names=F, sep="\t",
	file = paste(outbase, "heatmap.tsv", sep=".")
)

# mmin <- min(m)
# mmax <- max(m)
# breaks <- c(mmin,mmin+0.5,(mmax+mmin)/2,mmax-0.5,mmax)
# label <- c('',as.integer(mmin),labeldatatype,as.integer(mmax),'')
# pheatmap(m, color=rev(colorRampPalette(brewer.pal(9, 'RdBu'))(100)),
# 	cluster_cols = F, cluster_rows = F, filename = paste(outbase, "heatmap.pdf", sep="."),
# 	gaps_row = gr, annotation_row = ar, main = "Heatmap",
# 	annotation_names_row = F, show_rownames = F, border_color=NA,
# 	legend_breaks = breaks, legend_labels = label
# )

if(min(m)<0) {
	# use red (+) and blue (-)
	color = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
	# define breaks from min to 0 and 0 to max according to number colors/2
	breaks = c(seq(min(m), 0, length.out=51), seq(max(m)/100, max(m), length.out=50))
	# define breaks around 0 and labels
	legendbreaks = unique(c(seq(min(m), 0, length.out = round(-min(m)+1)) , seq(0, max(m), length.out = round(max(m)+1))))
	legendlabels = round(legendbreaks)
	legendlabels[legendlabels==0] = paste0(labeldatatype," ")
} else {
	# use simple gradient
	color = colorRampPalette(brewer.pal(9, "GnBu"))(100)
	# define breaks equal to heatmap default
	breaks = seq(min(m), max(m), length.out=101)
	# define breaks and labels
	legendbreaks = seq(min(m), max(m), length.out=7)
	legendlabels = round(legendbreaks)
	legendlabels[1] = paste0(labeldatatype," ")
}

pheatmap(m, color=color,
	cluster_cols = F, cluster_rows = F, filename = paste(outbase, "heatmap.pdf", sep="."),
	gaps_row = gr, annotation_row = ar, main = "Heatmap",
	annotation_names_row = F, show_rownames = F, border_color=NA,
	breaks = breaks, legend_breaks = legendbreaks, legend_labels = legendlabels
)


################# trajectories


# get cluster count
cluster_count <- length(unique(df$type))

# strip gene and cluster ID
foo <- df[,2:(ncol(df)-1)]
foo.colnames <- colnames(foo)

if (datatype == "TPM") {
	foo <- log2(foo+1)
	labeldatatype <- paste("log2",datatype,sep=" ")
} else {
	labeldatatype <- datatype
}

if (nrow(MEcor)>0 && type == "Cluster"){
	foo <- cbind(foo, df[,c(1,ncol(df))])
	res = data.frame()
	for (i in unique(MEcor$type)){
		cor = subset(MEcor[MEcor$type==i,],select=c("id",paste0(i,".cor")))
		cor[2] = abs(cor[2])
		colnames(cor)[2] = "me_corr"
		res = rbind(res, cor)
	}
	foo = merge(foo,res,by="id")
	molt <- reshape2::melt(foo, id.vars = c("type", "id", "me_corr"), measure.vars = foo.colnames)
	ggplot(molt, aes(x = variable, y = value, group = id, color = me_corr)) +
		theme_minimal() +
		theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(min(1,10/ncol(df))))) +
		labs(x = "", y = labeldatatype, color = "Abs. Correlation", linetype = "") +
		geom_line(alpha=0.3) +
		stat_summary(aes(group = 1, linetype = ''), fun.y = 'median', geom = 'line', size = 0.8, show.legend = TRUE, colour = 'green') +
		scale_color_gradient(low = "blue", high = "red") +
		scale_linetype_discrete(name = paste("Median",labeldatatype, sep=" ")) +
		facet_wrap( ~ type , ncol = (as.integer(sqrt(cluster_count)+1)))
	suppressMessages(ggsave(paste(outbase, "trajectories.pdf", sep=".")))

} else {
	if(ncol(m) < 3){
		foo <- cbind(foo, df[,c(1,ncol(df))])
		molt <- reshape2::melt(df, id.vars = c("type", "id"), measure.vars = foo.colnames)
		ggplot(molt, aes(x = variable, y = value, group = id)) +
			theme_minimal() +
			theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(min(1,10/ncol(df))))) +
			labs(x = "", y = labeldatatype, color = "Mean FC", linetype = "") +
			geom_line(alpha=0.3) +
			stat_summary(aes(group = 1, linetype = ''), fun.y = 'median', geom = 'line', size = 0.8, show.legend = TRUE, colour = 'green') +
			scale_linetype_discrete(name = paste("Median",labeldatatype, sep=" ")) +
			facet_wrap( ~ type , ncol = (as.integer(sqrt(cluster_count)+1)))
		suppressMessages(ggsave(paste(outbase, "trajectories.pdf", sep = ".")))
	} else {
		# get mean pairwise (fold) changes -> mean_pw_ch
		if (nrow(foo)>1){
			res <- as.data.frame(rowMeans(mapply(function(x,y) abs(x-y), foo[, 2:ncol(foo)], foo[, 1:(ncol(foo)-1)])))
		} else {
			res <- as.data.frame(mean(mapply(function(x,y) abs(x-y), foo[, 2:ncol(foo)], foo[, 1:(ncol(foo)-1)])))
		}
		colnames(res)[1] <- "mean_pw_ch"
		foo <- cbind(foo, df[,c(1,ncol(df))], res)
		molt <- reshape2::melt(foo, id.vars = c("type", "id", "mean_pw_ch"), measure.vars = foo.colnames)
		ggplot(molt, aes(x = variable, y = value, group = id, color = mean_pw_ch)) +
			theme_minimal() +
			theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = rel(min(1,10/ncol(df))))) +
			labs(x = "", y = labeldatatype, color = "Mean FC", linetype = "") +
			geom_line(alpha=0.3) +
			stat_summary(aes(group = 1, linetype = ''), fun.y = 'median', geom = 'line', size = 0.8, show.legend = TRUE, colour = 'green') +
			scale_color_gradient(low = "blue", high = "red") +
			scale_linetype_discrete(name = paste("Median",labeldatatype, sep=" ")) +
			facet_wrap( ~ type , ncol = (as.integer(sqrt(cluster_count)+1)))
		suppressMessages(ggsave(paste(outbase, "trajectories.pdf", sep=".")))
	}
}
