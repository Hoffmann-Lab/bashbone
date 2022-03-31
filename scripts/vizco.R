#! /usr/bin/env Rscript
# (c) Konstantin Riege, Jeanne Wilbrand
suppressMessages(library("pheatmap"))
suppressMessages(library("RColorBrewer"))

args <- commandArgs(TRUE);

type <- args[1] #e.g. Cluster or Module/Modules
datatype <- args[2] # 'TPM' (will be log transformed) or deseq 'VSC' else Z-Score
intsv <- args[3]
outbase <- args[4]

df = read.table(intsv, header=T, sep="\t", stringsAsFactors = F, quote="")
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
	rownames(tmp) <- dfg$id
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

write.table(data.frame(id=rownames(m),m,type=ac), quote=FALSE, row.names=F, sep="\t",
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


if(datatype != "TPM" && datatype != "VSC") quit()

suppressMessages(library(ggplot2))
suppressMessages(library(scales))
suppressMessages(library(reshape2))
suppressMessages(library(data.table))

# get cluster count
cluster_count <- length(unique(df$type))

# strip gene and cluster ID
foo <- df[,2:(ncol(df)-1)]

# log all data columns
if (datatype == "TPM") {
	flog <- log2(foo+1)
	labeldatatype=paste("log2",datatype,sep=" ")
} else {
	flog = foo
	labeldatatype=datatype
}
# get mean_pw_ch_log
if (nrow(flog)>1){
	res <- as.data.frame(rowMeans(mapply(function(x,y) abs(x-y), flog[, 2:ncol(flog)], flog[, 1:(ncol(flog)-1)])))
} else {
	res <- as.data.frame(mean(mapply(function(x,y) abs(x-y), flog[, 2:ncol(flog)], flog[, 1:(ncol(flog)-1)])))
}

colnames(res)[1] <- "mean_pw_ch_log"

# append gene name and cluster ID
flog2 <- cbind(flog,df[,c(1,ncol(df))], res)

molt <- reshape2::melt(flog2, id.vars = c("type", "id", "mean_pw_ch_log"), measure.vars = colnames(foo))

ggplot(molt, aes(x = variable, y = value, group = id, color = mean_pw_ch_log)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size = rel(min(1,15/ncol(df))))) +
	labs(x = "Type", y = labeldatatype, color = "Mean log2 FC", linetype = "") +
	geom_line(alpha=0.3) +
	stat_summary(aes(group = 1, linetype = ''), fun.y = 'median', geom = 'line', size = 1, show.legend = TRUE, colour = 'green') +
	scale_color_gradient(low = "blue", high = "red") +
	scale_linetype_discrete(name = paste("Median",labeldatatype, sep=" ")) +
	facet_wrap( ~ type , ncol = (as.integer(sqrt(cluster_count)+1)))
suppressMessages(ggsave(paste(outbase, "trajectories.pdf", sep=".")))
