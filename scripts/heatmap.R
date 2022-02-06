#! /usr/bin/env Rscript
# (c) Konstantin Riege

suppressMessages(library("pheatmap"))
suppressMessages(library("gplots"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("dplyr"))

args = commandArgs(TRUE)
cluster = args[1]
# ggplot default is 7x7
w = as.integer(args[2])
h = as.integer(args[3])
experiments = read.table(args[4], header = T, sep = ',', stringsAsFactors = F)
#experiments = data.frame(condition = df$condition, sample = df$sample)
io = args[5]
input = read.table(io, header = T, sep = '\t', stringsAsFactors = F)
keylabel = args[6]
# appended to "Heatmap of <nrow> "
title = args[7]

if(nrow(input)<2){
  quit("no")
}

# setup postscript device
setEPS(reset=T, width=w, height=h, onefile=T)
# use w=8.25, h=11.7 and convert via ps2pdf -sPAPERSIZE=a4 in.ps out.pdf
# or ps2pdf -g5950x8420 since 1 inch = 72 PostScript-points (pt) and a4 = 595x842pt = 5950x8420px
# or ps2pdf $(grep -m 1 -P '(BoundingBox|DocumentMedia)' in.ps | awk '{print "-g"$4*10"x"$5*10}') in.ps out.pdf
# cat(paste0("ps2pdf -g",round(w*72)*10,"x",round(h*72)*10," ",io,".localclust.ps ",io,".localclust.pdf\n"))

if(cluster){
  ##### do inner-condition clustering and store annotation and separation infos
  colclustlist = list()
  for (condition in unique(experiments$condition)) {
    m = as.matrix(input[ , colnames(input) %in% gsub("\\W",".",experiments$sample[experiments$condition %in% condition]) ])
    rownames(m) = input$id

    # dist(df, method = "euclidean")
    # hclust(df, method = "complete")
    colclust = hclust(dist(t(m)))
    m = m[ , colclust$order]

    if(length(colclustlist)==0){
      df = data.frame(id=rownames(m),m)
      colclustlist = list(colclust)
      colsep = ncol(m)
      colannotation = rep(condition,ncol(m))
    } else {
      df = full_join(df, data.frame(id=rownames(m),m), by = 'id')
      colclustlist = c(colclustlist, list(colclust))
      colsep = c(colsep,tail(colsep,1)+ncol(m))
      colannotation = c(colannotation,rep(condition,ncol(m)))
    }
  }

  # remove last separator and id-column
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
    legendlabels[legendlabels==0] = paste0(keylabel," ")
  } else {
    # use simple gradient
    color = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    # define breaks equal to heatmap default
    breaks = seq(min(df), max(df), length.out=101)
    # define breaks and labels
    legendbreaks = seq(min(df), max(df), length.out=7)
    legendlabels = round(legendbreaks)
    legendlabels[1] = paste0(keylabel," ")
  }

  # get number conditions colors and name them by conditions
  # colannotationcolor = colorRampPalette(brewer.pal(3, "Blues"))(length(unique(colannotation)))
  colannotationcolor = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(colannotation)))
  names(colannotationcolor) = unique(colannotation)
  # create dataframe of column annotations with rownames
  colannotations = data.frame(row.names = colnames(df), Group = colannotation)
  # create list of annotation colors according to dataframes of column and or row annotations
  annotationcolors = list(Group = colannotationcolor)

  p = pheatmap(df, color = color, cluster_rows = T, cluster_cols = colclust, cutree_cols = length(colclustlist),
               annotation_col = colannotations, annotation_colors = annotationcolors, annotation_names_col = F,
               border_color = NA,  fontsize_row = 7, fontsize_col = 7, angle_col = 45,
               breaks = breaks, legend_breaks = legendbreaks, legend_labels = legendlabels,
               main = paste0("Heatmap of ",nrow(df)," ",title,"\n")
  )
  # filename = "/../../png|pdf|tiff|bmp"
  # workaround to store figure as different file type - use ps to later replace row and column names + ps2pdf
  postscript(paste0(io,".localclust.ps"))
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  graphics.off()
  #cat(paste0("ps2pdf -g",round(w*72)*10,"x",round(h*72)*10," ",io,".localclust.ps ",io,".localclust.pdf\n"))


  p = pheatmap(df, color = color, cluster_rows = T, cluster_cols = T,
               annotation_col = colannotations, annotation_colors = annotationcolors, annotation_names_col = F,
               border_color = NA,  fontsize_row = 7, fontsize_col = 7, angle_col = 45,
               breaks = breaks, legend_breaks = legendbreaks, legend_labels = legendlabels,
               main = paste0("Heatmap of ",nrow(df)," ",title,"\n")
  )
  # filename = "/../../png|pdf|tiff|bmp"
  # workaround to store figure as different file type - use ps to later replace row and column names + ps2pdf
  postscript(paste0(io,".globalclust.ps"))
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  graphics.off()
 # cat(paste0("ps2pdf -g",round(w*72)*10,"x",round(h*72)*10," ",io,".globalclust.ps ",io,".globalclust.pdf\n"))

} else {

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
    legendlabels[legendlabels==0] = paste0(keylabel," ")
  } else {
    # use simple gradient
    color = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    # define breaks equal to heatmap default
    breaks = seq(min(df), max(df), length.out=101)
    # define breaks and labels
    legendbreaks = seq(min(df), max(df), length.out=7)
    legendlabels = round(legendbreaks)
    legendlabels[1] = paste0(keylabel," ")
  }

  p = pheatmap(df, color = color, cluster_rows = T, cluster_cols = F,
               border_color = NA,  fontsize_row = 7, fontsize_col = 7, angle_col = 45,
               breaks = breaks, legend_breaks = legendbreaks, legend_labels = legendlabels,
               main = paste0("Heatmap of ",nrow(df)," ",title,"\n")
  )
  # filename = "/../../png|pdf|tiff|bmp"
  # workaround to store figure as different file type - use ps to later replace row and column names + ps2pdf
  postscript(paste0(io,".ps"))
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  graphics.off()
  #cat(paste0("ps2pdf -g",round(w*72)*10,"x",round(h*72)*10," ",io,".ps ",io,".pdf\n"))
}

quit("no")

###### working alternative via gplots

if(min(df)<0) {
  # use red (+) and blue (-)
  color = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
  # define breaks from min to 0 and 0 to max according to number colors/2
  breaks = c(seq(min(df), 0, length.out=51), seq(max(df)/100, max(df), length.out=50))
} else {
  # use simple gradient
  color = colorRampPalette(brewer.pal(9, "GnBu"))(100)
  # define breaks equal to heatmap default
  breaks = seq(min(df), max(df), length.out=101)
}

# get number conditions colors and repeat them according to colannotation
# colannotationcolor = colorRampPalette(brewer.pal(3, "Blues"))(length(unique(colannotation)))
colannotationcolor = colorRampPalette(brewer.pal(11, "Spectral"))(length(unique(colannotation)))
i=0
colannotationcolor.full=c()
for (g in unique(colannotation)){
  i=i+1
  colannotationcolor.full=c(colannotationcolor.full,rep(colannotationcolor[i],length(colannotation[colannotation %in% g])))
}

m = as.matrix(df)
rownames(m) = rownames(df)

postscript("/u/people/kriege/heatmap.ps")
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
# coords <- locator(1)      + mouse click , then use values for x= and y= instead of keywords
graphics.off()
