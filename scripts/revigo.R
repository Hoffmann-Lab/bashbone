#! /usr/bin/env Rscript
# (c) Konstantin Riege, Arne Sahm
suppressMessages(library("httr"))

args = commandArgs(TRUE)
go = read.table(args[1], header=FALSE, sep="\t")
outprefix = args[2]
domain = args[3]

revigo = function(GOtable,destFolder,cutoff="0.70",isPValue="yes",whatIsBetter="higher",goSizes="0",measure="SIMREL",title=paste(sep="", "Gene Ontology treemap - ", gsub("_"," ",domain))){
	goList = paste0(apply(GOtable,1,function(x) paste0(x[1]," ",x[2])),collapse = "\n")
	x = POST(url="http://revigo.irb.hr/revigo.jsp", encode="form", body=list(goList=goList, cutoff=cutoff, isPValue=isPValue, whatIsBetter=whatIsBetter, goSizes=goSizes, measure=measure, startRevigo="Please wait"))
	if (x$status_code == 200){

		# in case of webgestalt gmt database contains all three domains fetch table 1 to 3
		# if (domain == "biological_process"){ 
			treemap_R = content(GET("http://revigo.irb.hr/toR_treemap.jsp?table=1"),"text")	
			write(content(GET("http://revigo.irb.hr/export_treemap.jsp?table=1"),"text"),paste(sep="",outprefix,".treemap.csv"))
			scatter_R = content(GET("http://revigo.irb.hr/toR.jsp?table=1"),"text")
			write(content(GET("http://revigo.irb.hr/export.jsp?table=1"),"text"),paste(sep="",outprefix,".scatterplot.csv"))
		# } else if (domain == "cellular_component") {
		# 	treemap_R = content(GET("http://revigo.irb.hr/toR_treemap.jsp?table=2"),"text")	
		# 	write(content(GET("http://revigo.irb.hr/export_treemap.jsp?table=2"),"text"),paste(sep="",outprefix,".treemap.csv"))
		# 	scatter_R = content(GET("http://revigo.irb.hr/toR.jsp?table=2"),"text")
		# 	write(content(GET("http://revigo.irb.hr/export.jsp?table=2"),"text"),paste(sep="",outprefix,".scatterplot.csv"))
		# } else {
		# 	treemap_R = content(GET("http://revigo.irb.hr/toR_treemap.jsp?table=3"),"text")	
		# 	write(content(GET("http://revigo.irb.hr/export_treemap.jsp?table=3"),"text"),paste(sep="",outprefix,".treemap.csv"))
		# 	scatter_R = content(GET("http://revigo.irb.hr/toR.jsp?table=3"),"text")
		# 	write(content(GET("http://revigo.irb.hr/export.jsp?table=3"),"text"),paste(sep="",outprefix,".scatterplot.csv"))
		# }
		
		treemap_R = gsub("\r","",treemap_R)
		treemap_R = gsub("REVIGO Gene Ontology treemap",title,treemap_R)
		treemap_R = gsub("file=[^,]+",paste(sep="","file='",outprefix,".treemap.pdf'"),treemap_R)
		save(treemap_R, file=paste(sep="",outprefix,".treemap.Rdata"))
		eval(parse(text=treemap_R))
		
		scatter_R = gsub("\r","",scatter_R)
		save(scatter_R, file=paste(sep="",outprefix,".scatterplot.Rdata"))
		eval(parse(text=paste0(scatter_R,paste(sep="","\nggsave(\"",outprefix,".scatterplot.pdf\")"))))	
	} else {
		stop()
	}
}

try = 0
while(try < 6){
	try = try + 1
	tryCatch({
		print(paste(sep="","start revigo try ",try," on ",outprefix))
		revigo(go,out)
		quit(status=0)
	}, error = function(e) {
		print(paste(sep="","catch revigo try ",try," on ",outprefix))
		Sys.sleep(1)
	})
}
print(paste(sep="","revigo failed on ",outprefix))
quit(status=1)
