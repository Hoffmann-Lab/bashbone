#! /usr/bin/env bash
# (c) Konstantin Riege

function enrichment::_ora(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-1 <cmds1>    | array of
			-d <title>    | prefix for plots. in case of GO may be domain
			-g <gmtfile>  | path to
			-r <refterms> | path to
			-i <idsfile>  | path to
			-j <tpmfile>  | tpm experiments file path to or -u
			-u <uiverse>  | path to or -j
			-o <outdir>   | path to

			legacy mode:
			-1 <cmds1>    | array of
			-d <domain>   | grepped for in gofile to generate gmt file and used as prefix for plot titles
			-g <gofile>   | path to
			-i <idsfile>  | path to
			-j <tpmfile>  | tpm experiments file path to
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory threads domain gmtfile refterms idsfile tpmtsv outdir universe
	declare -n _cmds1_ora
	while getopts '1:d:g:r:i:j:u:o:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_ora=$OPTARG;;
			d)	((++mandatory)); domain="$OPTARG";;
			g)	((++mandatory)); gmtfile="$OPTARG";;
			r)	refterms="$OPTARG";;
			i)	((++mandatory)); idsfile="$OPTARG";;
			j)	tpmtsv="$OPTARG";;
			u)	universe="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	if [[ -s "$gmtfile" && $(awk 'NF>1' "$gmtfile" | head -1 | cut -f 2) =~ GO:[0-9]+ ]]; then
		# for backwards compatibility
		declare -a cmdgmt
		if [[ $tpmtsv ]]; then
			commander::makecmd -a cmdgmt -s '|' -o "$outdir/reference.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD'
				grep -Fw "$domain" "$gmtfile" | tee -i >(cut -f 2,4 | sort -u > "$outdir/reference.terms")
			CMD
				grep -Fw -f <(perl -M'List::Util qw(max)' -lanE 'say \$F[0] if max(@F)>=1 && \$.>1' "$tpmtsv")
			CMD
				perl -F'\t' -lane '
					push @{$m{$F[1]}},$F[0];
					END{
						print join"\t",($_,"NA",@{$m{$_}}) for keys %m
					}
				'
			CMD
		else
			commander::makecmd -a cmdgmt -s '|' -o "$outdir/reference.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
				grep -Fw "$domain" "$gmtfile" | tee -i >(cut -f 2,4 | sort -u > "$outdir/reference.terms")
			CMD
				perl -F'\t' -lane '
					push @{$m{$F[1]}},$F[0];
					END{
						print join"\t",($_,"NA",@{$m{$_}}) for keys %m
					}
				'
			CMD
		fi
		gmtfile="$outdir/reference.gmt"
		refterms="$outdir/reference.terms"
		commander::runcmd -v -b -i 1 -a cmdgmt
	else
		[[ $refterms ]] || _usage
		[[ $universe ]] && tpmtsv="$universe"
	fi

	[[ $universe ]] || universe="$tpmtsv"
	[[ $universe ]] || universe="/dev/null"

	commander::makecmd -a _cmds1_ora -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			suppressMessages(library("clusterProfiler"));
			suppressMessages(library("DOSE"));
			suppressMessages(library("enrichplot"));
			suppressMessages(library("ggplot2"));
			options(warn=-1);
			args <- commandArgs(TRUE);
			gmt <- args[1];
			g2n <- args[2];
			idsfile <- args[3];
			unifile <- args[4];
			odir <- args[5];
			domain <- args[6];

			genes <- scan(idsfile, character(), quote="", quiet=T);
			universe <- scan(unifile, character(), quote="", quiet=T);
			ora <- data.frame(matrix(ncol = 6, nrow = 0));
			if(length(genes)>0){
				tg <- suppressMessages(read.gmt(gmt));
				tn <- read.table(g2n, sep="\t", stringsAsFactors=F, check.names=F, quote="");
				if(sum(genes %in% tg$gene)>0){
					if(length(universe)>1){
						ora <- enricher(genes, universe=universe, TERM2GENE=tg, TERM2NAME=tn, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500);
					} else {
						ora <- enricher(genes, TERM2GENE=tg, TERM2NAME=tn, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500);
					};
					if(!is.null(ora) && nrow(ora)>0){
						n <- min(20,nrow(ora));
						fs <- round(10-max(10,n)/10+5);
						dotplot(ora, title=paste0(domain," (Top ",n,")"), showCategory=n, label_format=max(nchar(ora$Description))+1) +
							theme(axis.text.y = element_text(size=fs), axis.text.x = element_text(size=fs), axis.title.x = element_text(size=fs), strip.background = element_rect(linetype = 0, fill=NA)) +
							guides(color = guide_colourbar("FDR"));
						suppressMessages(ggsave(file.path(odir,"dotplot.pdf"), width=10, height=round(log10(n+1)*8)));

						n <- nrow(ora);
						fs <- round(10-max(10,n)/10+5);
						dotplot(ora, title=paste0(domain," (",n,")"), showCategory=n, label_format=max(nchar(ora$Description))+1) +
							theme(axis.text.y = element_text(size=fs), axis.text.x = element_text(size=fs), axis.title.x = element_text(size=fs), strip.background = element_rect(linetype = 0, fill=NA)) +
							guides(color = guide_colourbar("FDR"));
						suppressMessages(ggsave(file.path(odir,"dotplot.full.pdf"), width=10, height=round(log10(n+1)*8)));

						df <- data.frame(fdr=-log(ora$p.adjust, base=10), description = paste0(ora$Description," (",ora$Count,")"));
						n <- min(20,nrow(df));
						dfp <- head(df,n);
						pdf(file.path(odir,"barplot.pdf"),width=max(nchar(df$description))/5,height=n/2+1.8);
						par(mar=c(5,max(nchar(df$description))/3+3,5,1));
						barplot(rev(dfp$fdr), main=paste0(domain," (Top ",n,")"), horiz=T, names.arg=rev(dfp$description),
							xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$fdr)+max(dfp$fdr)/10*2),
							cex.names=0.8, las=1);
						graphics.off();

						n <- nrow(df);
						dfp <- df;
						pdf(file.path(odir,"barplot.full.pdf"),width=max(nchar(df$description))/5,height=n/2+1.8);
						par(mar=c(5,max(nchar(df$description))/3+3,5,1));
						barplot(rev(dfp$fdr), main=paste0(domain," (",n,")"), horiz=T, names.arg=rev(dfp$description),
							xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$fdr)+max(dfp$fdr)/10*2),
							cex.names=0.8, las=1);
						graphics.off();

						ora <- as.data.frame(ora)[c("ID","Count","pvalue","p.adjust","Description","geneID")];
						ora$geneID <- gsub("/",";",ora$geneID);
					} else {
						ora <- data.frame(matrix(ncol = 6, nrow = 0));
					};
				};
			};
			colnames(ora) = c("id","count","pval","padj","description","enrichment");
			write.table(ora, row.names = F, file = file.path(odir,"goenrichment.tsv"), quote=F, sep="\t");
		'
	CMD
		"$gmtfile" "$refterms" "$idsfile" <(perl -M'List::Util qw(max)' -lanE 'say \$#F==0? \$F[0] : \$F[0] if max(@F)>=1 && \$.>1' "$universe") "$outdir" "$domain"
	CMD

	return 0
}

function enrichment::_gsea(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-1 <cmds1>    | array of
			-d <title>    | prefix for plots. in case of GO may be domain
			-g <gmtfile>  | path to
			-r <refterms> | path to
			-i <deseqtsv> | path to
			-o <outdir>   | path to

			legacy mode:
			-1 <cmds1>    | array of
			-d <domain>   | grepped for in gofile to generate gmt file and used as prefix for plot titles
			-g <gofile>   | path to
			-i <deseqtsv> | path to
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory threads domain gmtfile refterms deseqtsv outdir
	declare -n _cmds1_gsea
	while getopts '1:d:g:r:i:j:o:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_gsea=$OPTARG;;
			d)	((++mandatory)); domain="$OPTARG";;
			g)	((++mandatory)); gmtfile="$OPTARG";;
			r)	refterms="$OPTARG";;
			i)	((++mandatory)); deseqtsv="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	if [[ -s "$gmtfile" && $(awk 'NF>1' "$gmtfile" | head -1 | cut -f 2) =~ GO:[0-9]+ ]]; then
		# for backwards compatibility
		declare -a cmdgmt
		commander::makecmd -a cmdgmt -s '|' -o "$outdir/reference.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			grep -Fw "$domain" "$gmtfile" | tee -i >(cut -f 2,4 | sort -u > "$outdir/reference.terms")
		CMD
			perl -F'\t' -lane '
				push @{$m{$F[1]}},$F[0];
				END{
					print join"\t",($_,"NA",@{$m{$_}}) for keys %m
				}
			'
		CMD
		gmtfile="$outdir/reference.gmt"
		refterms="$outdir/reference.terms"
		commander::runcmd -v -b -i 1 -a cmdgmt
	else
		[[ $refterms ]] || _usage
	fi

	commander::makecmd -a _cmds1_gsea -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			suppressMessages(library("clusterProfiler"));
			suppressMessages(library("DOSE"));
			suppressMessages(library("enrichplot"));
			suppressMessages(library("ggplot2"));
			args <- commandArgs(TRUE);
			gmt <- args[1];
			g2n <- args[2];
			ddsr <- args[3];
			odir <- args[4];
			domain <- args[5];

			df <- read.table(ddsr, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
			df <- df[!is.na(df$padj) , ];
			df <- df[df$padj<=0.05 , ];

			gsea <- data.frame(matrix(ncol = 6, nrow = 0));
			if(nrow(df)>0){
				tg <- suppressMessages(read.gmt(gmt));
				tn <- read.table(g2n, sep="\t", stringsAsFactors=F, check.names=F, quote="");
				if(sum(df$id %in% tg$gene)>0){
					gl <- df$log2FoldChange;
					names(gl) <- df$id;
					gl <- sort(gl, decreasing = T);
					gsea <- GSEA(gl, TERM2GENE=tg, TERM2NAME=tn, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500, eps = 0, seed = T);

					if(!is.null(gsea) && nrow(gsea)>0){
						n <- min(20,nrow(gsea));
						fs <- round(10-max(10,n)/10+5);
						dotplot(gsea, title=paste0(domain," (Top ",n,")"), showCategory=n, split=".sign", label_format=max(nchar(gsea$Description))+1) +
							facet_grid(.~.sign) +
							theme(axis.text.y = element_text(size=fs), axis.text.x = element_text(size=fs), axis.title.x = element_text(size=fs), strip.background = element_rect(linetype = 0, fill=NA)) +
							guides(color = guide_colourbar("FDR"));
						suppressMessages(ggsave(file.path(odir,"dotplot.pdf"), width=10, height=round(log10(n+1)*8)));

						ridgeplot(gsea, showCategory=n, label_format=max(nchar(gsea$Description))+1) + xlab("log2 FC distribution") +
							guides(fill = guide_colourbar("FDR")) + ggtitle(paste0(domain," (Top ",n,")")) +
							theme(axis.text.y = element_text(size=fs), axis.text.x = element_text(size=fs), axis.title.x = element_text(size=fs));
						suppressMessages(ggsave(file.path(odir,"ridgeplot.pdf"), width=10, height=round(log10(n+1)*8)));

						n <- nrow(gsea);
						fs <- round(10-max(10,n)/10+5);
						dotplot(gsea, title=paste0(domain," (",n,")"), showCategory=n, split=".sign", label_format=max(nchar(gsea$Description))+1) +
							facet_grid(.~.sign) +
							theme(axis.text.y = element_text(size=fs), axis.text.x = element_text(size=fs), axis.title.x = element_text(size=fs), strip.background = element_rect(linetype = 0, fill=NA)) +
							guides(color = guide_colourbar("FDR"));
						suppressMessages(ggsave(file.path(odir,"dotplot.full.pdf"), width=10, height=round(log10(n+1)*8)));

						ridgeplot(gsea, showCategory=n, label_format=max(nchar(gsea$Description))+1) + xlab("log2 FC distribution") +
							guides(fill = guide_colourbar("FDR")) + ggtitle(paste0(domain," (",n,")")) +
							theme(axis.text.y = element_text(size=fs), axis.text.x = element_text(size=fs), axis.title.x = element_text(size=fs));
						suppressMessages(ggsave(file.path(odir,"ridgeplot.full.pdf"), width=10, height=round(log10(n+1)*8)));

						df <- data.frame(fdr=-log(gsea$p.adjust, base=10), description = paste0(gsea$Description," (",gsea$setSize,")"));
						n <- min(20,nrow(df));
						dfp <- head(df,n);
						pdf(file.path(odir,"barplot.pdf"),width=max(nchar(df$description))/5,height=n/2+1.8);
						par(mar=c(5,max(nchar(df$description))/3+3,5,1));
						barplot(rev(dfp$fdr), main=paste0(domain," (Top ",n,")"), horiz=T, names.arg=rev(dfp$description),
							xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$fdr)+max(dfp$fdr)/10*2),
							cex.names=0.8, las=1);
						graphics.off();

						df <- data.frame(fdr=-log(gsea$p.adjust, base=10), description = paste0(gsea$Description," (",gsea$setSize,")"));
						n <- nrow(df);
						dfp <- df;
						pdf(file.path(odir,"barplot.full.pdf"),width=max(nchar(df$description))/5,height=n/2+1.8);
						par(mar=c(5,max(nchar(df$description))/3+3,5,1));
						barplot(rev(dfp$fdr), main=paste0(domain," (",n,")"), horiz=T, names.arg=rev(dfp$description),
							xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$fdr)+max(dfp$fdr)/10*2),
							cex.names=0.8, las=1);
						graphics.off();

						dir.create(file.path(odir,"gsea_plots"), recursive = T, showWarnings = F);
						for (i in 1:length(gsea$Description)){
							if(! is.na(gsea$Description[i])){
								gseaplot2(gsea, title = gsea$Description[i], geneSetID = i);
								suppressMessages(ggsave(file.path(odir,"gsea_plots",paste0(gsub("[^A-Za-z0-9]+","_",gsea$Description[i]),".pdf"))));
							};
						};

						gsea <- as.data.frame(gsea)[c("ID","setSize","pvalue","p.adjust","Description","core_enrichment")];
						gsea$core_enrichment <- gsub("/",";",gsea$core_enrichment);
					} else {
						gsea <- data.frame(matrix(ncol = 6, nrow = 0));
					};
				};
			};
			colnames(gsea) = c("id","count","pval","padj","description","enrichment");
			write.table(gsea, row.names = F, file = file.path(odir,"goenrichment.tsv"), quote=F, sep="\t");
		'
	CMD
		"$gmtfile" "$refterms" "$deseqtsv" "$outdir" "$domain"
	CMD

	return 0
}

function enrichment::_reducego(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-1 <cmds1>    | array of
			-g <oRgdb>    | path to R library directory containing custom/"my" oRgdb (org.My.eg.db) - see genome::mkgodb
			-d <domain>   | biological_process or cellular_component or molecular_function
			-i <orafile>  | path to
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory threads domain orgdb orafile outdir
	declare -n _cmds1_reduce
	while getopts '1:2:d:i:g:o:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_reduce=$OPTARG;;
			d)	((++mandatory)); domain=$OPTARG;;
			i)	((++mandatory)); orafile="$OPTARG";;
			g)	((++mandatory)); orgdb="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	# sum up counts of clustered terms?
	# dfp$count <- sapply(unique(df$parent), function(x) sum(df$count[df$parent==x]));
	commander::makecmd -a _cmds1_reduce -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			args <- commandArgs(TRUE);
			libdir <- args[1];
			ora <- args[2];
			odir <- args[3];
			domain <- args[4];
			domshort <- args[5];

			.libPaths(c(libdir,.libPaths()));
			suppressMessages(library(rrvgo));
			suppressMessages(library(ggplot2));

			df <- read.table(ora, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
			dfc <- df[,1:2];
			colnames(dfc) <- c("go","count");

			if(nrow(df)>2){
				m <- calculateSimMatrix(df$id, orgdb="org.My.eg.db", ont=domshort, method="Wang", keytype="GID");
				if(! is.na(m) && nrow(m)>2){
					df <- reduceSimMatrix(m, setNames(-log10(df$padj), df$id), threshold=0.7, orgdb="org.My.eg.db", keytype="GID");
					df <- merge(dfc, df, by="go");
					df <- df[order(df$score, decreasing=T), ];

					write.table(df, row.names = F, file = file.path(odir,"goenrichment_reduced.tsv"), quote=F, sep="\t");

					pdf(file.path(odir,"treemap.pdf"), width=16, height=10);
					treemapPlot(df, force.print.labels=T);
					graphics.off();

					p <- scatterPlot(m, df, labelSize=2);
					ggsave(file.path(odir,"semantic_space.pdf"), plot=p, width=16, height=10);

					dfp <- df[df$go %in% unique(df$parent),];
					dfp$term <- paste0(dfp$term," (",dfp$count,")");
					df <- dfp;

					n <- min(10,nrow(df));
					dfp <- head(df,n);
					pdf(file.path(odir,"barplot_reduced.pdf"),width=max(nchar(dfp$term))/5,height=n/2+1.8);
					par(mar=c(5,max(nchar(dfp$term))/3+3,5,1));
					barplot(rev(dfp$score), main=paste0(domain," (Top ",n,")"), horiz=T, names.arg=rev(dfp$term),
						xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$score+5)),
						cex.names=0.8, las=1);
					graphics.off();

					n <- nrow(df);
					dfp <- df;
					pdf(file.path(odir,"barplot_reduced.full.pdf"),width=max(nchar(dfp$term))/5,height=n/2+1.8);
					par(mar=c(5,max(nchar(dfp$term))/3+3,5,1));
					barplot(rev(dfp$score), main=domain, horiz=T, names.arg=rev(dfp$term),
						xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$score+5)),
						cex.names=0.8, las=1);
					graphics.off();
				};
			};
		'
	CMD
		"$orgdb" "$orafile" "$outdir" "$(echo $domain | sed -E 's/([^_]{1,3})[^_]*_(\S+)/\u\1. \u\2/')" "$(echo $domain | sed -E 's/([^_]{1,1})[^_]*_(\S{1,1}).+/\u\1\u\2/')"
	CMD

	return 0
}

function enrichment::_revigo(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-d <domain>   | biological_process or cellular_component or molecular_function
			-i <orafile>  | path to
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory threads domain gofile orafile outdir
	declare -n _cmds1_revigo _cmds2_revigo
	while getopts '1:2:d:i:o:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_revigo=$OPTARG;;
			2)	((++mandatory)); _cmds2_revigo=$OPTARG;;
			d)	((++mandatory)); domain="$OPTARG";;
			i)	((++mandatory)); orafile="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	# for pvalue instead of padj do revigo <(awk 'NR>1 && \$3<=0.05' $odir/gsea.tsv) --stdout
	# need to remove [_'"\t*$#%^!]
	# use own buffered reader to handle corner case where orafile is empty or revigo db does not contain any go term and thus instead of an error returns userValue_2 to userValue_XXXXXX in a single line
	commander::makecmd -a _cmds1_revigo -s ' ' -o "$outdir/revigo.tsv" -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- 'CMD'
		l=""; while read -N 1 c; do
			if [[ "$c" == $'\n' ]]; then
				[[ ! $pgid ]] && pgid=$l || echo "$l";
				l="";
			else
				l+="$c";
				if [[ "$l" =~ userVal ]]; then
					cut -f 1-6 <<< "$l";
					sleep 1;
					env kill -TERM -- -$pgid &> /dev/null || true;
					break;
				fi;
			fi;
		done < <(setsid --wait bash -c "echo \$((\$(ps -o pgid= -p \$\$))); revigo -Xmx1024m -XX:ParallelGCThreads=1 -XX:ConcGCThreads=1 <(awk 'NR>1 && \$4<=0.05 {print \$1\"\\t\"\$4}'
	CMD
		'$orafile')
	CMD
		--stdout" & wait $! 2> /dev/null || { e=$?; [[ $e -eq 15 || $e -eq 0 ]] && exit 0 || exit $e; }
	CMD
		) | perl -lane '
			next if $.<3;
			if($.>3){
				$F[2]=~s/(^[\W_]+|[\W_]+$)//g;
				$F[2]=~s/([^\w,.;:=+-\\(\)\\[\]\\{\}]|_)+/ /g;
				$F[2]=~s/\s+/ /g;
			}
			shift @F;
			print join("\t",@F);
		'
	CMD

	# don't to 10**(-df$log10_p.value), values will get too small and treemap will fail
	# don't kick not representative values
	# df$count needs to be present in goenrichment.tsv
	commander::makecmd -a _cmds2_revigo -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			suppressMessages(library("treemap"));
			args <- commandArgs(TRUE);
			gsea <- args[1];

			revigo <- args[2];
			outdir <- args[3];
			domain <- args[4];

			c <- read.table(gsea, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
			colnames(c)[1] <- "term_ID";
			colnames(c)[5] <- "enricher_description";
			df <- read.table(revigo, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");

			if(nrow(df)>0){
				dfp <- df;
				dfp$log10_p.value <- -(dfp$log10_p.value);
				dfp$representative <- df[sapply(dfp$representative,function(x) which(dfp$term_ID==x)),"description"];
				pdf(file.path(outdir,"treemap.pdf"), width=16, height=10);
				treemap(dfp, index = "representative", vSize = "log10_p.value", type = "categorical",
					vColor = "representative", title = domain, inflate.labels = T,
					lowerbound.cex.labels = 0, force.print.labels = T, position.legend = "none");
				graphics.off();

				df <- merge(df,c,by="term_ID");
				df$description <- paste0(df$description," (",df$count,")");
				df <- df[df$term_ID %in% df$representative,];
				df <- df[order(df$log10_p.value),];
				bars <- min(10,nrow(df));
				dfp <- head(df,bars);
				pdf(file.path(outdir,"barplot_revigo.pdf"),width=max(nchar(dfp$description))/5,height=bars/2+1.8);
				par(mar=c(5,max(nchar(dfp$description))/3+3,5,1));
				barplot(rev(-(dfp$log10_p.value)), main=paste0(domain," (Top ",bars,")"), horiz=T, names.arg=rev(dfp$description),
					xlab="-log10 p-value", col=rainbow(9)[1], xlim=c(0,max(-(dfp$log10_p.value)+5)),
					cex.names=0.8, las=1);
				graphics.off();

				bars <- nrow(df);
				dfp <- df;
				pdf(file.path(outdir,"barplot_revigo.full.pdf"),width=max(nchar(dfp$description))/5,height=bars/2+1.8);
				par(mar=c(5,max(nchar(dfp$description))/3+3,5,1));
				barplot(rev(-(dfp$log10_p.value)), main=domain, horiz=T, names.arg=rev(dfp$description),
					xlab="-log10 p-value", col=rainbow(9)[1], xlim=c(0,max(-(dfp$log10_p.value)+5)),
					cex.names=0.8, las=1);
				graphics.off();
			};
		'
	CMD
		"$orafile" "$outdir/revigo.tsv" "$outdir" "$domain"
	CMD

	return 0
}

function enrichment::test(){
	enrichment::go "$@"
}

function enrichment::go(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-g <gtfile>   | path to 4-column, tab separated, gene transposed file to generate domain/collection specific (column 3) gene matrix transposed (gmt) and term files
			                ..
			                ENSG00000199065 GO:0005615 cellular_component extracellular space
			                ENSG00000199065 GO:1903231 molecular_function mRNA binding involved in posttranscriptional gene silencing
			                ENSG00000199065 GO:0035195 biological_process gene silencing by miRNA
			                ..
			-r <mapper>   | array of bams within array of
			-i <deseqdir> | for gsea path to
			-j <countsdir>| of joined tpms for ora gmt background creation
			-c <cmpfiles> | array of
			-l <idfiles>  | for ora array of (does not require -i -r -c)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads gofile deseqdir countsdir tmpdir=${TMPDIR:-/tmp}
    declare -n _mapper_go _cmpfiles_go _idfiles_go
	while getopts 'S:s:t:r:c:g:l:i:j:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			r) _mapper_go=$OPTARG;;
			c) _cmpfiles_go=$OPTARG;;
			g) ((++mandatory)); gofile="$OPTARG";;
			l) _idfiles_go="$OPTARG";;
			i) deseqdir="$OPTARG";;
			j) countsdir="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	commander::printinfo "calculating go enrichment"

	# https://www.bioconductor.org/help/course-materials/2017/BioC2017/Day1/Workshops/OmicsData/doc/enrichOmics.html
	# ORA tests the overlap between DE genes (typically DE p-value < 0.05) and genes in a gene set based on the hypergeometric distribution.
	# A major limitation of ORA is that it restricts analysis to DE genes, excluding genes not satisfying the chosen significance threshold (typically the vast majority).
	# This is resolved by gene set enrichment analysis (GSEA), which scores the tendency of gene set members to appear rather at the top or bottom of the ranked list of all measured genes.
	# The statistical significance of the enrichment score (ES) of a gene set is assessed via sample permutation, i.e. (1) sample labels (= group assignment) are shuffled, (2) per-gene DE statistics are recomputed, and (3) the enrichment score is recomputed. Repeating this procedure many times allows to determine the empirical distribution of the enrichment score and to compare the observed enrichment score against it.
	# As GSEA’s permutation procedure involves re-computation of per-gene DE statistics, adaptations are necessary for RNA-seq. The EnrichmentBrowser implements an accordingly adapted version of GSEA, which allows incorporation of limma/voom, edgeR, or DESeq2 for repeated DE re-computation within GSEA.
	# While it might be in some cases necessary to apply permutation-based GSEA for RNA-seq data, there are also alternatives avoiding permutation. Among them is ROtAtion gene Set Testing (ROAST), which uses rotation instead of permutation

	declare -a cmd1 cmd2 cmd3 mapdata domains
	local m f i c t title domain odir orgdb="$gofile.oRgdb"
	[[ -e "$orgdb" ]] || orgdb="$(dirname "$gofile")/oRgdb"
	local tdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.go)"
	mapfile -t domains < <(cut -f 3 "$gofile" | grep . | sort -u)

	for domain in "${domains[@]}"; do
		commander::makecmd -a cmd1 -s '|' -o "$tdir/$domain.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			grep -Fw $domain "$gofile" | tee -i >(cut -f 2,4 | sort -u > "$tdir/$domain.terms")
		CMD
			perl -F'\t' -lane '
				push @{$m{$F[1]}},$F[0];
				END{
					print join"\t",($_,"NA",@{$m{$_}}) for keys %m
				}
			'
		CMD
	done

	for f in "${_idfiles_go[@]}"; do
		for domain in "${domains[@]}"; do
			[[ $domain =~ ^(cellular_component|molecular_function|biological_process)$ ]] && title="$(echo $domain | sed -E 's/([^_]{1,3})[^_]*_(\S+)/\u\1. \u\2/')" || title="$(echo $domain | sed -E 's/_+/ /g')"
			odir="$(dirname "$f")/$domain"
			if [[ $countsdir ]]; then
				enrichment::_ora -1 cmd2 -d "$title" -g "$tdir/$domain.gmt" -r "$tdir/$domain.terms" -i "$f" -j "$(find -L "$countsdir" -name "experiments.tpm" -print -quit | grep .)" -o "$odir"
			else
				enrichment::_ora -1 cmd2 -d "$title" -g "$tdir/$domain.gmt" -r "$tdir/$domain.terms" -i "$f" -o "$odir"
			fi
			[[ -e "$orgdb" && $domain =~ ^(cellular_component|molecular_function|biological_process)$ ]] && enrichment::_reducego -1 cmd3 -d $domain -i "$odir/goenrichment.tsv" -g "$orgdb" -o "$odir"
		done
	done

	for m in "${_mapper_go[@]}"; do
		for f in "${_cmpfiles_go[@]}"; do
			mapfile -t mapdata < <(perl -F'\t' -lane 'next if exists $m{$F[1]}; $m{$F[1]}=1; print $F[1]' "$f")
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
					IFS=$'\n'
					for deseqtsv in $(find -L "$deseqdir/$m/$c-vs-$t/" -type f -name "deseq.full.fcshrunk.tsv" | grep .); do
						unset IFS
						if [[ $(head "$deseqtsv" | wc -l) -gt 1 ]]; then
							for domain in "${domains[@]}"; do
								[[ $domain =~ ^(cellular_component|molecular_function|biological_process)$ ]] && title="$(echo $domain | sed -E 's/([^_]{1,3})[^_]*_(\S+)/\u\1. \u\2/')" || title="$(echo $domain | sed -E 's/_+/ /g')"
								odir="$(dirname "$deseqtsv")/$domain"
								enrichment::_gsea -1 cmd2 -d "$title" -g "$tdir/$domain.gmt" -r "$tdir/$domain.terms" -i "$deseqtsv" -o "$odir"
								[[ -e "$orgdb" && $domain =~ ^(cellular_component|molecular_function|biological_process)$ ]] && enrichment::_reducego -1 cmd3 -d $domain -i "$odir/goenrichment.tsv" -g "$orgdb" -o "$odir"
							done
						fi
					done
				done
			done
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
	fi

	return 0
}
