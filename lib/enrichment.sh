#! /usr/bin/env bash
# (c) Konstantin Riege

function enrichment::_ora(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-d <domain>   | biological_process or cellular_component or molecular_function
			-g <gofile>   | path to
			-i <idsfile>  | path to
			-j <tpmfile>  | tpm experiments file path to
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory threads domain gofile idsfile tpmtsv outdir
	declare -n _cmds1_ora _cmds2_ora
	while getopts '1:2:d:g:i:j:o:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_ora=$OPTARG;;
			2)	((++mandatory)); _cmds2_ora=$OPTARG;;
			d)	((++mandatory)); domain=$OPTARG;;
			g)	((++mandatory)); gofile="$OPTARG";;
			i)	((++mandatory)); idsfile="$OPTARG";;
			j)	tpmtsv="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage

	if [[ $tpmtsv ]]; then
		commander::makecmd -a _cmds1_ora -s '|' -o "$outdir/reference.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD'
			grep -Fw $domain "$gofile" | tee -i >(cut -f 2,4 | grep '^GO:' | sort -u > "$outdir/reference.terms")
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
		commander::makecmd -a _cmds1_ora -s '|' -o "$outdir/reference.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			grep -Fw $domain "$gofile" | tee -i >(cut -f 2,4 | grep '^GO:' | sort -u > "$outdir/reference.terms")
		CMD
			perl -F'\t' -lane '
				push @{$m{$F[1]}},$F[0];
				END{
					print join"\t",($_,"NA",@{$m{$_}}) for keys %m
				}
			'
		CMD
	fi

	commander::makecmd -a _cmds2_ora -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			suppressMessages(library("clusterProfiler"));
			suppressMessages(library("DOSE"));
			suppressMessages(library("enrichplot"));
			suppressMessages(library("ggplot2"));
			args <- commandArgs(TRUE);
			gmt <- args[1];
			g2n <- args[2];
			idsfile <- args[3];
			odir <- args[4];
			domain <- args[5];

			genes <- scan(idsfile, character(), quote="", quiet=T);
			ora <- data.frame(matrix(ncol = 6, nrow = 0));
			if(length(genes)>0){
				tg <- suppressMessages(read.gmt(gmt));
				tn <- read.table(g2n, sep="\t", stringsAsFactors=F, check.names=F, quote="");
				if(sum(genes %in% tg$gene)>0){
					ora <- enricher(genes, TERM2GENE=tg, TERM2NAME=tn, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500);
					if(!is.null(ora) && nrow(ora)>0){
						dotplot(ora, showCategory=10, font.size=10) +
							theme(strip.background = element_rect(linetype = 0, fill=NA)) +
							guides(color = guide_colourbar("FDR"));
						suppressMessages(ggsave(file.path(odir,"dotplot.pdf")));

						df <- data.frame(fdr=-log(ora$p.adjust, base=10), description = paste0(ora$Description," (",ora$Count,")"));
						bars <- min(10,nrow(df));
						dfp <- head(df,bars);
						pdf(file.path(odir,"barplot.pdf"),width=max(nchar(df$description))/5,height=bars/2+1.8);
						par(mar=c(5,max(nchar(df$description))/3+3,5,1));
						barplot(rev(dfp$fdr), main=paste0(domain," (Top ",bars,")"), horiz=T, names.arg=rev(dfp$description),
							xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$fdr)+max(dfp$fdr)/10*2),
							cex.names=0.8, las=1);
						graphics.off();

						bars <- nrow(df);
						dfp <- df;
						pdf(file.path(odir,"barplot.full.pdf"),width=max(nchar(df$description))/5,height=bars/2+1.8);
						par(mar=c(5,max(nchar(df$description))/3+3,5,1));
						barplot(rev(dfp$fdr), main=paste0(domain," (",bars,")"), horiz=T, names.arg=rev(dfp$description),
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
		"$outdir/reference.gmt" "$outdir/reference.terms" "$idsfile" "$outdir" $(echo $domain | sed -E 's/([^_]{1,3})[^_]*_(\S+)/\u\1.\u\2/')
	CMD

	return 0
}

function enrichment::_gsea(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-d <domain>   | biological_process or cellular_component or molecular_function
			-g <gofile>   | path to
			-i <deseqtsv> | path to
			-j <tpmfile>  | tpm experiments file path to
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory threads domain gofile deseqtsv tpmtsv outdir
	declare -n _cmds1_gsea _cmds2_gsea
	while getopts '1:2:d:g:i:j:o:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_gsea=$OPTARG;;
			2)	((++mandatory)); _cmds2_gsea=$OPTARG;;
			d)	((++mandatory)); domain=$OPTARG;;
			g)	((++mandatory)); gofile="$OPTARG";;
			i)	((++mandatory)); deseqtsv="$OPTARG";;
			j)	tpmtsv="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage

	if [[ $tpmtsv ]]; then
		commander::makecmd -a _cmds1_gsea -s '|' -o "$outdir/reference.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD'
			grep -Fw $domain "$gofile" | tee -i >(cut -f 2,4 | grep '^GO:' | sort -u > "$outdir/reference.terms")
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
		commander::makecmd -a _cmds1_gsea -s '|' -o "$outdir/reference.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			grep -F $domain "$gofile" | tee -i >(cut -f 2,4 | grep '^GO:' | sort -u > "$outdir/reference.terms")
		CMD
			perl -F'\t' -lane '
				push @{$m{$F[1]}},$F[0];
				END{
					print join"\t",($_,"NA",@{$m{$_}}) for keys %m
				}
			'
		CMD
	fi

	commander::makecmd -a _cmds2_gsea -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
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
						dotplot(gsea, showCategory=10, split=".sign", font.size=10) +
							facet_grid(.~.sign) +
							theme(axis.text.y = element_text(size=8), strip.background = element_rect(linetype = 0, fill=NA)) +
							guides(color = guide_colourbar("FDR"));
						suppressMessages(ggsave(file.path(odir,"dotplot.pdf")));

						ridgeplot(gsea) + xlab("log2 FC distribution") +
							guides(fill = guide_colourbar("FDR")) +
							theme(axis.text.y = element_text(size=8), axis.text.x = element_text(size=10), axis.title.x = element_text(size=10));
						suppressMessages(ggsave(file.path(odir,"ridgeplot.pdf")));

						df <- data.frame(fdr=-log(gsea$p.adjust, base=10), description = paste0(gsea$Description," (",gsea$setSize,")"));
						bars <- min(10,nrow(df));
						dfp <- head(df,bars);
						pdf(file.path(odir,"barplot.pdf"),width=max(nchar(df$description))/5,height=bars/2+1.8);
						par(mar=c(5,max(nchar(df$description))/3+3,5,1));
						barplot(rev(dfp$fdr), main=paste0(domain," (Top ",bars,")"), horiz=T, names.arg=rev(dfp$description),
							xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$fdr)+max(dfp$fdr)/10*2),
							cex.names=0.8, las=1);
						graphics.off();

						df <- data.frame(fdr=-log(gsea$p.adjust, base=10), description = paste0(gsea$Description," (",gsea$setSize,")"));
						bars <- nrow(df);
						dfp <- df;
						pdf(file.path(odir,"barplot.full.pdf"),width=max(nchar(df$description))/5,height=bars/2+1.8);
						par(mar=c(5,max(nchar(df$description))/3+3,5,1));
						barplot(rev(dfp$fdr), main=paste0(domain," (",bars,")"), horiz=T, names.arg=rev(dfp$description),
							xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$fdr)+max(dfp$fdr)/10*2),
							cex.names=0.8, las=1);
						graphics.off();

						dir.create(file.path(odir,"gsea_plots"), recursive = T, showWarnings = F);
						for (i in 1:length(gsea$Description)){
							if(! is.na(gsea$Description[i])){
								gseaplot2(gsea, title = gsea$Description[i], geneSetID = i);
								suppressMessages(ggsave(file.path(odir,"gsea_plots",paste0(gsea$ID[i],".pdf"))));
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
		"$outdir/reference.gmt" "$outdir/reference.terms" "$deseqtsv" "$outdir" $(echo $domain | sed -E 's/([^_]{1,3})[^_]*_(\S+)/\u\1.\u\2/')
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

					pdf(file.path(odir,"treemap.pdf"), width=16, height=9);
					treemapPlot(df, force.print.labels=T);
					graphics.off();

					p <- scatterPlot(m, df, labelSize=2);
					ggsave(file.path(odir,"semantic_space.pdf"), plot=p, width=16, height=9);

					dfp <- df[df$go %in% unique(df$parent),];
					dfp$term <- paste0(dfp$term," (",dfp$count,")");
					df <- dfp;

					bars <- min(10,nrow(df));
					dfp <- head(df,bars);
					pdf(file.path(odir,"barplot_reduced.pdf"),width=max(nchar(dfp$term))/5,height=bars/2+1.8);
					par(mar=c(5,max(nchar(dfp$term))/3+3,5,1));
					barplot(rev(dfp$score), main=paste0(domain," (Top ",bars,")"), horiz=T, names.arg=rev(dfp$term),
						xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$score+5)),
						cex.names=0.8, las=1);
					graphics.off();

					bars <- nrow(df);
					dfp <- df;
					pdf(file.path(odir,"barplot_reduced.full.pdf"),width=max(nchar(dfp$term))/5,height=bars/2+1.8);
					par(mar=c(5,max(nchar(dfp$term))/3+3,5,1));
					barplot(rev(dfp$score), main=domain, horiz=T, names.arg=rev(dfp$term),
						xlab="-log10 q-value", col=rainbow(9)[1], xlim=c(0,max(dfp$score+5)),
						cex.names=0.8, las=1);
					graphics.off();
				};
			};
		'
	CMD
		"$orgdb" "$orafile" "$outdir" $(echo $domain | sed -E 's/([^_]{1,3})[^_]*_(\S+)/\u\1.\u\2/') $(echo $domain | sed -E 's/([^_]{1,1})[^_]*_(\S{1,1}).+/\u\1\u\2/')
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
			d)	((++mandatory)); domain=$OPTARG;;
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
				pdf(file.path(outdir,"treemap.pdf"), width=16, height=9);
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
		"$orafile" "$outdir/revigo.tsv" "$outdir" $(echo $domain | sed -E 's/([^_]{1,3})[^_]*_(\S+)/\u\1.\u\2/')
	CMD

	return 0
}

function enrichment::go(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-g <gofile>   | path to 4-column tab separated file
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

	local OPTIND arg mandatory skip=false threads gofile deseqdir countsdir
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

	declare -a cmd1 cmd2 enrichmentfiles mapdata
	local m f i c t domain odir

	for f in "${_idfiles_go[@]}"; do
		for domain in biological_process cellular_component molecular_function; do
			odir="$(dirname "$f")/$domain"
			if [[ $countsdir ]]; then
				enrichment::_ora -1 cmd1 -2 cmd2 -d $domain -g "$gofile" -i "$f" -j "$(find -L "$countsdir" -name "experiments.tpm" -print -quit | grep .)" -o "$odir"
			else
				enrichment::_ora -1 cmd1 -2 cmd2 -d $domain -g "$gofile" -i "$f" -o "$odir"
			fi
			enrichmentfiles+=("$odir/goenrichment.tsv")
		done
	done

	for m in "${_mapper_go[@]}"; do
		for f in "${_cmpfiles_go[@]}"; do
			mapfile -t mapdata < <(perl -F'\t' -lane 'next if exists $m{$F[1]}; $m{$F[1]}=1; print $F[1]' "$f")
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
					# deseqtsv="$deseqdir/$m/$c-vs-$t/deseq.tsv"
					IFS=$'\n'
					for deseqtsv in $(find -L "$deseqdir/$m/$c-vs-$t/" -type f -name "deseq.tsv" | grep .); do
						unset IFS
						if [[ $(wc -l < "$deseqtsv") -gt 1 ]]; then
							for domain in biological_process cellular_component molecular_function; do
								# odir="$deseqdir/$m/$c-vs-$t/$domain"
								odir="$(dirname "$deseqtsv")/$domain"
								if [[ -e "$deseqdir/$m/$c-vs-$t/experiments.tpm" ]]; then # ensure compatibility with previous version
									enrichment::_gsea -1 cmd1 -2 cmd2 -d $domain -g "$gofile" -i "$deseqtsv" -j "$deseqdir/$m/$c-vs-$t/experiments.tpm" -o "$odir"
								else
									enrichment::_gsea -1 cmd1 -2 cmd2 -d $domain -g "$gofile" -i "$deseqtsv" -o "$odir"
								fi
								enrichmentfiles+=("$odir/goenrichment.tsv")
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
	else
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	local orgdb="$gofile.oRgdb"
	[[ -e "$orgdb" ]] || orgdb="$(dirname "$gofile")/oRgdb"
	if [[ -e "$orgdb" ]]; then
		declare -a cmd3 cmd4
		local x
		for f in "${enrichmentfiles[@]}"; do
			x=$(head "$f" | wc -l)
			[[ $x -lt 2 ]] && commander::warn "no enriched go terms in $f"
			odir="$(dirname "$f")"
			domain="$(basename "$odir")"
			#enrichment::_revigo -1 cmd3 -2 cmd4 -d $domain -i "$f" -o "$odir"
			enrichment::_reducego -1 cmd3 -d $domain -i "$f" -g "$orgdb" -o "$odir"
		done

		if $skip; then
			commander::printcmd -a cmd3
			# commander::printcmd -a cmd4
		else
			# local instances ithreads
			# read -r instances ithreads < <(configure::instances_by_memory -m 1024 -T $threads)
			# commander::runcmd -v -b -i $instances -a cmd3
			# commander::runcmd -v -b -i $threads -a cmd4
			commander::runcmd -v -b -i $threads -a cmd3
		fi
	fi

	return 0
}
