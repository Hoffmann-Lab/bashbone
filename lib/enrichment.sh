#! /usr/bin/env bash
# (c) Konstantin Riege

enrichment::_ora(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage:
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-d <domain>   | biological_process or cellular_component or molecular_function
			-g <gofile>   | path to
			-i <idsfile>  | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory threads domain gofile idsfile outdir
	declare -n _cmds1_ora _cmds2_ora
	while getopts '1:2:d:g:i:o:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_ora=$OPTARG;;
			2) ((mandatory++)); _cmds2_ora=$OPTARG;;
			d) ((mandatory++)); domain=$OPTARG;;
			g) ((mandatory++)); gofile="$OPTARG";;
			i) ((mandatory++)); idsfile="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done

	[[ $mandatory -lt 6 ]] && _usage && return 1

	commander::makecmd -a _cmds1_ora -s '|' -o "$outdir/reference.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
		grep -F $domain $gofile
	CMD
		perl -F'\t' -lane '
			push @{$m{$F[1]}},$F[0];
			END{
				print join"\t",($_,"NA",@{$m{$_}}) for keys %m
			}
		'
	CMD

	commander::makecmd -a _cmds2_ora -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			suppressMessages(library("clusterProfiler"));
			args <- commandArgs(TRUE);
			gmt <- args[1];
			idsfile <- args[2];
			odir <- args[3];
			
			genes <- scan(idsfile, character(), quote = "", quiet = T);
			ora <- data.frame(matrix(ncol = 4, nrow = 0));
			if(length(genes)>0){
				tg <- suppressMessages(read.gmt(gmt));
				ora <- enricher(genes, TERM2GENE=tg, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500);
				if(!is.null(ora)){
					ora <- as.data.frame(ora)[c("ID","Count","pvalue","p.adjust")];
				};
			};
			colnames(ora) = c("id","count","pval","padj");
			write.table(ora, row.names = F, file = file.path(odir,"goenrichment.tsv"), quote=F, sep="\t");
		'
	CMD
		"$outdir/reference.gmt" "$idsfile" "$outdir"
	CMD

	return 0
}

enrichment::_gsea(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage:
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-d <domain>   | biological_process or cellular_component or molecular_function
			-g <gofile>   | path to
			-i <deseqtsv> | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory threads domain gofile deseqtsv outdir
	declare -n _cmds1_gsea _cmds2_gsea
	while getopts '1:2:d:g:i:o:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_gsea=$OPTARG;;
			2) ((mandatory++)); _cmds2_gsea=$OPTARG;;
			d) ((mandatory++)); domain=$OPTARG;;
			g) ((mandatory++)); gofile="$OPTARG";;
			i) ((mandatory++)); deseqtsv="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done

	[[ $mandatory -lt 6 ]] && _usage && return 1

	commander::makecmd -a _cmds1_gsea -s '|' -o "$outdir/reference.gmt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
		grep -F $domain $gofile
	CMD
		perl -F'\t' -lane '
			push @{$m{$F[1]}},$F[0];
			END{
				print join"\t",($_,"NA",@{$m{$_}}) for keys %m
			}
		'
	CMD

	commander::makecmd -a _cmds2_gsea -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			suppressMessages(library("clusterProfiler"));
			args <- commandArgs(TRUE);
			gmt <- args[1];
			ddsr <- args[2];
			odir <- args[3];

			df <- read.table(ddsr, header=T, sep="\t", stringsAsFactors=F);
			df <- df[!is.na(df$log2FoldChange) , ];
			df <- df[!is.na(df$padj) , ];
			df <- df[df$baseMean > 0 , ];
			df <- df[abs(df$log2FoldChange)>=0.5 , ];
			df <- df[df$padj<=0.05 , ];

			gsea <- data.frame(matrix(ncol = 4, nrow = 0));
			if(nrow(df)>0){
				tg <- suppressMessages(read.gmt(gmt));
				gl <- abs(df[, which(colnames(df)=="log2FoldChange") ]);
				names(gl) <- as.character(df[, which(colnames(df)=="id")]);
				gl <- sort(gl, decreasing = T);
				gsea <- GSEA(gl, TERM2GENE=tg, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500);
				if(!is.null(gsea)){
					gsea <- as.data.frame(gsea)[c("ID","setSize","pvalue","p.adjust")];
				};
			};
			colnames(gsea) = c("id","count","pval","padj");
			write.table(gsea, row.names = F, file = file.path(odir,"goenrichment.tsv"), quote=F, sep="\t");
		'
	CMD
		"$outdir/reference.gmt" "$deseqtsv" "$outdir"
	CMD

	return 0
}

enrichment::_revigo(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage:
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-d <domain>   | biological_process or cellular_component or molecular_function
			-i <orafile>  | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory threads domain gofile orafile outdir
	declare -n _cmds1_revigo _cmds2_revigo
	while getopts '1:2:d:i:o:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_revigo=$OPTARG;;
			2) ((mandatory++)); _cmds2_revigo=$OPTARG;;
			d) ((mandatory++)); domain=$(sed -r 's/([^_]{1,3})[^_]*_(\S+)/\u\1.\u\2/' <<< $OPTARG);;
			i) ((mandatory++)); orafile="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done

	[[ $mandatory -lt 5 ]] && _usage && return 1

	# for pvalue instead of padj do revigo <(awk 'NR>1 && \$(NF-1)<=0.05' $odir/gsea.tsv) --stdout
	# need to remove [_'"\t*$#%^!]
	commander::makecmd -a _cmds1_revigo -s '|' -o "$outdir/revigo.tsv" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
		revigo <(awk 'NR>1 && \$NF<=0.05 {print \$1"\t"\$NF}' $orafile) --stdout
	CMD
		perl -lane '
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
	commander::makecmd -a _cmds2_revigo -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			suppressMessages(library("treemap"));
			args <- commandArgs(TRUE);
			revigo <- args[1];
			outf <- args[2];
			domain <- args[3];
			df <- read.table(revigo, header=T, sep="\t");
			df$log10_p.value = -(df$log10_p.value);
			df <- df[df$dispensability<0.7,];
			pdf(outf, width=16, height=9);
			treemap(df, index = "description", vSize = "log10_p.value", type = "categorical",
				vColor = "representative", title = domain, inflate.labels = T,
				lowerbound.cex.labels = 0, force.print.labels = T, position.legend = "none");
			graphics.off();
		'
	CMD
		"$outdir/revigo.tsv" "$outdir/treemap.pdf" $domain
	CMD

	# df$count needs to be present in goenrichment.tsv
	commander::makecmd -a _cmds2_revigo -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			args <- commandArgs(TRUE);
			gsea <- args[1];
			revigo <- args[2];
			outf <- args[3];
			domain <- args[4];
			minbars=10;
			c <- read.table(gsea,header=T,sep="\t");
			colnames(c)[1] <- "term_ID";
			df <- read.table(revigo,header=T,sep="\t");
			df <- merge(df,c,by="term_ID");
			df$description <- paste(df$description," ","(",df$count,")", sep="");
			df=df[order(df$log10_p.value),];
			df=head(df,min(minbars,nrow(df)));
			pdf(outf,width=max(nchar(df$description))/5,height=min(minbars,nrow(df))/2+1.8);
			par(mar=c(5,max(nchar(df$description))/3+3,5,1));
			barplot(rev(-(df$log10_p.value)), main=domain, horiz=T, names.arg=rev(df$description),
				xlab="-log10 p-value", col=rainbow(9)[1], xlim=c(0,max(-(df$log10_p.value)+5)),
				cex.names=0.8, las=1);
			graphics.off();
		'
	CMD
		"$orafile" "$outdir/revigo.tsv" "$outdir/barplot.pdf" $(echo $domain | sed -r 's/([^_]{1,3})[^_]*_(\S+)/\u\1.\u\2/')
	CMD

	return 0
}

enrichment::go(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-g <gofile>   | path to
			                ..
			                ENSG00000199065 GO:0005615 cellular_component
			                ENSG00000199065 GO:1903231 molecular_function
			                ENSG00000199065 GO:0035195 biological_process
			                ..
			-i <deseqdir> | for gsea (logfc>=1 filtered) path to
			-c <cmpfiles> | array of
			-l <idfiles>  | for ora array of (does not require -i -c -o)
			
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads gofile deseqdir genelist best=false
    declare -n _mapper_go _cmpfiles_go _idfiles_go
	while getopts 'S:s:t:r:c:g:l:i:o:b:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			r) _mapper_go=$OPTARG;;
			c) _cmpfiles_go=$OPTARG;;
			g) ((mandatory++)); gofile="$OPTARG";;
			l) _idfiles_go="$OPTARG";;
			i) deseqdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	( [[ $mandatory -lt 2 ]] || ( [[ ! $_idfiles_go ]] && [[ ! "$deseqdir" || ! $_cmpfiles_go ]] ) ) && _usage && return 1

	commander::print "calculating go enrichment"

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
			odir="$(dirname $f)/$domain"
			mkdir -p "$odir"
			enrichment::_ora -1 cmd1 -2 cmd2 -d $domain -g $gofile -i $f -o $odir
			enrichmentfiles+=("$odir/goenrichment.tsv")
		done
	done
		
	for m in "${_mapper_go[@]}"; do
		for f in "${_cmpfiles_go[@]}"; do
			mapfile -t mapdata < <(cut -d $'\t' -f 2 $f | uniq)
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do 
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do 
					deseqtsv="$deseqdir/$m/$c-vs-$t/deseq.tsv"
					for domain in biological_process cellular_component molecular_function; do
						odir="$deseqdir/$m/$c-vs-$t/$domain"
						mkdir -p "$odir"
						enrichment::_gsea -1 cmd1 -2 cmd2 -d $domain -g "$gofile" -i "$deseqtsv" -o "$odir"
						enrichmentfiles+=("$odir/goenrichment.tsv")
					done
				done
			done
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	commander::runcmd -v -b -t $threads -a cmd1 && \
			conda activate py2r && \
			commander::runcmd -v -b -t $threads -a cmd2 && \
			conda activate py2
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	declare -a cmd3 cmd4
	for f in "${enrichmentfiles[@]}"; do
		[[ $(head $f | wc -l ) -lt 2 ]] && commander::warn "no enriched go terms in $f" && continue
		odir=$(dirname $f)
		domain=$(basename $odir)
		enrichment::_revigo -1 cmd3 -2 cmd4 -d $domain -i $f -o $odir
	done

	$skip && {
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	} || {
		{	commander::runcmd -v -b -t $threads -a cmd3 && \
			conda activate py2r && \
			commander::runcmd -v -b -t $threads -a cmd4 && \
			conda activate py2
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
