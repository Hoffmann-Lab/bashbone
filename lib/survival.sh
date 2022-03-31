#! /usr/bin/env bash
# (c) Konstantin Riege

survival::gettcga() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-i <tcga-ids> | array of
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir
	declare -a ids;
	while getopts 'S:s:t:i:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			i)	declare -n _ids_getfpkm=$OPTARG
				ids=("${_ids_getfpkm[@]}");;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir/FPKM";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	commander::printinfo "downloading tcga datasets"
	[[ ! $ids ]] && ids=($(Rscript - <<< 'library(TCGAbiolinks); tcga.cdr=TCGAbiolinks:::getGDCprojects()$id; cat(tcga.cdr[grep("TCGA-",tcga.cdr)])'))

	# gdc_manifest.TCGA-HNSC.2020-07-08.txt
	# gdc-client download --manifest gdc_manifest.TCGA-HNSC.2020-07-08.txt --log-file HNSC/gdc.log -d HNSC/ -n 8 --retry-amount 10 --wait-time 5

	# one may map file_id and barcodes/case_id via
	# library(TCGAutils)
	# file2case = UUIDtoUUID(file.ids, to_type="case_id")

	# legacy = T means hg19
	# TCGAbiolinks:::getProjectSummary(project="TCGA-HNSC", legacy=T)
	# TCGAbiolinks:::getProjectSummary(project="TCGA-HNSC")

	# GDCdownload adds clinical infos
	# fetch them independently via GDCquery_clinic("TCGA-HNSC", "clinical")

	# see https://www.bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html
	# data.category = "Gene Expression Quantification" | "Isoform Expression Quantification" | "miRNA Expression Quantification"
	# workflow.type = "HTSeq - FPKM-UQ" | "HTSeq - FPKM" | "HTSeq - Counts"
	# "The upper quartile FPKM (FPKM-UQ) is a modified FPKM calculation in which the total protein-coding read count is replaced by the 75th percentile read count value for the sample."

	declare -a cmd1 cmd2
	local i
	for i in "${ids[@]}"; do
		commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
			Rscript - <<< '
				suppressMessages({library(TCGAbiolinks);
					library(SummarizedExperiment);
				});
				args = commandArgs(TRUE);
				odir = args[1];
				q <- GDCquery(project = args[2],
					data.category = "Transcriptome Profiling",
					data.type = "Gene Expression Quantification",
					workflow.type = "HTSeq - FPKM"
				);
				GDCdownload(q, directory=file.path(odir,"FPKM"));
				data = GDCprepare(q, save = T, save.filename=file.path(odir,"FPKM",paste0(type,".Rdata")), directory=file.path(odir,"FPKM"));
				map = getResults(q)[c("cases","file_id")];
				colnames(map) = c("barcode","file_id");
				map$barcode = str_replace_all(map$barcode, "\\W", ".");
				write.table(map, row.names = F, file=file.path(odir,"FPKM",paste0(type,".case2id.tsv")), quote=F, sep="\t");
			'
		CMD
			"$outdir" "$i"
		CMD
	done

	commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			suppressMessages({library(TCGAbiolinks);
				library(SummarizedExperiment);
				library(stringr);
			});
			args = commandArgs(TRUE);
			odir = args[1];
			clinical=data.frame();
			fpkm=data.frame();
			for (type in args[2:length(args)]){
				cat(paste0("working on ",type,"\n"));
				load(file.path(odir,"FPKM",paste0(type,".Rdata")));
				df = as.data.frame(colData(data));
				df = df[,sapply(df,is.vector)];
				df$barcode = str_replace_all(df$barcode, "\\W", ".");
				write.table(df, row.names = F, file=file.path(odir,"FPKM",paste0(type,".CLINICAL.full.tsv")), quote=F, sep="\t");
				write.table(df[c("barcode","vital_status","days_to_death","days_to_last_follow_up")], row.names = F, file=file.path(odir,"FPKM",paste0(type,".CLINICAL.tsv")), quote=F, sep="\t");
				if(nrow(clinical)==0){
					clinical = df;
				} else {
					cols = intersect(colnames(clinical), colnames(df));
					clinical = rbind(clinical[,cols], df[,cols]);
				};
				df = assay(data);
				df = data.frame(gene_id=rownames(df),df);
				write.table(df, row.names = F, file=file.path(odir,"FPKM",paste0(type,".FPKM.tsv")), quote=F, sep="\t");
				if(nrow(fpkm)==0){
					fpkm=df;
				}else{
					fpkm = merge(fpkm, df, by="gene_id");
				};
			};
			cat("saving..\n");
			write.table(clinical, row.names = F, file=file.path(odir,"CLINICAL.full.tsv"), quote=F, sep="\t");
			clinical = clinical[c("barcode","vital_status","days_to_death","days_to_last_follow_up")];
			write.table(clinical, row.names = F, file=file.path(odir,"CLINICAL.tsv"), quote=F, sep="\t");
			write.table(fpkm, row.names = F, file=file.path(odir,"FPKM.tsv"), quote=F, sep="\t");
		'
	CMD
		"$outdir" $(printf '"%s" ' "${ids[@]}")
	CMD

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -v -b -t $threads -a cmd1
		commander::runcmd -v -b -t $threads -a cmd2
	fi

	return 0
}

survival::ssgsea(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-f <fpkms>    | path to tsv or gct
			-c <clinical> | path to clinical file
			-l <idfiles>  | array of gmt files or id lists
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir counts clinical
	declare -n _sets_ssgsea
	while getopts 'S:s:t:f:l:c:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			f)	((++mandatory)); counts="$OPTARG";;
			c)	clinical="$OPTARG";;
			l)	((++mandatory)); _sets_ssgsea=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "performing ssgsea"

	local f n setsize
	rm -f "$outdir/input.gmt" "$outdir/gmt.list"
	for f in "${_sets_ssgsea[@]}"; do
		if [[ $(awk '{if(NF>i){i=NF}}END{print i}' "$f") -eq 1 ]]; then
			n="$(basename "$f" | sed -E 's/\s+/_/g')"
			n="${n%.*}"
			sed -E -e ':a;N;$!ba' -e 's/(^\s+|\s+$)//g' -e's/\s+/\t/g' -e "s/^/$n\tna\t/" "$f" >> "$outdir/input.gmt"
		else
			realpath -s "$f" >> "$outdir/gmt.list"
		fi
	done
	[[ -s "$outdir/input.gmt" ]] && realpath -s "$outdir/input.gmt" >> "$outdir/gmt.list"

	local gct="$outdir/input.gct"
	if [[ "$(head -1 "$counts")" =~ ^# ]]; then
		gct="$counts"
	else
		nrow=$(wc -l < "$counts")
		ncol=$(head -1 "$counts" | awk '{print NF}')
		awk -v nrow=$((nrow-1)) -v ncol=$((ncol-1)) -v OFS='\t' 'BEGIN{print "#1.2"; print nrow,ncol} {if(NR==1){$1="Name\tDescription"}else{$1=$1"\tna"}; print}' "$counts" > "$gct"
	fi

	# ssgsea v2 defaults: -w 0.75 -m 10 -n rank
	# -> without permutations -p 0, results are equal to gpmodule implementation, failes with -l $threads
	# use (combine.replace) instead of seperate up/down scores (combine.all)

	# commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
	# 	ssgsea-cli.R -p 0 -i "$outdir/input.gct" -d "$outdir/input.gmt" -o $outdir/ssgsea
	# 	{ echo "barcode"; cut -f 1 RFX7sets.gmt; } | xargs echo | sed 's/ /\t/g' > ssgsea.tsv
	# 	datamash transpose < <(tail -n +3 ssgsea-scores.gct) | grep -w -F -m 1 -A 999999 'No.columns.scored' | tail -n +2 >> ssgsea.tsv
	# CMD

	declare -a cmd1 cmd2 cmd3
	local minsize=$(awk -v i=10 '{if(NF-2<i){i=NF-2}}END{print i}' $(cat "$outdir/gmt.list"))
	commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
		Rscript "$(which ssGSEA.R)"
			"-l$(dirname "$(which ssGSEA.R)")"
			"-i$gct" "-o$outdir/ssgsea"
			"-D$outdir/gmt.list"
			"-nrank"
			"-v$minsize"
			"-w0.75"
			"-Ccombine.replace"
	CMD

	commander::makecmd -a cmd2 -s '|' -o "$outdir/ssgsea.tsv" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
		datamash transpose < <(tail -n +3 "$outdir/ssgsea.gct")
	CMD
		awk -v OFS='\t' 'NR==1 || NR>2 {if(NR==1){$1="barcode"} print}'
	CMD

	if [[ -s "$clinical" ]]; then
		mkdir -p "$outdir/plots"
		commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
			survival.R "$clinical" "$outdir/ssgsea.tsv" "$outdir/plots"
		CMD
	fi

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -t $threads -a cmd1
		commander::runcmd -v -b -t $threads -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
	fi

	return 0
}
