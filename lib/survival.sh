#! /usr/bin/env bash
# (c) Konstantin Riege

function survival::gettcga(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of (requiered)
			-i <tcga-ids> | array of (e.g. TCGA-READ or TCGA-UCS)
			-o <outdir>   | path to (requiered)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}"
	declare -a ids;
	while getopts 'S:s:t:i:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			i)	declare -n _ids_getfpkm=$OPTARG
				ids=("${_ids_getfpkm[@]}")
			;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
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
	# new as of april 2022: harmonized tables: workflow.type = "STAR - counts" which contains tpm, fpkm, counts, fpkm_uq
	# old: workflow.type = "HTSeq - FPKM-UQ" | "HTSeq - FPKM" | "HTSeq - Counts"
	# "The upper quartile FPKM (FPKM-UQ) is a modified FPKM calculation in which the total protein-coding read count is replaced by the 75th percentile read count value for the sample."

	# ssm$Gene=="ENSG00000141510"
	# df <- as.data.frame(assays(counts)$fpkm_unstrand);
	# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
	# to get rid of ENGS.._PAR_Y use !duplicated, counts are always 0. otherwise sum or mean or max via dplyr: group_by(df,gene_id) %>% summarise_all(sum)

	declare -a cmd1
	local i
	local tdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.tcga)"
	for i in "${ids[@]}"; do
		commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
			Rscript - <<< '
				suppressMessages({library(TCGAbiolinks);
					library(TCGAutils);
					library(SummarizedExperiment);
				});
				args <- commandArgs(TRUE);
				odir <- file.path(args[2],"TCGA");
				p <- args[3];

				setwd(args[1]);
				q <- GDCquery(project=p,
					data.category="Transcriptome Profiling",
					data.type="Gene Expression Quantification",
					workflow.type="STAR - Counts"
				);
				GDCdownload(q, directory=odir);
				counts <- GDCprepare(q, directory=odir);
				metadata <- getResults(q)[c("cases","sample.submitter_id","file_id")];
				colnames(metadata) <- c("barcode","sample","file_id.star_counts");

				q <- GDCquery(project=p,
					data.category="Simple Nucleotide Variation",
					data.type="Masked Somatic Mutation"
				);
				GDCdownload(q, directory=odir);
				ssm <- GDCprepare(q, directory=odir);
				df <- UUIDtoBarcode(getResults(q)$file_id, from_type="file_id");
				colnames(df) <- c("file_id.somatic_mutations", "barcode.somatic_mutations");
				df$sample <- sapply(strsplit(df$barcode.somatic_mutations, "-"), function(x){paste(x[1],x[2],x[3],x[4],sep="-")});
				df <- df[!duplicated(df$sample),];
				metadata <- merge(metadata, df, by="sample", all.x=T);

				df <- data.frame();
				barcodes <- ssm[ssm$Hugo_Symbol=="TP53" & ssm$Variant_Classification!="Silent",]$Tumor_Sample_Barcode;
				if(length(barcodes)>0){
					df <- data.frame(sample=unique(sapply(strsplit(barcodes, "-"), function(x){paste(x[1],x[2],x[3],x[4],sep="-")})), TP53_status="mutated");
				};
				barcodes <- ssm$Tumor_Sample_Barcode[! ssm$Tumor_Sample_Barcode %in% barcodes];
				if(length(barcodes)>0){
					df <- rbind(df, data.frame(sample=unique(sapply(strsplit(barcodes, "-"), function(x){paste(x[1],x[2],x[3],x[4],sep="-")})), TP53_status="wildtype"));
				};
				df <- df[!duplicated(df$sample),];
				metadata <- merge(metadata, df, by="sample", all.x=T);

				save(counts, ssm, metadata, file=file.path(odir,paste0(p,".Rdata")));

				clin <- as.data.frame(colData(counts));
				clin <- clin[,!sapply(clin, is.list)];
				colnames(clin)[1] <- "barcode";
				clin <- merge(clin, metadata[,2:ncol(metadata)], by="barcode", all.x=T);

				clin <- clin[order(clin$patient, clin$sample_type_id),];

				clinT <- clin[!grepl("Normal",clin$definition),];

				clinTuniq <- clinT[!clinT$sample %in% clinT$sample[duplicated(clinT$sample)],];
				clinTuniq <- clinTuniq[!clinTuniq$patient %in% clinTuniq$patient[duplicated(clinTuniq$patient)],];

				clinTp53 <- clinT[!clinT$sample %in% clinT$sample[duplicated(clinT$sample)],];
				clinTp53 <- clinT[!is.na(clinT$TP53_status),];
				clinTp53 <- clinTp53[!clinTp53$patient %in% clinTp53$patient[duplicated(clinTp53$patient)],];

				write.table(clin, row.names=F, file=file.path(odir,paste0(p,".CLINICAL.full.tsv")), quote=F, sep="\t");
				write.table(clinT, row.names=F, file=file.path(odir,paste0(p,".CLINICAL.tumor.tsv")), quote=F, sep="\t");
				write.table(clinTuniq, row.names=F, file=file.path(odir,paste0(p,".CLINICAL.tsv")), quote=F, sep="\t");
				write.table(clinTp53, row.names=F, file=file.path(odir,paste0(p,".CLINICAL.p53.tsv")), quote=F, sep="\t");

				df <- as.data.frame(assays(counts)$tpm_unstrand);
				df$gene_id <- sapply(strsplit(rownames(df),split="\\."),function(x){x[1]});
				df <- df[!duplicated(df$gene_id),];

				write.table(df, row.names=F, file=file.path(odir,paste0(p,".TPM.full.tsv")), quote=F, sep="\t");
				write.table(df[,colnames(df) %in% c("gene_id",clinT$barcode)], row.names=F, file=file.path(odir,paste0(p,".TPM.tumor.tsv")), quote=F, sep="\t");
				write.table(df[,colnames(df) %in% c("gene_id",clinTuniq$barcode)], row.names=F, file=file.path(odir,paste0(p,".TPM.tsv")), quote=F, sep="\t");
				write.table(df[,colnames(df) %in% c("gene_id",clinTp53$barcode)], row.names=F, file=file.path(odir,paste0(p,".TPM.p53.tsv")), quote=F, sep="\t");
			'
		CMD
			"$tdir" "$(realpath -se "$outdir")" "$i"
		CMD
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $threads -a cmd1
	fi

	return 0
}

function survival::ssgsea(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-f <counts>    | path to tsv or gct
			-c <clinical>  | path to clinical file
			-r <right>     | censoring. default: 1800
			-v <variables> | array of variables from clinical file to use as covariates for cox regression. e.g. age_at_diagnosis or gender or TP53_status
			-l <idfiles>   | array of gmt files or id lists
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false outdir counts clinical censored=1800
	declare -n _sets_ssgsea
	declare -a covariates
	while getopts 'S:s:f:l:c:r:o:v:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			f)	((++mandatory)); counts="$OPTARG";;
			c)	clinical="$OPTARG";;
			r)	censored=$OPTARG;;
			l)	((++mandatory)); _sets_ssgsea=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			v)	declare -n _covariates_ssgsea=$OPTARG
				covariates=("${_covariates_ssgsea[@]}")
			;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "performing ssgsea"

	local f n setsize
	rm -f "$outdir/input.gmt" "$outdir/gmt.list"
	for f in "${_sets_ssgsea[@]}"; do
		if [[ $(awk '{if(NF>i){i=NF}}END{print i}' "$f") -eq 1 ]]; then
			n="$(basename "$f" | sed -E 's/\s+/_/g')"
			n="${n%.*}"
			sed -E -e ':a;N;$!ba' -e 's/(^\s+|\s+$)//g' -e's/\s+/\t/g' -e "s/^/$n\tna\t/" "$f" >> "$outdir/input.gmt"
		else
			realpath -se "$f" >> "$outdir/gmt.list"
		fi
	done
	[[ -s "$outdir/input.gmt" ]] && realpath -se "$outdir/input.gmt" >> "$outdir/gmt.list"

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
	# use (combine.replace) instead of separate up/down scores (combine.all)

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
			survival.R "$clinical" $censored "$outdir/ssgsea.tsv" "$outdir/plots" ${covariates[@]}
		CMD
	fi

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i 1 -a cmd2
		commander::runcmd -v -b -i 1 -a cmd3
	fi

	return 0
}
