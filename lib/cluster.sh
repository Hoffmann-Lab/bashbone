#! /usr/bin/env bash
# (c) Konstantin Riege

function cluster::wgcna_deseq(){
	cluster::coexpression_deseq "$@"
}

function cluster::coexpression_deseq(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-M <maxmemory>| amount of
			-c <cmpfiles> | array of
			-l <idfiles>  | empty array to be filled with returned cluster gene lists
			-r <mapper>   | array of bams within array of
			-n <value>    | [01234] filter input for padj < 0.05 (0) log2foldchange > 0.5 (1) basemean > 5 (2) lower 30% percentile (3) tpm > 5 (4)
			-g <gtf>      | path to (optional)
			-b <biotype>  | within gtf
			-i <countsdir>| path to joined counts directory
			-j <deseqdir> | path to
			-o <outdir>   | path to
			-f <feature>  | feature (default: gene)
			-x <idfile>   | path to file containing ids of interest to extract from experiments files at counts directory. replaces -c -n -g -b -j
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}" deseqdir countsdir clusterfilter=NA biotype gtf feature="gene" extractids
	declare -n _mapper_coexpression _cmpfiles_coexpression _idfiles_coexpression
	while getopts 'S:s:n:b:g:t:M:c:r:i:j:o:l:f:x:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			n)	clusterfilter=$OPTARG;;
			b)	biotype="$OPTARG";;
			g)	gtf="$OPTARG";;
			t)	((++mandatory)); threads=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			c)	_cmpfiles_coexpression=$OPTARG;;
			r)	((++mandatory)); _mapper_coexpression=$OPTARG;;
			i)	((++mandatory)); countsdir="$OPTARG";;
			j)	deseqdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			l)	_idfiles_coexpression=$OPTARG;;
			f)	feature=$OPTARG;;
			x)	extractids="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage
	[[ $clusterfilter =~ 5 && $biotype && ! $gtf ]] && _usage
	local instances=$((${#_mapper_coexpression[@]}*2)) imemory
	read -r instances imemory < <(configure::memory_by_instances -i $instances -M "$maxmemory")

	declare -p _idfiles_coexpression | grep -q '=' || {
		unset _idfiles_coexpression
		declare -a _idfiles_coexpression
	}

	commander::printinfo "inferring coexpression"

	declare -a mapdata cmd1 tpmtojoin
	declare -A visited
	local m f i c t e odir cdir ddir params tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.join)"
	local tojoin="$tmp.tojoin" joined="$tmp.joined"
	echo "rm -rf '$tmp'*" >> "$BASHBONE_CLEANUP"

	for m in "${_mapper_coexpression[@]}"; do
		odir="$outdir/$m"
		cdir="$countsdir/$m"
		ddir="$deseqdir/$m"
		mkdir -p "$odir"

		if [[ ! -e "$extractids" ]]; then
			visited=()
			tpmtojoin=()
			rm -f "$joined"
			for f in "${_cmpfiles_coexpression[@]}"; do
				mapfile -t mapdata < <(perl -F'\t' -lane 'next if exists $m{$F[1]}; $m{$F[1]}=1; print $F[1]' "$f")
				i=0
				for c in "${mapdata[@]::${#mapdata[@]}-1}"; do
					for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
						awk 'NR>1' "$ddir/$c-vs-$t/deseq.full.tsv" | sort -k1,1V | cut -f 1,2,3,7 | sed -E 's/\s+NA(\s+|$)/\t0\1/g' > "$tojoin"
						if [[ -s "$joined" ]]; then
							join -t $'\t' "$joined" "$tojoin" > "$tmp"
							mv "$tmp" "$joined"
						else
							mv "$tojoin" "$joined"
						fi
					done
				done
			done
			mv "$joined" "$odir/experiments.deseq.tsv"

			# filter joined deseq tables requiers padj >0 due to NA replacement
			perl -F'\t' -slane '
				$okmean=0;
				$okfc=0;
				$okpval=0;
				for (my $i=1; $i<$#F; $i+=3){
					$okmean=1 if $F[$i] >= 5;
					$okfc=1 if exists $F[$i+4] && abs(abs($F[$i+1])-abs($F[$i+4]))>=0.5;
					$okpval=1 if $F[$i+2] > 0 && $F[$i+2] <= 0.05;
				}
				if ($cf=~/[012]{3,3}/){
					print $F[0] if $okmean+$okfc+$okpval == 3;
				} elsif ($cf=~/[01]{2,2}/) {
					print $F[0] if $okfc+$okpval == 2;
				} elsif ($cf==0) {
					print $F[0] if $okpval == 1;
				} elsif	($cf==1) {
					print $F[0] if $okfc == 1;
				} else {
					print $F[0];
				}
			' -- -cf=$clusterfilter "$odir/experiments.deseq.tsv" > "$odir/experiments.filtered.genes"
			# 0 padj 1 fc>0.5 2 basemean > 5 3 30% see below

			if [[ $clusterfilter =~ 4 ]]; then
				awk '{i=2; ok=0; while(i<NF){if($i>=5){ok=1} i=i+1} if(ok){print $1}}' "$cdir/experiments.tpm" | grep -Fw -f "$odir/experiments.filtered.genes" > "$tmp.filtered.genes"
				mv "$tmp.filtered.genes" "$odir/experiments.filtered.genes"
			fi

			[[ $clusterfilter =~ 5 && "$biotype" != "." ]] && {
				perl -F'\t' -slane '
					next if /^#/;
					if($#F<5){
						print $F[0] if $F[2] eq $cb;
					} else {
						$F[-1]=~/${ft}_(bio)?type\s+"([^"]+)/;
						if ($2 =~ /$cb/ && $F[-1]=~/${ft}_id\s+"([^"]+)/){
							print $1 unless exists $m{$1};
							$m{$1}=1;
						}
					}
				' -- -cb="$biotype" -ft="$feature" "$({ realpath -se "$gtf.info" "$gtf" 2> /dev/null || true; } | head -1 | grep .)" > "$tmp.genes"
				grep -Fw -f "$tmp.genes" "$odir/experiments.filtered.genes" > "$tmp.filtered.genes"
				mv "$tmp.filtered.genes" "$odir/experiments.filtered.genes"
			}

			extractids="$odir/experiments.filtered.genes"
		fi

		for e in tpm vsc; do
			head -1 "$cdir/experiments.$e" > "$odir/experiments.filtered.$e"
			grep -Fw -f "$extractids" "$cdir/experiments.$e" >> "$odir/experiments.filtered.$e"

			mkdir -p "$odir/$e"
			[[ "$e" == "tpm" ]] && params=TRUE || params=FALSE
			[[ $clusterfilter =~ 3 ]] && params+=" TRUE" || params+=" FALSE"
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				ulimit -s $(ulimit -Hs)
			CMD
				wgcna.R $((imemory/1024)) $params "$odir/experiments.filtered.$e" "$odir/$e" $([[ $deseqdir ]] && echo "'$ddir/experiments.csv'")
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $instances -a cmd1
	fi

	declare -a cmd2
	local type wdir cdir wdir odir
	for m in "${_mapper_coexpression[@]}"; do
		cdir="$countsdir/$m"
		for e in tpm vsc; do
			wdir="$outdir/$m/$e"

			[[ -s "$wdir/wgcna.cluster2modules" ]] || continue

			echo "$(head -1 "$cdir/experiments.$e") type" | tr ' ' '\t' > "$wdir/cluster.$e"
			echo "$(head -1 "$cdir/experiments.mean.$e") type" | tr ' ' '\t' > "$wdir/cluster.mean.$e"
			cp "$wdir/cluster.$e" "$wdir/cluster.$e.zscores"
			cp "$wdir/cluster.mean.$e" "$wdir/cluster.mean.$e.zscores"

			cp "$wdir/cluster.$e" "$wdir/cluster.top.$e"
			cp "$wdir/cluster.$e" "$wdir/cluster.top.$e.zscores"
			cp "$wdir/cluster.mean.$e" "$wdir/cluster.top.mean.$e"
			cp "$wdir/cluster.mean.$e" "$wdir/cluster.top.mean.$e.zscores"

			for i in $(cut -d ' ' -f 1 "$wdir/wgcna.cluster2modules" | sort -n); do
				type="Cluster.$(printf '%03d' $i)"
				odir="$wdir/$type"
				mkdir -p "$odir"

				awk -F '\t' -v i=$type '$2==i{print $1}' "$wdir/wgcna.cluster.tsv" > "$odir/genes.list"
				awk -F '\t' -v i=$type '$2==i{print $1}' "$wdir/wgcna.cluster.top.tsv" > "$odir/genes.top.list"
				[[ -s "$odir/genes.list" ]] || continue

				# add to array for later go enrichment
				_idfiles_coexpression+=("$odir/genes.list")

				grep -Fw -f "$odir/genes.list" "$cdir/experiments.$e" | sed "s/$/\t$type/" >> "$wdir/cluster.$e"
				grep -Fw -f "$odir/genes.list" "$cdir/experiments.mean.$e" | sed "s/$/\t$type/" >> "$wdir/cluster.mean.$e"
				grep -Fw -f "$odir/genes.list" "$cdir/experiments.$e.zscores" | sed "s/$/\t$type/" >> "$wdir/cluster.$e.zscores"
				grep -Fw -f "$odir/genes.list" "$cdir/experiments.mean.$e.zscores" | sed "s/$/\t$type/" >> "$wdir/cluster.mean.$e.zscores"

				grep -Fw -f "$odir/genes.top.list" "$cdir/experiments.$e" | sed "s/$/\t$type/" >> "$wdir/cluster.top.$e"
				grep -Fw -f "$odir/genes.top.list" "$cdir/experiments.mean.$e" | sed "s/$/\t$type/" >> "$wdir/cluster.top.mean.$e"
				grep -Fw -f "$odir/genes.top.list" "$cdir/experiments.$e.zscores" | sed "s/$/\t$type/" >> "$wdir/cluster.top.$e.zscores"
				grep -Fw -f "$odir/genes.top.list" "$cdir/experiments.mean.$e.zscores" | sed "s/$/\t$type/" >> "$wdir/cluster.top.mean.$e.zscores"
			done

			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster ${e^^*} "$wdir/cluster.$e" "$wdir/cluster.$e" "$wdir/wgcna.cluster.tsv"
			CMD
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster ${e^^*} "$wdir/cluster.mean.$e" "$wdir/cluster.mean.$e" "$wdir/wgcna.cluster.tsv"
			CMD
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster Z-Score "$wdir/cluster.$e.zscores" "$wdir/cluster.$e.zscores" "$wdir/wgcna.cluster.tsv"
			CMD
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster Z-Score "$wdir/cluster.mean.$e.zscores" "$wdir/cluster.mean.$e.zscores" "$wdir/wgcna.cluster.tsv"
			CMD

			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster ${e^^*} "$wdir/cluster.top.$e" "$wdir/cluster.top.$e" "$wdir/wgcna.cluster.tsv"
			CMD
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster ${e^^*} "$wdir/cluster.top.mean.$e" "$wdir/cluster.top.mean.$e" "$wdir/wgcna.cluster.tsv"
			CMD
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster Z-Score "$wdir/cluster.top.$e.zscores" "$wdir/cluster.top.$e.zscores" "$wdir/wgcna.cluster.tsv"
			CMD
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster Z-Score "$wdir/cluster.top.mean.$e.zscores" "$wdir/cluster.top.mean.$e.zscores" "$wdir/wgcna.cluster.tsv"
			CMD

		done
	done

	if $skip; then
		commander::printcmd -a cmd2
	else
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function cluster::wgcna(){
	cluster::coexpression "$@"
}

function cluster::coexpression(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-M <maxmemory>| amount of
			-l <idfiles>  | empty array to be filled with returned cluster gene lists
			-r <mapper>   | array of bams within array of
			-n <value>    | [34] filter input for lower 30% percentile (3) tpm > 5 (4)
			-g <gtf>      | path to (optional)
			-b <biotype>  | within gtf
			-i <countsdir>| path to
			-o <outdir>   | path to
			-f <feature>  | feature (default: gene)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory outdir tmpdir="${TMPDIR:-/tmp}" countsdir biotype gtf clusterfilter=NA feature="gene"
	declare -n _mapper_coexpression _idfiles_coexpression
	while getopts 'S:s:n:b:g:t:M:c:r:i:j:o:l:f:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			n)	clusterfilter=$OPTARG;;
			b)	biotype="$OPTARG";;
			g)	gtf="$OPTARG";;
			t)	((++mandatory)); threads=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			r)	((++mandatory)); _mapper_coexpression=$OPTARG;;
			i)	((++mandatory)); countsdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			l)	_idfiles_coexpression=$OPTARG;;
			f)	feature=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage && return 1
	[[ $clusterfilter =~ 5 && $biotype && ! $gtf ]] && _usage
	local instances=$((${#_mapper_coexpression[@]}*2)) imemory
	read -r instances imemory < <(configure::memory_by_instances -i $instances -M "$maxmemory")

	declare -p _idfiles_coexpression | grep -q '=' || {
		unset _idfiles_coexpression
		declare -a _idfiles_coexpression
	}

	commander::printinfo "inferring coexpression"

	declare -a cmd1x cmd2 cmd3 tojoin tdirs
	local m f e odir suff header sample countfile csv c
	for m in "${_mapper_coexpression[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"

		declare -n _bams_coexpression=$m
		suff=$(cat <(echo -e "${_bams_coexpression[0]}\n${_bams_coexpression[1]}") | rev | paste - - | sed -E 's/(.+\.).+\t\1.+/\1/' | rev)

		# tpm needs to be first for csv!
		for e in tpm vsc; do
			if [[ "$e" == "tpm" ]]; then
				csv="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.csv)"
				header="id"
				tojoin=()
				for f in "${_bams_coexpression[@]}"; do
					sample=$(basename $f $suff)
					tojoin+=("$(find -L "$countsdir/$m" -maxdepth 1 -name "$sample*.${feature}counts.htsc.$e" -print -quit | grep .)")
					header+="\t$sample"
					echo "$sample,${tojoin[-1]}" >> "$csv"
				done

				for r in $(seq 0 $(echo ${#tojoin[@]} | awk '{h=log($1+1)/log(2); h=h>int(h)?int(h)+1:h; print h-1}')); do
					declare -a cmd1x$r
					cmd1x[$r]=cmd1x$r
				done
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.join)")
				helper::multijoin \
					-1 cmd1x \
					-p "${tdirs[-1]}" \
					-e 0 \
					-h "$(echo -e "$header")" \
					-o "$odir/experiments.$e" \
					"${tojoin[@]}"

			else
				commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
					Rscript - <<< '
						options(warn=-1);
						suppressMessages(library("DESeq2"));
						suppressMessages(library("argparser"));
						args = arg_parser("Calculate variance stabilized counts", hide.opts=T);
						args = add_argument(args, "--csv", short="-c", help="path to input csv", flag=F);
						args = add_argument(args, "--out", short="-o", help="path to output tsv", flag=F);
						args = parse_args(args);
						dds = DESeqDataSetFromHTSeqCount(sampleTable=read.table(args$csv, header=F, sep=",", stringsAsFactors=F, check.names=F, quote=""), directory="", design=~1);
						vsd = as.data.frame(assay(varianceStabilizingTransformation(dds, blind=TRUE)));
						write.table(data.frame(id=rownames(vsd), vsd, check.names=F), file=args$out, col.names=T, row.names=F, quote=F, sep = '\t');
					'
				CMD
					-c "$csv" -o "$odir/experiments.$e"
				CMD
			fi

			commander::makecmd -a cmd3 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				Rscript - <<< '
					args <- commandArgs(TRUE);
					intsv <- args[1];
					outf <- args[2];
					df <- read.table(intsv, row.names=1, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
					df <- log(df+1);
					df <- df-rowMeans(df);
					df <- df/apply(df,1,sd);
					df[is.na(df)] <- 0;
					write.table(data.frame(id=rownames(df),df,check.names=F), row.names = F, file = outf, quote=F, sep="\t");
				'
			CMD
				"$odir/experiments.$e" "$odir/experiments.$e.zscores"
			CMD
		done
	done

	if $skip; then
		for c in "${cmd1x[@]}"; do
			commander::printcmd -a $c
		done
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		for c in "${cmd1x[@]}"; do
			commander::runcmd -v -b -i $threads -a $c
		done
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
	fi

	declare -a cmd3
	local params tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.filter)"
	echo "rm -rf '$tmp'*" >> "$BASHBONE_CLEANUP"

	for m in "${_mapper_coexpression[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"

		if [[ $clusterfilter =~ 4 ]]; then
			awk 'NR>1{i=2; ok=0; while(i<NF){if($i>=5){ok=1} i=i+1} if(ok){print $1}}' "$odir/experiments.tpm" > "$odir/experiments.filtered.genes"
		else
			awk 'NR>1{print $1}' "$odir/experiments.tpm" > "$odir/experiments.filtered.genes"
		fi

		[[ $$clusterfilter =~ 5 && "$biotype" != "." ]] && {
			perl -F'\t' -slane '
				next if /^#/;
				if($#F<5){
					print $F[0] if $F[2] eq $cb;
				} else {
					$F[-1]=~/${ft}_(bio)?type\s+"([^"]+)/;
					if ($2 =~ /$cb/ && $F[-1]=~/${ft}_id\s+"([^"]+)/){
						print $1 unless exists $m{$1};
						$m{$1}=1;
					}
				}
			' -- -cb="$biotype" -ft="$feature" "$({ realpath -se "$gtf.info" "$gtf" 2> /dev/null || true; } | head -1 | grep .)" > "$tmp.genes"
			grep -Fw -f "$tmp.genes" "$odir/experiments.filtered.genes" > "$tmp.filtered.genes"
			mv "$tmp.filtered.genes" "$odir/experiments.filtered.genes"
		}

		for e in tpm vsc; do
			head -1 "$odir/experiments.e" > "$odir/experiments.filtered.e"
			grep -Fw -f "$odir/experiments.filtered.genes" "$odir/experiments.e" >> "$odir/experiments.filtered.e"

			mkdir -p "$odir/$e"
			[[ "$e" == "tpm" ]] && params=TRUE || params=FALSE
			[[ $clusterfilter =~ 3 ]] && params+=" TRUE" || params+=" FALSE"
			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				ulimit -s $(ulimit -Hs)
			CMD
				wgcna.R $((imemory/1024)) $params "$odir/experiments.filtered.$e" "$odir/$e"
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -i $instances -a cmd3
	fi

	declare -a cmd4
	local i j wdir odir type
	for m in "${_mapper_coexpression[@]}"; do
		wdir="$outdir/$m"
		for e in tpm vsc; do

			[[ -s "$wdir/wgcna.cluster2modules" ]] || continue

			echo "$(head -1 "$wdir/experiments.$e") type" | tr ' ' '\t' > "$wdir/cluster.$e"
			cp "$wdir/cluster.$e" "$wdir/cluster.$e.zscores"

			cp "$wdir/cluster.$e" "$wdir/cluster.top.$e"
			cp "$wdir/cluster.$e" "$wdir/cluster.top.$e.zscores"

			for i in $(cut -d ' ' -f 1 "$wdir/wgcna.cluster2modules" | sort -n); do
				type="Cluster.$(printf '%03d' $i)"
				odir="$wdir/$type"
				mkdir -p "$odir"

				awk -F '\t' -v i=$type '$2==i{print $1}' "$wdir/wgcna.cluster.tsv" > "$odir/genes.list"
				awk -F '\t' -v i=$type '$2==i{print $1}' "$wdir/wgcna.cluster.top.tsv" > "$odir/genes.top.list"
				[[ -s "$odir/genes.list" ]] || continue

				# add to array for later go enrichment
				_idfiles_coexpression+=("$odir/genes.list")

				grep -Fw -f "$odir/genes.list" "$wdir/experiments.$e" | sed "s/$/\t$type/" >> "$wdir/cluster.$e"
				grep -Fw -f "$odir/genes.list" "$wdir/experiments.$e.zscores" | sed "s/$/\t$type/" >> "$wdir/cluster.$e.zscores"

				grep -Fw -f "$odir/genes.top.list" "$wdir/experiments.$e" | sed "s/$/\t$type/" >> "$wdir/cluster.top.$e"
				grep -Fw -f "$odir/genes.top.list" "$wdir/experiments.$e.zscores" | sed "s/$/\t$type/" >> "$wdir/cluster.top.$e.zscores"
			done

			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster ${e^^*} "$wdir/cluster.$e" "$wdir/cluster.$e" "$wdir/wgcna.cluster.tsv"
			CMD
			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster Z-Score "$wdir/cluster.$e.zscores" "$wdir/cluster.$e.zscores" "$wdir/wgcna.cluster.tsv"
			CMD
			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster ${e^^*} "$wdir/cluster.top.$e" "$wdir/cluster.top.$e" "$wdir/wgcna.cluster.tsv"
			CMD
			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster Z-Score "$wdir/cluster.top.$e.zscores" "$wdir/cluster.top.$e.zscores" "$wdir/wgcna.cluster.tsv"
			CMD

		done
	done

	if $skip; then
		commander::printcmd -a cmd4
	else
		commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}
