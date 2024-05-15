#! /usr/bin/env bash
# (c) Konstantin Riege

function expression::diego(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-5 <skip>       | true/false md5sums, gtf prep respectively
			-t <threads>    | number of
			-r <mapper>     | array of bams within array of
			-g <gtf>        | path to
			-c <cmpfiles>   | array of
			-f <fcthreshold>| value of (default: 1)
			-o <outdir>     | absolute path to
			-e <exonmode>   | true/false in addition to splice junction mode (true)
			-x <strandness> | hash per bam of for exonmode
			-i <htscdir>    | absolute output path to for exonmode
			-F              | force index
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads countsdir tmpdir="${TMPDIR:-/tmp}" outdir gtf fcth exonmode=true forceidx=false
	declare -n _mapper_diego _cmpfiles_diego _strandness_diego
	while getopts 'S:s:5:t:r:x:g:c:f:e:i:o:F' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_diego=$OPTARG;;
			x)	_strandness_diego=$OPTARG;;
			g)	((++mandatory)); gtf="$OPTARG";;
			c)	((++mandatory)); _cmpfiles_diego=$OPTARG;;
			f)	fcth=$OPTARG;;
			e)	exonmode=$OPTARG;;
			i)	countsdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			F)	forceidx=true;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage
	$exonmode && [[ ${#_strandness_diego[@]} -eq 0 ]] && _usage
	$exonmode && [[ ! $countsdir ]] && _usage

	commander::printinfo "differential splice junction analyses"

	if $skipmd5; then
		commander::warn "skip checking md5 sums and thus annotation preparation"
	else
		commander::printinfo "checking md5 sums"
		local thismd5gtf tmp
		declare -a cmdprep1 cmdprep2
		[[ -s "$gtf" ]] && thismd5gtf=$(md5sum "$gtf" | cut -d ' ' -f 1)
		if $forceidx || [[ ! -s "${gtf%.*}.diego.bed" ]] || [[ "$thismd5gtf" && "$thismd5gtf" != "$md5gtf" ]]; then
			commander::printinfo "preparing annotation for differential splice junction analyses"
			commander::makecmd -a cmdprep1 -s ';' -c {COMMANDER[0]}<<- CMD
				gfftoDIEGObed.pl -g "$gtf" -o "${gtf%.*}.diego.bed"
			CMD
		fi

		if $forceidx || [[ ! -s "${gtf%.*}.aggregated.gtf" ]] || [[ "$thismd5gtf" && "$thismd5gtf" != "$md5gtf" ]]; then
			commander::printinfo "preparing annotation for exon_id tag based quantification"
			tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.gtf)"
			commander::makecmd -a cmdprep2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				dexseq_prepare_annotation2.py
					-r no
					-f "$tmp"
					"$gtf" "${gtf%.*}.dexseq.gtf"
			CMD
				sed -E 's/(.+gene_id\s+")([^"]+)(.+exon_number\s+")([^"]+)(.+)/\1\2\3\4\5; exon_id "\2:\4"/' "$tmp" > "${gtf%.*}.aggregated.gtf"
			CMD
		fi

		commander::runcmd -c diego -v -b -i $threads -a cmdprep1
		commander::runcmd -c htseq -v -b -i $threads -a cmdprep2
	fi

	if $exonmode; then
		quantify::featurecounts \
			-S false \
			-s true \
			-t $threads \
			-p "$tmpdir" \
			-g "${gtf%.*}.aggregated.gtf" \
			-l exon \
			-f exon \
			-o "$countsdir" \
			-r _mapper_diego \
			-x _strandness_diego
	fi

	declare -a cmd1 cmd2 tdirs mapdata mappeddirs
	local m f i c t odir sjfile countfile min sample condition library replicate factors
	for m in "${_mapper_diego[@]}"; do
		declare -n _bams_diego=$m
		mapfile -t mappeddirs < <(printf '%s\n' "${_bams_diego[@]}" | xargs -I {} dirname {} | sort -Vu)

		for f in "${_cmpfiles_diego[@]}"; do
			mapfile -t mapdata < <(perl -F'\t' -lane 'next if exists $m{$F[1]}; $m{$F[1]}=1; print $F[1]' "$f")
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
					odir="$outdir/$m/$c-vs-$t"
					mkdir -p "$odir"
					rm -f "$odir/groups.tsv" "$odir/list.sj.tsv" "$odir/list.ex.tsv"
					unset sample condition library replicate factors
					while read -r sample condition library replicate factors; do
						#echo -e "$condition\t$sample.$replicate" >> "$odir/groups.tsv"
						echo -e "$condition\t$condition.$replicate" >> "$odir/groups.tsv"

						sjfile="$(find -L "${mappeddirs[@]}" -maxdepth 1 -name "$sample*.sj" -print -quit | grep .)"
						#echo -e "$sample.$replicate\t$sjfile" >> "$odir/list.sj.tsv"
						echo -e "$condition.$replicate\t$(realpath -se "$sjfile")" >> "$odir/list.sj.tsv"

						if $exonmode; then
							countfile="$(find -L "$countsdir/$m" -maxdepth 1 -name "$sample*.exoncounts.htsc" -print -quit | grep .)"
							#echo -e "$sample.$replicate\t$countfile" >> "$odir/list.ex.tsv"
							echo -e "$condition.$replicate\t$(realpath -se "$countfile")" >> "$odir/list.ex.tsv"
						fi

					done < <(awk -v c=$c '$2==c' "$f" | sort -k4,4V && awk -v t=$t '$2==t' "$f" | sort -k4,4V)

					min=$(cut -f 1 "$odir/groups.tsv" | sort | uniq -c | column -t | cut -d ' ' -f 1 | sort -n | head -1)
					if [[ -s "$odir/list.sj.tsv" ]]; then
						if [[ $m == "segemehl" ]]; then
							tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.diego)")
							commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
								cd "${tdirs[-1]}"
							CMD
								pre_segemehl.pl
									-l "$(realpath -se "$odir/list.sj.tsv")"
									-a "$(realpath -se "${gtf%.*}.diego.bed")"
									-o "input.sj.tsv"
							CMD
								mv input.sj.tsv "$(realpath -s "$odir/input.sj.tsv")"
							CMD
						elif [[ $m == "star" ]]; then
							tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.diego)")
							commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
								cd "${tdirs[-1]}"
							CMD
								pre_STAR.py
									-l "$(realpath -se "$odir/list.sj.tsv")"
									-d "$(realpath -se "${gtf%.*}.diego.bed")"
							CMD
								mv junction_table.txt "$(realpath -s "$odir/input.sj.tsv")"
							CMD
						fi
						commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
							diego.py
								-d $((min<3?min:3))
								-a "$odir/input.sj.tsv"
								-b "$odir/groups.tsv"
								-q 0.05
								-z ${fcth:=1}
								-x $c
								> "$odir/diego.sj.tsv"
						CMD
						commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
							diego.py
								-d $((min<3?min:3))
								-a "$odir/input.sj.tsv"
								-b "$odir/groups.tsv"
								-q 0.05
								-z ${fcth:=1}
								-x $c
								-e
								-f "$odir/dendrogram.sj"
						CMD
					fi

					if $exonmode; then
						tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.diego)")
						commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
							cd "${tdirs[-1]}"
						CMD
							HTseq2DIEGO.pl
								-i "$(realpath -se "$odir/list.ex.tsv")"
								-o "$(realpath -s "$odir/input.ex.tsv")"
						CMD
						commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
							diego.py
								-d $((min<3?min:3))
								-a "$odir/input.ex.tsv"
								-b "$odir/groups.tsv"
								-q 0.05
								-z ${fcth:=1}
								-x $c
								> "$odir/diego.ex.tsv"
						CMD
						commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
							diego.py
								-d $((min<3?min:3))
								-a "$odir/input.ex.tsv"
								-b "$odir/groups.tsv"
								-q 0.05
								-z ${fcth:=1}
								-x $c
								-e
								-f "$odir/dendrogram.ex"
						CMD
					fi
				done
			done
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c diego -v -b -i $threads -a cmd1
		commander::runcmd -c diego -v -b -i $threads -a cmd2
	fi

	return 0
}

function expression::deseq(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-g <gtf>      | path to
			-c <cmpfiles> | array of
			-i <countsdir>| path to
			-o <outdir>   | path to
			-f <feature>  | feature (default: gene)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads countsdir outdir gtf gtfinfo feature="gene"
	declare -n _mapper_deseq _cmpfiles_deseq
	while getopts 'S:s:t:r:g:c:i:o:f:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_deseq=$OPTARG;;
			g)	gtf="$OPTARG"; gtfinfo="$([[ -s "$gtf.info" ]] && echo "$gtf.info" || true)";;
			c)	((++mandatory)); _cmpfiles_deseq=$OPTARG;;
			i)	((++mandatory)); countsdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			f)	feature="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "principal component and differential expression analyses"

	local instances ithreads i=${#_mapper_deseq[@]}
	read -r instances ithreads < <(configure::instances_by_threads -i $i -t 10 -T $threads)
	[[ $ithreads -gt 10 ]] && ithreads=10 # further speedup questionable and to avoid running out of memory due to too many forks.

	declare -a cmd1 cmd2 cmd3 cmdanno cmd4 cmd5 cmd6 cmps mapdata tojoin
	declare -A visited
	local m f h e i c t odir countfile sample condition library replicate factors header meanheader
	for m in "${_mapper_deseq[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"
		visited=()
		cmps=()

		head -1 ${_cmpfiles_deseq[0]} | perl -F'\t' -lane 'push @o,"sample,countfile,condition,replicate"; for (4..$#F){push @o,"factor".($_-3)}END{print join",",@o}' > "$odir/experiments.csv"
		for f in "${_cmpfiles_deseq[@]}"; do
			mapfile -t mapdata < <(perl -F'\t' -lane 'next if exists $m{$F[1]}; $m{$F[1]}=1; print $F[1]' "$f")
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
					cmps+=("$c $t")
					unset sample condition library replicate factors
					tojoin=()
					header="id"
					meanheader="id"

					mkdir -p "$odir/$c-vs-$t"
					head -1 "$odir/experiments.csv" > "$odir/$c-vs-$t/experiments.csv"
					while read -r sample condition library replicate factors; do
						countfile="$(find -L "$countsdir/$m" -maxdepth 1 -name "$sample*.${feature}counts.htsc" -print -quit | grep .)"
						[[ $factors ]] && factors=","$(echo $factors | sed -E 's/\s+/,/g')

						# echo "$sample,$countfile,$condition,$replicate$factors" >> "$odir/$c-vs-$t/experiments.csv"
						echo "$condition.$replicate,$countfile,$condition,$replicate$factors" >> "$odir/$c-vs-$t/experiments.csv"

						tojoin+=("$(find -L "$countsdir/$m" -maxdepth 1 -name "$sample*.${feature}counts.htsc.tpm" -print -quit | grep .)")
						#header+="\t$sample.$replicate"
						header+="\t$condition.$replicate"
						meanheader+="\t$condition"

						# for global experiments file
						#[[ ${visited["$sample.$replicate"]} ]] && continue || visited["$sample.$replicate"]=1
						[[ ${visited["$condition.$replicate"]} ]] && continue || visited["$condition.$replicate"]=1
						# echo "$sample,$countfile,$condition,$replicate$factors" >> "$odir/experiments.csv"
						echo "$condition.$replicate,$countfile,$condition,$replicate$factors" >> "$odir/experiments.csv"
					done < <(awk -v c=$c '$2==c' "$f" | sort -k4,4V && awk -v t=$t '$2==t' "$f" | sort -k4,4V)

					commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD
						helper::multijoin
							-h "$(echo -e "$header")"
							-o "$odir/$c-vs-$t/experiments.tpm"
							-f $(printf '"%s" ' "${tojoin[@]}");
					CMD
						echo -e "$meanheader" > "$odir/$c-vs-$t/experiments.mean.tpm";
					CMD
						tail -n +2 "$odir/$c-vs-$t/experiments.tpm" >> "$odir/$c-vs-$t/experiments.mean.tpm";
					CMD
						Rscript - <<< '
							args <- commandArgs(TRUE);
							tsv <- args[1];
							df <- read.table(tsv, row.names=1, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
							means <- t(apply(df, 1, function(x) tapply(x, colnames(df), mean)));
							write.table(data.frame(id=rownames(means),means[,unique(colnames(df))],check.names=F), row.names = F, file = tsv, quote=F, sep="\t");
						'
					CMD
						"$odir/$c-vs-$t/experiments.mean.tpm"
					CMD

					commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
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
						"$odir/$c-vs-$t/experiments.tpm" "$odir/$c-vs-$t/experiments.tpm.zscores"
					CMD

					commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
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
						"$odir/$c-vs-$t/experiments.mean.tpm" "$odir/$c-vs-$t/experiments.mean.tpm.zscores"
					CMD

					# run deseq cmd3

					commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						[[ \$(wc -l < "$odir/$c-vs-$t/deseq.fcshrunk.tsv") -le 2 ]] && exit 0
					CMD
						head -1 "$odir/$c-vs-$t/experiments.tpm" > "$odir/$c-vs-$t/heatmap.tpm"
					CMD
						grep -Fw -f <(head -51 "$odir/$c-vs-$t/deseq.fcshrunk.tsv" | cut -f 1 | tail -n +2) "$odir/$c-vs-$t/experiments.tpm" >> "$odir/$c-vs-$t/heatmap.tpm"
					CMD
						heatmap.R TRUE 8 8 "$odir/$c-vs-$t/experiments.csv" "$odir/$c-vs-$t/heatmap.tpm" "TPM" "most differentially expressed ${feature}s"
					CMD

					commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						[[ \$(wc -l < "$odir/$c-vs-$t/deseq.fcshrunk.tsv") -le 2 ]] && exit 0
					CMD
						head -1 "$odir/$c-vs-$t/experiments.tpm.zscores" > "$odir/$c-vs-$t/heatmap.tpm.zscores"
					CMD
						grep -Fw -f <(head -51 "$odir/$c-vs-$t/deseq.fcshrunk.tsv" | cut -f 1 | tail -n +2) "$odir/$c-vs-$t/experiments.tpm.zscores" >> "$odir/$c-vs-$t/heatmap.tpm.zscores"
					CMD
						heatmap.R TRUE 8 8 "$odir/$c-vs-$t/experiments.csv" "$odir/$c-vs-$t/heatmap.tpm.zscores" "Z-Score" "most differentially expressed ${feature}s"
					CMD

					commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						[[ \$(wc -l < "$odir/$c-vs-$t/deseq.fcshrunk.tsv") -le 2 ]] && exit 0
					CMD
						head -1 "$odir/$c-vs-$t/experiments.mean.tpm" > "$odir/$c-vs-$t/heatmap.mean.tpm"
					CMD
						grep -Fw -f <(head -51 "$odir/$c-vs-$t/deseq.fcshrunk.tsv" | cut -f 1 | tail -n +2) "$odir/$c-vs-$t/experiments.mean.tpm" >> "$odir/$c-vs-$t/heatmap.mean.tpm"
					CMD
						heatmap.R FALSE 8 8 "$odir/$c-vs-$t/experiments.csv" "$odir/$c-vs-$t/heatmap.mean.tpm" "TPM" "most differentially expressed ${feature}s"
					CMD

					commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						[[ \$(wc -l < "$odir/$c-vs-$t/deseq.fcshrunk.tsv") -le 2 ]] && exit 0
					CMD
						head -1 "$odir/$c-vs-$t/experiments.mean.tpm.zscores" > "$odir/$c-vs-$t/heatmap.mean.tpm.zscores"
					CMD
						grep -Fw -f <(head -51 "$odir/$c-vs-$t/deseq.fcshrunk.tsv" | cut -f 1 | tail -n +2) "$odir/$c-vs-$t/experiments.mean.tpm.zscores" >> "$odir/$c-vs-$t/heatmap.mean.tpm.zscores"
					CMD
						heatmap.R FALSE 8 8 "$odir/$c-vs-$t/experiments.csv" "$odir/$c-vs-$t/heatmap.mean.tpm.zscores" "Z-Score" "most differentially expressed ${feature}s"
					CMD

					if [[ $gtf ]]; then
						commander::makecmd -a cmd5 -s ';' -c {COMMANDER[0]}<<- CMD
							find -L "$odir/$c-vs-$t"
								-type f
								-name "deseq*.tsv" -and -not -name "*annotated*"
 								-exec annotate.pl "${gtfinfo:=0}" "$gtf" $feature "{}" \;
						CMD

						for e in tpm tpm.zscores mean.tpm mean.tpm.zscores vsc vsc.zscores mean.vsc mean.vsc.zscores; do
							commander::makecmd -a cmd5 -s ';' -c {COMMANDER[0]}<<- CMD
								find -L "$odir/$c-vs-$t"
									-type f
									-name "deseq*.tsv" -and -not -name "*annotated*"
	 								-exec annotate.pl "${gtfinfo:=0}" "$gtf" $feature "{}" "$odir/$c-vs-$t/experiments.$e" $e \;
							CMD
						done

						commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD
							find -L "$odir/$c-vs-$t"
								-type f
								-name "deseq.full.fcshrunk.annotated.tsv"
 								-exec volcano.R "{}" "volcano_plot.fcshrunk.pdf" \;
						CMD

						commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD
							find -L "$odir/$c-vs-$t"
								-type f
								-name "deseq.full.annotated.tsv"
 								-exec volcano.R "{}" "volcano_plot.pdf" \;
						CMD

						for e in vsc tpm; do
							for h in "$odir/$c-vs-$t/heatmap.$e.localclust.ps" "$odir/$c-vs-$t/heatmap.$e.globalclust.ps" \
									"$odir/$c-vs-$t/heatmap.$e.zscores.localclust.ps" "$odir/$c-vs-$t/heatmap.$e.zscores.globalclust.ps" \
									"$odir/$c-vs-$t/heatmap.mean.$e.ps" "$odir/$c-vs-$t/heatmap.mean.$e.zscores.ps"; do
								commander::makecmd -a cmd5 -s ';' -c {COMMANDER[0]}<<- CMD
									if [[ -e "$h" ]]; then
										annotate.pl "${gtfinfo:=0}" "$gtf" $feature "$h";
										ps2pdf \$(grep -m 1 -F BoundingBox ${h%.*}.annotated.ps | awk '{print "-g"\$4*10"x"\$5*10}') ${h%.*}.annotated.ps ${h%.*}.annotated.pdf;
									fi
								CMD
							done
						done
					else
						commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD
							find -L "$odir/$c-vs-$t"
								-type f
								-name "deseq.full.fcshrunk.tsv"
 								-exec volcano.R "{}" "volcano_plot.fcshrunk.pdf" \;
						CMD

						commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD
							find -L "$odir/$c-vs-$t"
								-type f
								-name "deseq.full.tsv"
 								-exec volcano.R "{}" "volcano_plot.pdf" \;
						CMD
					fi
				done
			done
		done

		commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
			deseq2.R $ithreads "$odir/experiments.csv" "$odir" ${cmps[*]}
		CMD

		# deprecated use since annotation of vsc heatmap ps files is caputured in loop above along with tpm
		# expression::_deseq \
		# 	-1 cmd3 \
		# 	-2 cmdanno \
		# 	-t $ithreads \
		# 	-i "$odir/experiments.csv" \
		# 	-g "$gtf" \
		# 	-c "${cmps[*]}" \
		# 	-o "$odir"
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		# commander::printcmd -a cmdanno
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
		commander::printcmd -a cmd6
	else
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -v -b -i $instances -a cmd3
		# commander::runcmd -v -b -i $threads -a cmdanno
		commander::runcmd -v -b -i $threads -a cmd4
		commander::runcmd -v -b -i $threads -a cmd5
		commander::runcmd -v -b -i $threads -a cmd6
	fi

	return 0
}

function expression::_deseq(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-1 <cmds1>   | array of
			-2 <cmds2>   | array of
			-t <threads> | number of
			-i <incsv>   | path to
			-g <gtf>     | path to
			-c <cmps>    | string of pairs
			-o <outdir>  | path to
			-f <feature> | feature (default: gene)
		EOF
		return 1
	}

	local OPTIND arg mandatory threads csvfile gtf outdir gtfinfo feature="gene"
	declare -n _cmds1_deseq _cmds2_deseq
	declare -a cmppairs
	while getopts '1:2:t:i:g:c:o:f:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_deseq=$OPTARG;;
			2)	((++mandatory)); _cmds2_deseq=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			i)	((++mandatory)); csvfile="$OPTARG";;
			g)	gtf="$OPTARG"; gtfinfo="$([[ -s "$gtf.info" ]] && echo "$gtf.info" || true)";;
			c)	((++mandatory)); mapfile -t -d ' ' cmppairs < <(printf '%s' "$OPTARG");;
			o)	((++mandatory)); outdir="$OPTARG";;
			f)	feature="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage

	commander::makecmd -a _cmds1_deseq -s ';' -c {COMMANDER[0]}<<- CMD
		deseq2.R $threads "$csvfile" "$outdir" ${cmppairs[*]}
	CMD

	if [[ $gtf ]]; then
		local c t odir
		for i in $(seq 0 2 $((${#cmppairs[@]}-1))); do
			c=${cmppairs[$i]}
			t=${cmppairs[$((i+1))]}
			odir="$outdir/$c-vs-$t"

			commander::makecmd -a _cmds2_deseq -s ';' -c {COMMANDER[0]}<<- CMD
				find -L "$odir"
					-type f
					-name "deseq*.tsv" -and -not -name "*annotated*"
					-exec annotate.pl "${gtfinfo:=0}" "$gtf" $feature "{}" \;
			CMD

			for f in "$odir/heatmap.vsc.ps" "$odir/heatmap.vsc.zscores.ps" \
					"$odir/heatmap.vsc.globalclust.ps" "$odir/heatmap.vsc.zscores.globalclust.ps" "$odir/heatmap.vsc.localclust.ps" "$odir/heatmap.vsc.zscores.localclust.ps" \
					"$odir/heatmap.mean.vsc.ps" "$odir/heatmap.mean.vsc.zscores.ps"; do
				commander::makecmd -a _cmds2_deseq -s ';' -c {COMMANDER[0]}<<- CMD
					if [[ -e "$f" ]]; then
						annotate.pl "${gtfinfo:=0}" "$gtf" $feature "$f";
						ps2pdf \$(grep -m 1 -F BoundingBox ${f%.*}.annotated.ps | awk '{print "-g"\$4*10"x"\$5*10}') ${f%.*}.annotated.ps ${f%.*}.annotated.pdf;
					fi
				CMD
			done
		done
	fi

	return 0
}

function expression::join_deseq(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-g <gtf>      | path to
			-r <mapper>   | array of bams within array of
			-c <cmpfiles> | array of
			-i <countsdir>| path to
			-j <deseqdir> | path to
			-o <outdir>   | path to
			-f <feature>  | feature (default: gene)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads countsdir deseqdir outdir tmpdir="${TMPDIR:-/tmp}" gtf gtfinfo joinheatmaps=true feature="gene"
    declare -n _mapper_join _cmpfiles_join
	while getopts 'S:s:t:g:r:c:i:j:o:f:n' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	gtf="$OPTARG"; gtfinfo="$([[ -s "$gtf.info" ]] && echo "$gtf.info" || true)";;
			r)	((++mandatory)); _mapper_join=$OPTARG;;
			c)	((++mandatory)); _cmpfiles_join=$OPTARG;;
			i)	((++mandatory)); countsdir="$OPTARG";;
			j)	((++mandatory)); deseqdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			f)	feature="$OPTARG";;
			n)	joinheatmaps=false;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "joining tpm, vsc, zscores and heatmaps"

	declare -a cmd1 cmd2 cmd3 cmd4 mapdata tojoin
	declare -A visited
	local m f i c t e header meanheader cf vsc sample condition library replicate factors tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.join)"
	local topids="$tmp.topids" height width
	echo "rm -f '$tmp'*" >> "$BASHBONE_CLEANUP"

	for m in "${_mapper_join[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"
		visited=()
		header="id"
		meanheader="id"
		tojoin=()
		for f in "${_cmpfiles_join[@]}"; do
			mapfile -t mapdata < <(perl -F'\t' -lane 'next if exists $m{$F[1]}; $m{$F[1]}=1; print $F[1]' "$f")
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
					vsc="$deseqdir/$m/$c-vs-$t/experiments.vsc"
					head -51 "$deseqdir/$m/$c-vs-$t/deseq.tsv" | perl -lane 'print $F[0]."\t".abs($F[2])' | tail -n +2 >> "$topids"

					unset sample condition library replicate factors
					while read -r sample condition library replicate factors; do
						[[ ${visited["$sample"]} ]] && continue || visited["$sample"]=1
						header+="\t$condition.$replicate"
						meanheader+="\t$condition"

						cf="$(find -L "$countsdir/$m" -maxdepth 1 -name "$sample*.${feature}counts.htsc.tpm" -print -quit | grep .)"
						cf="${cf/%\.htsc\.tpm/}"
						tojoin+=("$cf")
						# get column per sample and write it to counts dir for downstream joins
						perl -M'List::MoreUtils qw(first_index)' -slane 'if($.==1){$c = first_index {$_ eq $h} @F}else{print "$F[0]\t$F[$c]"}' -- -h="$condition.$replicate" "$vsc" > "$cf.htsc.vsc"
					done < <(awk -v c=$c '$2==c' "$f" | sort -k4,4V && awk -v t=$t '$2==t' "$f" | sort -k4,4V)
				done
			done
		done

		sort -k2,2gr "$topids" | perl -lane 'next if exists $m{$F[0]}; $m{$F[0]}=1; print $F[0]' > "$tmp"
		mv "$tmp" "$topids"

		height=$(cat "$topids" | wc -l | awk '{h=($1/50)*8; if(h<2){h=2}; printf "%0.f",h}')
		width=$(cat "$deseqdir/$m/experiments.csv" | wc -l | awk '{printf "%0.f",$1/10+6}')

		for e in htsc htsc.tpm htsc.vsc; do
			f="$(printf '"%s" ' "${tojoin[@]/%/.$e}")"
			e=${e#*.}

			commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD
				helper::multijoin
					-h "$(echo -e "$header")"
					-o "$odir/experiments.$e"
					-f $f;
			CMD
				echo -e "$meanheader" > "$odir/experiments.mean.$e";
			CMD
				tail -n +2 "$odir/experiments.$e" >> "$odir/experiments.mean.$e";
			CMD
				Rscript - <<< '
					args <- commandArgs(TRUE);
					tsv <- args[1];
					df <- read.table(tsv, row.names=1, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
					means <- t(apply(df, 1, function(x) tapply(x, colnames(df), mean)));
					write.table(data.frame(id=rownames(means),means[,unique(colnames(df))],check.names=F), row.names = F, file = tsv, quote=F, sep="\t");
				'
			CMD
				"$odir/experiments.mean.$e"
			CMD

			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				Rscript - <<< '
					args <- commandArgs(TRUE);
					intsv <- args[1];
					outf <- args[2];
					dolog <- args[3];
					df <- read.table(intsv, row.names=1, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
					if(dolog) df <- log(df+1);
					df <- df-rowMeans(df);
					df <- df/apply(df,1,sd);
					df[is.na(df)] <- 0;
					write.table(data.frame(id=rownames(df),df,check.names=F), row.names = F, file = outf, quote=F, sep="\t");
				'
			CMD
				"$odir/experiments.$e" "$odir/experiments.$e.zscores" $([[ ! $e =~ vsc ]] && echo TRUE || echo false)
			CMD

			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				Rscript - <<< '
					args <- commandArgs(TRUE);
					intsv <- args[1];
					outf <- args[2];
					dolog <- args[3];
					df <- read.table(intsv, row.names=1, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
					if(dolog) df <- log(df+1);
					df <- df-rowMeans(df);
					df <- df/apply(df,1,sd);
					df[is.na(df)] <- 0;
					write.table(data.frame(id=rownames(df),df,check.names=F), row.names = F, file = outf, quote=F, sep="\t");
				'
			CMD
				"$odir/experiments.mean.$e" "$odir/experiments.mean.$e.zscores" $([[ ! $e =~ vsc ]] && echo TRUE || echo false)
			CMD

			if $joinheatmaps && [[ "$e" != "htsc" && $(wc -l < "$topids") -gt 1 ]]; then
				commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					head -1 "$odir/experiments.$e" > "$deseqdir/$m/heatmap.$e"
				CMD
					grep -Fw -f "$topids" "$odir/experiments.$e" >> "$deseqdir/$m/heatmap.$e"
				CMD
					heatmap.R TRUE $width $height "$deseqdir/$m/experiments.csv" "$deseqdir/$m/heatmap.$e" "${e^^}" "most differentially expressed ${feature}s"
				CMD

				commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					head -1 "$odir/experiments.$e.zscores" > "$deseqdir/$m/heatmap.$e.zscores"
				CMD
					grep -Fw -f "$topids" "$odir/experiments.$e.zscores" >> "$deseqdir/$m/heatmap.$e.zscores"
				CMD
					heatmap.R TRUE $width $height "$deseqdir/$m/experiments.csv" "$deseqdir/$m/heatmap.$e.zscores" "Z-Score" "most differentially expressed ${feature}s"
				CMD

				commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					head -1 "$odir/experiments.mean.$e" > "$deseqdir/$m/heatmap.mean.$e"
				CMD
					grep -Fw -f "$topids" "$odir/experiments.mean.$e" >> "$deseqdir/$m/heatmap.mean.$e"
				CMD
					heatmap.R FALSE $width $height "$deseqdir/$m/experiments.csv" "$deseqdir/$m/heatmap.mean.$e" "${e^^}" "most differentially expressed ${feature}s"
				CMD

				commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					head -1 "$odir/experiments.mean.$e.zscores" > "$deseqdir/$m/heatmap.mean.$e.zscores"
				CMD
					grep -Fw -f "$topids" "$odir/experiments.mean.$e.zscores" >> "$deseqdir/$m/heatmap.mean.$e.zscores"
				CMD
					heatmap.R FALSE $width $height "$deseqdir/$m/experiments.csv" "$deseqdir/$m/heatmap.mean.$e.zscores" "Z-Score" "most differentially expressed ${feature}s"
				CMD

				if [[ $gtf ]]; then
					for f in "$deseqdir/$m/heatmap.$e.localclust.ps" "$deseqdir/$m/heatmap.$e.globalclust.ps" \
							"$deseqdir/$m/heatmap.$e.zscores.localclust.ps" "$deseqdir/$m/heatmap.$e.zscores.globalclust.ps" \
							"$deseqdir/$m/heatmap.mean.$e.ps" "$deseqdir/$m/heatmap.mean.$e.zscores.ps"; do
						commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
							annotate.pl "${gtfinfo:=0}" "$gtf" $feature "$f"
						CMD
							ps2pdf \$(grep -m 1 -F BoundingBox ${f%.*}.annotated.ps | awk '{print "-g"\$4*10"x"\$5*10}') ${f%.*}.annotated.ps ${f%.*}.annotated.pdf
						CMD
					done
				fi
			fi
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}

function expression::join(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-c <cmpfiles> | array of
			-i <countsdir>| path to
			-o <outdir>   | path to
			-f <feature>  | feature (default: gene)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads countsdir outdir feature="gene"
    declare -n _mapper_join _cmpfiles_join
	while getopts 'S:s:t:r:c:i:o:f:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_join=$OPTARG;;
			c)	((++mandatory)); _cmpfiles_join=$OPTARG;;
			i)	((++mandatory)); countsdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			f)	feature="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "joining htsc, zscores"

	declare -a cmd1 cmd2 cmd3 cmd4 mapdata tojoin
	declare -A visited
	local m f i c t e header meanheader cf sample condition library replicate factors
	for m in "${_mapper_join[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"
		visited=()
		header="id"
		meanheader="id"
		tojoin=()
		for f in "${_cmpfiles_join[@]}"; do
			mapfile -t mapdata < <(perl -F'\t' -lane 'next if exists $m{$F[1]}; $m{$F[1]}=1; print $F[1]' "$f")
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do

					unset sample condition library replicate factors
					while read -r sample condition library replicate factors; do
						[[ ${visited["$sample"]} ]] && continue || visited["$sample"]=1
						header+="\t$condition.$replicate"
						meanheader+="\t$condition"

						cf="$(find -L "$countsdir/$m" -maxdepth 1 -name "$sample*.${feature}counts.htsc.tpm" -print -quit | grep .)"
						tojoin+=("${cf/%\.htsc\.tpm/}")
					done < <(awk -v c=$c '$2==c' "$f" | sort -k4,4V && awk -v t=$t '$2==t' "$f" | sort -k4,4V)
				done
			done
		done

		for e in htsc htsc.tpm; do
			f="$(printf '"%s" ' "${tojoin[@]/%/.$e}")"
			e=${e#*.}

			commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD
				helper::multijoin
					-h "$(echo -e "$header")"
					-o "$odir/experiments.$e"
					-f $f;
			CMD
				echo -e "$meanheader" > "$odir/experiments.mean.$e";
			CMD
				tail -n +2 "$odir/experiments.$e" >> "$odir/experiments.mean.$e";
			CMD
				Rscript - <<< '
					args <- commandArgs(TRUE);
					tsv <- args[1];
					df <- read.table(tsv, row.names=1, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
					means <- t(apply(df, 1, function(x) tapply(x, colnames(df), mean)));
					write.table(data.frame(id=rownames(means),means[,unique(colnames(df))],check.names=F), row.names = F, file = tsv, quote=F, sep="\t");
				'
			CMD
				"$odir/experiments.mean.$e"
			CMD

			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
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

			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
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
				"$odir/experiments.mean.$e" "$odir/experiments.mean.$e.zscores"
			CMD

		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}
