#! /usr/bin/env bash
# (c) Konstantin Riege

expression::diego() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-5 <skip>     | true/false md5sums, gtf prep respectively
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-g <gtf>      | path to
			-c <cmpfiles> | array of
			-i <htscdir>  | path to
			-j <mapeddir> | path to
			-p <tmpdir>   | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads countsdir mappeddir tmpdir outdir gtf tmp
	declare -n _mapper_diego _cmpfiles_diego
	while getopts 'S:s:5:t:r:g:c:i:j:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			5) $OPTARG && skipmd5=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			r) ((mandatory++)); _mapper_diego=$OPTARG;;
			g) ((mandatory++)); gtf="$OPTARG";;
			c) ((mandatory++)); _cmpfiles_diego=$OPTARG;;
			i) ((mandatory++)); countsdir="$OPTARG";;
			j) ((mandatory++)); mappeddir="$OPTARG";;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 8 ]] && _usage && return 1

	commander::print "differential splice junction analyses"

	$skipmd5 && {
		commander::warn "skip checking md5 sums and thus annotation preparation"
	} || {
		commander::print "checking md5 sums"
		local thismd5gtf
		declare -a cmdprep
		[[ -s "$gtf" ]] && thismd5gtf=$(md5sum "$gtf" | cut -d ' ' -f 1)
		if [[ ! -s ${gtf%.*}.diego.bed ]] || [[ "$thismd5gtf" && "$thismd5gtf" != "$md5gtf" ]]; then
			commander::print "preparing annotation for differential splice junction analyses"
			commander::makecmd -a cmdprep -s '&&' -c {COMMANDER[0]}<<- CMD
				gfftoDIEGObed.pl -g "$gtf" -o "${gtf%.*}.diego.bed"
			CMD
		fi

		if [[ ! -s "${gtf%.*}.aggregated.gtf" ]] || [[ "$thismd5gtf" && "$thismd5gtf" != "$md5gtf" ]]; then
			commander::print "preparing annotation for exon_id tag based quantification"
			tmp=$(mktemp -p "$tmpdir" --suffix=".gtf" cleanup.XXXXXXXXXX)
			commander::makecmd -a cmdprep -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				dexseq_prepare_annotation2.py
					-r no
					-f "$tmp"
					"$gtf" "${gtf%.*}.dexseq.gtf"
			CMD
				sed -r 's/(.+gene_id\s+")([^"]+)(.+exon_number\s+")([^"]+)(.+)/\1\2\3\4\5; exon_id "\2:\4"/' "$tmp" > "${gtf%.*}.aggregated.gtf"
			CMD
				rm -f "$tmpdir/tmp.gtf"
			CMD
		fi

		{	conda activate py3 && \
			commander::runcmd -v -b -t $threads -a cmdprep && \
			conda activate py2
		} || {
			rm -f $tmp
			commander::printerr "$funcname failed at annotation preparation"
			return 1
		}
		rm -rf $tmp
	}

	quantify::featurecounts \
		-S false \
		-s "$skip" \
		-t $threads \
		-p $tmpdir \
		-g "${gtf%.*}.aggregated.gtf" \
		-l exon \
		-f exon_id \
		-o "$countsdir" \
		-r _mapper_diego || {

		commander::printerr "$funcname failed"
		return 1
	}

	declare -a cmd1 cmd2 mapdata tdirs
	local m f i c t odir tdir sjfile countfile min sample condition library replicate factors
	for m in "${_mapper_diego[@]}"; do
		for f in "${_cmpfiles_diego[@]}"; do
			mapfile -t mapdata < <(cut -d $'\t' -f 2 $f | uniq)
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do 
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
					odir="$outdir/$m/$c-vs-$t"
					tdir="$tmpdir/$m"
					mkdir -p "$odir" "$tdir"
					rm -f "$odir/groups.tsv" "$odir/list.sj.tsv" "$odir/list.ex.tsv"
					unset sample condition library replicate factors
					while read -r sample condition library replicate factors; do
						sjfile=$(readlink -e "$mappeddir/$m/$sample"*.sj | head -1)
						countfile=$(readlink -e "$countsdir/$m/$sample"*.exoncounts.htsc | head -1)
						[[ $sjfile ]] && echo -e "$sample.$replicate\t$sjfile" >> "$odir/list.sj.tsv"
						echo -e "$sample.$replicate\t$countfile" >> "$odir/list.ex.tsv"
						echo -e "$condition\t$sample.$replicate" >> "$odir/groups.tsv"
					done < <(awk -v c=$c '$2==c' $f | sort -k4,4V && awk -v t=$t '$2==t' $f | sort -k4,4V)

					min=$(cut -d $'\t' -f 1 "$odir/groups.tsv" | sort | uniq -c | column -t | cut -d ' ' -f 1 | sort -k1,1 | head -1)
					if [[ -s "$odir/list.sj.tsv" ]]; then
						if [[ $m == "segemehl" ]]; then
							tdirs+=("$(mktemp -d -p "$tdir" cleanup.XXXXXXXXXX)")
							commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
								cd "${tdirs[-1]}"
							CMD
								pre_segemehl.pl
									-l "$odir/list.sj.tsv"
									-a "${gtf%.*}.diego.bed"
									-o "input.sj.tsv"
							CMD
								mv input.sj.tsv "$odir/input.sj.tsv"
							CMD
						elif [[ $m == "star" ]]; then
							tdirs+=("$(mktemp -d -p "$tdir" cleanup.XXXXXXXXXX)")
							commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
								cd "${tdirs[-1]}"
							CMD
								pre_STAR.py
									-l "$odir/list.sj.tsv"
									-d "${gtf%.*}.diego.bed"
							CMD
								mv junction_table.txt "$odir/input.sj.tsv"
							CMD
						fi
						commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
							diego.py
								-d $min 
								-a "$odir/input.sj.tsv"
								-b "$odir/groups.tsv"
								-x $c
								> "$odir/diego.sj.tsv"
						CMD
						commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
							diego.py
								-d $min 
								-a "$odir/input.sj.tsv"
								-b "$odir/groups.tsv"
								-x $c
								-e
								-f "$odir/dendrogram.sj"
						CMD
					fi

					tdirs+=("$(mktemp -d -p "$tdir" cleanup.XXXXXXXXXX)")
					commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
						cd "${tdirs[-1]}"
					CMD
						HTseq2DIEGO.pl
							-i "$odir/list.ex.tsv"
							-o "$odir/input.ex.tsv"
					CMD
					commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
						diego.py
							-d $min 
							-a "$odir/input.ex.tsv"
							-b "$odir/groups.tsv"
							-x $c
							> "$odir/diego.ex.tsv"
					CMD
					commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
						diego.py
							-d $min 
							-a "$odir/input.ex.tsv"
							-b "$odir/groups.tsv"
							-x $c
							-e
							-f "$odir/dendrogram.ex"
					CMD
				done
			done
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	conda activate py3 && \
			commander::runcmd -v -b -t $threads -a cmd1 && \
			commander::runcmd -v -b -t $threads -a cmd2 && \
			conda activate py2
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}

expression::deseq() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-g <gtf>      | path to
			-c <cmpfiles> | array of
			-i <htscdir>  | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads countsdir outdir gtf
	declare -n _mapper_deseq _cmpfiles_deseq
	while getopts 'S:s:t:r:g:c:i:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			r) ((mandatory++)); _mapper_deseq=$OPTARG;;
			g) gtf="$OPTARG";;
			c) ((mandatory++)); _cmpfiles_deseq=$OPTARG;;
			i) ((mandatory++)); countsdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	commander::print "principal component and differential gene expression analyses"

	local instances=${#_mapper_deseq[@]} ithreads m
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 64 -T $threads)

	declare -a cmd1 cmd2 cmps mapdata
	declare -A visited
	local f i c t odir countfile sample condition library replicate factors
	for m in "${_mapper_deseq[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"
		visited=()
		cmps=()

		head -1 ${_cmpfiles_deseq[0]} | perl -lane 'push @o,"sample,countfile,condition,replicate"; for (4..$#F){push @o,"factor".($_-3)}END{print join",",@o}' > "$odir/experiments.csv"
		for f in "${_cmpfiles_deseq[@]}"; do
			mapfile -t mapdata < <(cut -d $'\t' -f 2 $f | uniq)
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do 
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
					cmps+=("$c $t")
					unset sample condition library replicate factors
					while read -r sample condition library replicate factors; do
						[[ ${visited["$sample.$replicate"]} ]] && continue || visited["$sample.$replicate"]=1
						countfile=$(readlink -e "$countsdir/$m/$sample"*.+(genecounts|counts).+(reduced|htsc) | head -1)
						[[ $factors ]] && factors=","$(echo $factors | sed -r 's/\s+/,/g')
						echo "$sample.$replicate,$countfile,$condition,$replicate$factors" >> "$odir/experiments.csv"
					done < <(awk -v c=$c '$2==c' $f | sort -k4,4V && awk -v t=$t '$2==t' $f | sort -k4,4V)
				done
			done
		done

		expression::_deseq \
			-1 cmd1 \
			-2 cmd2 \
			-t $ithreads \
			-i "$odir/experiments.csv" \
			-g $gtf \
			-c "${cmps[*]}" \
			-o "$odir"
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	conda activate py2r && \
			commander::runcmd -v -b -t $instances -a cmd1 && \
			conda activate py2 && \
			commander::runcmd -v -b -t $threads -a cmd2
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

expression::_deseq() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-1 <cmds1>   | array of
			-2 <cmds2>   | array of
			-t <threads> | number of
			-i <incsv>   | path to
			-g <gtf>     | path to
			-c <cmps>    | string of pairs
			-o <outdir>  | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory threads csvfile gtf outdir gtfinfo
	declare -n _cmds1_deseq _cmds2_deseq
	declare -a cmppairs
	while getopts '1:2:t:i:g:c:o:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_deseq=$OPTARG;;
			2) ((mandatory++)); _cmds2_deseq=$OPTARG;;
			t) ((mandatory++)); threads=$OPTARG;;
			i) ((mandatory++)); csvfile="$OPTARG";;
			g) gtf="$OPTARG";;
			c) ((mandatory++)); mapfile -t -d ' ' cmppairs < <(printf '%s' "$OPTARG");;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage && return 1

	commander::makecmd -a _cmds1_deseq -s '|' -c {COMMANDER[0]}<<- CMD
		deseq2.R $threads "$csvfile" "$outdir" ${cmppairs[*]}
	CMD

	if [[ $gtf ]]; then
		gtfinfo=$(readlink -e "$gtf"*.+(info|descr) | head -1);
		local c t odir
		for i in $(seq 0 2 $((${#cmppairs[@]}-1))); do
			c=${cmppairs[$i]}
			t=${cmppairs[$((i+1))]}
			odir="$outdir/$c-vs-$t"
			for f in "$odir/deseq.tsv" "$odir/deseq.full.tsv" "$odir/deseq.noNA.tsv" \
					"$odir/heatmap.vsc.ps" "$odir/heatmap.mean.vsc.ps" "$odir/heatmap.vsc.zscores.ps" "$odir/heatmap.mean.vsc.zscores.ps"; do
				commander::makecmd -a _cmds2_deseq -s '|' -c {COMMANDER[0]}<<- CMD
					[[ -e "$f" ]] && annotate.pl "${gtfinfo:=0}" "$gtf" "$f" || true
				CMD
			done
		done
	fi

	return 0
}

expression::joincounts() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-p <tmpdir>   | path to
			-r <mapper>   | array of bams within array of
			-c <cmpfiles> | array of
			-i <htscdir>  | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads countsdir deseqdir outdir tmpdir
    declare -n _mapper_join _cmpfiles_join
	while getopts 'S:s:t:r:p:c:i:j:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			r) ((mandatory++)); _mapper_join=$OPTARG;;
			c) ((mandatory++)); _cmpfiles_join=$OPTARG;;
			i) ((mandatory++)); countsdir="$OPTARG";;
			j) ((mandatory++)); deseqdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1

	commander::print "joining count values plus zscore calculation"

	declare -a cmd1 cmd2 mapdata
	local m f x i c t h mh sample condition library replicate factors cf e tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX)"
	local tojoin="$tmp.tojoin" joined="$tmp.joined"
	for m in "${_mapper_join[@]}"; do
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"
		declare -A countfiles
		declare -a header meanheader
		x=0
		for f in "${_cmpfiles_join[@]}"; do
			mapfile -t mapdata < <(cut -d $'\t' -f 2 $f | uniq)
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do 
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do 
					vsc="$deseqdir/$m/$c-vs-$t/experiments.vsc"
					
					unset sample condition library replicate factors
					while read -r sample condition library replicate factors; do
						[[ ${countfiles["$sample.$replicate"]} ]] && continue
						header[$x]="$sample.$replicate"
						meanheader[$x]="$condition"
						((x++))

						cf=$(readlink -e "$countsdir/$m/$sample"*.tpm | head -1)
						cf="${cf%.*}"
						countfiles["$sample.$replicate"]="$cf"

						perl -M'List::MoreUtils qw(first_index)' -slane 'if($.==1){$c = first_index {$_ eq $h} @F}else{print "$F[0]\t$F[$c]"}' -- -h="$sample.$replicate" "$vsc" > $cf.vsc

					done < <(awk -v c=$c '$2==c' $f | sort -k4,4V && awk -v t=$t '$2==t' $f | sort -k4,4V)
				done
			done
		done

		for e in tpm vsc; do
			h='id'
			mh='id'
			rm -f $joined
			for x in "${!header[@]}"; do
				h+="\t${header[$x]}"
				mh+="\t${meanheader[$x]}"
				sort -k 1,1V "${countfiles[${header[$x]}]}.$e" > "$tojoin"
				if [[ -s "$joined" ]]; then
					join -t $'\t' "$joined" "$tojoin" > "$tmp"
					mv "$tmp" "$joined"
				else
					mv "$tojoin" "$joined"
				fi
			done

			echo -e "$h" > "$odir/experiments.$e"
			cat "$joined" >> "$odir/experiments.$e"
			echo -e "$mh" > "$tmp.$e"
			cat "$joined" >> "$tmp.$e"

			commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				Rscript - <<< '
					args <- commandArgs(TRUE);
					intsv <- args[1];
					outf <- args[2];
					df <- read.table(intsv, row.names=1, header=T, sep="\t", stringsAsFactors=F, check.names = F);
					means <- t(apply(df, 1, function(x) tapply(x, colnames(df), mean)));
					write.table(data.frame(id=rownames(means),means[,unique(colnames(df))]), row.names = F, file = outf, quote=F, sep="\t");
				'
			CMD
				"$tmp.$e" "$odir/experiments.mean.$e"
			CMD

			#df <- scale(log(df+1)) supposed to calculate z-scores but returns different values
			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				Rscript - <<< '
					args <- commandArgs(TRUE);
					intsv <- args[1];
					outf <- args[2];
					df <- read.table(intsv, row.names=1, header=T, sep="\t", stringsAsFactors=F, check.names = F);
					df <- log(df+1);
					df <- df-rowMeans(df);
					df <- df/apply(df,1,sd);
					df[is.na(df)] <- 0;
					write.table(data.frame(id=rownames(df),df), row.names = F, file = outf, quote=F, sep="\t");
				' 
			CMD
				"$odir/experiments.$e" "$odir/experiments.$e.zscores"
			CMD

			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				Rscript - <<< '
					args <- commandArgs(TRUE);
					intsv <- args[1];
					outf <- args[2];
					df <- read.table(intsv, row.names=1, header=T, sep="\t", stringsAsFactors=F, check.names = F);
					df <- log(df+1);
					df <- df-rowMeans(df);
					df <- df/apply(df,1,sd);
					df[is.na(df)] <- 0;
					write.table(data.frame(id=rownames(df),df), row.names = F, file = outf, quote=F, sep="\t");
				' 
			CMD
				"$odir/experiments.mean.$e" "$odir/experiments.mean.$e.zscores"
			CMD
		done
	done

	$skip && {
		commander::printcmd -a cmd1
	} || {
		{	conda activate py2r && \
			commander::runcmd -v -b -t $threads -a cmd1 && \
			commander::runcmd -v -b -t $threads -a cmd2 && \
			conda activate py2
		} || {
			rm -f "$tmp".* 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -f "$tmp".*
	return 0
}
