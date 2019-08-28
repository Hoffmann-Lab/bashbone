#! /usr/bin/env bash
# (c) Konstantin Riege

expression::deseq() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip> | true/false return
			-s <softskip> | truefalse only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-g <gtf>      | path to
			-c <cmpfiles> | array of
			-i <htscdir>  | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads countsdir outdir gtf gtfinfo
	declare -n _mapper_deseq _cmpfiles_deseq
	while getopts 'S:s:t:r:g:c:i:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			r) ((mandatory++)); _mapper_deseq=$OPTARG;;
			g) gtf="$OPTARG"; [[ -s "$gtf.info" ]] && "gtfinfo=$gtf.info"; [[ -s "$gtf.descr" ]] && gtfinfo="$gtf.descr";;
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
	local f i c t o odir countfile sample condition library replicate pair
	for m in "${_mapper_deseq[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"
		visited=()
        cmps=()

		if [[ $(awk '{print $5}' ${_cmpfiles_deseq[0]}) ]]; then
			echo 'sample,countfile,condition,replicate,pairs' > "$odir/experiments.csv"
		else
			echo 'sample,countfile,condition,replicate' > "$odir/experiments.csv"
		fi

		for f in "${_cmpfiles_deseq[@]}"; do
			mapfile -t mapdata < <(cut -d $'\t' -f 2 $f | uniq)
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do 
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
					cmps+=("$c $t")
					unset sample condition library replicate pair
					while read -r sample condition library replicate pair; do
						[[ ${visited["$sample.$replicate"]} ]] && continue || visited["$sample.$replicate"]=1
						countfile=$(readlink -e "$countsdir/$m/$sample"*.+(reduced|htsc) | head -1)
						[[ $pair ]] && pair=",$pair"
						echo "$sample.$replicate,$countfile,$condition,$replicate$pair" >> "$odir/experiments.csv"
					done < <(awk -v c=$c '$2==c' $f | sort -k4,4V && awk -v t=$t '$2==t' $f | sort -k4,4V)

					if [[ $gtf ]]; then
						for o in "$odir/$m/$c-vs-$t/deseq.tsv" "$odir/$m/$c-vs-$t/deseq.full.tsv" "$odir/$m/$c-vs-$t/deseq.noNA.tsv" "$odir/$m/$c-vs-$t/heatmap.vsc.ps" "$odir/$m/$c-vs-$t/heatmap.mean.vsc.ps"; do
							commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
								annotate.pl "${gtfinfo:=0}" "$gtf" "$o"
							CMD
						done
					fi
				done
			done
		done

		commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
			deseq2.R $ithreads "$odir/experiments.csv" "$odir" ${cmps[*]}
		CMD
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
	local m f x i c t h mh sample condition library replicate pair cf e tmp="$(mktemp -p "$tmpdir")"
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
					
					unset sample condition library replicate pair
					while read -r sample condition library replicate pair; do
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
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
