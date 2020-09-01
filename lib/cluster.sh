#! /usr/bin/env bash
# (c) Konstantin Riege

cluster::coexpression_deseq(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-f <value>    | filter cluster for 0|1|2|20|21
			-t <threads>  | number of
			-m <memory>   | amount of
			-c <cmpfiles> | array of
			-l <idfiles>  | array of
			-r <mapper>   | array of bams within array of
			-p <tmpdir>   | path to
			-i <countsdir>| path to
			-j <deseqdir> | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory outdir tmpdir deseqdir countsdir clusterfilter=0 biotype gtf
	declare -n _mapper_coexpression _cmpfiles_coexpression _idfiles_coexpression
	while getopts 'S:s:f:b:g:t:m:c:r:p:i:j:o:l:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			f) clusterfilter=$OPTARG;;
			b) biotype=$OPTARG;;
			g) gtf="$OPTARG";;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			c) ((mandatory++)); _cmpfiles_coexpression=$OPTARG;;
			r) ((mandatory++)); _mapper_coexpression=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			i) ((mandatory++)); countsdir="$OPTARG";;
			j) ((mandatory++)); deseqdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			l) ((mandatory++)); _idfiles_coexpression=$OPTARG;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 9 ]] && _usage && return 1
	[[ $biotype && ! $gtf ]] && _usage && return 1

	commander::printinfo "inferring coexpression"

	declare -a mapdata cmd1
	declare -A visited
	local m f i c t e odir cdir ddir params tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.join)"
	local tojoin="$tmp.tojoin" joined="$tmp.joined"
	for m in "${_mapper_coexpression[@]}"; do
		odir="$outdir/$m"
		cdir="$countsdir/$m"
		ddir="$deseqdir/$m"
		mkdir -p "$odir"
		visited=()

		rm -f $joined
		for f in "${_cmpfiles_coexpression[@]}"; do
			mapfile -t mapdata < <(cut -d $'\t' -f 2 $f | uniq)
			i=0
			for c in "${mapdata[@]::${#mapdata[@]}-1}"; do
				for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
					[[ ${visited["$c-vs-$t"]} ]] && continue || visited["$c-vs-$t"]=1

					awk 'NR>1' "$ddir/$c-vs-$t/deseq.full.tsv" | sort -k1,1V | cut -f 1,2,3,7 | sed 's/NA/0/g' > "$tojoin"
					if [[ -s "$joined" ]]; then
						join -t $'\t' "$joined" "$tojoin" > "$tmp"
						mv "$tmp" "$joined"
					else
						mv "$tojoin" "$joined"
					fi
				done
			done
		done
		mv $joined $odir/experiments.deseq.tsv

		# filter joined deseq tables requiers padj >0 due to NA replacement
		perl -F'\t' -slane '
			$okmean=0;
			$okfc=0;
			$okpval=0;
			for (my $i=1; $i<$#F; $i+=3){
				$okmean=1 if $F[$i] > 5;
				$okfc=1 if exists $F[$i+4] && abs(abs($F[$i+1])-abs($F[$i+4]))>=0.5;
				$okpval=1 if $F[$i+2] > 0 && $F[$i+2] <= 0.05;
			}
			if ($cf=~/[012]{3,3}/){
				print $F[0] if $okmean+$okfc+$okpval == 3;
			} elsif ($cf=~/[01]{2,2}/) {
				print $F[0] if $okfc+$okpval == 2;
			} elsif ($sc==0) {
				print $F[0] if $okpval == 1;
			} elsif	($sc==1) {
				print $F[0] if $okfc == 1;
			} else {
				print $F[0];
			}
		' -- -cf=$clusterfilter "$odir/experiments.deseq.tsv" > "$odir/experiments.filtered.genes"
		# 0 padj 1 fc>0.5 2 basemean > 5 3 30% see below

		[[ $biotype ]] && {
			perl -lane -F'\t' '{
				if($#F==3){
					print $F[0] if $F[2]==$cb;
				} else {
					next unless $F[2]=="gene";
					$F[-1]=~/gene_biotype\s+"([^"]+)/;
					if ($1 eq $cb && $F[-1]=~/gene_id\s+"([^"]+)/){
						print $1;
					}
				}
			}' -- -cb=$biotype "$(readlink -e "$gtf"*.+(info|descr) "$gtf" | head -1)" > "$tmp.genes"
			grep -f "$tmp.genes" "$odir/experiments.filtered.genes" > "$tmp.filtered.genes"
			mv "$tmp.filtered.genes" "$odir/experiments.filtered.genes"
		}

		for e in tpm vsc; do
			[[ ! -s "$cdir/experiments.$e" ]] && continue
			head -1 "$cdir/experiments.$e" > "$odir/experiments.filtered.$e"
			grep -f "$odir/experiments.filtered.genes" "$cdir/experiments.$e" >> "$odir/experiments.filtered.$e"

			mkdir -p "$odir/$e"
			params=FALSE
			[[ $clusterfilter =~ 2 ]] && params=TRUE
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				wgcna.R $((memory/1024/2)) ${e^^*} $params "$odir/experiments.filtered.$e" "$odir/$e"
			CMD
		done
	done

	$skip && {
		commander::printcmd -a cmd1
	} || {
		{	conda activate py2r && \
			commander::runcmd -v -b -t $threads -a cmd1 && \
			conda activate py2
		} || {
			rm -f "$tmp"*
			commander::printerr "$funcname failed"
			return 1
		}
	}
	rm -f "$tmp"*

	declare -a cmd2
	local m e i type wdir cdir wdir odir
	for m in "${_mapper_coexpression[@]}"; do
		cdir="$countsdir/$m"
		for e in tpm vsc; do
			[[ ! -s "$cdir/experiments.$e" ]] && continue
			wdir="$outdir/$m/$e"

			# work on modules
			echo "$(head -1 $cdir/experiments.$e) type" | tr ' ' '\t' > "$wdir/modules.$e"
			echo "$(head -1 $cdir/experiments.mean.$e) type" | tr ' ' '\t' > "$wdir/modules.mean.$e"
			cp "$wdir/modules.$e" "$wdir/modules.$e.zscores"
			cp "$wdir/modules.mean.$e" "$wdir/modules.mean.$e.zscores"

			for i in $(seq 0 $(cut -d ' ' -f 2- "$wdir/wgcna.cluster2modules" | sed 's/ /\n/g' | sort -n | tail -1)); do
				type="Module.$(printf '%03d' $i)"
				odir="$wdir/$type"
				mkdir -p "$odir"
				awk -v i=$i '$NF==i{print $1}' "$wdir/wgcna.modules.tsv" > "$odir/genes.list"

				# add to array for later go enrichment
				_idfiles_coexpression+=("$odir/genes.list")

				echo "$(head -1 $cdir/experiments.$e) type" | tr ' ' '\t' > "$odir/experiments.$e"
				echo "$(head -1 $cdir/experiments.mean.$e) type" | tr ' ' '\t' > "$odir/experiments.mean.$e"
				cp "$odir/experiments.$e" "$odir/experiments.$e.zscores"
				cp "$odir/experiments.mean.$e" "$odir/experiments.mean.$e.zscores"
				grep -f "$odir/genes.list" "$cdir/experiments.$e" | sed "s/$/\t$type/" | tee -a "$wdir/modules.$e" >> "$odir/experiments.$e"
				grep -f "$odir/genes.list" "$cdir/experiments.mean.$e" | sed "s/$/\t$type/" | tee -a "$wdir/modules.mean.$e" >> "$odir/experiments.mean.$e"
				grep -f "$odir/genes.list" "$cdir/experiments.$e.zscores" | sed "s/$/\t$type/" | tee -a "$wdir/modules.$e.zscores" >> "$odir/experiments.$e.zscores"
				grep -f "$odir/genes.list" "$cdir/experiments.mean.$e.zscores" | sed "s/$/\t$type/" | tee -a "$wdir/modules.mean.$e.zscores" >> "$odir/experiments.mean.$e.zscores"

				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Module ${e^^*} "$odir/experiments.$e" "$odir/experiments.$e"
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Module ${e^^*} "$odir/experiments.mean.$e" "$odir/experiments.mean.$e"
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Module Z-Score "$odir/experiments.$e.zscores" "$odir/experiments.$e.zscores"
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Module Z-Score "$odir/experiments.mean.$e.zscores" "$odir/experiments.mean.$e.zscores"
				CMD
			done

			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Modules ${e^^*} "$wdir/modules.$e" "$wdir/modules.$e"
			CMD
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Modules ${e^^*} "$wdir/modules.mean.$e" "$wdir/modules.mean.$e"
			CMD
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Modules Z-Score "$wdir/modules.$e.zscores" "$wdir/modules.$e.zscores"
			CMD
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Modules Z-Score "$wdir/modules.mean.$e.zscores" "$wdir/modules.mean.$e.zscores"
			CMD

			# work on cluster
			echo "$(head -1 $cdir/experiments.$e) type" | tr ' ' '\t' > "$wdir/cluster.$e"
			echo "$(head -1 $cdir/experiments.mean.$e) type" | tr ' ' '\t' > "$wdir/cluster.mean.$e"
			cp "$wdir/cluster.$e" "$wdir/cluster.$e.zscores"
			cp "$wdir/cluster.mean.$e" "$wdir/cluster.mean.$e.zscores"
			for i in $(seq 0 $(cut -d ' ' -f 1 "$wdir/wgcna.cluster2modules" | sort -n | tail -1)); do
				type="Cluster.$(printf '%03d' $i)"
				odir="$wdir/$type"
				mkdir -p "$odir"
				awk -v i=$i '$NF==i{print $1}' "$wdir/wgcna.cluster.tsv" > "$odir/genes.list"

				# add to array for later go enrichment
				_idfiles_coexpression+=("$odir/genes.list")

				echo "$(head -1 $cdir/experiments.$e) type" | tr ' ' '\t' > "$odir/experiments.$e"
				echo "$(head -1 $cdir/experiments.mean.$e) type" | tr ' ' '\t' > "$odir/experiments.mean.$e"
				cp "$odir/experiments.$e" "$odir/experiments.$e.zscores"
				cp "$odir/experiments.mean.$e" "$odir/experiments.mean.$e.zscores"

				grep -f "$odir/genes.list" "$cdir/experiments.$e" | sed "s/$/\t$type/" >> "$wdir/cluster.$e"
				grep -f "$odir/genes.list" "$cdir/experiments.mean.$e" | sed "s/$/\t$type/" >> "$wdir/cluster.mean.$e"
				grep -f "$odir/genes.list" "$cdir/experiments.$e.zscores" | sed "s/$/\t$type/" >> "$wdir/cluster.$e.zscores"
				grep -f "$odir/genes.list" "$cdir/experiments.mean.$e.zscores" | sed "s/$/\t$type/" >> "$wdir/cluster.mean.$e.zscores"

				for j in $(awk -v i=$i '$1==i{$1=""; print}' "$wdir/wgcna.cluster2modules"); do
					type="Module.$(printf '%03d' $j)"
					awk 'NR>1' "$wdir/$type/experiments.$e" >> "$odir/experiments.$e"
					awk 'NR>1' "$wdir/$type/experiments.mean.$e" >> "$odir/experiments.mean.$e"
					awk 'NR>1' "$wdir/$type/experiments.$e.zscores" >> "$odir/experiments.$e.zscores"
					awk 'NR>1' "$wdir/$type/experiments.mean.$e.zscores" >> "$odir/experiments.mean.$e.zscores"
				done

				# Modules as legend title is correct here
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Modules ${e^^*} "$odir/experiments.$e" "$odir/experiments.$e"
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Modules ${e^^*} "$odir/experiments.mean.$e" "$odir/experiments.mean.$e"
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Modules Z-Score "$odir/experiments.$e.zscores" "$odir/experiments.$e.zscores"
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Modules Z-Score "$odir/experiments.mean.$e.zscores" "$odir/experiments.mean.$e.zscores"
				CMD
			done

			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster ${e^^*} "$wdir/cluster.$e" "$wdir/cluster.$e"
			CMD
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster ${e^^*} "$wdir/cluster.mean.$e" "$wdir/cluster.mean.$e"
			CMD
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster Z-Score "$wdir/cluster.$e.zscores" "$wdir/cluster.$e.zscores"
			CMD
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster Z-Score "$wdir/cluster.mean.$e.zscores" "$wdir/cluster.mean.$e.zscores"
			CMD

		done
	done

	$skip && {
		commander::printcmd -a cmd2
	} || {
		{	conda activate py2r && \
			commander::runcmd -v -b -t $threads -a cmd2 && \
			conda activate py2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

cluster::coexpression(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-f <value>    | filter cluster for 0|1|2|20|21
			-t <threads>  | number of
			-m <memory>   | amount of
			-l <idfiles>  | array of
			-r <mapper>   | array of bams within array of
			-p <tmpdir>   | path to
			-i <countsdir>| path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory outdir tmpdir countsdir biotype gtf clusterfilter=0
	declare -n _mapper_coexpression _idfiles_coexpression
	while getopts 'S:s:f:b:g:t:m:c:r:p:i:j:o:l:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			f) clusterfilter=$OPTARG;;
			b) biotype=$OPTARG;;
			g) gtf="$OPTARG";;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			r) ((mandatory++)); _mapper_coexpression=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			i) ((mandatory++)); countsdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			l) ((mandatory++)); _idfiles_coexpression=$OPTARG;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1
	[[ $biotype && ! $gtf ]] && _usage && return 1

	commander::printinfo "inferring coexpression"

	declare -a cmd1
	local m f cdir odir suff header sample countfile tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.join)"
	local tojoin="$tmp.tojoin" joined="$tmp.joined"

	for m in "${_mapper_coexpression[@]}"; do
		cdir="$countsdir/$m"
		odir="$outdir/$m"
		mkdir -p "$odir"

		declare -n _bams_coexpression=$m
		suff=$(cat <(echo -e "${_bams_coexpression[0]}\n${_bams_coexpression[1]}") | rev | paste - - | sed -E 's/(.+\.).+\t\1.+/\1/' | rev)
		header='id'
		rm -f $joined

		for f in "${_bams_coexpression[@]}"; do
			sample=$(basename $f $suff)
			countfile=$(readlink -e "$countsdir/$m/$sample"*.tpm | head -1)
			header+="\t$sample"

			sort -k 1,1V "$countfile" > "$tojoin"
			if [[ -s "$joined" ]]; then
				join -t $'\t' "$joined" "$tojoin" > "$tmp"
				mv "$tmp" "$joined"
			else
				mv "$tojoin" "$joined"
			fi
		done

		echo -e "$header" > "$odir/experiments.tpm"
		cat "$joined" >> "$odir/experiments.tpm"

		commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
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
			"$odir/experiments.tpm" "$odir/experiments.tpm.zscores"
		CMD
	done

	$skip && {
		commander::printcmd -a cmd1
	} || {
		{	conda activate py2r && \
			commander::runcmd -v -b -t $threads -a cmd1 && \
			conda activate py2
		} || {
			rm -f "$tmp"*
			commander::printerr "$funcname failed"
			return 1
		}
	}
	rm -f "$tmp"*

	declare -a cmd2
	local m odir params
	for m in "${_mapper_coexpression[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"

		[[ $clusterfilter =~ 2 ]] && {
			awk '{i=2; ok=0; while(i<NF){if($i>=5){ok=1} i=i+1} if(ok){print $1}}' "$odir/experiments.tpm" > "$odir/experiments.filtered.genes"
		} || {
			awk '{print $1}' "$odir/experiments.tpm" > "$odir/experiments.filtered.genes"
		}

		[[ $biotype ]] && {
			perl -F'\t' -slane '{
				if($#F==3){
					print $F[0] if $F[2] eq $cb;
				} else {
					next unless $F[2] eq "gene";
					$F[-1]=~/gene_biotype\s+"([^"]+)/;
					if ($1 eq $cb && $F[-1]=~/gene_id\s+"([^"]+)/){
						print $1;
					}
				}
			}' -- -cb=$biotype "$(readlink -e "$gtf"*.+(info|descr) "$gtf" | head -1)" > "$tmp.genes"
			grep -f "$tmp.genes" "$odir/experiments.filtered.genes" > "$tmp.filtered.genes"
			mv "$tmp.filtered.genes" "$odir/experiments.filtered.genes"
		}

		head -1 "$odir/experiments.tpm" > "$odir/experiments.filtered.tpm"
		grep -f "$odir/experiments.filtered.genes" "$odir/experiments.tpm" >> "$odir/experiments.filtered.tpm"

		params=FALSE
		[[ $clusterfilter =~ 2 ]] && params=TRUE
		commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
			wgcna.R $((memory/1024/2)) TPM $params "$odir/experiments.filtered.tpm" "$odir"
		CMD
	done

	$skip && {
		commander::printcmd -a cmd2
	} || {
		{	conda activate py2r && \
			commander::runcmd -v -b -t $threads -a cmd2 && \
			conda activate py2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	declare -a cmd3
	local m e i j wdir odir type
	for m in "${_mapper_coexpression[@]}"; do
		wdir="$outdir/$m"
		for e in tpm; do
			[[ ! -s "$wdir/experiments.$e" ]] && continue

			# work on modules
			echo "$(head -1 $wdir/experiments.$e) type" | tr ' ' '\t' > "$wdir/modules.$e"
			cp "$wdir/modules.$e" "$wdir/modules.$e.zscores"

			for i in $(seq 0 $(cut -d ' ' -f 2- "$wdir/wgcna.cluster2modules" | sed 's/ /\n/g' | sort -n | tail -1)); do
				type="Module.$(printf '%03d' $i)"
				odir="$wdir/$type"
				mkdir -p "$odir"
				awk -v i=$i '$NF==i{print $1}' "$wdir/wgcna.modules.tsv" > "$odir/genes.list"

				# add to array for later go enrichment
				_idfiles_coexpression+=("$odir/genes.list")

				echo "$(head -1 $wdir/experiments.$e) type" | tr ' ' '\t' > "$odir/experiments.$e"
				cp "$odir/experiments.$e" "$odir/experiments.$e.zscores"
				grep -f "$odir/genes.list" "$wdir/experiments.$e" | sed "s/$/\t$type/" | tee -a "$wdir/modules.$e" >> "$odir/experiments.$e"
				grep -f "$odir/genes.list" "$wdir/experiments.$e.zscores" | sed "s/$/\t$type/" | tee -a "$wdir/modules.$e.zscores" >> "$odir/experiments.$e.zscores"

				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Module ${e^^*} "$odir/experiments.$e" "$odir/experiments.$e"
				CMD
				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Module Z-Score "$odir/experiments.$e.zscores" "$odir/experiments.$e.zscores"
				CMD
			done

			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Modules ${e^^*} "$wdir/modules.$e" "$wdir/modules.$e"
			CMD
			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Modules Z-Score "$wdir/modules.$e.zscores" "$wdir/modules.$e.zscores"
			CMD

			# work on cluster
			echo "$(head -1 $wdir/experiments.$e) type" | tr ' ' '\t' > "$wdir/cluster.$e"
			cp "$wdir/cluster.$e" "$wdir/cluster.$e.zscores"

			for i in $(seq 0 $(cut -d ' ' -f 1 "$wdir/wgcna.cluster2modules" | sort -n | tail -1)); do
				type="Cluster.$(printf '%03d' $i)"
				odir="$wdir/$type"
				mkdir -p "$odir"
				awk -v i=$i '$NF==i{print $1}' "$wdir/wgcna.cluster.tsv" > "$odir/genes.list"

				# add to array for later go enrichment
				_idfiles_coexpression+=("$odir/genes.list")

				echo "$(head -1 $wdir/experiments.$e) type" | tr ' ' '\t' > "$odir/experiments.$e"
				cp "$odir/experiments.$e" "$odir/experiments.$e.zscores"

				grep -f "$odir/genes.list" "$wdir/experiments.$e" | sed "s/$/\t$type/" >> "$wdir/cluster.$e"
				grep -f "$odir/genes.list" "$wdir/experiments.$e.zscores" | sed "s/$/\t$type/" >> "$wdir/cluster.$e.zscores"

				for j in $(awk -v i=$i '$1==i{$1=""; print}' "$wdir/wgcna.cluster2modules"); do
					local type="Module.$(printf '%03d' $j)"
					awk 'NR>1' "$wdir/$type/experiments.$e" >> "$odir/experiments.$e"
					awk 'NR>1' "$wdir/$type/experiments.$e.zscores" >> "$odir/experiments.$e.zscores"
				done

				# Modules as legend title is correct here
				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Modules ${e^^*} "$odir/experiments.$e" "$odir/experiments.$e"
				CMD
				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
					vizco.R Modules Z-Score "$odir/experiments.$e.zscores" "$odir/experiments.$e.zscores"
				CMD
			done

			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster ${e^^*} "$wdir/cluster.$e" "$wdir/cluster.$e"
			CMD
			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
				vizco.R Cluster Z-Score "$wdir/cluster.$e.zscores" "$wdir/cluster.$e.zscores"
			CMD

		done
	done

	$skip && {
		commander::printcmd -a cmd3
	} || {
		{	conda activate py2r && \
			commander::runcmd -v -b -t $threads -a cmd3 && \
			conda activate py2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
