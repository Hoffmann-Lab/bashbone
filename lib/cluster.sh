#! /usr/bin/env bash
# (c) Konstantin Riege

cluster::coexpression(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-f <value>    | filter cluster for 0|1|2|20|21
			-t <threads>  | number of
			-m <memory>   | amount of
			-c <cmpfiles> | array of
			-z <idfiles>  | array of
			-r <mapper>   | array of bams within array of
			-p <tmpdir>   | path to
			-i <countsdir>| path to
			-j <deseqdir> | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory outdir tmpdir deseqdir countsdir clusterfilter=0
	declare -n _mapper_coexpression _cmpfiles_coexpression _idfiles_coexpression
	while getopts 'S:s:f:t:m:c:r:p:i:j:o:z:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			f) clusterfilter=$OPTARG;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			c) ((mandatory++)); _cmpfiles_coexpression=$OPTARG;;
			r) ((mandatory++)); _mapper_coexpression=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			i) ((mandatory++)); countsdir="$OPTARG";;
			j) ((mandatory++)); deseqdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			z) ((mandatory++)); _idfiles_coexpression=$OPTARG;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 9 ]] && _usage && return 1

	commander::print "inferring coexpression"

	declare -a mapdata cmd1
	declare -A visited
	local m f i j c t e tdir odir cdir ddir params tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX)"
	local tojoin="$tmp.tojoin" joined="$tmp.joined"
	for m in "${_mapper_coexpression[@]}"; do
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		cdir="$countsdir/$m"
		ddir="$deseqdir/$m"
		mkdir -p "$odir" "$tdir"
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
				$okmean=1 if $F[$i] > 0;
				$okfc=1 if exists $F[$i+4] && abs(abs($F[$i+1])-abs($F[$i+4]))>=0.5;
				$okpval=1 if $F[$i+2] > 0 && $F[$i+2] <= 0.05;
			}
			print $F[0] if ($cf =~ /2$/ ) || ($cf =~ /0$/ && $okpval == 1) || ($cf =~ /1$/ && $okmean+$okfc+$okpval == 3);
		' -- -cf=$clusterfilter "$odir/experiments.deseq.tsv" > "$odir/experiments.filtered.genes"

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
			commander::printerr "$funcname failed"
			return 1
		}
	}

	declare -a cmd2
	local type wdir
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
