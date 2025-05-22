#! /usr/bin/env bash
# (c) Konstantin Riege

function peaks::_bed2narrowpeak(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-1 <cmds1>    | array of. run with conda env macs
			-2 <cmds2>    | array of
			-i <file>     | treatment bed file
			-f <file>     | treatment bam file
			-b <file>     | bed file
			-o <outfile>  | narrowPeak file
			-p <tmpdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory n i f b o tmpdir
	declare -n _cmds1_bed2narrowpeak _cmds2_bed2narrowpeak
	while getopts '1:2:i:f:b:o:p:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_bed2narrowpeak=$OPTARG;;
			2)	((++mandatory)); _cmds2_bed2narrowpeak=$OPTARG;;
			i)	((++mandatory)); i="$OPTARG";;
			f)	((++mandatory)); f="$OPTARG";;
			b)	((++mandatory)); b="$OPTARG";;
			o)	((++mandatory)); o="$OPTARG";;
			p)	((++mandatory)); tmpdir="$OPTARG";; # needs to be given here, otherwise a mktemp dir will be deleted upon returning of this function
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 7 ]] && _usage

	# unless summit: refinepeak
	commander::makecmd -a _cmds1_bed2narrowpeak -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
		perl -lane '@o=("",0,0,"peak_$.",0,".",-1,-1,-1,-1); $o[$_]=$F[$_] for 0..$#F; print join"\t",@o'
	CMD
		"$b" | sort -k1,1 -k2,2n -k3,3n > "$tmpdir/peaks";
	CMD
		if [[ \$(head -1 "$tmpdir/peaks" | cut -f 10) == -1 ]]; then
			paste "$tmpdir/peaks" <(macs2 refinepeak
				-b "$tmpdir/peaks"
				-i "$i"
				-f BED
				-w 200
				-o /dev/stdout | sort -k1,1 -k2,2n -k3,3n)
			| awk -v OFS='\t' '{\$10=\$12-\$2; print}' | cut -f 1-10 > "$tmpdir/refinedpeaks";
		else
			ln -sfnr "$tmpdir/peaks" "$tmpdir/refinedpeaks";
		fi
	CMD

	# correct refinepeak results when not in region; afterwards get fold-enrichment unless already present
	# if multiple summits with identical hight, report first summit and not (summit1-summitN)/2
	# $F[9]=sprintf("%0.d",1+$x[-1]->[0]+($x[-1]->[-1]-$x[-1]->[0])/2);
	commander::makecmd -a _cmds2_bed2narrowpeak -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- 'CMD' {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- CMD {COMMANDER[7]}<<- CMD {COMMANDER[8]}<<- 'CMD' {COMMANDER[9]}<<- CMD {COMMANDER[10]}<<- CMD {COMMANDER[11]}<<- 'CMD' {COMMANDER[12]}<<- CMD
		samtools depth -aa -@ 1 -g 1796 -b <(awk '\$10<0 || \$2+\$10>\$3' "$tmpdir/refinedpeaks" | cut -f 1-3) "$f" |
	CMD
		awk -v OFS='\t' '{print $1,$2-1,$2,$3}' |
	CMD
		bedtools merge -i - -c 4 -o collapse 2> /dev/null |
	CMD
		sort -k1,1 -k2,2n -k3,3n |
	CMD
		perl -slane '
			BEGIN{open B,"<$b" or die $!}
			@c=split/,/,$F[3];
			@x=();
			push $x[$c[$_]]->@*,$_ for 0..$#c;
			$l=<B>; chomp($l); @F=split/\t/,$l;
			$F[9]=1+$x[-1]->[0];
			print join"\t",@F[0..9]
		'
	CMD
		-- -b=<(awk '\$10<0 || \$2+\$10>\$3' "$tmpdir/refinedpeaks") |
	CMD
		cat - <(awk '\$10<0 || \$2+\$10>\$3{next}{print}' "$tmpdir/refinedpeaks") |
	CMD
		sort -k1,1 -k2,2n -k3,3n |
	CMD
		awk -v OFS='\t' '{print $1,$2+$10,$2+$10+1,$0}'
	CMD
		> "$tmpdir/summits";
	CMD
		if [[ \$(head -1 "$tmpdir/summits" | cut -f 10) == -1 ]]; then
			paste
				<(bedtools intersect -loj -a "$tmpdir/summits" -b "$(dirname "$i")/model_fc.bedg" | awk -v OFS='\t' '{\$10=\$NF; print}' | cut -f 4-13)
				<(bedtools intersect -loj -a "$tmpdir/summits" -b "$(dirname "$i")/nomodel_fc.bedg" | awk -v OFS='\t' '{\$10=\$NF; print}' | cut -f 4-13);
		else
			paste <(cut -f 4-13 "$tmpdir/summits") <(cut -f 4-13 "$tmpdir/summits");
		fi |
	CMD
	 	awk -v OFS='\t' '{if($17>$7){$7=$17} print}' |
	CMD
		cut -f 1-10 > "$o"
	CMD
}

function peaks::macs(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-f <size>       | assumed mean fragment
			-g <genome>     | path to
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-q <RIP/ATAC>   | true/false
			-w <broad>      | true/false peak detection
			-y <CUT>        | true/false sparse CUT&TAG/CUT&RUN data or in RIP-seq/ATAC mode (see -q) find very narrow peaks
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false ripseq=false genome threads memory maxmemory fragmentsize outdir tmpdir="${TMPDIR:-/tmp}" strict=false pointy=false broad=false
	declare -n _mapper_macs _nidx_macs _tidx_macs
	while getopts 'S:s:t:m:M:g:f:q:r:a:i:o:w:z:y:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			r)	((++mandatory)); _mapper_macs=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			q)	ripseq=$OPTARG;;
			a)	_nidx_macs=$OPTARG;;
			i)	_tidx_macs=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			w)	broad=$OPTARG;;
			y)	pointy=$OPTARG;;
			z)	strict=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage
	[[ $_nidx_macs ]] && [[ ! $_tidx_macs ]] && _usage

	declare -n _bams_macs=${_mapper_macs[0]}
	if [[ ! $_nidx_macs && ! $_tidx_macs ]]; then
		declare -a tidx_macs=("${!_bams_macs[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_macs=tidx_macs
	fi

	commander::printinfo "peak calling macs"

	# get effective genome size
	# if multimapped reads: genome minus Ns , else genome minus Ns minus repetetive Elements
	local genomesize
	local x=$(samtools view -F 4 "${_bams_macs[0]}" | head -10000 | cat <(samtools view -H "${_bams_macs[0]}") - | samtools view -c -f 256)
	if [[ $x -gt 0 ]]; then
		# genomesize=$(faCount "$genome" | tail -1 | awk '{print $3+$4+$5+$6}')
		genomesize=$(cut -f 1 "$genome.fai" | helper::lapply -t $threads -d 1 -c "samtools faidx -r /dev/stdin '$genome' | faCount /dev/stdin" | awk '/^total/{c=c+$3+$4+$5+$6}END{print c}')
	else
		# genomesize=$(unique-kmers.py -k 100 $genome |& tail -1 | awk '{print $NF}')
		declare -a gscmd=("kmc -e -r -k120 -hp -fm -t$threads -m$((maxmemory/1024>1024?1024:maxmemory/1024)) "$genome" /dev/stdout "$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.kmc)" | head -1")
		read -r genomesize genomesize < <(commander::runcmd -c kmc -i 1 -a gscmd)
	fi

	local mult=5 # macs min mem in byte according to manual is numchr * buffer=100000 * mult=2 - observation max mem used ~ mult=5
	local numchr=$(samtools view -H "${_bams_macs[0]}" | grep -c '^@SQ')
	local buffer=$(( (1024*1024*memory)/(numchr*mult) ))

	local instances instances2 ithreads ithreads2 params='' params2 cutoff
	instances=$((${#_bams_macs[@]} * ${#_mapper_macs[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -T $threads -i $instances)
	if [[ $buffer -lt 100000 ]]; then
		params="--buffer-size $buffer"
	else
		memory=$(( (numchr*buffer*mult)/1024/1024 ))
	fi
	read -r instances2 ithreads2 < <(configure::instances_by_memory -T $threads -m $memory -M "$maxmemory")

	if ! $ripseq && $pointy; then
		# according to cut&tag benchmark dont use local background at all
		# https://doi.org/10.1101/2022.03.30.486382
		# highest precision at 10^-5 but for higher recall use 10^-4
		if $strict; then
			cutoff=0.00001
		 	params+=" --nolambda -q $cutoff"
		else
			cutoff=0.0001
			params+=" --nolambda -q $cutoff"
		fi
	else
		if $strict; then
			cutoff=0.01
		 	params+=" -q $cutoff"
		else
			# -p 0.01 is encode setting for downstream IDR, but precision will go down massively
			cutoff=0.1
			params+=" -q $cutoff"
		fi
	fi
	if $broad; then
		params+=" --min-length 150 --broad --broad-cutoff $cutoff"
		# --broad option uses threshold of 0.1 which needs to be set via --broad-cutoff == -p or -q
		# further --broad increases the maximum distance between significant sites by 4 times --max-gap
		# -> --max-gap 100 --broad <=> --max-gap 400 whereas the latter reports narrow peaks with much higher signals and scores under the summit
		# -> better use broad instead of --max-gap $((fragmentsize/2*mult))
		# in case of missing control, dont use local background to not miss broad peaks due to gapped peaks in local vicinity
		# https://github.com/macs3-project/MACS/issues/467
		# -> alternative increase by 4 times --slocal 1000 -> 4000 and increase --llocal 10000 -> 40000
		[[ $_nidx_macs ]] || $pointy || params+=' --nolambda' # for cutNtag i.e. braod and point with nidx params+=' --fe-cutoff 8'
		broad="broadPeak"
		mult=4
	else
		params+=' --min-length 100 --call-summits'
		broad="narrowPeak"
		mult=1
	fi
	# for macs2 in single-end mode, the maximum gap between significant regions to join them into a single peak is the average length of the first 10 reads in the BAM!!
	# and the minimum length of a peak is the predicted fragment size
	# Genrich uses a default of 100 for joining
	# https://github.com/jsh58/Genrich/issues/27
	# use options --max-gap and --min-length. minlength should be ~500 for broad peaks..
	if $ripseq; then
		if $pointy; then
			#params+=" --max-gap $((fragmentsize/4))"
			params+=" --max-gap 50"
			params2="$params --nomodel --shift -$((fragmentsize/4)) --extsize $((fragmentsize/2))"
		else
			# works also for atac seq: shift -75 (encode) to -100 and extsize 150 (encode) to 200
			#params+=" --max-gap $((fragmentsize/2))"
			params+=" --max-gap 100"
			params2="$params --nomodel --shift -$((fragmentsize/2)) --extsize $fragmentsize"
			# params2+=' --fe-cutoff 5' # for atac
		fi
	else
		#params+=" --max-gap $((fragmentsize/2))"
		params+=" --max-gap 100"
		# if $pointy; then
		# 	# nope!
		# 	params2="$params --nomodel --shift 0 --extsize $((fragmentsize/2))"
		# else
			params2="$params --nomodel --shift 0 --extsize $fragmentsize"
		# fi
		# mfold parameter (default: --mfold 5 50) is hard to estimate for small genome sizes and or huge coverages - as of chip-exo (-m 5 100) or RIP-seq
		# the first can be tackled by downsampling, but the latter also faces strong variation in base gene expression levels, which makes model estimation nearly impossible
		$pointy && params+=' -m 5 100' || params+=' -m 5 50'
	fi

	# infer SE or PE
	# x=$(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 1)
	# [[ $x -gt 0 ]] && params+=' -f BAMPE ' || params+=' -f BAM'
	# params+=' -f BAM'
	# do never use BAMPE! the only difference to BAM is the estimation of fragment sizes by the insert of both mates
	# in both cases only R1 is used to call peaks. one may circumvent this by splitNcigar strings (i.e change PNEXT in BAMs) or use custom tagAlign/BED format solutions in conjunction with bedtools -split

	local m i f o odir nf bam
	declare -a tdirs cmd1 cmd2 cmd3 cmd4 cmd5
	for m in "${_mapper_macs[@]}"; do
		declare -n _bams_macs=$m
		odir="$outdir/$m/macs"

		for i in "${!_tidx_macs[@]}"; do
			f="${_bams_macs[${_tidx_macs[$i]}]}"
			bam="$f"

			if [[ ${_nidx_macs[$i]} ]]; then
				nf="${_bams_macs[${_nidx_macs[$i]}]}"
				o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
				# encode replaces read name by static N and score by max value 1000 to spare disk space
				# SE: bedtools bamtobed -i $nf | awk -v OFS='\t' '{$4="N"; $5="1000"; print}' | gzip -nc > $nf.tagAlign.gz # bed3+3
				# PE: bedtools bamtobed -bedpe -mate1 -i <(samtools view -F 2028 -F 256 -F 4 -f 2 $nf -u | samtools sort -n -O BAM $nf) | awk -v OFS='\t' '{print $1,$2,$3,"N",1000,$9; print $4,$5,$6,"N",1000,$10}' | gzip -nc > $nf.tagAlign.gz

				# commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				# 	samtools view -u -e 'rname!~"^chrM" && rname!="MT"' "$nf"
				# CMD
				# 	bedtools bamtobed -split -i -
				# CMD
				# 	pigz -p $ithreads -k -c > "${tdirs[-1]}/$(basename "$nf").bed.gz"
				# CMD
				# nf="-c '${tdirs[-1]}/$(basename "$nf").bed.gz'"

				# add ctrl suffix in case files are named identically
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bedtools bamtobed -split -i "$nf"
				CMD
					helper::pgzip -t $ithreads -o "$odir/$o/$(basename "${nf%.*}").ctrl.bed.gz"
				CMD
				nf="-c '$odir/$o/$(basename "${nf%.*}").ctrl.bed.gz'"
			else
				unset nf
				o="$(basename "$f")"
				o="${o%.*}"
			fi

			mkdir -p "$odir/$o"

			# commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
			# 	samtools view -u -e 'rname!~"^chrM" && rname!="MT"' "$f"
			# CMD
			# 	bedtools bamtobed -split -i -
			# CMD
			# 	pigz -p $ithreads -k -c > "${tdirs[-1]}/$(basename "$f").bed.gz"
			# CMD
			# f="${tdirs[-1]}/$(basename "$f").bed.gz"
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bedtools bamtobed -split -i "$f"
			CMD
				helper::pgzip -t $ithreads -o "$odir/$o/$(basename "${f%.*}").bed.gz"
			CMD
			f="$odir/$o/$(basename "${f%.*}").bed.gz"

			if $ripseq; then
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					printf '' > "$odir/$o/$o.model_peaks.$broad"
				CMD
					ln -sfnr "$odir/$o/nomodel_fc.bedg" "$odir/$o/model_fc.bedg"
				CMD
			else
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.macs)")
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					macs2 callpeak
					-f BED
					-t "$f"
					$nf
					-g $genomesize
					--outdir "$odir/$o"
					-n "$o.model"
					--tempdir "${tdirs[-1]}"
					-B
					--SPMR
					--keep-dup all
					--verbose 2
					$params
				CMD
					macs2 bdgcmp
					-t "$odir/$o/$o.model_treat_pileup.bdg"
					-c "$odir/$o/$o.model_control_lambda.bdg"
					-p 1
					-m FE
					-o /dev/stdout | helper::sort -k1,1 -k2,2n -k3,3n -o "$odir/$o/model_fc.bedg" -t $ithreads2 -M $memory
				CMD
			fi
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.macs)")
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				macs2 callpeak
				-f BED
				-t "$f"
				$nf
				-g $genomesize
				--outdir "$odir/$o"
				-n "$o.nomodel"
				--tempdir "${tdirs[-1]}"
				-B
				--SPMR
				--keep-dup all
				--verbose 2
				$params2
			CMD
				macs2 bdgcmp
				-t "$odir/$o/$o.nomodel_treat_pileup.bdg"
				-c "$odir/$o/$o.nomodel_control_lambda.bdg"
				-p 1
				-m FE
				-o /dev/stdout | helper::sort -k1,1 -k2,2n -k3,3n -o "$odir/$o/nomodel_fc.bedg" -t $ithreads2 -M $memory
			CMD

			# if broad: summit becomes -1 and has to be estimated afterwards
			commander::makecmd -a cmd3 -s '|' -o "$odir/$o/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD'
				sort -k1,1 -k2,2n -k3,3n "$odir/$o/$o.nomodel_peaks.$broad" "$odir/$o/$o.model_peaks.$broad"
			CMD
				awk -v OFS='\t' '$2>=0 && $3>=0{if(!$10){$10=-1} print}'
			CMD
				bedtools merge -c 4,5,6,7,8,9,10 -o distinct,max,distinct,collapse,max,max,collapse
			CMD
				perl -F'\t' -lane '
					$i=0;
					@enrich=split/,/,$F[-4];
					@summits=split/,/,$F[-1];
					for (0..$#enrich){
						$i=$_ if $enrich[$_] > $enrich[$i]
					}
					$F[-4]=$enrich[$i];
					$F[-1]=$summits[$i];
					$F[3]="peak_$.";
					print join"\t",@F;
				'
			CMD

			if [[ "$broad" == "broadPeak" ]]; then
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed2narrowpeak)")
				peaks::_bed2narrowpeak \
					-1 cmd4 \
					-2 cmd5 \
					-i "$f" \
					-f "$bam" \
					-b "$odir/$o/$o.narrowPeak" \
					-p "${tdirs[-1]}" \
					-o "$odir/$o.narrowPeak"
			else
				commander::makecmd -a cmd5 -s ';' -c {COMMANDER[0]}<<- CMD
					ln -sfnr "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
				CMD
			fi
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
	else
		commander::runcmd -v -b -i $instances -a cmd1
		commander::runcmd -c macs -v -b -i $instances2 -a cmd2
		{ head -1 "$odir/$o/$o.nomodel_peaks.$broad" "$odir/$o/$o.model_peaks.$broad" | grep -q .; } &> /dev/null || {
			rm -f "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			touch "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			return 0
		}
		commander::runcmd -v -b -i $threads -a cmd3
		commander::runcmd -c macs -v -b -i $threads -a cmd4
		commander::runcmd -v -b -i $threads -a cmd5
	fi

	return 0
}

function peaks::seacr(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of. required unless macs directory is given by -c
			-M <maxmemory>  | amount of. unless macs directory is given by -c
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-f <size>       | assumed mean fragment. required unless macs directory is given by -c
			-g <genome>     | path to. required unless macs directory is given by -c
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}" strict=false macsdir genome fragmentsize memory maxmemory
	declare -n _mapper_seacr _strandness_seacr _nidx_seacr _tidx_seacr
	while getopts 'S:s:c:m:M:f:g:t:r:a:i:o:z:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			c)	macsdir=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_seacr=$OPTARG;;
			a)	_nidx_seacr=$OPTARG;;
			i)	_tidx_seacr=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			z)	strict=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 3 ]] && _usage
	else
		[[ $mandatory -lt 6 ]] && _usage
	fi
	[[ $_nidx_seacr ]] && [[ ! $_tidx_seacr ]] && _usage

	declare -n _bams_seacr=${_mapper_seacr[0]}
	if [[ ! $_nidx_seacr && ! $_tidx_seacr ]]; then
		declare -a tidx_seacr=("${!_bams_seacr[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_seacr=tidx_seacr
	fi

	commander::printinfo "peak calling seacr"

	if [[ ! $macsdir ]]; then
		macsdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.seacr)"
		peaks::macs \
			-S false \
			-s $skip \
			-q false \
			-f $fragmentsize \
			-g "$genome" \
			-a _nidx_seacr \
			-i _tidx_seacr \
			-r _mapper_seacr \
			-t $threads \
			-m $memory \
			-M "$maxmemory" \
			-o "$macsdir" \
			-w false \
			-y false \
			-z true
	fi

	local m i f o odir mkdir nf readsbed params
	declare -a cmd1 cmd2 cmd3 cmd4 cmd5 tdirs
	for m in "${_mapper_seacr[@]}"; do
		declare -n _bams_seacr=$m
		odir="$outdir/$m/seacr"
		mdir="$macsdir/$m/macs"

		for i in "${!_tidx_seacr[@]}"; do
			f="${_bams_seacr[${_tidx_seacr[$i]}]}"
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.seacr)")

			if [[ ${_nidx_seacr[$i]} ]]; then
				nf="${_bams_seacr[${_nidx_seacr[$i]}]}"
				o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"

				# commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				# 	samtools view -u -e 'rname!~"^chrM" && rname!="MT"' "$nf"
				# CMD
				# 	bedtools genomecov -bg -ibam - > "${tdirs[-1]}/$(basename "$nf").bedg"
				# CMD
				# nf="$(basename "$nf").bedg"
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					bedtools genomecov -bg -ibam "$nf" > "${tdirs[-1]}/$(basename "${nf%.*}").ctrl.bedg"
				CMD
				nf="$(basename "${nf%.*}").ctrl.bedg"
				$strict && params="stringent" || params="relaxed"
			else
				# value between 0 and 1 as top x% of regions by area under the curve (AUC)
				# https://doi.org/10.1101/2022.03.30.486382
				# highest precision at 0.01 but for better tradoff between precision and recall use 0.02
				# at 0.02 has 90% overlap with gopeaks, whereas 0.01 as 70% at half the number of peaks
				params="stringent"
				$strict && nf="0.02" || nf="0.03"
				o="$(basename "$f")"
				o="${o%.*}"
			fi

			readsbed="$(find -L "$mdir/$o" -maxdepth 1 -name "$(basename "${f%.*}").bed.gz" -print -quit | grep .)"
			mkdir -p "$odir/$o"

			# commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			# 	samtools view -u -e 'rname!~"^chrM" && rname!="MT"' "$f"
			# CMD
			# 	bedtools genomecov -bg -ibam - > "${tdirs[-1]}/$(basename "$f").bedg"
			# CMD
			# f="${tdirs[-1]}/$(basename "$f").bedg"
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				bedtools genomecov -bg -ibam "$f" > "${tdirs[-1]}/$(basename "${f%.*}").bedg"
			CMD

			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cd "${tdirs[-1]}"
			CMD
				SEACR.sh "$(basename "${f%.*}").bedg" "$nf" norm $params "$(realpath -s "$odir/$o/$o")"
			CMD
			# norm only applied on control input

			# likewise to own summit estimation, report first maximum summit and not middle of range in case of mutliple summits
			# awk -v OFS='\t' '{print $1,$2,$3,"peak_"NR,0,".",-1,-1,-1,sprintf("%0.f",($6+($7-$6)/2-$2))}'
			commander::makecmd -a cmd3 -s '|' -o "$odir/$o/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
				sed -E 's/(.+)\t\S+:([0-9]+)-([0-9]+)$/\1\t\2\t\3/' "$odir/$o/$o.$params.bed"
			CMD
				awk -v OFS='\t' '{print $1,$2,$3,"peak_"NR,0,".",-1,-1,-1,$6-$2}'
			CMD

			peaks::_bed2narrowpeak \
				-1 cmd4 \
				-2 cmd5 \
				-i "$readsbed" \
				-f "$f" \
				-b "$odir/$o/$o.narrowPeak" \
				-p "${tdirs[-1]}" \
				-o "$odir/$o.narrowPeak"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
	else
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -c seacr -v -b -i $threads -a cmd2
		{ head -1 "$odir/$o/$o.$params.bed" | grep -q .; } &> /dev/null || {
			rm -f "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			touch "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			return 0
		}
		commander::runcmd -v -b -i $threads -a cmd3
		commander::runcmd -c macs -v -b -i $threads -a cmd4
		commander::runcmd -v -b -i $threads -a cmd5
	fi

	return 0
}

function peaks::gopeaks(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-m <memory>     | amount of. required unless macs directory is given by -c
			-M <maxmemory>  | amount of. unless macs directory is given by -c
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-t <threads>    | number of
			-f <size>       | assumed mean fragment. required unless macs directory is given by -c
			-g <genome>     | path to. required unless macs directory is given by -c
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-w <broad>      | true/false peak detection
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}" strict=false broad=false macsdir genome fragmentsize memory maxmemory
	declare -n _mapper_gopeaks _strandness_gopeaks _nidx_gopeaks _tidx_gopeaks
	while getopts 'S:s:c:m:M:f:g:t:r:a:i:o:w:z:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			c)	macsdir=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_gopeaks=$OPTARG;;
			a)	_nidx_gopeaks=$OPTARG;;
			i)	_tidx_gopeaks=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			w)	broad=$OPTARG;;
			z)	strict=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 3 ]] && _usage
	else
		[[ $mandatory -lt 6 ]] && _usage
	fi
	[[ $_nidx_gopeaks ]] && [[ ! $_tidx_gopeaks ]] && _usage

	declare -n _bams_gopeaks=${_mapper_gopeaks[0]}
	if [[ ! $_nidx_gopeaks && ! $_tidx_gopeaks ]]; then
		declare -a tidx_gopeaks=("${!_bams_gopeaks[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_gopeaks=tidx_gopeaks
	fi

	local x=$(samtools view -F 4 "${_bams_gopeaks[0]}" | head -10000 | cat <(samtools view -H "${_bams_gopeaks[0]}") - | samtools view -c -f 1)
	[[ $x -eq 0 ]] && commander::warn "peak calling gopeaks not applied due to single-end data" && return 0

	# according to dev: memory requirement is ~17X per GB bam input for default -l 50
	# y=$($broad && echo 1 || echo 2)
	local minstances mthreads memory=$(du -k "${_bams_gopeaks[0]}" | awk -v x=$([[ $_nidx_gopeaks ]] && echo 1.5 || echo 1) -v y=2 '{printf "%.f",$1/1024*17*x*y}')
	read -r minstances mthreads < <(configure::instances_by_memory -T $threads -m $memory -M "$maxmemory")

	commander::printinfo "peak calling gopeaks"

	if [[ ! $macsdir ]]; then
		macsdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gopeaks)"
		peaks::macs \
			-S false \
			-s $skip \
			-q false \
			-f $fragmentsize \
			-g "$genome" \
			-a _nidx_gopeaks \
			-i _tidx_gopeaks \
			-r _mapper_gopeaks \
			-t $threads \
			-m $memory \
			-M "$maxmemory" \
			-o "$macsdir" \
			-w false \
			-y false \
			-z true
	fi

	local params
	$strict && params="-p 0.01" || params="-p 0.05"
	# $broad && params+=" -w 150 -t $((fragmentsize/2)) -l 25 -m $((fragmentsize/2*4))" || params+=" -w 50 -t 75 -l 25 -m $((fragmentsize/2))"
	$broad && params+=" -w 150 -t 100 -l 25 -m 400" || params+=" -w 100 -t 75 -l 25 -m 100"
	# warning!!: memory usage doubles with half sliding size. better dont touch and wait for summit.py stable release
	# $pointy && params+=" -w $((fragmentsize/2))" || params="-w $fragmentsize"
	# --broad          Run GoPeaks on broad marks (--step 5000 & --slide 1000) <- for seldom extreme broad marks by h3k27me3
	# -t  --step       Bin size for coverage bins. Default: 100 <- for standard histone chip signal detection
	# -l  --slide      Slide size for coverage bins. Default: 50
	# -w  --minwidth   Minimum width (bp) of a peak. Default: 150 <- minimum in output is actually -w + -l i.e. 200 for defaults
	# -m  --mdist      Merge peaks within <mdist> base pairs. Default: 1000
	# ==> -w 150 -l 25 -p 001 (i.e. min peask length 175) has highest overlap with macs in cutntag mode

	local m i f o odir mdir nf readsbed
	declare -a tdirs cmd1 cmd2 cmd3
	for m in "${_mapper_gopeaks[@]}"; do
		declare -n _bams_gopeaks=$m
		odir="$outdir/$m/gopeaks"
		mdir="$macsdir/$m/macs"

		for i in "${!_tidx_gopeaks[@]}"; do
			f="${_bams_gopeaks[${_tidx_gopeaks[$i]}]}"

			if [[ ${_nidx_gopeaks[$i]} ]]; then
				nf="${_bams_gopeaks[${_nidx_gopeaks[$i]}]}"
				o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
				nf="-c '$nf'"
			else
				unset nf
				o="$(basename "$f")"
				o="${o%.*}"
			fi

			readsbed="$(find -L "$mdir/$o" -maxdepth 1 -name "$(basename "${f%.*}").bed.gz" -print -quit | grep .)"
			mkdir -p "$odir/$o"
			# gopeaks is not parallelized. golang per default uses #cpu concurrent goroutines => != parallelism
			# on 85mb test bam: 40:3m40 vs 4:3m43 vs 2:3m50 vs 1:4m48
			# ==> set analogous to MALLOC_ARENA_MAX=4 for single-threads java applications
			# GOMAXPROCS=4 gopeaks <- old with heavy memory usage -> re-implementation in python called summit.py <- various bugs like if r1.is_forward & r2.is_reverse object has no attribute 'is_forward'
			# OMP_NUM_THREADS=$threads summit.py
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
				GOMAXPROCS=4 gopeaks
				$nf
				-b "$f"
				-o "$odir/$o/$o"
				$params
			CMD

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed2narrowpeak)")
			peaks::_bed2narrowpeak \
				-1 cmd2 \
				-2 cmd3 \
				-i "$readsbed" \
				-f "$f" \
				-b "$odir/$o/${o}_peaks.bed" \
				-p "${tdirs[-1]}" \
				-o "$odir/$o.narrowPeak"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -c gopeaks -v -b -i $minstances -a cmd1
		{ head -1 "$odir/$o/${o}_peaks.bed" | grep -q .; } &> /dev/null || {
			rm -f "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			touch "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			return 0
		}
		commander::runcmd -c macs -v -b -i $threads -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
	fi

	return 0
}

function peaks::gem(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-f <size>       | assumed mean fragment. required unless macs directory is given by -c
			-g <genome>     | path to
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-q <RIP/ATAC>   | true/false
			-x <strandness> | hash per bam of (if ripseq)
			-y <CUT>        | true/false sparse CUT&TAG/CUT&RUN data or in RIP-seq/ATAC mode (see -q) find very narrow peaks
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false ripseq=false threads memory maxmemory genome outdir tmpdir="${TMPDIR:-/tmp}" strict=false pointy=false macsdir fragmentsize
	declare -n _mapper_gem _strandness_gem _nidx_gem _tidx_gem
	while getopts 'S:s:c:f:t:m:M:g:q:r:x:a:i:o:z:y:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			c)	macsdir=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			r)	((++mandatory)); _mapper_gem=$OPTARG;;
			a)	_nidx_gem=$OPTARG;;
			i)	_tidx_gem=$OPTARG;;
			x)	_strandness_gem=$OPTARG;;
			q)	ripseq=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			z)	strict=$OPTARG;;
			y)	pointy=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 5 ]] && _usage
	else
		[[ $mandatory -lt 6 ]] && _usage
	fi
	[[ $_nidx_gem ]] && [[ ! $_tidx_gem ]] && _usage
	# $ripseq && [[ ${#_strandness_gem[@]} -eq 0 ]] && _usage

	declare -n _bams_gem=${_mapper_gem[0]}
	if [[ ! $_nidx_gem && ! $_tidx_gem ]]; then
		declare -a tidx_gem=("${!_bams_gem[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_gem=tidx_gem
	fi

	commander::printinfo "peak calling gem"

	if [[ ! $macsdir ]]; then
		macsdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gem)"
		peaks::macs \
			-S false \
			-s $skip \
			-q false \
			-f $fragmentsize \
			-g "$genome" \
			-a _nidx_gem \
			-i _tidx_gem \
			-r _mapper_gem \
			-t $threads \
			-m $memory \
			-M "$maxmemory" \
			-o "$macsdir" \
			-w false \
			-y false \
			-z true
	fi

	local minstances mthreads jmem jgct jcgct instances=$(( ${#_tidx_gem[@]} * ${#_mapper_gem[@]} ))
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i $instances -T $threads -m $memory -M "$maxmemory")

	# get effective genome size
	# if multimapped reads: genome minus Ns , else genome minus Ns minus repetetive Elements
	local genomesize
	local x=$(samtools view -F 4 "${_bams_gem[0]}" | head -10000 | cat <(samtools view -H "${_bams_gem[0]}") - | samtools view -c -f 256)
	if [[ $x -gt 0 ]]; then
		# genomesize=$(faCount "$genome" | tail -1 | awk '{print $3+$4+$5+$6}')
		genomesize=$(cut -f 1 "$genome.fai" | helper::lapply -t $threads -d 1 -c "samtools faidx -r /dev/stdin '$genome' | faCount /dev/stdin" | awk '/^total/{c=c+$3+$4+$5+$6}END{print c}')
	else
		# genomesize=$(unique-kmers.py -k 100 $genome |& tail -1 | awk '{print $NF}')
		declare -a gscmd=("kmc -e -r -k120 -hp -fm -t$threads -m$((maxmemory/1024>1024?1024:maxmemory/1024)) "$genome" /dev/stdout "$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.kmc)" | head -1")
		read -r genomesize genomesize < <(commander::runcmd -c kmc -i 1 -a gscmd)
	fi

	local tdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.genome)"
	samtools view -H "${_bams_gem[0]}" | sed -rn '/^@SQ/{s/.+\tSN:(\S+)\s+LN:(\S+).*/\1\t\2/p}' > "$tdir/chr.info"
	declare -a cmdg
	commander::makecmd -a cmdg -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		perl -slane '
			if($_=~/^>(\S+)/){
				close F;
				open F,">$t/$1.fa" or die $!;
			}
			print F $_;
			END{
				close F;
			}
		'
	CMD
		-- -t="$tdir" "$genome"
	CMD
	if $skip; then
		commander::printcmd -a cmdg
	else
		commander::runcmd -v -b -i $threads -a cmdg
	fi

	local params='' params2='' strandness=0
	if $ripseq; then
		local x=$(samtools view -F 4 "${_bams_gem[0]}" | head -10000 | cat <(samtools view -H "${_bams_gem[0]}") - | samtools view -c -f 1)
		[[ $x -eq 0 ]] && strandness=${_strandness_gem["${_bams_gem[0]}"]}
		if [[ $strandness -ne 0 ]]; then
			params+=' --strand_type 1'
			params2+=' -s' # for bedtools merge
		fi
		# --strand_type 1 disables search for asymmetry between sense and antisense mapped reads. instead leads to strand specific peak calls based on read orientation
		# if data is paired end, mates will be analyzed individually as single end reads
		# in case of PE this leads to duplicated peaks on both strands since orientation is not inferred from first in pair read
		# -> use this parameter only in case of strand specific SE reads
		[[ $_nidx_gem ]] && params+=' --nd 2' # --nd 2 is designed for RNA based approaches and uses the control data to model the noise

		params+=" --d '$(dirname "$(which gem)")/Read_Distribution_CLIP.txt'"
		$pointy && params+=' --smooth 5'

		# if very narrow peaks i.e. more pointy signals needs to be deteced (e.g. clip or chip-exo), then
		# decrease smoothing width for read distribution estimation from default 30 bp to e.g. 5 via --smooth 5
		# if enough initial peaks can be found, re-estimation of the read distribution takes place
		# see ${out}_outputs/*.Read_Distributions.png to compare distributions.
		# Ideally, the read distribution of later rounds should be smooth and similar to that of round 0 (i.e. the read distribution model file)
		# If learned distributions shifts too much from the initial distribution one can use option --constant_model_range to use the event predictions from round 0
	else
		if $pointy; then
			params+=" --d '$(dirname "$(which gem)")/Read_Distribution_ChIP-exo.txt'"
			params+=' --smooth 5'
		else
			params+=" --d '$(dirname "$(which gem)")/Read_Distribution_default.txt'"
		fi
	fi

	if $strict; then
		params+=" -q $(echo 0.05 | awk '{print -log($1)/log(10)}')"
		# 0.05 finds ~twice more good enriched peaks than the default 0.01
	else
		params+=' -q 0 --relax'
		# encode uses -q 0 (1)- which behaves very strange: its much! more stringend than -q 1 (0.1) ~ similar to -q 3 (0.001) unless --fold 1 is used
		# in the gem source there must be some relaxation which kinda adds a small number to q-vlaue cutoff which will lead to 10^-($q+0.001) = 0.001 instead of expected 1
		# when --fold 1 is used, it behaves as expected and -q 0 leads to loose results, despite much less (~20 fold) peaks than than macs -p 0.01
		# --relax leads to more peaks, but in contrast to --nf, the overlap with macs does not increase, thus --relax adds more low scored FP for IDR
	fi
	# summary:
	# use always --fold 1 cutoff (no big difference to --fold 2 at proper q-value cutoff), to make gem behaves as expected with respect to the -q 0 exceptional behaviour
	# in combination with --nf (disables fold and shape filters), in strict mode, gem reaches highest overlap with macs strict -q 0.05 results

	local m i f o odir mdir nf readsbed
	declare -a cmd1 cmd2 cmd3 cmd4 tdirs
	for m in "${_mapper_gem[@]}"; do
		declare -n _bams_gem=$m
		odir="$outdir/$m/gem"
		mdir="$macsdir/$m/macs"

		for i in "${!_tidx_gem[@]}"; do
			f="${_bams_gem[${_tidx_gem[$i]}]}"

			if [[ ${_nidx_gem[$i]} ]]; then
				nf="${_bams_gem[${_nidx_gem[$i]}]}"
				o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
				nf="--ctrl '$nf'"
			else
				unset nf
				o="$(basename "$f")"
				o="${o%.*}"
			fi

			readsbed="$(find -L "$mdir/$o" -maxdepth 1 -name "$(basename "${f%.*}").bed.gz" -print -quit | grep .)"
			mkdir -p "$odir/$o"

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				unset DISPLAY
			CMD
				gem
					-Xmx${jmem}m
					-XX:ParallelGCThreads=$jgct
					-XX:ConcGCThreads=$jcgct
					-Djava.io.TMPDIR="$tmpdir"
					--t $mthreads
					--genome "$tdir"
					--g "$tdir/chr.info"
					--out "$odir/$o"
					--expt "$f"
					$nf
					--f SAM
					--nrf
					--sl
					--s $genomesize
					--outNP
					--fold 1
					--nf
					$params
			CMD

			# gem does not predict a true summit of its peaks
			# col 4 (Fold) is NaN unless ctr is given. else IP base mean (col 2)
			# fold enrichment is IP/Control and 99999 if control==0 ...
			# paste "$odir/$o/$o.GPS_events.narrowPeak" <(tail -n +2 "$odir/$o/$o.GPS_events.txt") | perl -lane 'if($F[13] eq "NaN"){$F[6]=$F[11]}else{$F[6]=$F[13]}; print join"\t",@F[0..9]'
			commander::makecmd -a cmd2 -s '|' -o "$odir/$o/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- 'CMD' {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- 'CMD'
				cat "$odir/$o/$o.GPS_events.narrowPeak"
			CMD
				awk -v OFS='\t' '{$7=-1; $10=-1; print}'
			CMD
				sort -k1,1 -k2,2n -k3,3n
			CMD
				sed -E '/^chr(M|X|Y|[0-9]+)/!{s/^chr//}'
			CMD
				awk -v OFS='\t' '$2>=0 && $3>=0{$4="."; $10=-1; print}'
			CMD
				bedtools merge $params2 -c 4,5,6,7,8,9,10 -o distinct,max,distinct,max,max,max,distinct
			CMD
				awk -v OFS='\t' '{$4="peak_"NR; print}'
			CMD

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed2narrowpeak)")
			peaks::_bed2narrowpeak \
				-1 cmd3 \
				-2 cmd4 \
				-i "$readsbed" \
				-f "$f" \
				-b "$odir/$o/$o.narrowPeak" \
				-p "${tdirs[-1]}" \
				-o "$odir/$o.narrowPeak"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		commander::runcmd -v -b -i $minstances -a cmd1
		{ head -1 "$odir/$o/$o.GPS_events.narrowPeak" | grep -q .; } &> /dev/null || {
			rm -f "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			touch "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			return 0
		}
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -c macs -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}

function peaks::homer(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-f <size>       | assumed mean fragment
			-g <genome>     | path to
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-q <RIP/ATAC>   | true/false
			-w <broad>      | true/false peak detection
			-x <strandness> | hash per bam of (if RIP-seq mode)
			-y <CUT>        | true/false sparse CUT&TAG/CUT&RUN data or in RIP-seq/ATAC mode (see -q) find very narrow peaks
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false ripseq=false threads memory maxmemory genome outdir tmpdir="${TMPDIR:-/tmp}" strict=false broad=false pointy=false macsdir fragmentsize
	declare -n _mapper_homer _strandness_homer _nidx_homer _tidx_homer
	while getopts 'S:s:c:f:t:m:M:g:q:r:x:a:i:o:w:z:y:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			c)	macsdir=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			r)	((++mandatory)); _mapper_homer=$OPTARG;;
			a)	_nidx_homer=$OPTARG;;
			i)	_tidx_homer=$OPTARG;;
			x)	_strandness_homer=$OPTARG;;
			q)	ripseq=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			w)	broad=$OPTARG;;
			z)	strict=$OPTARG;;
			y)	pointy=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 5 ]] && _usage
	else
		[[ $mandatory -lt 6 ]] && _usage
	fi
	[[ $_nidx_homer ]] && [[ ! $_tidx_homer ]] && _usage

	declare -n _bams_homer=${_mapper_homer[0]}
	if [[ ! $_nidx_homer && ! $_tidx_homer ]]; then
		declare -a tidx_homer=("${!_bams_homer[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_homer=tidx_homer
	fi

	local instances ithreads imemory
	instances=$((${#_bams_homer[@]} * ${#_mapper_homer[@]}))
	read -r instances imemory < <(configure::memory_by_instances -i $instances -T $threads -M "$maxmemory")
	ithreads=$((threads/instances))

	commander::printinfo "peak calling homer"

	if [[ ! $macsdir ]]; then
		macsdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.homer)"
		peaks::macs \
			-S false \
			-s $skip \
			-q false \
			-f $fragmentsize \
			-g "$genome" \
			-a _nidx_homer \
			-i _tidx_homer \
			-r _mapper_homer \
			-t $threads \
			-m $memory \
			-M "$maxmemory" \
			-o "$macsdir" \
			-w false \
			-y false \
			-z true
	fi

	# get effective genome size
	# if multimapped reads: genome minus Ns , else genome minus Ns minus repetetive Elements
	local genomesize
	local x=$(samtools view -F 4 "${_bams_homer[0]}" | head -10000 | cat <(samtools view -H "${_bams_homer[0]}") - | samtools view -c -f 256)
	if [[ $x -gt 0 ]]; then
		# genomesize=$(faCount "$genome" | tail -1 | awk '{print $3+$4+$5+$6}')
		genomesize=$(cut -f 1 "$genome.fai" | helper::lapply -t $threads -d 1 -c "samtools faidx -r /dev/stdin '$genome' | faCount /dev/stdin" | awk '/^total/{c=c+$3+$4+$5+$6}END{print c}')
	else
		# genomesize=$(unique-kmers.py -k 100 $genome |& tail -1 | awk '{print $NF}')
		declare -a gscmd=("kmc -e -r -k120 -hp -fm -t$threads -m$((maxmemory/1024>1024?1024:maxmemory/1024)) "$genome" /dev/stdout "$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.kmc)" | head -1")
		read -r genomesize genomesize < <(commander::runcmd -c kmc -i 1 -a gscmd)
	fi

	# slocal default 10000 -> 2000 works better for narrow peaks
	# narrow: -L 2 (lower cutoff for smaller local sizes) -localSize 2000 (alternative: run with slocal and llocal and merge peaks)
	# broad: -minDist 400 (-localSize 40000 if no input to get at least one p-value)
	# atac: -fragLength 0 -> don't center reads -> actually not applied for -region, but results differ..
	# cutntag: -L 0 (or -localSize 40000)
	# score always in column 7
	# if -center -L 4 : merge -c 12,13
	# if -center -L 0 : merge -c 10,11
	# if -region -L 0 : merge -c 10,11
	# if -region -L 4 : merge -c 8,9 <- attention: column names out of order when mixing -region and local filter
	#	-> is: PeakID	chr	start	end	strand	Normalized Tag Count	region size	findPeaks Score	Total Tags (normalized to Control Experiment)	Control Tags	Fold Change vs Control	p-value vs Control	Fold Change vs Local	p-value vs Local
	#	-> must: PeakID	chr	start	end	strand	Normalized Tag Count	region size	findPeaks Score	Fold Change vs Local	p-value vs Local	Total Tags (normalized to Control Experiment)	Control Tags	Fold Change vs Control	p-value vs Control
	# if input is missing: -fdr 0.00001

	local params='' params2='' strandness=0 minsize=100 maxsize=400 mindist=100
	if [[ $(wc -l < "$genome.fai") -gt 100 ]]; then
		# create single tags.tsv to favor nfs at cost of higher memory consuption
		params+=' -single'
	fi
	if ! $ripseq && $pointy; then
		# cutntag
		params2+=' -L 4 -localSize 40000'
	else
		if $broad; then
			params2+=' -L 4 -localSize 10000'
			mult=4
		else
			params2+=' -L 2 -localSize 2000'
			mult=1
		fi
	fi
	if $ripseq; then
		local x=$(samtools view -F 4 "${_bams_homer[0]}" | head -10000 | cat <(samtools view -H "${_bams_homer[0]}") - | samtools view -c -f 1)
		[[ $x -eq 0 ]] && strandness=${_strandness_homer["${_bams_homer[0]}"]}
		if [[ $strandness -ne 0 ]]; then
			params+=' -ssep'
			if [[ $strandness -eq 2 ]]; then
				params+=' -flip'
			fi
			params2+=' -strand separate'
		fi

		if $pointy; then
			maxsize=200
			mindist=50
		fi
		# prevents shifting
		fragmentsize=0
	fi
	if $broad; then
		minsize=150
		maxsize=1000
		params2+=" -minDist 400"
	else
		params2+=" -minDist $mindist"
	fi

	if [[ $_nidx_homer ]]; then
		if $strict; then
			params2+=' -fdr 0.001' # default
		else
			params2+=" -P 0.1 -LP 0.1 -fdr 0.1"
		fi
	else
		if $strict; then
			params2+=' -fdr 0.00001'
		else
			params2+=" -P 0.01 -LP 0.01 -fdr 0.01"
		fi
	fi

	local m i s f o odir mdir nf readsbed
	declare -a cmd1 cmd2 cmd3 cmd4 cmd5 tdirs
	for m in "${_mapper_homer[@]}"; do
		declare -n _bams_homer=$m
		odir="$outdir/$m/homer"
		mdir="$macsdir/$m/macs"

		for i in "${!_tidx_homer[@]}"; do
			f="${_bams_homer[${_tidx_homer[$i]}]}"

			if [[ ${_nidx_homer[$i]} ]]; then
				nf="${_bams_homer[${_nidx_homer[$i]}]}"
				o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.tagdir)")
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					ln -sfn "$(realpath -se "$nf")" "${tdirs[-1]}/$(basename "$nf")"
				CMD
					makeTagDirectory
						"$odir/$o/normal"
						$params
						-keepAll
						-precision 3
						-totalReads all
						"${tdirs[-1]}/$(basename "$nf")"
				CMD
				nf="-i '$odir/$o/normal' -inputFragLength $fragmentsize"
			else
				unset nf
				o="$(basename "$f")"
				o="${o%.*}"
			fi

			readsbed="$(find -L "$mdir/$o" -maxdepth 1 -name "$(basename "${f%.*}").bed.gz" -print -quit | grep .)"
			mkdir -p "$odir/$o"

			# can also normalize for different gc levels
			# -ssep for PE strand specific data - if inverse FR, use also -flip -> allows strand specificv peak calling via findPeaks -strand separate
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.tagdir)")
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				ln -sfn "$(realpath -se "$f")" "${tdirs[-1]}/$(basename "$f")"
			CMD
				makeTagDirectory
					"$odir/$o/treatment"
					$params
					-keepAll
					-precision 3
					-totalReads all
					"${tdirs[-1]}/$(basename "$f")"
			CMD

			# or as one-liner: seq 100 10 400 | parallel --termseq INT,1000,TERM,0 --halt now,fail=1 --line-buffer -P 40 findPeaks TEST/ip/ -i TEST/input/ -gsize 2801039415 -region -L 2 -localSize 2000 -size {} -minDist 100 -fragLength 150 -inputFragLength 150 -C 0 | grep -v '^#' | cut -f 2- | awk -v OFS='\t' '$4="name\t"$7"\t"$4' | psort -k1,1 -k2,2n -k3,3n | bedtools merge -s -i - -c 4,5,6,10,11 -o distinct,max,distinct,max,min | awk -v s=$(${ss:-false} && echo 1 || echo 0) -v OFS='\t' '{$4="peak_"NR; if(!s){$6="."} print $0,-1,-1}' > homer.narrowPeak
			for s in $(seq $minsize 10 $maxsize); do
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					findPeaks
						"$odir/$o/treatment"
						$nf
						-fragLength $fragmentsize
						-gsize $genomesize
						-region
						-size $s
						-o "$odir/$o/size_$s.tsv"
						$params2
						-C 0
				CMD
			done

			if [[ $strandness -eq 0 ]]; then
				commander::makecmd -a cmd3 -s '|' -o "$odir/$o/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- 'CMD'
					grep -v '^#' "$odir/$o/size_"*.tsv
				CMD
					cut -f 2-
				CMD
					awk -v OFS='\t' '$4="name\t"$7"\t"$4'
				CMD
					helper::sort -t $ithreads -M $imemory -k1,1 -k2,2n -k3,3n
				CMD
					bedtools merge -s -i - -c 4,5,6,10,11 -o distinct,max,distinct,max,min
				CMD
					awk -v OFS='\t' '{$4="peak_"NR; $6="."; p=-log($8)/log(10); if(p=="inf"){$8=999}else{$8=p} print $0,-1,-1}'
				CMD
				# after awk:
				# (chr	start	end)	name	findPeaks_Score	strand	Normalized_Tag_Count	region_size	findPeaks_Score	Fold_Change_vs_Local	p-value_vs_Local	Total_Tags_(normalized_to_Control_Experiment)	Control_Tags	Fold_Change_vs_Control	p-value_vs_Control
				# after merge
				# (chr	start	end)	name	findPeaks_Score	strand	Fold_Change_vs_Local	p-value_vs_Local
			else
				commander::makecmd -a cmd3 -s '|' -o "$odir/$o/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- 'CMD'
					grep -v '^#' "$odir/$o/size_"*.tsv
				CMD
					cut -f 2-
				CMD
					awk -v OFS='\t' '$4="name\t"$7"\t"$4'
				CMD
					helper::sort -t $ithreads -M $imemory -k1,1 -k2,2n -k3,3n
				CMD
					bedtools merge -s -i - -c 4,5,6,10,11 -o distinct,max,distinct,max,min
				CMD
					awk -v OFS='\t' '{$4="peak_"NR; p=-log($8)/log(10); if(p=="inf"){$8=999}else{$8=p} print $0,-1,-1}'
				CMD
			fi

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed2narrowpeak)")
			peaks::_bed2narrowpeak \
				-1 cmd4 \
				-2 cmd5 \
				-i "$readsbed" \
				-f "$f" \
				-b "$odir/$o/$o.narrowPeak" \
				-p "${tdirs[-1]}" \
				-o "$odir/$o.narrowPeak"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
	else
		commander::runcmd -c homer -v -b -i $threads -a cmd1
		commander::runcmd -c homer -v -b -i $threads -a cmd2
		grep -q -v -m 1 '^#' "$odir/$o/size_$minsize.tsv" &> /dev/null || {
			rm -f "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			touch "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			return 0
		}
		commander::runcmd -v -b -i $instances -a cmd3
		commander::runcmd -c macs -v -b -i $threads -a cmd4
		commander::runcmd -v -b -i $threads -a cmd5
	fi

	return 0
}

function peaks::matk(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-f <size>       | assumed mean fragment. required unless macs directory is given by -c
			-g <genome>     | path to. required unless macs directory is given by -c
			-t <threads>    | number of
			-m <memory>     | amount of. required unless macs directory is given by -c
			-M <maxmemory>  | amount of
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-i <tidx>       | array of RIP* bam idices within -r
			-o <outdir>     | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory outdir tmpdir="${TMPDIR:-/tmp}" macsdir genome fragmentsize memory
	declare -n _mapper_matk _nidx_matk _tidx_matk _ridx_matk _pidx_matk
	while getopts 'S:s:c:f:m:g:t:M:r:x:a:i:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			c)	macsdir=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			t)	((++mandatory)); threads=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			r)	((++mandatory)); _mapper_matk=$OPTARG;;
			a)	((++mandatory)); _nidx_matk=$OPTARG;; # contains nridx by alignment::mkreplicates
			i)	((++mandatory)); _tidx_matk=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 5 ]] && _usage
	else
		[[ $mandatory -lt 8 ]] && _usage
	fi
	[[ $_nidx_matk ]] && [[ ! $_tidx_matk ]] && _usage

	declare -n _bams_matk=${_mapper_matk[0]}
	if [[ ! $_nidx_matk && ! $_tidx_matk ]]; then
		declare -a tidx_matk=("${!_bams_matk[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_matk=tidx_matk
	fi

	commander::printinfo "peak calling matk"

	if [[ ! $macsdir ]]; then
		macsdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.matk)"
		peaks::macs \
			-S false \
			-s $skip \
			-q false \
			-f $fragmentsize \
			-g "$genome" \
			-a _nidx_gopeaks \
			-i _tidx_gopeaks \
			-r _mapper_gopeaks \
			-t $threads \
			-m $memory \
			-M "$maxmemory" \
			-o "$macsdir" \
			-w false \
			-y false \
			-z true
	fi

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i 1 -T $threads -M "$maxmemory")

	local m i f o odir mdir nf readsbed
	declare -a cmd1 cmd2 cmd3 tdirs
	for m in "${_mapper_matk[@]}"; do
		declare -n _bams_matk=$m
		odir="$outdir/$m/matk"
		mdir="$macsdir/$m/macs"

		for i in "${!_tidx_matk[@]}"; do
			nf="${_bams_matk[${_nidx_matk[$i]}]}"
			f="${_bams_matk[${_tidx_matk[$i]}]}"
			o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"

			readsbed="$(find -L "$mdir/$o" -maxdepth 1 -name "$(basename "${f%.*}").bed.gz" -print -quit | grep .)"
			mkdir -p "$odir/$o"

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				OMP_NUM_THREADS=$threads matk
					-Xmx${jmem}m
					-XX:ParallelGCThreads=$jgct
					-XX:ConcGCThreads=$jcgct
					-Djava.io.TMPDIR="$tmpdir"
					-peakCalling
					-c 1
					-ip "$f"
					-input "$nf"
					-out "$odir/$o/$o.bed"
			CMD
				perl -lane 'print join"\t",(@F[0..3],0,".",-1,-1,-1*log(\$F[4])/log(10),-1)' "$odir/$o/$o.bed" | helper::sort -k1,1 -k2,2n -k3,3n -o "$odir/$o/$o.narrowPeak" -t $threads -M $jmem
			CMD

			# commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- 'CMD'
			# 	samtools bedcov -Q 0 -g 1796 -j "$odir/$o/$o.bed" "$nf" "$f"
			# CMD
			# 	awk -v c=$(samtools idxstats "$nf" | datamash sum 3) -v i=$(samtools idxstats "$f" | datamash sum 3) '\$NF>0 && \$NF>\$(NF-1){print \$0,(\$NF/i)/((\$(NF-1)+1)/c)}'
			# CMD
			# 	perl -lane 'print join"\t",(@F[0..3],0,".",$F[-1],-1,-1*log($F[4])/log(10),-1)'
			# CMD
			# 	sort -k1,1 -k2,2n -k3,3n
			# CMD
			# 	bedtools merge -c 4,5,6,7,8,9,10 -o distinct,max,distinct,max,max,max,distinct
			# CMD
			# 	awk -v OFS='\t' '{$4="peak_"NR; print}'
			# CMD

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed2narrowpeak)")
			peaks::_bed2narrowpeak \
				-1 cmd2 \
				-2 cmd3 \
				-i "$readsbed" \
				-f "$f" \
				-b "$odir/$o/$o.narrowPeak" \
				-p "${tdirs[-1]}" \
				-o "$odir/$o.narrowPeak"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -i 1 -a cmd1
		{ head -1 "$odir/$o/$o.narrowPeak" | grep -q .; } &> /dev/null || {
			rm -f "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			touch "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			return 0
		}
		commander::runcmd -c macs -v -b -i $threads -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
	fi

	return 0
}

function peaks::peakachu(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-g <genome>     | path to. required unless macs directory is given by -c
			-m <memory>     | amount of. required unless macs directory is given by -c
			-M <maxmemory>  | amount of
			-t <threads>    | number of
			-f <size>       | assumed mean fragment
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-z <strict>     | true/false peak filters
			-P              | run pair-wise i.e. ignores replicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads fragmentsize outdir strict=false genome macsdir memory maxmemory pairwise=false
	declare -n _mapper_peakachu _nidx_peakachu _tidx_peakachu
	while getopts 'S:s:c:g:m:M:t:f:r:a:i:o:z:P' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			c)	macsdir=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_peakachu=$OPTARG;;
			a)	_nidx_peakachu=$OPTARG;;
			i)	_tidx_peakachu=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			z)	strict=$OPTARG;;
			P)	pairwise=true;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 4 ]] && _usage
	else
		[[ $mandatory -lt 6 ]] && _usage
	fi
	[[ $_nidx_peakachu ]] && [[ ! $_tidx_peakachu ]] && _usage

	declare -n _bams_peakachu=${_mapper_peakachu[0]}
	if [[ ! $_nidx_peakachu && ! $_tidx_peakachu ]]; then
		declare -a tidx_peakachu=("${!_bams_peakachu[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_peakachu=tidx_peakachu
	fi

	commander::printinfo "peak calling peakachu"

	if $pairwise || [[ ${#_tidx_peakachu[@]} -eq 1 ]] && [[ ! $macsdir ]]; then
		macsdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.peakachu)"
		peaks::macs \
			-S false \
			-s $skip \
			-q false \
			-f $fragmentsize \
			-g "$genome" \
			-a _nidx_peakachu \
			-i _tidx_peakachu \
			-r _mapper_peakachu \
			-t $threads \
			-m $memory \
			-M "$maxmemory" \
			-o "$macsdir" \
			-w false \
			-y false \
			-z true
	fi
	# todo: else merge bams for single macs run to infer lambda track to get macs derived fc

	local x params
	# infer SE or PE
	x=$(samtools view -F 4 "${_bams_peakachu[0]}" | head -10000 | cat <(samtools view -H "${_bams_peakachu[0]}") - | samtools view -c -f 1)
	[[ $x -gt 0 ]] && params='--paired_end'

	local m f i nf o a b odir mdir readsbed
	declare -a cmd1 cmd2 cmd3 cmd4 tdirs
	for m in "${_mapper_peakachu[@]}"; do
		declare -n _bams_peakachu=$m
		odir="$outdir/$m/peakachu"
		mdir="$macsdir/$m/macs"

		if $pairwise || [[ ${#_tidx_peakachu[@]} -eq 1 ]]; then
			for i in "${!_tidx_peakachu[@]}"; do
				f="${_bams_peakachu[${_tidx_peakachu[$i]}]}"

				if [[ ${_nidx_peakachu[$i]} ]]; then
					nf="${_bams_peakachu[${_nidx_peakachu[$i]}]}"
					o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
					nf="--ctr_libs '$nf'"
				else
					unset nf
					o="$(basename "$f")"
					o="${o%.*}"
				fi

				readsbed="$(find -L "$mdir/$o" -maxdepth 1 -name "$(basename "${f%.*}").bed.gz" -print -quit | grep .)"

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					rm -rf "$odir/$o"
				CMD
					peakachu adaptive
						$params
						--norm_method none
						--exp_libs "$f"
						$nf
						--max_insert_size $((2*fragmentsize))
						--max_proc $threads
						--mad_multiplier 0.0
						--fc_cutoff 2
						--padj_threshold $($strict && echo 0.05 || echo 0.1)
						--output_folder "$odir/$o"
				CMD

				commander::makecmd -a cmd2 -s '|' -o "$odir/$o/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- 'CMD'
					awk -v OFS='\t' '!/^replicon/{print \$1,\$3,\$4,".",0,".",-1,-1,-1,-1}' "$odir/$o/peak_tables/"*.csv
				CMD
					sort -k1,1 -k2,2n -k3,3n
				CMD
					{ read -r l; [[ $l ]] && cat <(echo "$l") - | bedtools merge -c 4,5,6,7,8,9,10 -o distinct,distinct,distinct,max,distinct,distinct,distinct; }
				CMD
					awk -v OFS='\t' '{$4="peak_"NR; print}'
				CMD

				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed2narrowpeak)")
				peaks::_bed2narrowpeak \
					-1 cmd3 \
					-2 cmd4 \
					-i "$readsbed" \
					-f "$f" \
					-b "$odir/$o/$o.narrowPeak" \
					-p "${tdirs[-1]}" \
					-o "$odir/$o.narrowPeak"
			done
		else
			f='--exp_libs'
			for i in "${!_tidx_peakachu[@]}"; do
				f+="$(printf ' "%s"' "${_bams_peakachu[${_tidx_peakachu[$i]}]}")"
			done
			if [[ $_nidx_peakachu ]]; then
				nf="--ctr_libs"
				for i in "${!_nidx_peakachu[@]}"; do
					nf+="$(printf ' "%s"' "${_bams_peakachu[${_nidx_peakachu[$i]}]}")"
				done
				a="${_bams_peakachu[${_nidx_peakachu[0]}]}"
				b="${_bams_peakachu[${_tidx_peakachu[0]}]}"
				o="$(echo -e "$(basename "$a")\t$(basename "$b")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
			else
				unset nf
				o="$(basename "${_bams_peakachu[${_tidx_peakachu[0]}]}")"
				o="${o%.*}"
			fi

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				rm -rf "$odir/$o"
			CMD
				peakachu adaptive
					$params
					--norm_method deseq
					$f
					$nf
					--max_insert_size $((2*fragmentsize))
					--max_proc $threads
					--mad_multiplier 0.0
					--fc_cutoff 2
					--padj_threshold $($strict && echo 0.05 || echo 0.1)
					--output_folder "$odir/$o"
			CMD

			# NF-9 is either base_mean or fold_change, the latter if ctr_libs given
			if [[ $nf ]]; then
				commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- 'CMD'
					awk -v OFS='\t' '!/^replicon/ && \$(NF-9)!="inf"{print \$1,\$3,\$4,".",0,".",\$(NF-9),-1,-1,-1}' "$odir/$o/peak_tables/"*.csv
				CMD
					sort -k1,1 -k2,2n -k3,3n
				CMD
					{ read -r l; [[ $l ]] && cat <(echo "$l") - | bedtools merge -c 4,5,6,7,8,9,10 -o distinct,distinct,distinct,max,distinct,distinct,distinct; }
				CMD
					awk -v OFS='\t' '{$4="peak_"NR; print}'
				CMD
			else
				commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- 'CMD'
					awk -v OFS='\t' '!/^replicon/{print \$1,\$3,\$4,".",0,".",-1,-1,-1,-1}' "$odir/$o/peak_tables/"*.csv
				CMD
					sort -k1,1 -k2,2n -k3,3n
				CMD
					{ read -r l; [[ $l ]] && cat <(echo "$l") - | bedtools merge -c 4,5,6,7,8,9,10 -o distinct,distinct,distinct,max,distinct,distinct,distinct; }
				CMD
					awk -v OFS='\t' '{$4="peak_"NR; print}'
				CMD
			fi
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		commander::runcmd -c peakachu -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -c macs -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}

function peaks::genrich(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-f <size>       | assumed mean fragment
			-t <threads>    | number of
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-q <RIP/ATAC>   | true/false
			-w <broad>      | true/false peak detection
			-y <CUT>        | true/false sparse CUT&TAG/CUT&RUN data or in RIP-seq/ATAC mode (see -q) find very narrow peaks
			-z <strict>     | true/false peak filters
			-P              | run pair-wise i.e. ignores replicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false genome fragmentsize threads outdir tmpdir="${TMPDIR:-/tmp}" strict=false ripseq=false pairwise=false broad=false pointy=false
	declare -n _mapper_genrich _nidx_genrich _tidx_genrich
	while getopts 'S:s:f:t:r:a:i:o:q:z:y:w:P' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_genrich=$OPTARG;;
			a)	_nidx_genrich=$OPTARG;;
			i)	_tidx_genrich=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			q)	ripseq=$OPTARG;;
			z)	strict=$OPTARG;;
			y)	pointy=$OPTARG;;
			w)	broad=$OPTARG;;
			P)	pairwise=true;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage
	[[ $_nidx_genrich ]] && [[ ! $_tidx_genrich ]] && _usage

	declare -n _bams_genrich=${_mapper_genrich[0]}
	if [[ ! $_nidx_genrich && ! $_tidx_genrich ]]; then
		declare -a tidx_genrich=("${!_bams_genrich[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_genrich=tidx_genrich
	fi

	commander::printinfo "peak calling genrich"

	local params
	# f=200
	# for single sample ATAC, -a 50 or -a 0 (turns off filter for AUC) results are pretty much comparable with macs peaks when -d f/2 -g f/4
	# for ATAC relicates, -a set to -d (f/2) outcome is close to intersection of single runs with auc filter turned off (-a 0)
	# todo in case of RIP-seq on transcriptome, where multimapped reads have to be kept, prevent downsampling if #secondaries exceed 10 by manipulating sam? test it!
	# filter relaxed for poiny & ! ripsec i.e sparse cutntag data

	if $pairwise || [[ ${#_tidx_genrich[@]} -eq 1 ]]; then
		$ripseq && params="-a 100 -j -D -d 100" || params="-a 100"
	else
		$ripseq && params="-a 200 -j -D -d 100" || params="-a 400"
	fi
	if $broad; then
		params+=" -l 150 -g 400"
	else
		$ripseq && $pointy && params+=" -l 100 -g 50" || params+=" -l 100 -g 100"
	fi
	$strict && params+=' -p 0.01' || params+=' -p 0.05'

	local m f a b i nf o odir
	declare -a tdirs cmd1 cmd2 cmd3
	for m in "${_mapper_genrich[@]}"; do
		declare -n _bams_genrich=$m
		odir="$outdir/$m/genrich"
		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.genrich)")

		if $pairwise || [[ ${#_tidx_genrich[@]} -eq 1 ]]; then
			for i in "${!_tidx_genrich[@]}"; do
				f="${_bams_genrich[${_tidx_genrich[$i]}]}"

				if [[ ${_nidx_genrich[$i]} ]]; then
					nf="${_bams_genrich[${_nidx_genrich[$i]}]}"
					o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
					b="$(basename "${nf%.*}")"

					commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
						samtools sort -@ $threads -O BAM -T "${tdirs[-1]}/$b" -n -o "${tdirs[-1]}/$b.ctrl.bam" "$nf"
					CMD

					nf="-c '${tdirs[-1]}/$b.ctrl.bam'"
				else
					unset nf
					o="$(basename "$f")"
					o="${o%.*}"
				fi

				b="$(basename "${f%.*}")"
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					samtools sort -@ $threads -O BAM -T "${tdirs[-1]}/$b" -n -o "${tdirs[-1]}/$b.bam" "$f"
				CMD

				f="-t '${tdirs[-1]}/$b.bam'"
			done
		else
			if [[ $_nidx_genrich ]]; then
				nf='-c '
				for i in "${!_nidx_genrich[@]}"; do
					b="${_bams_genrich[${_nidx_genrich[$i]}]}"
					b="$(basename "${b%.*}")"

					commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
						samtools sort -@ $threads -O BAM -T "${tdirs[-1]}/$b" -n -o "${tdirs[-1]}/$b.ctrl.bam" "${_bams_genrich[${_nidx_genrich[$i]}]}"
					CMD

					nf+="$(printf '"%s",' "${tdirs[-1]}/$b.ctrl.bam")"
				done
				nf="${nf:0:-1}"
				a="${_bams_genrich[${_nidx_genrich[0]}]}"
				b="${_bams_genrich[${_tidx_genrich[0]}]}"
				o="$(echo -e "$(basename "$a")\t$(basename "$b")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
			else
				unset nf
				o="$(basename "${_bams_genrich[${_tidx_genrich[0]}]}")"
				o="${o%.*}"
			fi

			f='-t '
			for i in "${!_tidx_genrich[@]}"; do
				b="${_bams_genrich[${_tidx_genrich[$i]}]}"
				b="$(basename "${b%.*}")"

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					samtools sort -@ $threads -O BAM -T "${tdirs[-1]}/$b" -n -o "${tdirs[-1]}/$b.bam" "${_bams_genrich[${_tidx_genrich[$i]}]}"
				CMD

				f+="$(printf '"%s",' "${tdirs[-1]}/$b.bam")"
			done
			f="${f:0:-1}"
		fi

		mkdir -p "$odir/$o"

		commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
			Genrich
				$f
				$nf
				$params
				-y
				-v
				-o "$odir/$o/$o.narrowPeak"
		CMD

		commander::makecmd -a cmd3 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			sort -k1,1 -k2,2n -k3,3n "$odir/$o/$o.narrowPeak"
		CMD
			awk -v OFS='\t' '$2>=0 && $3>=0{$4="peak_"NR; print}'
		CMD
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -i 1 -a cmd1
		commander::runcmd -c genrich -v -b -i $threads -a cmd2
		{ head -1 "$odir/$o/$o.narrowPeak" | grep -q .; } &> /dev/null || {
			rm -f "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			touch "$odir/$o/$o.narrowPeak" "$odir/$o.narrowPeak"
			return 0
		}
		commander::runcmd -v -b -i $threads -a cmd3
	fi

	return 0
}

function peaks::m6aviewer(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-f <size>       | assumed mean fragment
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false genome threads memory maxmemory fragmentsize outdir tmpdir="${TMPDIR:-/tmp}"
	declare -n _mapper_m6aviewer _nidx_m6aviewer _tidx_m6aviewer
	while getopts 'S:s:t:m:M:f:r:a:b:i:j:k:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			r)	((++mandatory)); _mapper_m6aviewer=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			a)	_nidx_m6aviewer=$OPTARG;;
			i)	_tidx_m6aviewer=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage
	[[ $_nidx_m6aviewer ]] && [[ ! $_tidx_m6aviewer ]] && _usage

	declare -n _bams_m6aviewer=${_mapper_m6aviewer[0]}
	if [[ ! $_nidx_m6aviewer && ! $_tidx_m6aviewer ]]; then
		declare -a tidx_m6aviewer=("${!_bams_m6aviewer[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_m6aviewer=tidx_m6aviewer
	fi

	commander::printinfo "peak calling m6aviewer"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i 1 -T $threads -m $memory -M "$maxmemory")

	local m i f o odir nf
	declare -a cmd1 cmd2
	for m in "${_mapper_m6aviewer[@]}"; do
		declare -n _bams_m6aviewer=$m
		odir="$outdir/$m/m6aviewer"
		mkdir -p $odir

		commander::printinfo {COMMANDER[0]}<<- EOF
			HOWTO
			1a) settings - peak calling: minimum enrichment: 2, minimum peak height: 50, expected peak width=$((2*fragmentsize))
			1b) settings - other: threads: $threads
			2) load blocks of control and IP files as pairs into m6aviewer as printed below
			3) don't group anything, uncheck "limit to" and click "find peaks"
			4) file - save peaks to file "$odir/results"
		EOF

		for i in "${!_tidx_m6aviewer[@]}"; do
			f="${_bams_m6aviewer[${_tidx_m6aviewer[$i]}]}"

			if [[ $_nidx_m6aviewer ]]; then
				nf="${_bams_m6aviewer[${_nidx_m6aviewer[$i]}]}"
				commander::printinfo {COMMANDER[0]}<<- EOF
					load
					$nf
					versus
					$f
				EOF
				o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/').narrowPeak"
			else
				commander::printinfo {COMMANDER[0]}<<- EOF
					load
					$f
				EOF
				o="$(basename "$f")"
				o="${o%.*}.narrowPeak"
			fi

			commander::makecmd -a cmd2 -s '|' -o "$odir/$o" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- 'CMD'
				cat "$odir/results_$(basename "$f").txt"
			CMD
				awk -v OFS='\t' '!/^\s*$/ && NR>2{p=-log($5)/log(10); print $3,$4-75,$4+76,".",0,".",$6,$5,-1,75}'
			CMD
				sort -k1,1 -k2,2n -k3,3n
			CMD
				bedtools merge -c 4,5,6,7,8,9,10 -o distinct,max,distinct,collapse,max,max,collapse
			CMD
				perl -F'\t' -lane '
					$i=0;
					@enrich=split/,/,$F[-4];
					@summits=split/,/,$F[-1];
					for (0..$#enrich){
						$i=$_ if $enrich[$_] > $enrich[$i]
					}
					$F[-4]=$enrich[$i];
					$F[-1]=$summits[$i];
					$F[3]="peak_$.";
					print join"\t",@F;
				'
			CMD
		done

		#cmd1+=("m6aviewer -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.TMPDIR='$tmpdir'")
		m6aviewer -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.TMPDIR="$tmpdir" &> /dev/null
	done

	if $skip; then
		# commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		# commander::runcmd -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}
