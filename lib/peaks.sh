#! /usr/bin/env bash
# (c) Konstantin Riege

# https://sites.google.com/site/anshulkundaje/projects/idr
# https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit
# https://sites.google.com/site/anshulkundaje/projects/idr/deprecated

peaks::_idr(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-t <file>     | normal vs treatment peaks path to
			-r <file>     | normal vs replikate peaks path to
			-p <file>     | oracle: normal vs pool peaks path to
			-o <outfile>  | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory t r p o
	declare -n _cmds1_idr _cmds2_idr
	while getopts '1:2:t:r:p:o:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_idr=$OPTARG;;
			2)	((++mandatory)); _cmds2_idr=$OPTARG;;
			t)	((++mandatory)); t=$OPTARG;;
			r)	((++mandatory)); r=$OPTARG;;
			p)	((++mandatory)); p=$OPTARG;;
			o)	((++mandatory)); o="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage

	#p.value q.value signal.value
	commander::makecmd -a _cmds1_idr -s ';' -c {COMMANDER[0]}<<- CMD
		idr
		--samples "$t" "$r"
		--peak-list "$p"
		--input-file-type narrowPeak
		--output-file "$o"
		--rank signal.value
		--soft-idr-threshold 0.05
		--plot
		--use-best-multisummit-IDR
	CMD

	# Columns 1-10 are same as pooled common peaks narrowPeak columns
	# Col 11: -log10(local IDR value)
	# Col 12: -log10(global IDR value)
	# Col 15: ranking measure from Rep1
	# Col 19: ranking measure from Rep2

	local t=$(echo 0.05 | awk '{print -log($1)/log(10)}')
	commander::makecmd -a _cmds2_idr -s '|' -o "$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
		awk -v t=$t '\$12>=t' "$o"
	CMD
		cut -f 1-10
	CMD
		sort -k1,1 -k2,2n -k3,3n
	CMD
		uniq
	CMD

	commander::makecmd -a _cmds2_idr -s '|' -o "$o.full.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
		cut -f 1-10 "$o"
	CMD
		sort -k1,1 -k2,2n -k3,3n
	CMD
		uniq
	CMD

	return 0
}

peaks::macs(){
	declare -a tdirs
	_cleanup::peaks::macs(){
		rm -rf "${tdirs[@]}"
	}

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-f <size>       | assumed mean fragment
			-g <genome>     | path to
			-q <ripseq>     | true/false
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-p <tmpdir>     | path to
			-y <pointy>     | true/false call only pointy peaks
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false ripseq=false genome threads memory maxmemory fragmentsize outdir tmpdir strict=false pointy=false
	declare -n _mapper_macs _nidx_macs _tidx_macs
	while getopts 'S:s:t:m:M:g:f:q:r:a:i:o:p:z:y:' arg; do
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
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			y)	pointy=$OPTARG;;
			z)	strict=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage
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
		genomesize=$(faCount $genome | tail -1 | awk '{print $3+$4+$5+$6}')
	else
		genomesize=$(unique-kmers.py -k 100 $genome 2>&1 | tail -1 | awk '{print $NF}')
	fi

	local params='' mult=5 # macs min mem in byte according to manual is numchr * buffer=100000 * mult=2 - observation max mem used ~ mult=5
	local numchr=$(samtools view -H "${_bams_macs[0]}" | grep -c '^@SQ')
	local buffer=$(( (1024*1024*memory)/(numchr*mult) ))

	local instances instances2 ithreads ithreads2
	instances=$((${#_bams_macs[@]} * ${#_mapper_macs[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -T $threads -i $instances)
	if [[ $buffer -lt 100000 ]]; then
		params+=" --buffer-size $buffer"
		read -r instances2 ithreads2 < <(configure::instances_by_memory -T $threads -m $memory -M "$maxmemory")
	else
		read -r instances2 ithreads2 < <(configure::instances_by_memory -T $threads -m $(( (numchr*buffer*mult)/1024/1024 )) -M "$maxmemory")
	fi

	if $strict; then
		params+=' -q 0.05' # macs default. might be changed to 0.01
	else
		params+=' -p 0.01' # encode setting for downstream IDR. might be changed to -q 0.1
	fi

	# infer SE or PE
	# x=$(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 1)
	# [[ $x -gt 0 ]] && params+=' -f BAMPE ' || params+=' -f BAM'
	# params+=' -f BAM'
	# do never use BAMPE! the only difference to BAM is the estimation of fragment sizes by the insert of both mates
	# in both cases only R1 is used to call peaks. one may circumvent this by splitNcigar strings (i.e change PNEXT in BAMs) or use custom tagAlign/BED format solutions in conjunction with bedtools -split

	local m i f o odir nf
	declare -a cmd1 cmd2 cmd3
	for m in "${_mapper_macs[@]}"; do
		declare -n _bams_macs=$m
		odir="$outdir/$m/macs"

		for i in "${!_tidx_macs[@]}"; do
			f="${_bams_macs[${_tidx_macs[$i]}]}"

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.macs)")

			if [[ ${_nidx_macs[$i]} ]]; then
				nf="${_bams_macs[${_nidx_macs[$i]}]}"
				o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"

				# encode replaces read name by static N and score by max value 1000 to spare disk space
				# SE: bedtools bamtobed -i $nf | awk -v OFS='\t' '{$4="N"; $5="1000"; print}' | gzip -nc > $nf.tagAlign.gz # bed3+3
				# PE: bedtools bamtobed -bedpe -mate1 -i <(samtools view -F 2028 -F 256 -F 4 -f 2 $nf -u | samtools sort -n -O BAM $nf) | awk -v OFS='\t' '{print $1,$2,$3,"N",1000,$9; print $4,$5,$6,"N",1000,$10}' | gzip -nc > $nf.tagAlign.gz

				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bedtools bamtobed -split -i "$nf"
				CMD
					pigz -p $ithreads -k -c > "${tdirs[-1]}/$(basename "$nf").bed.gz"
				CMD
				nf="-c '${tdirs[-1]}/$(basename "$nf").bed.gz'"
			else
				unset nf
				o="$(basename "$f")"
				o="${o%.*}"
			fi

			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bedtools bamtobed -split -i "$f"
			CMD
				pigz -p $ithreads -k -c > "${tdirs[-1]}/$(basename "$f").bed.gz"
			CMD
			f="${tdirs[-1]}/$(basename "$f").bed.gz"

			mkdir -p "$odir/$o"
			printf '' > "$odir/$o/$o.model_peaks.narrowPeak"
			if $ripseq; then # works also for atac seq: shift -75 (encode) to -100 and extsize 150 (encode) to 200
				# mfold parameter (default: --mfold 5 50) is hard to estimate for small genome sizes and or huge coverages - as of chip-exo (-m 5 100) or RIP-seq
				# the first can be tackled by downsampling, but the latter also faces strong variation in base gene expression levels, which makes model estimation nearly impossible
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					macs2 callpeak
					-f BED
					-t "$f"
					$nf
					-g $genomesize
					--outdir "$odir/$o"
					-n "$o.nomodel"
					--tempdir "${tdirs[-1]}"
					--call-summits
					-B
					--SPMR
					--keep-dup all
					--verbose 2
					--nomodel
					--shift $($pointy && echo "-$((fragmentsize/8))" || echo "-$((fragmentsize/2))")
					--extsize $($pointy && echo "$((fragmentsize/2))" || echo $fragmentsize)
					$params
				CMD
			else
				# pointy for chip-exo: shift 0, extsize 30-50, bw 150, m 5 100-200
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					macs2 callpeak
					-f BED
					-t "$f"
					$nf
					-g $genomesize
					--outdir "$odir/$o"
					-n "$o.model"
					--tempdir "${tdirs[-1]}"
					--call-summits
					-B
					--SPMR
					--keep-dup all
					--verbose 2
					$params
					$($pointy && echo "-m 5 150 --bw 150")
				CMD

				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.macs)")
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					macs2 callpeak
					-f BED
					-t "$f"
					$nf
					-g $genomesize
					--outdir "$odir/$o"
					-n "$o.nomodel"
					--tempdir "${tdirs[-1]}"
					--call-summits
					-B
					--SPMR
					--keep-dup all
					--verbose 2
					--nomodel
					--shift 0
					--extsize $($pointy && echo "$((fragmentsize/4))" || echo $fragmentsize)
					$params
				CMD
			fi

			# keep summit of peak with highest summit enrichment
			commander::makecmd -a cmd3 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD'
				sort -k1,1 -k2,2n -k3,3n "$odir/$o/$o.nomodel_peaks.narrowPeak" "$odir/$o/$o.model_peaks.narrowPeak"
			CMD
				awk '$2>=0 && $3>=0'
			CMD
				bedtools merge -c 4,5,6,7,8,9,10 -o distinct,max,distinct,collapse,max,max,collapse
			CMD
				perl -F'\t' -lane '
					$i=0;
					@enrich=split/,/,$F[-4];
					@summits=split/,/,$F[-1];
					for (0..$#enrich){
						$i=$_ if $enrich[$_] < $enrich[$i]
					}
					$F[-4]=$enrich[$i];
					$F[-1]=$summits[$i];
					$F[3]="peak_$.";
					print join"\t",@F;
				'
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -i $instances -a cmd1
		commander::runcmd -c macs -v -b -i $instances2 -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
	fi

	return 0
}

peaks::macs_idr(){
	declare -a tdirs
	_cleanup::peaks::macs_idr(){
		rm -rf "${tdirs[@]}"
	}

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-f <size>       | assumed mean fragment
			-g <genome>     | path to
			-q <ripseq>     | true/false
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to
			-p <tmpdir>     | path to
			-y <pointy>     | true/false call only pointy peaks
			-z <strict>     | true/false peak filters

			requires prior execution of alignment::mkreplicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false ripseq=false genome threads memory maxmemory fragmentsize outdir tmpdir strict=false pointy=false
	declare -n _mapper_macs _nidx_macs _nridx_macs _tidx_macs _ridx_macs _pidx_macs
	while getopts 'S:s:t:m:M:g:f:q:r:a:b:i:j:k:o:p:z:y:' arg; do
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
			a)	((++mandatory)); _nidx_macs=$OPTARG;;
			b)	_nridx_macs=$OPTARG;;
			i)	((++mandatory)); _tidx_macs=$OPTARG;;
			j)	((++mandatory)); _ridx_macs=$OPTARG;;
			k)	((++mandatory)); _pidx_macs=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			y)	pointy=$OPTARG;;
			z)	strict=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 11 ]] && _usage

	commander::printinfo "peak calling macs"

	declare -n _bams_macs=${_mapper_macs[0]}
	local genomesize nf="${_bams_macs[${_nidx_macs[0]}]}"
	local x=$(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 256)
	if [[ $x -gt 0 ]]; then
		genomesize=$(faCount "$genome" | tail -1 | awk '{print $3+$4+$5+$6}')
	else
		genomesize=$(unique-kmers.py -k 100 "$genome" 2>&1 | tail -1 | awk '{print $NF}')
	fi

	local params='' mult=5
	local numchr=$(samtools view -H "$nf" | grep -c '^@SQ')
	local buffer=$(( (1024*1024*memory)/(numchr*mult) ))

	local instances instances2 ithreads ithreads2
	instances=$((${#_bams_macs[@]} * ${#_mapper_macs[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -T $threads -i $instances)
	if [[ $buffer -lt 100000 ]]; then
		params+=" --buffer-size $buffer"
		read -r instances2 ithreads2 < <(configure::instances_by_memory -T $threads -m $memory -M "$maxmemory")
	else
		read -r instances2 ithreads2 < <(configure::instances_by_memory -T $threads -m $(( (numchr*buffer*mult)/1024/1024 )) -M "$maxmemory")
	fi

	if $strict; then
		params+=' -q 0.05'
	else
		params+=' -p 0.01'
	fi

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a cmd1 cmd2 cmd3 cmd4 cmd5 toidr
	for m in "${_mapper_macs[@]}"; do
		declare -n _bams_macs=$m
		odir="$outdir/$m/macs"

		for i in "${!_nidx_macs[@]}"; do
			nf="${_bams_macs[${_nidx_macs[$i]}]}"
			tf="${_bams_macs[${_tidx_macs[$i]}]}"
			rf="${_bams_macs[${_ridx_macs[$i]}]}"
			pf="${_bams_macs[${_pidx_macs[$i]}]}"


			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.macs)")
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bedtools bamtobed -split -i "$nf"
			CMD
				pigz -p $ithreads -k -c > "${tdirs[-1]}/$(basename "$nf").bed.gz"
			CMD
			nf="${tdirs[-1]}/$(basename "$nf").bed.gz"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.macs)")
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bedtools bamtobed -split -i "$f"
				CMD
					pigz -p $ithreads -k -c > "${tdirs[-1]}/$(basename "$f").bed.gz"
				CMD
				f="${tdirs[-1]}/$(basename "$f").bed.gz"

				o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
				mkdir -p "$odir/$o"

				printf '' > "$odir/$o/$o.model_peaks.narrowPeak"
				if $ripseq; then
					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
						macs2 callpeak
						-f BED
						-t "$f"
						-c "$nf"
						-g $genomesize
						--outdir "$odir/$o"
						-n "$o.nomodel"
						--tempdir "${tdirs[-1]}"
						--call-summits
						-B
						--SPMR
						--keep-dup all
						--verbose 2
						--nomodel
						--shift $($pointy && echo "-$((fragmentsize/8))" || echo "-$((fragmentsize/2))")
						--extsize $($pointy && echo "$((fragmentsize/2))" || echo $fragmentsize)
						$params
					CMD
				else
					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
						macs2 callpeak
						-f BED
						-t "$f"
						-c "$nf"
						-g $genomesize
						--outdir "$odir/$o"
						-n "$o.model"
						--tempdir "${tdirs[-1]}"
						--call-summits
						-B
						--SPMR
						--keep-dup all
						--verbose 2
						$params
						$($pointy && echo "-m 5 150 --bw 150")
					CMD

					tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.macs)")
					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
						macs2 callpeak
						-f BED
						-t "$f"
						-c "$nf"
						-g $genomesize
						--outdir "$odir/$o"
						-n "$o.nomodel"
						--tempdir "${tdirs[-1]}"
						--call-summits
						-B
						--SPMR
						--keep-dup all
						--verbose 2
						--nomodel
						--shift 0
						--extsize $($pointy && echo "$((fragmentsize/4))" || echo $fragmentsize)
						$params
					CMD
				fi

				commander::makecmd -a cmd3 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD'
					sort -k1,1 -k2,2n -k3,3n "$odir/$o/$o.nomodel_peaks.narrowPeak" "$odir/$o/$o.model_peaks.narrowPeak"
				CMD
					awk '$2>=0 && $3>=0'
				CMD
					bedtools merge -c 4,5,6,7,8,9,10 -o distinct,max,distinct,collapse,max,max,collapse
				CMD
					perl -F'\t' -lane '
						$i=0;
						@enrich=split/,/,$F[-4];
						@summits=split/,/,$F[-1];
						for (0..$#enrich){
							$i=$_ if $enrich[$_] < $enrich[$i]
						}
						$F[-4]=$enrich[$i];
						$F[-1]=$summits[$i];
						$F[3]="peak_$.";
						print join"\t",@F;
					'
				CMD

				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd4 \
				-2 cmd5 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
		for i in "${!_nridx_macs[@]}"; do
			# if a normal replicate is given, then run idr on: n-pp (1/9) + nr-pp (3/9) vs nfp-fp (11/12)
			#                                           1  2   3   4  5  6  7  8  9   10   11  12   13  14
			# m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2] -> m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2 NFP1 FP1 NFP2 FP2]
			# n   1 2   3 4   11 13     # n   1  2  3  4     5  6  7  8    21 23 25 27
			# nr  3 4                   # nr  5  6  7  8
			# t   5 6   5 6   5  6      # t   9 10 11 12     9 10 11 12     9 10 11 12
			# r   7 8   7 8   7  8      # r  13 14 15 16    13 14 15 16    13 14 15 16
			# p   9 10  9 10  12 14     # p  17 18 19 20    17 18 19 20    22 24 26 28
			nf="${_bams_macs[${_nidx_macs[$i]}]}" # 1
			nrf="${_bams_macs[${_nridx_macs[$i]}]}" # 3
			pf="${_bams_macs[${_pidx_macs[$i]}]}" # 9

			x=$(( ${_nridx_macs[$i]} + ${_pidx_macs[$i]} )) # 12
			pff="${_bams_macs[$x]}"
			nff="${_bams_macs[$((--x))]}" # 11

			toidr=( "$odir/$(echo -e "$(basename "$nf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nrf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nff")\t$(basename "$pff")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd4 \
				-2 cmd5 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
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
		commander::runcmd -v -b -i $threads -a cmd3
		commander::runcmd -c idr -v -b -i $threads -a cmd4
		commander::runcmd -v -b -i $threads -a cmd5
	fi

	return 0
}

peaks::gem(){
	local tdir
	_cleanup::peaks::gem(){
		rm -rf "$tdir"
	}

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-g <genome>     | path to
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-p <tmpdir>     | path to
			-q <ripseq>     | true/false
			-x <strandness> | hash per bam of (if ripseq)
			-y <pointy>     | true/false call only pointy peaks
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false ripseq=false threads memory maxmemory genome outdir tmpdir strict=false pointy=false
	declare -n _mapper_gem _strandness_gem _nidx_gem _tidx_gem
	while getopts 'S:s:t:m:M:g:q:r:x:a:i:o:p:z:y:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			r)	((++mandatory)); _mapper_gem=$OPTARG;;
			a)	_nidx_gem=$OPTARG;;
			i)	_tidx_gem=$OPTARG;;
			x)	_strandness_gem=$OPTARG;;
			q)	ripseq=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			z)	strict=$OPTARG;;
			y)	pointy=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage

	$ripseq && [[ ${#_strandness_gem[@]} -eq 0 ]] && _usage
	[[ $_nidx_gem ]] && [[ ! $_tidx_gem ]] && _usage

	declare -n _bams_gem=${_mapper_gem[0]}
	if [[ ! $_nidx_gem && ! $_tidx_gem ]]; then
		declare -a tidx_gem=("${!_bams_gem[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_gem=tidx_gem
	fi

	commander::printinfo "peak calling gem"

	local minstances mthreads jmem jgct jcgct instances=$(( ${#_tidx_gem[@]} * ${#_mapper_gem[@]} ))
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i $instances -T $threads -m $memory -M "$maxmemory")

	# get effective genome size
	# if multimapped reads: genome minus Ns , else genome minus Ns minus repetetive Elements
	local genomesize
	local x=$(samtools view -F 4 "${_bams_gem[0]}" | head -10000 | cat <(samtools view -H "${_bams_gem[0]}") - | samtools view -c -f 256)
	if [[ $x -gt 0 ]]; then
		genomesize=$(faCount $genome | tail -1 | awk '{print $3+$4+$5+$6}')
	else
		genomesize=$(unique-kmers.py -k 100 $genome 2>&1 | tail -1 | awk '{print $NF}')
	fi
	tdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.genome)"
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
		if $pointy; then
			params+=' --smooth 5'
		fi
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

	local m i f o odir nf
	declare -a cmd1 cmd2
	for m in "${_mapper_gem[@]}"; do
		declare -n _bams_gem=$m
		odir="$outdir/$m/gem"

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
			commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- 'CMD' {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- 'CMD'
				paste "$odir/$o/$o.GPS_events.narrowPeak" <(tail -n +2 "$odir/$o/$o.GPS_events.txt")
			CMD
				perl -lane 'if($F[13] eq "NaN"){$F[6]=$F[11]}else{$F[6]=$F[13]}; print join"\t",@F[0..9]'
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
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -v -b -i $minstances -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

peaks::gem_idr(){
	local tdir
	_cleanup::peaks::gem_idr(){
		rm -rf "$tdir"
	}

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-g <genome>     | path to
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to
			-p <tmpdir>     | path to
			-q <ripseq>     | true/false
			-x <strandness> | hash per bam of (if ripseq)
			-y <pointy>     | true/false call only pointy peaks
			-z <strict>     | true/false peak filters

			requires prior execution of alignment::mkreplicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false ripseq=false threads memory maxmemory genome outdir tmpdir strict=false pointy=false
	declare -n _mapper_gem _strandness_gem _nidx_gem _nridx_gem _tidx_gem _ridx_gem _pidx_gem
	while getopts 'S:s:t:m:M:g:q:r:x:a:b:i:j:k:o:p:z:y:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			r)	((++mandatory)); _mapper_gem=$OPTARG;;
			x)	_strandness_gem=$OPTARG;;
			q)	ripseq=$OPTARG;;
			a)	((++mandatory)); _nidx_gem=$OPTARG;; # contains nridx by alignment::mkreplicates
			b)	_nridx_gem=$OPTARG;;
			i)	((++mandatory)); _tidx_gem=$OPTARG;;
			j)	((++mandatory)); _ridx_gem=$OPTARG;;
			k)	((++mandatory)); _pidx_gem=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			z)	$OPTARG && strict=true;;
			y)	pointy=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 10 ]] && _usage
	$ripseq && [[ ${#_strandness_gem[@]} -eq 0 ]] && _usage

	commander::printinfo "peak calling gem"

	local minstances mthreads jmem jgct jcgct instances=$(( (${#_tidx_gem[@]} + ${#_ridx_gem[@]} + ${#_pidx_gem[@]}) * ${#_mapper_gem[@]} ))
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i $instances -T $threads -m $memory -M "$maxmemory")

	declare -n _bams_gem=${_mapper_gem[0]}
	local nf="${_bams_gem[${_nidx_gem[0]}]}"

	local genomesize
	local x=$(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 256)
	if [[ $x -gt 0 ]]; then
		genomesize=$(faCount $genome | tail -1 | awk '{print $3+$4+$5+$6}')
	else
		genomesize=$(unique-kmers.py -k 100 $genome 2>&1 | tail -1 | awk '{print $NF}')
	fi
	tdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.genome)"
	samtools view -H $nf | sed -rn '/^@SQ/{s/.+\tSN:(\S+)\s+LN:(\S+).*/\1\t\2/p}' > "$tdir/chr.info"
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
		local x=$(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 1)
		[[ $x -eq 0 ]] && strandness=${_strandness_gem["$nf"]}
		if [[ $strandness -ne 0 ]]; then
			params+=' --strand_type 1'
			params2+=' -s' # for bedtools merge
		fi
		params+=' --nd 2'

		params+=" --d '$(dirname "$(which gem)")/Read_Distribution_CLIP.txt'"
		if $pointy; then
			params+=' --smooth 5'
		fi
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
	else
		params+=' -q 0 --relax'
	fi

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a cmd1 cmd2 cmd3 cmd4 toidr
	for m in "${_mapper_gem[@]}"; do
		declare -n _bams_gem=$m
		odir="$outdir/$m/gem"

		for i in "${!_nidx_gem[@]}"; do
			nf="${_bams_gem[${_nidx_gem[$i]}]}"
			tf="${_bams_gem[${_tidx_gem[$i]}]}"
			rf="${_bams_gem[${_ridx_gem[$i]}]}"
			pf="${_bams_gem[${_pidx_gem[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o=$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')

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
						--ctrl "$nf"
						--f SAM
						--nrf
						--sl
						--s $genomesize
						--outNP
						--fold 1
						--nf
						$params
				CMD

				commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- 'CMD' {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- 'CMD'
					paste "$odir/$o/$o.GPS_events.narrowPeak" <(tail -n +2 "$odir/$o/$o.GPS_events.txt")
				CMD
					perl -lane '$F[6]=$F[13]; print join"\t",@F[0..9]'
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

				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
		for i in "${!_nridx_gem[@]}"; do
			# if a normal replicate is given, then run idr on: n-pp (1/9) + nr-pp (3/9) vs nfp-fp (11/12)
			#                                           1  2   3   4  5  6  7  8  9   10   11  12   13  14
			# m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2] -> m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2 NFP1 FP1 NFP2 FP2]
			# n   1 2   3 4   11 13     # n   1  2  3  4     5  6  7  8    21 23 25 27
			# nr  3 4                   # nr  5  6  7  8
			# t   5 6   5 6   5  6      # t   9 10 11 12     9 10 11 12     9 10 11 12
			# r   7 8   7 8   7  8      # r  13 14 15 16    13 14 15 16    13 14 15 16
			# p   9 10  9 10  12 14     # p  17 18 19 20    17 18 19 20    22 24 26 28
			nf="${_bams_gem[${_nidx_gem[$i]}]}" # 1
			nrf="${_bams_gem[${_nridx_gem[$i]}]}" # 3
			pf="${_bams_gem[${_pidx_gem[$i]}]}" # 9

			x=$(( ${_nridx_gem[$i]} + ${_pidx_gem[$i]} )) # 12
			pff="${_bams_gem[$x]}"
			nff="${_bams_gem[$((x-1))]}" # 11

			toidr=( "$odir/$(echo -e "$(basename "$nf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nrf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nff")\t$(basename "$pff")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		commander::runcmd -v -b -i $minstances -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -c idr -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}

peaks::matk(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-M <maxmemory>  | amount of
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-i <tidx>       | array of RIP* bam idices within -r
			-o <outdir>     | path to
			-p <tmpdir>     | path to

			requires prior execution of alignment::mkreplicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory outdir tmpdir
	declare -n _mapper_matk _nidx_matk _nridx_matk _tidx_matk _ridx_matk _pidx_matk
	while getopts 'S:s:t:M:r:x:a:b:i:j:k:o:p:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			r)	((++mandatory)); _mapper_matk=$OPTARG;;
			a)	((++mandatory)); _nidx_matk=$OPTARG;; # contains nridx by alignment::mkreplicates
			b)	_nridx_matk=$OPTARG;;
			i)	((++mandatory)); _tidx_matk=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage

	declare -n _bams_matk=${_mapper_matk[0]}
	if [[ ! $_nidx_matk && ! $_tidx_matk ]]; then
		declare -a tidx_matk=("${!_bams_matk[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_matk=tidx_matk
	fi

	commander::printinfo "peak calling matk"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i 1 -T $threads -M "$maxmemory")

	local m i f o odir nf
	declare -a cmd1 cmd2
	for m in "${_mapper_matk[@]}"; do
		declare -n _bams_matk=$m
		odir="$outdir/$m/matk"

		for i in "${!_tidx_matk[@]}"; do
			nf="${_bams_matk[${_nidx_matk[$i]}]}"
			f="${_bams_matk[${_tidx_matk[$i]}]}"
			o="$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"

			mkdir -p "$odir/$o"

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
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

			commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- 'CMD'
				samtools bedcov -Q 0 -g 1796 -j "$odir/$o/$o.bed" "$nf" "$f"
			CMD
				awk -v c=$(samtools idxstats "$nf" | datamash sum 3) -v i=$(samtools idxstats "$f" | datamash sum 3) '\$NF>0 && \$NF>\$(NF-1){print \$0,(\$NF/i)/((\$(NF-1)+1)/c)}'
			CMD
				perl -lane 'print join"\t",(@F[0..3],0,".",$F[-1],-1,-1*log($F[4])/log(10),-1)'
			CMD
				sort -k1,1 -k2,2n -k3,3n
			CMD
				bedtools merge -c 4,5,6,7,8,9,10 -o distinct,max,distinct,max,max,max,distinct
			CMD
				awk -v OFS='\t' '{$4="peak_"NR; print}'
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

peaks::matk_idr(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-M <maxmemory>  | amount of
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to
			-p <tmpdir>     | path to

			requires prior execution of alignment::mkreplicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory outdir tmpdir
	declare -n _mapper_matk _nidx_matk _nridx_matk _tidx_matk _ridx_matk _pidx_matk
	while getopts 'S:s:t:M:r:x:a:b:i:j:k:o:p:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			r)	((++mandatory)); _mapper_matk=$OPTARG;;
			a)	((++mandatory)); _nidx_matk=$OPTARG;; # contains nridx by alignment::mkreplicates
			b)	_nridx_matk=$OPTARG;;
			i)	((++mandatory)); _tidx_matk=$OPTARG;;
			j)	((++mandatory)); _ridx_matk=$OPTARG;;
			k)	((++mandatory)); _pidx_matk=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 8 ]] && _usage

	commander::printinfo "peak calling matk"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i 1 -T $threads -M "$maxmemory")

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a cmd1 cmd2 cmd3 cmd4 toidr
	for m in "${_mapper_matk[@]}"; do
		declare -n _bams_matk=$m
		odir="$outdir/$m/matk"

		for i in "${!_nidx_matk[@]}"; do
			nf="${_bams_matk[${_nidx_matk[$i]}]}"
			tf="${_bams_matk[${_tidx_matk[$i]}]}"
			rf="${_bams_matk[${_ridx_matk[$i]}]}"
			pf="${_bams_matk[${_pidx_matk[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o=$(echo -e "$(basename "$nf")\t$(basename "$f")" | sed -E 's/(\..+)\t(.+)\1/-\2/')

				mkdir -p "$odir/$o"

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
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

				commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- 'CMD'
					samtools bedcov -Q 0 -g 1796 -j "$odir/$o/$o.bed" "$nf" "$f"
				CMD
					awk -v c=$(samtools idxstats "$nf" | datamash sum 3) -v i=$(samtools idxstats "$f" | datamash sum 3) '\$NF>0 && \$NF>\$(NF-1){print \$0,(\$NF/i)/((\$(NF-1)+1)/c)}'
				CMD
					perl -lane 'print join"\t",(@F[0..3],0,".",$F[-1],-1,-1*log($F[4])/log(10),-1)'
				CMD
					sort -k1,1 -k2,2n -k3,3n
				CMD
					bedtools merge -c 4,5,6,7,8,9,10 -o distinct,max,distinct,max,max,max,distinct
				CMD
					awk -v OFS='\t' '{$4="peak_"NR; print}'
				CMD

				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
		for i in "${!_nridx_matk[@]}"; do
			# if a normal replicate is given, then run idr on: n-pp (1/9) + nr-pp (3/9) vs nfp-fp (11/12)
			#                                           1  2   3   4  5  6  7  8  9   10   11  12   13  14
			# m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2] -> m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2 NFP1 FP1 NFP2 FP2]
			# n   1 2   3 4   11 13     # n   1  2  3  4     5  6  7  8    21 23 25 27
			# nr  3 4                   # nr  5  6  7  8
			# t   5 6   5 6   5  6      # t   9 10 11 12     9 10 11 12     9 10 11 12
			# r   7 8   7 8   7  8      # r  13 14 15 16    13 14 15 16    13 14 15 16
			# p   9 10  9 10  12 14     # p  17 18 19 20    17 18 19 20    22 24 26 28
			nf="${_bams_matk[${_nidx_matk[$i]}]}" # 1
			nrf="${_bams_matk[${_nridx_matk[$i]}]}" # 3
			pf="${_bams_matk[${_pidx_matk[$i]}]}" # 9

			x=$(( ${_nridx_matk[$i]} + ${_pidx_matk[$i]} )) # 12
			pff="${_bams_matk[$x]}"
			nff="${_bams_matk[$((x-1))]}" # 11

			toidr=( "$odir/$(echo -e "$(basename "$nf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nrf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nff")\t$(basename "$pff")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		commander::runcmd -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -c idr -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}

peaks::peakachu() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-f <size>       | assumed mean fragment
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r (if paired input, requires also -i)
			-i <tidx>       | array of IP* bam idices within -r (if paired input, requires also -a)
			-o <outdir>     | path to
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false genome threads fragmentsize outdir strict=false
	declare -n _mapper_peakachu _nidx_peakachu _tidx_peakachu
	while getopts 'S:s:t:f:r:a:i:o:z:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_peakachu=$OPTARG;;
			a)	_nidx_peakachu=$OPTARG;;
			i)	_tidx_peakachu=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			z)	strict=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage
	[[ $_nidx_peakachu ]] && [[ ! $_tidx_peakachu ]] && _usage

	declare -n _bams_peakachu=${_mapper_peakachu[0]}
	if [[ ! $_nidx_peakachu && ! $_tidx_peakachu ]]; then
		declare -a tidx_peakachu=("${!_bams_peakachu[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_peakachu=tidx_peakachu
	fi

	commander::printinfo "peak calling peakachu"

	local x params
	# infer SE or PE
	x=$(samtools view -F 4 "${_bams_peakachu[0]}" | head -10000 | cat <(samtools view -H "${_bams_peakachu[0]}") - | samtools view -c -f 1)
	[[ $x -gt 0 ]] && params='--paired_end'

	local m f i nf o odir
	declare -a cmd1 cmd2
	for m in "${_mapper_peakachu[@]}"; do
		declare -n _bams_peakachu=$m
		odir="$outdir/$m/peakachu"

		f='--exp_libs'
		for i in "${!_tidx_peakachu[@]}"; do
			f+="$(printf ' "%s"' "${_bams_peakachu[${_tidx_peakachu[$i]}]}")"
		done
		if [[ $_nidx_peakachu ]]; then
			nf="--ctr_libs"
			for i in "${!_nidx_peakachu[@]}"; do
				nf+="$(printf ' "%s"' "${_bams_peakachu[${_nidx_peakachu[$i]}]}")"
			done
			o="$(echo -e "$(basename "${_bams_peakachu[${_nidx_peakachu[0]}]}")\t$(basename "${_bams_peakachu[${_tidx_peakachu[0]}]}")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
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
				$f
				$nf
				--max_insert_size $((2*fragmentsize))
				--max_proc $threads
				--norm_method deseq
				--mad_multiplier 0.0
				--fc_cutoff 2
				--padj_threshold $($strict && echo 0.05 || echo 0.1)
				--output_folder "$odir/$o"
		CMD

		commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- 'CMD'
			awk '!/^replicon/' "$odir/$o/peak_tables/"*.csv
		CMD
			awk -v OFS='\t' '$(NF-9)!="inf"{print $1,$3,$4,".",0,".",$(NF-9),-1,-1,-1}'
		CMD
			sort -k1,1 -k2,2n -k3,3n
		CMD
			bedtools merge -c 4,5,6,7,8,9,10 -o distinct,distinct,distinct,max,distinct,distinct,distinct
		CMD
			awk -v OFS='\t' '{$4="peak_"NR; print}'
		CMD
		# NF-9 is either base_mean or if ctr_libs fold_change
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c peakachu -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

peaks::peakachu_idr() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-f <size>       | assumed mean fragment
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false genome threads fragmentsize outdir strict=false
	declare -n _mapper_peakachu _nidx_peakachu _nridx_peakachu _tidx_peakachu _ridx_peakachu _pidx_peakachu
	while getopts 'S:s:t:f:r:a:b:i:j:k:o:z:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_peakachu=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			a)	((++mandatory)); _nidx_peakachu=$OPTARG;;
			b)	_nridx_peakachu=$OPTARG;;
			i)	((++mandatory)); _tidx_peakachu=$OPTARG;;
			j)	((++mandatory)); _ridx_peakachu=$OPTARG;;
			k)	((++mandatory)); _pidx_peakachu=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			z)	strict=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 8 ]] && _usage

	commander::printinfo "peak calling peakachu"

	declare -n _bams_peakachu=${_mapper_peakachu[0]}
	local x params nf="${_bams_peakachu[${_nidx_peakachu[0]}]}"
	# infer SE or PE
	x=$(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 1)
	[[ $x -gt 0 ]] && params='--paired_end'

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a cmd1 cmd2 cmd3 cmd4 toidr
	for m in "${_mapper_peakachu[@]}"; do
		declare -n _bams_peakachu=$m
		odir="$outdir/$m/peakachu"

		for i in "${!_nidx_peakachu[@]}"; do
			nf="${_bams_peakachu[${_nidx_peakachu[$i]}]}"
			tf="${_bams_peakachu[${_tidx_peakachu[$i]}]}"
			rf="${_bams_peakachu[${_ridx_peakachu[$i]}]}"
			pf="${_bams_peakachu[${_pidx_peakachu[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o=$(echo -e "$(basename $nf)\t$(basename $f)" | sed -E 's/(\..+)\t(.+)\1/-\2/')

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					rm -rf "$odir/$o"
				CMD
					peakachu adaptive
						$params
						--exp_libs "$f"
						--ctr_libs "$nf"
						--max_insert_size $((2*fragmentsize))
						--max_proc $threads
						--norm_method deseq
						--mad_multiplier 0.0
						--fc_cutoff 2
						--padj_threshold $($strict && echo 0.05 || echo 0.1)
						--output_folder "$odir/$o"
				CMD

				commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- 'CMD'
					awk '!/^replicon/' "$odir/$o/peak_tables/"*.csv
				CMD
					awk -v OFS='\t' '$(NF-9)!="inf"{print $1,$3,$4,".",0,".",$(NF-9),-1,-1,-1}'
				CMD
					sort -k1,1 -k2,2n -k3,3n
				CMD
					bedtools merge -c 4,5,6,7,8,9,10 -o distinct,distinct,distinct,max,distinct,distinct,distinct
				CMD
					awk -v OFS='\t' '{$4="peak_"NR; print}'
				CMD
				# NF-9 is either base_mean or if ctr_libs fold_change

				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
		for i in "${!_nridx_peakachu[@]}"; do
			# if a normal replicate is given, then run idr on: n-pp (1/9) + nr-pp (3/9) vs nfp-fp (11/12)
			#                                           1  2   3   4  5  6  7  8  9   10   11  12   13  14
			# m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2] -> m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2 NFP1 FP1 NFP2 FP2]
			# n   1 2   3 4   11 13     # n   1  2  3  4     5  6  7  8    21 23 25 27
			# nr  3 4                   # nr  5  6  7  8
			# t   5 6   5 6   5  6      # t   9 10 11 12     9 10 11 12     9 10 11 12
			# r   7 8   7 8   7  8      # r  13 14 15 16    13 14 15 16    13 14 15 16
			# p   9 10  9 10  12 14     # p  17 18 19 20    17 18 19 20    22 24 26 28
			nf="${_bams_peakachu[${_nidx_peakachu[$i]}]}" # 1
			nrf="${_bams_peakachu[${_nridx_peakachu[$i]}]}" # 3
			pf="${_bams_peakachu[${_pidx_peakachu[$i]}]}" # 9

			x=$(( ${_nridx_peakachu[$i]} + ${_pidx_peakachu[$i]} )) # 12
			pff="${_bams_peakachu[$x]}"
			nff="${_bams_peakachu[$((--x))]}" # 11

			toidr=( "$odir/$(echo -e "$(basename "$nf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nrf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nff")\t$(basename "$pff")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		commander::runcmd -c peakachu -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -c idr -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}

peaks::m6aviewer() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
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
			-p <tmpdir>     | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false genome threads memory maxmemory fragmentsize outdir tmpdir
	declare -n _mapper_m6aviewer _nidx_m6aviewer _tidx_m6aviewer
	while getopts 'S:s:t:m:M:f:r:a:b:i:j:k:o:p:' arg; do
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
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage
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
				o="${o%.*}"
			fi

			commander::makecmd -a cmd2 -s '|' -o "$odir/$o" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- 'CMD'
				cat "$odir/results_$(basename "$f").txt"
			CMD
				awk -v OFS='\t' '!/^\s*$/ && NR>2{ if($5==0){logp=999}else{logp=log($5)/log(10)*-1}; print $3,$4-75,$4+76,".",0,".",$6,logp,-1,75}'
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
						$i=$_ if $enrich[$_] < $enrich[$i]
					}
					$F[-4]=$enrich[$i];
					$F[-1]=$summits[$i];
					$F[3]="peak_$.";
					print join"\t",@F;
				'
			CMD
		done

		#cmd1+=("m6aviewer -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.TMPDIR='$tmpdir'")
		m6aviewer -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.TMPDIR='$tmpdir' &> /dev/null
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

peaks::m6aviewer_idr() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-f <size>       | assumed mean fragment
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to
			-p <tmpdir>     | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false genome threads memory maxmemory fragmentsize outdir tmpdir
	declare -n _mapper_m6aviewer _nidx_m6aviewer _nridx_m6aviewer _tidx_m6aviewer _ridx_m6aviewer _pidx_m6aviewer
	while getopts 'S:s:t:m:M:f:r:a:b:i:j:k:o:p:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			r)	((++mandatory)); _mapper_m6aviewer=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			a)	((++mandatory)); _nidx_m6aviewer=$OPTARG;;
			b)	_nridx_m6aviewer=$OPTARG;;
			i)	((++mandatory)); _tidx_m6aviewer=$OPTARG;;
			j)	((++mandatory)); _ridx_m6aviewer=$OPTARG;;
			k)	((++mandatory)); _pidx_m6aviewer=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 10 ]] && _usage

	commander::printinfo "peak calling m6aviewer"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i 1 -T $threads -m $memory -M "$maxmemory")

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a cmd1 cmd2 cmd3 cmd4 toidr
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

		for i in "${!_nidx_m6aviewer[@]}"; do
			nf="${_bams_m6aviewer[${_nidx_m6aviewer[$i]}]}"
			tf="${_bams_m6aviewer[${_tidx_m6aviewer[$i]}]}"
			rf="${_bams_m6aviewer[${_ridx_m6aviewer[$i]}]}"
			pf="${_bams_m6aviewer[${_pidx_m6aviewer[$i]}]}"

			commander::printinfo {COMMANDER[0]}<<- EOF
				load
				$nf
				versus
				$tf
				$rf
				$pf
			EOF

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o=$(echo -e "$(basename $nf)\t$(basename $f)" | sed -E 's/(\..+)\t(.+)\1/-\2/').narrowPeak

				commander::makecmd -a cmd2 -s '|' -o "$odir/$o" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- 'CMD'
					cat "$odir/results_$(basename "$f").txt"
				CMD
					awk -v OFS='\t' '!/^\s*$/ && NR>2{ if($5==0){logp=999}else{logp=log($5)/log(10)*-1}; print $3,$4-75,$4+76,".",0,".",$6,logp,-1,75}'
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
							$i=$_ if $enrich[$_] < $enrich[$i]
						}
						$F[-4]=$enrich[$i];
						$F[-1]=$summits[$i];
						$F[3]="peak_$.";
						print join"\t",@F;
					'
				CMD

				toidr+=("$odir/$o")
			done

			peaks::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
		for i in "${!_nridx_m6aviewer[@]}"; do
			# if a normal replicate is given, then run idr on: n-pp (1/9) + nr-pp (3/9) vs nfp-fp (11/12)
			#                                           1  2   3   4  5  6  7  8  9   10   11  12   13  14
			# m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2] -> m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2 NFP1 FP1 NFP2 FP2]
			# n   1 2   3 4   11 13     # n   1  2  3  4     5  6  7  8    21 23 25 27
			# nr  3 4                   # nr  5  6  7  8
			# t   5 6   5 6   5  6      # t   9 10 11 12     9 10 11 12     9 10 11 12
			# r   7 8   7 8   7  8      # r  13 14 15 16    13 14 15 16    13 14 15 16
			# p   9 10  9 10  12 14     # p  17 18 19 20    17 18 19 20    22 24 26 28
			nf="${_bams_m6aviewer[${_nidx_m6aviewer[$i]}]}" # 1
			nrf="${_bams_m6aviewer[${_nridx_m6aviewer[$i]}]}" # 3
			pf="${_bams_m6aviewer[${_pidx_m6aviewer[$i]}]}" # 9

			x=$(( ${_nridx_m6aviewer[$i]} + ${_pidx_m6aviewer[$i]} )) # 12
			pff="${_bams_m6aviewer[$x]}"
			nff="${_bams_m6aviewer[$((--x))]}" # 11

			toidr=( "$odir/$(echo -e "$(basename "$nf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nrf")\t$(basename "$pf")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "$nff")\t$(basename "$pff")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done

		#cmd1+=("m6aviewer -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.TMPDIR='$tmpdir'")
		m6aviewer -Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.TMPDIR='$tmpdir' &> /dev/null
	done

	if $skip; then
		# commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		# commander::runcmd -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -c idr -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}
