#! /usr/bin/env bash
# (c) Konstantin Riege

# https://sites.google.com/site/anshulkundaje/projects/idr
# https://docs.google.com/document/d/1lG_Rd7fnYgRpSIqrIfuVlAz2dW1VaSQThzk836Db99c/edit
# https://sites.google.com/site/anshulkundaje/projects/idr/deprecated

function peaks::_idr(){
	function _usage(){
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

	local OPTIND arg mandatory t r p o y="narrowPeak"
	declare -n _cmds1_idr _cmds2_idr
	while getopts '1:2:t:r:p:o:y:' arg; do
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
	[[ "$y" == "narrowPeak" || "$y" == "bed" ]] || _usage

	# --rank signal.value (default) p.value q.value score (bed default) <columnidx>
	# bed does not work! bed_loader expects nine/ten columns i.e. narrowPeak format, but with score in column 5
	# no need for sorted input!
	commander::makecmd -a _cmds1_idr -s ';' -c {COMMANDER[0]}<<- CMD
		idr
		--samples "$t" "$r"
		--peak-list "$p"
		--input-file-type $y
		--rank signal.value
		--output-file "$o"
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

function peaks::macs_idr(){
	function _usage(){
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
			-w <broad>      | true/false peak detection
			-y <pointy>     | true/false call only pointy peaks
			-z <strict>     | true/false peak filters

			requires prior execution of alignment::mkreplicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false ripseq=false genome threads memory maxmemory fragmentsize outdir strict=false pointy=false broad=false
	declare -n _mapper_macs _nidx_macs _nridx_macs _tidx_macs _ridx_macs _pidx_macs
	while getopts 'S:s:t:m:M:g:f:q:r:a:b:i:j:k:o:w:z:y:' arg; do
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
			w)	broad=$OPTARG;;
			y)	pointy=$OPTARG;;
			z)	strict=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 10 ]] && _usage

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a nidx_macs tidx_macs cmd1 cmd2 toidr

	for m in "${_mapper_macs[@]}"; do
		declare -n _bams_macs=$m
		odir="$outdir/$m/macs"
		nidx_macs=()
		tidx_macs=()

		for i in "${!_nidx_macs[@]}"; do
			nidx_macs+=(${_nidx_macs[$i]} ${_nidx_macs[$i]} ${_nidx_macs[$i]})
			tidx_macs+=(${_tidx_macs[$i]} ${_ridx_macs[$i]} ${_pidx_macs[$i]})

			nf="${_bams_macs[${_nidx_macs[$i]}]}"
			tf="${_bams_macs[${_tidx_macs[$i]}]}"
			rf="${_bams_macs[${_ridx_macs[$i]}]}"
			pf="${_bams_macs[${_pidx_macs[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o="$(echo -e "$(basename "${nf%.*}")\t$(basename "${f%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
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

			toidr=( "$odir/$(echo -e "$(basename "${nf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nrf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nff%.*}")\t$(basename "${pff%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	peaks::macs \
		-S false \
		-s $skip \
		-q $ripseq \
		-f $fragmentsize \
		-g "$genome" \
		-a nidx_macs \
		-i tidx_macs \
		-r ${!_mapper_macs} \
		-t $threads \
		-m $memory \
		-M "$maxmemory" \
		-o "$outdir" \
		-w $broad \
		-y $pointy \
		-z $strict

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c idr -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function peaks::seacr_idr(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-m <memory>     | amount of. required if -c true
			-M <maxmemory>  | amount of. if -c true
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-f <size>       | assumed mean fragment. required if -c true
			-g <genome>     | path to. required if -c true
			-t <threads>    | number of
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to
			-z <strict>     | true/false peak filters

			requires prior execution of alignment::mkreplicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir strict=false macsdir memory maxmemory fragmentsize genome
	declare -n _mapper_seacr _nidx_seacr _nridx_seacr _tidx_seacr _ridx_seacr _pidx_seacr
	while getopts 'S:s:c:m:M:f:g:t:r:a:b:i:j:k:o:z:' arg; do
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
			a)	((++mandatory)); _nidx_seacr=$OPTARG;;
			b)	_nridx_seacr=$OPTARG;;
			i)	((++mandatory)); _tidx_seacr=$OPTARG;;
			j)	((++mandatory)); _ridx_seacr=$OPTARG;;
			k)	((++mandatory)); _pidx_seacr=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			z)	strict=$OPTARG;;
			*) _usage;;
		esac
	done
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 7 ]] && _usage
	else
		[[ $mandatory -lt 10 ]] && _usage
	fi

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a nidx_seacr tidx_seacr cmd1 cmd2 toidr

	for m in "${_mapper_seacr[@]}"; do
		declare -n _bams_seacr=$m
		odir="$outdir/$m/seacr"
		nidx_seacr=()
		tidx_seacr=()

		for i in "${!_nidx_seacr[@]}"; do
			nidx_seacr+=(${_nidx_seacr[$i]} ${_nidx_seacr[$i]} ${_nidx_seacr[$i]})
			tidx_seacr+=(${_tidx_seacr[$i]} ${_ridx_seacr[$i]} ${_pidx_seacr[$i]})

			nf="${_bams_seacr[${_nidx_seacr[$i]}]}"
			tf="${_bams_seacr[${_tidx_seacr[$i]}]}"
			rf="${_bams_seacr[${_ridx_seacr[$i]}]}"
			pf="${_bams_seacr[${_pidx_seacr[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o="$(echo -e "$(basename "${nf%.*}")\t$(basename "${f%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done

		for i in "${!_nridx_seacr[@]}"; do
			nf="${_bams_seacr[${_nidx_seacr[$i]}]}"
			nrf="${_bams_seacr[${_nridx_seacr[$i]}]}"
			pf="${_bams_seacr[${_pidx_seacr[$i]}]}"

			x=$(( ${_nridx_seacr[$i]} + ${_pidx_seacr[$i]} ))
			pff="${_bams_seacr[$x]}"
			nff="${_bams_seacr[$((--x))]}"

			toidr=( "$odir/$(echo -e "$(basename "${nf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nrf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nff%.*}")\t$(basename "${pff%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	peaks::seacr \
		-S false \
		-s $skip \
		-c "$macsdir" \
		-m $memory \
		-M "$maxmemory" \
		-f $fragmentsize \
		-g "$genome" \
		-a nidx_seacr \
		-i tidx_seacr \
		-r ${!_mapper_seacr} \
		-t $threads \
		-o "$outdir" \
		-z $strict

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c idr -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function peaks::gopeaks_idr(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-m <memory>     | amount of. required if -c true
			-M <maxmemory>  | amount of. if -c true
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-f <size>       | assumed mean fragment. required if -c true
			-g <genome>     | path to. required if -c true
			-t <threads>    | number of
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to
			-w <broad>      | true/false peak detection
			-z <strict>     | true/false peak filters

			requires prior execution of alignment::mkreplicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir strict=false broad=false macsdir genome fragmentsize memory maxmemory
	declare -n _mapper_gopeaks _nidx_gopeaks _nridx_gopeaks _tidx_gopeaks _ridx_gopeaks _pidx_gopeaks
	while getopts 'S:s:c:m:M:f:g:t:r:a:b:i:j:k:o:w:z:' arg; do
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
			a)	((++mandatory)); _nidx_gopeaks=$OPTARG;;
			b)	_nridx_gopeaks=$OPTARG;;
			i)	((++mandatory)); _tidx_gopeaks=$OPTARG;;
			j)	((++mandatory)); _ridx_gopeaks=$OPTARG;;
			k)	((++mandatory)); _pidx_gopeaks=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			w)	broad=$OPTARG;;
			z)	strict=$OPTARG;;
			*) _usage;;
		esac
	done
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 7 ]] && _usage
	else
		[[ $mandatory -lt 10 ]] && _usage
	fi

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a nidx_gopeaks tidx_gopeaks cmd1 cmd2 toidr

	for m in "${_mapper_gopeaks[@]}"; do
		declare -n _bams_gopeaks=$m
		odir="$outdir/$m/gopeaks"
		nidx_gopeaks=()
		tidx_gopeaks=()

		for i in "${!_nidx_gopeaks[@]}"; do
			nidx_gopeaks+=(${_nidx_gopeaks[$i]} ${_nidx_gopeaks[$i]} ${_nidx_gopeaks[$i]})
			tidx_gopeaks+=(${_tidx_gopeaks[$i]} ${_ridx_gopeaks[$i]} ${_pidx_gopeaks[$i]})

			nf="${_bams_gopeaks[${_nidx_gopeaks[$i]}]}"
			tf="${_bams_gopeaks[${_tidx_gopeaks[$i]}]}"
			rf="${_bams_gopeaks[${_ridx_gopeaks[$i]}]}"
			pf="${_bams_gopeaks[${_pidx_gopeaks[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o="$(echo -e "$(basename "${nf%.*}")\t$(basename "${f%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done

		for i in "${!_nridx_gopeaks[@]}"; do
			nf="${_bams_gopeaks[${_nidx_gopeaks[$i]}]}"
			nrf="${_bams_gopeaks[${_nridx_gopeaks[$i]}]}"
			pf="${_bams_gopeaks[${_pidx_gopeaks[$i]}]}"

			x=$(( ${_nridx_gopeaks[$i]} + ${_pidx_gopeaks[$i]} ))
			pff="${_bams_gopeaks[$x]}"
			nff="${_bams_gopeaks[$((--x))]}"

			toidr=( "$odir/$(echo -e "$(basename "${nf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nrf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nff%.*}")\t$(basename "${pff%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	peaks::gopeaks \
		-S false \
		-s $skip \
		-c "$macsdir" \
		-m $memory \
		-M "$maxmemory" \
		-f $fragmentsize \
		-g "$genome" \
		-a nidx_gopeaks \
		-i tidx_gopeaks \
		-r ${!_mapper_gopeaks} \
		-t $threads \
		-o "$outdir" \
		-w $broad \
		-z $strict

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c idr -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function peaks::gem_idr(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-M <maxmemory>  | amount of
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-f <size>       | assumed mean fragment. required unless macs directory is given by -c
			-g <genome>     | path to
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to
			-q <ripseq>     | true/false
			-x <strandness> | hash per bam of (if ripseq)
			-y <pointy>     | true/false call only pointy peaks
			-z <strict>     | true/false peak filters

			requires prior execution of alignment::mkreplicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false ripseq=false threads memory maxmemory genome outdir strict=false pointy=false macsdir fragmentsize
	declare -n _mapper_gem _strandness_gem _nidx_gem _nridx_gem _tidx_gem _ridx_gem _pidx_gem
	while getopts 'S:s:c:f:t:m:M:g:q:r:x:a:b:i:j:k:o:z:y:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			c)	macsdir=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
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
			z)	$OPTARG && strict=true;;
			y)	pointy=$OPTARG;;
			*)	_usage;;
		esac
	done
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 9 ]] && _usage
	else
		[[ $mandatory -lt 10 ]] && _usage
	fi
	$ripseq && [[ ${#_strandness_gem[@]} -eq 0 ]] && _usage


	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a nidx_gem tidx_gem cmd1 cmd2 toidr

	for m in "${_mapper_gem[@]}"; do
		declare -n _bams_gem=$m
		odir="$outdir/$m/gem"
		nidx_gem=()
		tidx_gem=()

		for i in "${!_nidx_gem[@]}"; do
			nidx_gem+=(${_nidx_gem[$i]} ${_nidx_gem[$i]} ${_nidx_gem[$i]})
			tidx_gem+=(${_tidx_gem[$i]} ${_ridx_gem[$i]} ${_pidx_gem[$i]})

			nf="${_bams_gem[${_nidx_gem[$i]}]}"
			tf="${_bams_gem[${_tidx_gem[$i]}]}"
			rf="${_bams_gem[${_ridx_gem[$i]}]}"
			pf="${_bams_gem[${_pidx_gem[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o=$(echo -e "$(basename "${nf%.*}")\t$(basename "${f%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2/')
				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done

		for i in "${!_nridx_gem[@]}"; do
			nf="${_bams_gem[${_nidx_gem[$i]}]}"
			nrf="${_bams_gem[${_nridx_gem[$i]}]}"
			pf="${_bams_gem[${_pidx_gem[$i]}]}"

			x=$(( ${_nridx_gem[$i]} + ${_pidx_gem[$i]} ))
			pff="${_bams_gem[$x]}"
			nff="${_bams_gem[$((x-1))]}"

			toidr=( "$odir/$(echo -e "$(basename "${nf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nrf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nff%.*}")\t$(basename "${pff%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	peaks::gem \
		-S false \
		-s $skip \
		-q $ripseq \
		-g "$genome" \
		-c "$macsdir" \
		-f $fragmentsize \
		-a nidx_gem \
		-i tidx_gem \
		-r ${!_mapper_gem} \
		-x ${!_strandness_gem} \
		-t $threads \
		-m "$memory" \
		-M "$maxmemory" \
		-o "$outdir" \
		-y $pointy \
		-z $strict

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c idr -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function peaks::matk_idr(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-f <size>       | assumed mean fragment. required if -c true
			-g <genome>     | path to. required if -c true
			-t <threads>    | number of
			-m <memory>     | amount of. required if -c true
			-M <maxmemory>  | amount of
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to

			requires prior execution of alignment::mkreplicates
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir macsdir memory maxmemory fragmentsize genome
	declare -n _mapper_matk _nidx_matk _nridx_matk _tidx_matk _ridx_matk _pidx_matk
	while getopts 'S:s:c:m:f:g:t:M:r:x:a:b:i:j:k:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			c)	macsdir=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_matk=$OPTARG;;
			a)	((++mandatory)); _nidx_matk=$OPTARG;; # contains nridx by alignment::mkreplicates
			b)	_nridx_matk=$OPTARG;;
			i)	((++mandatory)); _tidx_matk=$OPTARG;;
			j)	((++mandatory)); _ridx_matk=$OPTARG;;
			k)	((++mandatory)); _pidx_matk=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 7 ]] && _usage
	else
		[[ $mandatory -lt 10 ]] && _usage
	fi

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a nidx_matk tidx_matk cmd1 cmd2 toidr

	for m in "${_mapper_matk[@]}"; do
		declare -n _bams_matk=$m
		odir="$outdir/$m/matk"
		nidx_matk=()
		tidx_matk=()

		for i in "${!_nidx_matk[@]}"; do
			nidx_matk+=(${_nidx_matk[$i]} ${_nidx_matk[$i]} ${_nidx_matk[$i]})
			tidx_matk+=(${_tidx_matk[$i]} ${_ridx_matk[$i]} ${_pidx_matk[$i]})

			nf="${_bams_matk[${_nidx_matk[$i]}]}"
			tf="${_bams_matk[${_tidx_matk[$i]}]}"
			rf="${_bams_matk[${_ridx_matk[$i]}]}"
			pf="${_bams_matk[${_pidx_matk[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o=$(echo -e "$(basename "${nf%.*}")\t$(basename "${f%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2/')
				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done

		for i in "${!_nridx_matk[@]}"; do
			nf="${_bams_matk[${_nidx_matk[$i]}]}"
			nrf="${_bams_matk[${_nridx_matk[$i]}]}"
			pf="${_bams_matk[${_pidx_matk[$i]}]}"

			x=$(( ${_nridx_matk[$i]} + ${_pidx_matk[$i]} ))
			pff="${_bams_matk[$x]}"
			nff="${_bams_matk[$((x-1))]}"

			toidr=( "$odir/$(echo -e "$(basename "${nf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nrf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nff%.*}")\t$(basename "${pff%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	peaks::matk \
		-S false \
		-s $skip \
		-c "$macsdir" \
		-m "$memory" \
		-f "$fragmentsize" \
		-g "$genome" \
		-a nidx_matk \
		-i tidx_matk \
		-r ${!_mapper_matk} \
		-t $threads \
		-M "$maxmemory" \
		-o "$outdir"

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c idr -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function peaks::peakachu_idr(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-c <macsdir>    | base directory where to find "macs" sub-directory according to used mappers (see -r) i.e. <macsdir>/<mapper>/macs
			-g <genome>     | path to. required unless macs directory is given by -c
			-m <memory>     | amount of. required unless macs directory is given by -c
			-M <maxmemory>  | amount of
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

	local OPTIND arg mandatory skip=false skipmd5=false threads fragmentsize outdir strict=false macsdir genome memory maxmemory
	declare -n _mapper_peakachu _nidx_peakachu _nridx_peakachu _tidx_peakachu _ridx_peakachu _pidx_peakachu
	while getopts 'S:s:c:m:M:g:t:f:r:a:b:i:j:k:o:z:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			c)	macsdir=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
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
	if [[ $macsdir ]]; then
		[[ $mandatory -lt 8 ]] && _usage
	else
		[[ $mandatory -lt 10 ]] && _usage
	fi

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a nidx_peakachu tidx_peakachu cmd1 cmd2 toidr

	for m in "${_mapper_peakachu[@]}"; do
		declare -n _bams_peakachu=$m
		odir="$outdir/$m/peakachu"
		nidx_peakachu=()
		tidx_peakachu=()

		for i in "${!_nidx_peakachu[@]}"; do
			nidx_peakachu+=(${_nidx_peakachu[$i]} ${_nidx_peakachu[$i]} ${_nidx_peakachu[$i]})
			tidx_peakachu+=(${_tidx_peakachu[$i]} ${_ridx_peakachu[$i]} ${_pidx_peakachu[$i]})

			nf="${_bams_peakachu[${_nidx_peakachu[$i]}]}"
			tf="${_bams_peakachu[${_tidx_peakachu[$i]}]}"
			rf="${_bams_peakachu[${_ridx_peakachu[$i]}]}"
			pf="${_bams_peakachu[${_pidx_peakachu[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o=$(echo -e "$(basename "${nf%.*}")\t$(basename "${f%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2/')
				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done

		for i in "${!_nridx_peakachu[@]}"; do
			nf="${_bams_peakachu[${_nidx_peakachu[$i]}]}"
			nrf="${_bams_peakachu[${_nridx_peakachu[$i]}]}"
			pf="${_bams_peakachu[${_pidx_peakachu[$i]}]}"

			x=$(( ${_nridx_peakachu[$i]} + ${_pidx_peakachu[$i]} ))
			pff="${_bams_peakachu[$x]}"
			nff="${_bams_peakachu[$((--x))]}"

			toidr=( "$odir/$(echo -e "$(basename "${nf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nrf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nff%.*}")\t$(basename "${pff%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	peaks::peakachu \
		-S false \
		-s $skip \
		-c "$macsdir" \
		-g "$genome" \
		-m "$memory" \
		-M "$maxmemory" \
		-f "$fragmentsize" \
		-a nidx_peakachu \
		-i tidx_peakachu \
		-r ${!_mapper_peakachu} \
		-t $threads \
		-o "$outdir" \
		-z $strict \
		-P

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c idr -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function peaks::genrich_idr(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-f <size>       | assumed mean fragment
			-t <threads>    | number of
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal replicates bam idices within -r (optional)
			-i <tidx>       | array of IP* bam idices within -r
			-j <ridx>       | array of IP* bam replicates idices within -r
			-k <pidx>       | array of IP* bam pools idices within -r
			-o <outdir>     | path to
			-q <ripseq>     | true/false
			-z <strict>     | true/false peak filters
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false genome fragmentsize threads outdir strict=false ripseq=false
	declare -n _mapper_genrich _nidx_genrich _nridx_genrich _tidx_genrich _ridx_genrich _pidx_genrich
	while getopts 'S:s:t:f:r:a:b:i:j:k:o:q:z:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			f)	((++mandatory)); fragmentsize=$OPTARG;;
			r)	((++mandatory)); _mapper_genrich=$OPTARG;;
			a)	((++mandatory)); _nidx_genrich=$OPTARG;;
			b)	_nridx_genrich=$OPTARG;;
			i)	((++mandatory)); _tidx_genrich=$OPTARG;;
			j)	((++mandatory)); _ridx_genrich=$OPTARG;;
			k)	((++mandatory)); _pidx_genrich=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			q)	ripseq=$OPTARG;;
			z)	strict=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 8 ]] && _usage

	local m i f b o odir nf tf rf pf nrf pff nff x
	declare -a nidx_genrich tidx_genrich cmd1 cmd2 toidr

	for m in "${_mapper_genrich[@]}"; do
		declare -n _bams_genrich=$m
		odir="$outdir/$m/genrich"
		nidx_genrich=()
		tidx_genrich=()

		for i in "${!_nidx_genrich[@]}"; do
			nidx_genrich+=(${_nidx_genrich[$i]} ${_nidx_genrich[$i]} ${_nidx_genrich[$i]})
			tidx_genrich+=(${_tidx_genrich[$i]} ${_ridx_genrich[$i]} ${_pidx_genrich[$i]})

			nf="${_bams_genrich[${_nidx_genrich[$i]}]}"
			tf="${_bams_genrich[${_tidx_genrich[$i]}]}"
			rf="${_bams_genrich[${_ridx_genrich[$i]}]}"
			pf="${_bams_genrich[${_pidx_genrich[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o=$(echo -e "$(basename "${nf%.*}")\t$(basename "${f%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2/')
				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done

		for i in "${!_nridx_genrich[@]}"; do
			nf="${_bams_genrich[${_nidx_genrich[$i]}]}"
			nrf="${_bams_genrich[${_nridx_genrich[$i]}]}"
			pf="${_bams_genrich[${_pidx_genrich[$i]}]}"

			x=$(( ${_nridx_genrich[$i]} + ${_pidx_genrich[$i]} ))
			pff="${_bams_genrich[$x]}"
			nff="${_bams_genrich[$((--x))]}"

			toidr=( "$odir/$(echo -e "$(basename "${nf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nrf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nff%.*}")\t$(basename "${pff%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	peaks::genrich \
		-S false \
		-s $skip \
		-f $fragmentsize \
		-a nidx_genrich \
		-i tidx_genrich \
		-r ${!_mapper_genrich} \
		-t $threads \
		-o "$outdir" \
		-q $ripseq \
		-z $strict \
		-P

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c idr -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function peaks::m6aviewer_idr(){
	function _usage(){
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
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false genome threads memory maxmemory fragmentsize outdir tmpdir="${TMPDIR:-/tmp}"
	declare -n _mapper_m6aviewer _nidx_m6aviewer _nridx_m6aviewer _tidx_m6aviewer _ridx_m6aviewer _pidx_m6aviewer
	while getopts 'S:s:t:m:M:f:r:a:b:i:j:k:o:' arg; do
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
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 9 ]] && _usage

	local m i f o odir nf tf rf pf nrf pff nff x
	declare -a nidx_m6aviewer tidx_m6aviewer cmd1 cmd2 toidr
	for m in "${_mapper_m6aviewer[@]}"; do
		declare -n _bams_m6aviewer=$m
		odir="$outdir/$m/m6aviewer"
		nidx_m6aviewer=()
		tidx_m6aviewer=()

		for i in "${!_nidx_m6aviewer[@]}"; do
			nidx_m6aviewer+=(${_nidx_m6aviewer[$i]} ${_nidx_m6aviewer[$i]} ${_nidx_m6aviewer[$i]})
			tidx_m6aviewer+=(${_tidx_m6aviewer[$i]} ${_ridx_m6aviewer[$i]} ${_pidx_m6aviewer[$i]})

			nf="${_bams_m6aviewer[${_nidx_m6aviewer[$i]}]}"
			tf="${_bams_m6aviewer[${_tidx_m6aviewer[$i]}]}"
			rf="${_bams_m6aviewer[${_ridx_m6aviewer[$i]}]}"
			pf="${_bams_m6aviewer[${_pidx_m6aviewer[$i]}]}"

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o="$(echo -e "$(basename "${nf%.*}")\t$(basename "${f%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2/')"
				toidr+=("$odir/$o.narrowPeak")
			done

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done

		for i in "${!_nridx_m6aviewer[@]}"; do
			nf="${_bams_m6aviewer[${_nidx_m6aviewer[$i]}]}"
			nrf="${_bams_m6aviewer[${_nridx_m6aviewer[$i]}]}"
			pf="${_bams_m6aviewer[${_pidx_m6aviewer[$i]}]}"

			x=$(( ${_nridx_m6aviewer[$i]} + ${_pidx_m6aviewer[$i]} ))
			pff="${_bams_m6aviewer[$x]}"
			nff="${_bams_m6aviewer[$((--x))]}"

			toidr=( "$odir/$(echo -e "$(basename "${nf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nrf%.*}")\t$(basename "${pf%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )
			toidr+=( "$odir/$(echo -e "$(basename "${nff%.*}")\t$(basename "${pff%.*}")" | sed -E 's/(\..+)\t(.+)\1/-\2.narrowPeak/')" )

			peaks::_idr \
				-1 cmd1 \
				-2 cmd2 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.idr"
		done
	done

	peaks::m6aviewer \
		-S false \
		-s $skip \
		-f $fragmentsize \
		-a nidx_m6aviewer \
		-i tidx_m6aviewer \
		-r ${!_mapper_m6aviewer} \
		-t $threads \
		-m $memory \
		-M "$maxmemory" \
		-o "$outdir"

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c idr -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}
