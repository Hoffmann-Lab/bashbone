#! /usr/bin/env bash
# (c) Konstantin Riege

alignment::mkreplicates() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-r <mapper>     | array of sorted bams within array of
			-n <nidx>       | array of normal bam idices within -r
			-m <nridx>      | array of normal repl bam idices within -r
			-i <tidx>       | no treat repl ? array to store new pseudorepl bam : array of treat bam idices within -r
			-j <ridx>       | no treat repl ? array to store new pseudorepl bam : array of treat repl bam idices within -r
			-k <pidx>       | no treat repl ? array of treat aka pseudopool bam : array to store new fullpool bam idices within -r
			-o <outdir>     | path to
			-p <tmpdir>     | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads outdir tmpdir
	declare -n _mapper_mkreplicates _nidx_mkreplicates _nridx_mkreplicates _tidx_mkreplicates _ridx_mkreplicates _pidx_mkreplicates
	while getopts 'S:s:t:r:n:m:i:j:k:o:p:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			r) ((mandatory++)); _mapper_mkreplicates=$OPTARG;;
			n) ((mandatory++)); _nidx_mkreplicates=$OPTARG;;
			m) _nridx_mkreplicates=$OPTARG;;
			i) ((mandatory++)); _tidx_mkreplicates=$OPTARG;;
			j) ((mandatory++)); _ridx_mkreplicates=$OPTARG;;
			k) ((mandatory++)); _pidx_mkreplicates=$OPTARG;;
			o) ((mandatory++)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			p) ((mandatory++)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 8 ]] && _usage && return 1

	local m i odir o tmp nf nrf tf rf pf addindex=true ithreads1 ithreads2 instances1=1 instances2=1
	declare -a cmd1 cmd2 cmd3 tdirs
	if [[ $_ridx_mkreplicates ]]; then
		commander::printinfo "generating pseudo-pools"
		instances1=$((${#_mapper_mkreplicates[@]} * ${#_nidx_mkreplicates[@]} * 2 + ${#_mapper_mkreplicates[@]} * ${#_nridx_mkreplicates[@]} * 2))
		instances2=$((${#_mapper_mkreplicates[@]} * ${#_nidx_mkreplicates[@]}))
		read -r instances1 ithreads1 < <(configure::instances_by_threads -i $instances1 -t 10 -T $threads)
		read -r instances2 ithreads2 < <(configure::instances_by_threads -i $instances2 -t 10 -T $threads)

		# pool replicates:
		# m[N1 N2 T1 T2 R1 R2] -> m[N1 N2 T1 T2 R1 R2 PP1 PP2]
		# n   1 2
		# nr
		# t   3 4
		# r   5 6
		# p   7 8
		# -> idr: 1 vs 3 + 1 vs 5 + 1 vs 7
		#          2 vs 4 + 2 vs 6 + 2 vs 8
		for m in "${_mapper_mkreplicates[@]}"; do
			declare -n _bams_mkreplicates=$m
			odir=$outdir/$m
			mkdir -p "$odir"
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.mkreplicates)")
			for i in "${!_nidx_mkreplicates[@]}"; do
				tf=${_bams_mkreplicates[${_tidx_mkreplicates[$i]}]}
				rf=${_bams_mkreplicates[${_ridx_mkreplicates[$i]}]}
				o=$odir/$(echo -e "$(basename $tf)\t$(basename $rf)" | sed -E 's/(.+)\t(.+)\.\1/-\2.pseudopool.\1/')

				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					samtools view -@ $ithreads1 -b -s 0.5 $tf > "${tdirs[-1]}/$(basename "$tf")"
				CMD
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					samtools view -@ $ithreads1 -b -s 0.5 $rf > "${tdirs[-1]}/$(basename "$rf")"
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					samtools merge -f -@ $ithreads2 "$o" "${tdirs[-1]}/$(basename "$tf")" "${tdirs[-1]}/$(basename "$rf")"
				CMD

				_bams_mkreplicates+=("$o")
				_pidx_mkreplicates+=($((${#_bams_mkreplicates[@]}-1)))
			done
		done

		if [[ $_nridx_mkreplicates ]]; then
			# create fully pooled bams and extend _nidx_mkreplicates by indices for n/nr/nfp based peak calling
			# PP : pseudopool (2x0.5) , FP: fullpool (2x1)
			#                                           1  2   3   4  5  6  7  8  9   10   11  12   13  14
			# m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2] -> m[N1 N2 NR1 NR2 T1 T2 R1 R2 PP1 PP2 NFP1 FP1 NFP2 FP2]
			# n   1 2   3 4   11 13
			# nr  3 4
			# t   5 6   5 6   5  6
			# r   7 8   7 8   7  8
			# p   9 10  9 10  12 14
			# -> idr: 1 vs 7 + 1 vs 7 + 1 vs 9
			#          2 vs 6 + 2 vs 8 + 2 vs 10
			#         3 vs 5 + 3 vs 6 + 3 vs 9
			#          4 vs 6 + 4 vs 8 + 4 vs 10
			#         11 vs 5 + 11 vs 7 + 11 vs 12
			#          13 vs 6 + 13 vs 8 + 13 vs 14
			_nidx_mkreplicates=("${_nidx_mkreplicates[@]}" "${_nridx_mkreplicates[@]}")
			_tidx_mkreplicates+=("${_tidx_mkreplicates[@]}" "${_tidx_mkreplicates[@]}")
			_ridx_mkreplicates+=("${_ridx_mkreplicates[@]}" "${_ridx_mkreplicates[@]}")
			_pidx_mkreplicates+=("${_pidx_mkreplicates[@]}")
			for m in "${_mapper_mkreplicates[@]}"; do
				declare -n _bams_mkreplicates=$m
				odir=$outdir/$m
				mkdir -p $odir
				for i in "${!_nridx_mkreplicates[@]}"; do
					nf=${_bams_mkreplicates[${_nidx_mkreplicates[$i]}]}
					nrf=${_bams_mkreplicates[${_nridx_mkreplicates[$i]}]}
					o=$odir/$(echo -e "$(basename $nf)\t$(basename $nrf)" | sed -E 's/(.+)\t(.+)\.\1/-\2.fullpool.\1/')

					commander::makecmd -a cmd3 -s '&&' -c {COMMANDER[0]}<<- CMD
						samtools merge -f -@ $ithreads1 $o $nf $nrf
					CMD
					_bams_mkreplicates+=("$o")
					$addindex && _nidx_mkreplicates+=($((${#_bams_mkreplicates[@]}-1)))

					tf=${_bams_mkreplicates[${_tidx_mkreplicates[$i]}]}
					rf=${_bams_mkreplicates[${_ridx_mkreplicates[$i]}]}
					o=$odir/$(echo -e "$(basename $tf)\t$(basename $rf)" | sed -E 's/(.+)\t(.+)\.\1/-\2.fullpool.\1/')
					commander::makecmd -a cmd3 -s '&&' -c {COMMANDER[0]}<<- CMD
						samtools merge -f -@ $ithreads1 $o $tf $rf
					CMD
					_bams_mkreplicates+=("$o")
					$addindex && _pidx_mkreplicates+=($((${#_bams_mkreplicates[@]}-1)))
				done
				addindex=false
			done
		fi
	else
		commander::printinfo "generating pseudo-replicates"
		instances1=$((${#_mapper_mkreplicates[@]} * ${#_nidx_mkreplicates[@]}))
		read -r instances1 ithreads1 < <(configure::instances_by_threads -i $instances1 -t 10 -T $threads)
		# make pseudo-replicates from pseudo-pool:
		# m[N1 N2 P1 P2] -> m[N1 N2 P1 P2 T1 R1 T2 R2]
		# cmp 1 2
		# n   1 2
		# nr
		# t   5 7
		# r   6 8
		# p   3 4
		# -> idr: 1 vs 5 + 1 vs 6 + 1 vs 3
		#          2 vs 7 + 2 vs 8 + 2 vs 4
		for m in "${_mapper_mkreplicates[@]}"; do
			declare -n _bams_mkreplicates=$m
			odir=$outdir/$m
			mkdir -p $odir $tdir

			for i in "${!_nidx_mkreplicates[@]}"; do
				pf=${_bams_mkreplicates[${_pidx_mkreplicates[$i]}]}
				o=$odir/$(basename ${pf%.*}.pseudorep)
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.mkreplicates)")

				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					samtools collate
						-@ $ithreads1
						-n $((threads<64?64:threads))
						-u
						-O
						$pf ${tdirs[-1]}/$(basename ${pf%.*})
				CMD
					samtools view -@ $ithreads1
				CMD
					split
						--numeric-suffixes=1
						--additional-suffix=.bam
						-a 1
						-l \$(( (\$(samtools view -c -@ $ithreads1 $pf)+1)/2 ))
						--filter='cat <(samtools view -H $pf) - | samtools sort -@ $ithreads1 -O BAM -T ${tdirs[-1]}/\$(basename \${FILE%.*}) > \$FILE'
						- $o
				CMD
				# -n number of temporary files has default 64
				# -u not equals --output-fmt SAM => a small compression level (optimum: -l 3) reduces amount of data stream through pipe
				# which finally makes samtools view faster despite of decompression - tested with 650M bam (~1m30 vs ~1m)
				# samtools view necessary, because it removes header - grep -v '^@' may match CIGAR string
				_bams_mkreplicates+=("${o}1.bam")
				_tidx_mkreplicates+=($((${#_bams_mkreplicates[@]}-1)))
				_bams_mkreplicates+=("${o}2.bam")
				_ridx_mkreplicates+=($((${#_bams_mkreplicates[@]}-1)))
			done
		done
	fi

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	} || {
		{	commander::runcmd -v -b -t $instances1 -a cmd1 && \
			commander::runcmd -v -b -t $instances2 -a cmd2 && \
			commander::runcmd -v -b -t $instances1 -a cmd3
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}
