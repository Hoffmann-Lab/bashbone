#! /usr/bin/env bash
# (c) Konstantin Riege

callpeak::_idr() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-t <file>     | normal vs treatment peaks path to
			-r <file>     | normal vs replikate peaks path to
			-p <file>     | oracle: normal vs pool peaks path to
			-o <outfile>  | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory t r p o
	declare -n _cmds1_idr _cmds2_idr
	while getopts '1:2:t:r:p:o:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_idr=$OPTARG;;
			2) ((mandatory++)); _cmds2_idr=$OPTARG;;
			t) ((mandatory++)); t=$OPTARG;;
			r) ((mandatory++)); r=$OPTARG;;
			p) ((mandatory++)); p=$OPTARG;;
			o) ((mandatory++)); o="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage && return 1

	commander::makecmd -a _cmds1_idr -s '|' -c {COMMANDER[0]}<<- CMD
		idr
		--samples "$t" "$r"
		--peak-list "$p"
		--input-file-type narrowPeak
		--output-file "$o"
		--rank p.value
		--soft-idr-threshold 0.05
		--plot
		--use-best-multisummit-IDR
	CMD
	commander::makecmd -a _cmds2_idr -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
		cut -f 1-10 "$o"
	CMD
		sort -k1,1 -k2,2n -k3,3n
	CMD
		uniq > "$o.narrowPeak"
	CMD

	return 0
}

callpeak::macs() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-m <memory>     | amount of
			-f <size>       | assumed mean fragment
			-g <genome>     | path to
			-q <ripseq>     | true/false
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal repl bam idices within -r
			-i <tidx>       | no treat repl ? array to store new pseudorepl bam : array of treat bam idices within -r
			-j <ridx>       | no treat repl ? array to store new pseudorepl bam : array of treat repl bam idices within -r
			-k <pidx>       | no treat repl ? array of treat aka pseudopool bam : array to store new fullpool bam idices within -r
			-o <outdir>     | path to
			-p <tmpdir>     | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false ripseq=false genome threads memory fragmentsize outdir tmpdir
	declare -n _mapper_macs _macs _nidx_macs _nridx_macs _tidx_macs _ridx_macs _pidx_macs
	while getopts 'S:s:t:m:g:f:q:r:c:a:b:i:j:k:o:p:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			r) ((mandatory++)); _mapper_macs=$OPTARG;;
			f) ((mandatory++)); fragmentsize=$OPTARG;;
			q) ripseq=$OPTARG;;
			c) ((mandatory++))
				declare -n _caller_macs=$OPTARG
				_caller_macs+=(macs)
				_macs=macs
				;;
			a) ((mandatory++)); _nidx_macs=$OPTARG;;
			b) _nridx_macs=$OPTARG;;
			i) ((mandatory++)); _tidx_macs=$OPTARG;;
			j) ((mandatory++)); _ridx_macs=$OPTARG;;
			k) ((mandatory++)); _pidx_macs=$OPTARG;;
			o) ((mandatory++)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			p) ((mandatory++)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 12 ]] && _usage && return 1

	commander::printinfo "peak calling macs"
	commander::printinfo "preparing genome"

	# get effective genome size
	# if multimapped reads: genome minus Ns , else genome minus Ns minus repetetive Elements
	declare -n _bams_macs=${_mapper_macs[0]}
	local genomesize nf=${_bams_macs[${_nidx_macs[0]}]}
	if [[ $(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 256) -gt 0 ]]; then
		genomesize=$(faCount $genome | tail -1 | awk '{print $3+$4+$5+$6}')
	else
		genomesize=$(unique-kmers.py -k 100 $genome 2>&1 | tail -1 | awk '{print $NF}')
	fi

	local instances ithreads

	local params mult=5 # macs min mem in byte according to manual is numchr * buffer=100000 * mult=2 - observation max mem used ~ mult=5
	local numchr=$(samtools view -H results/mapped/segemehl/11q_HT.ip2.R1.final.sorted.bam | grep -c '^@SQ')
	local buffer=$(( (1024*1024*memory)/(numchr*mult) ))

	[[ $buffer -lt 100000 ]] && {
		params="--buffer-size $buffer"
		read -r instances ithreads < <(configure::instances_by_memory -i $threads -m $memory)
	} || {
		params=''
		read -r instances ithreads < <(configure::instances_by_memory -i $threads -m $(( (numchr*buffer*mult)/1024/1024 )))
	}

	local m i f o odir nf nrf tf rf x pff nff params2
	declare -a cmd1 cmd2 cmd3 cmd4 tdirs toidr
	for m in "${_mapper_macs[@]}"; do
		declare -n _bams_macs=$m
		odir=$outdir/$m/macs

		for i in "${!_nidx_macs[@]}"; do
			nf=${_bams_macs[${_nidx_macs[$i]}]}
			tf=${_bams_macs[${_tidx_macs[$i]}]}
			rf=${_bams_macs[${_ridx_macs[$i]}]}
			pf=${_bams_macs[${_pidx_macs[$i]}]}

			# infer SE or PE
			params+=' -f BAM'
			[[ $(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 1) -gt 0 ]] && params+=' -f BAMPE'

			toidr=()
			for f in $tf $rf $pf; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.macs)")
				o=$(echo -e "$(basename $nf)\t$(basename $f)" | sed -E 's/(.+)\t(.+)\1/-\2/')
				mkdir -p "$odir/$o"

				$ripseq && params2='--mfold 3 500'
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					macs2 callpeak
					$params
					-t "$f"
					-c "$nf"
					-g $genomesize
					--outdir "$odir/$o"
					-n "$o.model"
					--tempdir "${tdirs[-1]}"
					-B
					--SPMR
					--keep-dup all
					-q 0.05
					--verbose 1
					$params2
				CMD

				# shiftsize default 0 can not be altered for BAMPE and BEDPE
				$ripseq && params2="--shift $((-1*fragmentsize/2))" || params2="--shift 0"
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					macs2 callpeak
					$params
					-t "$f"
					-c "$nf"
					-g $genomesize
					--outdir "$odir/$o"
					-n "$o.nomodel"
					--tempdir "${tdirs[-1]}"
					-B
					--SPMR
					--keep-dup all
					-q 0.05
					--verbose 1
					--nomodel
					$params2
					--extsize $fragmentsize
				CMD

				commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
					bedtools merge
					-c 5,7,8,9
					-o collapse
					-i <(sort -k1,1 -k2,2n -k3,3n "$odir/$o/$o.model_peaks.narrowPeak" "$odir/$o/$o.nomodel_peaks.narrowPeak")
				CMD
					perl -M'List::Util qw(max)' -lane '
						print join "\t", (@F[0..2],"merged_peak_".(++$x),max(split/,/,$F[3]),".",max(split/,/,$F[4]),max(split/,/,$F[5]),max(split/,/,$F[6]),"-1")
					'
				CMD

				_macs+=("$odir/$o.narrowPeak")
				toidr+=("$odir/$o.narrowPeak")
			done

			callpeak::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.trp_idr"
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
			nf=${_bams_macs[${_nidx_macs[$i]}]} # 1
			nrf=${_bams_macs[${_nridx_macs[$i]}]} # 3
			pf=${_bams_macs[${_pidx_macs[$i]}]} # 9

			x=$(( ${_nridx_macs[$i]} + ${_pidx_macs[$i]} )) # 12
			pff=${_bams_macs[$x]}
			nff=${_bams_macs[$((--x))]} # 11

			toidr=( $odir/$(echo -e "$(basename $nf)\t$(basename $pf)" | sed -E 's/(.+)\t(.+)\1/-\2.narrowPeak/') )
			toidr+=( $odir/$(echo -e "$(basename $nrf)\t$(basename $pf)" | sed -E 's/(.+)\t(.+)\1/-\2.narrowPeak/') )
			toidr+=( $odir/$(echo -e "$(basename $nff)\t$(basename $pff)" | sed -E 's/(.+)\t(.+)\1/-\2.narrowPeak/') )

			callpeak::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.trp_idr"
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	} || {
		{	commander::runcmd -v -b -t $instances -a cmd1 && \
			commander::runcmd -v -b -t $threads -a cmd2 && \
			conda activate py3 && \
			commander::runcmd -v -b -t $threads -a cmd3 && \
			conda activate py2 && \
			commander::runcmd -v -b -t $threads -a cmd4
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}

callpeak::gem() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-g <genome>     | path to
			-f <gtf>        | path to
			-q <ripseq>     | true/false
			-r <mapper>     | array of sorted bams within array of
			-a <nidx>       | array of normal bam idices within -r
			-b <nridx>      | array of normal repl bam idices within -r
			-i <tidx>       | no treat repl ? array to store new pseudorepl bam : array of treat bam idices within -r
			-j <ridx>       | no treat repl ? array to store new pseudorepl bam : array of treat repl bam idices within -r
			-k <pidx>       | no treat repl ? array of treat aka pseudopool bam : array to store new fullpool bam idices within -r
			-o <outdir>     | path to
			-p <tmpdir>     | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false ripseq=false threads genome gtf outdir tmpdir
	declare -n _mapper_gem _gem _nidx_gem _nridx_gem _tidx_gem _ridx_gem _pidx_gem
	while getopts 'S:s:t:g:f:q:r:c:a:b:i:j:k:o:p:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			f) ((mandatory++)); gtf="$OPTARG";;
			r) ((mandatory++)); _mapper_gem=$OPTARG;;
			q) ripseq=$OPTARG;;
			c) ((mandatory++))
				declare -n _caller_gem=$OPTARG
				_caller_gem+=(gem)
				_gem=gem
				;;
			a) ((mandatory++)); _nidx_gem=$OPTARG;;
			b) _nridx_gem=$OPTARG;;
			i) ((mandatory++)); _tidx_gem=$OPTARG;;
			j) ((mandatory++)); _ridx_gem=$OPTARG;;
			k) ((mandatory++)); _pidx_gem=$OPTARG;;
			o) ((mandatory++)); outdir="$OPTARG";;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 11 ]] && _usage && return 1

	commander::printinfo "peak calling gem"
	commander::printinfo "preparing genome"

	local instances ithreads jmem jgct jcgct
	read -r instances ithreads jmem jgct jcgct < <(configure::jvm -t 20 -T $threads)

	# get effective genome size
	# if multimapped reads: genome minus Ns , else genome minus Ns minus repetetive Elements
	declare -n _bams_gem=${_mapper_gem[0]}
	local params genomesize nf=${_bams_gem[${_nidx_macs[0]}]}
	if [[ $(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 256) -gt 0 ]]; then
		genomesize=$(faCount $genome | tail -1 | awk '{print $3+$4+$5+$6}')
	else
		genomesize=$(unique-kmers.py -k 100 $genome 2>&1 | tail -1 | awk '{print $NF}')
	fi
	local tdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.genome)"
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
	$skip && {
		commander::printcmd -a cmdg
	} || {
		{	commander::runcmd -v -b -t $threads -a cmdg
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	local strandness=0
	$ripseq && {
		# --strand_type 1 disables search for asymmetry between sense and antisense mapped reads. instead leads to strand specific peak calls based on read orientation
		# if data is paired end, mates will be analyzed individually as single end reads
		# in case of PE this leads to duplicated peaks on both strands since orientation is not inferred from first in pair read
		# -> use this parameter only in case of strand specific SE reads
		[[ $(samtools view -F 4 "$nf" | head -10000 | cat <(samtools view -H "$nf") - | samtools view -c -f 1) -eq 0 ]] && {
			commander::printinfo "inferring library preparation method"

			local f tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.bed)"
			declare -a cmdi1 cmdi2
			alignment::_inferexperiment \
				-1 cmdi1 \
				-2 cmdi2 \
				-p "$tmp" \
				-i "$nf" \
				-g "$gtf"

			$skip && {
				commander::printcmd -a cmdi1
				commander::printcmd -a cmdi2
			} || {
				{	commander::runcmd -v -b -t $threads -a cmdi1 && \
					echo ":INFO: running commands of array cmdi2" && \
					commander::printcmd -a cmdi2 && \
					conda activate py3 && \
					read -r strandness f < <(commander::runcmd -t $threads -a cmdi2) && \
					conda activate py2
				} || {
					rm -f "$tmp"
					commander::printerr "$funcname failed"
					return 1
				}
			}
			rm -f "$tmp"
		}
		params='--relax --d $(dirname $(which gem))/Read_Distribution_CLIP.txt'
	} || {
		params="--d $(dirname $(which gem))/Read_Distribution_default.txt"
	}
	[[ $strandness -ne 0 ]] && params+=" --strand_type 1"

	local m i f o tdir odir nf nrf tf rf x pff nff
	declare -a cmd1 cmd2 cmd3 cmd4 toidr
	for m in "${_mapper_gem[@]}"; do
		declare -n _bams_gem=$m
		odir=$outdir/$m/gem

		for i in "${!_nidx_gem[@]}"; do
			nf=${_bams_gem[${_nidx_gem[$i]}]}
			tf=${_bams_gem[${_tidx_gem[$i]}]}
			rf=${_bams_gem[${_ridx_gem[$i]}]}
			pf=${_bams_gem[${_pidx_gem[$i]}]}

			toidr=()
			for f in "$tf" "$rf" "$pf"; do
				o=$(echo -e "$(basename $nf)\t$(basename $f)" | sed -E 's/(.+)\t(.+)\1/-\2/')
				mkdir -p "$odir/$o"

				commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					gem
						-Xmx${jmem}m
						-XX:ParallelGCThreads=$jgct
						-XX:ConcGCThreads=$jcgct
						-Djava.io.TMPDIR="$tmpdir"
						--t $ithreads
						--genome "$tdir"
						--g "$tdir/chr.info"
						--out "$odir/$o"
						--expt "$f"
						--ctrl "$nf"
						--f SAM
						--nrf
						--sl
						--s $genomesize
						--q $(echo 0.05 | awk '{print -log($1)/log(10)}')
						--outNP
						$params
				CMD
					cp "$odir/$o/$o.GPS_events.narrowPeak" "$odir/$o.narrowPeak"
				CMD

				[[ $strandness -eq 0 ]] && {
					commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
						bedtools merge
						-c 5,7,8,9
						-o collapse
						-i <(awk '$2>=0 && $3>=0' "$odir/$o/$o.GPS_events.narrowPeak" | sed -E '/^chr(M|X|Y|[0-9]+)/!{s/^chr//}')
					CMD
						perl -M'List::Util qw(max)' -lane '
							print join "\t", (@F[0..2],"merged_peak_".(++$x),max(split/,/,$F[3]),".",max(split/,/,$F[4]),max(split/,/,$F[5]),max(split/,/,$F[6]),"-1")
						'
					CMD
				} || {
					commander::makecmd -a cmd2 -s '|' -o "$odir/$o.narrowPeak" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
						bedtools merge
						-s
						-c 5,6,7,8,9
						-o collapse
						-i <(awk '$2>=0 && $3>=0' "$odir/$o/$o.GPS_events.narrowPeak" | sed -E '/^chr(M|X|Y|[0-9]+)/!{s/^chr//}')
					CMD
						perl -M'List::Util qw(max)' -lane '
							print join "\t", (@F[0..2],"merged_peak_".(++$x),max(split/,/,$F[3]),(split/,/,$F[4])[0],max(split/,/,$F[5]),max(split/,/,$F[6]),max(split/,/,$F[7]),"-1")
						'
					CMD
				}

				_gem+=("$odir/$o.narrowPeak")
				toidr+=("$odir/$o.narrowPeak")
			done

			callpeak::_idr \
				-1 cmd3 \
				-2 cmd4 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.trp_idr"
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
			nf=${_bams_gem[${_nidx_gem[$i]}]} # 1
			nrf=${_bams_gem[${_nridx_gem[$i]}]} # 3
			pf=${_bams_gem[${_pidx_gem[$i]}]} # 9

			x=$(( ${_nridx_gem[$i]} + ${_pidx_gem[$i]} )) # 12
			pff=${_bams_gem[$x]}
			nff=${_bams_gem[$((x-1))]} # 11

			toidr=( $odir/$(echo -e "$(basename $nf)\t$(basename $pf)" | sed -E 's/(.+)\t(.+)\1/-\2.narrowPeak/') )
			toidr+=( $odir/$(echo -e "$(basename $nrf)\t$(basename $pf)" | sed -E 's/(.+)\t(.+)\1/-\2.narrowPeak/') )
			toidr+=( $odir/$(echo -e "$(basename $nff)\t$(basename $pff)" | sed -E 's/(.+)\t(.+)\1/-\2.narrowPeak/') )

			callpeak::_idr \
				-1 cmd2 \
				-2 cmd3 \
				-t "${toidr[0]}" \
				-r "${toidr[1]}" \
				-p "${toidr[2]}" \
				-o "${toidr[2]%.*}.trp_idr"
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	} || {
		{	commander::runcmd -v -b -t $instances -a cmd1 && \
			commander::runcmd -v -b -t $threads -a cmd2 && \
			conda activate py3 && \
			commander::runcmd -v -b -t $threads -a cmd3 && \
			conda activate py2 && \
			commander::runcmd -v -b -t $threads -a cmd4
		} || {
			rm -rf "$tdir"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "$tdir"
	return 0
}
