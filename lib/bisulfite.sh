#! /usr/bin/env bash
# (c) Konstantin Riege

function bisulfite::mspicut(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-d <length>   | maximum of diversity adapters
			-c <number>   | bases to clip from R2 start (default: 2)
			-o <outdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir diversity=0 clip=2
	declare -n _fq1_mspicut _fq2_mspicut
	while getopts 'S:s:t:d:c:o:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			d) diversity=$OPTARG;;
			c) clip=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1) ((++mandatory)); _fq1_mspicut=$OPTARG;;
			2) _fq2_mspicut=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "selecting msp1 cut reads"

	declare -a cmd1
	local i o1 o2 e1 e2
	[[ $n -gt 2 ]] && n=2 # since only the best matching adapter is removed, run cutadapt twice
	for i in "${!_fq1_mspicut[@]}"; do
		helper::basename -f "${_fq1_mspicut[$i]}" -o o1 -e e1
		e1=$(echo $e1 | cut -d '.' -f 1)
		o1="$outdir/$o1.$e1.gz"


		if [[ "${_fq2_mspicut[$i]}" ]]; then
			helper::basename -f "${_fq2_mspicut[$i]}" -o o2 -e e2
			e2=$(echo $e2 | cut -d '.' -f 1)
			o2="$outdir/$o2.$e2.gz"

			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				rrbsMspIselection.sh
				-t $threads
				-d $diversity
				-i "${_fq1_mspicut[$i]}"
				-j "${_fq2_mspicut[$i]}"
				-o "$o1"
				-p "$o2"
				-c $clip
			CMD
			# do not clip off 2 nt from R2 starts, it will be done by cutadapt in case of nomspi=true (required for multi enzyme digestion)
			_fq1_mspicut[$i]="$o1"
			_fq2_mspicut[$i]="$o2"
		else
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				rrbsMspIselection.sh
				-t $threads
				-d $diversity
				-i "${_fq1_mspicut[$i]}"
				-o "$o1"
			CMD
			_fq1_mspicut[$i]="$o1"
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i 1 -a cmd1
	fi

	return 0
}

function bisulfite::segemehl(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-5 <skip>       | true/false md5sums, indexing respectively
			-t <threads>    | number of
			-a <accuracy>   | optional: 80 to 100
			-i <insertsize> | optional: 50 to 200000+
			-r <mapper>     | array of bams within array of
			-g <genome>     | path to
			-x <ctidx>      | path to
			-y <gaidx>      | path to
			-m <mode>       | lister (1), cokus (2) - default: 1
			-o <outdir>     | path to
			-1 <fastq1>     | array of
			-2 <fastq2>     | array of
			-F              | force index
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome gaidx ctidx outdir tmpdir="${TMPDIR:-/tmp}" mode=1 accuracy insertsize forceidx=false
	declare -n _fq1_segemehl _fq2_segemehl _mapper_segemehl
	declare -g -a segemehl=()
	while getopts 'S:s:5:t:a:i:g:x:y:m:r:o:1:2:F' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	mode=$OPTARG;;
			a)	accuracy=$OPTARG;;
			i)	insertsize=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			x)	((++mandatory)); ctidx="$OPTARG";;
			y)	((++mandatory)); gaidx="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/segemehl"; mkdir -p "$outdir";;
			r)	((++mandatory))
				_mapper_segemehl=$OPTARG
				_mapper_segemehl+=(segemehl)
			;;
			1)	((++mandatory)); _fq1_segemehl=$OPTARG;;
			2)	_fq2_segemehl=$OPTARG;;
			F)	forceidx=true;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 7 ]] && _usage

	commander::printinfo "bisufite mapping segemehl"

	if $skipmd5; then
		commander::warn "skip checking md5 sums and genome indexing respectively"
	else
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$BASHBONE_DIR/lib/md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		local thismd5genome thismd5segemehlbs
		thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
		[[ -s "$ctidx" ]] && thismd5segemehlbs=$(md5sum "$ctidx" | cut -d ' ' -f 1)
		if $forceidx || [[ "$thismd5genome" != "$md5genome" || ! "$thismd5segemehlbs" || "$thismd5segemehlbs" != "$md5segemehlbs" ]]; then
			commander::printinfo "indexing genome for segemehl bisulfite mode"
			declare -a cmdidx
			# lister/cokus mode use the same index, thus it does not matter if indexed with -F 1 or -F 2
			commander::makecmd -a cmdidx -s '|' -c {COMMANDER[0]}<<- CMD
				segemehl -F 1 -x "$ctidx" -y "$gaidx" -d "$genome"
			CMD
			commander::runcmd -v -b -i $threads -a cmdidx
			commander::printinfo "updating md5 sums"
			thismd5segemehlbs=$(md5sum "$ctidx" | cut -d ' ' -f 1)
			sed -i "s/md5segemehlbs=.*/md5segemehlbs=$thismd5segemehlbs/" $genome.md5.sh
			sed -i "s/md5genome=.*/md5genome=$thismd5genome/" "$genome.md5.sh"
		fi
	fi

	declare -a tdirs cmd1
	local o e params
	for i in "${!_fq1_segemehl[@]}"; do
		helper::basename -f "${_fq1_segemehl[$i]}" -o o -e e
		o="$outdir/$o"
		params=""
		[[ $accuracy ]] && params+=" -A $accuracy"
		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.segemehl)")
		if [[ ${_fq2_segemehl[$i]} ]]; then
			[[ $insertsize ]] && params+=" -I $insertsize"
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cd "${tdirs[-1]}"
			CMD
				segemehl
				$params
				-F $mode
				-i "$(realpath -se "$ctidx")"
				-j "$(realpath -se "$gaidx")"
				-d "$(realpath -se "$genome")"
				-q "$(realpath -se "${_fq1_segemehl[$i]}")"
				-p "$(realpath -se "${_fq2_segemehl[$i]}")"
				-t $threads
				-b
				-o "$(realpath -s "$o.bam")"
			CMD
		else
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cd "${tdirs[-1]}"
			CMD
				segemehl
				$params
				-F $mode
				-i "$(realpath -se "$ctidx")"
				-j "$(realpath -se "$gaidx")"
				-d "$(realpath -se "$genome")"
				-q "$(realpath -se "${_fq1_segemehl[$i]}")"
				-t $threads
				-b
				-o "$(realpath -s "$o.bam")"
			CMD
		fi
		segemehl+=("$o.bam")
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i 1 -a cmd1
	fi

	return 0
}

function bisulfite::bwa(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-5 <skip>       | true/false md5sums, indexing respectively
			-t <threads>    | number of
			-a <accuracy>   | optional: 80 to 100
			-r <mapper>     | array of bams within array of
			-g <genome>     | path to
			-o <outdir>     | path to
			-1 <fastq1>     | array of
			-2 <fastq2>     | array of
			-F              | force index
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome outdir forceidx=false accuracy
	declare -n _fq1_bwa _fq2_bwa _mapper_bwa
	declare -g -a bwa=()
	while getopts 'S:s:5:t:a:g:r:o:1:2:F' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/bwa"; mkdir -p "$outdir";;
			a)	accuracy=$OPTARG;;
			r)	((++mandatory))
				_mapper_bwa=$OPTARG
				_mapper_bwa+=(bwa)
			;;
			1)	((++mandatory)); _fq1_bwa=$OPTARG;;
			2)	_fq2_bwa=$OPTARG;;
			F)	forceidx=true;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "bisufite mapping bwa"
	if $skipmd5; then
		commander::warn "skip checking md5 sums and genome indexing respectively"
	else
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$BASHBONE_DIR/lib/md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		local thismd5genome thismd5bwameth
		thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
		[[ -s "$genome.bwameth.c2t.pac" ]] && thismd5bwameth=$(md5sum "$genome.bwameth.c2t.pac" | cut -d ' ' -f 1)
		if $forceidx || [[ "$thismd5genome" != "$md5genome" || ! "$thismd5bwameth" || "$thismd5bwameth" != "$md5bwameth" ]]; then
			commander::printinfo "indexing genome for bwameth"
			declare -a cmdidx
			commander::makecmd -a cmdidx -s ';' -c {COMMANDER[0]}<<- CMD
				bwameth.py index-mem2 "$genome"
			CMD
			commander::runcmd -c bwameth -v -b -i $threads -a cmdidx
			commander::printinfo "updating md5 sums"
			thismd5bwameth=$(md5sum "$genome.bwameth.c2t.pac" | cut -d ' ' -f 1)
			sed -i "s/md5bwameth=.*/md5bwameth=$thismd5bwameth/" "$genome.md5.sh"
			sed -i "s/md5genome=.*/md5genome=$thismd5genome/" "$genome.md5.sh"
		fi
	fi

	declare -a cmd1
	local o e readlength params=""
	for i in "${!_fq1_bwa[@]}"; do
		helper::basename -f "${_fq1_bwa[$i]}" -o o -e e
		o="$outdir/$o"

		helper::makecatcmd -c catcmd -f "${_fq1_bwa[$i]}"
		readlength=$($catcmd "${_fq1_bwa[$i]}" | head -4000 | awk 'NR%4==2{l+=length($0)}END{printf("%.f",l/(NR/4))}')
		[[ $accuracy ]] && params='--score '$(echo $accuracy | awk -v l=$readlength '{printf("%.f",l-l*(1-$1/100)*5)}')

		if [[ ${_fq2_bwa[$i]} ]]; then
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD
				bwameth.py
					$params
					--do-not-penalize-chimeras
					--read-group '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
					--threads $threads
					--reference "$genome"
					"${_fq1_bwa[$i]}" "${_fq2_bwa[$i]}"
			CMD
				sed 's/\t\t/\t*\t/'
			CMD
				sed -E -e ':a; s/(\s[SX]A:Z\S*[:;])[fr]/\1/; t a' -e 's/\tYC:Z:(\S+)/\tHI:i:0\tXD:i:0\tXF:i:0\tXB:Z:F1\/\1\tYZ:Z:0\tYC:Z:\1/'
			CMD
				samtools view -@ $threads -b > "$o.bam"
			CMD
			# tee >(samtools view -@ $((threads/2+1)) -b > "$o.raw.bam")
			# sed corrects SA and XA tags, which bwameth does not change back from c2t (f+r) chromosomes i.e. (f|r)chr to chr
		else
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD
				bwameth.py
					$params
					--do-not-penalize-chimeras
					--read-group '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
					--threads $threads
					--reference "$genome"
					"${_fq1_bwa[$i]}"
			CMD
				sed 's/\t\t/\t*\t/'
			CMD
				sed -E -e ':a; s/(\s[SX]A:Z\S*[:;])[fr]/\1/; t a' -e 's/\tYC:Z:(\S+)/\tHI:i:0\tXD:i:0\tXF:i:0\tXB:Z:F1\/\1\tYZ:Z:0\tYC:Z:\1/'
			CMD
				samtools view -@ $threads -b > "$o.bam"
			CMD
			# tee >(samtools view -@ $((threads/2+1)) -b > "$o.raw.bam")
		fi
		bwa+=("$o.bam")
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c bwameth -v -b -i 1 -a cmd1
	fi

	return 0
}

#todo gem: gem-mapper -p --bisulfite-mode -I c2t+g2a.fa.gem -s 1 -p -M 4 && bs_call -r genome.fa.gz -p -L5

function bisulfite::rmduplicates(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-k             | keep marked duplicates in bam
			-t <threads>   | number of
			-m <memory>    | amount of
			-M <maxmemory> | amount of
			-g <genome>    | mandatory
			-r <mapper>    | array of sorted, indexed bams within array of
			-3 <fastqUMI>  | array of
			-c <sliceinfo> | array of
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory maxmemory genome outdir remove=true tmpdir="${TMPDIR:-/tmp}"
	declare -n _mapper_rmduplicates _bamslices_rmduplicates _umi_rmduplicates
	while getopts 'S:s:t:m:M:g:r:3:c:o:k' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			k)	remove=false;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			r)	((++mandatory)); _mapper_rmduplicates=$OPTARG;;
			3)	_umi_rmduplicates=$OPTARG;;
			c)	((++mandatory)); _bamslices_rmduplicates=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "removing duplicates"

	local sinstances sthreads smemory minstances mthreads jmem jgct jcgct
	read -r sinstances sthreads smemory jgct jcgct < <(configure::jvm -i ${#_umi_rmduplicates[@]} -T $threads -m $memory -M "$maxmemory")
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory -M "$maxmemory")

	local m i o e slice instances ithreads odir params1 params2 x=0 catcmd oinstances othreads
	for m in "${_mapper_rmduplicates[@]}"; do
		declare -n _bams_rmduplicates=$m
		i=$(wc -l < "${_bamslices_rmduplicates[${_bams_rmduplicates[0]}]}")
		((instances+=i*${#_bams_rmduplicates[@]}))
		((oinstances+=${#_bams_rmduplicates[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)
	read -r oinstances othreads < <(configure::instances_by_threads -i $oinstances -t 10 -T $threads)

	declare -a cmdsort
	if [[ $_umi_rmduplicates ]]; then
		params2='-B'

		for i in "${!_umi_rmduplicates[@]}"; do
			helper::basename -f "${_umi_rmduplicates[$i]}" -o o -e e
			e=$(echo $e | cut -d '.' -f 1)
			o="$tmpdir/$o.$e.gz"

			helper::makecatcmd -c catcmd -f "${_umi_rmduplicates[$i]}"

			commander::makecmd -a cmdsort -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD
				$catcmd "${_umi_rmduplicates[$i]}" | paste - - - -
			CMD
				awk -v OFS='\t' '{print $1,$(NF-2),$(NF-1),$NF}'
			CMD
				LC_ALL=C sort --parallel=$sthreads -S ${smemory}M -T "$tmpdir" -k1,1V
			CMD
				tr '\t' '\n'
			CMD
				helper::pgzip -t $sthreads -o "$o"
			CMD
			_umi_rmduplicates[$i]="$o"
		done
	fi

	x=$(samtools view -F 4 "${_bams_rmduplicates[0]}" | head -10000 | cat <(samtools view -H "${_bams_rmduplicates[0]}") - | samtools view -c -f 1)
	[[ $x -gt 0 ]] && params1='--paired' || params2+=' -s'
	$remove && params2+=' -r'

	declare -a tdirs tomerge cmd1 cmd2 cmd3 cmd4 cmd5
	for m in "${_mapper_rmduplicates[@]}"; do
		declare -n _bams_rmduplicates=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_bams_rmduplicates[@]}"; do
			tomerge=()
			while read -r slice; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.rmduplicates)")

				if [[ $_umi_rmduplicates ]]; then
					commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- CMD
						samtools sort
							-n
							-@ $ithreads
							-u
							-T "${tdirs[-1]}/$(basename "${slice%.*}")"
							"$slice"
					CMD
						| samtools fixmate
							-@ $ithreads
							-r
							-O SAM
							- -
					CMD
						| perl -slane '
							BEGIN{
								open(F, "pigz -p 1 -cd \"$f\" | paste - - - - |")
							}
							if(/^@\S\S\s/){print; next}
							while($l[0] ne $F[0]){
								$l=<F>;
								exit unless defined $l;
								@l=split(/\t/,$l);
								$l[0]=~s/^.//
							}
							print join"\t",(@F,"RX:Z:$l[1]")
						'
					CMD
						-- -f="${_umi_rmduplicates[$i]}"
					CMD
						| samtools sort
							-@ $ithreads
							-O BAM
							-T "${tdirs[-1]}/$(basename "${slice%.*}").rx"
						> "$slice.rx";
					CMD
						mv "$slice.rx" "$slice";
					CMD
						samtools index -@ $ithreads "$slice" "${slice%.*}.bai"
					CMD

					# misuse corrected umi group tag (default often MI as corrected barcode tag CB for dupsifter -B == --has-barcodes - searches for CB tag)
					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
						umi_tools group
							-I "$slice"
							-S "$slice.grouped"
							--log2stderr
							--output-bam
							--extract-umi-method tag
							--umi-tag RX
							--umi-group-tag CB
							--temp-dir "${tdirs[-1]}"
							--method directional
							--edit-distance-threshold 1
							--random-seed 12345
							--no-sort-output
							$params1
					CMD
						mv "$slice.grouped" "$slice"
					CMD
				fi

				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD
					samtools sort
						-@ $ithreads
						-n
						-O SAM
						-T "${tdirs[-1]}/$(basename "${slice%.*}")"
						"$slice"
				CMD
					perl -F'\t' -lane '
						if(/^@\S\S\s/){print; next}
						$F[11]="EC:Z:$F[5]\t$F[11]";
						$F[5]=~s/[=X]/M/g;
						print join"\t",@F
					'
				CMD
					dupsifter
						-o "/dev/stdout"
						-O "/dev/stderr"
						$params2
						"$genome"
				CMD
					perl -F'\t' -lane '
						if(/^@\S\S\s/){print; next}
						$F[5]=substr($F[11],5);
						print join"\t",@F[0..10,12..$#F]
					'
				CMD
					samtools fixmate
						-@ $ithreads
						-r
						-u
						- -
				CMD
					samtools sort
						-@ $ithreads
						-O BAM
						-T "${tdirs[-1]}/$(basename "${slice%.*}").rmdup"
					> "$slice.rmdup"
				CMD
				# dupsifter handles only sam v 1.3 not extended cigar by X and =
				# solution: replace X and = by M, store extended cigar as EC tag and afterwards sqeeze back in

				commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					mv "$slice.rmdup" "$slice"
				CMD
					samtools index -@ $ithreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_rmduplicates[${_bams_rmduplicates[$i]}]}"

			o="$odir/$(basename "${_bams_rmduplicates[$i]}")"
			o="${o%.*}.rmdup.bam"

			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $othreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_rmduplicates["$o"]="${_bamslices_rmduplicates[${_bams_rmduplicates[$i]}]}"
			_bams_rmduplicates[$i]="$o"
		done
	done

	if $skip; then
		commander::printcmd -a cmdsort
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
	else
		commander::runcmd -v -b -i $sinstances -a cmdsort
		commander::runcmd -v -b -i $instances -a cmd1
		commander::runcmd -c umitools -v -b -i $minstances -a cmd2
		commander::runcmd -c dupsifter -v -b -i $instances -a cmd3
		commander::runcmd -v -b -i $instances -a cmd4
		commander::runcmd -v -b -i $oinstances -a cmd5
	fi

	return 0
}

function bisulfite::haarz(){
	bisulfite::mecall "$@"
}

function bisulfite::mecall(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-M <maxmemory>| amount of
			-g <genome>   | path to
			-x <context>  | 2-base - default: CG
			-r <mapper>   | array of bams within array of
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory genome outdir tmpdir="${TMPDIR:-/tmp}" context=CG
	declare -n _mapper_haarz
	while getopts 'S:s:t:M:r:g:x:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			x)	context="$OPTARG";;
			r)	((++mandatory)); _mapper_haarz=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "methylation calling haarz"

	declare -n _bams_haarz=${_mapper_haarz[0]}
	local ithreads imemory instances=$((${#_bams_haarz[@]}*${#_mapper_haarz[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -T $threads)
	read -r instances imemory < <(configure::memory_by_instances -i $instances -T $threads -M "$maxmemory")

	local m f o odir
	declare -a tdirs cmd1 cmd2 cmd3 cmd4
	for m in "${_mapper_haarz[@]}"; do
		declare -n _bams_haarz=$m
		odir="$outdir/$m/haarz"
		mkdir -p "$odir"

		for f in "${_bams_haarz[@]}"; do
			o=$(basename $f)
			o=${o%.*}

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.haarz)")

			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				haarz callmethyl -t $threads -d $genome -b "$f"
			CMD
				bgzip -k -c -@ $threads > "${tdirs[-1]}/$o.vcf.gz"
			CMD

			# pipe into bgzip is much faster than using bctools sort -O z
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bcftools sort -T "${tdirs[-1]}" -m ${imemory}M "${tdirs[-1]}/$o.vcf.gz" | bgzip -k -c -@ $ithreads > "$odir/$o.vcf.gz"
			CMD
				tabix -f -p vcf "$odir/$o.vcf.gz"
			CMD

			commander::makecmd -a cmd3 -s ' ' -o "$odir/$o.$context.full.bed" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- 'CMD'
				bcftools view -H "$odir/$o.vcf.gz" |
			CMD
				perl -slane '
					next unless $F[7]=~/;CC=$c;/;
					$F[1]-- if $F[7]=~/^CS=-;/;
					print join"\t",($F[0],$F[1]-1,$F[1],(split/:/,$F[-1])[-3,-4,0])
				'
			CMD
				-- -c="$context" |
			CMD
				bedtools merge -d -1 -c 4,5,6 -o sum,sum,max |
			CMD
				perl -lane 'print join"\t",(@F[0..2],$F[3]/$F[4],$F[5])'
			CMD

			commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				awk '\$NF>=10' "$odir/$o.$context.full.bed"
			CMD
				cut -f 1-4 > "$odir/$o.$context.bed"
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		commander::runcmd -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $instances -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}

function bisulfite::methyldackel(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-g <genome>   | path to
			-x <context>  | 2-base - default: CG
			-r <mapper>   | array of bams within array of
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads genome outdir context=CG tmpdir="${TMPDIR:-/tmp}"
	declare -n _mapper_methyldackel
	while getopts 'S:s:t:m:r:g:x:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			x)	context="$OPTARG";;
			r)	((++mandatory)); _mapper_methyldackel=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "methylation calling methyldackel"

	declare -a tdirs cmd1 cmd2
	local m f o odir
	for m in "${_mapper_methyldackel[@]}"; do
		declare -n _bams_methyldackel=$m
		odir="$outdir/$m/methyldackel"
		mkdir -p "$odir"

		for f in "${_bams_methyldackel[@]}"; do
			o=$(basename $f)
			o=${o%.*}

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.methyldackel)")

			commander::makecmd -a cmd1 -s ' ' -o "$odir/$o.$context.full.bed" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- 'CMD'
				ln -sfn "/dev/stdout" "${tdirs[-1]}/$o.cytosine_report.txt";
			CMD
				MethylDackel extract
					-q 0
					--keepDupes
					--keepSingleton
					--keepDiscordant
					--ignoreNH
					-F 0
					-@ $threads
					--CHH
					--cytosine_report
					-o "${tdirs[-1]}/$o"
					"$genome" "$f" |
			CMD
				tee >(helper::pgzip -t $threads -o "$odir/$o.cytosine_report.tsv.gz") |
			CMD
				perl -slane '
					next unless $F[-1]=~/^$c/;
					next if $F[3]+$F[4]==0;
					$F[1]-- if $F[2] eq "-";
					print join"\t",($F[0],$F[1]-1,$F[1],$F[3],$F[3]+$F[4])
				'
			CMD
				-- -c="$context" |
			CMD
				bedtools merge -d -1 -c 4,5 -o sum,sum |
			CMD
				perl -lane 'print join"\t",(@F[0..2],$F[3]/$F[4],$F[4])'
			CMD

			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				awk '\$NF>=10' "$odir/$o.$context.full.bed"
			CMD
				cut -f 1-4 > "$odir/$o.$context.bed"
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c methyldackel -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function bisulfite::metilene(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-c <cmpfiles>   | array of
			-m <minimum>    | values present in control and treatment samples as fraction or absolute number - default: 0.8
			-u <upperbound> | values present in control and treatment samples as absolute number (see -m) - default: not capped
			-x <context>    | bases - default: CG
			-r <mapper>     | array of bams within array of
			-i <methdir>    | path to
			-o <outdir>     | path to
			-d <tool>       | use if subdir in methdir. name identical to subdir. parameter can be used multiple times
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads genome mecalldir outdir min=0.8 cap=999999 context=CG
	declare -a tools
	declare -n _mapper_metilene _cmpfiles_metilene
	while getopts 'S:s:t:c:m:u:x:r:i:o:d:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			c)	((++mandatory)); _cmpfiles_metilene=$OPTARG;;
			m)	min=$OPTARG;;
			u)	cap=$OPTARG;;
			x)	context=$OPTARG;;
			r)	((++mandatory)); _mapper_metilene=$OPTARG;;
			i)	((++mandatory)); mecalldir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			d)	tools+=("/$OPTARG");;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage
	[[ $tools ]] || tools=("") # for backwards compatibility

	commander::printinfo "differential methylation analyses"

	declare -a cmd1 cmd2 cmd3 mapdata tojoin

	local m f i c t odir header sample condition library replicate factors crep trep tool
	for m in "${_mapper_metilene[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"

		for f in "${_cmpfiles_metilene[@]}"; do
			mapfile -t mapdata < <(perl -F'\t' -lane 'next if exists $m{$F[1]}; $m{$F[1]}=1; print $F[1]' "$f")
			for tool in "${tools[@]}"; do
				i=0
				for c in "${mapdata[@]::${#mapdata[@]}-1}"; do
					crep=$(awk -v s=$c -v n=$min '$2==s{i=i+1}END{if(n>1){if(i<n){n=i} print n}else{printf "%0.f", i*n}}' "$f")
					[[ $crep -gt $cap ]] && crep=$cap

					for t in "${mapdata[@]:$((++i)):${#mapdata[@]}}"; do
						odir="$outdir/$m$tool/$c-vs-$t"
						mkdir -p "$odir"
						tojoin=()
						header="chr\tpos"
						trep=$(awk -v s=$t -v n=$min '$2==s{i=i+1}END{if(n>1){if(i<n){n=i} print n}else{printf "%0.f", i*n}}' "$f")
						[[ $trep -gt $cap ]] && trep=$cap

						unset sample condition library replicate factors
						while read -r sample condition library replicate factors; do
							tojoin+=("$(find -L "$mecalldir/$m$tool" -maxdepth 1 -name "$sample*.$context.bed" -print -quit | grep .)")
							header+="\t${condition}_$replicate"
						done < <(awk -v c=$c '$2==c' "$f" | sort -k4,4V && awk -v t=$t '$2==t' "$f" | sort -k4,4V)

						commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
							echo -e "$header" > "$odir/merates.bedg"
						CMD
							bedtools unionbedg -filler . -i $(printf '"%s" ' "${tojoin[@]}") | cut -f 1,3- >> "$odir/merates.bedg"
						CMD

						bisulfite::_metilene \
							-1 cmd2 \
							-2 cmd3 \
							-t $threads \
							-i "$odir/merates.bedg" \
							-a $c \
							-b $t \
							-x $crep \
							-y $trep \
							-o "$odir"
					done
				done
			done
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -c metilene -v -b -i 1 -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
	fi

	return 0
}

function bisulfite::_metilene(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-1 <cmds1>     | array of
			-2 <cmds2>     | array of
			-t <threads>   | number of
			-i <inmerates> | path to bedgraph
			-a <condition> | name of control
			-b <condition> | name of treatment
			-x <minimum>   | of -a values
			-y <minimum>   | of -b values
			-d <distance>  | for CpGs within DMR maximum - default: 300
			-n <minimum>   | number of CpGs a DMR should consists of - default: 8
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory threads merates gtf outdir gtfinfo c t crep trep distance=300 ncpg=8
	declare -n _cmds1_metilene _cmds2_metilene
	while getopts '1:2:t:i:a:b:x:y:d:n:o:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_metilene=$OPTARG;;
			2)	((++mandatory)); _cmds2_metilene=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			i)	((++mandatory)); merates="$OPTARG";;
			a)	((++mandatory)); c=$OPTARG;;
			b)	((++mandatory)); t=$OPTARG;;
			x)	crep=$OPTARG;;
			y)	trep=$OPTARG;;
			d)	distance=$OPTARG;;
			n)	ncpg=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 7 ]] && _usage

	local params=''
	# this is the default behaviour of allowing for missing values
	[[ $crep ]] && params+=" -X $crep"
	[[ $trep ]] && params+=" -Y $trep"

	local header="chr start stop q-value mean_$c-mean_$t CpGs p-value_MWU p-value_2DKS mean_$c mean_$t"
	commander::makecmd -a _cmds1_metilene -s ';' -c {COMMANDER[0]}<<- CMD
		{	sed 's/ /\t/g' <<< "$header";
			metilene
				$params
				-t $threads
				-m $ncpg
				-M $distance
				-a $c
				-b $t
				"$merates" | sort -k1,1 -k2,2n -k3,3n;
		} > "$outdir/dmr.full.tsv"
	CMD

	commander::makecmd -a _cmds2_metilene -s ';' -c {COMMANDER[0]}<<- CMD
		awk '\$4=="q-value" || \$4<=0.05' "$outdir/dmr.full.tsv" > "$outdir/dmr.tsv"
	CMD

	return 0
}

function bisulfite::join(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-c <cmpfiles> | array of
			-x <context>  | Cp* base - default: CG
			-i <methdir>  | path to
			-o <outdir>   | path to
			-d <tool>     | use if subdir in methdir. name identical to subdir. parameter can be used multiple times
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads mecalldir outdir context=CG
	declare -a tools
    declare -n _mapper_join _cmpfiles_join
	while getopts 'S:s:t:r:c:x:i:o:f:d:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_join=$OPTARG;;
			c)	((++mandatory)); _cmpfiles_join=$OPTARG;;
			x)	context=$OPTARG;;
			i)	((++mandatory)); mecalldir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			d)	tools+=("/$OPTARG");;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage
	[[ $tools ]] || tools=("") # for backwards compatibility

	commander::printinfo "joining methylation rates, zscores"

	declare -a cmd1 cmd2 cmd3 cmd4 mapdata tojoin
	declare -A visited
	local m f i c t e header meanheader cf sample condition library replicate factors

	for tool in "${tools[@]}"; do
		for m in "${_mapper_join[@]}"; do
			odir="$outdir/$m$tool"
			mkdir -p "$odir"
			visited=()
			header="chr\tpos"
			meanheader="chr\tpos"
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
							tojoin+=("$(find -L "$mecalldir/$m$tool" -maxdepth 1 -name "$sample*.$context.bed" -print -quit | grep .)")
						done < <(awk -v c=$c '$2==c' "$f" | sort -k4,4V && awk -v t=$t '$2==t' "$f" | sort -k4,4V)
					done
				done
			done

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				echo -e "$header" > "$odir/merates.bedg"
			CMD
				bedtools unionbedg -filler . -i $(printf '"%s" ' "${tojoin[@]}") | cut -f 1,3- >> "$odir/merates.bedg"
			CMD

			# allow for Rscript | head without getting SIGPIPE error. redirection to file possible via sink()
			# println = function(F, sep="\t"){ tryCatch({cat(F,sep=sep); cat("\n");}, error=function(e){quit("no")}); };

			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				echo -e "$header" > "$odir/merates.bedg.zscores";
			CMD
				Rscript - <<< '
					options(warn=-1);
					args = commandArgs(TRUE);
					conin = file(args[1], open="r", raw=T);
					conz = file(args[2], open="a");
					conmean = file(args[3], open="w");
					conmeanz = file(args[4], open="w");

					readln = function(con, sep="\t"){ unlist(strsplit(readLines(con, n=1), split=sep)); };

					header = readln(conin);
					l=length(header);
					nfo=header[1:2];
					header=header[3:l];
					meanheader = unique(header);

					writeLines(paste(c(nfo,meanheader), collapse="\t"), conmean);
					writeLines(paste(c(nfo,meanheader), collapse="\t"), conmeanz);

					while (length(F <- readln(conin))>0){
						nfo = F[1:2];
						F = as.numeric(F[3:l]);

						Z = log(F+1);
						Z = Z-mean(Z, na.rm=T);
						Z = Z/sd(Z, na.rm=T);
						Z[is.nan(Z) | is.na(Z)] = ".";
						writeLines(paste(c(nfo,Z), collapse="\t"), conz);

						means = tapply(F, header, mean, na.rm=T)[meanheader];
						means[is.nan(means) | is.na(means)] = ".";
						writeLines(paste(c(nfo,means), collapse="\t"), conmean);

						F = as.numeric(means);
						Z = log(F+1);
						Z = Z-mean(Z, na.rm=T);
						Z = Z/sd(Z, na.rm=T);
						Z[is.nan(Z) | is.na(Z)] = ".";
						writeLines(paste(c(nfo,Z), collapse="\t"), conmeanz);
					};
					close(conin);
					close(conz);
					close(conmean);
					close(conmeanz);
				'
			CMD
				<(echo -e "$meanheader"; tail -n +2 "$odir/merates.bedg") "$odir/merates.bedg.zscores" "$odir/merates.mean.bedg" "$odir/merates.mean.bedg.zscores"
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
