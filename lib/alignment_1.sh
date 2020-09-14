#! /usr/bin/env bash
# (c) Konstantin Riege

alignment::segemehl() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-5 <skip>       | true/false md5sums, indexing respectively
			-t <threads>    | number of
			-a <accuracy>   | optional: 80 to 100
			-i <insertsize> | optional: 50 to 200000+
			-r <mapper>     | array of bams within array of
			-g <genome>     | path to
			-x <genomeidx>  | path to
			-n <nosplit>    | optional: true/false
			-o <outdir>     | path to
			-1 <fastq1>     | array of
			-2 <fastq2>     | array of
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome genomeidx outdir accuracy insertsize nosplitaln=false
	declare -n _fq1_segemehl _fq2_segemehl _mapper_segemehl
	declare -g -a segemehl=()
	while getopts 'S:s:5:t:g:x:a:n:i:r:o:1:2:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			x)	((++mandatory)); genomeidx="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/segemehl"; mkdir -p "$outdir" || return 1;;
			a)	accuracy=$OPTARG;;
			n)	nosplitaln=$OPTARG;;
			i)	insertsize=$OPTARG;;
			r)	((++mandatory))
				_mapper_segemehl=$OPTARG
				_mapper_segemehl+=(segemehl)
			;;
			1)	((++mandatory)); _fq1_segemehl=$OPTARG;;
			2)	_fq2_segemehl=$OPTARG;;
			*)	_usage; return 1;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage && return 1

	commander::printinfo "mapping segemehl"

	$skipmd5 && {
		commander::warn "skip checking md5 sums and genome indexing respectively"
	} || {
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")/md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		local thismd5genome thismd5segemehl
		thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
		[[ -s "$genomeidx" ]] && thismd5segemehl=$(md5sum "$genomeidx" | cut -d ' ' -f 1)
		if [[ "$thismd5genome" != "$md5genome" || ! "$thismd5segemehl" || "$thismd5segemehl" != "$md5segemehl" ]]; then
			commander::printinfo "indexing genome for segemehl"
			declare -a cmdidx
			commander::makecmd -a cmdidx -s '|' -c {COMMANDER[0]}<<- CMD
				segemehl -x "$genomeidx" -d "$genome"
			CMD
			commander::runcmd -v -b -t $threads -a cmdidx || {
				commander::printerr "$funcname failed at indexing"
				return 1
			}
			thismd5segemehl=$(md5sum "$genomeidx" | cut -d ' ' -f 1)
			sed -i "s/md5segemehl=.*/md5segemehl=$thismd5segemehl/" $genome.md5.sh
		fi
	}

	# read not properly paired - additional tag:
	# YI:i:0 (orientation)
	# YI:i:1 (insertsize)
	# YI:i:2 (orientation + insertsize)
	# YI:i:3 (chimeric)
	# for splice site detection use $o.sngl.bed
	declare -a cmd1
	local o e params
	for i in "${!_fq1_segemehl[@]}"; do
		helper::basename -f "${_fq1_segemehl[$i]}" -o o -e e
		o="$outdir/$o"
		$nosplitaln && params='' || params=" -S $o.sj" #segemehl trims suffix
		[[ $accuracy ]] && params+=" -A $accuracy"
		if [[ ${_fq2_segemehl[$i]} ]]; then
			[[ $insertsize ]] && params+=" -I $insertsize"
			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				segemehl
				$params
				-i "$genomeidx"
				-d "$genome"
				-q "${_fq1_segemehl[$i]}"
				-p "${_fq2_segemehl[$i]}"
				-t $threads
				-b
				-o "$o.bam"
			CMD
				ln -sfnr $o.sngl.bed $o.sj
			CMD
		else
			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				segemehl
				$params
				-i "$genomeidx"
				-d "$genome"
				-q "${_fq1_segemehl[$i]}"
				-t $threads
				-b
				-o "$o.bam"
			CMD
				ln -sfnr $o.sngl.bed $o.sj
			CMD
		fi
		segemehl+=("$o.bam")
	done

	$skip && {
		commander::printcmd -a cmd1
	} || {
		{	commander::runcmd -v -b -t 1 -a cmd1
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

alignment::star() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip>     | true/false return
			-s <softskip>     | true/false only print commands
			-5 <skip>         | true/false md5sums, indexing respectively
			-t <threads>      | number of
			-a <accuracy>     | optional: 80 to 100
			-i <insertsize>   | optional: 50 to 200000+
			-r <mapper>       | array of bams within array of
			-g <genome>       | path to
			-f <gtf>          | path to
			-x <genomeidxdir> | path to
			-n <nosplit>      | optional: true/false
			-o <outdir>       | path to
			-p <tmpdir>       | path to
			-1 <fastq1>       | array of
			-2 <fastq2>       | array of
			-c <opts>         | passed to star
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome gtf genomeidxdir outdir accuracy insertsize nosplitaln=false params=''
	declare -n _fq1_star _fq2_star _mapper_star
	declare -g -a star=()
	while getopts 'S:s:5:t:g:f:x:a:n:i:r:o:p:1:2:c:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			f)	gtf="$OPTARG";;
			x)	((++mandatory)); genomeidxdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/star"; mkdir -p "$outdir" || return 1;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			a)	accuracy=$OPTARG;;
			n)	nosplitaln=$OPTARG;;
			i)	insertsize=$OPTARG;;
			r)	((++mandatory))
				_mapper_star=$OPTARG
				_mapper_star+=(star)
			;;
			1)	((++mandatory)); _fq1_star=$OPTARG;;
			2)	_fq2_star=$OPTARG;;
			c)	params="$OPTARG";;
			*)	_usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1
	commander::printinfo "mapping star"

	$skipmd5 && {
		commander::warn "skip checking md5 sums and genome indexing respectively"
	} || {
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")/md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		local thismd5genome thismd5star thismd5gtf
		thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
		[[ -s "$genomeidxdir/SA" ]] && thismd5star=$(md5sum "$genomeidxdir/SA" | cut -d ' ' -f 1)
		[[ -s $gtf ]] && thismd5gtf=$(md5sum "$gtf" | cut -d ' ' -f 1)
		if [[ "$thismd5genome" != "$md5genome" || ! "$thismd5star" || "$thismd5star" != "$md5star" ]] || [[ "$thismd5gtf" && "$thismd5gtf" != "$md5gtf" ]]; then
			commander::printinfo "indexing genome for star"
			#100 = assumend usual read length
			[[ "$thismd5gtf" ]] && params+=" --sjdbGTFfile '$gtf' --sjdbOverhang 200"
			local genomesize=$(du -sb "$genome" | cut -f 1)
			params+=' --genomeSAindexNbases '$(echo $genomesize | perl -M'List::Util qw(min)' -lane 'printf("%d",min(14, log($_)/log(2)/2 - 1))')
			genomeseqs=$(grep -c '^>' "$genome")
			[[ $genomeseqs -gt 5000 ]] && params+=' --genomeChrBinNbits '$(echo "$genomesize $genomeseqs" | perl -M'List::Util qw(min)' -lane 'printf("%d",min(18, log($F[0]/$F[1])/log(2)))')

			declare -a cmdidx
			commander::makecmd -a cmdidx -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				mkdir -p "$genomeidxdir"
			CMD
				STAR
				--runMode genomeGenerate
				$params
				--runThreadN $threads
				--genomeDir "$genomeidxdir"
				--genomeFastaFiles "$genome"
				--outFileNamePrefix "$genomeidxdir/$(basename "$genome")."
			CMD
			commander::runcmd -c star -v -b -t $threads -a cmdidx || {
				commander::printerr "$funcname failed at indexing"
				return 1
			}

			thismd5star=$(md5sum "$genomeidxdir/SA" | cut -d ' ' -f 1)
			sed -i "s/md5star=.*/md5star=$thismd5star/" $genome.md5.sh
		fi
	}

	declare -a cmd1 tdirs
	local a o e extractcmd
	for i in "${!_fq1_star[@]}"; do
		helper::basename -f "${_fq1_star[$i]}" -o o -e e
		o="$outdir/$o"
		tdirs+=("$(mktemp -u -d -p "$tmpdir" cleanup.XXXXXXXXXX.star)")

		params+=' --outSAMmapqUnique 60' #use 60 instead of default 255 - necessary for gatk implemented MappingQualityAvailableReadFilter
		helper::makecatcmd -c extractcmd -f "${_fq1_star[$i]}"
		[[ $extractcmd != "cat" ]] && params+=" --readFilesCommand '$extractcmd'"
		[[ $accuracy ]] && params+=' --outFilterMismatchNoverReadLmax '$(echo $accuracy | awk '{print $1/100}')

		if [[ ${_fq2_star[$i]} ]]; then
			$nosplitaln && params+=' --alignIntronMax 1 --alignSJDBoverhangMin=999999' || {
				[[ $insertsize ]] || insertsize=200000
				params+=" --alignMatesGapMax $insertsize --alignIntronMax $insertsize --alignSJDBoverhangMin 10"
			}
			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				STAR
				--runMode alignReads
				$params
				--runThreadN $threads
				--genomeDir "$genomeidxdir"
				--outTmpDir "${tdirs[-1]}"
				--readFilesIn "${_fq1_star[$i]}" "${_fq2_star[$i]}"
				--outFileNamePrefix "$o."
				--runRNGseed 12345
				--genomeLoad NoSharedMemory
				--outSAMtype BAM Unsorted
				--outMultimapperOrder Random
				--outReadsUnmapped None
				--outSAMunmapped Within
				--twopassMode Basic
				--outSAMstrandField intronMotif
				--alignInsertionFlush Right
				--outSAMattrRGline ID:A1 SM:sample1 LB:library1 PU:unit1 PL:illumina
			CMD
				mv $o.Aligned.out.bam $o.bam
			CMD
				ln -sfnr $o.SJ.out.tab $o.sj
			CMD
		else
			$nosplitaln && params+=' --alignIntronMax 1 --alignSJDBoverhangMin=999999' || {
				[[ $insertsize ]] || insertsize=200000
				params+=" --alignIntronMax $insertsize --alignSJDBoverhangMin 10"
			}
			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				STAR
				--runMode alignReads
				$params
				--runThreadN $threads
				--genomeDir "$genomeidxdir"
				--outTmpDir "${tdirs[-1]}"
				--readFilesIn "${_fq1_star[$i]}"
				--outFileNamePrefix "$o."
				--runRNGseed 12345
				--genomeLoad NoSharedMemory
				--outSAMtype BAM Unsorted
				--outMultimapperOrder Random
				--outReadsUnmapped None
				--outSAMunmapped Within
				--twopassMode Basic
				--outSAMstrandField intronMotif
				--alignInsertionFlush Right
				--outSAMattrRGline ID:A1 SM:sample1 LB:library1 PU:unit1 PL:illumina
			CMD
				mv $o.Aligned.out.bam $o.bam
			CMD
				ln -sfnr $o.SJ.out.tab $o.sj
			CMD
		fi
		star+=("$o.bam")
	done

	$skip && {
		commander::printcmd -a cmd1
	} || {
		{	commander::runcmd -c star -v -b -t 1 -a cmd1
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}

alignment::bwa() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-5 <skip>       | true/false md5sums, indexing respectively
			-t <threads>    | number of
			-a <accuracy>   | optional: 80 to 100
			-r <mapper>     | array of bams within array of
			-g <genome>     | path to
			-x <idxprefix>  | path to
			-f <memalog>    | true/false force mem algorithm independend of readlength
			-o <outdir>     | path to
			-1 <fastq1>     | array of
			-2 <fastq2>     | array of
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome idxprefix outdir accuracy forcemem=true
	declare -n _fq1_bwa _fq2_bwa _mapper_bwa
	declare -g -a bwa=()
	while getopts 'S:s:5:t:g:x:a:f:i:r:o:1:2:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			x)	((++mandatory)); idxprefix="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/bwa"; mkdir -p "$outdir" || return 1;;
			a)	accuracy=$((100-$OPTARG));;
			f)	forcemem=$OPTARG;;
			r)	((++mandatory))
				_mapper_bwa=$OPTARG
				_mapper_bwa+=(bwa)
			;;
			1)	((++mandatory)); _fq1_bwa=$OPTARG;;
			2)	_fq2_bwa=$OPTARG;;
			*)	_usage; return 1;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage && return 1

	commander::printinfo "mapping bwa"

	$skipmd5 && {
		commander::warn "skip checking md5 sums and genome indexing respectively"
	} || {
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")/md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		local thismd5genome thismd5bwa
		thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
		[[ -s "$idxprefix.bwt" ]] && thismd5bwa=$(md5sum "$idxprefix.bwt" | cut -d ' ' -f 1)
		if [[ "$thismd5genome" != "$md5genome" || ! "$thismd5bwa" || "$thismd5bwa" != "$md5bwa" ]]; then
			commander::printinfo "indexing genome for bwa"
			declare -a cmdidx
			commander::makecmd -a cmdidx -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				mkdir -p "$(dirname "$idxprefix")"
			CMD
				bwa index -p "$idxprefix" "$genome"
			CMD
			commander::runcmd -c bwa -v -b -t $threads -a cmdidx || {
				commander::printerr "$funcname failed at indexing"
				return 1
			}
			thismd5bwa=$(md5sum "$idxprefix.bwt" | cut -d ' ' -f 1)
			sed -i "s/md5bwa=.*/md5bwa=$thismd5bwa/" $genome.md5.sh
		fi
	}

	local instances=${#_fq1_bwa[@]} ithreads
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -T $threads)

	declare -a cmd1 cmd2
	local i o1 e1 o2 e2 readlength catcmd params
	for i in "${!_fq1_bwa[@]}"; do

		helper::basename -f "${_fq1_bwa[$i]}" -o o1 -e e1
		helper::basename -f "${_fq2_bwa[$i]}" -o o2 -e e2
		o1="$outdir/$o1"
		o2="$outdir/$o2"

		helper::makecatcmd -c catcmd -f ${_fq1_bwa[$i]}
		readlength=$($catcmd $f | head -800 | awk 'NR%4==2{l+=length($0)}END{printf("%d",l/(NR/4))}')

		if $forcemem || [[ $readlength -gt 70 ]]; then
			# minOUTscore:30 @ MM/indelpenalty:4/6 -> (100-30)/4=~17% errors -> increase minOUTscore
			# solve (readlength-x)/5=accuracy -> r/5 - x/5 = a -> x/5 = r/5 - a -> x = (r/5-a)*5
			# => minOUTscore = $(( readlength - (readlength/5-accuracy)*5 ))
			[[ $accuracy ]] && params="-T $(( readlength - (readlength/5-accuracy)*5 ))"
			if [[ ${_fq2_bwa[$i]} ]]; then
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa mem
						$params
						-R '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						-a
						-t $threads
						"$idxprefix"
						${_fq1_bwa[$i]} ${_fq2_bwa[$i]}
				CMD
					samtools view -@ $threads -b > $o1.bam
				CMD
			else
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa mem
						$params
						-R '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						-a
						-t $threads
						"$idxprefix"
						${_fq1_bwa[$i]}
				CMD
					samtools view -@ $threads -b > $o1.bam
				CMD
			fi
		else
			if [[ ${_fq2_bwa[$i]} ]]; then
				commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa	aln
						-t $threads
						"$idxprefix"
						${_fq1_bwa[$i]}
					> $o1.sai
				CMD
					bwa	aln
						-t $threads
						"$idxprefix"
						${_fq2_bwa[$i]}
					> $o2.sai
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa sampe
						-r '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						"$idxprefix"
						$o1.sai $o2.sai
						${fastq1[$i]} ${fastq2[$i]}
				CMD
					samtools view -@ $ithreads -b > $o1.bam
				CMD
			else
				commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
					bwa	aln
						-t $threads
						"$idxprefix"
						${_fq1_bwa[$i]}
					> $o1.sai
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa samse
						-r '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						"$idxprefix"
						$o1.sai
						${fastq1[$i]}
				CMD
					samtools view -@ $ithreads -b > $o1.bam
				CMD
			fi
		fi
		bwa+=("$o1.bam")
	done


	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	commander::runcmd -c bwa -v -b -t 1 -a cmd1 && \
			commander::runcmd -c bwa -v -b -t $instances -a cmd2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

alignment::postprocess() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-j <job>      | [uniqify|sort|index]
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-p <tmpdir>   | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir job
	declare -n _mapper_process
	while getopts 'S:s:t:j:r:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			j) ((++mandatory)); job=${OPTARG,,*};;
			r) ((++mandatory)); _mapper_process=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	local instances ithreads m i outbase newbam
	for m in "${_mapper_process[@]}"; do
		declare -n _bams_process=$m
		((instances+=${#_bams_process[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	commander::printinfo "$job alignments"

	declare -a cmd1 cmd2 tdirs
	for m in "${_mapper_process[@]}"; do
		declare -n _bams_process=$m
		mkdir -p "$outdir/$m"
		for i in "${!_bams_process[@]}"; do
			outbase=$outdir/$m/$(basename "${_bams_process[$i]}")
			outbase="${outbase%.*}"
			case $job in
				uniqify)
					alignment::_uniqify \
						-1 cmd1 \
						-2 cmd2 \
						-t $ithreads \
						-m $m \
						-i "${_bams_process[$i]}" \
						-o "$outbase" \
						-r newbam
					_bams_process[$i]="$newbam"
				;;
				sort)
					instances=1
					tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.samtools)")
					alignment::_sort \
						-1 cmd1 \
						-t $threads \
						-i "${_bams_process[$i]}" \
						-o "$outbase" \
						-p "${tdirs[-1]}" \
						-r newbam
					_bams_process[$i]="$newbam"
				;;
				index)
					alignment::_index \
						-1 cmd1 \
						-t $ithreads \
						-i "${_bams_process[$i]}" \
				;;
				*) _usage; return 1;;
			esac
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	commander::runcmd -v -b -t $instances -a cmd1 && \
			commander::runcmd -v -b -t $instances -a cmd2
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}

alignment::_uniqify() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-t <threads>  | number of
			-m <mapper>   | name of
			-i <sam|bam>  | alignment file
			-o <outbase>  | path to
			-r <var>      | returned alignment file
		EOF
		return 0
	}

	local OPTIND arg mandatory threads sambam outbase m
	declare -n _cmds1_uniqify _cmds2_uniqify _returnfile_uniqify
	while getopts '1:2:t:i:o:r:m:' arg; do
		case $arg in
			1) ((++mandatory)); _cmds1_uniqify=$OPTARG;;
			2) ((++mandatory)); _cmds2_uniqify=$OPTARG;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); sambam="$OPTARG";;
			m) m=${OPTARG,,*};;
			o) ((++mandatory)); outbase="$OPTARG";;
			r) _returnfile_uniqify=$OPTARG; ;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	_returnfile_uniqify="$outbase.unique.bam"

	readlink -e "$sambam" | file -f - | grep -qF compressed || {
		commander::makecmd -a _cmds1_uniqify -s '|' -c {COMMANDER[0]}<<- CMD
			samtools view
				-@ $threads
				-b
				-F 4
				"$sambam"
				> "$outbase.bam"
		CMD
	}

	commander::makecmd -a _cmds2_uniqify -s '|' -c {COMMANDER[0]}<<- CMD
		samtools view
			-@ $threads
			-b
			-f 4
			"$sambam"
			> "$outbase.unmapped.bam"
	CMD

	# infer SE or PE filter
	local params=''
	[[ $(samtools view -F 4 "$sambam" | head -10000 | cat <(samtools view -H "$sambam") - | samtools view -c -f 1) -gt 0 ]] && params+='-f 2 '

	if [[ "$m" =~ bwa || "$m" =~ bowtie || "$m" =~ tophat || "$m" =~ hisat || "$m" =~ star || \
		$(samtools view -F 4 "$sambam" | head -10000 | grep -cE '\s+NH:i:[0-9]+\s+' ) -eq 0 ]]; then

		#extract uniques just by MAPQ
		[[ "$m" == "star" ]] && params+='-q 60' || params+='-q 1'
		commander::makecmd -a _cmds2_uniqify -s '&&' -c {COMMANDER[0]}<<- CMD
			samtools view
				$params
				-@ $ithreads
				-F 4
				-F 256
				-b
				"$sambam"
				> "$_returnfile_uniqify"
		CMD
	else
		# sed is faster than grep here
		commander::makecmd -a _cmds2_uniqify -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
			LC_ALL=C samtools view
				$params
				-h
				-@ $ithreads
				-F 4
				-F 256
				"$sambam"
		CMD
			sed -n '/^@/p; /\tNH:i:1\t/p'
		CMD
			samtools view
				-@ $ithreads
				-b
				> "$_returnfile_uniqify"
		CMD
	fi

	return 0
}

alignment::_sort() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-1 <cmds1>    | array of
			-t <threads>  | number of
			-i <bam>      | binary alignment file
			-o <outbase>  | path to
			-p <tmpdir>   | path to
			-r <var>      | returned alignment file
		EOF
		return 0
	}

	local OPTIND arg mandatory threads bam outbase tmpdir
	declare -n _cmds1_sort _returnfile_sort
	while getopts '1:t:i:o:p:r:' arg; do
		case $arg in
			1) ((++mandatory)); _cmds1_sort=$OPTARG;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); bam="$OPTARG";;
			o) ((++mandatory)); outbase="$OPTARG";;
			p) ((++mandatory)); tmpdir="$OPTARG";;
			r) _returnfile_sort=$OPTARG; ;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	_returnfile_sort="$outbase.sorted.bam"

	commander::makecmd -a _cmds1_sort -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
		cd "$tmpdir"
	CMD
		samtools sort
			-@ $threads
			-O BAM
			-T "$(basename "$outbase")"
			"$bam"
			> "$_returnfile_sort"
	CMD

	return 0
}

alignment::_index() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-1 <cmds1>    | array of
			-t <threads>  | number of
			-i <bam>      | sorted alignment file
		EOF
		return 0
	}

	local OPTIND arg mandatory threads bam
	declare -n _cmds1_index
	while getopts '1:t:i:' arg; do
		case $arg in
			1) ((++mandatory)); _cmds1_index=$OPTARG;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); bam="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage && return 1

	commander::makecmd -a _cmds1_index -s '|' -c {COMMANDER[0]}<<- CMD
		samtools index
			-@ $threads
			"$bam"
			"${bam%.*}.bai"
	CMD

	return 0
}

alignment::inferstrandness(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-r <mapper>     | array of sorted, indexed bams within array of
			-x <strandness> | hash per bam of
			-g <gtf>        | path to
			-p <tmpdir>     | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads outdir tmpdir gtf level="exon" featuretag="gene_id"
	declare -n _mapper_inferstrandness _strandness_inferstrandness
	while getopts 'S:s:t:r:x:g:l:f:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			r) ((++mandatory)); _mapper_inferstrandness=$OPTARG;;
			x) ((++mandatory)); _strandness_inferstrandness=$OPTARG;;
			g) ((++mandatory)); gtf="$OPTARG";;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			*) _usage; return 1;;
		esac
	done

	[[ $mandatory -lt 5 ]] && _usage && return 1

	commander::printinfo "inferring library preparation method"

	declare -a cmd1 cmd2 tfiles
	local m f
	for m in "${_mapper_inferstrandness[@]}"; do
		declare -n _bams_inferstrandness=$m
		for f in "${_bams_inferstrandness[@]}"; do
			tfiles+=("$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.bed)")

			commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				perl -lane '
					next if $F[0] =~ /^(MT|chrM)/i;
					next unless $F[2] eq "exon";
					if($F[6] eq "+"){
						++$plus;
						print join"\t",($F[0],$F[3],$F[4],$F[1],0,$F[6]);
					} else {
						++$minus;
						print join"\t",($F[0],$F[3],$F[4],$F[1],0,$F[6]);
					}
					exit if $plus>2000 && $minus>2000;
				'
			CMD
				"$gtf" > "${tfiles[-1]}"
			CMD
			# requires, sorted, indexed bam
			# 0 - unstranded
			# 1 - dUTP ++,-+ (FR stranded)
			# 2 - dUTP +-,++ (FR, reversely stranded)
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
				{ echo "$f" &&
					infer_experiment.py
					-q 0
					-s 100000
					-i "$f"
					-r "${tfiles[-1]}";
				}
			CMD
				perl -lane '
					BEGIN{
						$p=1;
					}
					$f=$_ if $.==1;
					$p=0 if /SingleEnd/i;
					$s1=$F[-1] if $F[-2]=~/^\"\d?\+\+/;
					$s2=$F[-1] if $F[-2]=~/^\"\d?\+-/;
					END{
						print $s1 > 0.7 ? "1 $f" : $s2 > 0.7 ? "2 $f" : "0 $f";
					}
				'
			CMD

			# skip case
			_strandness_inferstrandness["$f"]="?"
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	local l
			declare -a a mapdata
			commander::runcmd -v -b -t $threads -a cmd1 && \
			echo ":INFO: running commands of array cmd2" && \
			commander::printcmd -a cmd2 && \
			mapfile -t mapdata < <(commander::runcmd -c rseqc -t $threads -a cmd2)
			for l in "${mapdata[@]}"; do
				a=($l)
				_strandness_inferstrandness["${a[@]:1}"]="${a[0]}"
			done
		} || {
			rm -f "${tfiles[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -f "${tfiles[@]}"
	return 0
}

alignment::add4stats(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-r <mapper>    | array of bams within array of
		EOF
		return 0
	}

	local OPTIND arg mandatory
	declare -n _mapper_add4stats
	while getopts 'r:' arg; do
		case $arg in
			r) ((++mandatory)); _mapper_add4stats=$OPTARG;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage && return 1

	local m
	for m in "${_mapper_add4stats[@]}"; do
		declare -n _bams_add4stats=$m
		for i in "${!_bams_add4stats[@]}"; do
			declare -g -a $m$i #optional (see below), declare can be used with $var! but then without assignment
			declare -n _mi_add4stats=$m$i # non existing reference (here $m$i) will be always globally declared ("declare -g $m$i")
			_mi_add4stats+=("${_bams_add4stats[$i]}")
		done
	done
	# idea for datastructure mapper(segemehl,star)
	# 1: segemehl(bam1,bam2) -> segemehl0+=(bam1), segemehl1+=(bam2) -> segemehl0(bam1), segemehl1(bam2)
	# 2: segemehl(bam1.uniq,bam2.uniq) -> segemehl0+=(bam1.uniq), segemehl1+=(bam2.uniq) -> segemehl0(bam,bam.uniq) segemehl1(bam,bam.uniq)
	# 3: [..] -> segemehl0(bam,bam.uniq,bam.uniq.rmdup) segemehl1(bam,bam.uniq,bam.uniq.rmdup)

	return 0
}

alignment::bamstats(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of sorted, indexed bams within array of
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads outdir
	declare -n _mapper_bamstats
	while getopts 'S:s:r:t:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			r) ((++mandatory)); _mapper_bamstats=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage && return 1

	commander::printinfo "summarizing mapping stats"

	local m instances ithreads
	for m in "${_mapper_bamstats[@]}"; do
		declare -n _bams_bamstats=$m
		for i in "${!_bams_bamstats[@]}"; do
			declare -n _mi_bamstats=$m$i # reference declaration in alignment::add4stats
			((instances+=${#_mi_bamstats[@]}))
		done
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a cmd1
	for m in "${_mapper_bamstats[@]}"; do
		declare -n _bams_bamstats=$m
		for i in "${!_bams_bamstats[@]}"; do
			declare -n _mi_bamstats=$m$i # reference declaration in alignment::add4stats
			for bam in "${_mi_bamstats[@]}"; do
				[[ "$bam" =~ (fullpool|pseudopool|pseudorep|pseudoreplicate) ]] && continue
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					samtools flagstat
						-@ $ithreads
						"$bam"
						> "${bam%.*}.flagstat"
				CMD
			done
		done
	done

	$skip && {
		commander::printcmd -a cmd1
	} || {
		{	commander::runcmd -v -b -t $instances -a cmd1
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	local filter b o all a s c
	declare -a cmd2
	for m in "${_mapper_bamstats[@]}"; do
		declare -n _bams_bamstats=$m
		odir="$outdir/$m"
		mkdir -p "$odir"
		echo -e "sample\ttype\tcount" > "$odir/mapping.barplot.tsv"
		for i in "${!_bams_bamstats[@]}"; do
			[[ "${_bams_bamstats[$i]}" =~ (fullpool|pseudopool|pseudorep|pseudoreplicate) ]] && continue
			declare -n _mi_bamstats=$m$i # reference declaration in alignment::add4stats
			filter=''
			b="$(basename "${_mi_bamstats[0]}")"
			b="${b%.*}"
			o="$odir/$b.stats"
			all=''
			if [[ -s "$outdir/$b.stats" ]]; then # check if there is a preprocessing fastq stats file
				all=$(tail -1 "$outdir/$b.stats" | cut -f 3)
				tail -1 "$outdir/$b.stats" > $o
			else
				rm -f $o
			fi
			for bam in "${_mi_bamstats[@]}"; do
				[[ ! $filter ]] && filter='mapped' || filter=$(echo $bam | rev | cut -d '.' -f 2 | rev)
				a=$(grep mapped -m 1 ${bam%.*}.flagstat | cut -d ' ' -f 1)
				s=$(grep secondary -m 1 ${bam%.*}.flagstat | cut -d ' ' -f 1) # get rid of multicounts - 0 if unique reads in bam only
				c=$((a-s))
				[[ ! $all ]] && all=$c # set all to what was mapped in first file unless preprocessing fastq stats file was found
				echo -e "$b\t$filter reads\t$c" >> $o
			done
			perl -F'\t' -lane '$all=$F[2] unless $all; $F[0].=" ($all)"; $F[2]=(100*$F[2]/$all); print join"\t",@F' $o | tac | awk -F '\t' '{OFS="\t"; if(c){$NF=$NF-c} c=c+$NF; print}' | tac >> "$odir/mapping.barplot.tsv"
		done

		commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
			Rscript - <<< '
				suppressMessages(library("ggplot2"));
				suppressMessages(library("scales"));
				args <- commandArgs(TRUE);
				intsv <- args[1];
				outfile <- args[2];
				m <- read.table(intsv, header=T, sep="\t");
				l <- length(m$type)/length(unique(m$sample));
				l <- m$type[1:l];
				m$type = factor(m$type, levels=l);
				pdf(outfile);
				ggplot(m, aes(x = sample, y = count, fill = type)) +
					ggtitle("Mapping") + xlab("Sample") + ylab("Readcount in %") +
					theme_bw() + guides(fill=guide_legend(title=NULL)) +
					theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
					geom_bar(position = "fill", stat = "identity") +
					scale_y_continuous(labels = percent_format());
				graphics.off();
			'
		CMD
			"$odir/mapping.barplot.tsv"  "$odir/mapping.barplot.pdf"
		CMD
	done

	$skip && {
		commander::printcmd -a cmd2
	} || {
		{	commander::runcmd -c r -v -b -t $threads -a cmd2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
