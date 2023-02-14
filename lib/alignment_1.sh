#! /usr/bin/env bash
# (c) Konstantin Riege

alignment::segemehl() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} alignes pre-processed read data in fastq(.gz) format utilizing the mapping software segemehl

			-S <hardskip>   | optional
			                | [true|false] do nothing and return
			-s <softskip>   | optional
			                | [true|false] do nothing but check for files and print commands
			-5 <skip>       | optional
			                | true/false skip md5sum check, indexing respectively
			-t <threads>    | mandatory
			                | number of threads
			-a <accuracy>   | optional
			                | 80 to 100 (%) of required matching bases per read
			-i <insertsize> | optional
			                | 50 to 200000+
			-r <mapper>     | mandatory
			                | array of array names which contain alignment paths. segemehl will be added.
			                | mapper+=(segemehl); segemehl=(/outdir/segemehl/1.bam /outdir/segemehl/2.bam ..)
			-g <genome>     | mandatory
			                | path to genome in fasta format
			-x <genomeidx>  | mandatory
			                | path to segemehl genome index
			-n <nosplit>    | optional
			                | true/false disable split read alignments. use e.g. for DNA-Seq derived data
			-o <outdir>     | mandatory
			                | path to output directory. subdirectory segemehl will be created according to array of array names (see -r)
			-1 <fastq1>     | mandatory
			                | array which contains single or first mate fastq(.gz) paths
			-2 <fastq2>     | optional
			                | array which contains mate pair fastq(.gz) paths
			-F              | optional
			                | force indexing even if md5sums match. ignored upon -5

			example:
			    declare -a fastq1=(/path/to/1.fq.gz /path/to/2.fq.gz ..)
			    declare -a mapper=()
			    ${FUNCNAME[1]} -5 true -t 16 -r mapper -g /path/to/genome.fa -x /path/to/genome.segemehl.idx -o /path/to/outdir -1 fastq1

			access bam paths directly:
			    printf '%s\n' "\${segemehl[@]}"

			access bam paths via array of arrays:
			    declare -n _bams=\${mapper[-1]}
			    printf '%s\n' "\${_bams[@]}"
		EOF
		BASHBONE_ERROR="false"
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome genomeidx outdir accuracy insertsize nosplitaln=false forceidx=false
	declare -n _fq1_segemehl _fq2_segemehl _mapper_segemehl
	declare -g -a segemehl=()
	while getopts 'S:s:5:t:g:x:a:n:i:r:o:1:2:F' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			x)	((++mandatory)); genomeidx="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/segemehl"; mkdir -p "$outdir";;
			a)	accuracy=$OPTARG;;
			n)	nosplitaln=$OPTARG;;
			i)	insertsize=$OPTARG;;
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
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "mapping segemehl"

	if $skipmd5; then
		commander::warn "skip checking md5 sums and genome indexing respectively"
	else
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")/md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		local thismd5genome thismd5segemehl
		thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
		[[ -s "$genomeidx" ]] && thismd5segemehl=$(md5sum "$genomeidx" | cut -d ' ' -f 1)
		if $forceidx || [[ "$thismd5genome" != "$md5genome" || ! "$thismd5segemehl" || "$thismd5segemehl" != "$md5segemehl" ]]; then
			commander::printinfo "indexing genome for segemehl"
			declare -a cmdidx
			commander::makecmd -a cmdidx -s '|' -c {COMMANDER[0]}<<- CMD
				segemehl -x "$genomeidx" -d "$genome"
			CMD
			commander::runcmd -v -b -i $threads -a cmdidx
			commander::printinfo "updating md5 sums"
			thismd5segemehl=$(md5sum "$genomeidx" | cut -d ' ' -f 1)
			sed -i "s/md5segemehl=.*/md5segemehl=$thismd5segemehl/" "$genome.md5.sh"
		fi
	fi

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
		$nosplitaln && params='' || params=" -S '$o.sj'" #segemehl trims suffix
		[[ $accuracy ]] && params+=" -A $accuracy"
		if [[ ${_fq2_segemehl[$i]} ]]; then
			[[ $insertsize ]] && params+=" -I $insertsize"
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
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
				touch "$o.sngl.bed"
			CMD
				ln -sfnr "$o.sngl.bed" "$o.sj"
			CMD
		else
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				segemehl
				$params
				-i "$genomeidx"
				-d "$genome"
				-q "${_fq1_segemehl[$i]}"
				-t $threads
				-b
				-o "$o.bam"
			CMD
				touch "$o.sngl.bed"
			CMD
				ln -sfnr "$o.sngl.bed" "$o.sj"
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

alignment::star() {
	declare -a tdirs
	_cleanup::alignment::star(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} alignes pre-processed read data in fastq(.gz) format utilizing the mapping software STAR

			-S <hardskip>     | optional
			                  | [true|false] do nothing and return
			-s <softskip>     | optional
			                  | [true|false] do nothing but check for files and print commands
			-5 <skip>         | optional
			                  | true/false skip md5sum check, indexing respectively
			-t <threads>      | mandatory
			                  | number of threads
			-a <accuracy>     | optional
			                  | 80 to 100 (%) of required matching bases per read
			-i <insertsize>   | optional
			                  | 50 to 200000+
			-r <mapper>       | mandatory
			                  | array of array names which contain alignment paths. star will be added.
			                  | mapper+=(star); star=(/outdir/star/1.bam /outdir/star/2.bam ..)
			-g <genome>       | mandatory
			                  | path to genome in fasta format
			-f <gtf>          | optional
			                  | path to genome annotation in gtf format. triggers creation of splice junction database during indexing
			-x <genomeidxdir> | mandatory
			                  | path to star genome index directory
			-n <nosplit>      | optional
			                  | true/false disable split read alignments. use e.g. for DNA-Seq derived data
			-o <outdir>       | mandatory
			                  | path to output directory. subdirectory star will be created according to array of array names (see -r)
			-p <tmpdir>       | mandatory
			                  | path to temporary directory
			-1 <fastq1>       | mandatory
			                  | array which contains single or first mate fastq(.gz) paths
			-2 <fastq2>       | optional
			                  | array which contains mate pair fastq(.gz) paths
			-F                | optional
			                  | force indexing even if md5sums match. ignored upon -5

			example:
			    declare -a fastq1=(/path/to/1.fq.gz /path/to/2.fq.gz ..)
			    declare -a mapper=()
			    ${FUNCNAME[1]} -5 true -t 16 -r mapper -g /path/to/genome.fa -x /path/to/genome.star.idx/ -o /path/to/outdir -p /path/to/tmpdir -1 fastq1

			access bam paths directly:
			    printf '%s\n' "\${star[@]}"

			access bam paths via array of arrays:
			    declare -n _bams=\${mapper[-1]}
			    printf '%s\n' "\${_bams[@]}"
		EOF
		BASHBONE_ERROR="false"
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome gtf genomeidxdir outdir tmpdir accuracy insertsize nosplitaln=false inparams='' forceidx=false
	declare -n _fq1_star _fq2_star _mapper_star
	declare -g -a star=()
	while getopts 'S:s:5:t:g:f:x:a:n:i:r:o:p:1:2:c:F' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			f)	gtf="$OPTARG";;
			x)	((++mandatory)); genomeidxdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/star"; mkdir -p "$outdir";;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			a)	accuracy=$OPTARG;;
			n)	nosplitaln=$OPTARG;;
			i)	insertsize=$OPTARG;;
			r)	((++mandatory))
				_mapper_star=$OPTARG
				_mapper_star+=(star)
			;;
			1)	((++mandatory)); _fq1_star=$OPTARG;;
			2)	_fq2_star=$OPTARG;;
			c)	inparams="$OPTARG";;
			F)	forceidx=true;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage
	commander::printinfo "mapping star"

	if $skipmd5; then
		commander::warn "skip checking md5 sums and genome indexing respectively"
	else
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")/md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		declare -a cmdidx=("STAR --help | grep -m 1 -F versionGenome | sed -E 's/.+\s+(.+)/\1/'")
		local doindex=$forceidx idxparams thisidxversion idxversion=$(commander::runcmd -c star -i $threads -a cmdidx)
		[[ -s "$genomeidxdir/genomeParameters.txt" ]] && thisidxversion=$(grep -m 1 -F versionGenome "$genomeidxdir/genomeParameters.txt" | sed -E 's/.+\s+(.+)/\1/')

		if [[ "$thisidxversion" != "$idxversion" ]]; then
			doindex=true
		else
			local thismd5genome thismd5star thismd5gtf
			thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
			[[ -s "$genomeidxdir/SA" ]] && thismd5star=$(md5sum "$genomeidxdir/SA" | cut -d ' ' -f 1)
			[[ -s "$gtf" ]] && thismd5gtf=$(md5sum "$gtf" | cut -d ' ' -f 1)
			if [[ "$thismd5genome" != "$md5genome" || ! "$thismd5star" || "$thismd5star" != "$md5star" ]] || [[ "$thismd5gtf" && "$thismd5gtf" != "$md5gtf" ]]; then
				doindex=true
			fi
		fi
		if $doindex; then
			commander::printinfo "indexing genome for star"
			#100 = assumend usual read length. tools like arriba use 250 in case of longer reads
			[[ -s "$gtf" ]] && idxparams+=" --sjdbGTFfile '$gtf' --sjdbOverhang 250"
			local genomesize=$(du -sb "$genome" | cut -f 1)
			idxparams+=' --genomeSAindexNbases '$(echo $genomesize | perl -M'List::Util qw(min)' -lane 'printf("%d",min(14, log($_)/log(2)/2 - 1))')
			local genomeseqs=$(grep -c '^>' "$genome")
			[[ $genomeseqs -gt 5000 ]] && idxparams+=' --genomeChrBinNbits '$(echo "$genomesize $genomeseqs" | perl -M'List::Util qw(min)' -lane 'printf("%d",min(18, log($F[0]/$F[1])/log(2)))')

			cmdidx=()
			commander::makecmd -a cmdidx -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				mkdir -p "$genomeidxdir"
			CMD
				STAR
				--runMode genomeGenerate
				$idxparams
				--runThreadN $threads
				--genomeDir "$genomeidxdir"
				--genomeFastaFiles "$genome"
				--outFileNamePrefix "$genomeidxdir/$(basename "$genome")."
			CMD
			commander::runcmd -c star -v -b -i $threads -a cmdidx
			commander::printinfo "updating md5 sums"
			thismd5star=$(md5sum "$genomeidxdir/SA" | cut -d ' ' -f 1)
			sed -i "s/md5star=.*/md5star=$thismd5star/" "$genome.md5.sh"
		fi
	fi

	declare -a cmd1
	local params a o e extractcmd
	for i in "${!_fq1_star[@]}"; do
		helper::basename -f "${_fq1_star[$i]}" -o o -e e
		o="$outdir/$o"
		tdirs+=("$(mktemp -u -d -p "$tmpdir" cleanup.XXXXXXXXXX.star)")

		params="$inparams --outSAMmapqUnique 60" #use 60 instead of default 255 - necessary for gatk implemented MappingQualityAvailableReadFilter
		helper::makecatcmd -c extractcmd -f "${_fq1_star[$i]}"
		[[ $extractcmd != "cat" ]] && params+=" --readFilesCommand '$extractcmd'"
		#[[ $accuracy ]] && params+=' --outFilterMismatchNoverReadLmax '$(echo $accuracy | awk '{print 1-$1/100}')
		[[ $accuracy ]] && params+=' --outFilterMatchNminOverLread '$(echo $accuracy | awk '{print $1/100}')

		if [[ ${_fq2_star[$i]} ]]; then
			$nosplitaln && params+=' --alignIntronMax 1 --alignSJDBoverhangMin=999999' || {
				[[ $insertsize ]] || insertsize=200000
				params+=" --alignMatesGapMax $insertsize --alignIntronMax $insertsize --alignSJDBoverhangMin 10"
			}
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
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
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
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

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c star -v -b -i 1 -a cmd1
	fi

	return 0
}

alignment::bwa() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} alignes pre-processed read data in fastq(.gz) format utilizing the mapping software BWA

			-S <hardskip>  | optional
			               | [true|false] do nothing and return
			-s <softskip>  | optional
			               | [true|false] do nothing but check for files and print commands
			-5 <skip>      | optional
			               | true/false skip md5sum check, indexing respectively
			-t <threads>   | mandatory
			               | number of threads
			-a <accuracy>  | optional
			               | 80 to 100 (%) of required matching bases per read via adjusted bwa minOUTscore parameter
			-r <mapper>    | mandatory
			               | array of array names which contain alignment paths. bwa will be added.
			               | mapper+=(bwa); bwa=(/outdir/bwa/1.bam /outdir/bwa/2.bam ..)
			-g <genome>    | mandatory
			               | path to genome in fasta format
			-x <idxprefix> | mandatory
			               | path to bwa genome index directory plus common index files prefix
			-f <memalgo>   | optional. default: true
			               | true/false force to use mem algorithm. else determine algorithm by first 1000 reads mean length (>70 triggers mem)
			-o <outdir>    | mandatory
			               | path to output directory. subdirectory bwa will be created according to array of array names (see -r)
			-1 <fastq1>    | mandatory
			               | array which contains single or first mate fastq(.gz) paths
			-2 <fastq2>    | optional
			               | array which contains mate pair fastq(.gz) paths
			-F             | optional
			               | force indexing even if md5sums match. ignored upon -5

			example:
			    declare -a fastq1=(/path/to/1.fq.gz /path/to/2.fq.gz ..)
			    declare -a mapper=()
			    ${FUNCNAME[1]} -5 true -t 16 -r mapper -g /path/to/genome.fa -x /path/to/genome.bwa.idx/bwa -o /path/to/outdir -1 fastq1

			access bam paths directly:
			    printf '%s\n' "\${bwa[@]}"

			access bam paths via array of arrays:
			    declare -n _bams=\${mapper[-1]}
			    printf '%s\n' "\${_bams[@]}"
		EOF
		BASHBONE_ERROR="false"
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome idxprefix outdir accuracy forcemem=true forceidx=false
	declare -n _fq1_bwa _fq2_bwa _mapper_bwa
	declare -g -a bwa=()
	while getopts 'S:s:5:t:g:x:a:f:i:r:o:1:2:f:F' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			x)	((++mandatory)); idxprefix="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/bwa"; mkdir -p "$outdir";;
			a)	accuracy=$OPTARG;;
			f)	forcemem=$OPTARG;;
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
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "mapping bwa"
	declare -a cmdchk=("which bwa-mem2 &> /dev/null && echo bwa-mem2 || echo bwa")
	local bwacmd=$(commander::runcmd -c bwa -a cmdchk)

	if $skipmd5; then
		commander::warn "skip checking md5 sums and genome indexing respectively"
	else
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")/md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		local thismd5genome thismd5bwa
		thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
		[[ -s "$idxprefix.pac" ]] && thismd5bwa=$(md5sum "$idxprefix.pac" | cut -d ' ' -f 1)
		if $forceidx || [[ "$thismd5genome" != "$md5genome" || ! "$thismd5bwa" || "$thismd5bwa" != "$md5bwa" ]]; then
			commander::printinfo "indexing genome for bwa"
			declare -a cmdidx
			commander::makecmd -a cmdidx -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				mkdir -p "$(dirname "$idxprefix")"
			CMD
				$bwacmd index -p "$idxprefix" "$genome"
			CMD
			commander::runcmd -c bwa -v -b -i $threads -a cmdidx
			commander::printinfo "updating md5 sums"
			thismd5bwa=$(md5sum "$idxprefix.pac" | cut -d ' ' -f 1)
			sed -i "s/md5bwa=.*/md5bwa=$thismd5bwa/" "$genome.md5.sh"
		fi
	fi

	local instances=${#_fq1_bwa[@]} ithreads
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -T $threads)

	declare -a cmd1 cmd2
	local i o1 e1 o2 e2 readlength catcmd params
	for i in "${!_fq1_bwa[@]}"; do
		helper::basename -f "${_fq1_bwa[$i]}" -o o1 -e e1
		o1="$outdir/$o1"

		helper::makecatcmd -c catcmd -f "${_fq1_bwa[$i]}"
		readlength=$($catcmd "${_fq1_bwa[$i]}" | head -4000 | awk 'NR%4==2{l+=length($0)}END{printf("%.f",l/(NR/4))}')

		if $forcemem || [[ $readlength -gt 70 ]]; then
			# minOUTscore:30 @ MM/indelpenalty:4/6 -> (100-30)/5=~14% errors -> increase minOUTscore
			# 100*(1-95/100)*6 = 25 allowed penalties -> minOUTscore = 70
			# => minOUTscore = readlength − readlength*(1−accuracy/100)*5
			[[ $accuracy ]] && params='-T '$(echo $accuracy | awk -v l=$readlength '{printf("%.f",l-l*(1-$1/100)*5)}')
			if [[ ${_fq2_bwa[$i]} ]]; then
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					$bwacmd mem
						$params
						-R '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						-a
						-Y
						-t $threads
						"$idxprefix"
						"${_fq1_bwa[$i]}" "${_fq2_bwa[$i]}"
				CMD
					samtools view -@ $threads -b > "$o1.bam"
				CMD
			else
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					$bwacmd mem
						$params
						-R '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						-a
						-Y
						-t $threads
						"$idxprefix"
						"${_fq1_bwa[$i]}"
				CMD
					samtools view -@ $threads -b > "$o1.bam"
				CMD
			fi
		else
			[[ $accuracy ]] && params='-n '$(echo $accuracy | awk -v l=$readlength '{printf("%.f",l*(1-$1/100))}')
			if [[ ${_fq2_bwa[$i]} ]]; then
				helper::basename -f "${_fq2_bwa[$i]}" -o o2 -e e2
				o2="$outdir/$o2"

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa	aln
						$params
						-t $threads
						"$idxprefix"
						"${_fq1_bwa[$i]}"
					> "$o1.sai"
				CMD
					bwa	aln
						$params
						-t $threads
						"$idxprefix"
						"${_fq2_bwa[$i]}"
					> "$o2.sai"
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa sampe
						-r '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						"$idxprefix"
						"$o1.sai" "$o2.sai"
						"${fastq1[$i]}" "${fastq2[$i]}"
				CMD
					samtools view -@ $ithreads -b > "$o1.bam"
				CMD
			else
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					bwa	aln
						-t $threads
						"$idxprefix"
						"${_fq1_bwa[$i]}"
					> "$o1.sai"
				CMD
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa samse
						-r '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						"$idxprefix"
						"$o1.sai"
						"${_fq1_bwa[$i]}"
				CMD
					samtools view -@ $ithreads -b > "$o1.bam"
				CMD
			fi
		fi
		bwa+=("$o1.bam")
	done


	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c bwa -v -b -i 1 -a cmd1
		commander::runcmd -c bwa -v -b -i $instances -a cmd2
	fi

	return 0
}

alignment::postprocess() {
	declare -a tdirs
	_cleanup::alignment::postprocess(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} (converts sam to bam and) either filteres alignments for uniqueness (and properly aligned mate pairs), sorts or index them.

			-S <hardskip> | optional
			              | [true|false] do nothing and return
			-s <softskip> | optional
			              | [true|false] do nothing but check for files and print commands
			-j <job>      | mandatory
			              | [uniqify|sort|index] to be applied on alignments (see -r)
			-t <threads>  | mandatory
			              | number of threads
			-r <mapper>   | mandatory
			              | array of array names which contain alignment paths. resulting files will be suffixed according to job (see -j). alignment paths will be updated.
			              | mapper=(segemehl star); [segemehl|star]=(/outdir/[segemehl|star]/1.[unique|sorted].bam /outdir/[segemehl|star]/2.[unique|sorted].bam ..);
			-o <outdir>   | mandatory
			              | path to output directory. subdirectories will be created according to array of array names (see -r)
			-p <tmpdir>   | mandatory
			              | path to temporary directory

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			    ${FUNCNAME[1]} -t 16 -r mapper -j uniqify -o /path/to/outdir -p /path/to/tmpdir

			access bam paths directly:
			    printf '%s\n' "\${segemehl[@]}"
			    printf '%s\n' "\${star[@]}"

			access bam paths via array of arrays:
			    for tool in \${mapper[@]}; do
			        declare -n _bams=\$tool
			        printf '%s\n' "\${_bams[@]}"
			    done
		EOF
		BASHBONE_ERROR="false"
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir job
	declare -n _mapper_process
	while getopts 'S:s:t:j:r:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			j)	((++mandatory)); job=${OPTARG,,*};;
			r)	((++mandatory)); _mapper_process=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	local instances ithreads m i outbase newbam
	for m in "${_mapper_process[@]}"; do
		declare -n _bams_process=$m
		((instances+=${#_bams_process[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	commander::printinfo "$job alignments"

	declare -a cmd1 cmd2
	for m in "${_mapper_process[@]}"; do
		declare -n _bams_process=$m
		mkdir -p "$outdir/$m"
		for i in "${!_bams_process[@]}"; do
			outbase="$outdir/$m/$(basename "${_bams_process[$i]}")"
			outbase="${outbase%.*}"
			case $job in
				uniqify)
					alignment::_uniqify \
						-1 cmd1 \
						-2 cmd2 \
						-t $ithreads \
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
				*) _usage;;
			esac
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -v -b -i $instances -a cmd1
		commander::runcmd -v -b -i $instances -a cmd2
	fi
	return 0
}

alignment::_uniqify() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} crafts command to (convert SAM to BAM and) filter alignments for uniqueness (and properly aligned mate pairs).

			-1 <cmds1>   | mandatory
			             | first array to append commands to
			-2 <cmds2>   | mandatory
			             | second array to append commands to which depend on first array commands
			-t <threads> | mandatory
			             | number of threads
			-i <sam|bam> | mandatory
			             | path to alignment in SAM or BAM format
			-o <outbase> | mandatory
			             | path to output directory plus alignment file prefix
			-r <var>     | optional
			             | variable to store resulting alignment file path

			example:
			    declare -a cmds1 cmds2
			    declare outfile=""
			    ${FUNCNAME[1]} -1 cmds1 -2 cmds2 -t 16 -i /path/to/alignment.[sam|bam] -o /path/to/outdir/alignment -r outfile

			create uniqified bam:
			    commander::runcmd [..] -a cmds1
			    commander::runcmd [..] -a cmds2

			access uniqified bam:
			    echo "\$outfile"
		EOF
		BASHBONE_ERROR="false"
		return 1
	}

	local OPTIND arg mandatory threads sambam outbase _returnfile_uniqify
	declare -n _cmds1_uniqify _cmds2_uniqify
	while getopts '1:2:t:i:o:r:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_uniqify=$OPTARG;;
			2)	((++mandatory)); _cmds2_uniqify=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			i)	((++mandatory)); sambam="$OPTARG";;
			o)	((++mandatory)); outbase="$OPTARG";;
			r)	declare -n _returnfile_uniqify=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	_returnfile_uniqify="$outbase.unique.bam"

	readlink -e "$sambam" | file -f - | grep -qF compressed || {
		commander::makecmd -a _cmds1_uniqify -s '|' -c {COMMANDER[0]}<<- CMD
			samtools view
				-@ $threads
				-b
				"$sambam"
				> "$outbase.bam"
		CMD
	}

	# infer SE or PE filter
	local params=''
	local x=$(samtools view -F 4 "$sambam" | head -10000 | cat <(samtools view -H "$sambam") - | samtools view -c -f 1)
	[[ $x -gt 0 ]] && params+='-f 2 '

	if [[ $(samtools view -F 4 "$sambam" | head -10000 | grep -cE '\s+NH:i:[0-9]+\s+' ) -eq 0 ]]; then
		#extract uniques just by MAPQ
		# commander::makecmd -a _cmds2_uniqify -s ';' -c {COMMANDER[0]}<<- CMD
		# 	samtools view
		# 		-q 1
		# 		$params
		# 		-@ $ithreads
		# 		-F 4
		# 		-F 256
		# 		-F 2048
		# 		-b
		# 		"$sambam"
		# 		> "$_returnfile_uniqify"
		# CMD
		commander::makecmd -a _cmds2_uniqify -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
			samtools view
				-h
				-q 1
				$params
				-@ $ithreads
				-F 4
				-F 256
				-F 2048
				"$sambam"
		CMD
			sed '/^@\S\S\s/!{s/$/\tNH:i:1/}'
		CMD
			samtools view
				-@ $ithreads
				-b
				> "$_returnfile_uniqify"
		CMD
	else
		declare -a cmdchk=("samtools --version | head -1 | cut -d '.' -f 2")
		local version=$(commander::runcmd -a cmdchk)
		if [[ $version -lt 12 ]]; then
			# sed is faster than grep here
			commander::makecmd -a _cmds2_uniqify -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				samtools view
					$params
					-h
					-@ $ithreads
					-F 4
					-F 256
					-F 2048
					"$sambam"
			CMD
				sed -n '/^@\S\S\s/p; /\tNH:i:1\t/p'
			CMD
				samtools view
					-@ $ithreads
					-b
					> "$_returnfile_uniqify"
			CMD
		else
			commander::makecmd -a _cmds2_uniqify -s ';' -c {COMMANDER[0]}<<- CMD
				samtools view
					$params
					-b
					-@ $ithreads
					-F 4
					-F 256
					-F 2048
					-d NH:1
					"$sambam"
				> "$_returnfile_uniqify"
			CMD
		fi
	fi

	return 0
}

alignment::_sort() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} crafts command to sort an alignment file by coordinate.

			-1 <cmds>    | mandatory
			             | array to append commands to
			-t <threads> | mandatory
			             | number of threads
			-i <sam|bam> | mandatory
			             | path to alignment in SAM or BAM format
			-o <outbase> | mandatory
			             | path to output directory plus alignment file prefix
			-r <var>     | optional
			             | variable to store resulting alignment file path

			example:
			    declare -a cmds
			    declare outfile=""
			    ${FUNCNAME[1]} -1 cmds -t 16 -i /path/to/alignment.[sam|bam] -o /path/to/outdir/alignment -r outfile

			create sorted bam:
			    commander::runcmd [..] -a cmds

			access sorted bam:
			    echo "\$outfile"
		EOF
		BASHBONE_ERROR="false"
		return 1
	}

	local OPTIND arg mandatory threads bam outbase _returnfile_sort tmpdir
	declare -n _cmds1_sort
	while getopts '1:t:i:o:r:p:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_sort=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			i)	((++mandatory)); bam="$OPTARG";;
			o)	((++mandatory)); outbase="$OPTARG";;
			p)	((++mandatory)); tmpdir="$OPTARG";;
			r)	declare -n _returnfile_sort=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	_returnfile_sort="$outbase.sorted.bam"

	commander::makecmd -a _cmds1_sort -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
		rm -f "$tmpdir/$(basename "$outbase")"*
	CMD
		samtools sort
			-@ $threads
			-O BAM
			-T "$tmpdir/$(basename "$outbase")"
			"$bam"
			> "$_returnfile_sort"
	CMD

	return 0
}

alignment::_index() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} crafts command to index a coordinate sorted alignment file.

			-1 <cmds>    | mandatory
			             | array to append commands to
			-t <threads> | mandatory
			             | number of threads
			-i <bam>     | mandatory
			             | path to coordinate sorted alignment in BAM format
			-o <outbase> | optional
			             | path to output directory plus alignment file prefix
			-r <var>     | optional
			             | variable to store resulting alignment index file path

			example:
			    declare -a cmds
			    declare idxfile=""
			    ${FUNCNAME[1]} -1 cmds -t 16 -i /path/to/alignment.bam -r idxfile

			create bam index:
			    commander::runcmd [..] -a cmds

			access bam index:
			    echo "\$idxfile"
		EOF
		BASHBONE_ERROR="false"
		return 1
	}

	local OPTIND arg mandatory threads bam bai _returnfile_index
	declare -n _cmds1_index
	while getopts '1:t:i:o:r:' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_index=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			i)	((++mandatory)); bam="$OPTARG";;
			o)	bai="$OPTARG.bai";;
			r)	declare -n _returnfile_index=$OPTARG; ;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage
	[[ $bai ]] || bai="${bam%.*}.bai"

	_returnfile_index="$bai"

	commander::makecmd -a _cmds1_index -s ';' -c {COMMANDER[0]}<<- CMD
		samtools index
			-@ $threads
			"$bam"
			"$bai"
	CMD

	return 0
}


alignment::tobed() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} converts alignments from BAM to BED format splitting N cigar strings.

			-S <hardskip> | optional
			              | [true|false] do nothing and return
			-s <softskip> | optional
			              | [true|false] do nothing but check for files and print commands
			-t <threads>  | mandatory
			              | number of threads
			-r <mapper>   | mandatory
			              | array of array names which contain alignment paths. resulting, compressed BED files will be placed next to them and suffixed with bed.gz.
			              | mapper=(segemehl star); [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			    ${FUNCNAME[1]} -t 16 -r mapper

			access bam paths directly:
			    printf '%s.bed.gz\n' "\${segemehl[@]%.*}"
			    printf '%s.bed.gz\n' "\${star[@]%.*}"

			access bam paths via array of arrays:
			    for tool in \${mapper[@]}; do
			        declare -n _bams=\$tool
			        printf '%s.bed.gz\n' "\${_bams[@]%.*}"
			    done
		EOF
		BASHBONE_ERROR="false"
		return 1
	}

	local OPTIND arg mandatory skip=false threads
	declare -n _mapper_tobed
	while getopts 'S:s:t:r:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_tobed=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	local instances ithreads m f
	for m in "${_mapper_tobam[@]}"; do
		declare -n _bams_tobam=$m
		((instances+=${#_bams_tobam[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	commander::printinfo "convertig alignments to bed.gz"

	declare -a cmd1
	for m in "${_mapper_tobed[@]}"; do
		declare -n _bams_tobed=$m
		for f in "${_bams_process[@]}"; do
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bedtools bamtobed -split -i "$f"
			CMD
				helper::pgzip -t $ithreads -o "${f%.*}.bed.gz"
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $instances -a cmd1
	fi

	return 0
}

alignment::inferstrandness(){
	local tmpfile
	_cleanup::alignment::inferstrandness(){
		rm -f "$tmpfile"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands. use with -d
			-d <default>    | do not infer. assign 0 (unstranded), 1 (stranded/fr second strand) or 2 (reversely stranded /fr first strand)
			-t <threads>    | number of
			-r <mapper>     | array of sorted, indexed bams within array of
			-x <strandness> | hash per bam of
			-g <gtf>        | path to
			-l <level>      | feature (default: exon)
			-p <tmpdir>     | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads outdir gtf level="exon" default tmpdir
	declare -n _mapper_inferstrandness _strandness_inferstrandness
	while getopts 'S:s:t:r:x:g:d:l:p:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			d)	default=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_inferstrandness=$OPTARG;;
			x)	((++mandatory)); _strandness_inferstrandness=$OPTARG;;
			g)	gtf="$OPTARG";;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			l)	level="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage
	[[ ! $default && ! $gtf ]] && _usage

	local m f
	if [[ $default ]]; then
		commander::printinfo "assigning default library preparation method"
		for m in "${_mapper_inferstrandness[@]}"; do
			declare -n _bams_inferstrandness=$m
			for f in "${_bams_inferstrandness[@]}"; do
				_strandness_inferstrandness["$f"]=$default
			done
		done
		return 0
	else
		commander::printinfo "inferring library preparation method"
	fi

	tmpfile="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.bed)"
	declare -a cmd1
	commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		perl -slane '
			next if $F[0] eq "MT" || $F[0] eq "chrM";
			next unless $F[2] eq $lvl;
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
		-- -lvl="$level" "$gtf" > "$tmpfile"
	CMD

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $threads -a cmd1
	fi

	declare -a cmd2
	local m f
	for m in "${_mapper_inferstrandness[@]}"; do
		declare -n _bams_inferstrandness=$m
		for f in "${_bams_inferstrandness[@]}"; do
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
					-r "$tmpfile";
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

			_strandness_inferstrandness["$f"]='?' # skip case
		done
	done

	if $skip; then
		commander::printcmd -a cmd2
	else
		local l
		declare -a a mapdata
		commander::printinfo "running commands of array cmd2"
		commander::printcmd -a cmd2
		mapfile -t mapdata < <(commander::runcmd -c rseqc -i $threads -a cmd2)
		for l in "${mapdata[@]}"; do
			a=($l)
			_strandness_inferstrandness["${a[@]:1}"]="${a[0]}"
		done
	fi

	return 0
}

alignment::add4stats(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-r <mapper>   | array of bams within array of
			-f <force>    | true/false start from scratch
		EOF
		return 1
	}

	local OPTIND arg mandatory force=false
	declare -n _mapper_add4stats
	while getopts 'r:f:' arg; do
		case $arg in
			r)	((++mandatory)); _mapper_add4stats=$OPTARG;;
			f)	force=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	local m
	for m in "${_mapper_add4stats[@]}"; do
		declare -n _bams_add4stats=$m
		for i in "${!_bams_add4stats[@]}"; do
			declare -g -a $m$i #optional (see below), declare can be used with $var! but then without assignment
			declare -n _mi_add4stats=$m$i # non existing reference (here $m$i) will be always globally declared ("declare -g $m$i")
			$force && _mi_add4stats=("${_bams_add4stats[$i]}") || _mi_add4stats+=("${_bams_add4stats[$i]}")
		done
	done
	# idea for datastructure mapper(segemehl,star)
	# 1: segemehl(bam1,bam2) -> segemehl0+=(bam1), segemehl1+=(bam2) -> segemehl0(bam1), segemehl1(bam2)
	# 2: segemehl(bam1.uniq,bam2.uniq) -> segemehl0+=(bam1.uniq), segemehl1+=(bam2.uniq) -> segemehl0(bam,bam.uniq) segemehl1(bam,bam.uniq)
	# 3: [..] -> segemehl0(bam,bam.uniq,bam.uniq.rmdup) segemehl1(bam,bam.uniq,bam.uniq.rmdup)

	return 0
}

alignment::bamqc(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of

			to create mapping statistics afterwards, please run alignment::add4stats hereafter
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads
	declare -n _mapper_bamqc
	while getopts 'S:s:r:t:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_bamqc=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	commander::printinfo "counting primary and supplementary alignments"

	declare -n _bams_bamqc=${_mapper_bamqc[0]}
	local ithreads instances=$((${#_mapper_bamqc[@]}*${#_bams_bamqc[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a cmd1
	local m bam
	for m in "${_mapper_bamqc[@]}"; do
		declare -n _bams_bamqc=$m
		for bam in "${_bams_bamqc[@]}"; do
			[[ "$bam" =~ (fullpool|pseudopool|pseudorep|pseudoreplicate) ]] && continue
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools flagstat
					-@ $ithreads
					"$bam"
					> "${bam%.*}.flagstat"
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $instances -a cmd1
	fi

	return 0
}

alignment::qcstats(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-f <force>    | true/false rerun bamqc
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir force=false
	declare -n _mapper_bamstats
	while getopts 'S:s:f:r:t:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			f)	force=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_bamstats=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	local m i bam
	declare -a mapper_qcstats
	for m in "${_mapper_bamstats[@]}"; do
		declare -n _bams_bamstats=$m
		for i in "${!_bams_bamstats[@]}"; do
			declare -n _mi_bamstats=$m$i
			[[ ${#_mi_bamstats[@]} -eq 0 || "${_bams_bamstats[$i]}" =~ (fullpool|pseudopool|pseudorep|pseudoreplicate) ]] && continue
			if $force || [[ ! -s "${_mi_bamstats[0]%.*}.flagstat" ]]; then
		 		mapper_qcstats+=("$m$i")
			fi
			# rescue enables to call add4stats without bamqc after each bam processing step - instead run flagstat on all files in parallel
		done
	done
	if [[ ${#mapper_qcstats[@]} -gt 0 ]]; then
		alignment::bamqc -S false -s $skip -t $threads -r mapper_qcstats
	fi

	commander::printinfo "plotting mapping stats"

	local filter b o all a s c x odir
	declare -a cmd1
	for m in "${_mapper_bamstats[@]}"; do
		declare -n _bams_bamstats=$m
		odir="$outdir/$m"
		mkdir -p "$odir"
		echo -e "sample\ttype\tcount" > "$odir/mapping.barplot.tsv"
		for i in "${!_bams_bamstats[@]}"; do
			[[ "${_bams_bamstats[$i]}" =~ (fullpool|pseudopool|pseudorep|pseudoreplicate) ]] && continue
			declare -n _mi_bamstats=$m$i # reference declaration in alignment::add4stats

			x=${#_mi_bamstats[@]}
			[[ $x -eq 0 ]] && continue


			b="$(basename "${_mi_bamstats[0]}")"
			b="${b%.*}"
			o="$odir/$b.stats"
			unset all filter
			for bam in "${_mi_bamstats[@]}"; do
				if [[ ! $all ]]; then
					if [[ -s "$outdir/$b.stats" ]]; then # check if there is a preprocessing fastq stats file
						all=$(tail -1 "$outdir/$b.stats" | cut -f 3)
						tail -1 "$outdir/$b.stats" > "$o"
					else
						all=$(grep -m 1 -F primary "${bam%.*}.flagstat" | cut -d ' ' -f 1) # get all primaries ie. primary mapped + unmapped reads
						echo -e "$b\tinput reads\t$all" > "$o"
					fi
				fi

				# total = primary mapped + secondary + supplementary (primary mapped may include unampped reads if any i.e. 0 = total - secondary - supplementary)
				# primary mapped = mapped - secondary - supplementary
				[[ ! $filter ]] && filter='mapped' || filter=$(echo "${bam/\.sorted\./.}" | rev | cut -d '.' -f 2 | rev)
				a=$(grep -m 1 -F mapped "${bam%.*}.flagstat" | cut -d ' ' -f 1)
				s=$(grep -m 1 -F secondary "${bam%.*}.flagstat" | cut -d ' ' -f 1)
				c=$((a-s))
				s=$(grep -m 1 -F supplementary "${bam%.*}.flagstat" | cut -d ' ' -f 1)
				c=$((c-s))
				echo -e "$b\t$filter reads\t$c" >> "$o"
			done
			perl -F'\t' -lane '$all=$F[2] unless $all; $F[0].=" ($all)"; $F[2]=(100*$F[2]/$all); print join"\t",@F' $o | tac | awk -F '\t' '{OFS="\t"; if(c){$NF=$NF-c} c=c+$NF; print}' | tac >> "$odir/mapping.barplot.tsv"
		done

		commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
			Rscript - <<< '
				suppressMessages(library("ggplot2"));
				suppressMessages(library("scales"));
				args <- commandArgs(TRUE);
				intsv <- args[1];
				outfile <- args[2];
				m <- read.table(intsv, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
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

	if [[ $x -lt 2 ]]; then
		commander::warn "too few postprocessing steps applied for plotting"
	else
		if $skip; then
			commander::printcmd -a cmd1
		else
			commander::runcmd -v -b -i $threads -a cmd1
		fi
	fi

	return 0
}
