#! /usr/bin/env bash
# (c) Konstantin Riege

function alignment::segemehl(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} alignes pre-processed read data in fastq(.gz) format utilizing the mapping software segemehl

			-S <hardskip>   | optional. default: false
			                | [true|false] do nothing and return
			-s <softskip>   | optional. default: false
			                | [true|false] do nothing but check for files and print commands
			-5 <skip>       | optional, default: false
			                | true/false skip md5sum check, indexing respectively
			-t <threads>    | mandatory
			                | number of threads
			-a <accuracy>   | optional. default: 95
			                | 80 to 100 (%) of required matching bases per read
			-i <insertsize> | optional. default: 200000
			                | 50 to 200000+ of theoretical maximum aligned mate pair distance
			-r <mapper>     | mandatory
			                | array of array names which contain alignment paths. segemehl will be added.
			                | mapper+=(segemehl); segemehl=(/outdir/segemehl/1.bam /outdir/segemehl/2.bam ..)
			-g <genome>     | mandatory
			                | path to genome in fasta format
			-x <genomeidx>  | mandatory
			                | path to segemehl genome index
			-n <nosplit>    | optional. default: false
			                | true/false disable split read alignments. use e.g. for DNA-Seq derived data
			-o <outdir>     | mandatory
			                | path to output directory. subdirectory segemehl will be created according to array of array names (see -r)
			-1 <fastq1>     | mandatory
			                | array which contains single or first mate fastq(.gz) paths
			-2 <fastq2>     | optional
			                | array which contains mate pair fastq(.gz) paths
			-F              | optional
			                | force indexing even if md5sums match. ignored upon -5
			-P <parameter>  | optional
			                | additional segemehl parameter
			example:
			    declare -a fastq1=(/path/to/1.fq.gz /path/to/2.fq.gz ..)
			    declare -a mapper=()
			    ${FUNCNAME[1]} -5 true -t 16 -r mapper -g /path/to/genome.fa -x /path/to/genome.segemehl.idx -o /path/to/outdir -1 fastq1

			access bam paths directly:
			    printf '%s\n' "\${segemehl[@]}" # /path/to/outdir/segemehl/1.bam /path/to/outdir/segemehl/2.bam ..

			access bam paths via array of arrays:
			    echo \${mapper[-1]} # segemehl
			    declare -n _bams=\${mapper[-1]}
			    printf '%s\n' "\${_bams[@]}" # /path/to/outdir/segemehl/1.bam /path/to/outdir/segemehl/2.bam ..
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome genomeidx outdir accuracy=95 insertsize=200000 nosplitaln=false forceidx=false inparams
	declare -n _fq1_segemehl _fq2_segemehl _mapper_segemehl
	declare -g -a segemehl=()
	while getopts 'S:s:5:t:g:x:a:n:i:r:o:1:2:P:Fh' arg; do
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
			P)	inparams="$OPTARG";;
			F)	forceidx=true;;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
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
		params="$inparams"
		$nosplitaln || params+=" -S '$o.sj'" #segemehl trims suffix
		if [[ ${_fq2_segemehl[$i]} ]]; then
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				segemehl
				$params
				-A $accuracy
				-I $insertsize
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
				-A $accuracy
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

function alignment::star(){
	declare -a tdirs
	function _cleanup::alignment::star(){
		rm -rf "${tdirs[@]}"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} alignes pre-processed read data in fastq(.gz) format utilizing the mapping software STAR

			-S <hardskip>     | optional. default: false
			                  | [true|false] do nothing and return
			-s <softskip>     | optional. default: false
			                  | [true|false] do nothing but check for files and print commands
			-5 <skip>         | optional. default: false
			                  | true/false skip md5sum check, indexing respectively
			-t <threads>      | mandatory
			                  | number of threads
			-a <accuracy>     | optional. default: 95
			                  | 80 to 100 (%) of required matching bases per read
			-i <insertsize>   | optional. default: 200000
			                  | 50 to 200000+ of maximum aligned mate pair or split-read distance
			-r <mapper>       | mandatory
			                  | array of array names which contain alignment paths. star will be added.
			                  | mapper+=(star); star=(/outdir/star/1.bam /outdir/star/2.bam ..)
			-g <genome>       | mandatory
			                  | path to genome in fasta format
			-f <gtf>          | optional
			                  | path to genome annotation in gtf format. triggers creation of splice junction database during indexing
			-x <genomeidxdir> | mandatory
			                  | path to star genome index directory
			-n <nosplit>      | optional. default: false
			                  | true/false disable split read alignments. use e.g. for DNA-Seq derived data
			-o <outdir>       | mandatory
			                  | path to output directory. subdirectory star will be created according to array of array names (see -r)
			-1 <fastq1>       | mandatory
			                  | array which contains single or first mate fastq(.gz) paths
			-2 <fastq2>       | optional
			                  | array which contains mate pair fastq(.gz) paths
			-F                | optional
			                  | force indexing even if md5sums match. ignored upon -5
			-P <parameter>    | optional
			                  | additional star parameter

			example:
			    declare -a fastq1=(/path/to/1.fq.gz /path/to/2.fq.gz ..)
			    declare -a mapper=()
			    ${FUNCNAME[1]} -5 true -t 16 -r mapper -g /path/to/genome.fa -x /path/to/genome.star.idx/ -o /path/to/outdir -1 fastq1

			access bam paths directly:
			    printf '%s\n' "\${star[@]}" # /path/to/outdir/star/1.bam /path/to/outdir/star/2.bam ..

			access bam paths via array of arrays:
			    echo \${mapper[-1]} # star
			    declare -n _bams=\${mapper[-1]}
			    printf '%s\n' "\${_bams[@]}" # /path/to/outdir/star/1.bam /path/to/outdir/star/2.bam ..
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome gtf genomeidxdir outdir accuracy=95 insertsize=200000 nosplitaln=false inparams forceidx=false tmpdir="${TMPDIR:-/tmp}"
	declare -n _fq1_star _fq2_star _mapper_star
	declare -g -a star=()
	while getopts 'S:s:5:t:g:f:x:a:n:i:r:o:1:2:P:Fh' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			f)	gtf="$OPTARG";;
			x)	((++mandatory)); genomeidxdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/star"; mkdir -p "$outdir";;
			a)	accuracy=$OPTARG;;
			n)	nosplitaln=$OPTARG;;
			i)	insertsize=$OPTARG;;
			r)	((++mandatory))
				_mapper_star=$OPTARG
				_mapper_star+=(star)
			;;
			1)	((++mandatory)); _fq1_star=$OPTARG;;
			2)	_fq2_star=$OPTARG;;
			P)	inparams="$OPTARG";;
			F)	forceidx=true;;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage

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

		params="$inparams"
		helper::makecatcmd -c extractcmd -f "${_fq1_star[$i]}"
		[[ $extractcmd != "cat" ]] && params+=" --readFilesCommand '$extractcmd'"

		if [[ ${_fq2_star[$i]} ]]; then
			$nosplitaln && params+=' --alignIntronMax 1 --alignSJDBoverhangMin=999999' || params+=" --alignMatesGapMax $insertsize --alignIntronMax $insertsize --alignSJDBoverhangMin 10"
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				STAR
				--runMode alignReads
				$params
				--outSAMmapqUnique 60
				--outFilterMatchNminOverLread $(echo $accuracy | awk '{print $1/100}')
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
			#use unique score 60 instead of default 255 - necessary for gatk implemented MappingQualityAvailableReadFilter
		else
			$nosplitaln && params+=' --alignIntronMax 1 --alignSJDBoverhangMin=999999' || params+=" --alignIntronMax $insertsize --alignSJDBoverhangMin 10"
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				STAR
				--runMode alignReads
				$params
				--outSAMmapqUnique 60
				--outFilterMatchNminOverLread $(echo $accuracy | awk '{print $1/100}')
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

function alignment::bwa(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} alignes pre-processed read data in fastq(.gz) format utilizing the mapping software BWA

			-S <hardskip>  | optional. default: false
			               | [true|false] do nothing and return
			-s <softskip>  | optional. default: false
			               | [true|false] do nothing but check for files and print commands
			-5 <skip>      | optional. default: false
			               | true/false skip md5sum check, indexing respectively
			-t <threads>   | mandatory
			               | number of threads
			-a <accuracy>  | optional. default: 95
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
			-P <parameter> | optional
			               | additional bwa parameter

			example:
			    declare -a fastq1=(/path/to/1.fq.gz /path/to/2.fq.gz ..)
			    declare -a mapper=()
			    ${FUNCNAME[1]} -5 true -t 16 -r mapper -g /path/to/genome.fa -x /path/to/genome.bwa.idx/bwa -o /path/to/outdir -1 fastq1

			access bam paths directly:
			    printf '%s\n' "\${bwa[@]}" # /path/to/outdir/bwa/1.bam /path/to/outdir/bwa/2.bam ..

			access bam paths via array of arrays:
			    echo \${mapper[-1]} # bwa
			    declare -n _bams=\${mapper[-1]}
			    printf '%s\n' "\${_bams[@]}" # /path/to/outdir/bwa/1.bam /path/to/outdir/bwa/2.bam ..
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome idxprefix outdir accuracy=95 forcemem=true forceidx=false inparams
	declare -n _fq1_bwa _fq2_bwa _mapper_bwa
	declare -g -a bwa=()
	while getopts 'S:s:5:t:g:x:a:f:i:r:o:1:2:f:P:Fh' arg; do
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
			P)	inparams="$OPTARG";;
			F)	forceidx=true;;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
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
		params="$inparams"
		if $forcemem || [[ $readlength -gt 70 ]]; then
			# minOUTscore:30 @ MM/indelpenalty:4/6 -> (100-30)/5=~14% errors -> increase minOUTscore
			# 100*(1-95/100)*6 = 25 allowed penalties -> minOUTscore = 70
			# => minOUTscore = readlength − readlength*(1−accuracy/100)*5

			if [[ ${_fq2_bwa[$i]} ]]; then
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					$bwacmd mem
						$params
						-T $(echo $accuracy | awk -v l=$readlength '{printf("%.f",l-l*(1-$1/100)*5)}')
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
						-T $(echo $accuracy | awk -v l=$readlength '{printf("%.f",l-l*(1-$1/100)*5)}')
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
			if [[ ${_fq2_bwa[$i]} ]]; then
				helper::basename -f "${_fq2_bwa[$i]}" -o o2 -e e2
				o2="$outdir/$o2"

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa	aln
						$params
						-n $(echo $accuracy | awk -v l=$readlength '{printf("%.f",l*(1-$1/100))}')
						-t $threads
						"$idxprefix"
						"${_fq1_bwa[$i]}"
					> "$o1.sai"
				CMD
					bwa	aln
						$params
						-n $(echo $accuracy | awk -v l=$readlength '{printf("%.f",l*(1-$1/100))}')
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
						$params
						-n $(echo $accuracy | awk -v l=$readlength '{printf("%.f",l*(1-$1/100))}')
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

function alignment::postprocess(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} (converts sam to bam and) either filteres alignments for uniqueness (and properly aligned mate pairs), sorts by coordinate or index them

			-S <hardskip>  | optional. default: false
			               | [true|false] do nothing and return
			-s <softskip>  | optional. default: false
			               | [true|false] do nothing but check for files and print commands
			-j <job>       | mandatory
			               | [uniqify|sort|index] to be applied on alignments (see -r). index requires coordinate sorted alignment files
			-t <threads>   | mandatory
			               | number of threads
			-r <mapper>    | mandatory
			               | array of array names which contain alignment paths. will be updated by suffixes according to job (see -j)
			               | mapper=(segemehl star); [segemehl|star]=(/outdir/[segemehl|star]/1.[unique|sorted].bam /outdir/[segemehl|star]/2.[unique|sorted].bam ..);
			-o <outdir>    | mandatory
			               | path to output directory. subdirectories will be created according to array of array names (see -r)
			-P <parameter> | optional
			               | additional samtools [view|sort|index] parameter

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			    ${FUNCNAME[1]} -t 16 -r mapper -j uniqify -o /path/to/outdir

			access bam paths directly:
			    printf '%s\n' "\${segemehl[@]}" # /path/to/outdir/segemehl/1.unique.bam /path/to/outdir/segemehl/2.unique.bam ..
			    printf '%s\n' "\${star[@]}" # /path/to/outdir/star/1.unique.bam /path/to/outdir/star/2.unique.bam ..

			access bam paths via array of arrays:
			    for tool in \${mapper[@]}; do # segemehl star
			        declare -n _bams=\$tool
			        printf '%s\n' "\${_bams[@]}" # /path/to/outdir/[segemehl|star]/1.unique.bam /path/to/outdir/[segemehl|star]/2.unique.bam ..
			    done
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir job tmpdir="${TMPDIR:-/tmp}" inparams
	declare -n _mapper_process
	while getopts 'S:s:t:j:r:o:P:h' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			j)	((++mandatory)); job=${OPTARG,,*};;
			r)	((++mandatory)); _mapper_process=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage

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
						-f "${_bams_process[$i]}" \
						-o "$outbase" \
						-r newbam \
						-P "$inparams"
					_bams_process[$i]="$newbam"
				;;
				sort)
					instances=1
					alignment::_sort \
						-1 cmd1 \
						-t $threads \
						-f "${_bams_process[$i]}" \
						-o "$outbase" \
						-r newbam \
						-P "$inparams"
					_bams_process[$i]="$newbam"
				;;
				index)
					alignment::_index \
						-1 cmd1 \
						-t $ithreads \
						-f "${_bams_process[$i]}" \
						-P "$inparams"
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

function alignment::_uniqify(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} crafts command to (convert SAM to BAM and) filter alignments for uniqueness (and properly aligned mate pairs)

			-1 <cmds1>     | mandatory
			               | first array to store commands
			-2 <cmds2>     | optional. use if samtools version < 1.12
			               | second array to store commands which depend on -1
			-t <threads>   | mandatory
			               | number of threads
			-i <sam|bam>   | mandatory
			               | path to alignment in SAM or BAM format
			-o <outbase>   | mandatory
			               | path to output directory plus alignment file prefix
			-r <var>       | optional
			               | variable to store resulting alignment file path
			-P <parameter> | optional
			               | additional samtools view parameter

			example:
			    declare -a cmds1 cmds2
			    declare outfile=""
			    ${FUNCNAME[1]} -1 cmds1 -2 cmds2 -t 16 -i /path/to/alignment.[sam|bam] -o /path/to/outdir/alignment -r outfile

			create uniqified bam:
			    commander::printcmd -a cmds1 # samtools [..] /path/to/alignment.[sam|bam] > /path/to/alignment.unique.bam
			    commander::runcmd [..] -a cmds1
			    commander::runcmd [..] -a cmds2

			access uniqified bam:
			    ls "\$outfile" # /path/to/outdir/alignment.unique.bam
		EOF
		return 1
	}

	local OPTIND arg mandatory threads sambam outbase _returnfile_uniqify inparams
	declare -n _cmds1_uniqify _cmds2_uniqify
	while getopts '1:2:t:f:o:r:P:h' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_uniqify=$OPTARG;;
			2)	_cmds2_uniqify=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			f)	((++mandatory)); sambam="$OPTARG";;
			o)	((++mandatory)); outbase="$OPTARG";;
			r)	declare -n _returnfile_uniqify=$OPTARG;;
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage

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
	local params="$inparams"
	local x=$(samtools view -F 4 "$sambam" | head -10000 | cat <(samtools view -H "$sambam") - | samtools view -c -f 1)
	[[ $x -gt 0 ]] && params+='-f 2 '

	if [[ $(samtools view -F 4 "$sambam" | head -10000 | grep -cE '\s+NH:i:[0-9]+\s+' ) -eq 0 ]]; then
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

function alignment::_sort(){
	declare -a tdirs
	function _cleanup::alignment::_sort(){
		rm -rf "${tdirs[@]}"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} crafts command to sort an alignment file by coordinate

			-1 <cmds>      | mandatory
			               | array to append commands to
			-t <threads>   | mandatory
			               | number of threads
			-f <sam|bam>   | mandatory
			               | path to alignment in SAM or BAM format
			-o <outbase>   | mandatory
			               | path to output directory plus alignment file prefix
			-r <var>       | optional
			               | variable to store resulting alignment file path
			-P <parameter> | optional
			               | additional samtools view parameter

			example:
			    declare -a cmds
			    declare outfile=""
			    ${FUNCNAME[1]} -1 cmds -t 16 -i /path/to/alignment.[sam|bam] -o /path/to/outdir/alignment -r outfile

			create sorted bam:
			    commander::printcmd -a cmds # samtools [..] /path/to/alignment.[sam|bam] > /path/to/alignment.sorted.bam
			    commander::runcmd [..] -a cmds

			access sorted bam:
			    ls "\$outfile" # /path/to/outdir/alignment.sorted.bam
		EOF
		return 1
	}

	local OPTIND arg mandatory threads bam outbase _returnfile_sort tmpdir="${TMPDIR:-/tmp}" inparams
	declare -n _cmds1_sort
	while getopts '1:t:f:o:r:P:h' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_sort=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			f)	((++mandatory)); bam="$OPTARG";;
			o)	((++mandatory)); outbase="$OPTARG";;
			r)	declare -n _returnfile_sort=$OPTARG;;
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage

	_returnfile_sort="$outbase.sorted.bam"

	local params="$inparams"
	tdirs+=("$(mktemp -u -d -p "$tmpdir" cleanup.XXXXXXXXXX.samtools)")

	commander::makecmd -a _cmds1_sort -s ';' -c {COMMANDER[0]}<<- CMD
		samtools sort
			$params
			-@ $threads
			-O BAM
			-T "${tdirs[-1]}/$(basename "$outbase")"
			"$bam"
			> "$_returnfile_sort"
	CMD

	return 0
}

function alignment::_index(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} crafts command to index a coordinate sorted alignment file

			-1 <cmds>      | mandatory
			               | array to append commands to
			-t <threads>   | mandatory
			               | number of threads
			-f <bam>       | mandatory
			               | path to coordinate sorted alignment in BAM format
			-o <outbase>   | optional
			               | path to output directory plus alignment file prefix
			-r <var>       | optional
			               | variable to store resulting alignment index file path
			-P <parameter> | optional
			               | additional samtools view parameter

			example:
			    declare -a cmds
			    declare idxfile=""
			    ${FUNCNAME[1]} -1 cmds -t 16 -i /path/to/alignment.bam -r idxfile

			create bam index:
			    commander::printcmd -a cmds # samtools [..] /path/to/alignment.sorted.bam /path/to/alignment.sorted.bai
			    commander::runcmd [..] -a cmds

			access bam index:
			    ls "\$idxfile"
		EOF
		return 1
	}

	local OPTIND arg mandatory threads bam bai _returnfile_index inparams
	declare -n _cmds1_index
	while getopts '1:t:f:o:r:P:h' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_index=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			f)	((++mandatory)); bam="$OPTARG";;
			o)	bai="$OPTARG.bai";;
			r)	declare -n _returnfile_index=$OPTARG;;
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage
	[[ $bai ]] || bai="${bam%.*}.bai"

	_returnfile_index="$bai"
	local params="$inparams"

	commander::makecmd -a _cmds1_index -s ';' -c {COMMANDER[0]}<<- CMD
		samtools index
			$params
			-@ $threads
			"$bam"
			"$bai"
	CMD

	return 0
}


function alignment::tobed(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} converts alignments from BAM to BED format splitting N cigar strings and thus decoupling mate pairs

			-S <hardskip>  | optional. default: false
			               | [true|false] do nothing and return
			-s <softskip>  | optional. default: false
			               | [true|false] do nothing but check for files and print commands
			-t <threads>   | mandatory
			               | number of threads
			-r <mapper>    | mandatory
			               | array of array names which contain alignment paths. resulting, compressed BED files will be placed next to them and suffixed with bed.gz
			               | mapper=(segemehl star); [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			-P <parameter> | optional
			               | additional bedtools bamtobed parameter

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			    ${FUNCNAME[1]} -t 16 -r mapper

			access bed paths directly:
			    printf '%s.bed.gz\n' "\${segemehl[@]%.*}" # /path/to/outdir/segemehl/1.bed.gz /path/to/outdir/segemehl/2.bed.gz ..
			    printf '%s.bed.gz\n' "\${star[@]%.*}"

			access bed paths via array of arrays:
			    for tool in \${mapper[@]}; do # segemehl star
			        declare -n _bams=\$tool
			        printf '%s.bed.gz\n' "\${_bams[@]%.*}" # /path/to/outdir/segemehl/1.bed.gz /path/to/outdir/segemehl/2.bed.gz ..
			    done
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads inparams
	declare -n _mapper_tobed
	while getopts 'S:s:t:r:P:h' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_tobed=$OPTARG;;
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
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

	commander::printinfo "converting alignments to bed.gz"

	local params="$inparams"
	declare -a cmd1
	for m in "${_mapper_tobed[@]}"; do
		declare -n _bams_tobed=$m
		for f in "${_bams_process[@]}"; do
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bedtools bamtobed $params -split -i "$f"
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

function alignment::inferstrandness(){
	local tmpfile
	function _cleanup::alignment::inferstrandness(){
		rm -f "$tmpfile"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} infers library strandness from alignment files in BAM format given a reference annotation

			-S <hardskip>   | optional. default: false
			                | [true|false] do nothing and return
			-s <softskip>   | optional. default: false
			                | [true|false] do nothing but check for files and print commands (see -d)
			-d <method>     | optional. default when used with -s: ?
			                | [0|1|2] to define method instead of inference. 0 = unstranded, 1 = stranded/fr second strand or 2 = reversely stranded /fr first strand
			-t <threads>    | mandatory
			                | number of threads
			-r <mapper>     | mandatory
			                | array of array names which contain alignment paths
			                | mapper=(segemehl star); [segemehl|star]=(/outdir/[segemehl|star]/1.bam /outdir/[segemehl|star]/2.bam ..);
			-x <strandness> | mandatory
			                | associative array to store strandness information (see -d)
			                | strandness=([/outdir/[segemehl|star]/1.bam]=[0|1|2] /outdir/[segemehl|star]/2.bam]=[0|1|2] ..)
			-g <gtf>        | mandatory unless -d
			                | path to genome annotation in gtf format
			-l <level>      | optional. default: exon
			                | feature level to use for inference from mapped reads (3rd column in gtf)

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			    declare -A strandness
			    ${FUNCNAME[1]} -t 16 -r mapper -x strandness -x 2

			access strandness directly:
			    for file in "\${segemehl[@]}"; do
			        echo "\$file \${strandness[\$file]}" # /path/to/outdir/segemehl/1.bam [0|1|2] /path/to/outdir/segemehl/2.bam [0|1|2] ..
			    done
			    for file in "\${star[@]}"; do
			        echo "\$file \${strandness[\$file]}" # /path/to/outdir/star/1.bam [0|1|2] /path/to/outdir/star/2.bam [0|1|2] ..
			    done

			access strandness via array of arrays:
			    for tool in \${mapper[@]}; do # segemehl star
			        declare -n _bams=\$tool
			        for file in "\${_bams[@]}"; do
			            echo "\$file \${strandness[\$file]}" # /path/to/outdir/[segemehl|star]/1.bam [0|1|2] /path/to/outdir/[segemehl|star]/2.bam [0|1|2] ..
			        done
			    done
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads outdir gtf level="exon" default tmpdir="${TMPDIR:-/tmp}"
	declare -n _mapper_inferstrandness _strandness_inferstrandness
	while getopts 'S:s:t:r:x:g:d:l:h' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			d)	default=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_inferstrandness=$OPTARG;;
			x)	((++mandatory)); _strandness_inferstrandness=$OPTARG;;
			g)	gtf="$OPTARG";;
			l)	level="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
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

function alignment::add4stats(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} appends current paths of postprocessed alignmet files to collection data structures

			-r <mapper>   | mandatory
			              | array of array names which contain alignment paths
			              | mapper=(segemehl star); [segemehl|star]=(/outdir/[segemehl|star]/1.bam /outdir/[segemehl|star]/2.bam ..);
			-f <force>    | optional. default: false
			              | [true|false] unset collections and start from scratch

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			    ${FUNCNAME[1]} -r mapper
			    [segemehl|star]=(/path/to/[segemehl|star]/1.uniq.bam /path/to/[segemehl|star]/2.uniq.bam ..)
			    ${FUNCNAME[1]} -r mapper

			access collection
			    echo "\${segemehl0[@]}" # /path/to/segemehl/1.bam /path/to/segemehl/1.unique.bam ..
			    echo "\${segemehl1[@]}" # /path/to/segemehl/2.bam /path/to/segemehl/2.unique.bam ..
			    ..
			    echo "\${star0[@]}" # /path/to/star/1.bam /path/to/star/1.unique.bam ..
			    echo "\${star1[@]}" # /path/to/star/2.bam /path/to/star/2.unique.bam ..
			    ..

			access collection via array of arrays:
			    for tool in \${mapper[@]}; do # segemehl star
			        declare -n _bams=\$tool
			        for idx in "\${!_bams[@]}"; do
			            declare -n _collection="\$tool\$idx"
			            printf '%s\n' "\${_collection[@]}" # /path/to/outdir/[segemehl|star]/1.bam /path/to/outdir/[segemehl|star]/1.unique.bam
			        done
			    done
		EOF
		return 1
	}

	local OPTIND arg mandatory force=false
	declare -n _mapper_add4stats
	while getopts 'r:f:h' arg; do
		case $arg in
			r)	((++mandatory)); _mapper_add4stats=$OPTARG;;
			f)	force=$OPTARG;;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
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

function alignment::bamqc(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} summarizes numbers of primary, secondary and supplementary mapped reads from alignment files

			-S <hardskip> | optional. default: false
			              | [true|false] do nothing and return
			-s <softskip> | optional. default: false
			              | [true|false] do nothing but check for files and print commands
			-t <threads>  | mandatory
			              | number of threads
			-r <mapper>   | mandatory
			              | array of array names which contain alignment paths
			              | mapper=(segemehl star); [segemehl|star]=(/outdir/[segemehl|star]/1.bam /outdir/[segemehl|star]/2.bam ..);

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			    ${FUNCNAME[1]} -t 16 -r mapper

			access summaries directly:
			    ls "\${segemehl[@]/%.*/.flagstat}" # /path/to/outdir/segemehl/1.flagstat /path/to/outdir/segemehl/2.flagstat ..
			    ls "\${star[@]/%.*/.flagstat}" # /path/to/outdir/star/1.flagstat /path/to/outdir/star/2.flagstat ..

		    access summaries via array of arrays:
			    for tool in \${mapper[@]}; do # segemehl star
			        declare -n _bams=\$tool
			        ls "\${_bams[@]/%.*/.flagstat}" # /path/to/outdir/[segemehl|star]/1.flagstat /path/to/outdir/[segemehl|star]/2.flagstat ..
			    done

		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads
	declare -n _mapper_bamqc
	while getopts 'S:s:r:t:h' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_bamqc=$OPTARG;;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
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

function alignment::bulkindex(){
	declare -a tdirs
	function _cleanup::alignment::bulkindex(){
		rm -rf "${tdirs[@]}"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} index alignment files from collections generated via alignment::add4stats

			-S <hardskip>  | optional. default: false
			               | [true|false] do nothing and return
			-s <softskip>  | optional. default: false
			               | [true|false] do nothing but check for files and print commands
			-t <threads>   | mandatory
			               | number of threads
			-r <mapper>    | mandatory
			               | array of array names which contain alignment paths
			               | mapper=(segemehl star); [segemehl|star]=(/outdir/[segemehl|star]/1.bam /outdir/[segemehl|star]/2.bam ..);
			-P <parameter> | optional
			               | additional samtools index parameter

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.sorted.bam /path/to/[segemehl|star]/2.sorted.bam ..)
			    alignment::add4stats -r mapper
			    [segemehl|star]=(/path/to/[segemehl|star]/1.sorted.uniq.bam /path/to/[segemehl|star]/2.sorted.uniq.bam ..)
			    alignment::add4stats -r mapper
			    ${FUNCNAME[1]} -t 16 -r mapper
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads inparams tmpdir="${TMPDIR:-/tmp}"
	declare -n _mapper_bulkindex
	while getopts 'S:s:t:j:r:o:P:h' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_bulkindex=$OPTARG;;
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage

	local m i
	declare -a mapper_index
	for m in "${_mapper_bulkindex[@]}"; do
		declare -n _bams_bulkindex=$m
		for i in "${!_bams_bulkindex[@]}"; do
			declare -n _mi_bulkindex=$m$i
		 	mapper_index+=("$m$i")
		done
	done

	tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.samtools)")
	alignment::postprocess \
		-S false \
		-s $skip \
		-j index \
		-t $threads \
		-o "${tdirs[-1]}" \
		-r mapper_index \
		-P "$inparams"
}

function alignment::qcstats(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} generates bar plots from read count summaries of alignment files stored in collection data structures prior generated by alignment::add4stats.
			if non existing or forced, summaries will be generated first (see also alignment::bamqc).

			-S <hardskip> | optional. default: false
			              | [true|false] do nothing and return
			-s <softskip> | optional. default: false
			              | [true|false] do nothing but check for files and print commands
			-f <force>    | optional. default: false
			              | [true|false] (re-)generate summaries (see alignment::bamqc) from alignment files stored in collection data structures
			-t <threads>  | mandatory
			              | number of threads
			-r <mapper>   | mandatory
			              | array of array names which contain alignment paths
			              | mapper=(segemehl star); [segemehl|star]=(/outdir/[segemehl|star]/1.bam /outdir/[segemehl|star]/2.bam ..);
			-o <outdir>   | mandatory
			              | path to output directory. subdirectory star will be created according to array of array names (see -r)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir force=false
	declare -n _mapper_bamstats
	while getopts 'S:s:f:r:t:o:h' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			f)	force=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_bamstats=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
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
				# total = primary mapped + secondary + supplementary (primary mapped may include unampped reads if any i.e. 0 = total - secondary - supplementary)
				# primary mapped = mapped - secondary - supplementary
				[[ ! $filter ]] && filter='mapped' || filter=$(echo "${bam/\.sorted\./.}" | rev | cut -d '.' -f 2 | rev)
				a=$(grep -m 1 -F mapped "${bam%.*}.flagstat" | cut -d ' ' -f 1)
				s=$(grep -m 1 -F secondary "${bam%.*}.flagstat" | cut -d ' ' -f 1)
				c=$((a-s))
				s=$(grep -m 1 -F supplementary "${bam%.*}.flagstat" | cut -d ' ' -f 1)
				c=$((c-s))

				if [[ ! $all ]]; then
					if [[ -s "$outdir/$b.stats" ]]; then # check if there is a preprocessing fastq stats file
						all=$(tail -1 "$outdir/$b.stats" | cut -f 3)
						tail -1 "$outdir/$b.stats" > "$o"
					else
						all=$(grep -m 1 -F primary "${bam%.*}.flagstat" | cut -d ' ' -f 1) # get all primaries ie. primary mapped + unmapped reads
						if [[ $all -eq $c ]]; then
							rm -f "$o"
						else
							echo -e "$b\tinput reads\t$all" > "$o"
						fi
					fi
				fi
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
