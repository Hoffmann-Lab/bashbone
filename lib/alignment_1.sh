#! /usr/bin/env bash
# (c) Konstantin Riege

function alignment::segemehl(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} alignes pre-processed read data in fastq(.gz) format utilizing the mapping software segemehl

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
			    ${FUNCNAME[-2]} -5 true -t 16 -r mapper -g /path/to/genome.fa -x /path/to/genome.segemehl.idx -o /path/to/outdir -1 fastq1

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
		[[ ! -s "$genome.md5.sh" ]] && cp "$BASHBONE_DIR/lib/md5.sh" "$genome.md5.sh"
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
			sed -i "s/md5genome=.*/md5genome=$thismd5genome/" "$genome.md5.sh"
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
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} alignes pre-processed read data in fastq(.gz) format utilizing the mapping software STAR

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
			-c <convert>      | optional
			                  | [true|false] convert genomic into transcriptomic alignments
			-F                | optional
			                  | force indexing even if md5sums match. ignored upon -5
			-P <parameter>    | optional
			                  | additional star parameter

			example:
			    declare -a fastq1=(/path/to/1.fq.gz /path/to/2.fq.gz ..)
			    declare -a mapper=()
			    ${FUNCNAME[-2]} -5 true -t 16 -r mapper -g /path/to/genome.fa -x /path/to/genome.star.idx/ -o /path/to/outdir -1 fastq1

			access bam paths directly:
			    printf '%s\n' "\${star[@]}" # /path/to/outdir/star/1.bam /path/to/outdir/star/2.bam ..

			access bam paths via array of arrays:
			    echo \${mapper[-1]} # star
			    declare -n _bams=\${mapper[-1]}
			    printf '%s\n' "\${_bams[@]}" # /path/to/outdir/star/1.bam /path/to/outdir/star/2.bam ..
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome gtf convert=false genomeidxdir outdir accuracy=95 insertsize=200000 nosplitaln=false inparams forceidx=false tmpdir="${TMPDIR:-/tmp}"
	declare -n _fq1_star _fq2_star _mapper_star
	declare -g -a star=()
	while getopts 'S:s:5:t:g:f:c:x:a:n:i:r:o:1:2:P:Fh' arg; do
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
			c)	convert="$OPTARG";;
			P)	inparams="$OPTARG";;
			F)	forceidx=true;;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage
	$convert && nosplitaln=false

	commander::printinfo "mapping star"

	if $skipmd5; then
		commander::warn "skip checking md5 sums and genome indexing respectively"
	else
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$BASHBONE_DIR/lib/md5.sh" "$genome.md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		declare -a cmdidx=("STAR --help | grep -m 1 -F versionGenome | sed -E 's/.+\s+(.+)/\1/'")
		local idxparams thisidxversion idxversion=$(commander::runcmd -c star -i $threads -a cmdidx)
		[[ -s "$genomeidxdir/genomeParameters.txt" ]] && thisidxversion=$(grep -m 1 -F versionGenome "$genomeidxdir/genomeParameters.txt" | sed -E 's/.+\s+(.+)/\1/')
		local thismd5genome thismd5star thismd5gtf
		thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
		[[ -s "$genomeidxdir/SA" ]] && thismd5star=$(md5sum "$genomeidxdir/SA" | cut -d ' ' -f 1)
		[[ -s "$gtf" ]] && thismd5gtf=$(md5sum "$gtf" | cut -d ' ' -f 1)

		if $forceidx || [[ "$thisidxversion" != "$idxversion" || "$thismd5genome" != "$md5genome" || ! "$thismd5star" || "$thismd5star" != "$md5star" ]] || [[ "$thismd5gtf" && "$thismd5gtf" != "$md5gtf" ]]; then
			commander::printinfo "indexing genome for star"
			#100 = assumend usual read length. tools like arriba use 250 in case of longer reads
			#TODO ensure --sjdbGTFfeatureExon exon --sjdbGTFtagExonParentTranscript transcript_id --sjdbGTFtagExonParentGene gene_id
			[[ -s "$gtf" ]] && idxparams+=" --sjdbGTFfile '$gtf' --sjdbOverhang 250"
			# local genomesize=$(du -sb "$(readlink -e "$genome")" | cut -f 1)
			local genomesize=$(( $(wc -c "$genome" | cut -d " " -f 1) - $(wc -l "$genome" | cut -d " " -f 1) ))

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
			sed -i "s/md5genome=.*/md5genome=$thismd5genome/" "$genome.md5.sh"
			[[ -s "$gtf" ]] && sed -i "s/md5gtf=.*/md5gtf=$thismd5gtf/" "$genome.md5.sh"
		fi
	fi

	declare -a tdirs cmd1 catcmd
	local params a o e
	for i in "${!_fq1_star[@]}"; do
		helper::basename -f "${_fq1_star[$i]}" -o o -e e
		o="$outdir/$o"
		tdirs+=("$(mktemp -u -d -p "$tmpdir" cleanup.XXXXXXXXXX.star)")

		params="$inparams"
		if $convert; then
			grep -m 1 '^sjdbGTFfile' "$genomeidxdir/genomeParameters.txt" | cut -f 2 | grep -qvFw -- -
			# for salmon or eXpress quantification. switch to IndelSoftclipSingleend for RSEM
			params+=" --quantMode TranscriptomeSAM --quantTranscriptomeBan Singleend"
		fi

		helper::makecatcmd -l -v catcmd -f "${_fq1_star[$i]}"
		[[ $extractcmd != "cat" ]] && params+=" --readFilesCommand '${catcmd[*]}'"

		if [[ ${_fq2_star[$i]} ]]; then
			$nosplitaln && params+=" --alignMatesGapMax $insertsize --alignIntronMax 1 --alignSJDBoverhangMin 999999" || params+=" --alignMatesGapMax $insertsize --alignIntronMax $insertsize --alignSJDBoverhangMin 10"

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
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
				--outSAMattributes NH HI AS nM NM MD XS
			CMD
			#use unique score 60 instead of default 255 - necessary for gatk implemented MappingQualityAvailableReadFilter
		else
			$nosplitaln && params+=' --alignIntronMax 1 --alignSJDBoverhangMin=999999' || params+=" --alignIntronMax $insertsize --alignSJDBoverhangMin 10"
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
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
				--outSAMattributes NH HI AS nM NM MD XS MC
			CMD
		fi

		if $convert; then
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				mv "$o.Aligned.out.bam" "$o.bam"
			CMD
				mv "$o.Aligned.toTranscriptome.out.bam" "$o.transcriptomic.bam"
			CMD
				ln -sfnr "$o.SJ.out.tab" "$o.sj"
			CMD
			star+=("$o.transcriptomic.bam")
		else
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				mv "$o.Aligned.out.bam" "$o.bam"
			CMD
				ln -sfnr "$o.SJ.out.tab" "$o.sj"
			CMD
			star+=("$o.bam")
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c star -v -b -i 1 -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
	fi

	return 0
}

function alignment::bwa(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} alignes pre-processed read data in fastq(.gz) format utilizing the mapping software BWA

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
			    ${FUNCNAME[-2]} -5 true -t 16 -r mapper -g /path/to/genome.fa -x /path/to/genome.bwa.idx/bwa -o /path/to/outdir -1 fastq1

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

	# use absolute path, because on some maschines sometimes bwa-mem2 aborts with error: prefix is too long
	declare -a cmdchk=("which bwa-mem2 &> /dev/null && which bwa-mem2 || echo bwa")
	local bwacmd=$(commander::runcmd -c bwa -a cmdchk)

	if $skipmd5; then
		commander::warn "skip checking md5 sums and genome indexing respectively"
	else
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$BASHBONE_DIR/lib/md5.sh" "$genome.md5.sh"
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
			sed -i "s/md5genome=.*/md5genome=$thismd5genome/" "$genome.md5.sh"
		fi
	fi

	local instances=${#_fq1_bwa[@]} ithreads
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -T $threads)

	declare -a cmd1 cmd2
	local i o1 e1 o2 e2 readlength reflength=$(cat "$genome.fai" | datamash min 2) params
	for i in "${!_fq1_bwa[@]}"; do
		helper::basename -f "${_fq1_bwa[$i]}" -o o1 -e e1
		o1="$outdir/$o1"

		readlength=$(helper::cat -f "${_fq1_bwa[$i]}" | head -4000 | awk 'NR%4==2{l+=length($0)}END{printf("%.f",l/(NR/4))}')
		# [[ $readlength -lt $reflength ]] || readlength=$reflength # option to force heavily soft-cliped read alignments on short references or to detect adapter content
		params="$inparams"
		if $forcemem || [[ $readlength -gt 70 ]]; then
			# minOUTscore:30 @ MM/indelpenalty:4/6 -> (100-30)/5=~14% errors -> increase minOUTscore
			# 100*(1-95/100)*6 = 25 allowed penalties -> minOUTscore = 70
			# => minOUTscore = readlength − readlength*(1−accuracy/100)*5
			# round: -T $(echo $accuracy | awk -v l=$readlength '{printf("%.f",l-l*(1-$1/100)*5)}')
			# ceil: -T $(echo $accuracy | awk -v l=$readlength '{print l-sprintf("%.0d",(1-$1/100)*l+1)*5}')

			if [[ ${_fq2_bwa[$i]} ]]; then
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					$bwacmd mem
						$params
						-T $(echo $accuracy | awk -v l=$readlength '{print l-sprintf("%.0d",(1-$1/100)*l+1)*5}')
						-R '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						-a
						-Y
						-t $threads
						"$idxprefix"
						"${_fq1_bwa[$i]}" "${_fq2_bwa[$i]}"
				CMD
					samtools view --no-PG -@ $threads -b -o "$o1.bam"
				CMD
			else
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					$bwacmd mem
						$params
						-T $(echo $accuracy | awk -v l=$readlength '{print l-sprintf("%.0d",(1-$1/100)*l+1)*5}')
						-R '@RG\tID:A1\tSM:sample1\tLB:library1\tPU:unit1\tPL:illumina'
						-a
						-Y
						-t $threads
						"$idxprefix"
						"${_fq1_bwa[$i]}"
				CMD
					samtools view --no-PG -@ $threads -b -o "$o1.bam"
				CMD
			fi
		else
			if [[ ${_fq2_bwa[$i]} ]]; then
				helper::basename -f "${_fq2_bwa[$i]}" -o o2 -e e2
				o2="$outdir/$o2"

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bwa aln
						$params
						-n $(echo $accuracy | awk -v l=$readlength '{printf("%.f",l*(1-$1/100))}')
						-t $threads
						"$idxprefix"
						"${_fq1_bwa[$i]}"
					> "$o1.sai"
				CMD
					bwa aln
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
					samtools view --no-PG -@ $ithreads -b -o "$o1.bam"
				CMD
			else
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					bwa aln
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
					samtools view --no-PG -@ $ithreads -b -o "$o1.bam"
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
			${FUNCNAME[-2]} (converts sam to bam and) either filteres alignments for uniqueness (and properly aligned mate pairs), sorts by coordinate or index them

			-S <hardskip>    | optional. default: false
			                 | [true|false] do nothing and return
			-s <softskip>    | optional. default: false
			                 | [true|false] do nothing but check for files and print commands
			-j <job>         | mandatory
			                 | [uniqify|blacklist|sizeselect|sort|collate|index] to be applied on alignments (see -r). index requires coordinate sorted alignment files
			-f <path/string> | bedfile of regions or reference/chromosome name to remove alignments from (-j blacklist) or a range of insert sizes to keep (-j sizeselect). default: 0:1000
			-t <threads>     | mandatory
			                 | number of threads
			-r <mapper>      | mandatory
			                 | array of array names which contain alignment paths. will be updated by suffixes according to job (see -j)
			                 | mapper=(segemehl star); [segemehl|star]=(/outdir/[segemehl|star]/1.[unique|sorted].bam /outdir/[segemehl|star]/2.[unique|sorted].bam ..);
			-o <outdir>      | mandatory
			                 | path to output directory. subdirectories will be created according to array of array names (see -r)
			-M <maxmemory>   | optional. default: all available
			                 | amount of memory to allocate
			-P <parameter>   | optional
			                 | additional samtools [view|sort|index] parameter

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			    ${FUNCNAME[-2]} -t 16 -r mapper -j uniqify -o /path/to/outdir

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

	local OPTIND arg mandatory skip=false threads outdir job tmpdir="${TMPDIR:-/tmp}" inparams misc maxmemory
	declare -n _mapper_process
	while getopts 'S:s:t:j:f:r:o:M:P:h' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			j)	((++mandatory)); job=${OPTARG,,*};;
			f)	misc="$OPTARG";;
			M)	maxmemory=$OPTARG;;
			r)	((++mandatory)); _mapper_process=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage
	[[ "$job" == "exclude" && ! $misc ]] && BASHBONE_ERROR="blacklist file missing" && _usage

	declare -n _bams_process="${_mapper_process[0]}"
	local ithreads instances=$((${#_mapper_process[@]}*${#_bams_process[@]})) memory
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)
	read -r instances imemory < <(configure::memory_by_instances -i $instances -T $threads -M "$maxmemory")

	commander::printinfo "$job alignments"

	local m i outbase newbam params
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
				blacklist)
					instances=1 # use ithreads and imemory otherwise
					alignment::_blacklist \
						-1 cmd1 \
						-2 cmd2 \
						-t $threads \
						-f "${_bams_process[$i]}" \
						-b "$misc" \
						-o "$outbase" \
						-m "$maxmemory" \
						-p "$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.samtools)" \
						-r newbam \
						-P "$inparams"
					_bams_process[$i]="$newbam"
				;;
				sizeselect)
					params="-c picard"
					alignment::_sizeselect \
						-1 cmd1 \
						-2 cmd2 \
						-t $ithreads \
						-f "${_bams_process[$i]}" \
						-m $(cut -d ':' -f 1 <<< ${misc:-"0:1000"}) \
						-x $(cut -d ':' -f 2 <<< ${misc:-"0:1000"}) \
						-o "$outbase" \
						-p "$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.picard)" \
						-r newbam \
						-P "$inparams"
					_bams_process[$i]="$newbam"
				;;
				sort|collate)
					instances=1
					alignment::_$job \
						-1 cmd1 \
						-t $threads \
						-f "${_bams_process[$i]}" \
						-o "$outbase" \
						-p "$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.samtools)" \
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
		commander::runcmd $params -v -b -i $instances -a cmd1
		commander::runcmd -v -b -i $instances -a cmd2
	fi
	return 0
}

function alignment::_uniqify(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} crafts command to (convert SAM to BAM and) filter alignments for uniqueness (and properly aligned mate pairs)

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
			    ${FUNCNAME[-2]} -1 cmds1 -2 cmds2 -t 16 -i /path/to/alignment.[sam|bam] -o /path/to/outdir/alignment -r outfile

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

	readlink -e "$sambam" | file -b --mime-type -f - | grep -qF -e 'gzip' || {
		commander::makecmd -a _cmds1_uniqify -s '|' -c {COMMANDER[0]}<<- CMD
			samtools view
				--no-PG
				-@ $threads
				-b
				-o "$outbase.bam"
				"$sambam"
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
				-@ $threads
				-F 4
				-F 256
				-F 2048
				"$sambam"
		CMD
			sed '/^@\S\S\s/!{s/$/\tNH:i:1/}'
		CMD
			samtools view
				--no-PG
				-@ $threads
				-b
				-o "$_returnfile_uniqify"
		CMD
	else
		declare -a cmdchk=("samtools --version | head -1 | cut -d '.' -f 2")
		local version=$(commander::runcmd -a cmdchk)
		if [[ $version -lt 10 ]]; then
			# sed is faster than grep here
			commander::makecmd -a _cmds2_uniqify -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				samtools view
					$params
					-h
					-@ $threads
					-F 4
					-F 256
					-F 2048
					"$sambam"
			CMD
				sed -n '/^@\S\S\s/p; /\tNH:i:1\t/p'
			CMD
				samtools view
					--no-PG
					-@ $threads
					-b
					-o "$_returnfile_uniqify"
			CMD
		else
			commander::makecmd -a _cmds2_uniqify -s ';' -c {COMMANDER[0]}<<- CMD
				samtools view
					$params
					-b
					-@ $threads
					-F 4
					-F 256
					-F 2048
					-d NH:1
					-o "$_returnfile_uniqify"
					"$sambam"
			CMD
		fi
	fi

	return 0
}

function alignment::_blacklist(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} crafts command to filter out alignments in blacklisted regions

			-1 <cmds1>     | mandatory
			               | first array to store commands
			-2 <cmds2>     | mandatory
			               | second array to store commands
			-t <threads>   | mandatory
			               | number of threads
			-b <blacklist> | mandatory
			               | bed file or reference/chromosome name
			-f <sam|bam>   | mandatory
			               | path to alignment in SAM or BAM format
			-o <outbase>   | mandatory
			               | path to output directory plus alignment file prefix
			-p <tmpdir>    | mandatory
			               | path to temporary directory
			-m <memory>    | mandatory
			               | amount of memory to allocate
			-r <var>       | optional
			               | variable to store resulting alignment file path
			-P <parameter> | optional
			               | additional samtools view parameter

			example:
			    declare -a cmds1 cmds2
			    declare outfile=""
			    ${FUNCNAME[-2]} -1 cmds1 -2 cmds2 -t 16 -b /path/to/blacklist.bed -f /path/to/alignment.[sam|bam] -o /path/to/outdir/alignment -p /path/to/tmpdir -r outfile

			create blacklisted bam:
			    commander::printcmd -a cmds1 # samtools [..] /path/to/alignment.[sam|bam] > /path/to/alignment.blacklisted.bam
			    commander::runcmd [..] -a cmds1
			    commander::runcmd [..] -a cmds2

			access blacklisted bam:
			    ls "\$outfile" # /path/to/outdir/alignment.blacklisted.bam
		EOF
		return 1
	}

	local OPTIND arg mandatory threads bam outbase _returnfile_blacklist blacklist inparams tmpdir memory
	declare -n _cmds1_blacklist _cmds2_blacklist
	while getopts '1:2:t:b:f:o:p:r:m:P:h' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_blacklist=$OPTARG;;
			2)	((++mandatory)); _cmds2_blacklist=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			b)	((++mandatory)); blacklist="$OPTARG";;
			f)	((++mandatory)); bam="$OPTARG";;
			o)	((++mandatory)); outbase="$OPTARG";;
			m)	((++mandatory)); memory=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG";; # needs to be given here, otherwise a mktemp dir will be deleted upon returning of this function
			r)	declare -n _returnfile_blacklist=$OPTARG;;
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 8 ]] && _usage

	_returnfile_blacklist="$outbase.blacklisted.bam"

	if [[ -s "$blacklist" ]]; then

		# or bedtools complement
		commander::makecmd -a _cmds1_blacklist -s ';' -c {COMMANDER[0]}<<- CMD
			bedtools subtract
				-a <(samtools view -H "$bam" | awk '/^@SQ/{\$1=""; print}' | cut -d ':' -f 2,3 | sed -r 's/(\S+)\s+LN:(.+)/\1\t0\t\2/' | helper::sort -k1,1 -k2,2n -k3,3n -t $threads -M $((memory/2)))
				-b <(helper::sort -k1,1 -k2,2n -k3,3n -t $threads -f "$blacklist" -M $((memory/2)))
			> "$tmpdir/whitelist.bed"
		CMD

		# infer SE or PE
		local params params2 x=$(samtools view -F 4 "$bam" | head -10000 | cat <(samtools view -H "$bam") - | samtools view -c -f 1)
		samtools view "$bam" $(samtools view -H "$bam" | grep -m 1 -F @SQ | cut -f 2 | cut -d : -f 2) 2> /dev/null | head -1 &> /dev/null && params="$inparams -M" || params="$inparams"

		if [[ $x -gt 0 ]]; then
			# to remove both mates, needs workaraound, because
			# filter for a region only allows to keep second mate and set either to
			# 	- unmapped by -p switch (cigar to * from version 1.16) -> does not work for multi-region iterator (-M)
			# 	- or keep it by -P switch. from v1.14 iterator initialization bug still existing in v1.17
			# even if second mate is not kept, all flags, including proper paired flag kept in first mate
			# => solution1: mark as unmapped, fixmate and filter by proper pair
			# ==> solution2: fast remove via -M, fixmate and filter by pair
			# attention: disable FR check by fixmate to not loose proper pair flag for kept pairs
			commander::makecmd -a _cmds2_blacklist -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD
				samtools view
					$params
					-@ $threads
					-u
					-L "$tmpdir/whitelist.bed"
					"$bam"
			CMD
				samtools sort
					-@ $threads
					-u
					-n
					-T "$tmpdir/nsrt.$(basename "$outbase")"
			CMD
				samtools fixmate
					-@ $threads
					-p
					-u
					-r
					- -
			CMD
				samtools view
					-@ $threads
					-u
					-f 1
			CMD
				samtools sort
					-@ $threads
					-O BAM
					-T "$tmpdir/psrt.$(basename "$outbase")"
					--write-index
					-o "$_returnfile_blacklist##idx##${_returnfile_blacklist%.*}.bai"
			CMD
		else
			commander::makecmd -a _cmds2_blacklist -s ';' -c {COMMANDER[0]}<<- CMD
				samtools view
					$params
					-@ $threads
					-b
					-L "$tmpdir/whitelist.bed"
					-o "$_returnfile_blacklist"
					"$bam"
			CMD
		fi
	else
		commander::makecmd -a _cmds2_blacklist -s ';' -c {COMMANDER[0]}<<- CMD
			samtools view
				$inparams
				-@ $threads
				-b
				-e 'rname!="$blacklist"'
				-o "$_returnfile_blacklist"
				"$bam"
		CMD
	fi

	return 0
}

function alignment::_sizeselect(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} crafts command to filter properly aligned mate pairs by fragment size

			-1 <cmds1>     | mandatory
			               | first array to store commands
			-2 <cmds2>     | mandatory
			               | second array to store commands
			-t <threads>   | mandatory
			               | number of threads
			-m <minsize>   | mandatory
			               | minimum fragment size
			-x <maxsize>   | mandatory
			               | maximum fragment size
			-f <sam|bam>   | mandatory
			               | path to alignment in SAM or BAM format
			-o <outbase>   | mandatory
			               | path to output directory plus alignment file prefix
			-p <tmpdir>    | mandatory
			               | path to temporary directory
			-r <var>       | optional
			               | variable to store resulting alignment file path
			-P <parameter> | optional
			               | additional samtools view parameter

			example:
			    declare -a cmds1 cmds2
			    declare outfile=""
			    ${FUNCNAME[-2]} -1 cmds1 -2 cmds2 -t 16 -m 50 -x 1000 -f /path/to/alignment.[sam|bam] -o /path/to/outdir/alignment -p /path/to/tmpdir -r outfile

			create sizeselected bam:
			    commander::printcmd -a cmds1 # samtools [..] /path/to/alignment.[sam|bam] > /path/to/alignment.sizeselected.bam
			    commander::runcmd [..] -a cmds1
			    commander::runcmd [..] -a cmds2

			access sizeselected bam:
			    ls "\$outfile" # /path/to/outdir/alignment.sizeselected.bam
		EOF
		return 1
	}

	local OPTIND arg mandatory threads bam outbase _returnfile_sizeselect minsize maxsize inparams tmpdir
	declare -n _cmds1_sizeselect _cmds2_sizeselect
	while getopts '1:2:t:m:x:f:o:p:r:P:h' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_sizeselect=$OPTARG;;
			2)	((++mandatory)); _cmds2_sizeselect=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); minsize=$OPTARG;;
			x)	((++mandatory)); maxsize=$OPTARG;;
			f)	((++mandatory)); bam="$OPTARG";;
			o)	((++mandatory)); outbase="$OPTARG";;
			p)	((++mandatory)); tmpdir="$OPTARG";; # needs to be given here, otherwise a mktemp dir will be deleted upon returning of this function
			r)	declare -n _returnfile_sizeselect=$OPTARG;;
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 8 ]] && _usage

	_returnfile_sizeselect="$outbase.sizeselected.bam"

	# infer SE or PE
	local params="$inparams"
	local x=$(samtools view -F 4 "$bam" | head -10000 | cat <(samtools view -H "$bam") - | samtools view -c -f 1)
	BASHBONE_ERROR="no paird-end data given"
	[[ $x -gt 0 ]]
	BASHBONE_ERROR="fragment size range $minsize > $maxsize"
	[[ $minsize -lt $maxsize ]]
	unset BASHBONE_ERROR

	commander::makecmd -a _cmds1_sizeselect -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
		MALLOC_ARENA_MAX=4 picard
			-Xmx250m
			-XX:ParallelGCThreads=1
			-XX:ConcGCThreads=1
			-Djava.io.tmpdir="$tmpdir"
			CollectInsertSizeMetrics
			I="$bam"
			O="$outbase.sizemetrics.txt"
			H="$outbase.sizemetrics.pdf"
			MW=800
			INCLUDE_DUPLICATES=true
			TMP_DIR="$tmpdir"
			ASSUME_SORTED=true
			VALIDATION_STRINGENCY=SILENT
			VERBOSITY=WARNING
	CMD
		awk -v OFS='\t' '!/^\S*$/{if(p){print \$1,\$2}if(/insert_size/){p=1}}' "$outbase.sizemetrics.txt" > "$outbase.sizemetrics.tsv"
	CMD

	# declare -a cmdchk=("samtools --version | head -1 | cut -d '.' -f 2")
	# local version=$(commander::runcmd -a cmdchk)
	# if [[ $version -lt 12 ]]; then
	# 	commander::makecmd -a _cmds2_sizeselect -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
	# 		samtools view
	# 			$params
	# 			-h
	# 			-f 2
	# 			-@ $threads
	# 			"$bam"
	# 	CMD
	# 		awk -F '\t' -v m=$minsize -v x=$maxsize '/^@\S\S\s/ || (\$9 >= 0 && \$9 >= m && \$9 <= x) || (\$9 < 0 && \$9 <= -m && \$9 >= -x)'
	# 	CMD
	# 		samtools view
	# 			-@ $threads
	# 			-b
	# 		> "$_returnfile_sizeselect"
	# 	CMD
	# else
		commander::makecmd -a _cmds2_sizeselect -s ';' -c {COMMANDER[0]}<<- CMD
			samtools view
				$params
				-b
				-@ $threads
				-e "(tlen >= 0 && tlen >= $minsize && tlen <= $maxsize) || (tlen < 0 && tlen <= -$minsize && tlen >= -$maxsize)"
				-o "$_returnfile_sizeselect"
				"$bam"
		CMD
	# fi

	return 0
}

function alignment::_sort(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} crafts command to sort an alignment file by coordinate

			-1 <cmds>      | mandatory
			               | array to append commands to
			-t <threads>   | mandatory
			               | number of threads
			-f <sam|bam>   | mandatory
			               | path to alignment in SAM or BAM format
			-o <outbase>   | mandatory
			               | path to output directory plus alignment file prefix
			-p <tmpdir>    | mandatory
			               | path to temporary directory
			-r <var>       | optional
			               | variable to store resulting alignment file path
			-P <parameter> | optional
			               | additional samtools view parameter

			example:
			    declare -a cmds
			    declare outfile=""
			    ${FUNCNAME[-2]} -1 cmds -t 16 -i /path/to/alignment.[sam|bam] -o /path/to/outdir/alignment -p /path/to/tmpdir -r outfile

			create sorted bam:
			    commander::printcmd -a cmds # samtools [..] /path/to/alignment.[sam|bam] > /path/to/alignment.sorted.bam
			    commander::runcmd [..] -a cmds

			access sorted bam:
			    ls "\$outfile" # /path/to/outdir/alignment.sorted.bam
		EOF
		return 1
	}

	local OPTIND arg mandatory threads bam outbase _returnfile_sort inparams tmpdir
	declare -n _cmds1_sort
	while getopts '1:t:f:o:p:r:P:h' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_sort=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			f)	((++mandatory)); bam="$OPTARG";;
			o)	((++mandatory)); outbase="$OPTARG";;
			r)	declare -n _returnfile_sort=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG";; # needs to be given here, otherwise a mktemp dir will be deleted upon returning of this function
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage

	_returnfile_sort="$outbase.sorted.bam"
	local params="$inparams"

	commander::makecmd -a _cmds1_sort -s ';' -c {COMMANDER[0]}<<- CMD
		samtools sort
			$params
			-@ $threads
			-O BAM
			-T "$tmpdir/$(basename "$outbase")"
			--write-index
			-o "$_returnfile_sort##idx##${_returnfile_sort%.*}.bai"
			"$bam"
	CMD

	return 0
}

function alignment::_collate(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} crafts command to collate alignments by read name/id while re-pairing/coupling mate pairs mapped on same reference sequence

			-1 <cmds>      | mandatory
			               | array to append commands to
			-t <threads>   | mandatory
			               | number of threads
			-f <sam|bam>   | mandatory
			               | path to alignment in SAM or BAM format
			-o <outbase>   | mandatory
			               | path to output directory plus alignment file prefix
			-p <tmpdir>    | mandatory
			               | path to temporary directory
			-r <var>       | optional
			               | variable to store resulting alignment file path
			-P <parameter> | optional
			               | additional samtools view parameter

			example:
			    declare -a cmds
			    declare outfile=""
			    ${FUNCNAME[-2]} -1 cmds -t 16 -i /path/to/alignment.[sam|bam] -o /path/to/outdir/alignment -p /path/to/tmpdir -r outfile

			create sorted bam:
			    commander::printcmd -a cmds # samtools [..] /path/to/alignment.[sam|bam] > /path/to/alignment.sorted.bam
			    commander::runcmd [..] -a cmds

			access sorted bam:
			    ls "\$outfile" # /path/to/outdir/alignment.sorted.bam
		EOF
		return 1
	}

	local OPTIND arg mandatory threads bam outbase _returnfile_namesort inparams tmpdir
	declare -n _cmds1_namesort
	while getopts '1:t:f:o:p:r:P:h' arg; do
		case $arg in
			1)	((++mandatory)); _cmds1_namesort=$OPTARG;;
			t)	((++mandatory)); threads=$OPTARG;;
			f)	((++mandatory)); bam="$OPTARG";;
			o)	((++mandatory)); outbase="$OPTARG";;
			r)	declare -n _returnfile_namesort=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG";; # needs to be given here, otherwise a mktemp dir will be deleted upon returning of this function
			P)	inparams="$OPTARG";;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage

	_returnfile_namesort="$outbase.namesorted.bam"
	local params="$inparams"
	local x=$(samtools view -F 4 "$bam" | head -10000 | cat <(samtools view -H "$bam") - | samtools view -c -f 1)
	local y=$(samtools view -F 4 "$bam" | head -1 | cat <(samtools view -H "$bam") - | samtools view -c -d HI)

	if [[ $y -eq 1 && $x -gt 0 ]]; then
		# using HI is save but requires collate
		commander::makecmd -a _cmds1_namesort -s '|' -c {COMMANDER[0]}<<- CMD -c {COMMANDER[1]}<<- CMD
			{
				samtools view --no-PG -H "$bam";
				paste -d '\n'
					<(samtools view -@ $(((threads+1)/2)) -u -e 'rnext==rname' -f 2 -F 4 -f 64 "$bam" | samtools sort -t HI -n -@ $(((threads+1)/2)) -T "$tmpdir/$(basename "$outbase").nsrtR1" -u | samtools collate -@ $(((threads+1)/2)) -u -O -n 1 - "$tmpdir/$(basename "$outbase").collateR1" | samtools view -@ $(((threads+1)/2)))
					<(samtools view -@ $(((threads+1)/2)) -u -e 'rnext==rname' -f 2 -F 4 -F 64 "$bam" | samtools sort -t HI -n -@ $(((threads+1)/2)) -T "$tmpdir/$(basename "$outbase").nsrtR2" -u | samtools collate -@ $(((threads+1)/2)) -u -O -n 1 - "$tmpdir/$(basename "$outbase").collateR2" | samtools view -@ $(((threads+1)/2)));
			}
		CMD
			samtools view --no-PG -@ 100 -b -o "$_returnfile_namesort"
		CMD
	elif [[ $x -gt 0 ]]; then
		# name sort seems to work, too - shuffling via collate not necessary i.e. beeing the better name sorting for PE data
		# -e 'rnext==rname' is mainly a bwa fix for quantification with salmon
		commander::makecmd -a _cmds1_namesort -s '|' -c {COMMANDER[0]}<<- CMD -c {COMMANDER[1]}<<- CMD
			{
				samtools view --no-PG -H "$bam";
				paste -d '\n'
					<(samtools view -@ $(((threads+1)/2)) -u -e 'rnext==rname' -f 2 -F 4 -f 64 "$bam" | samtools sort -n -@ $(((threads+1)/2)) -T "$tmpdir/$(basename "$outbase").nsrtR1" -u | samtools collate -@ $(((threads+1)/2)) -u -O -n 1 - "$tmpdir/$(basename "$outbase").collateR1" | samtools view -@ $(((threads+1)/2)))
					<(samtools view -@ $(((threads+1)/2)) -u -e 'rnext==rname' -f 2 -F 4 -F 64 "$bam" | samtools sort -n -@ $(((threads+1)/2)) -T "$tmpdir/$(basename "$outbase").nsrtR2" -u | samtools collate -@ $(((threads+1)/2)) -u -O -n 1 - "$tmpdir/$(basename "$outbase").collateR2" | samtools view -@ $(((threads+1)/2)));
			}
		CMD
			samtools view --no-PG -@ 100 -b -o "$_returnfile_namesort"
		CMD
	else
		commander::makecmd -a _cmds1_namesort -s ';' -c {COMMANDER[0]}<<- CMD
			samtools collate
				--no-PG
				-l 6
				-@ $threads
				--output-fmt BAM
				-n 1
				-o "$_returnfile_namesort"
				"$bam" "$tmpdir/$(basename "$outbase")"
		CMD
	fi

	return 0
}

function alignment::_index(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} crafts command to index a coordinate sorted alignment file

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
			    ${FUNCNAME[-2]} -1 cmds -t 16 -i /path/to/alignment.bam -r idxfile

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
			"$bam" "$bai"
	CMD

	return 0
}

function alignment::tobed(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} converts alignments from BAM to BED format splitting N cigar strings and thus decoupling mate pairs

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
			    ${FUNCNAME[-2]} -t 16 -r mapper

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
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 2 ]] && _usage

	declare -n _bams_tobed="${_mapper_tobed[0]}"
	local ithreads instances=$((${#_mapper_tobed[@]}*${#_bams_tobed[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	commander::printinfo "converting alignments to bed.gz"

	local m f params="$inparams"
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
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} infers library strandness from alignment files in BAM format given a reference annotation

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
			    ${FUNCNAME[-2]} -t 16 -r mapper -x strandness -x 2

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

	declare -n _bams_inferstrandness="${_mapper_inferstrandness[0]}"
	local ithreads instances=$((${#_mapper_inferstrandness[@]}*${#_bams_inferstrandness[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -T $threads)

	declare -a cmd1 cmd2 tdirs=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.inferstrandness)")

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
		-- -lvl="$level" "$gtf" > "${tdirs[0]}/features.bed"
	CMD

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $threads -a cmd1
	fi

	local m f params
	for m in "${_mapper_inferstrandness[@]}"; do
		declare -n _bams_inferstrandness=$m
		for f in "${_bams_inferstrandness[@]}"; do

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.inferstrandness)")

			samtools view "$f" $(samtools view -H "$f" | grep -m 1 -F @SQ | cut -f 2 | cut -d : -f 2) 2> /dev/null | head -1 &> /dev/null && params="-M" || params=""

			# requires, sorted, indexed bam
			# 0 - unstranded
			# 1 - dUTP ++,-+ (FR stranded)
			# 2 - dUTP +-,++ (FR, reversely stranded)
			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD'
				samtools view
					$params
					-L "${tdirs[0]}/features.bed"
					-@ $ithreads
					-u
					"$f" |
			CMD
				samtools sort
					-O BAM
					-@ $ithreads
					-T "${tdirs[-1]}/$(basename "${f%.*}")"
					--write-index
					-o "${tdirs[-1]}/$(basename "${f%.*}").bam##idx##${tdirs[-1]}/$(basename "${f%.*}").bai";
			CMD
				{ echo "$f" &&
					infer_experiment.py
					-q 0
					-s 100000
					-i "${tdirs[-1]}/$(basename "${f%.*}").bam"
					-r "${tdirs[0]}/features.bed";
				} |
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
			${FUNCNAME[-2]} appends current paths of postprocessed alignmet files to collection data structures

			-r <mapper>   | mandatory
			              | array of array names which contain alignment paths
			              | mapper=(segemehl star); [segemehl|star]=(/outdir/[segemehl|star]/1.bam /outdir/[segemehl|star]/2.bam ..);
			-f <force>    | optional. default: false
			              | [true|false] unset collections and start from scratch

			example:
			    mapper=(segemehl star)
			    [segemehl|star]=(/path/to/[segemehl|star]/1.bam /path/to/[segemehl|star]/2.bam ..)
			    ${FUNCNAME[-2]} -r mapper
			    [segemehl|star]=(/path/to/[segemehl|star]/1.uniq.bam /path/to/[segemehl|star]/2.uniq.bam ..)
			    ${FUNCNAME[-2]} -r mapper

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
			${FUNCNAME[-2]} summarizes numbers of primary, secondary and supplementary mapped reads from alignment files

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
			    ${FUNCNAME[-2]} -t 16 -r mapper

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

	declare -n _bams_bamqc="${_mapper_bamqc[0]}"
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
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} index alignment files from collections generated via alignment::add4stats

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
			    ${FUNCNAME[-2]} -t 16 -r mapper
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

	local tdir+="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.samtools)"
	alignment::postprocess \
		-S false \
		-s $skip \
		-j index \
		-t $threads \
		-o "$tdir" \
		-r mapper_index \
		-P "$inparams"
}

function alignment::qcstats(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} generates bar plots from read count summaries of alignment files stored in collection data structures prior generated by alignment::add4stats.
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

	local filter b o all a s c odir
	declare -a cmd1 tojoin header
	for m in "${_mapper_bamstats[@]}"; do
		declare -n _bams_bamstats=$m
		odir="$outdir/$m"
		mkdir -p "$odir"
		header="size"
		tojoin=()
		echo -e "sample\ttype\tcount" > "$odir/mapping.barplot.tsv"
		for i in "${!_bams_bamstats[@]}"; do
			declare -n _mi_bamstats=$m$i # reference declaration in alignment::add4stats
			[[ ${#_mi_bamstats[@]} -eq 0 || "${_bams_bamstats[$i]}" =~ (fullpool|pseudopool|pseudorep|pseudoreplicate) ]] && continue

			b="$(basename "${_mi_bamstats[0]}")"
			b="${b%.*}"
			o="$odir/$b.stats"
			unset all filter
			for bam in "${_mi_bamstats[@]}"; do
				[[ -s "${bam%.*}.sizemetrics.tsv" ]] && header+="\t$b" && tojoin+=("${bam%.*}.sizemetrics.tsv")

				# total = primary + secondary + supplementary (primary includes unmapped reads)
				# primary mapped = mapped - secondary - supplementary
				# unmapped = primary - primary mapped
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

		if [[ $(wc -l < "$odir/mapping.barplot.tsv") -gt 2 ]]; then
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
		fi

		if [[ ${#tojoin[@]} -eq 1 ]]; then
			commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				{	echo -e "$header";
					cat "$tojoin";
				} >  "$odir/insertsizes.histogram.tsv";
			CMD
				Rscript - <<< '
					suppressMessages(library("ggpubr"));
					suppressMessages(library("tidyr"));
					args <- commandArgs(TRUE);
					intsv <- args[1];
					outfile <- args[2];
					df <- read.table(intsv, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
					df <- as.data.frame(pivot_longer(df,2:ncol(df),names_to="sample",values_to="count"));
					ggline(df, x = "size", y = "count", color = "sample", plot_type = "l");
					ggsave(outfile, width=16, height=10);
				'
			CMD
				"$odir/insertsizes.histogram.tsv" "$odir/insertsizes.histogram.pdf"
			CMD
		elif [[ ${#tojoin[@]} -gt 1 ]]; then
			commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				helper::multijoin
					-e 'NA'
					-h "$(echo -e "$header")"
					-o "$odir/insertsizes.histogram.tsv"
					-f $(printf '"%s" ' "${tojoin[@]}");
			CMD
				Rscript - <<< '
					suppressMessages(library("ggpubr"));
					suppressMessages(library("tidyr"));
					args <- commandArgs(TRUE);
					intsv <- args[1];
					outfile <- args[2];
					df <- read.table(intsv, header=T, sep="\t", stringsAsFactors=F, check.names=F, quote="");
					df <- as.data.frame(pivot_longer(df,2:ncol(df),names_to="sample",values_to="count"));
					ggline(df, x = "size", y = "count", color = "sample", plot_type = "l");
					ggsave(outfile, width=16, height=10);
				'
			CMD
				"$odir/insertsizes.histogram.tsv" "$odir/insertsizes.histogram.pdf"
			CMD
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $threads -a cmd1
	fi

	return 0
}
