#! /usr/bin/env bash
# (c) Konstantin Riege

function preprocess::dedup(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-o <outdir>    | path to
			-M <maxmemory> | amount of
			-1 <fastq1>    | array of
			-2 <fastq2>    | array of
			-3 <fastqUMI>  | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory tmpdir="${TMPDIR:-/tmp}" outdir
	declare -n _fq1_dedup _fq2_dedup _umi_dedup
	while getopts 'S:s:t:M:o:1:2:3:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1)	((++mandatory)); _fq1_dedup=$OPTARG;;
			2)	_fq2_dedup=$OPTARG;;
			3)	((++mandatory)); _umi_dedup=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "umi based de-duplication"

	declare -a cmd1
	local i o1 e1 o2 e2 instances memory catcmd
	read -r instances memory < <(configure::memory_by_instances -i 1 -M "$maxmemory")

	for i in "${!_fq1_dedup[@]}"; do

		helper::basename -f "${_fq1_dedup[$i]}" -o o1 -e e1
		e1=$(echo $e1 | cut -d '.' -f 1)
		o1="$outdir/$o1.$e1.gz"

		helper::makecatcmd -c catcmd -f "${_fq1_dedup[$i]}"

		if [[ "${_fq2_dedup[$i]}" ]]; then
			helper::basename -f "${_fq2_dedup[$i]}" -o o2 -e e2
			e2=$(echo $e2 | cut -d '.' -f 1)
			o2="$outdir/$o2.$e2.gz"

			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD
				paste <($catcmd "${_fq1_dedup[$i]}" | paste - - - -) <($catcmd "${_fq2_dedup[$i]}" | paste - - - -) <($catcmd "${_umi_dedup[$i]}" | paste - - - -)
			CMD
				awk -F '\t' -v OFS='\t' '{print $2$6$10,$0}'
			CMD
				LC_ALL=C sort --parallel=$threads -S ${memory}M -T "$tmpdir" -k1,1
			CMD
				awk '{if($1!=s){print}; s=$1}'
			CMD
				tee -i >(awk -F '\t' -v OFS='\n' '{print \$2,\$3,\$4,\$5}' | helper::pgzip -t $(((threads+1)/2)) -o "$o1") >(awk -F '\t' -v OFS='\n' '{print \$6,\$7,\$8,\$9}' | helper::pgzip -t $(((threads+1)/2)) -o "$o2") > /dev/null
			CMD
				cat
			CMD

			_fq1_dedup[$i]="$o1"
			_fq2_dedup[$i]="$o2"
		else
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD
				paste <($catcmd "${_fq1_dedup[$i]}" | paste - - - -) <($catcmd "${_umi_dedup[$i]}" | paste - - - -)
			CMD
				awk -F '\t' -v OFS='\t' '{print $2$6,$0}'
			CMD
				LC_ALL=C sort --parallel=$threads -S ${memory}M -T "$tmpdir" -k1,1
			CMD
				awk -F '\t' -v OFS='\n' '{if($1!=s){print $2,$3,$4,$5}; s=$1}'
			CMD
				helper::pgzip -t $threads -o "$o1"
			CMD
			_fq1_dedup[$i]="$o1"
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i 1 -a cmd1
	fi

	return 0
}

function preprocess::fastqc(){
	declare -a tdirs
	function _cleanup::preprocess::fastqc(){
		rm -rf "${tdirs[@]}"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-M <maxmemory>| amount of
			-o <outdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
			-a <adapter1> | array of
			-A <adapter2> | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}" maxmemory
	declare -n _fq1_fastqc _fq2_fastqc _adapter1_fastqc _adapter2_fastqc
	while getopts 'S:s:t:M:o:1:2:a:A:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1)	((++mandatory)); _fq1_fastqc=$OPTARG;;
			2)	_fq2_fastqc=$OPTARG;;
			a)	_adapter1_fastqc=$OPTARG;;
			A)	_adapter2_fastqc=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "calculating qualities"

	local instances=$((${#_fq1_fastqc[@]}+${#_fq2_fastqc[@]})) ithreads jmem jgct jcgct
	read -r instances ithreads jmem jgct jcgct < <(configure::jvm -i $instances -T $threads -m 250 -M "$maxmemory")

	declare -a cmd1 cmd2 cmd3
	local f b e i=0
	for f in {"${_fq1_fastqc[@]}","${_fq2_fastqc[@]}"}; do
		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.fastqc)")
		commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			MALLOC_ARENA_MAX=4
			JAVA_OPTS="-Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.tmpdir='$tmpdir'"
			fastqc
			-d "${tdirs[-1]}"
			-outdir "$outdir"
			"$f" 2>&1
		CMD
			sed -u '/Exception/{q 1};${/Analysis complete/!{q 1}}'
		CMD

		helper::basename -f "$f" -o b -e e
		e=$(echo $e | cut -d '.' -f 1) # if e == fastq or fq : check for ${b}_fastqc.zip else $b.${e}_fastqc.zip
		[[ $e == "fastq" || $e == "fq" ]] && f="${b}_fastqc.zip" || f="$b.${e}_fastqc.zip"
		commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			unzip -c "$outdir/$f" "${f%.*}/fastqc_data.txt" | tac
		CMD
			perl -M'List::Util q(max)' -M'Switch' -lane '{
				next unless /^\d+/;
				shift @F;
				$m=max(@F);
				exit if $m<0.001;
				$i=(grep {$F[$_]==$m} 0..$#F)[0];
				switch($i){
					case 0 {print "AGATCGGAAGAGC"}
					case 1 {print "TGGAATTCTCGGGTGCCAAGG"}
					case 2 {print "GTTCAGAGTTCTACAGTCCGACGATC"}
					case 3 {print "CTGTCTCTTATACACATCT"}
					case 4 {print "CGCCTTGGCCGT"}
				}
				exit
			}'
		CMD
		# uiversal, srna3' srna5', nextera, solexa
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c fastqc -v -b -i $threads -a cmd1
	fi

	declare -p _adapter1_fastqc | grep -q '=' && {
		mapfile -t _adapter1_fastqc < <(commander::runcmd -i $threads -a cmd2 -s 1:${#_fq1_fastqc[@]} | sort -u)
		if [[ $_fq2_fastqc ]]; then
			mapfile -t _adapter2_fastqc < <(commander::runcmd -i $threads -a cmd2 -s $((${#_fq1_fastqc[@]}+1)):$((${#_fq1_fastqc[@]}+${#_fq2_fastqc[@]})) | sort -u)
			[[ $_adapter2_fastqc ]] || _adapter2_fastqc=("${_adapter1_fastqc[@]}")
			[[ $_adapter1_fastqc ]] || _adapter1_fastqc=("${_adapter2_fastqc[@]}")
		fi
		[[ $_adapter1_fastqc ]] && commander::printinfo "Inferred adapter sequences: ${_adapter1_fastqc[*]}" || commander::warn "No adapter sequence inferred"
	}

	return 0
}

function preprocess::rmpolynt(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-o <outdir>   | path to
			-d <dinucs>   | true/false trim di-nucleotide ends
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false dinuc=true threads outdir
	declare -n _fq1_rmpolynuc _fq2_rmpolynuc
	while getopts 'S:s:d:t:o:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			d) dinuc=$OPTARG;;
			t) ((++mandatory)); threads=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1) ((++mandatory)); _fq1_rmpolynuc=$OPTARG;;
			2) _fq2_rmpolynuc=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "clipping poly N-, mono- and di-nucleotide ends"

	# -a ADAPTERX : allows partial matches, but disallow internal matches
	# -a ADAPTER$ : adapter will only be found if it is a true suffix of the read
	declare -a poly
	local i
	for i in A C G T; do
		poly+=("$(printf "$i%.0s" {1..200})X")
	done
	if $dinuc; then
		for i in AB CD GH TV; do #iupac
			poly+=("$(printf "$i%.0s" {1..200})X")
		done
	fi

	preprocess::cutadapt \
		-S false \
		-s $skip \
		-a poly \
		-A poly \
		-b false \
		-t $threads \
		-o "$outdir" \
		-1 _fq1_rmpolynuc \
		-2 _fq2_rmpolynuc

	return 0
}

function preprocess::cutadapt(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-a <adapter1> | array of
			-A <adapter2> | array of
			-b <rrbs>     | true/false
			-t <threads>  | number of
			-o <outdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir rrbs=false
	declare -n _adaptera_cutadapt _adapterA_cutadapt _fq1_cutadapt _fq2_cutadapt
	while getopts 'S:s:a:A:b:t:o:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			a) ((++mandatory)); _adaptera_cutadapt=$OPTARG;;
			A) _adapterA_cutadapt=$OPTARG;;
			b) rrbs=$OPTARG;;
			t) ((++mandatory)); threads=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1) ((++mandatory)); _fq1_cutadapt=$OPTARG;;
			2) _fq2_cutadapt=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "adapter clipping"

	declare -a cmd1
	local i o1 e1 o2 e2 params n=$((${#_adaptera_cutadapt[@]}))
	[[ $n -gt 2 ]] && n=2 # since only the best matching adapter is removed, one may run cutadapt twice with -n $n
	if $rrbs; then
		_adaptera_cutadapt=(${_adaptera_cutadapt[@]/#/NN} ${_adaptera_cutadapt[@]})
		# add instead of replace - necessary due to NN mismatch if read == ^adapter
	fi

	for i in "${!_fq1_cutadapt[@]}"; do
		helper::basename -f "${_fq1_cutadapt[$i]}" -o o1 -e e1
		e1=$(echo $e1 | cut -d '.' -f 1)
		o1="$outdir/$o1.$e1.gz"

		$rrbs && params='-O 3' || params='-O 5'
		if [[ "${_fq2_cutadapt[$i]}" ]]; then
			$rrbs && params+=' -U 2' # r=NAATT a=TT -> A, r=NNNAATT a=TT -> AA

			helper::basename -f "${_fq2_cutadapt[$i]}" -o o2 -e e2
			e2=$(echo $e2 | cut -d '.' -f 1)
			o2="$outdir/$o2.$e2.gz"

			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				cutadapt
				${_adaptera_cutadapt[@]/#/-a }
				${_adapterA_cutadapt[@]/#/-A }
				$params
				-q 20
				--trim-n
				-j $threads
				-m 18
				-o >(helper::pgzip -t $(((threads+1)/2)) -o "$o1")
				-p >(helper::pgzip -t $(((threads+1)/2)) -o "$o2")
				"${_fq1_cutadapt[$i]}" "${_fq2_cutadapt[$i]}"
				| cat
			CMD
			_fq1_cutadapt[$i]="$o1"
			_fq2_cutadapt[$i]="$o2"
		else
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				cutadapt
				${_adaptera_cutadapt[@]/#/-a }
				$params
				-q 20
				--trim-n
				-j $threads
				-m 18
				-o >(helper::pgzip -t $threads -o "$o1")
				"${_fq1_cutadapt[$i]}"
				| cat
			CMD
			_fq1_cutadapt[$i]="$o1"
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c cutadapt -v -b -i 1 -a cmd1
	fi

	return 0
}

function preprocess::trimmomatic(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-b <rrbs>     | true/false if true does not trim read starts
			-t <threads>  | number of
			-M <maxmemory>| amount of
			-o <outdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory outdir tmpdir="${TMPDIR:-/tmp}" rrbs=false
	declare -n _fq1_trimmomatic _fq2_trimmomatic
	while getopts 'S:s:t:M:b:m:o:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			M) maxmemory=$OPTARG;;
			b) rrbs=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1) ((++mandatory)); _fq1_trimmomatic=$OPTARG;;
			2) _fq2_trimmomatic=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "trimming"

	#offset 64: ASCII 64 to 106 (solexa: 59 to 106)
	#offset 33: ASCII 33 to 75
	#64 to 33: ord(char)-33+2
	#theoretical max range is 126 for all encodings, thus more reliable detection would be just min based
	#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/
	#https://www.drive5.com/usearch/manual/quality_score.html
	#od -v -A n -t u1
	declare -a cmd1
	local f catcmd params
	$rrbs || params='LEADING:20'

	for f in "${_fq1_trimmomatic[@]}"; do
		helper::makecatcmd -c catcmd -f "$f"
		commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
			$catcmd "$f" | head -4000
		CMD
			| perl -M'List::Util qw(min max)' -slne '
				BEGIN{
					$min=106
				}
				if($.%4==0){
					@x=unpack"C*";
					$min=min($min,@x);
					$max=max($max,@x);
				}
				END{
					if($min>=33 && $max<=75){
						print "phred33 $f";
					}elsif($min>=64 && $max>75 && $max<=106){
						print "phred64 $f";
					}elsif($min>=59 && $min<64 && $max>75 && $max<=106){
						print "solexa64 $f";
					}else{
						print "unknown $f";
					}
				}
			'
		CMD
			-- -f="$f"
		CMD
	done

	declare -A phred
	local l
	declare -a a mapdata
	commander::printcmd -a cmd1
	mapfile -t mapdata < <(commander::runcmd -i $threads -a cmd1)
	for l in "${mapdata[@]}"; do
		a=($l)
		phred["${a[@]:1}"]="${a[0]}"
	done

	local instances ithreads jmem jgct jcgct
	read -r instances ithreads jmem jgct jcgct < <(configure::jvm -i 1 -T $threads -M "$maxmemory")

	# trimmomatic bottleneck are number of used compression threads (4)
	declare -a cmd2 cmd3
	local i o1 o2 e1 e2
	for i in "${!_fq1_trimmomatic[@]}"; do
		helper::basename -f "${_fq1_trimmomatic[$i]}" -o o1 -e e1
		e1=$(echo $e1 | cut -d '.' -f 1)
		os1="$outdir/singletons.$o1.$e1.gz"
		o1="$outdir/$o1.$e1.gz"

		if [[ ${_fq2_trimmomatic[$i]} ]]; then
			helper::basename -f "${_fq2_trimmomatic[$i]}" -o o2 -e e2
			e2=$(echo $e2 | cut -d '.' -f 1)
			os2="$outdir/singletons.$o2.$e2.gz"
			o2="$outdir/$o2.$e2.gz"

			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				trimmomatic
				-Xmx${jmem}m
				-XX:ParallelGCThreads=$jgct
				-XX:ConcGCThreads=$jcgct
				-Djava.io.tmpdir="$tmpdir"
				PE
				-threads $threads
				-${phred["${_fq1_trimmomatic[$i]}"]}
				"${_fq1_trimmomatic[$i]}" "${_fq2_trimmomatic[$i]}"
				>(helper::pgzip -t $(((threads+1)/2)) -o "$o1") >(helper::pgzip -t $(((threads+1)/2)) -o "$os1")
				>(helper::pgzip -t $(((threads+1)/2)) -o "$o2") >(helper::pgzip -t $(((threads+1)/2)) -o "$os2")
				$params
				SLIDINGWINDOW:5:20
				MINLEN:18
				TOPHRED33
				| cat
			CMD
			_fq1_trimmomatic[$i]="$o1"
			_fq2_trimmomatic[$i]="$o2"
		else
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				trimmomatic
				-Xmx${jmem}m
				-XX:ParallelGCThreads=$jgct
				-XX:ConcGCThreads=$jcgct
				-Djava.io.tmpdir="$tmpdir"
				SE
				-threads $threads
				-${phred["${_fq1_trimmomatic[$i]}"]}
				"${_fq1_trimmomatic[$i]}"
				>(helper::pgzip -t $threads -o "$o1")
				$params
				SLIDINGWINDOW:5:20
				MINLEN:18
				TOPHRED33
				| cat
			CMD
			_fq1_trimmomatic[$i]="$o1"
		fi
	done

	if $skip; then
		commander::printcmd -a cmd2
	else
		commander::runcmd -v -b -i 1 -a cmd2
	fi

	return 0
}

function preprocess::rcorrector(){
	declare -a tdirs
	function _cleanup::preprocess::rcorrector(){
		rm -rf "${tdirs[@]}"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-o <outdir>   | absolute path to
			-1 <fastq1>   | array of absolute paths
			-2 <fastq2>   | array of absolute paths
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}"
	declare -n _fq1_rcorrector _fq2_rcorrector
	while getopts 'S:s:t:o:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir"; outdir=$(realpath -se "$outdir");;
			1) ((++mandatory)); _fq1_rcorrector=$OPTARG;;
			2) _fq2_rcorrector=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "correcting read errors"

	declare -a cmd1 cmd2
	local i o1 e1 o2 e2 r1 r2
	for i in "${!_fq1_rcorrector[@]}"; do
		helper::basename -f "${_fq1_rcorrector[$i]}" -o o1 -e e1
		e1=$(echo $e1 | cut -d '.' -f 1) # if e1 == fastq or fq : mv $o1.cor.fq.gz $o1.$e1.gz else mv $o1.$e1.cor.fq.gz $o1.$e1.gz
		o1="$outdir/$o1"
		[[ $e1 == "fastq" || $e1 == "fq" ]] && r1="cor.fq" || r1="$e1.cor.fq"

		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.rcorrector)")
		if [[ ${_fq2_rcorrector[$i]} ]]; then
			helper::basename -f "${_fq2_rcorrector[$i]}" -o o2 -e e2
			e2=$(echo $e2 | cut -d '.' -f 1)
			o2="$outdir/$o2"
			[[ $e2 == "fastq" || $e2 == "fq" ]] && r2="cor.fq" || r2="$e2.cor.fq"

			readlink -e "${_fq1_rcorrector[$i]}" | file -f - | grep -qF 'compressed' && {
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					cd "${tdirs[-1]}"
				CMD
					run_rcorrector.pl
					-1 "$(realpath -se "${_fq1_rcorrector[$i]}")"
					-2 "$(realpath -se "${_fq2_rcorrector[$i]}")"
					-od "$outdir"
					-t $threads
				CMD
					mv "$o1.$r1.gz" "$o1.$e1.gz"
				CMD
					mv "$o2.$r2.gz" "$o2.$e2.gz"
				CMD
			} || {
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD
					cd "${tdirs[-1]}"
				CMD
					exec 11>&1; exec 12>&1
				CMD
					ln -sfn "/dev/fd/11" "$(basename "$o1").$r1"
				CMD
					ln -sfn "/dev/fd/12" "$(basename "$o2").$r2"
				CMD
					run_rcorrector.pl
					-1 "$(realpath -se "${_fq1_rcorrector[$i]}")"
					-2 "$(realpath -se "${_fq2_rcorrector[$i]}")"
					-od "${tdirs[-1]}"
					-t $threads
					11> >(helper::pgzip -t $(((threads+1)/2)) -o "$o1.$e1.gz")
					12> >(helper::pgzip -t $(((threads+1)/2)) -o "$o2.$e2.gz")
					| cat
				CMD
					exec 11>&-; exec 12>&-
				CMD
			}
			# -stdout is broken for PE data. reports R1 twice instead of R1 and R2

			_fq1_rcorrector[$i]="$o1.$e1.gz"
			_fq2_rcorrector[$i]="$o2.$e2.gz"
		else
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cd "${tdirs[-1]}"
			CMD
				run_rcorrector.pl
				-s "$(realpath -se "${_fq1_rcorrector[$i]}")"
				-stdout
				-t $threads
				> >(helper::pgzip -t $threads -o "$o1.$e1.gz")
				| cat
			CMD
			_fq1_rcorrector[$i]="$o1.$e1.gz"
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c rcorrector -v -b -i 1 -a cmd1
	fi

	return 0
}


function preprocess::sortmerna_new(){
	declare -a tdirs
	function _cleanup::preprocess::sortmerna(){
		rm -rf "${tdirs[@]}"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-o <outdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}"
	declare -n _fq1_sortmerna _fq2_sortmerna
	while getopts 'S:s:t:i:o:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1) ((++mandatory)); _fq1_sortmerna=$OPTARG;;
			2) _fq2_sortmerna=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "filtering rRNA fragments"

	declare -a cmd1 cmd2
	local i f catcmd o1 o2 e1 e2
	for i in "${!_fq1_sortmerna[@]}"; do
		helper::basename -f "${_fq1_sortmerna[$i]}" -o o1 -e e1
		e1=$(echo $e1 | cut -d '.' -f 1)

		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.sortmerna)")

		helper::makecatcmd -c catcmd -f "${_fq1_sortmerna[$i]}"

		if [[ ${_fq2_sortmerna[$i]} ]]; then
			helper::basename -f "${_fq2_sortmerna[$i]}" -o o2 -e e2
			e2=$(echo $e2 | cut -d '.' -f 1)

			if [[ "$catcmd" == "cat" ]]; then
				# exec/ln trick does not work any longer. outfile will be overridden
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					sortmerna
						--index 0
						\$(ls \$CONDA_PREFIX/rRNA_databases/*.fasta | xargs -I {} printf "--ref %q\n" "{}")
						--idx-dir "\$CONDA_PREFIX/rRNA_databases/index"
						--workdir "${tdirs[-1]}"
						--threads $threads
						--reads "${_fq1_sortmerna[$i]}"
						--reads "${_fq2_sortmerna[$i]}"
						--fastx
						--paired_out
						--out2
						--no-best
						--num_alignments 1
						--aligned "$outdir/rRNA.$o1"
						--other "$outdir/$o1"
				CMD

				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					helper::pgzip -t $threads -f "${tdirs[-1]}/${o1}_fwd.fq" -o "$outdir/$o1.$e1.gz"
				CMD
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					helper::pgzip -t $threads -f "${tdirs[-1]}/${o1}_rev.fq" -o "$outdir/$o2.$e2.gz"
				CMD
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					helper::pgzip -t $threads -f "${tdirs[-1]}/rRNA.${o1}_fwd.fq" -o "$outdir/rRNA.$o1.$e1.gz"
				CMD
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					helper::pgzip -t $threads -f "${tdirs[-1]}/rRNA.${o1}_rev.fq" - "$outdir/rRNA.$o2.$e2.gz"
				CMD
			else
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					sortmerna
						--index 0
						\$(ls \$CONDA_PREFIX/rRNA_databases/*.fasta | xargs -I {} printf "--ref %q\n" "{}")
						--idx-dir "\$CONDA_PREFIX/rRNA_databases/index"
						--workdir "${tdirs[-1]}"
						--threads $threads
						--reads "${_fq1_sortmerna[$i]}"
						--reads "${_fq2_sortmerna[$i]}"
						--fastx
						--paired_out
						--out2
						--no-best
						--num_alignments 1
						--aligned "$outdir/rRNA.$o1"
						--other "$outdir/$o1"
				CMD
					mv "$outdir/${o1}_fwd.fq.gz" "$outdir/$o1.$e1.gz";
					mv "$outdir/${o1}_rev.fq.gz" "$outdir/$o2.$e2.gz";
					mv "$outdir/rRNA.${o1}_fwd.fq.gz" "$outdir/rRNA.$o1.$e1.gz";
					mv "$outdir/rRNA.${o1}_rev.fq.gz" "$outdir/rRNA.$o2.$e2.gz"
				CMD
			fi

			_fq1_sortmerna[$i]="$outdir/$o1.$e1.gz"
			_fq2_sortmerna[$i]="$outdir/$o2.$e2.gz"
		else
			if [[ "$catcmd" == "cat" ]]; then
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					sortmerna
						--index 0
						\$(ls \$CONDA_PREFIX/rRNA_databases/*.fasta | xargs -I {} printf "--ref %q\n" "{}")
						--idx-dir "\$CONDA_PREFIX/rRNA_databases/index"
						--workdir "${tdirs[-1]}"
						--threads $threads
						--reads "${_fq1_sortmerna[$i]}"
						--fastx
						--no-best
						--num_alignments 1
						--aligned "${tdirs[-1]}/rRNA.$o1"
						--other "${tdirs[-1]}/$o1"
				CMD

				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					helper::pgzip -t $threads -f "${tdirs[-1]}/$o1.fq" -o "$outdir/$o1.$e1.gz"
				CMD
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					helper::pgzip -t $threads -f "${tdirs[-1]}/rRNA.$o1.fq" -o "$outdir/rRNA.$o1.$e1.gz"
				CMD
			else
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					sortmerna
						--index 0
						\$(ls \$CONDA_PREFIX/rRNA_databases/*.fasta | xargs -I {} printf "--ref %q\n" "{}")
						--idx-dir "\$CONDA_PREFIX/rRNA_databases/index"
						--workdir "${tdirs[-1]}"
						--threads $threads
						--reads "${_fq1_sortmerna[$i]}"
						--fastx
						--no-best
						--num_alignments 1
						--aligned "$outdir/rRNA.$o1"
						--other "$outdir/$o1"
				CMD
			fi

			_fq1_sortmerna[$i]="$outdir/$o1.$e1.gz"
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c sortmerna -v -b -i 1 -a cmd1
		commander::runcmd -c sortmerna -v -b -i 1 -a cmd2
	fi

	return 0
}

function preprocess::sortmerna(){
	declare -a tdirs
	function _cleanup::preprocess::sortmerna(){
		rm -rf "${tdirs[@]}"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-o <outdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}"
	declare -n _fq1_sortmerna _fq2_sortmerna
	while getopts 'S:s:t:i:o:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1) ((++mandatory)); _fq1_sortmerna=$OPTARG;;
			2) _fq2_sortmerna=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "filtering rRNA fragments"

	declare -a cmd1 cmd2 cmd3
	local i catcmd tmp o1 o2 or1 or2 e1 e2
	for i in "${!_fq1_sortmerna[@]}"; do
		helper::basename -f "${_fq1_sortmerna[$i]}" -o o1 -e e1
		e1=$(echo $e1 | cut -d '.' -f 1)

		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.sortmerna)")
		tmp="${tdirs[-1]}/$o1"

		or1="$outdir/rRNA.$o1.$e1.gz"
		o1="$outdir/$o1.$e1.gz"

		helper::makecatcmd -c catcmd -f "${_fq1_sortmerna[$i]}"

		# sortmerna v2.1 input must not be compressed (v.3.* creates empty files)
		# outfile gets extension from input file
		# in.fq.bz2 > in.fq + rRNA.out|out -> rRNA.out.fq|out.fq -> rRNA.out.fq.gz|out.fq.gz
		if [[ ${_fq2_sortmerna[$i]} ]]; then
			helper::basename -f "${_fq2_sortmerna[$i]}" -o o2 -e e2
			e2=$(echo $e2 | cut -d '.' -f 1)
			or2="$outdir/rRNA.$o2.$e2.gz"
			o2="$outdir/$o2.$e2.gz"

			commander::makecmd -a cmd1 -s '|' -o "$tmp.merged.$e1" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD'
				paste <($catcmd "${_fq1_sortmerna[$i]}") <($catcmd "${_fq2_sortmerna[$i]}")
			CMD
				paste - - - -
			CMD
				awk -F '\t' -v OFS='\n' '{print $1,$3,$5,$7; print $2,$4,$6,$8}'
			CMD

			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD
				exec 11>&1; exec 12>&1
			CMD
				ln -sfn /dev/fd/11 "$tmp.ok.$e1"
			CMD
				ln -sfn /dev/fd/12 "$tmp.rRNA.$e1"
			CMD
				sortmerna
				--ref \$(ls \$CONDA_PREFIX/rRNA_databases/*.fasta | xargs -I {} bash -c 'printf ":%q,%q" "\$1" "\$(dirname "\$1")/index/\$(basename "\$1" .fasta)-L18"' bash {} | sed 's/://')
				--reads "$tmp.merged.$e1"
				--num_alignments 1
				--fastx
				--paired_out
				--aligned "$tmp.rRNA"
				--other "$tmp.ok"
				-a $threads
				11> >(sed -E '/^\s*$/d' | paste - - - - | tee -i >(sed -n '1~2{s/\t/\n/gp}' | helper::pgzip -t $(((threads+1)/2)) -o "$o1") >(sed -n '2~2{s/\t/\n/gp}' | helper::pgzip -t $(((threads+1)/2)) -o "$o2") > /dev/null)
				12> >(sed -E '/^\s*$/d' | paste - - - - | tee -i >(sed -n '1~2{s/\t/\n/gp}' | helper::pgzip -t $(((threads+1)/2)) -o "$or1") >(sed -n '2~2{s/\t/\n/gp}' | helper::pgzip -t $(((threads+1)/2)) -o "$or2") > /dev/null)
				| cat
			CMD
				rm -f "$tmp.merged.$e1"; exec 11>&-; exec 12>&-
			CMD
			# attention: sometimes sortmerna inserts empty lines - use sed /^\s*$/d

			_fq1_sortmerna[$i]="$o1"
			_fq2_sortmerna[$i]="$o2"
		else
			helper::makecatcmd -c catcmd -f "${_fq1_sortmerna[$i]}"
			[[ $catcmd == "cat" ]] && {
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					ln -sfn "${_fq1_sortmerna[$i]}" "$tmp.$e1"
				CMD
			} || {
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					$catcmd "${_fq1_sortmerna[$i]}" > "$tmp.$e1"
				CMD
			}

			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD
				exec 11>&1; exec 12>&1
			CMD
				ln -sfn /dev/fd/11 "$tmp.ok.$e1"
			CMD
				ln -sfn /dev/fd/12 "$tmp.rRNA.$e1"
			CMD
				sortmerna
				--ref \$(ls \$CONDA_PREFIX/rRNA_databases/*.fasta | xargs -I {} bash -c 'printf ":%q,%q" "\$1" "\$(dirname "\$1")/index/\$(basename "\$1" .fasta)-L18"' bash {} | sed 's/://')
				--reads "$tmp.$e1"
				--fastx
				--aligned "$tmp.rRNA"
				--other "$tmp.ok"
				-a $threads
				11> >(sed -E '/^\s*$/d' | helper::pgzip -t $threads -o "$o1")
				12> >(sed -E '/^\s*$/d' | helper::pgzip -t $threads -o "$or1")
				| cat
			CMD
				rm -f "$tmp.$e1"; exec 11>&-; exec 12>&-
			CMD

			_fq1_sortmerna[$i]="$o1"
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -c sortmerna -v -b -i 1 -a cmd2
	fi

	return 0
}

function preprocess::add4stats(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-r <qualdirs>  | array of
			-a <path>      | to add
			-1 <fastq1>    | array of
			-2 <fastq2>    | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory
	declare -n _qdirs_add4stats _fq1_add4stats _fq2_add4stats
	while getopts 'r:a:1:2:' arg; do
		case $arg in
			r)	((++mandatory)); _qdirs_add4stats=$OPTARG;;
			a)	((++mandatory)); qdir="$OPTARG";;
			1)	((++mandatory)); _fq1_add4stats=$OPTARG;;
			2)	_fq2_add4stats=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	local f i=${#_qdirs_add4stats[@]}
	_qdirs_add4stats+=("$qdir")
	declare -g -a preprocess$i
	declare -n _pi_add4stats=preprocess$i

	for f in {"${_fq1_add4stats[@]}","${_fq2_add4stats[@]}"}; do
		_pi_add4stats+=("$f")
	done

	return 0
}

function preprocess::qcstats(){
	local tmp
	function _cleanup::preprocess::qcstats(){
		rm -f "$tmp"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-f <force>    | true/false rerun fastqc
			-t <threads>  | number of
			-i <qualdirs> | array of
			-o <outdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory threads skip=false outdir tmpdir="${TMPDIR:-/tmp}" force=false
	declare -n _qualdirs_qcstats _fq1_qcstats _fq2_qcstats
	while getopts 'S:s:t:f:i:o:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			f) force=$OPTARG;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); _qualdirs_qcstats=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1) ((++mandatory)); _fq1_qcstats=$OPTARG;;
			2) _fq2_qcstats=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	declare -a cmdqc
	local qdir i f b e qczip
	for i in "${!_qualdirs_qcstats[@]}"; do
		qdir="${_qualdirs_qcstats[$i]}"
		declare -n _pi_qcstats=preprocess$i
		for f in "${_pi_qcstats[@]}"; do
			helper::basename -f "$f" -o b -e e
			e=$(echo $e | cut -d '.' -f 1) # if e == fastq or fq : check for ${b}_fastqc.zip else $b.${e}_fastqc.zip
			[[ $e == "fastq" || $e == "fq" ]] && qczip="${b}_fastqc.zip" || qczip="$b.${e}_fastqc.zip"
			if $force || [[ ! -s "$qdir/$qczip" ]]; then
				commander::makecmd -a cmdqc -s ' ' -c {COMMANDER[0]}<<- CMD
					fq=("$f");
					preprocess::fastqc -S false -s $skip -t 1 -p "$tmpdir" -o "$qdir" -1 fq
				CMD
				# rescue enables to store quality directories without running fastqc after each fastq processing step - instead run fastqc on all files in parallel
			fi
		done
	done
	if [[ ${#cmdqc[@]} -gt 0 ]]; then
		commander::runcmd -b -i $threads -a cmdqc
	fi

	commander::printinfo "plotting preprocessing stats"

	tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.tsv)"
	local o c multiplier tool
	declare -a counts
	echo -e "sample\ttype\tcount" > "$outdir/preprocessing.barplot.tsv"
	for i in "${!_fq1_qcstats[@]}"; do
		helper::basename -f "${_fq1_qcstats[$i]}" -o b -e e
		e=$(echo $e | cut -d '.' -f 1) # if e == fastq or fq : check for ${b}_fastqc.zip else $b.${e}_fastqc.zip
		[[ $e == "fastq" || $e == "fq" ]] && qczip="${b}_fastqc.zip" || qczip="$b.${e}_fastqc.zip"
		o="$outdir/$b.stats"
		multiplier=1
		[[ "${_fq2_qcstats[$i]}" ]] && multiplier=2
		rm -f $o
		for qdir in "${_qualdirs_qcstats[@]}"; do
			tool=$(basename "$qdir")
			c=$(unzip -c "$qdir/$qczip" "${qczip%.*}/fastqc_data.txt" | grep -m 1 -F Total | awk -v mult=$multiplier '{print $3*mult}')
			counts+=($c)
			echo -e "$b\t$tool reads\t$c" >> $o
			perl -sle 'print join"\t",("$sample ($all)","$tool reads",(100*$c/$all))' -- -all=${counts[$((i*${#_qualdirs_qcstats[@]}))]} -c=$c -sample=$b -tool=$tool
		done > "$tmp" # strange!!! if piped directly into tac - tac's awk implementation fails - not a shournal raceexception bug!
		tac "$tmp" | awk -F '\t' '{OFS="\t"; if(c){$NF=$NF-c} c=c+$NF; print}' | tac >> "$outdir/preprocessing.barplot.tsv"
	done

	declare -a cmd1
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
				ggtitle("Preprocessing") + xlab("Sample") + ylab("Readcount in %") +
				theme_bw() + guides(fill=guide_legend(title=NULL)) +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
				geom_bar(position = "fill", stat = "identity") +
				scale_y_continuous(labels = percent_format());
			graphics.off();
		'
	CMD
		"$outdir/preprocessing.barplot.tsv"  "$outdir/preprocessing.barplot.pdf"
	CMD

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $threads -a cmd1
	fi

	return 0
}
