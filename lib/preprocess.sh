#! /usr/bin/env bash
# (c) Konstantin Riege

function preprocess::dedup(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "removing read duplicates"

	declare -a cmd1
	local i o1 e1 o2 e2 instances memory
	read -r instances memory < <(configure::memory_by_instances -i 1 -M "$maxmemory")

	for i in "${!_fq1_dedup[@]}"; do

		helper::basename -f "${_fq1_dedup[$i]}" -o o1 -e e1
		e1=$(echo $e1 | cut -d '.' -f 1)
		o1="$outdir/$o1.$e1.gz"

		if [[ "${_fq2_dedup[$i]}" ]]; then
			helper::basename -f "${_fq2_dedup[$i]}" -o o2 -e e2
			e2=$(echo $e2 | cut -d '.' -f 1)
			o2="$outdir/$o2.$e2.gz"

			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD
				paste <(helper::cat -f "${_fq1_dedup[$i]}" | paste - - - -) <(helper::cat -f "${_fq2_dedup[$i]}" | paste - - - -) <(helper::cat -f "${_umi_dedup[$i]}" | paste - - - -)
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
			# paste can not be used if files differ in length and join needs alnum sorting. and pay attention to added info not added to umi read ids
			# -> solution works, but fastq deduplication should actually not applied on trimmed/clipped/... data
			# commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD
			# 	awk -v f=<(helper::cat -f "${_umi_dedup[$i]}" | paste - - - - | awk -v OFS='\t' '{print \$1,\$0}') -F '\t' '{getline l < f; split(l,a,FS); while(\$1!=a[1]){getline l < f; split(l,a,FS)} print \$3""a[3],\$2,\$3,\$4,\$5}' <(helper::cat -f "${_fq1_dedup[$i]}" | paste - - - - | awk -v OFS='\t' '{print \$1,\$0}')
			# CMD
			# 	LC_ALL=C sort --parallel=$threads -S ${memory}M -T "$tmpdir" -k1,1
			# CMD
			# 	awk -F '\t' -v OFS='\n' '{if($1!=s){print $2,$3,$4,$5}; s=$1}'
			# CMD
			# 	helper::pgzip -t $threads -o "$o1"
			# CMD

			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD
				paste <(helper::cat -f "${_fq1_dedup[$i]}" | paste - - - -) <(helper::cat -f "${_umi_dedup[$i]}" | paste - - - -)
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
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "calculating qualities"

	local instances=$((${#_fq1_fastqc[@]}+${#_fq2_fastqc[@]})) ithreads
	read -r instances ithreads < <(configure::instances_by_memory -T $threads -m 1024 -M "$maxmemory")

	declare -a cmdchk=("fastqc -v | sed -E 's/.*v([0-9]+\.[0-9]+).*/\1/'")
	local version=$(commander::runcmd -c fastqc -a cmdchk)
	# increase default v0.12:512m|v0.11:250m java xmx via threads parameter trick. (workaround for fastqc freezing at 95%)
	ithreads=$(awk '{if ($0<=0.11){print 4}else{print 2}}' <<< $version)

	declare -a tdirs cmd1 cmd2
	local f b e
	for f in {"${_fq1_fastqc[@]}","${_fq2_fastqc[@]}"}; do
		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.fastqc)")

		# also possible: cat fastq | fastqc -o <outdir> stdin:<basename>
		# Illumina Small RNA 5' Adapter	GATCGTCGGACT
		commander::makecmd -m -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
			cat <<- EOF | MALLOC_ARENA_MAX=4 fastqc -a /dev/stdin -t $ithreads -d "${tdirs[-1]}" -outdir "$outdir" "$f" 2>&1 | sed -u '/Exception/{q 1};\${/Analysis complete/!{q 1}}'
				Illumina Universal Adapter	AGATCGGAAGAGC
				Nextera Transposase Sequence	CTGTCTCTTATA
				SOLID Small RNA Adapter	CGCCTTGGCCGT
				10x Genomics TSO	CCCATGTACTCTGCGTTGATACCACTGCTT
				Illumina Small RNA 3' Adapter	TGGAATTCTCGG
				Tecan NuGEN Methyl-Seq Adapter	AAATCAAAAAAAC
			EOF
		CMD

		helper::basename -f "$f" -o b -e e
		e=$(echo $e | cut -d '.' -f 1) # if e == fastq or fq : check for ${b}_fastqc.zip else $b.${e}_fastqc.zip
		[[ "$e" == "fastq" || "$e" == "fq" ]] && f="${b}_fastqc.zip" || f="$b.${e}_fastqc.zip"
		# attention: nugen adapter search causes low false positive rates in R2 seqeunces (<0.1) likewise Small RNA adapter in R1 (<0.01)
		commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
			unzip -c "$outdir/$f" "${f%.*}/fastqc_data.txt" | tac
		CMD
			perl -F'\t' -M'List::Util q(max)' -M'Switch' -lane '{
				next unless /^\d+/;
				$m=max(@F[1..4]);
				if($m>=0.001){
					$i=(grep {$F[$_]==$m} 1..4)[0];
				} else {
					$m=max(@F[5..$#F]);
					$i=(grep {$F[$_]==$m} 5..$#F)[0];
					if ($i==5){
						exit if $m<0.01;
					} else {
						exit if $m<0.1;
					}
				}
				switch($i){
					case 1 {print "AGATCGGAAGAGC"}
					case 2 {print "CTGTCTCTTATACACATCT"}
					case 3 {print "CGCCTTGGCCGT"}
					case 4 {print "CCCATGTACTCTGCGTTGATACCACTGCTT"}
					case 5 {print "TGGAATTCTCGGGTGCCAAGG"}
					case 6 {print "AAATCAAAAAAAC"}
				}
				exit
			}'
		CMD
			tee "$($skip && echo "/dev/null" || echo "$outdir/${f%_fastqc.zip}.adapter")"
		CMD
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c fastqc -v -b -i $instances -a cmd1
	fi

	declare -p _adapter1_fastqc | grep -q '=' && {
		commander::printinfo "inferring adapter"
		commander::printcmd -a cmd2

		mapfile -t _adapter1_fastqc < <(commander::runcmd -i $threads -a cmd2 -s 1:${#_fq1_fastqc[@]} | sort -u)

		if [[ $_fq2_fastqc ]]; then

			mapfile -t _adapter2_fastqc < <(commander::runcmd -i $threads -a cmd2 -s $((${#_fq1_fastqc[@]}+1)):$((${#_fq1_fastqc[@]}+${#_fq2_fastqc[@]})) | sort -u)

			if ! [[ $_adapter1_fastqc && $_adapter2_fastqc ]]; then
				read -r instances maxmemory < <(configure::memory_by_instances -i 1 -T $threads -M "$maxmemory")
				local f1 f2 b1 b2 e1 e2 l=20

				[[ $_adapter1_fastqc ]] && l=$(printf '%s\n' "${_adapter1_fastqc[@]}" | awk '{print length($1)}' | sort -nr | head -1)
				[[ $_adapter2_fastqc ]] && l=$(printf '%s\n' "${_adapter2_fastqc[@]}" | awk '{print length($1)}' | sort -nr | head -1)

				declare -a tomerge1 tomerge2 cmd3

				for i in ${!_fq1_fastqc[@]}; do
					f1="${_fq1_fastqc[$i]}"
					f2="${_fq2_fastqc[$i]}"

					helper::basename -f "$f1" -o b1 -e e1
					helper::basename -f "$f2" -o b2 -e e2

					if [[ "$b1" == "$b2" ]]; then
						b1+="$e1"
						b2+="$e2"
					fi

					tomerge1+=("$outdir/$b1.adapter")
					tomerge2+=("$outdir/$b2.adapter")

					commander::makecmd -a cmd3 -s ' ' -o "${tdirs[$i]}/adapter.tsv" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD
						bbmerge.sh
							-Xmx${maxmemory}m
							reads=1000000
							t=$threads
							in1="${_fq1_fastqc[$i]}"
							in2="${_fq2_fastqc[$i]}"
							outa=/dev/stdout |
					CMD
						awk -v l=$l
					CMD
						'!/^>/{sub(/^N$/,"",$1); o[i++]=substr($1,0,l)}END{print o[0]"\t"o[1]}' |
					CMD
						tee >(cut -f 1 > "$outdir/$b1.adapter") >(cut -f 2 > "$outdir/$b2.adapter") > /dev/null | cat
					CMD
				done

				if $skip; then
					commander::printcmd -a cmd3
				else
					commander::runcmd -c bbmap -v -b -i 1 -a cmd3
				fi
				[[ $_adapter1_fastqc ]] || _adapter1_fastqc=($(sort -u "${tomerge1[@]}"))
				[[ $_adapter2_fastqc ]] || _adapter2_fastqc=($(sort -u "${tomerge2[@]}"))
			fi
			[[ $_adapter1_fastqc ]] && commander::printinfo "Inferred first mate adapter sequence(s): ${_adapter1_fastqc[*]}" || commander::warn "No first mate adapter sequence inferred"
			[[ $_adapter2_fastqc ]] && commander::printinfo "Inferred mate pair adapter sequence(s): ${_adapter2_fastqc[*]}" || commander::warn "No mate pair adapter sequence inferred"
		else
			[[ $_adapter1_fastqc ]] && commander::printinfo "Inferred adapter sequence(s): ${_adapter1_fastqc[*]}" || commander::warn "No adapter sequence inferred"
		fi
	}

	return 0
}

function preprocess::rmpolynt(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
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

	local OPTIND arg mandatory skip=false dinuc=false threads outdir
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "clipping poly mono- and di-nucleotide ends"

	# -a ADAPTERX : allows partial matches, but disallow internal matches
	# -a ADAPTER$ : adapter will only be found if it is a true suffix of the read
	declare -a poly
	local i
	for i in A C G T; do
		poly+=("$(printf "$i%.0s" {1..200})X")
	done
	# or -a {A}{100}X
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
			${FUNCNAME[-2]} usage:
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
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
				<(helper::cat -f "${_fq1_cutadapt[$i]}") <(helper::cat -f "${_fq2_cutadapt[$i]}")
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
				<(helper::cat -f "${_fq1_cutadapt[$i]}")
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
			${FUNCNAME[-2]} usage:
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "trimming"

	#offset 64: ASCII 64 to 106 (solexa: 59 to 106)
	#offset 33: ASCII 33 to 88
	#64 to 33: ord(char)-33+2
	#theoretical max range is 126 for all encodings, thus more reliable detection would be just min based
	#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2847217/
	#https://www.drive5.com/usearch/manual/quality_score.html
	#od -v -A n -t u1
	declare -a cmd1
	local f params
	$rrbs || params='LEADING:20'

	for f in "${_fq1_trimmomatic[@]}"; do
		commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
			helper::cat -f "$f" | head -4000
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
					if($min>=33 && $max<=88){
						print "phred33 $f";
					}elsif($min>=64 && $max>88 && $max<=106){
						print "phred64 $f";
					}elsif($min>=59 && $min<64 && $max>88 && $max<=106){
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

			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
				trimmomatic
				-Xmx${jmem}m
				-XX:ParallelGCThreads=$jgct
				-XX:ConcGCThreads=$jcgct
				-Djava.io.tmpdir="$tmpdir"
				PE
				-threads $threads
				-${phred["${_fq1_trimmomatic[$i]}"]}
				<(helper::cat -f "${_fq1_trimmomatic[$i]}") <(helper::cat -f "${_fq2_trimmomatic[$i]}")
				>(helper::pgzip -t $(((threads+1)/2)) -o "$o1") >(helper::pgzip -t $(((threads+1)/2)) -o "$os1")
				>(helper::pgzip -t $(((threads+1)/2)) -o "$o2") >(helper::pgzip -t $(((threads+1)/2)) -o "$os2")
				$params
				SLIDINGWINDOW:5:20
				MINLEN:18
				TOPHRED33
			CMD
				sed -u '/Exception/{q 1};${/Completed successfully/!{q 1}}'
			CMD
			_fq1_trimmomatic[$i]="$o1"
			_fq2_trimmomatic[$i]="$o2"
		else
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
				trimmomatic
				-Xmx${jmem}m
				-XX:ParallelGCThreads=$jgct
				-XX:ConcGCThreads=$jcgct
				-Djava.io.tmpdir="$tmpdir"
				SE
				-threads $threads
				-${phred["${_fq1_trimmomatic[$i]}"]}
				<(helper::cat -f "${_fq1_trimmomatic[$i]}")
				>(helper::pgzip -t $threads -o "$o1")
				$params
				SLIDINGWINDOW:5:20
				MINLEN:18
				TOPHRED33
			CMD
				sed -u '/Exception/{q 1};${/Completed successfully/!{q 1}}'
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
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "correcting read errors"

	declare -a tdirs cmd1 cmd2
	local i o1 e1 o2 e2 r1 r2

	declare -a cmdchk=("conda list | grep -F rcorrector | awk '{print \$2}' | sed 's/\./\t/' | awk '\$1>1||\$2>=0.7{print true; exit}{print false}'")
	local workaround=$(commander::runcmd -c rcorrector -a cmdchk)

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

			if $workaround; then
				# -stdout is broken for PE data. reports R1 twice instead of R1 and R2
				# does not work with dynamic FD
				readlink -e "${_fq1_rcorrector[$i]}" | file -b --mime-type -f - | grep -qF -e 'gzip' -e 'bzip2' && {
					commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD
						cd "${tdirs[-1]}"
					CMD
						exec 11>&1; exec 12>&1
					CMD
						ln -sfn "/dev/fd/11" "$(basename "$o1").$r1.gz"
					CMD
						ln -sfn "/dev/fd/12" "$(basename "$o2").$r2.gz"
					CMD
						run_rcorrector.pl
							-1 "$(realpath -se "${_fq1_rcorrector[$i]}")"
							-2 "$(realpath -se "${_fq2_rcorrector[$i]}")"
							-od "${tdirs[-1]}"
							-t $threads
							11> >(helper::index -o "$o1.$e1.gz")
							12> >(helper::index -o "$o2.$e2.gz")
						| cat
					CMD
						exec 11>&-; exec 12>&-
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
			else
				# stdout fixed in v1.0.6 and -tmpd added in v1.0.7
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					run_rcorrector.pl
					-tmpd "${tdirs[-1]}"
					-1 "$(realpath -se "${_fq1_rcorrector[$i]}")"
					-2 "$(realpath -se "${_fq2_rcorrector[$i]}")"
					-t $threads
					-stdout
				CMD
					paste - - - -
				CMD
					tee -i
					>(sed -n '1~2{s/\t/\n/gp}' | helper::pgzip -t $(((threads+1)/2)) -o "$o1.$e1.gz")
					>(sed -n '2~2{s/\t/\n/gp}' | helper::pgzip -t $(((threads+1)/2)) -o "$o2.$e2.gz") > /dev/null)
				CMD
					cat
				CMD
			fi

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

function preprocess::sortmerna(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-o <outdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	# sortmerna --version |& grep version | tail -1 | grep -oE '[0-9.]+' | head -1

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
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "filtering rRNA fragments"

	declare -a cmdchk=("sortmerna --version |& grep version | tail -1 | grep -oE '[0-9]' | head -1")
	local version=$(commander::runcmd -c sortmerna -a cmdchk)

	declare -a tdirs cmd1 cmd2
	local i tmp o1 o2 or1 or2 e1 e2

	if [[ $version -eq 2 ]]; then
		for i in "${!_fq1_sortmerna[@]}"; do
			helper::basename -f "${_fq1_sortmerna[$i]}" -o o1 -e e1
			e1=$(echo $e1 | cut -d '.' -f 1)

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.sortmerna)")
			tmp="${tdirs[-1]}/$o1"

			or1="$outdir/rRNA.$o1.$e1.gz"
			o1="$outdir/$o1.$e1.gz"

			# sortmerna v2.1 input must not be compressed (v.3.* creates empty files)
			# outfile gets extension from input file
			# in.fq.bz2 > in.fq + rRNA.out|out -> rRNA.out.fq|out.fq -> rRNA.out.fq.gz|out.fq.gz
			if [[ ${_fq2_sortmerna[$i]} ]]; then
				helper::basename -f "${_fq2_sortmerna[$i]}" -o o2 -e e2
				e2=$(echo $e2 | cut -d '.' -f 1)
				or2="$outdir/rRNA.$o2.$e2.gz"
				o2="$outdir/$o2.$e2.gz"

				commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- CMD
					paste <(helper::cat -f "${_fq1_sortmerna[$i]}") <(helper::cat -f "${_fq2_sortmerna[$i]}") | paste - - - - |
				CMD
					awk -F '\t' -v OFS='\n' '{print $1,$3,$5,$7; print $2,$4,$6,$8}'
				CMD
					> "$tmp.merged.$e1";
				CMD
					exec 11>&1; exec 12>&1;
				CMD
					ln -sfn /dev/fd/11 "$tmp.ok.$e1"; ln -sfn /dev/fd/12 "$tmp.rRNA.$e1";
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
					| cat;
				CMD
					rm -f "$tmp.merged.$e1"; exec 11>&-; exec 12>&-
				CMD
				# attention: sometimes sortmerna inserts empty lines - use sed /^\s*$/d

				_fq1_sortmerna[$i]="$o1"
				_fq2_sortmerna[$i]="$o2"
			else
				helper::makecatcmd -a catcmd -f "${_fq1_sortmerna[$i]}"

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD
					if [[ "$catcmd" == "cat" ]]; then ln -sfn "${_fq1_sortmerna[$i]}" "$tmp.$e1"; else helper::cat -f "${_fq1_sortmerna[$i]}" > "$tmp.$e1"; fi
				CMD
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
		else
			commander::runcmd -c sortmerna -v -b -i 1 -a cmd1
		fi
	else
		# NEW
		# comparison for 150m 51bp reads @ 56 cores:
		# old: 60 min
		# new without tweak and fast ref: 70 min
		# new with tweak and fast ref: 20 min (65m rRNA)
		# new with tweak and default or old ref: 30 min (65m or 67m rRNA)
		for i in "${!_fq1_sortmerna[@]}"; do
			helper::basename -f "${_fq1_sortmerna[$i]}" -o o1 -e e1
			e1=$(echo $e1 | cut -d '.' -f 1)
			or1="$outdir/rRNA.$o1.$e1.gz"
			o1="$outdir/$o1.$e1.gz"

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.sortmerna)")

			if [[ ${_fq2_sortmerna[$i]} ]]; then
				helper::basename -f "${_fq2_sortmerna[$i]}" -o o2 -e e2
				e2=$(echo $e2 | cut -d '.' -f 1)
				or2="$outdir/rRNA.$o2.$e2.gz"
				o2="$outdir/$o2.$e2.gz"

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					helper::index -f "${_fq1_sortmerna[$i]}"
				CMD
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					helper::index -f "${_fq2_sortmerna[$i]}"
				CMD

				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- CMD {COMMANDER[7]}<<- CMD {COMMANDER[8]}<<- CMD {COMMANDER[9]}<<- CMD {COMMANDER[10]}<<- CMD {COMMANDER[11]}<<- CMD {COMMANDER[12]}<<- CMD {COMMANDER[13]}<<- CMD {COMMANDER[14]}<<- CMD
					helper::cat -f "${_fq1_sortmerna[$i]}" | head -n $((threads*4)) | helper::pgzip -t $threads -o "${tdirs[-1]}/input.R1.fastq.gz"
				CMD
					helper::cat -f "${_fq2_sortmerna[$i]}" | head -n $((threads*4)) | helper::pgzip -t $threads -o "${tdirs[-1]}/input.R2.fastq.gz"
				CMD
					sortmerna
						--idx-dir "\$CONDA_PREFIX/rRNA_databases/index"
						--ref "\$CONDA_PREFIX/rRNA_databases/smr_v4.3_fast_db.fasta"
						--workdir "${tdirs[-1]}"
						--threads $threads
						--reads "${tdirs[-1]}/input.R1.fastq.gz"
						--reads "${tdirs[-1]}/input.R2.fastq.gz"
						--fastx
						--no-best
						--num_alignments 1
						--aligned
						--other
						--zip-out true
						--paired_out
						--out2
						--task 1
				CMD
					rm -rf "${tdirs[-1]}/kvdb" "${tdirs[-1]}/out/"*
				CMD
					helper::capply -i $threads -f "${_fq1_sortmerna[$i]}" -c 'bgzip -@ 1 -kc > "${tdirs[-1]}/readb/fwd_\$((JOB_ID-1)).fq.gz"'
				CMD
					helper::capply -i $threads -f "${_fq2_sortmerna[$i]}" -c 'bgzip -@ 1 -kc > "${tdirs[-1]}/readb/rev_\$((JOB_ID-1)).fq.gz"'
				CMD
					mkfifo "${tdirs[-1]}/out/aligned_fwd_0.fq.gz" "${tdirs[-1]}/out/aligned_rev_0.fq.gz" "${tdirs[-1]}/out/other_fwd_0.fq.gz" "${tdirs[-1]}/out/other_rev_0.fq.gz"
				CMD
					{ cat "${tdirs[-1]}/out/aligned_fwd_0.fq.gz" | helper::index -o "$or1" & } 2> /dev/null
				CMD
					{ cat "${tdirs[-1]}/out/aligned_rev_0.fq.gz" | helper::index -o "$or2" & } 2> /dev/null
				CMD
					{ cat "${tdirs[-1]}/out/other_fwd_0.fq.gz" | helper::index -o "$o1" & } 2> /dev/null
				CMD
					{ cat "${tdirs[-1]}/out/other_rev_0.fq.gz" | helper::index -o "$o2" & } 2> /dev/null
				CMD
					exec {FDAF}<> "${tdirs[-1]}/out/aligned_fwd_0.fq.gz"; exec {FDAR}<> "${tdirs[-1]}/out/aligned_rev_0.fq.gz"; exec {FDOF}<> "${tdirs[-1]}/out/other_fwd_0.fq.gz"; exec {FDOR}<> "${tdirs[-1]}/out/other_rev_0.fq.gz"
				CMD
					sortmerna
						--idx-dir "\$CONDA_PREFIX/rRNA_databases/index"
						--ref "\$CONDA_PREFIX/rRNA_databases/smr_v4.3_fast_db.fasta"
						--workdir "${tdirs[-1]}"
						--threads $threads
						--reads "${tdirs[-1]}/input.R1.fastq.gz"
						--reads "${tdirs[-1]}/input.R2.fastq.gz"
						--fastx
						--no-best
						--num_alignments 1
						--aligned
						--other
						--zip-out true
						--paired_out
						--out2
				CMD
					exec {FDAF}>&-; exec {FDAR}>&-; exec {FDOF}>&-; exec {FDOR}>&-
				CMD
					rm -rf "${tdirs[-1]}"
				CMD
				# needs wait after closing FDs?
				# { cat "${tdirs[-1]}/out/aligned_fwd_0.fq.gz" | helper::index -o "$or1" & } 2> /dev/null
				# -> unforunately, we have to de -and re-compress the data to remove intermediate gzip headers, that cause segemehl to process only the first chunk

				_fq1_sortmerna[$i]="$o1"
				_fq2_sortmerna[$i]="$o2"
			else
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					helper::index -f "${_fq1_sortmerna[$i]}"
				CMD

				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD {COMMANDER[6]}<<- CMD {COMMANDER[7]}<<- CMD {COMMANDER[8]}<<- CMD {COMMANDER[9]}<<- CMD {COMMANDER[10]}<<- CMD
					helper::cat -f "${_fq1_sortmerna[$i]}" | head -n $((threads*4)) | helper::pgzip -t $threads -o "${tdirs[-1]}/input.fastq.gz"
				CMD
					sortmerna
						--idx-dir "\$CONDA_PREFIX/rRNA_databases/index"
						--ref "\$CONDA_PREFIX/rRNA_databases/smr_v4.3_fast_db.fasta"
						--workdir "${tdirs[-1]}"
						--threads $threads
						--reads "${tdirs[-1]}/input.fastq.gz"
						--fastx
						--no-best
						--num_alignments 1
						--aligned
						--other
						--zip-out true
						--task 1
				CMD
					rm -rf "${tdirs[-1]}/kvdb" "${tdirs[-1]}/out/"*
				CMD
					helper::capply -i $threads -f "${_fq1_sortmerna[$i]}" -c 'bgzip -@ 1 -kc > "${tdirs[-1]}/readb/fwd_\$((JOB_ID-1)).fq.gz"'
				CMD
					mkfifo "${tdirs[-1]}/out/aligned_0.fq.gz" "${tdirs[-1]}/out/other_0.fq.gz"
				CMD
					{ cat "${tdirs[-1]}/out/aligned_0.fq.gz" | helper::index -o "$or1" & } 2> /dev/null
				CMD
					{ cat "${tdirs[-1]}/out/other_0.fq.gz" | helper::index -o "$o1" & } 2> /dev/null
				CMD
					exec {FDA}<> "${tdirs[-1]}/out/aligned_0.fq.gz"; exec {FDO}<> "${tdirs[-1]}/out/other_0.fq.gz"
				CMD
					sortmerna
						--idx-dir "\$CONDA_PREFIX/rRNA_databases/index"
						--ref "\$CONDA_PREFIX/rRNA_databases/smr_v4.3_fast_db.fasta"
						--workdir "${tdirs[-1]}"
						--threads $threads
						--reads "${tdirs[-1]}/input.fastq.gz"
						--fastx
						--no-best
						--num_alignments 1
						--aligned
						--other
						--zip-out true
				CMD
					exec {FDA}>&-; exec {FDO}>&-
				CMD
					rm -rf "${tdirs[-1]}"
				CMD
				# needs wait after closing FDs?

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
	fi

	return 0
}

function preprocess::add4stats(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
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
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
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
					preprocess::fastqc -S false -s $skip -t 1 -o "$qdir" -1 fq
				CMD
				# rescue enables to store quality directories without running fastqc after each fastq processing step - instead run fastqc on all files in parallel
			fi
		done
	done
	if [[ ${#cmdqc[@]} -gt 0 ]]; then
		commander::runcmd -b -i $threads -a cmdqc
	fi

	commander::printinfo "plotting preprocessing stats"

	local o c multiplier tool tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.tsv)"
	declare -a counts
	echo -e "sample\ttype\tcount" > "$outdir/preprocessing.barplot.tsv"
	echo -e "sample\ttype\tcount" > "$outdir/preprocessing.barplot.overlayed.tsv"
	for i in "${!_fq1_qcstats[@]}"; do
		helper::basename -f "${_fq1_qcstats[$i]}" -o b -e e
		e=$(echo $e | cut -d '.' -f 1) # if e == fastq or fq : check for ${b}_fastqc.zip else $b.${e}_fastqc.zip
		[[ $e == "fastq" || $e == "fq" ]] && qczip="${b}_fastqc.zip" || qczip="$b.${e}_fastqc.zip"
		o="$outdir/$b.stats"
		multiplier=1
		[[ "${_fq2_qcstats[$i]}" ]] && multiplier=2
		rm -f "$o" "$tmp"
		for qdir in "${_qualdirs_qcstats[@]}"; do
			tool=$(basename "$qdir")
			c=$(unzip -c "$qdir/$qczip" "${qczip%.*}/fastqc_data.txt" | grep -m 1 -F Total | awk -v mult=$multiplier '{print $3*mult}')
			counts+=($c)
			echo -e "$b\t$tool reads\t$c" >> $o
			perl -sle 'print join"\t",("$sample ($all)","$tool reads",(100*$c/$all))' -- -all=${counts[$((i*${#_qualdirs_qcstats[@]}))]} -c=$c -sample=$b -tool=$tool >> "$tmp"
			perl -sle 'print join"\t",("$sample ($all)","$tool reads",(100*$c/$all))' -- -all=${counts[$((i*${#_qualdirs_qcstats[@]}))]} -c=$c -sample=$b -tool=$tool
		done >> "$outdir/preprocessing.barplot.tsv"
		tac "$tmp" | awk -F '\t' '{OFS="\t"; if(c){$NF=$NF-c} c=c+$NF; print}' | tac >> "$outdir/preprocessing.barplot.overlayed.tsv"
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
			ggplot(m, aes(x = sample, y = count, fill = type)) +
				ggtitle("Preprocessing") + xlab("Sample") + ylab("Readcount") +
				theme_bw() + guides(fill=guide_legend(title=NULL)) +
				theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
				geom_bar(position = "fill", stat = "identity") +
				scale_y_continuous(breaks = pretty_breaks(), labels = percent_format());
			suppressMessages(ggsave(outfile));
		'
	CMD
		"$outdir/preprocessing.barplot.overlayed.tsv"  "$outdir/preprocessing.barplot.overlayed.pdf"
	CMD

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
			ggplot(m, aes(x = sample, y = count, fill = type)) +
				ggtitle("Preprocessing") + xlab("Sample") + ylab("Readcount") +
				theme_bw() + guides(fill=guide_legend(title=NULL)) +
				theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5, size = 8)) +
				geom_bar(stat = "identity", position = "dodge") +
				scale_y_continuous(breaks = pretty_breaks(), labels = percent_format(scale = 1));
			suppressMessages(ggsave(outfile));
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
