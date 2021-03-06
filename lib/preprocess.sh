#! /usr/bin/env bash
# (c) Konstantin Riege

preprocess::fastqc() {
	declare -a tdirs
	_cleanup::preprocess::fastqc(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-o <outdir>   | path to
			-p <tmpdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir
	declare -n _fq1_fastqc _fq2_fastqc
	while getopts 'S:s:t:p:o:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1) ((++mandatory)); _fq1_fastqc=$OPTARG;;
			2) _fq2_fastqc=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "calculating qualities"

	declare -a cmd1
	local f
	for f in {"${_fq1_fastqc[@]}","${_fq2_fastqc[@]}"}; do
		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.fastqc)")
		commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			fastqc
			-d "${tdirs[-1]}"
			-outdir "$outdir"
			"$f" 2>&1
		CMD
			sed -u '/Exception/{q 1};${/Analysis complete/!{q 1}}'
		CMD
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c fastqc -v -b -t $threads -a cmd1
	fi

	return 0
}

preprocess::rmpolynt(){
	_usage() {
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
		poly+=("$(printf "$i%.0s" {1..100})X")
	done
	if $dinuc; then
		for i in AB CD GH TV; do #iupac
			poly+=("$(printf "$i%.0s" {1..100})X")
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

preprocess::cutadapt(){
	_usage() {
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
				-o >(bgzip -@ $(((threads+1)/2)) -c > "$o1")
				-p >(bgzip -@ $(((threads+1)/2)) -c > "$o2")
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
				-o >(bgzip -@ $threads -c > "$o1")
				"${_fq1_cutadapt[$i]}"
				| cat
			CMD
			_fq1_cutadapt[$i]="$o1"
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c cutadapt -v -b -t 1 -a cmd1
	fi

	return 0
}

preprocess::trimmomatic() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-b <rrbs>     | true/false
			-t <threads>  | number of
			-o <outdir>   | path to
			-p <tmpdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir rrbs=false
	declare -n _fq1_trimmomatic _fq2_trimmomatic
	while getopts 'S:s:t:b:m:o:p:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			b) rrbs=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			1) ((++mandatory)); _fq1_trimmomatic=$OPTARG;;
			2) _fq2_trimmomatic=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

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
		helper::makecatcmd -c catcmd -f $f
		commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
			$catcmd $f | head -4000
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
	mapfile -t mapdata < <(commander::runcmd -t $threads -a cmd1)
	for l in "${mapdata[@]}"; do
		a=($l)
		phred["${a[@]:1}"]="${a[0]}"
	done

	local instances ithreads jmem jgct jcgct
	read -r instances ithreads jmem jgct jcgct < <(configure::jvm -i 1 -T $threads)

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
				>(bgzip -@ $(((threads+1)/2)) -c > "$o1") >(bgzip -@ $(((threads+1)/2)) -c > "$os1")
				>(bgzip -@ $(((threads+1)/2)) -c > "$o2") >(bgzip -@ $(((threads+1)/2)) -c > "$os2")
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
				>(bgzip -@ $threads -c > "$o1")
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
		commander::runcmd -v -b -t 1 -a cmd2
	fi

	return 0
}

preprocess::rcorrector(){
	declare -a tdirs
	_cleanup::preprocess::rcorrector(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-o <outdir>   | path to
			-p <tmpdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir
	declare -n _fq1_rcorrector _fq2_rcorrector
	while getopts 'S:s:t:o:p:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			1) ((++mandatory)); _fq1_rcorrector=$OPTARG;;
			2) _fq2_rcorrector=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "correcting read errors"

	declare -a cmd1 cmd2
	local i o1 e1 o2 e2
	for i in "${!_fq1_rcorrector[@]}"; do
		helper::basename -f "${_fq1_rcorrector[$i]}" -o o1 -e e1
		e1=$(echo $e1 | cut -d '.' -f 1)
		o1="$outdir/$o1"
		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.rcorrector)")
		if [[ ${_fq2_rcorrector[$i]} ]]; then
			helper::basename -f "${_fq2_rcorrector[$i]}" -o o2 -e e2
			e2=$(echo $e2 | cut -d '.' -f 1)
			o2="$outdir/$o2"

			readlink -e "${_fq1_rcorrector[$i]}" | file -f - | grep -qF 'compressed' && {
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					cd "${tdirs[-1]}"
				CMD
					run_rcorrector.pl
					-1 "${_fq1_rcorrector[$i]}"
					-2 "${_fq2_rcorrector[$i]}"
					-od "$outdir"
					-t $threads
				CMD
					mv "$o1.cor.fq.gz" "$o1.$e1.gz"
				CMD
					mv "$o2.cor.fq.gz" "$o2.$e2.gz"
				CMD
			} || {
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD {COMMANDER[5]}<<- CMD
					cd "${tdirs[-1]}"
				CMD
					exec 11>&1; exec 12>&1
				CMD
					ln -sfn "/dev/fd/11" "$(basename "$o1").cor.fq"
				CMD
					ln -sfn "/dev/fd/12" "$(basename "$o2").cor.fq"
				CMD
					run_rcorrector.pl
					-1 "${_fq1_rcorrector[$i]}"
					-2 "${_fq2_rcorrector[$i]}"
					-od "${tdirs[-1]}"
					-t $threads
					11> >(bgzip -@ $(((threads+1)/2)) -c > "$o1.$e1.gz")
					12> >(bgzip -@ $(((threads+1)/2)) -c > "$o2.$e2.gz")
					| cat
				CMD
					exec 11>&-; exec 12>&-
				CMD
			}

			_fq1_rcorrector[$i]="$o1.$e1.gz"
			_fq2_rcorrector[$i]="$o2.$e2.gz"
		else
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cd "${tdirs[-1]}"
			CMD
				run_rcorrector.pl
				-s "${_fq1_rcorrector[$i]}"
				-stdout
				-t $threads
				> >(bgzip -@ $threads -c > "$o1.$e1.gz")
				| cat
			CMD
			_fq1_rcorrector[$i]="$o1.$e1.gz"
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c rcorrector -v -b -t 1 -a cmd1
	fi

	return 0
}

preprocess::sortmerna(){
	declare -a tdirs
	_cleanup::preprocess::sortmerna(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-o <outdir>   | path to
			-p <tmpdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir
	declare -n _fq1_sortmerna _fq2_sortmerna
	while getopts 'S:s:t:i:o:p:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			1) ((++mandatory)); _fq1_sortmerna=$OPTARG;;
			2) _fq2_sortmerna=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "filtering rRNA fragments"

	local insdir=$(dirname $(dirname $(which sortmerna)))
	local sortmernaref=$(for i in $insdir/rRNA_databases/*.fasta; do echo $i,$insdir/index/$(basename $i .fasta)-L18; done | xargs -echo | sed 's/ /:/g')

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
				paste - - - - - - - -
			CMD
				awk -F '\t' -v OFS='\n' '{print $1,$3,$5,$7; print $2,$4,$6,$8}'
			CMD

			# alternative:
			# 11> >(sed -E '/^\s*$/d' | paste - - - - | awk -F '\t' -v OFS='\n' '{if(NR%2==1){print \$1,\$2,\$3,\$4 > "/dev/fd/1"}else{print \$1,\$2,\$3,\$4 > "/dev/fd/2"}}' > >(bgzip -@ $(((threads+1)/2)) -c > "$o1") 2> >(bgzip -@ $(((threads+1)/2)) -c > "$o2"))
			# 12> >(sed -E '/^\s*$/d' | paste - - - - | awk -F '\t' -v OFS='\n' '{if(NR%2==1){print \$1,\$2,\$3,\$4 > "/dev/fd/1"}else{print \$1,\$2,\$3,\$4 > "/dev/fd/2"}}' > >(bgzip -@ $(((threads+1)/2)) -c > "$or1") 2> >(bgzip -@ $(((threads+1)/2)) -c > "$or2"))
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD
				exec 11>&1; exec 12>&1
			CMD
				ln -sfn /dev/fd/11 "$tmp.ok.$e1"
			CMD
				ln -sfn /dev/fd/12 "$tmp.rRNA.$e1"
			CMD
				sortmerna
				--ref "$sortmernaref"
				--reads "$tmp.merged.$e1"
				--fastx
				--paired_out
				--aligned "$tmp.rRNA"
				--other "$tmp.ok"
				-a $threads
				11> >(sed -E '/^\s*$/d' | paste - - - - | tee >(sed -n '1~2{s/\t/\n/gp}' | bgzip -@ $(((threads+1)/2)) -c > "$o1") >(sed -n '2~2{s/\t/\n/gp}' | bgzip -@ $(((threads+1)/2)) -c > "$o2") > /dev/null)
				12> >(sed -E '/^\s*$/d' | paste - - - - | tee >(sed -n '1~2{s/\t/\n/gp}' | bgzip -@ $(((threads+1)/2)) -c > "$or1") >(sed -n '2~2{s/\t/\n/gp}' | bgzip -@ $(((threads+1)/2)) -c > "$or2") > /dev/null)
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
				--ref "$sortmernaref"
				--reads "$tmp.$e1"
				--fastx
				--aligned "$tmp.rRNA"
				--other "$tmp.ok"
				-a $threads
				11> >(sed -E '/^\s*$/d' | bgzip -@ $threads -c > "$o1")
				12> >(sed -E '/^\s*$/d' | bgzip -@ $threads -c > "$or1")
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
		commander::runcmd -v -b -t $threads -a cmd1
		commander::runcmd -v -b -t 1 -a cmd2
	fi

	return 0
}

preprocess::add4stats(){
	_usage() {
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

preprocess::qcstats(){
	local tmp
	_cleanup::preprocess::qcstats(){
		rm -f "$tmp"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-f <force>    | true/false rerun fastqc
			-t <threads>  | number of
			-i <qualdirs> | array of
			-o <outdir>   | path to
			-p <tmpdir>   | path to
			-1 <fastq1>   | array of
			-2 <fastq2>   | array of
		EOF
		return 1
	}

	local OPTIND arg mandatory threads skip=false outdir tmpdir force=false
	declare -n _qualdirs_qcstats _fq1_qcstats _fq2_qcstats
	while getopts 'S:s:t:f:i:o:p:1:2:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			f) force=$OPTARG;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); _qualdirs_qcstats=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			1) ((++mandatory)); _fq1_qcstats=$OPTARG;;
			2) _fq2_qcstats=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	declare -a cmdqc
	local qdir i f b e
	for i in "${!_qualdirs_qcstats[@]}"; do
		qdir="${_qualdirs_qcstats[$i]}"
		declare -n _pi_qcstats=preprocess$i
		for f in "${_pi_qcstats[@]}"; do
			helper::basename -f "$f" -o b -e e
			if $force || [[ ! -s "$qdir/${b}_fastqc.zip" ]]; then
				commander::makecmd -a cmdqc -s ' ' -c {COMMANDER[0]}<<- CMD
					fq=("$f");
					preprocess::fastqc -S false -s $skip -t 1 -p "$tmpdir" -o "$qdir" -1 fq
				CMD
				# rescue enables to store quality directories without running fastqc after each fastq processing step - instead run fastqc on all files in parallel
			fi
		done
	done
	if [[ ${#cmdqc[@]} -gt 0 ]]; then
		commander::runcmd -b -t $threads -a cmdqc
	fi

	commander::printinfo "plotting preprocessing stats"

	tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.tsv)"
	local o c multiplier tool
	declare -a counts
	echo -e "sample\ttype\tcount" > "$outdir/preprocessing.barplot.tsv"
	for i in "${!_fq1_qcstats[@]}"; do
		helper::basename -f "${_fq1_qcstats[$i]}" -o b -e e
		o="$outdir/$b.stats"
		multiplier=1
		[[ "${_fq2_qcstats[$i]}" ]] && multiplier=2
		rm -f $o
		for qdir in "${_qualdirs_qcstats[@]}"; do
			tool=$(basename "$qdir")
			c=$(unzip -p "$qdir/${b}_fastqc.zip" "${b}_fastqc/fastqc_data.txt" | grep -m 1 -F Total | awk -v mult=$multiplier '{print $3*mult}')
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
			m <- read.table(intsv, header=T, sep="\t");
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
		commander::runcmd -v -b -t $threads -a cmd1
	fi

	return 0
}
