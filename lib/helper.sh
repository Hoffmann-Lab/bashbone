#!/usr/bin/env bash
# (c) Konstantin Riege

function helper::index(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::index indexes a flat file or gzip compressed file for random access data retrieval by line number (see helper::cat, helper::makecatcmd)

			-f <infile>  | somewhat optional (see -o). default: stdin
			             | path to input file. unless given, -o is mandatory
			-o <outfile> | somewhat optional (see -f). default: stdout
			             | path to output file. unless given, -f is mandatory
		EOF
		return 1
	}

	local OPTIND arg f o
	while getopts 'f:o:' arg; do
		case $arg in
			f)	f="$OPTARG";;
			o)	o="$OPTARG"; mkdir -p "$(dirname "$o")";;
			*)	_usage;;
		esac
	done
	[[ ! $f && ! $o ]] && { _usage || return 0; }

	if [[ $f ]]; then
		if readlink -e "$f" | file -b --mime-type -f - | grep -qF 'gzip'; then
			if rapidgzip --version 2> /dev/null | sed -E 's/.+\s([0-9]+)\.([0-9]+)\.([0-9]+)$/\1.\2\t\2.\3/' | awk '$1>0.14 || $2>=14.3{exit 0}{exit 1}'; then
				# test if available index is okay
				# if [[ -e "${f%.*}.gzi" ]]; then
				# 	rapidgzip -qkcd -P 1 --import-index "${f%.*}.gzi" "$f" 2> /dev/null | head -1 | grep -q . && return 0 || true
				# fi
				rapidgzip -qf -P 4 --index-format gztool-with-lines --export-index "${f%.*}.gzi" "$f"
			else
				# test if available index is okay
				# if [[ -e "${f%.*}.rgzi" ]]; then
				# 	gztool -v 0 -l "$f" && rapidgzip -kcd -P 2 --import-index "${f%.*}.rgzi" "$f" 2> /dev/null | head -1 | grep -q . && return 0
				# fi
				cat "$f" | tee -i >(gztool -v 0 -f -i -x -C -I "${f%.*}.gzi") >(rapidgzip -qf -P 2 --export-index "${f%.*}.rgzi") > /dev/null | cat
			fi
		else
			fftool.sh -i -f "$f"
		fi
	else
		local magic recordsize=1000
		read -n 2 magic
		[[ $(wc -c <<< $magic) -eq 3 ]] || magic+=$'\n';
		# diff <(echo -n $magic | od -N 2 -t x1) <(echo | gzip -kc | od -N 2 -t x1) &> /dev/null # 0x1f8b
		if echo -n "$magic" | file -b --mime-type - | grep -qF 'gzip'; then
			if rapidgzip --version 2> /dev/null | sed -E 's/.+\s([0-9]+)\.([0-9]+)\.([0-9]+)$/\1.\2\t\2.\3/' | awk '$1>0.14 || $2>=14.3{exit 0}{exit 1}'; then
				# as of version 14.3 indexing from stdin still does not work for gztool index - requires seekable file
				# cat <(echo -n "$magic") - | tee -i >(rapidgzip -qf -P 1 --index-format gztool-with-lines --export-index "${f%.*}.gzi") > "$o" | cat
				cat <(echo -n "$magic") - | tee -i >(gztool -v 0 -f -i -x -C -I "${o%.*}.gzi") > "$o" | cat
			else
				cat <(echo -n "$magic") - | tee -i >(gztool -v 0 -f -i -x -C -I "${o%.*}.gzi") >(rapidgzip -f -P 2 --export-index "${o%.*}.rgzi") > "$o" | cat
			fi
		else
			cat <(echo -n "$magic") - | fftool.sh -i -o "$o"
		fi
	fi

	return 0
}

function helper::pgzip(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::pgzip compresses a file into gzip format while indexing it on the fly for random access data retrieval by line number (see helper::cat, helper::makecatcmd)

			-f <infile>  | somewhat optional (see -o). default: stdin
			             | path to input file. unless given, -o is mandatory
			-o <outfile> | somewhat optional (see -f). default: stdout
			             | path to output file. unless given, -f is mandatory
			-t <threads> | optional.
			             | number of threads
			-p           | optional.
			             | switch from bgzip to pigz to not compress block wise at the cost of higher runtime, but lower disk footprint (5%)
		EOF
		return 1
	}

	local OPTIND arg threads=1 f o tool="bgzip"
	while getopts 'f:t:o:p' arg; do
		case $arg in
			f)	f="$OPTARG";;
			o)	o="$OPTARG"; mkdir -p "$(dirname "$o")";;
			t)	threads=$OPTARG;;
			p)	tool="pigz";;
			*)	_usage;;
		esac
	done
	[[ ! $f && ! $o ]] && { _usage || return 0; }
	[[ ! $f ]] && f=/dev/stdin
	[[ ! $o ]] && o="$f.gz"

	if [[ $tool == "bgzip" && -x "$(command -v bgzip)" ]]; then
		# compression level changed from v1.16 to ~match gzip size at cost of throughput (still being faster than pigz)
		tool+=" -l $(bgzip --version | head -1 | awk '$NF>1.15{print 5; exit}{print 6}') -@"
	else
		tool="pigz -p"
	fi

	$tool $threads -k -c "$f" | helper::index -o "$o"

	return 0
}

function helper::apply(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::apply executes commands in parallel

			-t <threads> | optional. default: 1
			             | number of threads
			-k           | optional.
			             | keep ordered output by fully buffering output stream instead of immediate line bufferd output
			-c <cmd>     | optional.
			             | ALWAYS LAST OPTION. command or exported function and parameters to apply

			example 1:
			{ echo "echo foo"; echo "echo bar"; } | helper::apply -t 10

			example 2:
			fun(){
				echo $JOB_ID $*
			}
			export -f fun
			{ echo foo; echo bar; } | helper::apply -t 10 -c fun {}

			example 3:
			{ echo foo; echo bar; } | helper::apply -t 10 -c echo \\\$JOB_ID {}

			example 4:
			helper::apply -t 10 -c echo \$JOB_ID {} ::: foo bar

			example 5:
			helper::apply -t 10 -c 'echo $JOB_ID {} ::: foo bar'
		EOF
		return 1
	}

	local OPTIND arg mandatory threads=1 params="--line-buffer" tmpdir="${TMPDIR:-/tmp}"
	while getopts 't:kc' arg; do
		case $arg in
			t)	threads=$OPTARG;;
			k)	params="--keep-order";;
			c)	break;;
			*)	_usage;;
		esac
	done
	shift $((OPTIND-1))
	[[ $1 ]] || [[ ! -t 0 ]] || _usage

	local shell="$(mktemp -p "$tmpdir" shell.XXXXXXXXXX)"
	cat <<- 'EOF' > "$shell"
		#!/usr/bin/env bash
		shift
		JOB_ID=$PARALLEL_SEQ
		cenv="$(basename "$CONDA_PREFIX")"
		source "$BASHBONE_DIR/activate.sh" -l ${BASHBONE_LEGACY:-true} -s "$BASHBONE_EXTENSIONDIR" -c ${BASHBONE_CONDA:-false} -r false -i "$BASHBONE_TOOLSDIR"
		if ${BASHBONE_CONDA:-false} && [[ $cenv != "$(basename "$CONDA_PREFIX")" ]]; then
		    conda activate --no-stack $cenv
		fi
		eval "$*"
	EOF
	chmod 755 "$shell"
	PARALLEL_SHELL="$shell" parallel --tmpdir "$tmpdir" --termseq INT,1000,TERM,0 --halt now,fail=1 $params -P $threads $*

	return 0
}

function helper::lapply(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::lapply applies a command or function on each line of input in parallel

			-l <lines>   | optional. default: 1
			             | number of lines that make a record
			-d <records> | optional. default: 5000
			             | number of records processed as chunk in a thread
			-t <threads> | optional. default: 1
			             | number of threads
			-o <odir>    | optional. default: TMPDIR
			             | path to output directory of temporary files in order to keep them (see -f)
			-f           | optional.
			             | create seekable temporary files. file name is $FILE
			-r           | optional. default: initiate new job per chunk (job id is $JOB_ID)
			             | switch to round robin mode i.e. distribute all chunks evenly across once only spawned threads
			-k           | optional.
			             | keep ordered output by fully buffering output stream instead of immediate line bufferd output. somewhat mutally exclusive to -r
			-c <cmd>     | mandatory.
			             | ALWAYS LAST OPTION. command or exported function and parameters to apply

			example 1:
			cat file | helper::lapply -t 10 -c cat

			example 2:
			fun(){
				cat > $JOB_ID.out
			}
			export -f fun
			cat file | helper::lapply -t 10 -c fun

			example 3:
			cat file | helper::lapply -f -t 10 -c "cat \$FILE > \$FILE.out"

			example 4:
			cat file | helper::lapply -f -t 10 -c 'cat > $FILE.out'
		EOF
		return 1
	}

	local OPTIND arg mandatory threads=1 l=1 d=5000 k=false; params="" odir="${TMPDIR:-/tmp}"
	while getopts 'l:d:t:o:rfkc' arg; do
		case $arg in
			l)	l=$OPTARG;;
			d)	d=$OPTARG;;
			t)	threads=$OPTARG;;
			o)	odir="$OPTARG"; mkdir -p "$odir";;
			r)	params+=" --round-robin";;
			f)	params+=" --cat";;
			k)	k=true; params+=" --keep-order";;
			c)	((++mandatory)); shift $((OPTIND-1)); break;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage
	$k || params+=" --line-buffer"

	# attention: this solution requires escape sequences within single quotes. (example: cat file | papply -t 10 -c 'cat \> \$JOB_ID.out')
	# function helper::_papply(){
	# 	FILE="${1:-/dev/stdin}"
	# 	JOB_ID="$2" # or $PARALLEL_SEQ
	# 	shift 2
	# 	if [[ "$FILE" == "/dev/stdin" ]]; then
	# 		IFS=$'\n'
	# 		read -r l < "$FILE" || exit 0
	# 		cat <(echo "$l") - | eval "$*"
	# 	else
	# 		[[ -s "$FILE" ]] || exit 0
	# 		eval "$*"
	# 	fi
	# }
	# export -f helper::_papply
	# parallel --tmpdir "$odir" --termseq INT,1000,TERM,0 --halt now,fail=1 --line-buffer --pipe -L $l -N $((5000/l)) $params -P $threads helper::_papply "{}" "{#}" "$@"

	local shell="$(mktemp -p "${TMPDIR:-/tmp}" shell.XXXXXXXXXX)"
	cat <<- 'EOF' > "$shell"
		#!/usr/bin/env bash
		shift
		cenv="$(basename "$CONDA_PREFIX")"
		source "$BASHBONE_DIR/activate.sh" -l ${BASHBONE_LEGACY:-true} -s "$BASHBONE_EXTENSIONDIR" -c ${BASHBONE_CONDA:-false} -r false -i "$BASHBONE_TOOLSDIR"
		if ${BASHBONE_CONDA:-false} && [[ $cenv != "$(basename "$CONDA_PREFIX")" ]]; then
		    conda activate --no-stack $cenv
		fi
		eval "$*"
	EOF
	chmod 755 "$shell"
	PARALLEL_SHELL="$shell" parallel --tmpdir "$odir" --termseq INT,1000,TERM,0 --halt now,fail=1 --pipe -L $l -N $((d/l)) $params -P $threads "JOB_ID={#}; FILENO={#}; FILE=\"{}\"; if [[ -s \"\$FILE\" ]]; then $*; else read -r l || exit 0; cat <(echo \"\$l\") - | $*; fi"

	return 0
}

function helper::capply(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::capply applies a command or function on equally sized chunks of an input file in parallel by random access (see helper::index, helper::pgzip)

			-i <instances> | optional. default: 1 (job id is $JOB_ID)
			               | number of
			-t <threads>   | optional. default: 1
			               | number of decompressison threads per instance
			-f <infile>    | mandatory.
			               | path to indexed flat-file or indexed gzip compressed input file (see helper::index and helper::pgzip)
			-m <modulo>    | optional. default: 4 (to not break fastq file structure. number of records is $JOB_NR)
			               | on which per chunk record number modulo becomes zero (1 line is 1 record)
			-l             | optional.
			               | switch to immediate line bufferd output instead of keeping ordered output by fully buffering output stream
   			-c <cmd>       | mandatory.
			               | ALWAYS LAST OPTION. command or exported function and parameters to apply

			example 1:
			helper::capply -i 10 -f file[.gz] -c cat

			example 2:
			fun(){
				cat > $JOB_ID.out
			}
			export -f fun
			helper::capply -i 10 -f file[.gz] -c fun

			example 3:
			helper::capply -i 10 -f file[.gz] -c "cat > \$JOB_ID.out"

			example 4:
			helper::capply -i 10 -f file[.gz] -c 'cat > $JOB_ID.out'
		EOF
		return 1
	}

	local OPTIND arg mandatory instances=1 m=4 f tdir="${TMPDIR:-/tmp}" params="--keep-order" threads=1
	while getopts 'i:f:t:m:lc' arg; do
		case $arg in
			i)	instances=$OPTARG;;
			t)	threads=$OPTARG;;
			m)	m=$OPTARG;;
			l)	params="--line-buffer";;
			f)	((++mandatory)); f="$OPTARG";;
			c)	((++mandatory)); shift $((OPTIND-1)); break;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	local l c x
	if readlink -e "$f" | file -b --mime-type -f - | grep -qF 'gzip'; then
		# count-lines does not use index yet
		# if rapidgzip --version 2> /dev/null | sed -E 's/.+\s([0-9]+)\.([0-9]+)\.([0-9]+)$/\1.\2\t\2.\3/' | awk '$1>0.14 || $2>=14.3{exit 0}{exit 1}'; then
		# 	l=$(rapidgzip -P 1 -q --count-lines --import-index "${f%.*}.gzi" "$f")
		# else
			l=$(gztool -l "$f" |& sed -nE 's/.*\s+lines\s+:\s+([0-9]+).*/\1/p')
		# fi
	else
		l=$(tail -1 "${f%.*}.ffi" | cut -f 1)
	fi
	c=$((l%instances==0 ? l/instances : l/instances+1))
	c=$((c%m == 0 ? c : c+m-c%m))

	# attention: this solution requires escape sequences within single quotes. (example: papply_gzip -i 10 -f file.gz -c 'cat \> \$JOB_ID.out')
	# function helper::_papply_gzip(){
	# 	JOB_ID="$1"
	# 	shift
	# 	eval "$*"
	# }
	# export -f helper::_papply_gzip
	# for i in $(seq 0 $((instances-1))); do
	# 	echo "gztool -v 0 -L $((i*c+1)) -R $c '$f' |"
	# done | parallel --tmpdir "$tdir" --termseq INT,1000,TERM,0 --halt now,fail=1 $params -P $instances helper::_papply_gzip "{#}" "{}" "$@"

	# for i in $(seq 0 $((instances-1))); do
	# 	echo "gztool -v 0 -L $((i*c+1)) -R $c '$f'"
	# done | parallel --tmpdir "$tdir" --termseq INT,1000,TERM,0 --halt now,fail=1 --keep-order -P $instances '{= $ENV{JOB_ID}=seq() =}' :::: - ::: " | $*"
	# for i in $(seq 0 $((instances-1))); do
	# 	echo "FILENO=\$PARALLEL_SEQ; JOB_ID=\$PARALLEL_SEQ; gztool -v 0 -L $((i*c+1)) -R $c '$f'"
	# done | parallel --tmpdir "$tdir" --termseq INT,1000,TERM,0 --halt now,fail=1 --keep-order -P $instances {} :::: - ::: " | $*"
	# v0.14+:
	# for i in $(seq 0 $((instances-1))); do
	# 	echo "FILENO=\$PARALLEL_SEQ; JOB_ID=\$PARALLEL_SEQ; rapidgzip -P 1 -kcd --ranges ${c}L@$((i*c))L --import-index '${f%.*}.gzi' '$f'"
	# done | parallel --tmpdir "$tdir" --termseq INT,1000,TERM,0 --halt now,fail=1 --keep-order -P $instances {} :::: - ::: " | $*"

	local shell="$(mktemp -p "${TMPDIR:-/tmp}" shell.XXXXXXXXXX)"
	cat <<- 'EOF' > "$shell"
		#!/usr/bin/env bash
		shift
		cenv="$(basename "$CONDA_PREFIX")"
		source "$BASHBONE_DIR/activate.sh" -l ${BASHBONE_LEGACY:-true} -s "$BASHBONE_EXTENSIONDIR" -c ${BASHBONE_CONDA:-false} -r false -i "$BASHBONE_TOOLSDIR"
		if ${BASHBONE_CONDA:-false} && [[ $cenv != "$(basename "$CONDA_PREFIX")" ]]; then
		    conda activate --no-stack $cenv
		fi
		eval "$*"
	EOF
	chmod 755 "$shell"

	declare -a catcmd
	for i in $(seq 0 $((instances-1))); do
		x=$((l-i*c > c ? c : l-i*c)) # -r "${c}L@$((i*c))L" also works, but sometimes it is necessary to know the exact number of records in the last chunk
		helper::makecatcmd -a catcmd -t $threads -r "${x}L@$((i*c))L" -f "$f"
		echo "JOB_NR=$x; FILENO=\$PARALLEL_SEQ; JOB_ID=\$PARALLEL_SEQ; ${catcmd[*]} '$f'"
	done | PARALLEL_SHELL="$shell" parallel --tmpdir "$tdir" --termseq INT,1000,TERM,0 --halt now,fail=1 $params -P $instances {} :::: - ::: " | $*"

	return 0
}

function helper::sort(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::sort wraps GNU sort in a parallelized fashion as drop-in replacement

			-f <infile>    | optional. default: stdin
			               | path to input file
			-o <outfile>   | optional. default: stdout
			               | path to output file
			-t <threads>   | mandatory
			               | number of threads
			-M <maxmemory> | optional. default: all available
			               | amount of memory to allocate
		EOF
		return 1
	}

	local OPTIND threads=1 f=/dev/stdin o=/dev/stdout maxmemory tmpdir="${TMPDIR:-/tmp}"
	declare -a args=();
	while [[ $# -gt 0 ]]; do
		case "$1" in
			-M)	maxmemory=$2; shift 2;;
			-f)	f=$2; shift 2;;
			-o)	o=$2; shift 2; mkdir -p "$(dirname "$o")";;
			-t)	threads=$2; shift 2;;
			*)	args+=("$1"); shift;;
		esac
	done
	[[ $args ]] || [[ ! -t 0 ]] || _usage

	local instances tdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.sort)"
	read -r instances maxmemory < <(configure::memory_by_instances -i 1 -M "$maxmemory")
	LC_ALL=C sort --parallel="$threads" -S "${maxmemory}M" -T "$tdir" "${args[@]}" "$f" > "$o" | cat
	# fun(){ echo foo > /dev/stdout; }; echo bar > tmp; fun >> tmp. cat tmp # foo
	# fun(){ echo foo > /dev/stdout | cat; }; echo bar > tmp; fun >> tmp. cat tmp # bar foo

	return 0
}

function helper::vcfsort(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::vcfsort wraps GNU sort in a parallelized fashion to properly sort VCF files

			-f <infile>    | optional. default: stdin
			               | path to input file
			-o <outfile>   | optional. default: stdout
			               | path to output file
			-t <threads>   | mandatory
			               | number of threads
			-M <maxmemory> | optional. default: all available
			               | amount of memory to allocate
			-z             | optional.
			               | switch to compress (and index) output
		EOF
		return 1
	}

	local OPTIND threads=1 f=/dev/stdin o=/dev/stdout maxmemory tmpdir="${TMPDIR:-/tmp}" zip=false
	while getopts 'f:t:o:M:z' arg; do
		case $arg in
			f) f="$OPTARG";;
			o) o="$OPTARG"; mkdir -p "$(dirname "$o")";;
			t) threads=$OPTARG;;
			M) maxmemory=$OPTARG;;
			z) zip=true;;
			*) _usage;;
		esac
	done

	local instances maxmemory tdir="$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.vcfsort)"
	read -r instances maxmemory < <(configure::memory_by_instances -i 1 -M "$maxmemory")
	if $zip; then
		bcftools view "$f" 2> >(sed -u '/vcf_parse/{q 1}' >&2) | awk -F'\t' -v t=$threads -v m=$maxmemory -v p="$tdir" -v OFS='\t' '/^#/{if(match($0,/^##contig=<ID=(\S+),length.*/,a)){i++; c[a[1]]=i}; print; next} {$1=c[$1]; print | "LC_ALL=C sort -k1,1n -k2,2n -k4,4 -k5,5 --parallel="t" -S "m"M -T \""p"\""}' | awk -F'\t' -v OFS='\t' '/^#/{if(match($0,/^##contig=<ID=(\S+),length.*/,a)){i++; c[i]=a[1]}; print; next} {$1=c[$1]; print}' | helper::pgzip -t $threads -o "$o" | cat
		[[ -f "$o" ]] && tabix -f -p vcf "$o"
	else
		bcftools view "$f" 2> >(sed -u '/vcf_parse/{q 1}' >&2) | awk -F'\t' -v t=$threads -v m=$maxmemory -v p="$tdir" -v OFS='\t' '/^#/{if(match($0,/^##contig=<ID=(\S+),length.*/,a)){i++; c[a[1]]=i}; print; next} {$1=c[$1]; print | "LC_ALL=C sort -k1,1n -k2,2n -k4,4 -k5,5 --parallel="t" -S "m"M -T \""p"\""}' | awk -F'\t' -v OFS='\t' '/^#/{if(match($0,/.*#contig=<ID=(\S+),length.*/,a)){i++; c[i]=a[1]}; print; next} {$1=c[$1]; print}' > "$o" | cat
	fi

	return 0
}

function helper::makecatcmd(){
	function _usage(){
		commander::print <<- 'EOF'
			description
			helper::cat / helper::makecatcmd identifies the fastest way on fully decompressing gzip or bzip files.
			Indexed flat files or indexed gzip files can be accessed randomly by line number (see helper::index, helper::pgzip)

			synopsis
			helper::makecatcmd <OPTIONS> -f <file>

			options
			-f <file>    | mandatory. path to input file
			-t <threads> | optional. number of threads to use. default: 4
			-r <range>   | optional. to access indexed flat or gzip data given line number ranges.
			             | range format: <m|inf>L@<n>L means to retrieve m or infinit number of lines after n lines

			example 1
			helper::cat -f file.[txt|bz2|gz]

			example 2
			helper::cat -r 1L@100L -f file.[txt|bz2|gz]

			developer options
			-a <cmds>    | array to append crafted commands to instead of execution
			-l           | switch to legacy mode i.e. decompression without index

			example 1
			declare -a cmds
			helper::cat -1 cmds -f file.[txt|bz2|gz]
			commander::runcmd -a cmds
		EOF
		return 1
	}

	local OPTIND arg mandatory f v legacy=false threads=4 range
	while getopts 'f:a:t:r:l' arg; do
		case $arg in
			a)	v=$OPTARG;;
			t)	threads=$OPTARG;;
			r)	range="$OPTARG";;
			l)	legacy=true;;
			f)	((++mandatory)); f="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 1 ]] && _usage

	if [[ $v ]]; then
		declare -n _makecatcmd=$v
	else
		declare -a _makecatcmd
	fi

	### SSD
	# pigz p=1  : 400MB/s
	# pigz p=2  : 500MB/s
	# bgzip @=1 : 400MB/s
	# bgzip @=2 : 400MB/s
	# rapidgzip P=4 : index or bgzip : 1300MB/s | no index : 800MB/s  <-
	# rapidgzip P=8 : index or bgzip : 1500MB/s | no index : 1100MB/s
	### NFS
	# pigz p=1  : 250MB/s
	# pigz p=2  : 200MB/s
	# bgzip @=1 : 550MB/s
	# bgzip @=4 : 2500MB/s
	# bgzip @=8 : 3800MB/s
	# rapidgzip P=4               : index or bgzip : 2500MB/s | no index : 550MB/s (P=1!)
	# rapidgzip P=8               : index or bgzip : 3800MB/s | no index : 550MB/s (P=1!)
	# rapidgzip P=12              : index or bgzip : 5500MB/s | no index : 550MB/s (P=1!)
	# rapidgzip (sequential) P=4  : index or bgzip : 2500MB/s | no index : 2200MB/s       <-
	# rapidgzip (sequential) P=8  : index or bgzip : 3800MB/s | no index : 3500MB/s
	# rapidgzip (sequential) P=12 : index or bgzip : 5500MB/s | no index : 5000MB/s

	# _makecatcmd=$({ readlink -e "$f" | file -b --mime-type -f - | grep -oF -e 'gzip' -e 'bzip2' || echo cat; } | sed '/cat/!{s/gzip/pigz -p 2/; s/$/ -kcd/}')
	# _makecatcmd=$({ readlink -e "$f" | file -b --mime-type -f - | grep -oF -e 'gzip' -e 'bzip2' || echo cat; } | sed '/cat/!{s/gzip/bgzip -@ 1/; s/$/ -kcd/}')
	# _makecatcmd=$({ readlink -e "$f" | file -b --mime-type -f - | grep -oF -e 'gzip' -e 'bzip2' || echo cat; } | sed '/cat/!{s/gzip/rapidgzip -P 4/; s/$/ -kcd/}')

	# rapidgzip supports: gzip zlib raw-deflate bzip2

	case "$(readlink -e "$f" | file -b --mime-type -f - | grep -oF -e 'gzip' -e 'bzip2' || echo cat)" in
		gzip)
			if rapidgzip --version 2> /dev/null | sed -E 's/.+\s([0-9]+)\.([0-9]+)\.([0-9]+)$/\1.\2\t\2.\3/' | awk '$1>0.14 || $2>=14.3{exit 0}{exit 1}'; then
				# from v0.14.1+ rapidgzip supports gztool indices and likewise to its default indices gives all the benefits (ie. ISA-L based decoding)
				# -> no rgzi indices necesseray anymore
				# from v0.14.3+ rapidgzip supports byte or line offset based extraction
				if [[ $range ]]; then
					_makecatcmd=("rapidgzip" "-qkcd" "-P" "$threads" "--import-index" "${f%.*}.gzi" "--ranges" "$range")
				elif ! $legacy && [[ -e "${f%.*}.gzi" ]]; then
					_makecatcmd=("rapidgzip" "-qkcd" "-P" "$threads" "--import-index")
					[[ $v ]] && _makecatcmd+=("'${f%.*}.gzi'") || _makecatcmd+=("${f%.*}.gzi")
				else
					local rota=$(
						rota="$(df --output=source "$(readlink -e "$f")" | tail -1)"
						lsblk -n -o ROTA "$rota" 2> /dev/null || {
							rota="$(realpath "/dev/disk/by-label/$(rev <<< "$rota" | xargs basename | rev)")"
							lsblk -n -o ROTA "$rota" 2> /dev/null || echo 1
						}
					)
					if [[ $rota -gt 0 ]]; then
						# hdd/nfs/overlay/...
						_makecatcmd=("rapidgzip" "-qkcd" "-P" "$threads" "--io-read-method" "sequential")
					else
						# bgzip || gzip
						# params="-P $(gzip -l "$f" | tail -1 | awk '$1>0 && $2==0{print 4; exit}{print 8}')"
						_makecatcmd=("rapidgzip" "-qkcd" "-P" "$threads")
					fi
				fi
			else
				# wo avoid non mutable "[Warning] The index only has an effect for parallel decoding"
				[[ $threads -eq 1 ]] && threads=2
				if [[ $range ]]; then
					local l r
					read -r r l < <(sed -E 's/L@?/ /g' <<< "$range")
					# attention: sometimes complains about corrupt file, altough it is not?! rapidgzip --ranges works with same gztool index
					# https://github.com/circulosmeos/gztool/issues/20 solved as of version 1.7.0 i.e. actually no requirement for -p switch any more
					_makecatcmd=("gztool" "-v" "0" "-p" "-L" "$((l+1))" "-R" "$r")
				elif ! $legacy && [[ -e "${f%.*}.rgzi" ]]; then
					_makecatcmd=("rapidgzip" "-kcd" "-P" "$threads" "--import-index")
					[[ $v ]] && _makecatcmd+=("'${f%.*}.rgzi'") || _makecatcmd+=("${f%.*}.rgzi")
				else
					local rota=$(
						rota="$(df --output=source "$(readlink -e "$f")" | tail -1)"
						lsblk -n -o ROTA "$rota" 2> /dev/null || {
							rota="$(realpath "/dev/disk/by-label/$(rev <<< "$rota" | xargs basename | rev)")"
							lsblk -n -o ROTA "$rota" 2> /dev/null || echo 1
						}
					)
					if [[ $rota -gt 0 ]]; then
						# hdd/nfs/overlay/...
						_makecatcmd=("rapidgzip" "-kcd" "-P" "$threads" "--io-read-method" "sequential")
					else
						# bgzip || gzip
						# params="-P $(gzip -l "$f" | tail -1 | awk '$1>0 && $2==0{print 4; exit}{print 8}')"
						_makecatcmd=("rapidgzip" "-kcd" "-P" "$threads")
					fi
				fi
			fi
		;;
		bzip2)
			# bzip2 replaced by faster indexed_bzip2
			# rapidgzip can be used as well, but is slower
			if ibzip2 --version 2> /dev/null; then
				_makecatcmd=("ibzip2" "-qkcd" "-P" "$threads")
			else
				_makecatcmd=("rapidgzip" "-qkcd" "-P" "$threads")
			fi
			;;
		*)
			if [[ $range ]]; then
				_makecatcmd=("fftool.sh" "-r" "$range" "-f")
			else
				_makecatcmd=("cat")
			fi
	esac

	[[ $v ]] || "${_makecatcmd[@]}" "$f"
	return 0
}

function helper::cat(){
	helper::makecatcmd "$@"
}

function helper::basename(){
	function _usage(){
		commander::print <<- 'EOF'
			description
			helper::basename returns a truncated basename of a given file by one extension or two in case of compressed data

			synopsis
			helper::basename <OPTIONS>

			options
			-f <file>     | mandatory. path to input file
			-o <variable> | mandatory. variable to store the truncated basename in
			-e <variable> | mandatory. variable to store the truncation i.e. extension in

			example
			helper::basename -f foo.txt.gz -o base -e ex
			echo $base
			echo $ex
		EOF
		return 1
	}

	local OPTIND arg f mandatory
	declare -n _basename _basenamex
	while getopts 'f:o:e:' arg; do
		case $arg in
			f) ((++mandatory)); f="$OPTARG";;
			o) ((++mandatory)); _basename=$OPTARG;;
			e) ((++mandatory)); _basenamex=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage

	if readlink -e "$f" | file -b --mime-type -f - | grep -qF -e 'gzip' -e 'bzip2'; then
		_basename="$(basename "$f" | rev | cut -d '.' -f 3- | rev)"
		_basenamex="$(basename "$f" | rev | cut -d '.' -f 1-2 | rev)"
	else
		_basename="$(basename "$f" | rev | cut -d '.' -f 2- | rev)"
		_basenamex="$(basename "$f" | rev | cut -d '.' -f 1 | rev)"
	fi

	return 0
}

function helper::ps2pdf(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			helper::ps2pdf converts a ps file into pdf while preserving dimensions

			-f <file> | mandatory. path to input ps file
			-o <file> | optional. path to output pdf file
		EOF
		return 1
	}

	local OPTIND arg mandatory f o
	while getopts 'f:o:' arg; do
		case $arg in
			f) ((++mandatory)); f="$OPTARG";;
			o) o="$OPTARG"; mkdir -p "$(dirname "$o")";;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 1 ]] && _usage
	[[ $o ]] || o="${f%.*}.pdf"

	ps2pdf -g$(grep -m 1 -F BoundingBox "$f" | sed -E 's/.+\s+([0-9]+)\s+([0-9]+)$/\10x\20/') "$f" "$o"

	return 0
}

function helper::join(){
	helper::multijoin "$@"
}

function helper::multijoin(){
	function _usage(){
		commander::print <<- 'EOF'
			description
			helper::multijoin, given multiple lexicographically sorted input files, returns a fully joined table, while executing pairwise full-join commands in parallel

			synopsis
			helper::multijoin <OPTIONS> [<files>]

			options
			-t <threads>   | optional. number of threads to use (see -e). default: 1
			-s <separator> | optional. character to separate columns. default: '\t'
			-n <columns>   | optional. number of leading columns to be used as id group. default: 1
			-r <range>     | optional. range of columns to be reported (i[,-][j]) of columns to be reported. default: 1-
			-h <header>    | optional. header string to be reported
			-e <empty>     | optional. string for empty/null fields. default: NA
			-o <outfile>   | optional. path to output table. default: stdout
			-b <seqidsfile>| optional. data is positional sorted bed/bedgraph format with sequence identifiers ordered according to first column of given file
			-d             | optional. all columns to be used as id group contain numerically, increasingly sorted digits

			example
			helper::multijoin <path> <path> [<path> ..]

			developer options
			-1 <cmds>     | array to append crafted commands to instead of execution
			-p <tmpdir>   | path to unique temporary directory

			example - pre-declare command arrays to stay in local scope. array number equals the levels/height of a complete binary tree i.e. ceil[log2(#leaves)]
			declare -a cmds
			for h in $(seq 0 $(echo <number_files> | awk '{h=log($1+1)/log(2); h=h>int(h)?int(h)+1:h; print h-1}')); do
				declare -a cmd_$h
				cmds[$h]=cmd_$h
			done
			helper::multijoin -1 cmds -p <tmpdir> <path> <path> [<path> ..]
			for c in ${cmds[@]}; do
				commander::runcmd -a $c
			done
		EOF
		return 1
	}

	local OPTIND arg mandatory outfile=/dev/stdout header empty=NA sep=$'\t' tmpdir group=1 range=1- bed=false execute=true threads=1 digits=false seqids
	declare -n _cmds_multijoin _files_multijoin
	while getopts '1:t:s:h:e:o:n:r:p:f:b:d' arg; do
		case $arg in
			1)	execute=false; _cmds_multijoin=$OPTARG;;
			t)	threads="$OPTARG";;
			s)	sep=$(echo -e "$OPTARG");;
			e)	empty="$OPTARG";;
			n)	group="$OPTARG";;
			r)	range="$OPTARG";;
			h)	header="$OPTARG";;
			o)	outfile="$OPTARG"; mkdir -p "$(dirname "$outfile")";;
			p)	tmpdir="$OPTARG";;
			f)	_files_multijoin="$OPTARG";;
			b)	seqids="$OPTARG"; bed=true;;
			d)	digits=true;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }

	if $execute; then
		unset _cmds_multijoin
		declare -a _cmds_multijoin
		tmpdir=$(mktemp -d -p "${TMPDIR:-/tmp}" cleanup.XXXXXXXXXX.multijoin)
	else
		[[ $tmpdir ]] || _usage
	fi

	shift $((OPTIND-1))
	if [[ $1 ]]; then
		declare -a tojoin=("$@") joined
	else
		[[ $_files_multijoin ]] || _usage
		declare -a tojoin=("${_files_multijoin[@]}") joined
	fi

	local i j f1 f2 o r=0
	if $bed; then
		while [[ ${#tojoin[@]} -gt 1 ]]; do
			joined=()
			j=0
			if [[ ${_cmds_multijoin[$r]} ]]; then
				declare -n __cmds_multijoin=${_cmds_multijoin[$r]}
			else
				declare -g -a multijoincmd$r
				declare -n __cmds_multijoin=multijoincmd$r
				__cmds_multijoin=()
				_cmds_multijoin[$r]=multijoincmd$r
			fi

			for i in $(seq 0 2 $((${#tojoin[@]}-1))); do
				if [[ ${tojoin[$((i+1))]} ]]; then
					f1="${tojoin[$i]}"
					f2="${tojoin[$((i+1))]}"
					if [[ ${#tojoin[@]} -eq 2 ]]; then
						o="$outfile"
						# padding positions by leading zeroes
						# and use padded rank instead of fasta seq ids.
						commander::makecmd -s ' ' -a __cmds_multijoin -c {COMMANDER[0]}<<-CMD {COMMANDER[1]}<<-CMD {COMMANDER[2]}<<-'CMD' {COMMANDER[3]}<<-CMD {COMMANDER[4]}<<-'CMD' {COMMANDER[5]}<<-CMD {COMMANDER[6]}<<-CMD {COMMANDER[7]}<<-CMD
							$([[ $header ]] && echo "echo -e '$header' > '$o';" || rm -f "$o" &> /dev/null || true)
						CMD
							join -o 0\$({ head -1 "$f1"; head -1 "$f2"; } | awk -F '\t' '{for(i=2;i<=NF-2;i++){printf ",%s",NR"."i}}') -1 1 -2 1 -a 1 -a 2 -e "$empty" -t \$'\t'
						CMD
							<(awk -F '\t' -v OFS='\t' '{f+=FNR==1?1:0} f==1{h[$1]=NR} f==2{$3=sprintf("%010d",h[$1])sprintf("%010d",$2)sprintf("%010d",$3)"#:::#"$1"#:::#"$2"#:::#"$3; print}'
						CMD
								"$seqids" "$f1" | cut -f 3-)
						CMD
							<(awk -F '\t' -v OFS='\t' '{f+=FNR==1?1:0} f==1{h[$1]=NR} f==2{$3=sprintf("%010d",h[$1])sprintf("%010d",$2)sprintf("%010d",$3)"#:::#"$1"#:::#"$2"#:::#"$3; print}'
						CMD
								"$seqids" "$f2" | cut -f 3-)
						CMD
							| sed 's/#:::#/\t/g' | cut -f 2- | cut -f $range >> "$o";
						CMD
							$([[ "$f1" == *"$tmpdir/tojoin"* ]] && echo "rm '$f1';"; [[ "$f2" == *"$tmpdir/tojoin"* ]] && echo "rm '$f2'")
						CMD
					else
						((++j))
						o="$tmpdir/tojoin.$r.$j"
						commander::makecmd -s ' ' -a __cmds_multijoin -c {COMMANDER[0]}<<-CMD {COMMANDER[1]}<<-'CMD' {COMMANDER[2]}<<-CMD {COMMANDER[3]}<<-'CMD' {COMMANDER[4]}<<-CMD {COMMANDER[5]}<<-CMD {COMMANDER[6]}<<-CMD
							join -o 0\$({ head -1 "$f1"; head -1 "$f2"; } | awk -F '\t' '{for(i=2;i<=NF-2;i++){printf ",%s",NR"."i}}') -1 1 -2 1 -a 1 -a 2 -e "$empty" -t \$'\t'
						CMD
							<(awk -F '\t' -v OFS='\t' '{f+=FNR==1?1:0} f==1{h[$1]=NR} f==2{$3=sprintf("%010d",h[$1])sprintf("%010d",$2)sprintf("%010d",$3)"#:::#"$1"#:::#"$2"#:::#"$3; print}'
						CMD
								"$seqids" "$f1" | cut -f 3-)
						CMD
							<(awk -F '\t' -v OFS='\t' '{f+=FNR==1?1:0} f==1{h[$1]=NR} f==2{$3=sprintf("%010d",h[$1])sprintf("%010d",$2)sprintf("%010d",$3)"#:::#"$1"#:::#"$2"#:::#"$3; print}'
						CMD
								"$seqids" "$f2" | cut -f 3-)
						CMD
							| sed 's/#:::#/\t/g' | cut -f 2- > "$o";
						CMD
							$([[ "$f1" == *"$tmpdir/tojoin"* ]] && echo "rm '$f1';"; [[ "$f2" == *"$tmpdir/tojoin"* ]] && echo "rm '$f2'";)
						CMD
						joined+=("$o")
					fi
				else
					joined+=("${tojoin[$i]}")
				fi
			done
			tojoin=("${joined[@]}")
			((++r))
		done
	elif $digits; then
		local zeroes=$(tail -q -n 1 "${tojoin[@]}" | cut -f 1-$group | tr "$sep" '\n' | awk '{x=length($0);m=m>x?m:x}END{print m}')
		while [[ ${#tojoin[@]} -gt 1 ]]; do
			joined=()
			j=0
			if [[ ${_cmds_multijoin[$r]} ]]; then
				declare -n __cmds_multijoin=${_cmds_multijoin[$r]}
			else
				declare -g -a multijoincmd$r
				declare -n __cmds_multijoin=multijoincmd$r
				__cmds_multijoin=()
				_cmds_multijoin[$r]=multijoincmd$r
			fi

			for i in $(seq 0 2 $((${#tojoin[@]}-1))); do
				if [[ ${tojoin[$((i+1))]} ]]; then
					f1="${tojoin[$i]}"
					f2="${tojoin[$((i+1))]}"
					if [[ ${#tojoin[@]} -eq 2 ]]; then
						o="$outfile"
						# padding positions by leading zeroes
						# and use padded rank instead of fasta seq ids.
						commander::makecmd -s ' ' -a __cmds_multijoin -c {COMMANDER[0]}<<-CMD {COMMANDER[1]}<<-CMD {COMMANDER[2]}<<-CMD {COMMANDER[3]}<<-'CMD' {COMMANDER[4]}<<-CMD {COMMANDER[5]}<<-'CMD' {COMMANDER[6]}<<-CMD {COMMANDER[7]}<<-CMD
							$([[ $header ]] && echo "echo -e '$header' > '$o';" || rm -f "$o" &> /dev/null || true)
						CMD
							join -o 0\$({ head -1 "$f1"; head -1 "$f2"; } | awk -F "$sep" -v r=$group '{for(i=2;i<=NF-r+1;i++){printf ",%s",NR"."i}}') -1 1 -2 1 -a 1 -a 2 -e "$empty" -t \$'\t'
						CMD
							<(cat "$f1" | awk -F '\t' -v OFS='\t' -v r=$group -v n=$zeroes
						CMD
								'{x=$r; $r=sprintf("%010d",$r); for(i=1;i<r;i++){x=$i"#:::#"x; $r=sprintf("%0"n"d",$i)$r; $i=""} $r=$r"#:::#"x; print}' | sed 's/^\s*//')
						CMD
							<(cat "$f2" | awk -F '\t' -v OFS='\t' -v r=$group -v n=$zeroes
						CMD
								'{x=$r; $r=sprintf("%010d",$r); for(i=1;i<r;i++){x=$i"#:::#"x; $r=sprintf("%0"n"d",$i)$r; $i=""} $r=$r"#:::#"x; print}' | sed 's/^\s*//')
						CMD
							| sed 's/#:::#/\t/g' | cut -f 2- | cut -f $range >> "$o";
						CMD
							$([[ "$f1" == *"$tmpdir/tojoin"* ]] && echo "rm '$f1';"; [[ "$f2" == *"$tmpdir/tojoin"* ]] && echo "rm '$f2'")
						CMD
					else
						((++j))
						o="$tmpdir/tojoin.$r.$j"
						commander::makecmd -s ' ' -a __cmds_multijoin -c {COMMANDER[0]}<<-CMD {COMMANDER[1]}<<-CMD {COMMANDER[2]}<<-'CMD' {COMMANDER[3]}<<-CMD {COMMANDER[4]}<<-'CMD' {COMMANDER[5]}<<-CMD {COMMANDER[6]}<<-CMD
							join -o 0\$({ head -1 "$f1"; head -1 "$f2"; } | awk -F "$sep" -v r=$group '{for(i=2;i<=NF-r+1;i++){printf ",%s",NR"."i}}') -1 1 -2 1 -a 1 -a 2 -e "$empty" -t \$'\t'
						CMD
							<(cat "$f1" | awk -F '\t' -v OFS='\t' -v r=$group -v n=$zeroes
						CMD
								'{x=$r; $r=sprintf("%010d",$r); for(i=1;i<r;i++){x=$i"#:::#"x; $r=sprintf("%0"n"d",$i)$r; $i=""} $r=$r"#:::#"x; print}' | sed 's/^\s*//')
						CMD
							<(cat "$f2" | awk -F '\t' -v OFS='\t' -v r=$group -v n=$zeroes
						CMD
								'{x=$r; $r=sprintf("%010d",$r); for(i=1;i<r;i++){x=$i"#:::#"x; $r=sprintf("%0"n"d",$i)$r; $i=""} $r=$r"#:::#"x; print}' | sed 's/^\s*//')
						CMD
							| sed 's/#:::#/\t/g' | cut -f 2- > "$o";
						CMD
							$([[ "$f1" == *"$tmpdir/tojoin"* ]] && echo "rm '$f1';"; [[ "$f2" == *"$tmpdir/tojoin"* ]] && echo "rm '$f2'";)
						CMD
						joined+=("$o")
					fi
				else
					joined+=("${tojoin[$i]}")
				fi
			done
			tojoin=("${joined[@]}")
			((++r))
		done
	elif [[ $group -eq 1 ]]; then
		while [[ ${#tojoin[@]} -gt 1 ]]; do
			joined=()
			j=0
			if [[ ${_cmds_multijoin[$r]} ]]; then
				declare -n __cmds_multijoin=${_cmds_multijoin[$r]}
			else
				declare -g -a multijoincmd$r
				declare -n __cmds_multijoin=multijoincmd$r
				__cmds_multijoin=()
				_cmds_multijoin[$r]=multijoincmd$r
			fi

			for i in $(seq 0 2 $((${#tojoin[@]}-1))); do
				if [[ ${tojoin[$((i+1))]} ]]; then
					f1="${tojoin[$i]}"
					f2="${tojoin[$((i+1))]}"
					if [[ ${#tojoin[@]} -eq 2 ]]; then
						o="$outfile"
						# needs true in case of /dev/stdout
						commander::makecmd -m -a __cmds_multijoin -c <<-CMD
							$([[ $header ]] && echo "echo -e '$header' > '$o'" || rm -f "$o" &> /dev/null || true)
							join -o 0\$({ head -1 "$f1"; head -1 "$f2"; } | awk -F "$sep" '{for(i=2;i<=NF;i++){printf ",%s",NR"."i}}') -1 1 -2 1 -a 1 -a 2 -e "$empty" -t "$sep" "$f1" "$f2" | cut -d "$sep" -f $range >> "$o"
							$([[ "$f1" == *"$tmpdir/tojoin"* ]] && echo "rm '$f1';"; [[ "$f2" == *"$tmpdir/tojoin"* ]] && echo "rm '$f2'")
						CMD
					else
						((++j))
						o="$tmpdir/tojoin.$r.$j"
						commander::makecmd -m -a __cmds_multijoin -c <<-CMD
							join -o 0\$({ head -1 "$f1"; head -1 "$f2"; } | awk -F "$sep" '{for(i=2;i<=NF;i++){printf ",%s",NR"."i}}') -1 1 -2 1 -a 1 -a 2 -e "$empty" -t "$sep" "$f1" "$f2" > "$o"
							$([[ "$f1" == *"$tmpdir/tojoin"* ]] && echo "rm '$f1';"; [[ "$f2" == *"$tmpdir/tojoin"* ]] && echo "rm '$f2'")
						CMD
						joined+=("$o")
					fi
				else
					joined+=("${tojoin[$i]}")
				fi
			done
			tojoin=("${joined[@]}")
			((++r))
		done
	else
		while [[ ${#tojoin[@]} -gt 1 ]]; do
			joined=()
			j=0
			if [[ ${_cmds_multijoin[$r]} ]]; then
				declare -n __cmds_multijoin=${_cmds_multijoin[$r]}
			else
				declare -g -a multijoincmd$r
				declare -n __cmds_multijoin=multijoincmd$r
				__cmds_multijoin=()
				_cmds_multijoin[$r]=multijoincmd$r
			fi

			for i in $(seq 0 2 $((${#tojoin[@]}-1))); do
				if [[ ${tojoin[$((i+1))]} ]]; then
					f1="${tojoin[$i]}"
					f2="${tojoin[$((i+1))]}"
					if [[ ${#tojoin[@]} -eq 2 ]]; then
						o="$outfile"
						# needs true in case of /dev/stdout
						commander::makecmd -m -a __cmds_multijoin -c <<-CMD
							$([[ $header ]] && echo "echo -e '$header' > '$o'" || rm -f "$o" &> /dev/null || true)
							join -o 0\$({ head -1 "$f1"; head -1 "$f2"; } | awk -F "$sep" -v r=$group '{for(i=2;i<=NF-r+1;i++){printf ",%s",NR"."i}}') -1 1 -2 1 -a 1 -a 2 -e "$empty" -t "$sep" <(paste -d "$sep" <(cut -d "$sep" -f 1-$group "$f1" | sed "s/$sep/#:::#/g") <(cut -d "$sep" -f $((group+1))- "$f1")) <(paste -d "$sep" <(cut -d "$sep" -f 1-$group "$f2" | sed "s/$sep/#:::#/g") <(cut -d "$sep" -f $((group+1))- "$f2")) | sed "s/#:::#/$sep/g" | cut -d "$sep" -f $range >> "$o"
							$([[ "$f1" == *"$tmpdir/tojoin"* ]] && echo "rm '$f1';"; [[ "$f2" == *"$tmpdir/tojoin"* ]] && echo "rm '$f2'")
						CMD
					else
						((++j))
						o="$tmpdir/tojoin.$r.$j"
						commander::makecmd -m -a __cmds_multijoin -c <<-CMD
							join -o 0\$({ head -1 "$f1"; head -1 "$f2"; } | awk -F "$sep" -v r=$group '{for(i=2;i<=NF-r+1;i++){printf ",%s",NR"."i}}') -1 1 -2 1 -a 1 -a 2 -e "$empty" -t "$sep" <(paste -d "$sep" <(cut -d "$sep" -f 1-$group "$f1" | sed "s/$sep/#:::#/g") <(cut -d "$sep" -f $((group+1))- "$f1")) <(paste -d "$sep" <(cut -d "$sep" -f 1-$group "$f2" | sed "s/$sep/#:::#/g") <(cut -d "$sep" -f $((group+1))- "$f2")) | sed "s/#:::#/$sep/g" > "$o"
							$([[ "$f1" == *"$tmpdir/tojoin"* ]] && echo "rm '$f1';"; [[ "$f2" == *"$tmpdir/tojoin"* ]] && echo "rm '$f2'")
						CMD
						joined+=("$o")
					fi
				else
					joined+=("${tojoin[$i]}")
				fi
			done
			tojoin=("${joined[@]}")
			((++r))
		done
	fi

	if $execute; then
		local c
		for c in "${_cmds_multijoin[@]}"; do
			commander::runcmd -i $threads -a $c
		done
	fi

	return 0
}

function helper::ishash(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::ishash returns without error if given variable or refers to an associative array i.e. hash

			-v <var> | variable
		EOF
		return 1
	}

	local OPTIND arg
	while getopts 'v:' arg; do
		case $arg in
			v)	declare -n __="$OPTARG"
				[[ ${!__@a} == *A* ]] || return 1
			;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }

	return 0
}

function helper::isarray(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::isarray returns without error if given variable is or refers to an array

			-v <var> | variable
		EOF
		return 1
	}

	local OPTIND arg
	while getopts 'v:' arg; do
		case $arg in
			v)	declare -n __="$OPTARG"
				[[ ${!__@a} == *a* ]] || return 1
			;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }

	return 0
}

function helper::addmemberfunctions(){
	function _usage(){
		commander::print <<- 'EOF'
			helper::addmemberfunctions adds bashbone builtin OOP like syntax to variables

			-v <var> | variable
		EOF
		return 1
	}

	local OPTIND arg mandatory var
	declare -a vars
	while getopts 'v:' arg; do
		case $arg in
			v) ((++mandatory)); vars+=("$(printf '%q' "$OPTARG")");;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	shopt -s expand_aliases
	for var in "${vars[@]}"; do
		eval "alias $var.get='helper::_get $var'"
		eval "alias $var.push='helper::_push $var'"
		eval "alias $var.pop='helper::_pop $var'"
		eval "alias $var.slice='helper::_slice $var'"
		eval "alias $var.join='helper::_join $var'"
		eval "alias $var.print='helper::_print $var'"
		eval "alias $var.println='helper::_println $var'"
		eval "alias $var.shift='helper::_shift $var'"
		eval "alias $var.length='helper::_length $var'"
		eval "alias $var.lastidx='helper::_lastidx $var'"
		eval "alias $var.idxs='helper::_idxs $var'"
		eval "alias $var.uc='helper::_uc $var'"
		eval "alias $var.ucfist='helper::_ucfirst $var'"
		eval "alias $var.lc='helper::_lc $var'"
		eval "alias $var.lcfist='helper::_lcfirst $var'"
		eval "alias $var.sum='helper::_sum $var'"
		eval "alias $var.trimprefix='helper::_prefix $var'"
		eval "alias $var.trimprefixfirst='helper::_prefixfirst $var'"
		eval "alias $var.trimsuffix='helper::_trimsuffix $var'"
		eval "alias $var.trimsuffixfirst='helper::_suffixfirst $var'"
		eval "alias $var.substring='helper::_substring $var'"
		eval "alias $var.replace='helper::_replace $var'"
		eval "alias $var.replaceprefix='helper::_replaceprefix $var'"
		eval "alias $var.replacesuffix='helper::_replacesuffix $var'"
		eval "alias $var.uniq='helper::_uniq $var'"
		eval "alias $var.sort='helper::_sort $var'"
		eval "alias $var.basename='helper::_basename $var'"
		eval "alias $var.dirname='helper::_dirname $var'"
	done

	return 0
}

function helper::_get(){
	declare -n __="$1"
	[[ $2 ]] && {
		local j=1
		[[ $3 ]] && {
			[[ $3 -lt 0 ]] && j=$((${#__[@]}-$2+$3)) || {
				[[ $3 -eq 0 ]] && j=${#__[@]} || j=$(($3-$2+1))
			}
		}
		echo "${__[@]:$2:$j}"
	} || echo ${__[*]}
}

function helper::_push(){
	declare -n __="$1"
	shift
	__+=("$@")
}

function helper::_pop(){
	declare -n __="$1"
	__=("${__[@]:0:$((${#__[@]}-1))}")
}

function helper::_slice(){
	declare -n __="$1"
	local j=$(($3-$2+1))
	__=("${__[@]:$2:$3}")
}

function helper::_join(){
	declare -n __="$1" ___
	shift
	for ___ in "$@"; do
		__+=("${___[@]}")
	done
}

function helper::_shift(){
	declare -n __="$1"
	__=("${__[@]:1}")
}

function helper::_idxs(){
	declare -n __="$1"
	echo "${!__[@]}"
}

function helper::_lastidx(){
	declare -n __="$1"
	echo $((${#__[@]}-1))
}

function helper::_length(){
	declare -n __="$1"
	echo ${#__[@]}
}

function helper::_print(){
	declare -n __="$1"
	echo ${__[*]}
}

function helper::_println(){
	declare -n __="$1"
	printf '%s\n' "${__[@]}"
}

function helper::_uc(){
	declare -n __="$1"
	__=("${__[@]^^${2:-*}}")
}

function helper::_ucfirst(){
	declare -n __="$1"
	__=("${__[@]^${2:-*}}")
}

function helper::_lc(){
	declare -n __="$1"
	__=("${__[@],,${2:-*}}")
}

function helper::_lcfirst(){
	declare -n __="$1"
	__=("${__[@],${2:-*}}")
}

function helper::_sum(){
	declare -n __="$1"
	__=$(("${__[@]/%/+}"0))
}

function helper::_trimsuffixfirst(){
	declare -n __="$1"
	__=("${__[@]%"$2"*}")
}

function helper::_trimsuffix(){
	declare -n __="$1"
	__=("${__[@]%%"$2"*}")
}

function helper::_trimprefixfirst(){
	declare -n __="$1"
	__=("${__[@]#*"$2"}")
}

function helper::_trimprefix(){
	declare -n __="$1"
	__=("${__[@]##*"$2"}")
}

function helper::_substring(){
	declare -n __="$1"
	local i
	if [[ $3 ]]; then
		for i in "${!__[@]}"; do
			__[$i]="${__[$i]:$2:$3}"
		done
	else
		for i in "${!__[@]}"; do
			__[$i]="${__[$i]:$2}"
		done
	fi
}

function helper::_replace(){
	declare -n __="$1"
	__=("${__[@]/${2:-*}/"$3"}")
}

function helper::_replaceprefix(){
	declare -n __="$1"
	__=("${__[@]/#${2:-*}/"$3"}")
}

function helper::_replacesuffix(){
	declare -n __="$1"
	__=("${__[@]/%${2:-*}/"$3"}")
}

function helper::_uniq(){
	declare -n __="$1"
	declare -A ___
	local e
	for e in "${__[@]}"; do
		___["$e"]=1
	done
	__=("${!___[@]}")
}

function helper::_sort(){
	declare -n __="$1"
	mapfile -t __ < <(printf '%s\n' "${__[@]}" | sort -V)
}

function helper::_basename(){
	declare -n __="$1"
	local i
	for i in "${!__[@]}"; do
		${__[$i]}="$(basename "${__[$i]}" "$2")"
	done
}

function helper::_dirname(){
	declare -n __="$1"
	local i
	for i in "${!__[@]}"; do
		${__[$i]}="$(dirname "${__[$i]}")"
	done
}

function helper::_test(){
	declare -a arr
	helper::addmemberfunctions -v arr
	arr.push "f.o.o f.o.o"
	arr.push bar
	arr.push zar
	arr.print
	arr.get 0 -1
	arr.get 1 0
	arr.sort
	arr.print
	arr.uc [fa]
	arr.print
	arr.replace A X
	arr.print
	arr.lc
	arr.print
	arr.slice 1 2
	arr.print
	arr.substring 1 2
	arr.print
	arr.uniq
	arr.print
}
