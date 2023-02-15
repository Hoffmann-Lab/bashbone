#!/usr/bin/env bash
# (c) Konstantin Riege

function helper::pgzip(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-f <infile>  | path to. else stdin
			-o <outfile> | path to
			-t <threads> | number of
			-b           | block zip via bgzip
		EOF
		return 1
	}

	local OPTIND arg threads=1 f=/dev/stdin o tool="pigz -p"
	while getopts 'f:t:o:b' arg; do
		case $arg in
			f) f="$OPTARG";;
			o) o="$OPTARG"; mkdir -p "$(dirname "$o")";;
			t) threads=$OPTARG;;
			b) tool="bgzip -@";;
			*) _usage;;
		esac
	done
	[[ ! $f && ! $o ]] && _usage
	[[ ! $o ]] && o="$f.gz"

	$tool $threads -k -c "$f" | tee -i "$o" | gztool -v 0 -f -i -x -C -I "${o%.*}.gzi"

	return 0
}

function helper::sort(){
	declare -a tdirs
	function _cleanup::helper::sort(){
		rm -rf "${tdirs[@]}"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-f <infile>    | path to. else stdin
			-o <outfile>   | path to. else stdout
			-t <threads>   | number of
			-M <maxmemory> | amount of
		EOF
		return 1
	}

	local OPTIND threads=1 f=/dev/stdin o=/dev/stdout maxmemory tmpdir="${TMPDIR:-/tmp}"
	declare -a args=();
	while [[ $# -gt 0 ]]; do
		case "$1" in
			-m)	maxmemory=$2; shift 2;;
			-f)	f=$2; shift 2;;
			-o)	o=$2; shift 2; mkdir -p "$(dirname "$o")";;
			-t)	threads=$2; shift 2;;
			*)	args+=("$1"); shift;;
		esac
	done

	local instances maxmemory
	read -r instances maxmemory < <(configure::memory_by_instances -i 1 -M "$maxmemory")
	tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.sort)")
	sort --parallel="$threads" -S "$maxmemory" -T "${tdirs[-1]}" "${args[@]}" "$f" > "$o"

	return 0
}

function helper::vcfsort(){
	declare -a tdirs
	function _cleanup::helper::vcfsort(){
		rm -rf "${tdirs[@]}"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-f <infile>    | path to. else stdin
			-o <outfile>   | path to. else stdout
			-t <threads>   | number of
			-M <maxmemory> | amount of
			-z             | compress output
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

	local instances maxmemory
	read -r instances maxmemory < <(configure::memory_by_instances -i 1 -M "$maxmemory")
	tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.vcfsort)")
	if $zip; then
		bcftools view "$f" | awk -F'\t' -v t=$threads -v m=$maxmemory -v p="${tdirs[-1]}" -v OFS='\t' 'match($0,/.*#contig=<ID=(\S+),length.*/,a){i++; c[a[1]]=i} /^#/{print; next} {$1=c[$1]; print | "sort -k1,1n -k2,2n -k4,4 -k5,5 --parallel="t" -S "m"M -T \""p"\""}' | awk -F'\t' -v OFS='\t' 'match($0,/.*#contig=<ID=(\S+),length.*/,a){i++; c[i]=a[1]} /^#/{print; next} {$1=c[$1]; print}' | bgzip -k -c -@ $threads /dev/stdin > "$o"
		tabix -f -p vcf "$o"
	else
		bcftools view "$f" | awk -F'\t' -v t=$threads -v m=$maxmemory -v p="${tdirs[-1]}" -v OFS='\t' 'match($0,/.*#contig=<ID=(\S+),length.*/,a){i++; c[a[1]]=i} /^#/{print; next} {$1=c[$1]; print | "sort -k1,1n -k2,2n -k4,4 -k5,5 --parallel="t" -S "m"M -T \""p"\""}' | awk -F'\t' -v OFS='\t' 'match($0,/.*#contig=<ID=(\S+),length.*/,a){i++; c[i]=a[1]} /^#/{print; next} {$1=c[$1]; print}' > "$o"
	fi

	return 0
}

function helper::makecatcmd(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-c <var>  | cmd
			-f <file> | to ascii
			example:
			${FUNCNAME[1]} -c cmd -f [txt|bz2|gz]
		EOF
		return 1
	}

	local OPTIND arg mandatory f
	declare -n _makecatcmd
	while getopts 'f:c:' arg; do
		case $arg in
			c) ((++mandatory)); _makecatcmd=$OPTARG;;
			f) ((++mandatory)); f="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	_makecatcmd=$({ readlink -e "$f" | file -f - | grep -Eo '(gzip|bzip)' || echo cat; } | sed '/gzip/{s/gzip/pigz -p 1/; s/$/ -cd/}')

	return 0
}

function helper::basename(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-f <file> | path to
			-o <var>  | basename
			-e <var>  | extension

			example:
			${FUNCNAME[1]} -f foo.txt.gz -o base -e ex
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
	[[ $mandatory -lt 3 ]] && _usage

	if readlink -e "$f" | file -f - | grep -qE '(gzip|bzip)'; then
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
			${FUNCNAME[1]} usage:
			-f <file> | to ps
		EOF
		return 1
	}

	local OPTIND arg mandatory f
	while getopts 'f:' arg; do
		case $arg in
			f) ((++mandatory)); f="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	ps2pdf -g$(grep -m 1 -F BoundingBox "$f" | sed -E 's/.+\s+([0-9]+)\s+([0-9]+)$/\10x\20/') "$f" "${f%.*}.pdf"

	return 0
}

function helper::join(){
	helper::multijoin "$@"
}

function helper::multijoin(){
	local tmp joined
	function _cleanup::helper::multijoin(){
		rm -f "$tmp" "$joined"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-s <separator> | i/o character - default: tab
			-h <header>    | string of
			-e <empty>     | cell character - default: 0
			-o <outfile>   | path to
			-f <files>     | ALWAYS LAST OPTION
			                 tab-separated, to join by unique id in first column
			example:
			${FUNCNAME[1]} -f <path> <path> [<path> ..]
		EOF
		return 1
	}

	local OPTIND arg mandatory outfile=/dev/stdout empty=0 sep=$'\t' tmpdir="${TMPDIR:-/tmp}"
	declare -n _makepdfcmd
	while getopts 's:h:e:o:f:' arg; do
		case $arg in
			s)	sep=$(echo -e "$OPTARG");;
			e)	empty="$OPTARG";;
			h)	header="$OPTARG";;
			o)	outfile="$OPTARG"; mkdir -p "$(dirname "$outfile")";;
			f)	((++mandatory)); shift $((OPTIND-2)); break;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	local format i
	tmp="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.join)"
	joined="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.joined)"
	join -t "$sep" -1 1 -2 1 -a 1 -a 2 -e "$empty" -o '0,1.2,2.2' "$1" "$2" > "$joined"
	for i in $(seq 3 $#); do
		format="0"
		for j in $(seq 2 $i); do format+=",1.$j"; done
		format+=",2.2"
		join -t "$sep" -1 1 -2 1 -a 1 -a 2 -e "$empty" -o "$format" "$joined" "${!i}" > "$tmp"
		mv "$tmp" "$joined"
	done

	if [[ $header ]]; then
		cat <(echo -e "$header") "$joined" > "$outfile"
	else
		cat "$joined" > "$outfile"
	fi

	return 0
}

function helper::ishash(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v <var> | variable
		EOF
		return 1
	}

	local OPTIND arg
	while getopts 'v:' arg; do
		case $arg in
			v)	{	declare -p "$OPTARG" &> /dev/null
					declare -n __="$OPTARG"
					[[ "$(declare -p ${!__})" =~ ^declare\ \-[Ab-z]+ ]]
				} || return 1
			;;
			*)	_usage;;
		esac
	done

	return 0
}

function helper::isarray(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v <var> | variable
		EOF
		return 1
	}

	local OPTIND arg mandatory var
	declare -a vars
	while getopts 'v:' arg; do
		case $arg in
			v)	{	declare -p "$OPTARG" &> /dev/null
					declare -n __="$OPTARG"
					[[ "$(declare -p ${!__})" =~ ^declare\ \-[a-zB-Z]+ ]]
				} || return 1
			;;
			*)	_usage;;
		esac
	done

	return 0
}

function helper::addmemberfunctions(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
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
	__=("${__[@]/${2:-*}%/"$3"}")
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
