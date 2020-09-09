#!/usr/bin/env bash
# (c) Konstantin Riege

helper::makezipcmd(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-a <cmds>    | array of
			-t <threads> | number of
			-c <file>    | compress if not compressed
			-z <var>     | of path to file
			example:
			$funcname -a cmds -c f1.txt -c f2.txt -z o1 -z o2
		EOF
		return 0
	}

	local OPTIND arg mandatory threads=1
	declare -a check_makezipcmd tozip_makezipcmd
	declare -n _cmds_makezipcmd
	while getopts 'a:t:c:z:' arg; do
		case $arg in
			a) ((mandatory++)); _cmds_makezipcmd=$OPTARG;;
			t) threads=$OPTARG;;
			c) ((mandatory++)); check_makezipcmd+=("$OPTARG");;
			z) ((mandatory++)); tozip_makezipcmd+=("$OPTARG");;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 3 || ${#check_makezipcmd[@]} -ne ${#tozip_makezipcmd[@]} ]] && _usage && return 1

	local i
	for i in "${!check_makezipcmd[@]}"; do
		readlink -e "${check_makezipcmd[$i]}" | file -f - | grep -qE '(gzip|bzip)' || {
			declare -n _f_makezipcmd=${tozip_makezipcmd[$i]}
			commander::makecmd -a _cmds_makezipcmd -s '|' -c {COMMANDER[0]}<<- CMD
				pigz -p $threads -k -c "$_f_makezipcmd" > "$_f_makezipcmd.gz"
			CMD
			_f_makezipcmd="$_f_makezipcmd.gz"
		}
	done

	return 0
}

helper::makecatcmd(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-c <var>  | cmd
			-f <file> | to ascii
			example:
			$funcname -f [txt|bz2|gz]
		EOF
		return 0
	}

	local OPTIND arg mandatory f
	declare -n _makecatcmd
	while getopts 'f:c:' arg; do
		case $arg in
			c) ((mandatory++)); _makecatcmd=$OPTARG;;
			f) ((mandatory++)); f="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage && return 1

	_makecatcmd=$(readlink -e "$f" | file -f - | grep -Eo '(gzip|bzip)' && echo -cd || echo cat)

	return 0
}

helper::basename(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-f <file> | path to
			-o <var>  | basename
			-e <var>  | extension

			example:
			$funcname -f foo.txt.gz -o base -e ex
		EOF
		return 0
	}

	local OPTIND arg f mandatory
	declare -n _basename _basenamex
	while getopts 'f:o:e:' arg; do
		case $arg in
			f) ((mandatory++)); f="$OPTARG";;
			o) ((mandatory++)); _basename=$OPTARG;;
			e) ((mandatory++)); _basenamex=$OPTARG;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage && return 1

	readlink -e "$f" | file -f - | grep -qE '(gzip|bzip)' && {
		_basename=$(basename $f | rev | cut -d '.' -f 3- | rev)
		_basenamex=$(basename $f | rev | cut -d '.' -f 1-2 | rev)
	} || {
		_basename=$(basename $f | rev | cut -d '.' -f 2- | rev)
		_basenamex=$(basename $f | rev | cut -d '.' -f 1 | rev)
	}

	return 0
}

helper::ishash(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v <var> | variable
		EOF
		return 0
	}
	local OPTIND arg mandatory var
	declare -a vars
	while getopts 'v:' arg; do
		case $arg in
			v)	((mandatory++)); vars+=("$(printf '%q' "$OPTARG")")
				{	declare -p "$OPTARG" &> /dev/null && \
					declare -n __="$OPTARG"
					[[ "$(declare -p ${!__})" =~ ^declare\ \-[Ab-z]+ ]]
				} || return 1
			;;
			*)	_usage; return 1;;
		esac
	done

	return 0
}

helper::isarray(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v <var> | variable
		EOF
		return 0
	}
	local OPTIND arg mandatory var
	declare -a vars
	while getopts 'v:' arg; do
		case $arg in
			v)	((mandatory++)); vars+=("$(printf '%q' "$OPTARG")")
				{	declare -p "$OPTARG" &> /dev/null && \
					declare -n __="$OPTARG"
					[[ "$(declare -p ${!__})" =~ ^declare\ \-[a-zB-Z]+ ]]
				} || return 1
			;;
			*)	_usage; return 1;;
		esac
	done

	return 0
}

helper::addmemberfunctions(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v <var> | variable
		EOF
		return 0
	}
	local OPTIND arg mandatory var
	declare -a vars
	while getopts 'v:' arg; do
		case $arg in
			v) ((mandatory++)); vars+=("$(printf '%q' "$OPTARG")");;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage && return 1

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
		eval "alias $var.trimsuffix='helper::_suffix $var'"
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

helper::_get(){
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

helper::_push(){
	declare -n __="$1"
	shift
	__+=("$@")
}

helper::_pop(){
	declare -n __="$1"
	__=("${__[@]:0:$((${#__[@]}-1))}")
}

helper::_slice(){
	declare -n __="$1"
	local j=$(($3-$2+1))
	__=("${__[@]:$2:$3}")
}

helper::_join(){
	declare -n __="$1" ___
	shift
	for ___ in "$@"; do
		__+=("${___[@]}")
	done
}

helper::_shift(){
	declare -n __="$1"
	__=("${__[@]:1}")
}

helper::_idxs(){
	declare -n __="$1"
	echo "${!__[@]}"
}

helper::_lastidx(){
	declare -n __="$1"
	echo $((${#__[@]}-1))
}

helper::_length(){
	declare -n __="$1"
	echo ${#__[@]}
}

helper::_print(){
	declare -n __="$1"
	echo ${__[*]}
}

helper::_println(){
	declare -n __="$1"
	printf '%s\n' "${__[@]}"
}

helper::_uc(){
	declare -n __="$1"
	__=("${__[@]^^${2:-*}}")
}

helper::_ucfirst(){
	declare -n __="$1"
	__=("${__[@]^${2:-*}}")
}

helper::_lc(){
	declare -n __="$1"
	__=("${__[@],,${2:-*}}")
}

helper::_lcfirst(){
	declare -n __="$1"
	__=("${__[@],${2:-*}}")
}

helper::_sum(){
	declare -n __="$1"
	__=$(("${__[@]/%/+}"0))
}

helper::_trimsuffixfirst(){
	declare -n __="$1"
	__=("${__[@]%"$2"*}")
}

helper::_trimsuffix(){
	declare -n __="$1"
	__=("${__[@]%%"$2"*}")
}

helper::_trimprefixfirst(){
	declare -n __="$1"
	__=("${__[@]#*"$2"}")
}

helper::_trimprefix(){
	declare -n __="$1"
	__=("${__[@]##*"$2"}")
}

helper::_substring(){
	declare -n __="$1"
	local i
	for i in "${!__[@]}"; do
		__[$i]="${__[$i]:$2:$3}"
	done
}

helper::_replace(){
	declare -n __="$1"
	__=("${__[@]/${2:-*}/"$3"}")
}

helper::_replaceprefix(){
	declare -n __="$1"
	__=("${__[@]/#${2:-*}/"$3"}")
}

helper::_replacesuffix(){
	declare -n __="$1"
	__=("${__[@]/${2:-*}%/"$3"}")
}

helper::_uniq(){
	declare -n __="$1"
	declare -A ___
	local e
	for e in "${__[@]}"; do
		___["$e"]=1
	done
	__=("${!___[@]}")
}

helper::_sort(){
	declare -n __="$1"
	mapfile -t __ < <(printf '%s\n' "${__[@]}" | sort -V)
}

helper::_basename(){
	declare -n __="$1"
	local i
	for i in "${!__[@]}"; do
		${__[$i]}="$(basename "${__[$i]}" "$2")"
	done
}

helper::_dirname(){
	declare -n __="$1"
	local i
	for i in "${!__[@]}"; do
		${__[$i]}="$(dirname "${__[$i]}")"
	done
}

helper::_test(){
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
