#!/usr/bin/env bash
# (c) Konstantin Riege

helper::makezipcmd(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage:
			-a <cmds>    | array of
			-t <threads> | number of
			-c <file>    | compress if not compressed
			-z <var>     | to compress file
			example: 
			$funcname -a cmds -c f1.txt -c f2.txt -z o1 -z o2
		EOF
		return 0
	}

	local OPTIND arg mandatory threads=1
	declare -a check_makezipcmd
	declare -a tozip_makezipcmd # be very careful with references name space
	declare -n _cmds_makezipcmd # be very careful with circular name reference
	while getopts 'a:t:c:z:' arg; do
		case $arg in
			a)	mandatory=1; _cmds_makezipcmd=$OPTARG;;
			t) 	threads=$OPTARG;;
			c)	check_makezipcmd+=("$OPTARG");;
			z)	tozip_makezipcmd+=("$OPTARG");;
			*)	_usage;	return 1;;
		esac
	done
	[[ ! $mandatory || ${#check_makezipcmd[@]} -ne ${#tozip_makezipcmd[@]} ]] && { _usage; return 1; }

    local i
	for i in "${!check_makezipcmd[@]}"; do
		readlink -e "${check_makezipcmd[$i]}" | file -f - | grep -qF compressed || {
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
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage:
			-c <var>  | cmd
			-f <file> | to ascii
			example: 
			$funcname -f [txt|bz2|gz]
		EOF
		return 0
	}

	local OPTIND arg f
	declare -n _ref_makecatcmd
	while getopts 'f:c:' arg; do
		case $arg in
			c)	_ref_makecatcmd=$OPTARG;;
			f)	f="$OPTARG";;
			*)	_usage; return 1;;
		esac
	done

	_ref_makecatcmd=$(readlink -e "$f" | file -f - | grep -Eo '(gzip|bzip)' && echo -cd || echo cat)

	return 0
}

helper::basename(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage:
			-f <file> | path to
			-o <var>  | basename
			-e <var>  | extension

			example: 
			$funcname -f foo.txt.gz -o base -e ex
		EOF
		return 0
	}

	local OPTIND arg f
	declare -n _ref_basename _ref_basenamex
	while getopts 'f:o:e:' arg; do
		case $arg in
			f)	f="$OPTARG";;
			o)	_ref_basename=$OPTARG;;
			e)	_ref_basenamex=$OPTARG;;
			*)	_usage;	return 1;;
		esac
	done

	readlink -e "$f" | file -f - | grep -qE '(gzip|bzip)' && {
		_ref_basename=$(basename $f | rev | cut -d '.' -f 3- | rev)
		_ref_basenamex=$(basename $f | rev | cut -d '.' -f 1-2 | rev)
	} || {
		_ref_basename=$(basename $f | rev | cut -d '.' -f 2- | rev)
		_ref_basenamex=$(basename $f | rev | cut -d '.' -f 1 | rev)
	}

	return 0
}

helper::loadaliases(){
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v <var> | variable
		EOF
		return 0
	}
	local OPTIND arg mandatory var
	while getopts 'v:' arg; do
		case $arg in
			v)	mandatory=1; var=$OPTARG;;
			*)	_usage;	return 1;;
		esac
	done
	[[ ! $mandatory ]] && _usage && return 1

	__=$(printf '%q' ${1})
	eval "alias $var.push='helper::_push $var'"
	eval "alias $var.pop='helper::_pop $var'"
	eval "alias $var.slice='helper::_slice $var'"
	eval "alias $var.print='helper::_print $var'"
	eval "alias $var.println='helper::_println $var'"
	eval "alias $var.shift='helper::_shift $var'"
	eval "alias $var.uc='helper::_uc $var'"
	eval "alias $var.ucfist='helper::_ucfirst $var'"
	eval "alias $var.lc='helper::_lc $var'"
	eval "alias $var.lcfist='helper::_lcfirst $var'"
	eval "alias $var.sum='helper::_sum $var'"
	eval "alias $var.prefix='helper::_prefix $var'"
	eval "alias $var.prefixfirst='helper::_prefixfirst $var'"
	eval "alias $var.suffix='helper::_suffix $var'"
	eval "alias $var.suffixfirst='helper::_suffixfirst $var'"
	eval "alias $var.substring='helper::_substring $var'"
	eval "alias $var.replace='helper::_replace $var'"
	eval "alias $var.replaceprefix='helper::_replaceprefix $var'"
	eval "alias $var.replacesuffix='helper::_replacesuffix $var'"
	eval "alias $var.uniq='helper::_uniq $var'"
	eval "alias $var.sort='helper::_sort $var'"

	return 0
}

helper::_push(){
	declare -n __=$1
	__+=("$2")
}

helper::_pop(){
	declare -n __=$1
	__=("${__[@]:0:$((${#__[@]}-1))}")
}

helper::_slice(){
	declare -n __=$1
	__=("${__[@]:$2:$3}")
}

helper::_shift(){
	declare -n __=$1
	__=("${__[@]:1}")
}

helper::_print(){
	declare -n __=$1
	echo ${__[*]}
}

helper::_println(){
	declare -n __=$1
	local e
	printf '%s\n' "${__[@]}"
}

helper::_uc(){
	declare -n __=$1
	__=("${__[@]^^${2:-*}}")
}

helper::_ucfirst(){
	declare -n __=$1
	__=("${__[@]^${2:-*}}")
}

helper::_lc(){
	declare -n __=$1
	__=("${__[@],,${2:-*}}")
}

helper::_lcfirst(){
	declare -n __=$1
	__=("${__[@],${2:-*}}")
}

helper::_sum(){
	declare -n __=$1
	echo $((${__[@]/%/+}0))
}

helper::_prefix(){
	declare -n __=$1
	__=("${__[@]%$2*}")
}

helper::_prefixfirst(){
	declare -n __=$1
	__=("${__[@]%%$2*}")
}

helper::_suffix(){
	declare -n __=$1
	__=("${__[@]#*$2}")
}

helper::_suffixfirst(){
	declare -n __=$1
	__=("${__[@]##*$2}")
}

helper::_substring(){
	declare -n __=$1
	local i
	for i in "${!__[@]}"; do
		__[$i]=${__[$i]:$2:$3}
	done
}

helper::_replace(){
	declare -n __=$1
	__=("${__[@]/${2:-*}/$3}")
}

helper::_replaceprefix(){
	declare -n __=$1
	__=("${__[@]/#${2:-*}/$3}")
}

helper::_replacesuffix(){
	declare -n __=$1
	__=("${__[@]/${2:-*}%/$3}")
}

helper::_uniq(){
	declare -n __=$1
	declare -A ___
	local i
	for i in "${!__[@]}"; do
		___[${__[$i]}]=1
	done
	__=("${!___[@]}")
}

helper::_sort(){
	declare -n __=$1
	mapfile -t __ < <(printf '%s\n' "${__[@]}" | sort -V)
}

helper::_test(){
	declare -a arr
	helper::loadaliases -v arr
	arr.push "foo foo"
	arr.push bar
	arr.push zar
	arr.print
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
