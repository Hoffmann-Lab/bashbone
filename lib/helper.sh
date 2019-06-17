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
