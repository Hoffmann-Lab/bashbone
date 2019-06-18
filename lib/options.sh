#! /usr/bin/env bash
# (c) Konstantin Riege

options::usage(){
	cat <<- EOF
		SYNOPSIS
		$(basename $0) -i [all|upgrade|[<tool>,..]] -d [path]

		OPTIONS
		-i | --install [all|upgrade|[<tool>,..]] : install into given directory
		-d | --directory [path]                  : installation path 
		-t | --threads [value]                   : threads - predicted default: $THREADS
		-l | --log [path]                        : log file - default: [-d]/install.log
		-v | --verbose                           : enable verbose mode
		-h | --help                              : prints this message

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 0
}


options::checkopt(){
	local arg=false
	case $1 in
		-h | --help) usage;;
		-v | --verbose) VERBOSITY=2;;
		-d | --directory) arg=true; INSDIR="$2";;
		-t | --threads) arg=true; THREADS=$2;;
		-l | --log) arg=true; LOG=$2;;
		-i | --install) arg=true; mapfile -d ',' -t INSTALL <<< $2;;
		-*) echo ":ERROR: illegal option $1"; return 1;; 
		*) echo ":ERROR: illegal option $2"; return 1;;
	esac
	$arg && {
		[[ ! $2 ]] && echo ":ERROR: argument missing for option $1" && return 1
		[[ "$2" =~ ^- ]] && echo ":ERROR: illegal argument $2 for option $1" && return 1
		return 0
	} || {
		[[ $2 ]] && [[ ! "$2" =~ ^- ]] && echo ":ERROR: illegal argument $2 for option $1" && return 1
		return 0
	}
}

options::parse(){
	[[ $# -eq 0 ]] && options::usage
	[[ $# -eq 1 ]] && [[ ! $1 =~ ^- ]] && echo "illegal option $1" >&2 && return 1
	for i in $(seq 1 $#); do
		if [[ ${!i} =~ ^- ]]; then
			j=$((i+1))
			options::checkopt "${!i}" "${!j}" || return 1
		else 
			((++i))
		fi
	done

	return 0
}
