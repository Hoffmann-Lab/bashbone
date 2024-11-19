#! /usr/bin/env bash
# (c) Konstantin Riege

function options::usage(){
	commander::print {COMMANDER[0]}<<- EOF
		DESCRIPTION
		Bashbone setup routine

		SYNOPSIS
		setup.sh -i [all|upgrade] -d [path]

		OPTIONS
		-i | --install [lite|all|upgrade] : install into given directory
		-g | --use-config                 : use supplied yaml files and URLs instead of cutting edge tools
		-d | --directory [path]           : installation path
		-t | --threads [value]            : threads - predicted default: $THREADS
		-l | --log [path]                 : log file - default: [-d]/install.log
		-v | --verbose                    : enable verbose mode
		-h | --help                       : prints this message

		DEVELOPER OPTIONS
		-s | --source [path,..]           : source file(s) to overload compile::[lite|all|upgrade|<tool>] functions
		-i | --install [<tool>,..]        : install into given directory

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	return 1
}

function options::checkopt(){
	local arg=false
	case $1 in
		-h | --help) options::usage || exit 0;;
		-v | --verbose) VERBOSITY=2;;
		-d | --directory) arg=true; INSDIR="$2";;
		-t | --threads) arg=true; THREADS=$2;;
		-l | --log) arg=true; LOG="$2";;
		-g | --use-config) USECONFIG=true;;
		-s | --source)
			arg=true
			declare -a mapdata
			mapfile -t -d ',' mapdata < <(printf '%s' "$2")
			for f in "${mapdata[@]}"; do
				source $f || {
					commander::printerr "unexpected error in source code - please contact developer"
					return 1
				}
			done
		;;
		-i | --install) arg=true; mapfile -t -d ',' INSTALL < <(printf '%s' "$2");;
		-*) commander::printerr "illegal option $1"; return 1;;
		*) commander::printerr "illegal option $2"; return 1;;
	esac
	if $arg; then
		[[ ! $2 ]] && commander::printerr "argument missing for option $1" && return 1
		[[ "$2" =~ ^- ]] && commander::printerr "illegal argument $2 for option $1" && return 1
	else
		[[ $2 ]] && [[ ! "$2" =~ ^- ]] && commander::printerr "illegal argument $2 for option $1" && return 1
	fi

	return 0
}

function options::parse(){
	[[ $# -eq 0 ]] && { options::usage || exit 0; }
	[[ $# -eq 1 ]] && [[ ! $1 =~ ^- ]] && commander::printerr "illegal option $1" && return 1
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
