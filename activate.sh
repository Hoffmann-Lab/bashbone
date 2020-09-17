#! /usr/bin/env bash
# (c) Konstantin Riege
set -o pipefail

unset ERROR
trap 'e=$?; echo ":ERROR: ${ERROR:-"..an unexpected one"} (exit $e) @ $(basename $0) (line: $LINENO) $BASH_COMMAND" >&2; exit $e' ERR
ERROR="$(basename "$0") script needs to be sourced"
[[ "${BASH_SOURCE[0]}" != "$0" ]]

trap 'trap - ERR; trap - RETURN' RETURN
trap 'e=$?; echo ":ERROR: ${ERROR:-"..an unexpected one"} (exit $e) @ $(basename ${BASH_SOURCE[0]}) (line: $LINENO) $BASH_COMMAND" >&2; return $e' ERR

ERROR="requires a bash shell"
[[ "$(ps -p $$ -o command= | cut -d ' ' -f 1)" == "bash" ]]
ERROR="unsupported operating system"
[[ $OSTYPE == "linux" ]]
ERROR="requieres bash >= v4.4"
[[ ${BASH_VERSINFO[0]} -gt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -ge 4) ]]

unset ERROR
INSDIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
INSDIR_TOOLS="$(dirname "$INSDIR")"
unset OPTIND activate
while getopts :i:c: arg; do
	case $arg in
		i) INSDIR_TOOLS="$OPTARG";;
		c) activate="$OPTARG";;
		:) ERROR="argument missing"; false;;
	esac
done

_IFS=$IFS
IFS=$'\n'
for f in "$INSDIR/lib/"*.sh; do
	source "$f"
done
IFS=$_IFS

ERROR="environment setup failed (use -c false to disable tools and conda activation)"
configure::environment -i "$INSDIR_TOOLS" -b "$INSDIR" -c ${activate:-false}
[[ $activate ]] || {
	commander::printinfo {COMMANDER[0]}<<- EOF
		to activate conda environment do
		source $(basename "${BASH_SOURCE[0]}") -c true
	EOF
}

bashbone(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			Bashbone

			A bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses.

			Usage:
			-h | help
			-l | list functions for users
			-e | list functions for experienced users
			-d | list functions for developers
			-a | list all, including non-standalone, functions
		EOF
		return 0
	}

	local OPTIND arg
	while getopts 'hleda' arg; do
		case $arg in
		h) _usage; return 0;;
		l) declare -F | grep -oE '\S+::\S+' | grep -vF -e ::_ -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		e) declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		d) declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper::_ | grep -F -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		a) declare -F | grep -oE '\S+::\S+' | sort -t ':' -k1,1 -k3,3V; return 0;;
		*) _usage; return 1;;
		esac
	done
	_usage
	return 0
}

return 0
