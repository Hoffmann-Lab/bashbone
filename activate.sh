#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$(basename "${BASH_SOURCE[0]}") script needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && return 1
[[ $OSTYPE != "linux" ]] && echo "unsupported operating system" >&2 && return 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && return 1

declare -F bashbone &> /dev/null && bashbone -x
mapfile -t BASCHBONE_BAK_SHOPT < <(shopt | sed -E '/off$/d;{s/^(\S+).+/shopt -s \1/}')
mapfile -t BASCHBONE_BAK_ERR < <(trap -p ERR)
mapfile -t BASCHBONE_BAK_RET < <(trap -p RETURN)

BASHBONE_DIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
toolsdir="$(dirname "$BASHBONE_DIR")"
unset OPTIND activate
while getopts :i:c: arg; do
	case $arg in
		i) toolsdir="$OPTARG";;
		c) activate="$OPTARG";;
		:) echo "argument missing" >&2; return 1;;
	esac
done

_IFS=$IFS
IFS=$'\n'
for f in "$BASHBONE_DIR/lib/"*.sh; do
	source "$f" || {
		echo "file not found $f" >&2
		return 1
	}
done
IFS=$_IFS
BASHBONE_VERSION=$version

BASHBONE_ERROR="environment setup failed. use -c false to disable tools and conda activation"
configure::environment -i "$toolsdir" -c ${activate:-false}

 [[ $activate ]] || {
	commander::printinfo {COMMANDER[0]}<<- EOF
		to activate conda environment do
		source $(basename "${BASH_SOURCE[0]}") -c true
	EOF
}

unset BASHBONE_ERROR

bashbone(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			Bashbone v$BASHBONE_VERSION

			A bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses.

			Usage:
			-h | help
			-l | list functions for users
			-e | list functions for experienced users
			-d | list functions for developers
			-a | list all, including non-standalone, functions
			-x | exit bashbone and revert changes to environment
		EOF
		return 0
	}

	local OPTIND arg
	while getopts 'hledax' arg; do
		case $arg in
		h)	_usage; return 0;;
		l)	declare -F | grep -oE '\S+::\S+' | grep -vF -e ::_ -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		e)	declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		d)	declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper::_ | grep -F -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		a)	declare -F | grep -oE '\S+::\S+' | sort -t ':' -k1,1 -k3,3V; return 0;;
		x)	shopt -u extdebug
			set +E +o pipefail +o functrace
			trap - ERR
			trap - RETURN
			source <(printf '%s\n' "${BASCHBONE_BAK_RET[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_ERR[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_SHOPT[@]}")

			local x
			for x in $(declare -p  | cut -d ' ' -f 3 | grep -oE '^BASHBONE[^=]+'); do
				unset $x
			done
			for x in $(declare -F | grep -oE '\S+::\S+') bashbone; do
				unset -f $x
			done
			return 0
		;;
		*) _usage; return 1;;
		esac
	done
	_usage
	return 0
}

return 0
