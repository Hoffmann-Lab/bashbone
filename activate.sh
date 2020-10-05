#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$(basename "${BASH_SOURCE[0]}") script needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && return 1
[[ $OSTYPE != "linux" ]] && echo "unsupported operating system" >&2 && return 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && return 1

declare -F bashbone &> /dev/null && bashbone -x
BASHBONE_BAK_PATH="$PATH"
mapfile -t BASCHBONE_BAK_SHOPT < <(shopt | sed -E '/off$/d;{s/^(\S+).+/shopt -s \1/}')
mapfile -t BASCHBONE_BAK_ERR < <(trap -p ERR)
mapfile -t BASCHBONE_BAK_RET < <(trap -p RETURN)

BASHBONE_WORKDIR="$PWD"
BASHBONE_DIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
BASHBONE_TOOLSDIR="$(dirname "$BASHBONE_DIR")"
unset OPTIND
while getopts ':i:c:x:h' arg; do
	case $arg in
		i)	BASHBONE_TOOLSDIR="$OPTARG";;
		c)	BASHBONE_CONDA="$OPTARG";;
		x)	BASHBONE_EXITFUN="$OPTARG";;
		:)	echo "argument missing" >&2; return 1;;
		h)	cat <<- EOF
				This is bashbone activation script.

				Afterwards bashbone functions can be called. To see lists of available functions, execute bashbone -h
				To revert changes made to the current shell environment, execute bashbone -x

				Usage:
				-h            | this help
				-i <path>     | to installation root <path>/latest/<tools>/<bins>
				                default: inferred from script location
				                hint: run activation from source code, indeed enables basic functions, but will fail on executing tools
				-c <activate> | true/false conda from [-i]/conda/bin
				                default: false
				-x <fun>      | a function or command to be called upon EXIT signal with the exit code appended to its argument list
				                default: none

				Example:
				source activate.sh -i <path> -c true -x "<fun> [<arg>..]"
				bashbone -h
				bashbone -x
			EOF
			return 0
		;;
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

${BASHBONE_CONDA:-false} && {
	BASHBONE_ERROR="conda activation failed. use -c false to disable conda activation"
	source $BASHBONE_TOOLSDIR/conda/bin/activate base &> /dev/null
	commander::printinfo "utilizing $(conda --version)"
}
set -E -o pipefail -o functrace # -E allows simple trap bubbeling and -o functrace enables inheritance of RETURN and DEBUG trap
shopt -s extdebug # do not shopt -u extdebug, otherwise set -E -o pipefail configuration will be nuked
shopt -s extglob
shopt -s expand_aliases
ulimit -n $(ulimit -Hn)
export MALLOC_ARENA_MAX=4
BASHBONE_ERROR="environment setup failed. use -i [path] to point towards a custom installation directory"
td="$(readlink -e "$BASHBONE_TOOLSDIR/latest")"
# better stay off custom java path to avoid conflicts with conda openjdk and pre-compiled jars requesting specific versions (IncompatibleClassChangeError)
# [[ $td && -e $td/java ]] && export JAVA_HOME=$(dirname $(readlink -e $td/java))
[[ $td ]] && export PATH="$(readlink -e "$td/"!(java) | xargs -echo | sed 's/ /:/g'):$PATH"
export PATH="$(readlink -e "$BASHBONE_DIR/scripts" | xargs -echo | sed 's/ /:/g'):$PATH"
if [[ $BASHBONE_EXITFUN ]]; then
	trap 'configure::exit -p $$ -f $BASHBONE_EXITFUN $?' EXIT
else
	trap 'configure::exit -p $$' EXIT
fi
# make use of local scope during trace to trigger tmp file deletion etc.
trap 'declare -F _cleanup::${FUNCNAME[0]} &> /dev/null && _cleanup::${FUNCNAME[0]}' RETURN
# check execution from terminal (use return only) or a script/subshell (use return and exit)
# error traps must not be splited into muliple lines to hold correct lineno
if ! [[ ${BASH_EXECUTION_STRING} ]]; then # prefer over [[ $- =~ i ]] which is true in bash -i -c 'echo $-'
	# since trap needs to persist in shell, make sure return is triggerd only from sourced bashbone functions. otherwise there will be issues with bash completion, vte and other functions
	trap 'e=$?; if [[ $e -ne 141 ]]; then if [[ ${BASH_SOURCE[0]} && "$(cd "$BASHBONE_WORKDIR"; readlink -e "${BASH_SOURCE[0]}")" =~ "$BASHBONE_DIR" ]]; then configure::err -x $e -e "$BASHBONE_ERROR" -l $LINENO -f ${FUNCNAME[0]} -w "$BASHBONE_WORKDIR"; return $e; fi; fi' ERR
else
	# if activate used within script or subshell, traps will be nuked anyways, thus allow tracing for all functions
	# dont call exit directly. allow for back trace through all functions. local scopes are available
	trap 'e=$?; if [[ $e -ne 141 ]]; then if [[ "${BASH_SOURCE[0]}" == "$0" && "${FUNCNAME[0]}" == "main" ]]; then configure::err -x $e -e "$BASHBONE_ERROR" -l $LINENO -s "$0" -w "$BASHBONE_WORKDIR"; exit $e; else configure::err -x $e -e "$BASHBONE_ERROR" -l $LINENO -f ${FUNCNAME[0]} -w "$BASHBONE_WORKDIR"; return $e; fi; fi' ERR
fi
trap 'BASHBONE_ERROR="killed"' INT TERM

[[ $BASHBONE_CONDA ]] || {
	commander::printinfo {COMMANDER[0]}<<- EOF
		to activate conda environment execute
		source $(basename "${BASH_SOURCE[0]}") -c true
	EOF
}

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
			PATH="$BASHBONE_BAK_PATH"

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

unset BASHBONE_ERROR
return 0
