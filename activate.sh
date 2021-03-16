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
mapfile -t BASCHBONE_BAK_EXIT < <(trap -p EXIT)

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

set -o pipefail -o errtrace -o functrace # traces enable trap to know local scope of functions and shubshells that inherit ERR trap (-o errtrace) and RETURN and DEBUG trap (-o functrace)
shopt -s extdebug # do not shopt -u extdebug, otherwise set -E -o pipefail configuration will be nuked
shopt -s extglob
shopt -s expand_aliases
ulimit -n $(ulimit -Hn)
export MALLOC_ARENA_MAX=4
BASHBONE_ERROR="environment setup failed. use -i [path] to point towards a custom installation directory"


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


if ${BASHBONE_CONDA:-false}; then
	BASHBONE_ERROR="conda activation failed. use -c false to disable conda activation"
	source $BASHBONE_TOOLSDIR/conda/bin/activate base &> /dev/null
else
	[[ $BASHBONE_CONDA ]] || {
		commander::printinfo {COMMANDER[0]}<<- EOF
			to activate conda environment either execute: bashbone -c
			or execute: source $(basename "${BASH_SOURCE[0]}") -c true
		EOF
	}
fi


# better stay off custom java path to avoid conflicts with conda openjdk and pre-compiled jars requesting specific versions (IncompatibleClassChangeError)
# [[ $td && -e $td/java ]] && export JAVA_HOME=$(dirname $(readlink -e $td/java))
td="$(readlink -e "$BASHBONE_TOOLSDIR/latest")"
[[ $td ]] && export PATH="$(readlink -e "$td/"!(java) | xargs -echo | sed 's/ /:/g'):$PATH"
export PATH="$(readlink -e "$BASHBONE_DIR/scripts" | xargs -echo | sed 's/ /:/g'):$PATH"


if [[ $BASHBONE_EXITFUN ]]; then
	trap 'e=$?; sleep $((++BASHBONE_TRAPPED==1?1:0)); configure::exit -p $$ -f $BASHBONE_EXITFUN $e' EXIT
else
	trap 'sleep $((++BASHBONE_TRAPPED==1?1:0)); configure::exit -p $$' EXIT
fi
# make use of local scope during functrace to trigger tmp file deletion etc.
trap 'declare -F _cleanup::${FUNCNAME[0]} &> /dev/null && _cleanup::${FUNCNAME[0]}' RETURN
# error traps must not be splited into muliple lines to hold correct lineno
if [[ $- =~ i ]]; then
	# since trap needs to persist in shell, make sure return is triggerd only from sourced bashbone functions. otherwise there will be issues with bash completion, vte and other functions
	# although errtrace ignores errors upon 'while until for if || &&' in interactive shell it still triggers ERR trap on function definition line i.e. at LINENO==0 where error occured
	trap 'e=$?; if [[ $e -ne 141 ]]; then if [[ $LINENO -gt 0 && ${FUNCNAME[0]} && "$(cd "$BASHBONE_WORKDIR"; readlink -e "${BASH_SOURCE[0]}")" =~ "$BASHBONE_DIR" ]]; then sleep $((++BASHBONE_TRAPPED==1?1:0)); configure::err -x $e -e "$BASHBONE_ERROR" -l $LINENO -f "${FUNCNAME[0]}" -w "$BASHBONE_WORKDIR"; return $e; else unset BASHBONE_TRAPPED; fi; fi' ERR
else
	# xargs creates funcname main for executing script. to prevent xargs in commander to continue despite of an error, set exit code to 255 so that xargs will stop immediately
	trap 'e=$?; if [[ $e -ne 141 ]]; then e=255; if [[ ${FUNCNAME[0]} && "${FUNCNAME[0]}" != "main" ]]; then sleep $((++BASHBONE_TRAPPED==1?1:0)); configure::err -x $e -e "$BASHBONE_ERROR" -l $LINENO -f "${FUNCNAME[0]}" -w "$BASHBONE_WORKDIR"; return $e; else configure::err -x $e -e "$BASHBONE_ERROR" -l $LINENO -s "${BASH_SOURCE[0]}" -w "$BASHBONE_WORKDIR"; exit $e; fi; fi' ERR
fi
trap 'BASHBONE_ERROR="killed"' INT TERM


bashbone(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			Bashbone v$BASHBONE_VERSION

			A bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses.

			Usage:
			-h | help
			-r | readme
			-c | activate conda if installed
			-s | stop conda if activated
			-u | list scripts for users
			-v | list scripts for experienced users
			-l | list functions for users
			-e | list functions for experienced users
			-d | list functions for developers
			-a | list all, including non-standalone, functions
			-x | exit bashbone and revert changes to environment
		EOF
		return 0
	}

	local OPTIND arg
	while getopts 'hrcsuvledax' arg; do
		case $arg in
		h)	_usage; return 0;;
		r)	mdless -P $BASHBONE_DIR/README.md | less; return 0;;
		c)	source $BASHBONE_TOOLSDIR/conda/bin/activate base &> /dev/null;	commander::printinfo "utilizing $(conda --version)"; return 0;;
		s)	while [[ -n $CONDA_PREFIX ]]; do conda deactivate &> /dev/null; done; return 0;;
		u)	find -L $BASHBONE_TOOLSDIR/latest -maxdepth 2 -name "*.sh" -not -name "activate.sh" -not -name "setup.sh" -printf "%f\n"; find "$BASHBONE_DIR/scripts/" -type f -name "*.pl" -printf "%f\n" -o -name "*.sh" -printf "%f\n" | rev | sort | rev; return 0;;
		v)	find -L $BASHBONE_TOOLSDIR/latest -maxdepth 2 -name "*.sh" -not -name "activate.sh" -not -name "setup.sh" -printf "%f\n"; find "$BASHBONE_DIR/scripts/" -type f -printf "%f\n" | rev | sort | rev; return 0;;
		l)	declare -F | grep -oE '\S+::\S+' | grep -vF -e ::_ -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		e)	declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		d)	declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper::_ | grep -F -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		a)	declare -F | grep -oE '\S+::\S+' | sort -t ':' -k1,1 -k3,3V; return 0;;
		x)	shopt -u extdebug
			set +E +o pipefail +o functrace
			trap - ERR
			trap - RETURN
			trap - EXIT
			trap - INT
			trap - TERM
			source <(printf '%s\n' "${BASCHBONE_BAK_RET[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_ERR[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_EXIT[@]}")
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