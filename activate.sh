#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$(basename "${BASH_SOURCE[0]}") script needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && return 1
[[ $OSTYPE != "linux" ]] && echo "unsupported operating system" >&2 && return 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && return 1


declare -F bashbone &> /dev/null && bashbone -x
BASHBONE_BAK_PATH="$PATH"
mapfile -t BASCHBONE_BAK_SHOPT < <(shopt | sed -E '/off$/d;{s/^(\S+).+/shopt -s \1/}')
# mapfile -t BASCHBONE_BAK_ALIAS < <(alias)
# unalias -a
mapfile -t BASCHBONE_BAK_ERR < <(trap -p ERR)
mapfile -t BASCHBONE_BAK_RETURN < <(trap -p RETURN)
mapfile -t BASCHBONE_BAK_EXIT < <(trap -p EXIT)
mapfile -t BASCHBONE_BAK_INT < <(trap -p INT)
mapfile -t BASCHBONE_BAK_TERM < <(trap -p TERM)


BASHBONE_WORKDIR="$PWD"
export BASHBONE_DIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
export BASHBONE_TOOLSDIR="$(dirname "$BASHBONE_DIR")"
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


# +m turns off job control to avoid "done" aor "terminated" messages of subshells
set +m -o pipefail -o errtrace -o functrace # traces enable trap to know local scope of functions and shubshells that inherit ERR trap (-o errtrace) and RETURN and DEBUG trap (-o functrace)
shopt -s extdebug # do not shopt -u extdebug, otherwise set -E -o pipefail configuration will be nuked
shopt -s extglob
ulimit -n $(ulimit -Hn)
export MALLOC_ARENA_MAX=4


# trap '_trap_e=$?; pid=$BASHPID; echo error ${FUNCNAME[0]} $_trap_e $LINENO; if [[ ${FUNCNAME[0]} && "${FUNCNAME[0]}" != "main" ]]; then return $_trap_e; else kill -USR1 $$; fi; kill -PIPE $(pstree -Alps $pid | grep -oE "\($$).+\($pid\)" | grep -oE "[0-9]+" | tail -n +2)' ERR
# bubble-up/propagate exit code > 0 via return to continue ERR trap recursion on higher level function.
# use _trap_VAR to derease risk of interference with dynamically loaded variable name space of current function recursion level
# BASHPID set to last nested subshell pid keeps constant during trace -> contrast to LINENO which referes to recent BASH_LINENO[0] (in line to BASH_SOURCE[@]/FUNCNAME[@])
# send PIPE signal to all nested subshells (parent pids between $BASHPID and $$) to kill them silently -> no conflict with 141/PIPE catch in ERR trap, since kill alone does not invoke ERR, despite of trap being inherited, without calling wait $! to receive exit code e.g. cat <(kill -PIPE $BASHPID) && echo mustNotShow vs cat <(kill -PIPE $BASHPID) && wait $! && echo mustNotShow
# main as function name is not set from interactive shells or script execution as an other process like xargs -I {} bash -c '{}'
_bashbone_errtrap(){
	local err=$1 pid=$2 ppid=$3 line=$4 fun=$5 src=$6
	[[ $err -eq 141 ]] && return 0
	if [[ $fun && "$fun" != "main" ]]; then
		[[ "$BASHBONE_ERROR" != "false" && $line -gt 1 ]] && configure::err -x $err -e "$BASHBONE_ERROR" -l $line -f "$fun"
		if [[ $pid -ne $ppid ]]; then
			kill -SIGRTMIN+2 $ppid
			kill -PIPE $(pstree -Alps $pid | grep -oE "\($ppid).+\($pid\)" | grep -oE "[0-9]+" | tail -n +2) || true
		fi
		return $err
	else
		[[ "$BASHBONE_ERROR" != "false" && $line -gt 1 ]] && configure::err -x $err -e "$BASHBONE_ERROR" -l $line -s "$src" -w "$BASHBONE_WORKDIR"
		kill -SIGRTMIN+1 $ppid
		if [[ $pid -ne $ppid ]]; then
			kill -PIPE $(pstree -Alps $pid | grep -oE "\($ppid).+\($pid\)" | grep -oE "[0-9]+" | tail -n +2) || true
		fi
		return 0
	fi
}
trap '_trap_e=$?; _bashbone_errtrap $_trap_e $BASHPID $$ $LINENO "${FUNCNAME[0]}" "${BASH_SOURCE[0]}"; if [[ $? -gt 0 ]]; then return $_trap_e; fi' ERR
trap '_trap_cmd=$BASH_COMMAND; if [[ ! "$_trap_cmd" =~ ^source[[:space:]] ]]; then declare -f _cleanup::${FUNCNAME[0]} > /dev/null && _cleanup::${FUNCNAME[0]}; fi' RETURN # fire cleanup function from local name space upon any return (success or failure). do not fire cleanup if returning from a source command
trap '_trap_e=$?; BASHBONE_ERROR="false"; (exit $_trap_e)' INT TERM # silence error messages which go out of sync, when e.g. TERM signal is received from failed sibling xargs job. INT signal is often not correctly received and thus needs to trigger ERR manually. note that INT/ctr-c by user kills subshells immediately -> no cleanup for inner subshell functions possible
trap 'configure::exit -p $$ -f "$BASHBONE_EXITFUN" $?' EXIT # fire last, global cleanup function and kill all remaining subshells
trap 'exit 255' SIGRTMIN+1 # define user signal, that can be received by $$ from any subshell or other process. use 255, in case xargs job failed, to prevent xargs to load and execute further, queued commands
trap 'return 1' SIGRTMIN+2
trap 'BASHBONE_ERROR="$BASH_COMMAND: command not found"; (exit 127)' SIGRTMIN+3
if [[ $- =~ i ]]; then
	# do not exit on interactive shells ;)
	trap '_trap_e=$?; if [[ $LINENO -gt 0 && ${FUNCNAME[0]} && "$(cd "$BASHBONE_WORKDIR"; readlink -e "${BASH_SOURCE[0]}")" =~ "$BASHBONE_DIR" ]]; then _bashbone_errtrap $_trap_e $BASHPID $$ $LINENO "${FUNCNAME[0]}" "${BASH_SOURCE[0]}"; if [[ $? -gt 0 ]]; then return $_trap_e; fi; fi' ERR
	trap 'echo "$BASH_COMMAND: command not found" >&2; (exit 127)' SIGRTMIN+3
fi
command_not_found_handle(){ # overload bash builtin
	# echo "$*: command not found" >&2
	kill -SIGRTMIN+3 $$
	return 127
}



_IFS=$IFS
IFS=$'\n'
for f in "$BASHBONE_DIR/lib/"*.sh; do
	BASHBONE_ERROR="bashbone library file not found $f"
	source "$f"
done
IFS=$_IFS
BASHBONE_VERSION=$version


if ${BASHBONE_CONDA:-false}; then
	BASHBONE_ERROR="conda activation failed. use -c false to disable conda activation or use -i [path] to point towards a custom installation directory"
	source $BASHBONE_TOOLSDIR/conda/bin/activate bashbone
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
td="$(readlink -e "$BASHBONE_TOOLSDIR/latest")" && PATH="$(readlink -e "$td/"!(java) | xargs -echo | sed 's/ /:/g'):$PATH"
PATH="$(readlink -e "$BASHBONE_DIR/scripts" | xargs -echo | sed 's/ /:/g'):$PATH"
_bashbone_pathdedub() {
	local P PATHNEW
	declare -A PATHUNIQ
	declare -a PATHARR
	mapfile -t PATHARR < <(sed 's/::*/\n/g' <<< "$PATH")
	for P in "${PATHARR[@]}"; do
		${PATHUNIQ["$P"]:-false} || {
			 PATHUNIQ["$P"]=true
			 PATHNEW+=":$P"
		}
	done
	PATH="${PATHNEW:1}"
}
_bashbone_pathdedub


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
			-t | list tools and versions in setupped environment
			-x | exit bashbone and revert changes to environment
		EOF
		return 0
	}

	local OPTIND arg
	while getopts 'hrcsuvledtax' arg; do
		case $arg in
		h)	_usage; return 0;;
		r)	mdless -P $BASHBONE_DIR/README.md | less; return 0;;
		c)	source $BASHBONE_TOOLSDIR/conda/bin/activate bashbone &> /dev/null; return 0;;
		s)	while [[ -n $CONDA_PREFIX ]]; do conda deactivate &> /dev/null; done; return 0;;
		u)	find -L $BASHBONE_TOOLSDIR/latest -maxdepth 2 -name "*.sh" -not -name "activate.sh" -not -name "setup.sh" -printf "%f\n"
			find "$BASHBONE_DIR/scripts/" -type f -name "*.pl" -printf "%f\n" -o -name "*.sh" -printf "%f\n" | rev | sort | rev
			return 0
		;;
		v)	find -L $BASHBONE_TOOLSDIR/latest -maxdepth 2 -name "*.sh" -not -name "activate.sh" -not -name "setup.sh" -printf "%f\n"
			find "$BASHBONE_DIR/scripts/" -type f -printf "%f\n" | rev | sort | rev
			return 0
		;;
		l)	declare -F | grep -oE '\S+::\S+' | grep -vF -e ::_ -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		e)	declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		d)	declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper::_ | grep -F -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		a)	declare -F | grep -oE '\S+::\S+' | sort -t ':' -k1,1 -k3,3V; return 0;;
		t)	(	source $BASHBONE_TOOLSDIR/conda/bin/activate base
				mapfile -t mapdata < <({ readlink -e "$BASHBONE_TOOLSDIR"/latest/* | sed -nE 's@.*\/([^/]+)-([0-9][^/]+)\/*.*@\L\1\t\2@p;'; conda list | grep -vE '^(#|lib|perl-|xorg-|r-|python-|font)' | tr -s ' ' '\t'; } | sort -k1,1 -k2,2Vr | cut -f 1,2)
				for e in $(conda env list | grep -F "$BASHBONE_TOOLSDIR" | grep -v '^base' | cut -f 1 -d ' '); do
					conda list -n $e | grep -vE '^(#|lib|perl-|xorg-|r-|python-|font)'
				done | tr -s ' ' '\t' | sort -k1,1 -k2,2Vr | cut -f 1,2 | rev | uniq -f 1 | rev | grep -v -F -f <(printf "%s\n" "${mapdata[@]}" | cut -f 1) | sort -k1,1 - <(printf "%s\n" "${mapdata[@]}")
			)
			return 0
		;;
		x)	shopt -u extdebug
			set -m +E +o pipefail +o functrace
			trap - RETURN
			trap - ERR
			trap - EXIT
			trap - INT
			trap - TERM
			source <(printf '%s\n' "${BASCHBONE_BAK_RETURN[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_ERR[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_EXIT[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_INT[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_TERM[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_SHOPT[@]}")
			# source <(printf '%s\n' "${BASCHBONE_BAK_ALIAS[@]}")
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
