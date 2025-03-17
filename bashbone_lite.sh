#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$0 needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && exit 1
[[ "$(uname)" != "Linux" ]] && echo "unsupported operating system" >&2 && exit 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && exit 1

export BASHBONE_WORKINGDIR="$PWD"
export BASHBONE_DIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
export BASHBONE_EXTENSIONDIR="${BASHBONE_EXTENSIONDIR:-$BASHBONE_DIR}/lib"
export BASHBONE_LEGACY=${BASHBONE_LEGACY:-false}

unset OPTIND
while getopts 'p:x:s:ah' arg; do
	case $arg in
		p)	TMPDIR="$OPTARG";;
		x)	BASHBONE_EXITFUN="$OPTARG";;
		s)	BASHBONE_EXTENSIONDIR="$OPTARG";;
		a)	shift $((OPTIND-1)); break;;
		h)	cat <<- 'EOF'
				This is bashbone-lite.

				Designed to be used in a bash script, it allows for
				  - error tracing
				  - proper subprocess termination upon exit or error
				  - cleanup of temporary files upon function return, exit or error

				Usage:
				-h            | this help
				-p <path>     | to temporary directory. default: $TMPDIR or /tmp
				-x <string>   | a function or command string to be called upon EXIT signal. this function or command will receive the exit code as last positional argument
				-a <optarg>   | use as last option in case bashbone is sourced in a script wich makes use of optargs

				Howto:
				1) Activate
				  # source at the very top of your script! activate.sh re-runs the script under a new process group
				    #! /usr/bin/env bash
				    source bashbone_lite.sh -a "$@"

				2) Add function specific cleanup commands to be executed upon return or error
				2a) Wrap function manually
				  # wrap your function using an alias and thereby assign a cleanup script to it (access via $BASHBONE_CLEANUP)
				    alias myfun="_bashbone_wrapper myfun"

				  # now define your function (function keyword is required!)
				    function myfun(){
					  # and either append custom commands to $BASHBONE_CLEANUP script
					    echo "echo 'cleanup myfun'" >> "$BASHBONE_CLEANUP"
					    touch /tmp/tmp.123456.remove.me1
					    echo "rm -rf '/tmp/tmp.123456.remove.me1'" >> "$BASHBONE_CLEANUP"
					  # or use mktemp
					    mktemp -p /tmp --suffix=".remove.me2"
				    }
				2b) Wrap functions from a custom library to be sourced from a <path/to/library>/lib directory
				  # source bashbone_lite and hand over the library parent directory
				    #! /usr/bin/env bash
				    source bashbone_lite.sh -s <path/to/library> -a "$@"

				  # see 2a) for how to use function specific cleanups

				3) Run a custom cleanup function upon exit or error
				  # source bashbone_lite and hand over a cleanup function
				    #! /usr/bin/env bash
				    source bashbone_lite.sh -x cleanup -a "$@"

				    function cleanup(){
				    	echo "exit code is: $1"
				        rm -f "$mytempfile"
				    }
			EOF
			unset OPTIND
			return 0
		;;
	esac
done
unset OPTIND
[[ $- == *i* ]] && source "${BASH_SOURCE[0]}" -h && return 0

if ${BASHBONE_SETSID:-false}; then
	export BASHBONE_SETSID=false
	if [[ -t 0 ]]; then
		{ trap 'exit 130' INT; exec setsid --wait bash "$(realpath -s "$0")" "$@"; } &
	else
		cat | { trap 'exit 130' INT; exec setsid --wait bash "$(realpath -s "$0")" "$@"; } &
	fi
	BASHBONE_PGID=$!
	trap 'trap "" INT TERM; env kill -INT -- -$BASHBONE_PGID' INT TERM
	wait $BASHBONE_PGID
	wait $BASHBONE_PGID
	BASHBONE_EX=$?
	trap - EXIT
	exit $BASHBONE_EX
fi
export -n BASHBONE_SETSID

mapfile -t BASHBONE_BAK_SHOPT < <(shopt | awk '$2=="off"{print "shopt -u "$1}'; shopt | awk '$2=="on"{print "shopt -s "$1}')
BASHBONE_BAK_TRAPS=$(trap -p)
mapfile -t BASHBONE_BAK_SET < <(printf "%s" $- | sed 's/[isc]//g' | sed -E 's/(.)/set -\1\n/g')
mapfile -t BASHBONE_BAK_ALIASES < <(declare -p BASH_ALIASES | sed 's/^declare/declare -g/') # squeeze in global paramater otherwise _bashbone_reset function call declares BASH_ALIASES locally
mapfile -t BASHBONE_BAK_CMD_NOT_FOUND < <(declare -f command_not_found_handle)

set -o pipefail -o errtrace -o functrace
shopt -s expand_aliases extdebug extglob globstar
[[ $BASH_VERSINFO -gt 4 ]] && shopt -u localvar_inherit
ulimit -n $(ulimit -Hn)

function command_not_found_handle(){
	echo "$1: command not found" >&2
	(exit 127)
	return 127
}

function mktemp(){
	local tmp="$(command mktemp "$@")"
	echo "rm -rf '$tmp'" >> "${BASHBONE_CLEANUP:-/dev/null}"
	echo "$tmp"
	return 0
}

BASHBONE_PGID=$(($(ps -o pgid= -p $BASHPID)))
export TMPDIR="$(command mktemp -d -p "${TMPDIR:-/tmp}" bashbone.XXXXXXXXXX)"
trap '_bashbone_on_error 130 $LINENO; exit 130 &> /dev/null' INT
trap '_bashbone_on_error $? $LINENO || { [[ $BASHPID -ne $BASHBONE_PGID ]] && { return 130 &> /dev/null || exit 130; } || exit 130 &> /dev/null; }' ERR
trap '_bashbone_on_exit $?' EXIT
trap '_bashbone_on_return $? "${BASH_COMMAND%% *}"' RETURN

function _bashbone_on_exit(){
	trap "" INT TERM ERR

	if [[ $1 -eq 0 ]] && ! [[ -t 1 && -t 2 ]]; then
		:
	else
		if [[ $$ -eq $BASHBONE_PGID ]]; then
			env kill -INT -- -$BASHBONE_PGID
			sleep 0.1
			env kill -TERM -- -$BASHBONE_PGID
		fi
	fi
	find -L "$TMPDIR" -type f -name "cleanup.*.sh" -exec bash -c 'tac "$1" | bash || true; rm -f "$1"' bash {} \;
	"${BASHBONE_EXITFUN:-:}" $1
	rm -rf "$TMPDIR"

	set +o pipefail +o errtrace +o functrace
	trap - INT TERM RETURN ERR EXIT
	source <(IFS=; echo "$BASHBONE_BAK_TRAPS")
	source <(printf '%s\n' "${BASHBONE_BAK_CMD_NOT_FOUND[@]}")
	source <(printf '%s\n' "${BASHBONE_BAK_SHOPT[@]}")
	source <(printf '%s\n' "${BASHBONE_BAK_SET[@]}")
	source <(printf '%s\n' "${BASHBONE_BAK_ALIASES[@]}")
	eval "$(trap -p EXIT | sed 's/^trap -- .//;s/. EXIT$//')"
	return 0
}

function _bashbone_wrapper(){
	local BASHBONE_FUNCNAME=$1
	shift
	local f_bashbone_wrapper l_bashbone_wrapper s_bashbone_wrapper BASHBONE_LEGACY=$BASHBONE_LEGACY
	read -r f_bashbone_wrapper l_bashbone_wrapper s_bashbone_wrapper < <(declare -F $BASHBONE_FUNCNAME) || true
	if [[ "$s_bashbone_wrapper" == "$BASHBONE_DIR/lib/"* ]]; then
		[[ "$BASHBONE_FUNCNAME" == commander::* || "$BASHBONE_FUNCNAME" == progress::* ]] || BASHBONE_LEGACY=true
	else
		BASHBONE_LEGACY=false
	fi
	local BASHBONE_CLEANUP="$(command mktemp -p "$TMPDIR" cleanup.XXXXXXXXXX.sh)"
	$BASHBONE_FUNCNAME "$@"
	return $?
}

function _bashbone_trace(){
	local error=$1 line l fun f src frame last=$((${#BASH_LINENO[@]}-2)) cmd
	for frame in $(seq 1 $last); do
		fun=${FUNCNAME[$((frame+1))]}
		[[ "$fun" == "_bashbone_wrapper" || "$fun" == "command_not_found_handle" ]] && continue
		src="${BASH_SOURCE[$((frame+1))]}"
		line=${BASH_LINENO[$frame]}
		if [[ $line -eq 1 ]]; then
			echo ":ERROR: ${BASHBONE_ERROR:+$BASHBONE_ERROR }in ${src:-shell} (function: ${fun:-main})" >&2
		else
			if [[ $fun && "$fun" != "main" ]]; then
				read -r f l src < <(declare -F "$fun")
			fi
			if [[ "$src" == "environment" ]]; then
				line=$((line+3))
				cmd=$(declare -f $fun | awk -v l=$line '{ if(NR>=l){if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{print o$0; exit}}else{if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{o=""}}}' | sed -E -e 's/\s+/ /g' -e 's/(^\s+|\s+$)//g')
			else
				cmd=$(awk -v l=$line '{ if(NR>=l){if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{print o$0; exit}}else{if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{o=""}}}' "$BASHBONE_WORKINGDIR/$src" | sed -E -e 's/\s+/ /g' -e 's/(^\s+|\s+$)//g')
			fi
			echo ":ERROR: ${BASHBONE_ERROR:+$BASHBONE_ERROR }in ${src:-shell} (function: ${fun:-main}) @ line $line: $cmd" >&2
		fi
	done
	return 0
}

function _bashbone_on_error(){
	[[ $1 -eq 141 ]] && return 0
	trap "" INT
	mkdir "$TMPDIR/.killed" &> /dev/null && _bashbone_trace $1 && env kill -INT -- -$BASHBONE_PGID
	return 130
}

function _bashbone_on_return(){
	if [[ $1 -eq 0 && "$2" != "source" && "$2" != "." && "$BASHBONE_FUNCNAME" == "${FUNCNAME[1]}" ]]; then
		COMMANDER=()
		[[ -e "$BASHBONE_CLEANUP" ]] && {
			tac "$BASHBONE_CLEANUP" | bash || true
			rm -f "$BASHBONE_CLEANUP"
		}
	fi
	return 0
}

if [[ -e "$BASHBONE_EXTENSIONDIR" ]]; then
	declare -x -a BASHBONE_FUNCNAMES
	mapfile -t BASHBONE_FUNCNAMES < <(
		shopt -s extdebug
		trap - RETURN
		while read -r f; do
			unset -f $f
		done < <(compgen -A function)
		while read -r f; do
			source "$f"
		done < <(find -L "$BASHBONE_EXTENSIONDIR" -name "*.sh" -not -name "#*")
		compgen -A function
	)
	for f in "${BASHBONE_FUNCNAMES[@]}"; do
		alias $f="_bashbone_wrapper $f"
	done

	while read -r f; do
		source "$f"
	done < <(find -L "$BASHBONE_EXTENSIONDIR" -name "*.sh" -not -name "#*")
fi
