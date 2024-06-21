#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$0 needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && exit 1
[[ "$(uname)" != "Linux" ]] && echo "unsupported operating system" >&2 && exit 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && exit 1

export BASHBONE_DIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
export BASHBONE_EXTENSIONDIR="${BASHBONE_EXTENSIONDIR:-$BASHBONE_DIR}/lib"

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
			return 0
		;;
	esac
done
unset OPTIND
[[ $- =~ i ]] && source "${BASH_SOURCE[0]}" -h && return 0

if [[ ! $- =~ i ]] && ${BASHBONE_SETSID:-true}; then
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
BASHBONE_PGID=$(($(ps -o pgid= -p $BASHPID)))

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

export TMPDIR="$(command mktemp -d -p "${TMPDIR:-/tmp}" XXXXXXXXXX)"
function _bashbone_reset(){
	trap "" INT TERM ERR

	{ env kill -INT -- -$BASHBONE_PGID & wait $!; } &> /dev/null
	sleep 0.1
	{ env kill -TERM -- -$BASHBONE_PGID & wait $!; } &> /dev/null

	mkdir "$TMPDIR/.lock.cleanup" &> /dev/null && {
		find "$TMPDIR" -type f -name "cleanup.*.sh" -exec bash -c 'tac "$1" | bash || true; rm -f "$1"' bash {} \;
		"${BASHBONE_EXITFUN:-:}" $1
		rm -rf "$TMPDIR"
	}
	return 0
}

function _bashbone_wrapper(){
	local BASHBONE_FUNCNAME=$1
	shift
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
		if [[ $fun && "$fun" != "main" ]]; then
			read -r f l src < <(declare -F "$fun")
		fi
		cmd=$(awk -v l=$line '{ if(NR>=l){if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{print o$0; exit}}else{if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{o=""}}}' "$src" | sed -E -e 's/\s+/ /g' -e 's/(^\s+|\s+$)//g')
		echo ":ERROR:$BASHBONE_ERROR in ${src:-shell} (function: ${fun:-main}) @ line $line: $cmd" >&2
	done
	echo ":ERROR: exit code 143" >&2
	return 0
}

function _bashbone_on_error(){
	[[ $1 -eq 141 ]] && return 0
	trap "" INT
	[[ $1 -ne 130 && $1 -ne 143 ]] && mkdir "$TMPDIR/.lock.error" &> /dev/null && _bashbone_trace $1
	set -m
	{ env kill -TERM -- -$BASHBONE_PGID; sleep 0.1; _bashbone_reset 143; } &
	set +m
	wait $!
	return 143
}

function _bashbone_on_return(){
	if [[ $1 -eq 0 && "$2" != "source" && "$2" != "." && "$BASHBONE_FUNCNAME" == "${FUNCNAME[1]}" ]]; then
		if [[ -e "$BASHBONE_CLEANUP" ]]; then
			tac "$BASHBONE_CLEANUP" | bash || true
			rm -f "$BASHBONE_CLEANUP"
		fi
	fi
	return 0
}

trap '_bashbone_reset $?' EXIT
trap '_bashbone_on_error $? || return 143' ERR
trap '(exit 130)' INT
trap 'exit 143' TERM
trap '_bashbone_on_return $? "${BASH_COMMAND%% *}"' RETURN

if [[ -e "$BASHBONE_EXTENSIONDIR" ]]; then
	while read -r f; do
		alias $f="_bashbone_wrapper $f"
	done < <(
		shopt -s extdebug
		while read -r f; do
			source "$f"
		done < <(find -L "$BASHBONE_EXTENSIONDIR" -name "*.sh")
		for f in $(declare -F | awk '{print $3}'); do
			read -r f l s < <(declare -F $f)
			[[ "$s" == "$BASHBONE_EXTENSIONDIR"* ]] && echo $f
		done
	)
	while read -r f; do
		source "$f"
	done < <(find -L "$BASHBONE_EXTENSIONDIR" -name "*.sh")
fi
