#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$(basename "${BASH_SOURCE[0]}") script needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && return 1
[[ "$(uname)" != "Linux" ]] && echo "unsupported operating system" >&2 && return 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && return 1

declare -F bashbone &> /dev/null && bashbone -y
mapfile -t BASCHBONE_BAK_SHOPT < <(shopt | awk '$2=="off"{print "shopt -u "$1}'; shopt | awk '$2=="on"{print "shopt -s "$1}')
mapfile -t BASCHBONE_BAK_TRAPS < <(trap -p)
mapfile -t BASCHBONE_BAK_SET < <(printf "%s" $- | sed 's/[is]//g' | sed -E 's/(.)/set -\1\n/g')
mapfile -t BASCHBONE_BAK_ALIASES < <(declare -p BASH_ALIASES | sed 's/^declare/declare -g/') # squeeze in global paramater otherwise bashbone -x function call declares BASH_ALIASES locally
BASHBONE_BAK_TMPDIR="$TMPDIR"

export BASHBONE_DIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
export BASHBONE_TOOLSDIR="${BASHBONE_TOOLSDIR:-$(dirname "$BASHBONE_DIR")}"
unset OPTIND
# while getopts ':i:p:c:x:a:hj' arg; do
while getopts 'i:p:c:x:ahj' arg; do
	case $arg in
		i)	BASHBONE_TOOLSDIR="$OPTARG";;
		c)	BASHBONE_CONDA="$OPTARG";;
		x)	BASHBONE_EXITFUN="$OPTARG";;
		j)  BASHBONE_JOBCONTROL=true;;
		p)	TMPDIR="$OPTARG";;
		# a)	shift $((OPTIND-2)); break;; # check shift IND-1 or IND-2 via optarg if activate is called with -a "$*"
		# :)	if [[ "$OPTARG" == "a" ]]; then # use shift IND-1 if no (:optional i.e OPTARG becomes arg) arg given for -a if activate is called with -a "$@"
		# 		shift $((OPTIND-1))
		# 		break
		# 	else
		# 		echo "argument missing" >&2
		# 	fi
		# 	return 1
		# 	;;
		a) 	shift $((OPTIND-1)); break;;
		h)	cat <<- 'EOF'
				This is bashbone activation script.

				Afterwards bashbone functions can be called. To see lists of available functions, execute bashbone -h
				To revert changes made to the current shell environment, execute bashbone -x

				Usage:
				-h            | this help
				-i <path>     | to installation root <path>/latest/<tools>/<bins>
				                default: inferred from script location
				                hint: run activation from source code, indeed enables basic functions, but will fail on executing tools
				-p <path>     | to temporary directory. default: $TMPDIR (fallback: /tmp)
				-c <activate> | true/false conda from [-i]/conda/bin
				                default: false
				-j            | experimental: enable job control when interactive
				-x <fun>      | an optional function or command to be called upon EXIT signal. this function or command will receive the exit code as last argument
				-a <optarg>   | use as last option in case bashbone is sourced in a script wich makes use of optargs

				Example:
				source activate.sh -i <path> -c true -x "<fun> [<arg>..]" -a "$@"
				bashbone -h
				bashbone -x
			EOF
			return 0
		;;
	esac
done
export TMPDIR="${TMPDIR:-/tmp}" && mkdir -p "$TMPDIR" || return 1
export BASHBONE_PGID=$(($(ps -o pgid= -p $$)))
if ! [[ $- =~ i ]] && ${BASHBONE_SETSID:-true}; then
	export BASHBONE_SETSID=false
	exec bash -c 'trap "kill -TERM -- -\$BASHBONE_PGID;" INT TERM; setsid --wait env --default-signal=INT,QUIT bash "$0" "$@" & BASHBONE_PGID=$!; wait $BASHBONE_PGID' "$0" "$@"
fi


# https://www.gnu.org/software/bash/manual/html_node/The-Set-Builtin.html
set +a +b +f +k +n +t +u +v +x -B +C -H +P +o posix
# +m turns off job control so that a) SIGINT is ignored in asynchronous commands and b) subshells will not run in own process groups anymore i.e. setsid ensures all subprocesses and forks to run under same process group id
# thus, killing the process group $$ upon exit a) terminates all subshells and b) waiting for asynchronous subshells does not print "done" or "terminated" messages anymore
set +m -o pipefail -o errtrace -o functrace # traces enable trap to know local scope of functions and shubshells that inherit ERR trap (-E|-o errtrace) and RETURN and DEBUG trap (-T|-o functrace)
# https://www.gnu.org/software/bash/manual/html_node/The-Shopt-Builtin.html
shopt -u nocaseglob nocasematch xpg_echo
[[ $BASH_VERSINFO -gt 4 ]] && shopt -u localvar_inherit
# extdebug required by declare -F to get src and line; note that this option implicitly set -E -T and thus shopt -u extdebug nukes those $SHELLOPTS
shopt -s expand_aliases extdebug extglob extquote promptvars sourcepath

ulimit -n $(ulimit -Hn)
# export MALLOC_ARENA_MAX=4 # limit number of parallel memory pools by malloc() to e.g. reduce jvm memory allocation. use for single-threaded and parallellized java instances only to avoid side effects of drastically lowered performance of other multithreaded applications which make use of malloc() e.g. segemehl
# [[ $BASH_VERSINFO -lt 5 ]] && enable -n command # always use env - command in bash 4 is buggy : trap '/bin/echo > /dev/null' DEBUG; command cat <(echo foo)

if [[ $- =~ i ]]; then
	_bashbone_traps(){
		# do not fire on non-bashbone functions
		trap '_trap_e=$?; if [[ $_trap_e -eq 254 ]] || [[ ${FUNCNAME[0]} && "$(readlink -e "${BASH_SOURCE[0]}")" =~ "$BASHBONE_DIR" ]]; then _bashbone_errtrap $_trap_e || return 254; else trap - ERR; (exit $_trap_e); fi' ERR
		# do not exit programmatically but cleanup processes
		# trap 'trap "" INT TERM; rm -f "/dev/shm/BASHBONE_CLEANED.$$"; unset BASHBONE_ERROR BASHBONE_CLEANED; configure::exit -x 0 -p $$ -f "$BASHBONE_EXITFUN"' SIGRTMIN+1
		trap 'trap "" INT TERM ERR; rm -f "/dev/shm/BASHBONE_CLEANED.$$"; unset BASHBONE_ERROR BASHBONE_CLEANED; "${BASHBONE_EXITFUN:-:}" 0; { env kill -INT -- -$$; sleep 0.2; env kill -TERM -- -$$ & wait $!; } &> /dev/null' SIGRTMIN+1
		trap 'if [[ ${FUNCNAME[0]} && ${FUNCNAME[0]} != "main" ]]; then read -r BASHBONE_CLEANED < "/dev/shm/BASHBONE_CLEANED.$$"; return 254; else env kill -SIGRTMIN+1 $$; fi' SIGRTMIN+2
		# in order to gain back job control when interactive, detect functions via DEBUG trap, disable DEBUG to not recurse (re-enable upon return) and set +/-m prior to execution
		# might be extended by a more general pre-env setup to e.g. spare func src check in err trap above
		# attention: if code within DEBUG trap throws error, DEBUG returns without executing bash_command
		if ${BASHBONE_JOBCONTROL:=false}; then
			trap '_trap_cmd="${BASH_COMMAND%% *}"; set -m; if declare -f "$_trap_cmd" &> /dev/null; then trap - DEBUG; read -r _trap_fun _trap_line _trap_src < <(declare -F "$_trap_cmd"); if [[ "$(readlink -e "$_trap_src")" =~ "$BASHBONE_DIR" ]]; then set +m; fi; fi' DEBUG
			trap '_trap_cmd="${BASH_COMMAND%% *}"; if [[ "$_trap_cmd" != "source" && "$BASHBONE_CLEANED" != "${FUNCNAME[0]}" ]]; then declare -f _cleanup::${FUNCNAME[0]} > /dev/null && _cleanup::${FUNCNAME[0]}; fi; if [[ ${#FUNCNAME[@]} -eq 1 && ! $COMP_LINE ]]; then _bashbone_traps; fi' RETURN
		else
			trap '_trap_cmd="${BASH_COMMAND%% *}"; if [[ "$_trap_cmd" != "source" && "$BASHBONE_CLEANED" != "${FUNCNAME[0]}" ]]; then declare -f _cleanup::${FUNCNAME[0]} > /dev/null && _cleanup::${FUNCNAME[0]}; fi' RETURN
		fi
		# due to set -m no need to check for bashbone function
		trap '_trap_e=$?; trap "" INT; BASHBONE_ERROR="false"; (exit $_trap_e)' INT TERM
	}
	_bashbone_traps
	PROMPT_COMMAND="${PROMPT_COMMAND:+$PROMPT_COMMAND$'\n'}_bashbone_traps"
else
	trap '_trap_e=$?; _bashbone_errtrap $_trap_e || return 254' ERR
	# define user signal, that can be received by $$ from any subshell or other process. use 255, in case xargs job failed, to prevent xargs to load and execute further, queued commands
	trap 'exit 255' SIGRTMIN+1
	# ignore TERM on exit to prevent ERR/EXIT loop upon kill pgid
	# trap '_trap_e=$?; trap "" INT TERM ERR; rm -f "/dev/shm/BASHBONE_CLEANED.$$"; configure::exit -x $_trap_e -p $$ -f "$BASHBONE_EXITFUN"' EXIT
	trap '_trap_e=$?; trap "" INT TERM ERR; rm -f "/dev/shm/BASHBONE_CLEANED.$$"; unset BASHBONE_ERROR BASHBONE_CLEANED; "${BASHBONE_EXITFUN:-:}" $_trap_e; { env kill -INT -- -$$; sleep 0.2; env kill -TERM -- -$$ & wait $!; } &> /dev/null' EXIT

	trap 'if [[ ${FUNCNAME[0]} && ${FUNCNAME[0]} != "main" ]]; then read -r BASHBONE_CLEANED < "/dev/shm/BASHBONE_CLEANED.$$"; return 254; else exit 255; fi' SIGRTMIN+2
	# bubble-up/propagate exit code > 0 via return to continue ERR trap recursion on higher level function.
	# fire cleanup function from local name space upon any return (success or failure). do not fire cleanup if returning from a source command. cleanup may use TERM
	trap '_trap_cmd="${BASH_COMMAND%% *}"; if [[ "$_trap_cmd" != "source" && "$BASHBONE_CLEANED" != "${FUNCNAME[0]}" ]]; then declare -f _cleanup::${FUNCNAME[0]} > /dev/null && _cleanup::${FUNCNAME[0]}; fi' RETURN
	# silence error messages which go out of sync, when e.g. TERM signal is received from failed sibling xargs job.
	trap '_trap_e=$?; trap "" INT; BASHBONE_ERROR="false"; (exit $_trap_e)' INT TERM
fi

_bashbone_errtrap(){
	# process substitutions are asynchronouse subshells, and thus needs wait for exit code or programmatically terminated upon ERR
	# send PIPE signal to all nested subshells (parent pids between $BASHPID and $$) to kill them silently. PIPE is ignored in ERR trap due to e.g. samtools view <bam> | head ... whre samtools sends signal 141 to head
	# cat <(kill -PIPE $BASHPID) && echo mustNotShow vs cat <(kill -PIPE $BASHPID) && wait $! && echo mustNotShow
	# main as function name is not set from interactive shells or script execution as an other process like xargs -I {} bash -c '{}'
	# us underscore variables to avoid clash of this local scope mixed up with local scope of function to be returned from upon error!!!
	local _trap_e=$1 _trap_pid=$BASHPID _trap_ppid=$$
	[[ $_trap_e -eq 141 ]] && return 0

	declare -a _trap_pids=($(pstree -Alps $_trap_pid | grep -oE "\($_trap_ppid\).+\($_trap_pid\)" | grep -oE "[0-9]+" | grep -vFw $_trap_ppid))

	if [[ $_trap_e -eq 127 ]]; then # from command_not_found_handle, which is always a subshell

		_bashbone_trace 2 $_trap_e # start from 2 since 1 is command_not_found_handle and 0 this trap function

		echo "${FUNCNAME[2]}" > "/dev/shm/BASHBONE_CLEANED.$_trap_ppid" # if function to be terminated is not part of this (nested) subshell(s), prevents return at $$ to trigger cleanup again under different scope
		env kill -SIGRTMIN+2 $_trap_ppid
		declare -f _cleanup::${FUNCNAME[2]} > /dev/null && _cleanup::${FUNCNAME[2]} >&2 # clean under $BASHPID to have variables in local scope
		env kill -PIPE ${_trap_pids[@]} >& /dev/null || true
	else
		if [[ $_trap_e -ne 254 && "$BASHBONE_ERROR" != "false" ]]; then # first execution of trap determined by 254
			_bashbone_trace 1 $_trap_e
		fi

		if [[ ${FUNCNAME[1]} && ${FUNCNAME[1]} != "main" ]]; then
			if [[ $_trap_ppid -ne $_trap_pid ]]; then # in subshell
				echo "${FUNCNAME[1]}" > "/dev/shm/BASHBONE_CLEANED.$_trap_ppid"
				env kill -SIGRTMIN+2 $_trap_ppid
				declare -f _cleanup::${FUNCNAME[1]} > /dev/null && _cleanup::${FUNCNAME[1]} >&2 # clean under $BASHPID
				env kill -PIPE ${_trap_pids[@]} >& /dev/null || true
			fi
			return 254
		else
			env kill -SIGRTMIN+1 $_trap_ppid # exit
			if [[ $_trap_ppid -ne $_trap_pid ]]; then # in subshell
				env kill -PIPE ${_trap_pids[@]} >& /dev/null || true
			fi
		fi
	fi

	return 0
}

_bashbone_trace(){
	local line fun src frame=$1 e=$2 last=$((${#BASH_LINENO[@]}-2))
	# f=0
	# while caller $f >&2; do
	# 	((++f))
	# done
	#while read -r line fun src < <(caller $frame || true); do
	for frame in $(seq $frame $last); do
		[[ $frame -eq $last ]] && e=255
		line=${BASH_LINENO[$frame]}
		fun=${FUNCNAME[$((frame+1))]}
		src=${BASH_SOURCE[$((frame+1))]}
		[[ $fun == "command_not_found_handle" ]] && continue
		#echo $line $fun $src $$ $BASHPID >&2
		if [[ $fun && "$fun" != "main" ]]; then
			configure::err -x $e -e "$BASHBONE_ERROR" -l $line -f "$fun"
		else
			configure::err -x $e -e "$BASHBONE_ERROR" -l $line -s "$src"
		fi
		BASHBONE_ERROR="traceback"
		((++frame))
	done

	return 0
}

command_not_found_handle(){ # override bash builtin. runs in own subshell
	echo "$*: command not found" >&2
	BASHBONE_ERROR="command not found"
	(exit 127)
	return 127 # wont be reached
}


# ensure $f to be absolute so is BASH_SOURCE
for f in "$BASHBONE_DIR/lib/"*.sh; do
	BASHBONE_ERROR="bashbone library file not found $f"
	source "$f"
done
BASHBONE_VERSION=$version


if ${BASHBONE_CONDA:-false}; then
	BASHBONE_ERROR="conda activation failed. use -c false to disable conda activation or use -i [path] to point towards a custom installation directory"
	source "$BASHBONE_TOOLSDIR/conda/bin/activate" bashbone
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

	return 0
}

if [[ ! $BASHBONE_PATH || ! "$PATH" == *"$BASHBONE_PATH"* ]]; then
	if [[ -e "$BASHBONE_TOOLSDIR/latest" ]]; then
		export BASHBONE_PATH="$(realpath -s "$BASHBONE_TOOLSDIR/latest/"!(java|bashbone) | xargs -echo | sed 's/ /:/g')"
	fi
	if [[ -e "$BASHBONE_DIR/scripts" ]]; then
		export BASHBONE_PATH="$(realpath -s "$BASHBONE_DIR/scripts" | xargs -echo | sed 's/ /:/g')${BASHBONE_PATH:+":$BASHBONE_PATH"}"
	fi
	if [[ -e "$BASHBONE_DIR/tools" ]]; then
		export BASHBONE_PATH="$(realpath -s "$BASHBONE_DIR"/tools/*/bin | xargs -echo | sed 's/ /:/g')${BASHBONE_PATH:+":$BASHBONE_PATH"}"
	fi
	PATH="${BASHBONE_PATH:+"$BASHBONE_PATH:"}$PATH"
	_bashbone_pathdedub
fi

bashbone(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			Bashbone v$BASHBONE_VERSION

			A bash library for workflow and pipeline design within but not restricted to the scope of Next Generation Sequencing (NGS) data analyses.

			Usage:
			-h | help
			-r | readme
			-j | experimental: enables or disables job control when interactive
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
	while getopts 'hrjcsuvledtayx' arg; do
		case $arg in
		h)	_usage; return 0;;
		r)	mdless -P "$BASHBONE_DIR/README.md" | less; return 0;;
		j)	$BASHBONE_JOBCONTROL && BASHBONE_JOBCONTROL=false || BASHBONE_JOBCONTROL=true; return 0;;
		c)	source "$BASHBONE_TOOLSDIR/conda/bin/activate" bashbone &> /dev/null; return 0;;
		s)	while [[ -n $CONDA_PREFIX ]]; do conda deactivate &> /dev/null; done; return 0;;
		u)	[[ -e "$BASHBONE_TOOLSDIR/latest" ]] && find -L "$BASHBONE_TOOLSDIR/latest" -maxdepth 2 -name "*.sh" -not -name "activate.sh" -not -name "setup.sh" -printf "%f\n"
			find "$BASHBONE_DIR/scripts/" -type f -name "*.pl" -printf "%f\n" -o -name "*.sh" -printf "%f\n" | rev | sort | rev
			return 0
		;;
		v)	[[ -e "$BASHBONE_TOOLSDIR/latest" ]] && find -L "$BASHBONE_TOOLSDIR/latest" -maxdepth 2 -name "*.sh" -not -name "activate.sh" -not -name "setup.sh" -printf "%f\n" || true
			find "$BASHBONE_DIR/scripts/" -type f -printf "%f\n" | rev | sort | rev
			return 0
		;;
		l)	declare -F | grep -oE '\S+::\S+' | grep -vF -e ::_ -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		e)	declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		d)	declare -F | grep -oE '\S+::\S+' | grep -vF -e compile:: -e helper::_ | grep -F -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V; return 0;;
		a)	declare -F | grep -oE '\S+::\S+' | sort -t ':' -k1,1 -k3,3V; return 0;;
		t)	(	source "$BASHBONE_TOOLSDIR/conda/bin/activate" base
				mapfile -t mapdata < <({ readlink -e "$BASHBONE_TOOLSDIR/latest/"* | sed -nE 's@.*\/([^/]+)-([0-9][^/]+)\/*.*@\L\1\t\2@p;'; conda list | grep -vP '^(#|_|lib|perl-|xorg-|r-(?!(base|wgcna))|r\s|python-|font|gcc_|gxx_|gfortran_|ca-certificates|pkg-config|pthread)' | tr -s ' ' '\t'; } | sort -k1,1 -k2,2Vr | cut -f 1,2)
				for e in $(conda env list | grep -F "$BASHBONE_TOOLSDIR" | grep -v '^base' | cut -f 1 -d ' '); do
					conda list -n $e | grep -vP '^(#|_|lib|perl-|xorg-|r-(?!(base|wgcna))|r\s|python-|font|gcc_|gxx_|gfortran_|ca-certificates|pkg-config|pthread)'
				done | tr -s ' ' '\t' | sort -k1,1 -k2,2Vr | cut -f 1,2 | rev | uniq -f 1 | rev | grep -v -Fw -f <(printf "%s\n" "${mapdata[@]}" | cut -f 1) | sort -k1,1 - <(printf "%s\n" "${mapdata[@]}")
			)
			return 0
		;;
		y)	set +m +o pipefail +o errtrace +o functrace
			trap - RETURN ERR EXIT INT TERM DEBUG SIGRTMIN+1 SIGRTMIN+2
			source <(printf '%s\n' "${BASCHBONE_BAK_TRAPS[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_SHOPT[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_SET[@]}")
			source <(printf '%s\n' "${BASCHBONE_BAK_ALIASES[@]}")
			[[ $BASHBONE_BAK_TMPDIR ]] && TMPDIR="$BASHBONE_BAK_TMPDIR"
			[[ $BASHBONE_PATH ]] && PATH="${PATH/$BASHBONE_PATH:/}"
			[[ "$PROMPT_COMMAND" == "_bashbone_traps" ]] && unset PROMPT_COMMAND || PROMPT_COMMAND="${PROMPT_COMMAND/$'\n'_bashbone_traps/}"
			return 0
		;;
		x)	local f l s x
			for f in $(declare -F | cut -d ' ' -f 3 | grep -vFx bashbone); do
				read -r f l s < <(declare -F $f)
				[[ "$(readlink -e "$s")" =~ "$BASHBONE_DIR" ]] && unset -f $f
			done

			bashbone -y
			unset -f bashbone

			for x in $(declare -p | cut -d ' ' -f 3 | grep -oE '^BASHBONE[^=]+'); do
				unset $x
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


################### test1

source ./activate.sh -c false

fun(){
	local x=ffun
	_cleanup::fun(){
		echo cleaning fun $x
	}
	## must not fail. all nested functions will be cleaned
	echo foo
	## simple fail
	cat qwert
	## cmd not found handler in subshell
	qwert
	## cmd in subshell + cmd not found handler in subshell
	if qwert; then
		echo if fin
	fi
	echo fun fin
}

inter(){
	local x=finter
	_cleanup::inter(){
		echo cleaning inter $x
	}
	cat $1
}

## error in main
#fun; echo wrap fin1
#x=$(fun; echo wrap fin1 >&2)
#cat <(inter <(fun); echo wrap fin1)
#echo main fin
#exit

## error in function
wrap(){
	local x=fwrap
	_cleanup::wrap(){
		echo cleaning wrap $x
	}
	## simple case
	fun; echo wrap fin1
	## subsehll case 0
	#(fun; echo wrap fin1 >&2)
	## subshell case 1
	#y=$(fun; echo wrap fin1 >&2)
	## subsehll case 2
	#cat <(cat <(fun); echo wrap fin1)
	## subshell case 3. cannot catch/cleanup inter upon error. also caller/FUNCNAME[*] does not know about inter
	#cat <(inter <(fun); echo wrap fin1)
	echo wrap fin2
}

outa(){
	wrap
}

outa
echo main fin


################### test2

source ./activate.sh -c false

cat <<-'EOF' > funs.sh
fun(){
	declare -a cmds
	cmds=('{ echo running fun under $BASHPID $$ $PPID $(($(ps -o pgid= -p $$))); qwert; echo fun1; } & wait $!')
	cmds+=('{ echo running fun under $BASHPID $$ $PPID $(($(ps -o pgid= -p $$))); sleep 2; qwert; echo fun2; } & wait $!')
	cmds+=('{ echo running fun under $BASHPID $$ $PPID $(($(ps -o pgid= -p $$))); sleep 2; qwert; echo fun3; } & wait $!')
	cmds+=('{ echo running fun under $BASHPID $$ $PPID $(($(ps -o pgid= -p $$))); sleep 2; qwert; echo fun4; } & wait $!')
	commander::runcmd -i 4 -a cmds
}
EOF

wrap(){
	declare -a cmds
	cmds=('source ./funs.sh; echo running wrap under $BASHPID $$ $PPID $(($(ps -o pgid= -p $$))); sleep 1; fun; echo HIER1')
	cmds+=('source ./funs.sh; echo running wrap under $BASHPID $$ $PPID $(($(ps -o pgid= -p $$))); sleep 5; fun; echo HIER2')
	commander::runcmd -i 2 -b -a cmds
}
echo starting $(($(ps -o pgid= -p $$)))
wrap
