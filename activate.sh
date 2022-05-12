#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$(basename "${BASH_SOURCE[0]}") script needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && return 1
[[ "$(uname)" != "Linux" ]] && echo "unsupported operating system" >&2 && return 1
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

export BASHBONE_DIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
export BASHBONE_TOOLSDIR="${BASHBONE_TOOLSDIR:-$(dirname "$BASHBONE_DIR")}"
unset OPTIND
while getopts ':i:c:x:a:h' arg; do
	case $arg in
		i)	BASHBONE_TOOLSDIR="$OPTARG";;
		c)	BASHBONE_CONDA="$OPTARG";;
		x)	BASHBONE_EXITFUN="$OPTARG";;
		a)  shift $((OPTIND-2)); break;;
		h)	cat <<- 'EOF'
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

# +m turns off job control to avoid "done" or "terminated" messages when waiting for asynchronouse subshells because subshells will not run in own process groups anymore
set +m -o pipefail -o errtrace -o functrace # traces enable trap to know local scope of functions and shubshells that inherit ERR trap (-o errtrace) and RETURN and DEBUG trap (-o functrace)
shopt -s extdebug # do not shopt -u extdebug, otherwise set -E -o pipefail configuration will be nuked
shopt -s extglob
ulimit -n $(ulimit -Hn)
export MALLOC_ARENA_MAX=4
enable -n kill # either disable bash builtin kill here or use always "env kill"


# bubble-up/propagate exit code > 0 via return to continue ERR trap recursion on higher level function.
# use _trap_VAR to derease risk of interference with dynamically loaded variable name space of current function recursion level
# BASHPID set to last nested subshell pid keeps constant during trace -> contrast to LINENO which referes to recent BASH_LINENO[0] (in line to BASH_SOURCE[@]/FUNCNAME[@])
# send PIPE signal to all nested subshells (parent pids between $BASHPID and $$) to kill them silently -> no conflict with 141/PIPE catched in ERR trap, since kill alone does not invoke ERR, despite of trap being inherited, without calling wait $! to receive exit code e.g. cat <(kill -PIPE $BASHPID) && echo mustNotShow vs cat <(kill -PIPE $BASHPID) && wait $! && echo mustNotShow
# main as function name is not set from interactive shells or script execution as an other process like xargs -I {} bash -c '{}'
if [[ $BASHBONE_PGID ]]; then
	if [[ $$ -eq $BASHBONE_PGID ]]; then
		# fire last, global cleanup function and kill all remaining subshells
		trap '_trap_e=$?; trap "" INT TERM; rm -f "/dev/shm/BASHBONE_CLEANED.$$"; configure::exit -x $_trap_e -p $$ -f "$BASHBONE_EXITFUN"' EXIT
	else
		trap '_trap_e=$?; trap "" INT TERM; rm -f "/dev/shm/BASHBONE_CLEANED.$$"; configure::exit_job -x $_trap_e -p $$ -f "$BASHBONE_EXITFUN"' EXIT
		# runcmd without setsid: don't kill $$ == pgid, in case of error simply exit with 255 for further traceback
		# runcmd with setsid: kill -- -$$ does not work because $$ is not process group leader
	fi
else
	if [[ $$ -eq $(($(ps -o pgid= -p $$))) ]]; then
		export BASHBONE_PGID=$$
	else
		exec setsid --wait env bash "$0" "$@"
	fi
fi

_bashbone_settrap(){
	# silence error messages which go out of sync, when e.g. TERM signal is received from failed sibling xargs job. INT signal is often not correctly received and thus needs to trigger ERR manually. note that INT/ctr-c by user kills subshells immediately -> no cleanup for inner subshell functions possible
	trap '_trap_e=$?; trap "" INT TERM; BASHBONE_ERROR="false"; (exit $_trap_e)' INT TERM
}
_bashbone_settrap

if [[ $- =~ i ]]; then
	# do not run on non-bashbone functions. $LINENO -gt 0 && ... necessary?
	# || [[ ${FUNCNAME[0]} == "command_not_found_handle" && ${FUNCNAME[1]} && "${BASH_SOURCE[1]}" =~ "$BASHBONE_DIR" ]]
	trap '_trap_e=$?; if [[ $_trap_e -eq 254 ]] || [[ ${FUNCNAME[0]} && "${BASH_SOURCE[0]}" =~ "$BASHBONE_DIR" ]]; then _bashbone_errtrap $_trap_e || return 254; fi' ERR
	# restore INT TERM traps upon return
	trap '_trap_cmd=$BASH_COMMAND; if [[ ! "$_trap_cmd" =~ ^source[[:space:]] && "$BASHBONE_CLEANED" != "${FUNCNAME[0]}" ]]; then declare -f _cleanup::${FUNCNAME[0]} > /dev/null && _cleanup::${FUNCNAME[0]}; fi; _bashbone_settrap' RETURN
	# do not exit programmatically
	trap '_bashbone_settrap' SIGRTMIN+1
else
	trap '_trap_e=$?; _bashbone_errtrap $_trap_e || return 254' ERR
	trap '_trap_cmd=$BASH_COMMAND; if [[ ! "$_trap_cmd" =~ ^source[[:space:]] && "$BASHBONE_CLEANED" != "${FUNCNAME[0]}" ]]; then declare -f _cleanup::${FUNCNAME[0]} > /dev/null && _cleanup::${FUNCNAME[0]}; fi; _bashbone_settrap' RETURN # fire cleanup function from local name space upon any return (success or failure). do not fire cleanup if returning from a source command
	trap 'exit 255' SIGRTMIN+1 # define user signal, that can be received by $$ from any subshell or other process. use 255, in case xargs job failed, to prevent xargs to load and execute further, queued commands
fi
trap 'read -r BASHBONE_CLEANED < "/dev/shm/BASHBONE_CLEANED.$$"; return 254' SIGRTMIN+2

_bashbone_errtrap(){
	local _trap_e=$1 _trap_pid=$BASHPID _trap_ppid=$$ # avoid clash of this local scope mixed up with local scope of function to be returned from upon error!!!
	# e.g. samtools view <bam> | head ... samtools sends signal 141 to head, thus $? gets 141
	[[ $_trap_e -eq 141 ]] && return 0

	declare -a pids=($(pstree -Alps $_trap_pid | grep -oE "\($_trap_ppid\).+\($_trap_pid\)" | grep -oE "[0-9]+" | grep -vFw $_trap_ppid))

	if [[ $_trap_e -eq 127 ]]; then # from command_not_found_handle, which is always a subshell
		_bashbone_trace 2 $_trap_e # start from 2 since 1 is command_not_found_handle and 0 this trap function

		echo "${FUNCNAME[2]}" > "/dev/shm/BASHBONE_CLEANED.$_trap_ppid" # if function to be terminated is not part of this (nested) subshell(s), return at $$ triggers cleanup again
		if [[ ${FUNCNAME[3]} && ${FUNCNAME[3]} != "main" ]]; then
			env kill -SIGRTMIN+2 $_trap_ppid # return
		else
			env kill -SIGRTMIN+1 $_trap_ppid # exit
		fi
		declare -f _cleanup::${FUNCNAME[2]} > /dev/null && _cleanup::${FUNCNAME[2]} >&2 # clean under $BASHPID
		env kill -PIPE ${pids[@]} >& /dev/null || true
	else
		if [[ $_trap_e -ne 254 && "$BASHBONE_ERROR" != "false" ]]; then # first execution of trap determined by 254
			_bashbone_trace 1 $_trap_e
		fi

		if [[ ${FUNCNAME[1]} && ${FUNCNAME[1]} != "main" ]]; then
			if [[ $_trap_ppid -ne $_trap_pid ]]; then # in subshell
				echo "${FUNCNAME[1]}" > "/dev/shm/BASHBONE_CLEANED.$_trap_ppid"
				if [[ ${FUNCNAME[2]} && ${FUNCNAME[2]} != "main" ]]; then
					env kill -SIGRTMIN+2 $_trap_ppid # return
				else
					env kill -SIGRTMIN+1 $_trap_ppid # exit
				fi
				declare -f _cleanup::${FUNCNAME[1]} > /dev/null && _cleanup::${FUNCNAME[1]} >&2 # clean under $BASHPID
				env kill -PIPE ${pids[@]} >& /dev/null || true
			fi
			return 254
		else
			env kill -SIGRTMIN+1 $_trap_ppid # exit
			if [[ $_trap_ppid -ne $_trap_pid ]]; then # in subshell
				env kill -PIPE ${pids[@]} >& /dev/null || true
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
		[[ $frame -eq  $last ]] && e=255
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
	unset BASHBONE_ERROR

	return 0
}

command_not_found_handle(){ # override bash builtin. runs in own subshell
	BASHBONE_ERROR="$*: command not found"
	(exit 127)
	return 127 # wont be reached
}


# ensure $f to be absolute so is BASH_SOURCE
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
if [[ ! $PATH == *"$BASHBONE_TOOLSDIR/latest"* && -e "$BASHBONE_TOOLSDIR/latest" ]]; then
	PATH="$(realpath -s "$BASHBONE_TOOLSDIR/latest/"!(java|bashbone) | xargs -echo | sed 's/ /:/g'):$PATH"
fi
if [[ ! $PATH == *"$BASHBONE_DIR/scripts"* && -e "$BASHBONE_DIR/scripts" ]]; then
	PATH="$(realpath -s "$BASHBONE_DIR/scripts" | xargs -echo | sed 's/ /:/g'):$PATH"
fi
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
				mapfile -t mapdata < <({ readlink -e "$BASHBONE_TOOLSDIR"/latest/* | sed -nE 's@.*\/([^/]+)-([0-9][^/]+)\/*.*@\L\1\t\2@p;'; conda list | grep -vP '^(#|_|lib|perl-|xorg-|r-(?!(base|wgcna))|r\s|python-|font|gcc_|gxx_|gfortran_|ca-certificates|pkg-config|pthread)' | tr -s ' ' '\t'; } | sort -k1,1 -k2,2Vr | cut -f 1,2)
				for e in $(conda env list | grep -F "$BASHBONE_TOOLSDIR" | grep -v '^base' | cut -f 1 -d ' '); do
					conda list -n $e | grep -vP '^(#|_|lib|perl-|xorg-|r-(?!(base|wgcna))|r\s|python-|font|gcc_|gxx_|gfortran_|ca-certificates|pkg-config|pthread)'
				done | tr -s ' ' '\t' | sort -k1,1 -k2,2Vr | cut -f 1,2 | rev | uniq -f 1 | rev | grep -v -Fw -f <(printf "%s\n" "${mapdata[@]}" | cut -f 1) | sort -k1,1 - <(printf "%s\n" "${mapdata[@]}")
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
			# do not unset to keep pgid
			# for x in $(declare -p  | cut -d ' ' -f 3 | grep -oE '^BASHBONE[^=]+'); do
			# 	unset $x
			# done
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
	cmds=('{ echo running fun under $BASHPID $$ $PPID $(($(ps -o pgid= -p $$))); qwert; echo HIER3; } & wait $!')
	commander::runcmd -t 1 -a cmds
}
EOF

wrap(){
	declare -a cmds
	cmds=('source ./funs.sh; running wrap under $BASHPID $$ $PPID $(($(ps -o pgid= -p $$))); sleep 1; fun; echo HIER1')
	cmds+=('source ./funs.sh; running wrap under $BASHPID $$ $PPID $(($(ps -o pgid= -p $$))); sleep 5; fun; echo HIER2')

	cmds=("echo running wrap under \$((\$(ps -o pgid= -p \$\$))); source ./funs.sh; fun")
	cmds+=("echo running wrap under \$((\$(ps -o pgid= -p \$\$))); source ./funs.sh; sleep 5; fun")
	commander::runcmd -t 2 -b -a cmds
}
echo starting $(($(ps -o pgid= -p $$)))
wrap
