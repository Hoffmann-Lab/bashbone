#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$0 needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && return 1
[[ "$(uname)" != "Linux" ]] && echo "unsupported operating system" >&2 && return 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && return 1

#####################################################################################

export BASHBONE_WORKINGDIR="$PWD"
export BASHBONE_DIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
export BASHBONE_VERSION=$(source "$BASHBONE_DIR/lib/version.sh"; echo $version)
export BASHBONE_TOOLSDIR="${BASHBONE_TOOLSDIR:-$(dirname "$BASHBONE_DIR")}" # export to not override if previously defined
export BASHBONE_EXTENSIONDIR="${BASHBONE_EXTENSIONDIR:-$BASHBONE_DIR}" # directory of sh files to be sourced (functions will be wrapped) after bashbone
export BASHBONE_CONDA=${BASHBONE_CONDA:-false}
export BASHBONE_LEGACY=${BASHBONE_LEGACY:-false}

while getopts 'i:p:c:x:s:r:l:ah' arg; do
	case $arg in
		i)	BASHBONE_TOOLSDIR="$OPTARG";;
		p)	TMPDIR="$OPTARG";;
		c)	BASHBONE_CONDA="$OPTARG";;
		x)	BASHBONE_EXITFUN="$OPTARG";;
		s)	BASHBONE_EXTENSIONDIR="$OPTARG";;
		r)	BASHBONE_SETSID=${BASHBONE_SETSID:-$OPTARG};; # use exported state upon restart to not recurse
		l)	BASHBONE_LEGACY="$OPTARG";;
		a)	shift $((OPTIND-1)); break;;
		h)	cat <<- 'EOF'
				This is bashbone activation script.

				To see lists of available options and functions, source me and execute bashbone -h

				Usage:
				-h              | this help
				-l <legacymode> | true/false let commander inerts line breaks, thus crafts one-liners from makecmd here-documents
				                  default: false
				-i <path>       | to installation root <path>/latest
				                  default: inferred from script location
				                  hint: run activation from source code, indeed enables basic functions, but will fail on executing tools
				-p <path>       | to temporary directory. default: $TMPDIR or /tmp
				-c <activate>   | true/false conda from [-i]/conda/bin and prepend bashbone installed tools to PATH
				                  default: false
				-x <string>     | a function or command string to be called upon EXIT signal. this function or command will receive the exit code as last positional argument
				-s <path>       | extend/overload bashbone functions by libraries at <path>/lib and via -c option also add <path>/scripts to PATH
				-r <re-execute> | true/false the script that implements activate.sh under a new process group to not bubble up INT/TERM or kill a pipeline upon exit
				-a <optarg>     | use as last option in case bashbone is sourced in a script wich makes use of positional arguments (optargs)

				Example:
				source activate.sh -i <path> -c true -x "<fun> [<arg>..]" -a "$@"
				bashbone -h
			EOF
			unset OPTIND
			return 0
		;;
	esac
done
unset OPTIND

#####################################################################################

if [[ ! $- == *i* ]] && ${BASHBONE_SETSID:-false}; then
	# re-run under own pgid via setsid, to be able to kill process group when kill signal is send to pid only e.g. by SGE/PBE qdel or xargs
	# do not set sid when executed with gnu parallel!

	# exec bash -c 'export PGID=$$; trap "trap \"\" INT TERM; env kill -INT -- -\$PGID" INT; trap "trap \"\" INT TERM; env kill -TERM -- -\$PGID" TERM; setsid --wait env --default-signal=INT,QUIT bash "$0" "$@" & PGID=$!; wait $PGID' "$(realpath -s "$0")" "$@"
	export BASHBONE_SETSID=false
	if [[ -t 0 ]]; then
		{ trap 'exit 130' INT; exec setsid --wait bash "$(realpath -s "$0")" "$@"; } &
	else
		cat | { trap 'exit 130' INT; exec setsid --wait bash "$(realpath -s "$0")" "$@"; } &
	fi
	BASHBONE_PGID=$!
	# always trap INT to trigger ERR that triggers TERM to implement termination sequence with INT coming first to avoid most termination messages
	# termination sequence: INT,1000,TERM,0 plus error tracing and cleanup
	trap 'trap "" INT TERM; env kill -INT -- -$BASHBONE_PGID' INT TERM
	wait $BASHBONE_PGID # would be 130 upon INT and not exit code of script if script defines trap 'exit 42' INT
	wait $BASHBONE_PGID # would be 42 or again 130 unless script defines trap 'exit 42' INT
	BASHBONE_EX=$?
	trap - EXIT # due to eval of potential user exit trap after bashbone_reset in on_exit function
	#former without eval: [[ $BASHBONE_EX -gt 0 ]] && trap - EXIT # in case user defined EXIT trap, which must not be called here again
	exit $BASHBONE_EX
fi
export -n BASHBONE_SETSID # unset or remove export property, so that a truely desired restart of a script e.g. via qsub is realized

#####################################################################################

shopt -s expand_aliases extglob progcomp # inevitable for defining wrapper aliases, to source this file which makes use of extglob and use bash completions

export TMPDIR="${TMPDIR:-/tmp}" && mkdir -p "$TMPDIR" || return 1

function _bashbone_setpath(){
	export BASHBONE_PATH=""
	if [[ -e "$BASHBONE_DIR/scripts" ]]; then
		BASHBONE_PATH+="$(realpath -s "$BASHBONE_DIR/scripts" | xargs -echo | sed 's/ /:/g'):"
	fi
	if [[ "$BASHBONE_EXTENSIONDIR" != "$BASHBONE_DIR" && -e "$BASHBONE_EXTENSIONDIR/scripts" ]]; then
		BASHBONE_PATH+="$BASHBONE_EXTENSIONDIR/scripts:"
	fi
	# since tools are shipped compressed, this does not work. ie. gnu parallel needs to be present on target system e.g due to bashbone setup
	# if [[ -e "$BASHBONE_DIR/tools" ]]; then
	# 	BASHBONE_PATH+="$(realpath -s "$BASHBONE_DIR"/tools/*/bin | xargs -echo | sed 's/ /:/g'):"
	# fi
	if [[ -e "$BASHBONE_TOOLSDIR/latest" ]]; then
		BASHBONE_PATH+="$(realpath -s "$BASHBONE_TOOLSDIR/latest/"!(java|bashbone) | xargs -echo | sed 's/ /:/g'):"
	fi
	[[ ":$PATH:" == *":$BASHBONE_PATH"* ]] || PATH="$BASHBONE_PATH$PATH"
}

if ${BASHBONE_CONDA:-false}; then
	source "$BASHBONE_TOOLSDIR/conda/bin/activate" bashbone || return 1
	_bashbone_setpath
fi

#####################################################################################

function _bashbone_reset(){
	trap "" INT ERR

	# disabled for now -> see wait based solution in _wrapper
	if [[ $- == *i* ]]; then
		local key
		# close unused FD and try to read <- which is blocking until EOF comes in
		exec {BASHBONE_FDW}<&-
		read key <&$BASHBONE_FDR
		exec {BASHBONE_FDR}<&-

		# when in a script, daemons will be slayed within on_exit, taking care of if being part of pipeline
		# when interactive, being part of a pipeline (i.e. subprocess) and exitcode is 0, never reached, because PROMPT_COMMAND of shell PID not set
		# when interactive, being part of a pipeline and exitcode is not 0, whole pipeline is going to be killed, which may lead to dataloss
		# test::error 2>&1 | cat <- cat gets killed by INT, so no error messages will be displayed
		[[ $1 -eq 0 ]] && env kill -INT -- -$p &> /dev/null
		sleep 0.1
		env kill -TERM -- -$p &> /dev/null
	fi

	[[ -e "$BASHBONE_TMPDIR/trace" ]] && tac "$BASHBONE_TMPDIR/trace" | awk '!a[$0]{a[$0]=1;print}' | tac >&2
	find -L "$BASHBONE_TMPDIR" -type f -name "cleanup.*.sh" -exec bash -c 'tac "$1" | bash || true; rm -f "$1"' bash {} \;
	"${BASHBONE_EXITFUN:-:}" $1
	rm -rf "$BASHBONE_TMPDIR"

	set -m +o pipefail +o errtrace +o functrace
	# trap - INT RETURN ERR EXIT
	# trap bugs
	# 1) set +m; { trap -p; trap 'echo got int >&2' INT; sleep 2; } & sleep 1; kill -INT -- -$$ -> initial trap -p permanently sets SIGIGN to ingore INT
	# -> reported Jan 2023
	# 2) prompt command cannot reset traps in bash v4
	# trap ':' TERM; PROMPT_COMMAND='trap - TERM'; kill $$ -> kills shell!
	# trap ':' INT; PROMPT_COMMAND='trap - INT'; kill -INT $$ -> kills shell!
	# compared to v5, in bash v4:
	# - only INT can be resetted: PROMPT_COMMAND='trap "trap - INT" INT; kill -INT $$'
	# - this script, needs an other signal. INT gets ignored ¯\_(ツ)_/¯
	trap 'trap - INT USR1 RETURN ERR EXIT' USR1
	env kill -USR1 $BASHPID

	source <(IFS=; echo "$BASHBONE_BAK_TRAPS")
	source <(printf '%s\n' "${BASHBONE_BAK_CMD_NOT_FOUND[@]}")
	source <(printf '%s\n' "${BASHBONE_BAK_SHOPT[@]}")
	source <(printf '%s\n' "${BASHBONE_BAK_SET[@]}")
	source <(printf '%s\n' "${BASHBONE_BAK_ALIASES[@]}")
	${BASHBONE_CONDA:-false} || PATH="${PATH/$BASHBONE_PATH/}"
	if compgen -A arrayvar BASHBONE_BAK_PROMPT_CMD > /dev/null; then
		mapfile -t PROMPT_COMMAND < <(printf "%s\n" "${BASHBONE_BAK_PROMPT_CMD[@]}")
	else
		PROMPT_COMMAND="$BASHBONE_BAK_PROMPT_CMD"
	fi

	# remove all non-exported variables
	local v p=$BASHBONE_PGID
	for v in COMMANDER $(compgen -A variable BASHBONE_ | grep -vFx -f <(compgen -A export BASHBONE_)); do
		unset $v
	done

	return 0
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

function _bashbone_on_exit(){
	trap "" INT TERM ERR
	# when stacked i.e. bash -c 'activate; echo fin lvl1; bash -c "activate; echo fin lvl2";' kill at lvl1 only
	# check also if any stdin/out/err fd is still connected, to not kill a pipe i.e. script_with_bashbone.sh | wc -l . wc runs under same pgid unless script is re-executed via setsid
	# && -t 0 && -t 1 && -t 2 does not work when re-directed into file and when script_with_bashbone.sh is used in a wrapper script like a job script for a workload manager (which may not use tty at all)

	if [[ $1 -eq 0 ]] && ! [[ -t 1 && -t 2 ]]; then
		if [[ ! $BASHBONE_SETSID ]]; then
			commander::warn "Bashbone detected a pipeline and will not sent termination signals to PGID $BASHBONE_PGID" >&2
			commander::warn "Consider to run your script under an exclusive PGID via activate.sh re-execution option '-r true' to avoid orphan processes" >&2
		fi
	else
		if [[ $$ -eq $BASHBONE_PGID ]]; then
			env kill -INT -- -$BASHBONE_PGID
			sleep 0.1
			env kill -TERM -- -$BASHBONE_PGID
		fi
	fi
	_bashbone_reset $1
	eval "$(trap -p EXIT | sed 's/^trap -- .//;s/. EXIT$//')"
	return 0
}

function _bashbone_on_error(){
	# ignore PIPE e.e.g thrown by samtools view | head
	[[ $1 -eq 141 ]] && return 0

	trap "" INT
	if [[ $- == *i* ]]; then
		_bashbone_trace_interactive $1 $2 "$BASHBONE_TMPDIR/trace" &> /dev/null
		mkdir "$BASHBONE_TMPDIR/.killed" &> /dev/null && env kill -INT -- -$BASHBONE_PGID
	else
		mkdir "$BASHBONE_TMPDIR/.killed" &> /dev/null && _bashbone_trace $1 && env kill -INT -- -$BASHBONE_PGID
	fi

	return 130
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

function _bashbone_trace_interactive(){
	local error=$1 line=$2 l fun f src frame last=$((${#BASH_LINENO[@]}-2)) o="${3:-/dev/stderr}" cmd sum add=0
	[[ $BASH_VERSINFO -eq 4 ]] && add=1 && ((++line))
	for frame in $(seq 1 $last); do
		fun=${FUNCNAME[$((frame+1))]}
		[[ "$fun" == "_bashbone_wrapper" || "$fun" == "command_not_found_handle" ]] && continue
		src="${BASH_SOURCE[$((frame+1))]}"

		# in first iteration unless function is command_not_found_handle, line should be set to correct LINENO passed via ERR trap
		[[ $frame -gt 1 || $((line+add)) -eq 1 ]] && line=$((${BASH_LINENO[$frame]}+add))
		[[ $line -eq 1 ]] && continue

		if [[ $fun && "$fun" != "main" ]]; then
			read -r f l src < <(declare -F "$fun")
			l=$((line+l-1))
		fi
		if [[ $line -eq 1 ]]; then
			echo ":ERROR: ${BASHBONE_ERROR:+$BASHBONE_ERROR }in ${src:-shell} (function: ${fun:-main} @ line $l)" >> "$o"
		else
			cmd=$(awk -v l=$l '{ if(NR>=l){if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{print o$0; exit}}else{if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{o=""}}}' "$BASHBONE_WORKINGDIR/$src" | sed -E -e 's/\s+/ /g' -e 's/(^\s+|\s+$)//g')
			echo ":ERROR: ${BASHBONE_ERROR:+$BASHBONE_ERROR }in ${src:-shell} (function: ${fun:-main}) @ line $l: $cmd" >> "$o"
		fi
	done

	return 0
}

#####################################################################################

function _bashbone_wrapper(){
	local BASHBONE_FUNCNAME=$1
	shift

	if ${BASHBONE_SET_ENV:-true}; then
		mapfile -t BASHBONE_BAK_SHOPT < <(shopt | awk '$2=="off"{print "shopt -u "$1}'; shopt | awk '$2=="on"{print "shopt -s "$1}')
		BASHBONE_BAK_TRAPS=$(trap -p)
		mapfile -t BASHBONE_BAK_SET < <(printf "%s" $- | sed 's/[isc]//g' | sed -E 's/(.)/set -\1\n/g')
		mapfile -t BASHBONE_BAK_ALIASES < <(declare -p BASH_ALIASES | sed 's/^declare/declare -g/') # squeeze in global paramater otherwise _bashbone_reset function call declares BASH_ALIASES locally
		mapfile -t BASHBONE_BAK_CMD_NOT_FOUND < <(declare -f command_not_found_handle)
		if compgen -A arrayvar PROMPT_COMMAND > /dev/null; then
			mapfile -t BASHBONE_BAK_PROMPT_CMD < <(printf "%s\n" "${PROMPT_COMMAND[@]}")
			# use to reset temporary env changes when interactive
			# in bash v4 there is a bug where PROMPT_COMMAND is triggered under BASHPID of an async process substitution (test::error @ cat <(cat noexist))
			PROMPT_COMMAND=('[[ $BASHPID -eq $$ ]] || exit; _bashbone_reset $?',"${PROMPT_COMMAND[@]}")
		else
			BASHBONE_BAK_PROMPT_CMD="$PROMPT_COMMAND"
			PROMPT_COMMAND='[[ $BASHPID -eq $$ ]] || exit; _bashbone_reset $?; '"$PROMPT_COMMAND"
		fi

		# +m turns off job control i.e. subshells will not run under own pgid anymore and thus do not ignore int signal (see also setsid)
		set +m
		# traces enable trap to know local scope of functions and shubshells that inherit ERR trap (-E|-o errtrace) and RETURN and DEBUG trap (-T|-o functrace)
		set -o pipefail -o errtrace -o functrace
		# extdebug required by declare -F to get src and line. note that this option implicitly set -E -T. therefore shopt -u extdebug nukes those $SHELLOPTS
		shopt -s expand_aliases extdebug extglob globstar
		[[ $BASH_VERSINFO -gt 4 ]] && shopt -u localvar_inherit
		ulimit -n $(ulimit -Hn)

		BASHBONE_PID=$BASHPID
		BASHBONE_PGID=$(($(ps -o pgid= -p $BASHPID)))
		BASHBONE_TMPDIR="$(command mktemp -d -p "${TMPDIR:-/tmp}" bashbone.XXXXXXXXXX)"

		_bashbone_setpath

		function command_not_found_handle(){
			echo "$1: command not found" >&2
			# necessary to capture contructs like if failcmd; then ... !!!
			# drawback when interactive: incomplete error trace e.g. fun(){ cat <(failcmd); }; cat <(fun)
			(exit 127)
			return 127
		}

		# override mktemp for ease usability of BASHBONE_CLEANUP script per bashbone function upon return
		function mktemp(){
			local tmp="$(command mktemp "$@")"
			echo "rm -rf '$tmp'" >> "${BASHBONE_CLEANUP:-/dev/null}"
			echo "$tmp"
			return 0
		}

		if [[ $- == *i* ]]; then
			# hack to let reset, triggered in parent shell only, to wait for termination of all started async. processes including double forks
			# open fifo read-writable, writebale only and read only. FDs will be inherited by tasks.
			# upon reading from FD in reset function, it will block until last task terminated to receive EOF.
			# -> see also on_reset
			# wait only waits for all direct child processes, but this may be sufficient, because less chances to get into deadlock where also a child waits for parent
			# -> see bash v4 bug below, where async. process can trigger PROMPT_COMMAND execution so parent function and PROMPT_COMMAND both running on_reset waiting for file descriptors
			mkfifo "$BASHBONE_TMPDIR/lock"
			exec {BASHBONE_FD}<>"$BASHBONE_TMPDIR/lock"
			exec {BASHBONE_FDR}<"$BASHBONE_TMPDIR/lock"
			exec {BASHBONE_FD}<&-
			exec {BASHBONE_FDW}>"$BASHBONE_TMPDIR/lock"
			trap '_bashbone_on_error 130 $LINENO; _bashbone_reset 130; return 130 &> /dev/null' INT
			# trap 'wait; _bashbone_on_error 130 $LINENO; _bashbone_reset 130; return 130 &> /dev/null' INT
			trap '_bashbone_on_error $? $LINENO || { [[ $BASHPID -ne $BASHBONE_PID ]] && { return 130 &> /dev/null || exit 130; } || return 130 &> /dev/null; }' ERR
		else
			trap '_bashbone_on_error 130 $LINENO; exit 130 &> /dev/null' INT
			trap '_bashbone_on_error $? $LINENO || { [[ $BASHPID -ne $BASHBONE_PGID ]] && { return 130 &> /dev/null || exit 130; } || exit 130 &> /dev/null; }' ERR
			trap '_bashbone_on_exit $?' EXIT
		fi
		trap '_bashbone_on_return $? "${BASH_COMMAND%% *}"' RETURN
	else
		# necessary to be implemented in order to trigger ERR->TERM and not EXIT directly in case of hidden exit bit thrown e.g. by sleep
		# INT is further triggered in subshells to return and cleanup under PGID in case of e.g. fun(){ cat <(error); echo end; }
		[[ $- == *i* ]] && trap '(exit 130)' INT
	fi
	BASHBONE_SET_ENV=false
	# reset COMMANDER and execute function specific/local $BASHBONE_CLEANUP script

	if [[ $BASHBONE_FUNCNAME ]]; then
		local f_bashbone_wrapper l_bashbone_wrapper s_bashbone_wrapper BASHBONE_LEGACY=$BASHBONE_LEGACY
		read -r f_bashbone_wrapper l_bashbone_wrapper s_bashbone_wrapper < <(declare -F $BASHBONE_FUNCNAME) || true
		if [[ "$s_bashbone_wrapper" == "$BASHBONE_DIR/lib/"* ]]; then
			[[ "$BASHBONE_FUNCNAME" == commander::* || "$BASHBONE_FUNCNAME" == progress::* ]] || BASHBONE_LEGACY=true
		# else
		# 	BASHBONE_LEGACY=false
		fi
		local BASHBONE_BPID=$BASHPID
		local BASHBONE_CLEANUP="$(command mktemp -p "$BASHBONE_TMPDIR" cleanup.XXXXXXXXXX.sh)"
		$BASHBONE_FUNCNAME "$@"
	else
		# allow user to use bashbone cleanup ie auto-cleanup upon mktemp usage and to append their cleanup commands to $BASHBONE_CLEANUP
		# further, users can override _on_exit() when defining qsub/runcmd jobs
		BASHBONE_CLEANUP="$(command mktemp -p "$BASHBONE_TMPDIR" cleanup.XXXXXXXXXX.sh)"
	fi
	return $?
}
[[ $- == *i* ]] || _bashbone_wrapper

#####################################################################################

# ensure expand_aliases shopt is enabled and BASHBONE_DIR to be absolute so is $f and $BASH_SOURCE
# further, ensure to unset all bashbone functions to not wrap nested functions like _usage which leads to unexpected error trace when interactive
# also unset RETURN which tries to call already unset _on_return triggered by source

declare -x -a BASHBONE_FUNCNAMES
mapfile -t BASHBONE_FUNCNAMES < <(
	shopt -s extdebug
	trap - RETURN
	while read -r f; do
		unset -f $f
	done < <(compgen -A function)
	while read -r f; do
		source "$f"
	done < <(find -L "$BASHBONE_DIR/lib/" "$BASHBONE_EXTENSIONDIR/lib/" -name "*.sh" -not -name "#*")
	compgen -A function
)
for f in "${BASHBONE_FUNCNAMES[@]}"; do
	alias $f="_bashbone_wrapper $f"
done

while read -r f; do
	source "$f"
done < <(find -L "$BASHBONE_DIR/lib/" "$BASHBONE_EXTENSIONDIR/lib/" -name "*.sh" -not -name "#*")

#####################################################################################

function bashbone(){
	function _usage(){
		cat <<- EOF
			Bashbone ${BASHBONE_VERSION:+v$BASHBONE_VERSION}

			Is a bash/biobash library for workflow and pipeline design within, but not restricted to, the scope of Next Generation Sequencing (NGS) data analyses.

			Usage:
			-h        | this help
			-l        | enable or disable legacy mode
			-r        | open readme
			-c        | activate bashbone conda environment or deactivate conda
			-s        | list bashbone scripts
			-f <name> | list all functions or if <name> is supplied, execute it
			-d        | list bashbone developer functions
			-e        | list installed tools and versions
			-x        | exit bashbone i.e. remove functions and disable conda. within scripts also restores environment i.e. traps and shell options
		EOF
		return 0
	}

	local OPTIND arg f l s v e i
	while getopts 'hlrcsfdex' arg; do
		case $arg in
		h)	_usage
			;;
		l)	${BASHBONE_LEGACY:-false} && BASHBONE_LEGACY=false || BASHBONE_LEGACY=true
			;;
		r)	if mdless -h &> /dev/null; then
				# or mdless -P <file> | less
				mdless "$BASHBONE_DIR/README.md"
			else
				less "$BASHBONE_DIR/README.md"
			fi
			;;
		c)	if ${BASHBONE_CONDA:-false} && [[ $CONDA_PREFIX ]]; then
				while [[ $CONDA_PREFIX ]]; do
					conda deactivate &> /dev/null
				done
				BASHBONE_CONDA=false
				PATH="${PATH/$BASHBONE_PATH/}"
			else
				source "$BASHBONE_TOOLSDIR/conda/bin/activate" bashbone &> /dev/null || {
					echo ":ERROR: no bashbone installation found" >&2
					return 1
				}
				BASHBONE_CONDA=true
				_bashbone_setpath
			fi
			;;
		s)	find -L "$BASHBONE_DIR/scripts/" "$BASHBONE_EXTENSIONDIR/scripts/" -type f -not -name "test.sh" -not -name "setup.sh" -not -name "#*" | rev | sort -t '/' -u | rev
			;;
		f)	shift $((OPTIND-1));
			if [[ $1 ]]; then
				"$@"
			else
				printf '%s\n' "${BASHBONE_FUNCNAMES[@]}" | grep -vF -e ::_ -e test:: -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V
			fi
			;;
		d)	printf '%s\n' "${BASHBONE_FUNCNAMES[@]}" | grep -F -e helper:: -e progress:: -e commander:: -e configure:: | sort -t ':' -k1,1 -k3,3V
			;;
		e)	(	source "$BASHBONE_TOOLSDIR/conda/bin/activate" bashbone &> /dev/null || {
					echo ":ERROR: no bashbone installation found" >&2
					exit 1
				}
				mapfile -t mapdata < <({ readlink -e "$BASHBONE_TOOLSDIR/latest/"* | sed -nE 's@.*\/([^/]+)-([0-9][^/]+)\/*.*@\L\1\t\2@p;'; mamba list | grep -vP '^(#|_|lib|perl-|xorg-|r-(?!(base|wgcna))|r\s|python-|font|gcc_|gxx_|gfortran_|ca-certificates|pkg-config|pthread)' | tr -s ' ' '\t'; } | sort -k1,1 -k2,2Vr | cut -f 1,2)
				for e in $(mamba env list | grep -F "$BASHBONE_TOOLSDIR" | grep -v '^base' | cut -f 1 -d ' '); do
					mamba list -n $e | grep -vP '^(#|_|lib|perl-|xorg-|r-(?!(base|wgcna))|r\s|python-|font|gcc_|gxx_|gfortran_|ca-certificates|pkg-config|pthread)'
				done | tr -s ' ' '\t' | sort -k1,1 -k2,2Vr | cut -f 1,2 | rev | uniq -f 1 | rev | grep -v -Fw -f <(printf "%s\n" "${mapdata[@]}" | cut -f 1) | sort -k1,1 - <(printf "%s\n" "${mapdata[@]}")
			) || return 1
			;;
		x)	while ${BASHBONE_CONDA:-false} && [[ $CONDA_PREFIX ]]; do
				conda deactivate &> /dev/null
			done
			PATH="${PATH/$BASHBONE_PATH/}"

			if [[ $- == *i* ]]; then
				if [[ ${BASH_VERSINFO[0]} -ge 5 && -r /usr/share/bash-completion/bash_completion ]]; then
					complete -r bashbone samtools bcftools bedtools picard conda mamba datamash bam
					complete -r -I
					_bashbone_completion_default
					if compgen -A arrayvar PROMPT_COMMAND > /dev/null; then
						l=${#PROMPT_COMMAND[@]}
						for i in "${!PROMPT_COMMAND[@]}"; do
							[[ ${PROMPT_COMMAND[$i]} == "_bashbone_completion_default" ]] || PROMPT_COMMAND+=("${PROMPT_COMMAND[$i]}")
						done
						PROMPT_COMMAND=("${PROMPT_COMMAND[@]:$l}")
					else
						PROMPT_COMMAND="${PROMPT_COMMAND/$'\n'_bashbone_completion_default/}"
					fi
				fi
			else
				# already done in prompt_command done when interactive
				_bashbone_reset 0
			fi

			for f in "${BASHBONE_FUNCNAMES[@]}" bashbone $(compgen -A function _bashbone); do
				unset -f $f
				unalias $f &> /dev/null || true
			done

			for v in $(compgen -A variable BASHBONE_); do
				unset $v
			done
			;;
		*)	_usage
			return 1
			;;
		esac
		return 0
	done

	_usage
	return 0
}

# activate colon aware _InitialWorD_ completion
if [[ $- == *i* && ${BASH_VERSINFO[0]} -ge 5 && -r /usr/share/bash-completion/bash_completion ]]; then
	[[ $BASH_COMPLETION_VERSINFO ]] || source /usr/share/bash-completion/bash_completion
	function _bashbone_completion(){
		declare -a args
		mapfile -d ' ' -t args < <(sed 's/\s*$//'< <(printf '%s' "$COMP_LINE"))

		local cur # cur prev words cword
		_get_comp_words_by_ref -n : cur	# _init_completion -n :
		[[ "$3" == "-f" ]] || [[ ${#args[@]} -ge 2 && ${args[-2]} == "-f" && $3 && $cur ]] || return
		COMPREPLY=($(compgen -W "${BASHBONE_FUNCNAMES[*]}" -- "$cur"))
		__ltrim_colon_completions "$cur"

		return 0
	}

	# override to ensure COMP_WORDBREAKS reset in loaded completion functions like _scp, _make or _git ..
	function _completion_loader(){
	    local cmd="${1:-_EmptycmD_}"
	    __load_completion "$cmd" || complete -F _minimal -- "$cmd" && {
	    	local fun=$(complete -p $cmd | grep -oE ' -F [^[:space:]]+' | awk '{print $NF}')
	    	source <(
	    		echo "$fun (){"
	    		echo 'COMP_WORDBREAKS="${COMP_WORDBREAKS/:/}:"'
	    		declare -f $fun | tail -n +3
	    	)
	    	return 124
	    }
	}

	function _bashbone_completion_init(){
		COMP_WORDBREAKS=${COMP_WORDBREAKS/:/}
		# necessary detour to reset in current readline
		complete -F _bashbone_completion_default -D
	}

	function _bashbone_completion_default(){
		COMP_WORDBREAKS="${COMP_WORDBREAKS/:/}:"
		complete -F _completion_loader -D
	}

	if [[ ! ${PROMPT_COMMAND[*]} == *_bashbone_completion_default* ]]; then
		if compgen -A arrayvar PROMPT_COMMAND > /dev/null; then
			PROMPT_COMMAND+=("_bashbone_completion_default")
		else
			PROMPT_COMMAND="${PROMPT_COMMAND:+$PROMPT_COMMAND$'\n'}_bashbone_completion_default"
		fi
	fi

	complete -F _bashbone_completion_init -c -W "${BASHBONE_FUNCNAMES[*]}" -I
	complete -F _bashbone_completion bashbone

	_bashbone_completion_samtools(){
		local cur prev words cword
		_init_completion
		if [[ $cword -eq 1 && -x "$(which samtools)" ]]; then
			COMPREPLY=($(compgen -W "$(samtools --help | sed -nE 's/^\s{5,5}(\S+).*/\1/p' | xargs echo)" -- "$cur"))
		fi
	}

	_bashbone_completion_bcftools(){
		local cur prev words cword
		_init_completion
		if [[ $cword -eq 1 && -x "$(which bcftools)" ]]; then
			COMPREPLY=($(compgen -W "$(bcftools --help | sed -nE '1,/Plugins/{s/^\s{4,4}(\S+).*/\1/p}' | xargs echo)" -- "$cur"))
		fi
	}

	_bashbone_completion_bedtools(){
		local cur prev words cword
		_init_completion
		if [[ $cword -eq 1 && -x "$(which bedtools)" ]]; then
			COMPREPLY=($(compgen -W "$(bedtools --help | sed -nE '/--/!{s/^\s{4,4}(\S+).*/\1/p}' | xargs echo)" -- "$cur"))
		fi
	}

	_bashbone_completion_picard(){
		local cur prev words cword
		_init_completion
		if [[ $cword -eq 1 && -x "$(which picard)" ]]; then
			COMPREPLY=($(compgen -W "$({ picard || true; } |& sed -nE 's/[[:cntrl:]]\[[0-9]{1,3}m//g;s/^\s{4,4}(\S+).*/\1/p' | xargs echo)" -- "$cur"))
		fi
	}

	_bashbone_completion_conda(){
		local cur prev words cword
		_init_completion
		if [[ $cword -eq 1 && -x "$(which conda)" ]]; then
			COMPREPLY=($(compgen -W "$(conda --help | sed -nE 's/^\s{4,4}(\S+).*/\1/p' | xargs echo)" -- "$cur"))
		fi
	}

	_bashbone_completion_mamba(){
		local cur prev words cword
		_init_completion
		if [[ $cword -eq 1 && -x "$(which mamba)" ]]; then
			COMPREPLY=($(compgen -W "$(mamba --help | sed -nE 's/^\s{4,4}(\S+).*/\1/p' | xargs echo)" -- "$cur"))
		fi
	}

	_bashbone_completion_datamash(){
		local cur prev words cword
		_init_completion
		if [[ $cword -eq 1 && -x "$(which datamash)" ]]; then
			COMPREPLY=($(compgen -W "$(datamash --help | sed -nE '1,/Options/{/^\s{2,2}/{s/,//gp}}' | xargs echo)" -- "$cur"))
		fi
	}

	_bashbone_completion_bam(){
		local cur prev words cword
		_init_completion
		if [[ $cword -eq 1 && -x "$(which bam)" ]]; then
			COMPREPLY=($(compgen -W "$(bam help | sed -nE 's/^\s(\S+).*/\1/p' | xargs echo)" -- "$cur"))
		fi
	}

	complete -o default -F _bashbone_completion_samtools samtools
	complete -o default -F _bashbone_completion_bcftools bcftools
	complete -o default -F _bashbone_completion_bedtools bedtools
	complete -o default -F _bashbone_completion_picard picard
	complete -o default -F _bashbone_completion_conda conda
	complete -o default -F _bashbone_completion_mamba mamba
	complete -o default -F _bashbone_completion_datamash datamash
	complete -o default -F _bashbone_completion_bam bam
fi
