#! /usr/bin/env bash
# (c) Konstantin Riege

[[ "${BASH_SOURCE[0]}" == "$0" ]] && echo "$(basename "${BASH_SOURCE[0]}") script needs to be sourced" >&2 && exit 1
[[ ! "$(ps -p $$ -o command= | cut -d ' ' -f 1)" =~ bash ]] && echo "requires a bash shell" >&2 && return 1
[[ "$(uname)" != "Linux" ]] && echo "unsupported operating system" >&2 && return 1
[[ ${BASH_VERSINFO[0]} -lt 4 || (${BASH_VERSINFO[0]} -eq 4 && ${BASH_VERSINFO[1]} -lt 4) ]] && echo "requieres bash >= v4.4" >&2 && return 1

#####################################################################################

if [[ ! $- =~ i ]] && ${BASHBONE_SETSID:-true}; then
	# re-run under own pgid via setsid, to be able to kill process group when kill signal is send to pid only e.g. by SGE/PBE qdel or xargs
	# do not set sid when executed with gnu parallel!

	# exec bash -c 'export PGID=$$; trap "trap \"\" INT TERM; env kill -INT -- -\$PGID" INT; trap "trap \"\" INT TERM; env kill -TERM -- -\$PGID" TERM; setsid --wait env --default-signal=INT,QUIT bash "$0" "$@" & PGID=$!; wait $PGID' "$(realpath -s "$0")" "$@"
	export BASHBONE_SETSID=false
	{ trap 'exit 130' INT; exec setsid --wait bash "$(realpath -s "$0")" "$@"; } &
	BASHBONE_PGID=$!
	# always trap INT to trigger ERR that triggers TERM to implement termination sequence with INT coming first to avoid most termination messages
	# termination sequence: INT,1000,TERM,0 plus error tracing and cleanup
	trap 'trap "" INT TERM; env kill -INT -- -$BASHBONE_PGID' INT TERM
	wait $BASHBONE_PGID # would be 130 upon INT and not exit code of script if script defines trap 'exit 42' INT
	wait $BASHBONE_PGID # would be 42 or again 130 unless script defines trap 'exit 42' INT
	trap - EXIT # in case user defined EXIT trap, which must not be called here again
	exit $?
fi

#####################################################################################

export BASHBONE_DIR="$(dirname "$(readlink -e "${BASH_SOURCE[0]}")")"
export BASHBONE_VERSION=$(source "$BASHBONE_DIR/lib/version.sh"; echo $version)
export BASHBONE_TOOLSDIR="${BASHBONE_TOOLSDIR:-$(dirname "$BASHBONE_DIR")}" # export to not override if previously defined
shopt -s expand_aliases extglob # inevitable for defining wrapper aliases and to source this file which makes use of extglob

unset OPTIND
while getopts 'i:p:c:x:ah' arg; do
	case $arg in
		i)	BASHBONE_TOOLSDIR="$OPTARG";;
		p)	TMPDIR="$OPTARG";;
		c)	BASHBONE_CONDA="$OPTARG";;
		x)	BASHBONE_EXITFUN="$OPTARG";;
		a)	shift $((OPTIND-1)); break;;
		h)	cat <<- 'EOF'
				This is bashbone activation script.

				To see lists of available options and functions, source me and execute bashbone -h

				Usage:
				-h            | this help
				-i <path>     | to installation root <path>/latest/<tools>/<bins>
				                default: inferred from script location
				                hint: run activation from source code, indeed enables basic functions, but will fail on executing tools
				-p <path>     | to temporary directory. default: $TMPDIR or /tmp
				-c <activate> | true/false conda from [-i]/conda/bin
				                default: false
				-x <string>   | a function or command string to be called upon EXIT signal. this function or command will receive the exit code as last positional argument
				-a <optarg>   | use as last option in case bashbone is sourced in a script wich makes use of optargs

				Example:
				source activate.sh -i <path> -c true -x "<fun> [<arg>..]" -a "$@"
				bashbone -h
			EOF
			return 0
		;;
	esac
done

export TMPDIR="${TMPDIR:-/tmp}" && mkdir -p "$TMPDIR" || return 1

if ${BASHBONE_CONDA:-false}; then
	source "$BASHBONE_TOOLSDIR/conda/bin/activate" bashbone || return 1
fi

#####################################################################################

function _bashbone_reset(){
	# for sake of interactive execution, stay away from touching/defining TERM !!!
	# trap ':' TERM; PROMPT_COMMAND='trap - TERM'; kill $$ -> kills shell! TERM can not be resetted via PROMPT_COMMAND nor USR traps and propably more. trap -p lies!
	# see also { trap -p; ..} & bug where trap -p permanently sets SIGIGN to ingore INT
	trap "" INT ERR
	# when in PROMT_COMMAND: in case cleanup was already triggered by ERR, use atomic mkdir which will block if lock dir exists or if tmpdir already removed
	mkdir "$BASHBONE_TMPDIR/.lock4" &> /dev/null && {
		# echo ":INFO: cleanup $BASHBONE_TMPDIR in progress" >&2
		find "$BASHBONE_TMPDIR" -type f -name "cleanup.*.sh" -exec bash -c 'bash "$1"; rm -f "$1"' bash {} \;
		rm -rf "$BASHBONE_TMPDIR"
	}

	# get rid of progress::log
	{ env kill -INT -- -$BASHBONE_PGID & wait $!; } &> /dev/null

	set -m +o pipefail +o errtrace +o functrace
	trap - RETURN ERR INT EXIT

	source <(IFS=; echo "$BASHBONE_BAK_TRAPS")
	source <(printf '%s\n' "${BASHBONE_BAK_CMD_NOT_FOUND[@]}")
	source <(printf '%s\n' "${BASHBONE_BAK_SHOPT[@]}")
	source <(printf '%s\n' "${BASHBONE_BAK_SET[@]}")
	source <(printf '%s\n' "${BASHBONE_BAK_ALIASES[@]}")
	PATH="${PATH/$BASHBONE_PATH:/}"
	PROMPT_COMMAND=("${BASHBONE_BAK_PROMPT_CMD[@]}")

	# remove all non-exported variables
	local v
	for v in COMMANDER $(declare -p | grep -E '^declare -[^x] BASHBONE_' | cut -d ' ' -f 3 | cut -d '=' -f 1); do
		unset $v
	done
	return 0
}

function _bashbone_on_return(){
	if [[ $1 -eq 0 && "$2" != "source" && "$2" != "." && "$BASHBONE_FUNCNAME" == "${FUNCNAME[1]}" ]]; then
		COMMANDER=()
		if [[ -e "$BASHBONE_CLEANUP" ]]; then
			bash "$BASHBONE_CLEANUP" || true
			rm -f "$BASHBONE_CLEANUP"
		fi
	fi
	return 0
}

function _bashbone_on_exit(){
	trap "" INT TERM ERR
	{ env kill -INT -- -$BASHBONE_PGID & wait $!; } &> /dev/null
	sleep 0.1
	{ env kill -TERM -- -$BASHBONE_PGID & wait $!; } &> /dev/null
	_bashbone_reset
	"${BASHBONE_EXITFUN:-:}" $1
	return 0
}

function _bashbone_on_err(){
	[[ $1 -eq 141 ]] && return 0

	trap "" INT
	if [[ $1 -eq 130 ]]; then
		_bashbone_trace(){ :; }
		_bashbone_trace_interactive(){ :; }
	fi
	if [[ $- =~ i ]]; then
		# only recent subshell (if not further stacked in command_not_found subshell) and PGID shell can trace without wrong line numbers
		# process substitutions and async subshells need LINENO passthough vie ERR trap for proper tracing
		# use atomic mkdir to block other subshells error trace upon killing process group and upon return function stack within current shell
		if [[ $BASHPID -ne $BASHBONE_PGID ]]; then
			# _bashbone_trace_interactive may sets $BASHBONE_TMPDIR/.lock2
			# further use INT trap to trigger error->return->cleanup via under PGID (see wrapper) in case of fun() { cat <(error); echo end; }
			mkdir "$BASHBONE_TMPDIR/.lock1" &> /dev/null && _bashbone_trace_interactive $1 $2 && env kill -INT $BASHBONE_PGID
		else
			mkdir "$BASHBONE_TMPDIR/.lock2" &> /dev/null && _bashbone_trace_interactive $1 $2
		fi
	else
		[[ $1 -eq 143 ]] && _bashbone_trace(){ :; }
		mkdir "$BASHBONE_TMPDIR/.lock1" &> /dev/null && _bashbone_trace $1
	fi

	# { env kill -TERM -- -$BASHBONE_PGID & wait $!; } &> /dev/null
	# better spawn suicide job with new PGID, which also cleans up
	# this allows calling functions interactively via process substitution, where PROMPT_COMMAND is not set to perform cleanup and reset under PGID e.g. cat <(fun)
	# _bashbone_reset uses $BASHBONE_TMPDIR/.lock4 in case error happens, wehn interactive, cleanup starts and then prompt_command tries to cleanup again
	mkdir "$BASHBONE_TMPDIR/.lock3" &> /dev/null && {
		set -m
		{ env kill -TERM -- -$BASHBONE_PGID; _bashbone_reset; } &
		set +m
		wait $!
	}
	return 143
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
	echo ":ERROR: exit code 143" >&2;
	return 0
}

function _bashbone_trace_interactive(){
	local error=$1 line=$2 l fun f src frame last=$((${#BASH_LINENO[@]}-2)) cmd sum
	for frame in $(seq 1 $last); do
		fun=${FUNCNAME[$((frame+1))]}
		[[ "$fun" == "_bashbone_wrapper" || "$fun" == "command_not_found_handle" ]] && continue
		src="${BASH_SOURCE[$((frame+1))]}"
		# in first iteration unless function is command_not_found_handle, line may be set to proper LINENO passed via ERR trap
		# if $line > 1 && $line == {BASH_LINENO[1]} : error & or { error; } & or (error) & plus wait $!
		# if $line > 1 : error or { error; } or (error)
		# if $line == 1 && {BASH_LINENO[1]} > 1 : cat <(error)
		[[ $frame -eq 1 && $line -eq ${BASH_LINENO[$frame]} ]] && sum=false || sum=true
		[[ $frame -gt 1 || $line -eq 1 ]] && line=${BASH_LINENO[$frame]}
		# rest is done under PGID
		[[ $line -eq 1 ]] && return 0

		if [[ $fun && "$fun" != "main" ]]; then
			read -r f l src < <(declare -F "$fun")
			$sum && ((line+=l-1))
		fi
		cmd=$(awk -v l=$line '{ if(NR>=l){if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{print o$0; exit}}else{if($0~/\s\\\s*$/){o=o""gensub(/\\\s*$/,"",1,$0)}else{o=""}}}' "$src" | sed -E -e 's/\s+/ /g' -e 's/(^\s+|\s+$)//g')
		echo ":ERROR:$BASHBONE_ERROR in ${src:-shell} (function: ${fun:-main}) @ line $line: $cmd" >&2
	done
	echo ":ERROR: exit code 143" >&2
	# no further need to do anything under PGID
	mkdir "$BASHBONE_TMPDIR/.lock2" &> /dev/null || true
	return 0
}

#####################################################################################

# override mktemp for ease usability of BASHBONE_CLEANUP script per bashbone function upon return
function mktemp(){
	local tmp="$(command mktemp "$@")"
	echo "rm -rf '$tmp'" >> "${BASHBONE_CLEANUP:-/dev/null}"
	echo "$tmp"
	return 0
}

function _bashbone_wrapper(){
	local BASHBONE_FUNCNAME=$1
	shift

	if ${BASHBONE_SET_ENV:-true}; then
		mapfile -t BASHBONE_BAK_SHOPT < <(shopt | awk '$2=="off"{print "shopt -u "$1}'; shopt | awk '$2=="on"{print "shopt -s "$1}')
		BASHBONE_BAK_TRAPS=$(trap -p)
		mapfile -t BASHBONE_BAK_SET < <(printf "%s" $- | sed 's/[is]//g' | sed -E 's/(.)/set -\1\n/g')
		mapfile -t BASHBONE_BAK_ALIASES < <(declare -p BASH_ALIASES | sed 's/^declare/declare -g/') # squeeze in global paramater otherwise _bashbone_reset function call declares BASH_ALIASES locally
		mapfile -t BASHBONE_BAK_CMD_NOT_FOUND < <(declare -f command_not_found_handle)
		mapfile -t BASHBONE_BAK_PROMPT_CMD < <(printf "%s\n" "${PROMPT_COMMAND[@]}")

		# +m turns off job control i.e. subshells will not run under own pgid anymore and thus do not ignore int signal (see also setsid)
		set +m
		# traces enable trap to know local scope of functions and shubshells that inherit ERR trap (-E|-o errtrace) and RETURN and DEBUG trap (-T|-o functrace)
		set -o pipefail -o errtrace -o functrace
		# extdebug required by declare -F to get src and line. note that this option implicitly set -E -T. therefore shopt -u extdebug nukes those $SHELLOPTS
		shopt -s expand_aliases extdebug extglob globstar
		[[ $BASH_VERSINFO -gt 4 ]] && shopt -u localvar_inherit
		ulimit -n $(ulimit -Hn)

		# use to reset temporary env changes when interactive
		PROMPT_COMMAND="_bashbone_reset"

		BASHBONE_PGID=$(($(ps -o pgid= -p $BASHPID)))
		BASHBONE_TMPDIR="$(command mktemp -d -p "${TMPDIR:-/tmp}" bashbone.XXXXXXXXXX)"

		if [[ -e "$BASHBONE_TOOLSDIR/latest" ]]; then
			BASHBONE_PATH="$(realpath -s "$BASHBONE_TOOLSDIR/latest/"!(java|bashbone) | xargs -echo | sed 's/ /:/g')"
		fi
		if [[ -e "$BASHBONE_DIR/scripts" ]]; then
			BASHBONE_PATH="$(realpath -s "$BASHBONE_DIR/scripts" | xargs -echo | sed 's/ /:/g')${BASHBONE_PATH:+":$BASHBONE_PATH"}"
		fi
		if [[ -e "$BASHBONE_DIR/tools" ]]; then
			BASHBONE_PATH="$(realpath -s "$BASHBONE_DIR"/tools/*/bin | xargs -echo | sed 's/ /:/g')${BASHBONE_PATH:+":$BASHBONE_PATH"}"
		fi
		PATH="${BASHBONE_PATH:+"$BASHBONE_PATH:"}$PATH"

		function command_not_found_handle(){
			echo "$1: command not found" >&2
			# necessary to capture contructs like if failcmd; then ... !!!
			# drawback when interactive: incomplete error trace e.g. fun(){ cat <(failcmd); }; cat <(fun)
			(exit 127)
			return 127
		}
	fi
	BASHBONE_SET_ENV=false

	# reset COMMANDER and execute function specific/local $BASHBONE_CLEANUP script
	trap '_bashbone_on_return $? "${BASH_COMMAND%% *}"' RETURN
	# pass though LINENO which is more often correct than ${BASH_LINENO[0]} when interactive
	trap '_bashbone_on_err $? $LINENO || return 143' ERR
	# necessary to be implemented in order to trigger ERR->TERM and not EXIT directly in case of hidden exit bit thrown e.g. by sleep
	# INT is further triggered in subshells to return and cleanup under PGID in case of e.g. fun(){ cat <(error); echo end; }
	trap '(exit 130)' INT
	# just in case defined by user ...
	trap - TERM

	local BASHBONE_CLEANUP="$(command mktemp -p "$BASHBONE_TMPDIR" cleanup.XXXXXXXXXX.sh)"
	$BASHBONE_FUNCNAME "$@"
	return 0
}

# to capture non-function error
[[ $- =~ i ]] || {
	_bashbone_wrapper true
	trap '_bashbone_on_exit $?' EXIT
}

#####################################################################################

# ensure expand_aliases shopt is enabled and BASHBONE_DIR to be absolute so is $f and $BASH_SOURCE
while read -r f; do
	alias $f="_bashbone_wrapper $f"
done < <(
	shopt -s extdebug
	for f in "$BASHBONE_DIR/lib/"*.sh; do
		source $f
	done
	for f in $(declare -F | awk '{print $3}'); do
		read -r f l s < <(declare -F $f)
		[[ "$s" == "$BASHBONE_DIR/lib/"* ]] && echo $f
	done
)
for f in "$BASHBONE_DIR/lib/"*.sh; do
	source $f
done

#####################################################################################

function bashbone(){
	function _usage(){
		cat <<- EOF
			Bashbone ${BASHBONE_VERSION:+v$BASHBONE_VERSION}

			Is a bash/biobash library for workflow and pipeline design within, but not restricted to, the scope of Next Generation Sequencing (NGS) data analyses.

			Usage:
			-h | this help
			-r | open readme
			-c | activate bashbone conda environment or deactivate current conda environment
			-p | list installed pipelines that utilize bashbone
			-s | list bashbone scripts
			-f | list bashbone function
			-d | list bashbone developer functions
			-e | list installed tools and versions
			-x | remove bashbone functions and within scripts, restore environment (traps and shell options)
		EOF
		return 0
	}

	local OPTIND arg f l s
	while getopts 'hrcpsfdex' arg; do
		case $arg in
		h)	_usage
			;;
		r)	if mdless -h &> /dev/null; then
				# or mdless -P <file> | less
				mdless "$BASHBONE_DIR/README.md"
			else
				less "$BASHBONE_DIR/README.md"
			fi
			;;
		c)	if [[ $CONDA_PREFIX ]]; then
				while [[ $CONDA_PREFIX ]]; do
					conda deactivate &> /dev/null
				done
			else
				source "$BASHBONE_TOOLSDIR/conda/bin/activate" bashbone &> /dev/null || {
					echo ":ERROR: no bashbone installation found" >&2
					return 1
				}
			fi
			;;
		p)	[[ -e "$BASHBONE_TOOLSDIR/latest" ]] && find -L "$BASHBONE_TOOLSDIR/latest" -maxdepth 2 -name "*.sh" -not -name "activate.sh" -not -name "setup.sh" -printf "%f\n" || {
				echo ":ERROR: no bashbone installation found" >&2
				return 1
			}
			;;
		s)	find "$BASHBONE_DIR/scripts/" -type f -printf "%f\n" | rev | sort | rev
			;;
		f)	(	shopt -s extdebug
				for f in $(declare -F | awk '{print $3}'); do
					read -r f l s < <(declare -F $f)
					[[ "$s" == "$BASHBONE_DIR/lib/"* ]] && echo $f
				done | grep -vF -e ::_ -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V
			)
			;;
		d)	(	shopt -s extdebug
				for f in $(declare -F | awk '{print $3}'); do
					read -r f l s < <(declare -F $f)
					[[ "$s" == "$BASHBONE_DIR/lib/"* ]] && echo $f
				done | grep -F -e ::_ -e compile:: -e helper:: -e progress:: -e commander:: -e configure:: -e options:: | sort -t ':' -k1,1 -k3,3V
			)
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
		x)	[[ $- =~ i ]] || _bashbone_reset # already done when interactive
			# do in subshell to not alter enviroment again by extdebug
			while read -r f; do
				unset -f $f
				unalias $f &> /dev/null || true
			done < <(
				shopt -s extdebug
				for f in $(declare -F | awk '{print $3}'); do
					read -r f l s < <(declare -F $f)
					[[ "$s" == "$BASHBONE_DIR/lib/"* || "$f" == *"bashbone"* ]] && echo $f
				done
			)
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
