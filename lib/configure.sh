#!/usr/bin/env bash
# (c) Konstantin Riege

configure::exit(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-p <pid>      | process id
			-f <function> | to call
		EOF
		return 1
	}

	local OPTIND arg mandatory pid fun
	while getopts 'p:f:' arg; do
		case $arg in
			p) ((++mandatory)); pid=$OPTARG;;
			f) fun=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	shift $((OPTIND-1))
	$(declare -F $fun &> /dev/null) && $fun "$@"

	sleep 1 # to get very last entry of logifle by tail -f before being killed
	declare -a pids=($(pstree -p $pid | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+" | tail -n +2))
	{ kill -KILL "${pids[@]}" && wait "${pids[@]}"; } &> /dev/null
	printf "\r"

	return 0
}

configure::err(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-x <exit>     | better use as first option. code to print and return
			-f <function> | name to print
			-s <source>   | filename of
			-l <lineno>   | LINENO
			-e <error>    | message
			-c <cmd>      | command

		EOF
		return 1
	}

	local OPTIND arg mandatory fun src lineno error cmd ex
	while getopts 'f:s:l:e:c:x:' arg; do
		case $arg in
			f) fun="$OPTARG";;
			s) src="$(basename "$OPTARG")";;
			l) ((++mandatory)); lineno=$OPTARG;;
			e) error="$OPTARG";;
			c) cmd="$OPTARG";;
			x) ((++mandatory)); ex=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	[[ $fun ]] && {
		local line
		shopt -s extdebug
		read -r fun line src < <(declare -F "$fun")
		shopt -u extdebug
		src="$(basename "$src")"
		src="$fun / $src"
		[[ -t 1 ]] && ((lineno+=line))
	}

	unset BASH_COMMAND
	commander::printerr "${error:-"..an unexpected one"} (exit $ex) @ $src (line $lineno) $cmd"
	return $ex
}

configure::environment(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-i <insdir> | root path to tools
			-b <insdir> | root path to bashbone
			-c <conda>  | true/false activate
		EOF
		return 1
	}

	local OPTIND arg mandatory insdir_tools insdir_bashbone activate_conda=true
	while getopts 'i:b:c:' arg; do
		case $arg in
			i) ((++mandatory)); insdir_tools="$OPTARG";;
			b) ((++mandatory)); insdir_bashbone="$OPTARG";;
			c) activate_conda="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	set -o pipefail
	shopt -s extglob
	shopt -s expand_aliases
	ulimit -n $(ulimit -Hn)
	export MALLOC_ARENA_MAX=4

	$activate_conda && {
		commander::printinfo "setting up environment"
		source $insdir_tools/conda/bin/activate base &> /dev/null
	}

	local tp=$(readlink -e $insdir_tools/latest)
	# better stay off custom java path to avoid conflicts with conda openjdk and pre-compiled jars requesting specific versions (IncompatibleClassChangeError)
	# [[ $tp && -e $tp/java ]] && export JAVA_HOME=$(dirname $(readlink -e $tp/java))
	[[ $tp ]] && export PATH=$(readlink -e $tp/!(java) | xargs -echo | sed 's/ /:/g'):$PATH
	export PATH=$(readlink -e $insdir_bashbone/scripts | xargs -echo | sed 's/ /:/g'):$PATH

	return 0
}

configure::instances_by_threads(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-i <instances> | number of all
			-t <threads>   | per instance targeted
			-T <threads>   | available
		EOF
		return 1
	}

	local OPTIND arg mandatory instances ithreads=1 maxthreads
	while getopts 'i:t:T:m:' arg; do
		case $arg in
			i) ((++mandatory)); instances=$OPTARG;;
			t) ithreads=$OPTARG;;
			T) ((++mandatory)); maxthreads=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	local maxinstances=$maxthreads
	[[ $maxinstances -gt $(( (maxthreads+10)/ithreads==0?1:(maxthreads+10)/ithreads )) ]] && maxinstances=$(( (maxthreads+10)/ithreads==0?1:(maxthreads+10)/ithreads )) #+10 for better approximation
	[[ $instances -gt $maxthreads ]] && instances=$maxthreads
	[[ $instances -gt $maxinstances ]] && instances=$maxinstances
	ithreads=$((maxthreads/instances))

	echo "$instances $ithreads"
	return 0
}

configure::instances_by_memory(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-t <threads> | available
			-m <memory>  | per instance maximum
		EOF
		return 1
	}

	local OPTIND arg mandatory threads memory
	while getopts 't:m:' arg; do
		case $arg in
			t) ((++mandatory)); threads=$OPTARG;;
			m) ((++mandatory)); memory=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	local maxmemory=$(grep -F -i memavailable /proc/meminfo | awk '{printf("%d",$2*0.9/1024)}')
	[[ $memory -gt $maxmemory ]] && memory=$maxmemory
	local instances=$((maxmemory/memory))
	[[ $instances -gt $threads ]] && instances=$threads
	local ithreads=$((threads/instances))

	echo "$instances $ithreads"
	return 0
}

configure::jvm(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-i <instances> | number of all
			-t <threads>   | per instance targeted
			-T <threads>   | available
			-m <memory>    | per instance maximum
		EOF
		return 1
	}

	local OPTIND arg mandatory instances ithreads=1 maxthreads memory=1
	while getopts 'i:t:T:m:' arg; do
		case $arg in
			i)	instances=$OPTARG;;
			t)	ithreads=$OPTARG;;
			T)	((++mandatory)); maxthreads=$OPTARG;;
			m)	memory=$OPTARG;;
			*)	_usage
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage
	[[ ! $instances ]] && instances=$maxthreads

	local jmem jgct jcgct maxmemory=$(grep -F -i memavailable /proc/meminfo | awk '{printf("%d",$2*0.9/1024)}')
	local maxinstances=$((maxmemory/memory))
	[[ $maxinstances -gt $(( (maxthreads+10)/ithreads==0?1:(maxthreads+10)/ithreads )) ]] && maxinstances=$(( (maxthreads+10)/ithreads==0?1:(maxthreads+10)/ithreads )) #+10 for better approximation
	[[ $instances -gt $maxthreads ]] && instances=$maxthreads
	[[ $instances -gt $maxinstances ]] && instances=$maxinstances
	ithreads=$((maxthreads/instances))

	jmem=$((maxmemory/instances))
	[[ $memory -gt 1 ]] && [[ $jmem -gt $memory ]] && jmem=$memory
	jgct=$(((3+5*ithreads/8)/instances))
	[[ $jgct -eq 0 ]] && jgct=1
	jcgct=$((jgct/4))
	[[ $jcgct -eq 0 ]] && jcgct=1

	echo "$instances $ithreads $jmem $jgct $jcgct"
	return 0
}
