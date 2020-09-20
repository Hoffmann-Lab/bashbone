#!/usr/bin/env bash
# (c) Konstantin Riege

_configure::test(){
	_cleanup::_configure::test(){
		echo "${FUNCNAME[0]} of ${FUNCNAME[1]}"
	}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			test usage
		EOF
	}
	echo "${FUNCNAME[0]}"
	failcmd -a -b
	echo "this message must not be shown"
}

configure::exit(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-p <pid>      | process id
			-f <function> | to call
		EOF
		return 1
	}

	local OPTIND arg mandatory pid fun
	while getopts 'p:f:' arg; do
		case $arg in
			p)	((++mandatory)); pid=$OPTARG;;
			f)	fun=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	shift $((OPTIND-1))
	declare -F $fun &> /dev/null && $fun "$@"

	sleep 1 # to get very last entry of logifle by tail -f before being killed
	declare -a pids=($(pstree -p $pid | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+" | tail -n +2))
	{ kill -KILL "${pids[@]}" && wait "${pids[@]}"; } &> /dev/null || true # includes pids of pstree parser pipeline above, thus throws errors necessary to be catched
	printf "\r"
	return 0
}

configure::err(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-x <exit>     | better use as first option. code to print and return
			-f <function> | name to print
			-s <source>   | filename of
			-l <lineno>   | LINENO
			-e <error>    | message
		EOF
		return 1
	}

	local OPTIND arg mandatory fun src lineno error cmd ex
	while getopts 'f:s:l:e:x:' arg; do
		case $arg in
			f) fun="$OPTARG";;
			s) src="$(basename "$OPTARG")";;
			l) ((++mandatory)); lineno=$OPTARG;;
			e) error="$OPTARG";;
			x) ((++mandatory)); ex=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	[[ $fun ]] && {
		local line
		read -r fun line src < <(declare -F "$fun") # requires shopt -s extdebug
		[[ $- =~ i ]] && ((lineno+=line))
	}
	local cmd=$(awk -v l=$lineno '{ if(NR>=l){if($0~/\s\\\s*$/){o=o$0}else{print o$0; exit}}else{if($0~/\s\\\s*$/){o=o$0}else{o=""}}}' $src | sed -E -e 's/\s+/ /g' -e 's/(^\s+|\s+$)//g')
	[[ $fun ]] && src="$src ($fun)"

	sleep 1 # to be very last entry of a logifle
	commander::printerr "${error:-"..an unexpected one"} (exit $ex) @ $src @ line $lineno @ $cmd"
	return $ex
}

configure::environment(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-i <insdir> | root path to tools
			-b <insdir> | root path to bashbone
			-c <conda>  | true/false activate
		EOF
		return 1
	}

	local OPTIND arg mandatory insdir_tools activate_conda=true
	while getopts 'i:b:c:' arg; do
		case $arg in
			i) ((++mandatory)); insdir_tools="$OPTARG";;
			c) activate_conda="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	$activate_conda && {
		source $insdir_tools/conda/bin/activate base &> /dev/null && {
			commander::printinfo "utilizing $(conda --version)"
		}|| {
			commander::printerr "conda activation failed. use -i [path] to point towards installation directory"
			return 1
		}
	}

	set -E -o pipefail -o functrace # -E allows simple trap bubbeling and -o functrace enables inheritance of RETURN and DEBUG trap
	shopt -s extdebug # do not shopt -u extdebug, otherwise set -E -o pipefail configuration will be nuked
	shopt -s extglob
	shopt -s expand_aliases
	ulimit -n $(ulimit -Hn)
	export MALLOC_ARENA_MAX=4

	local tp="$(readlink -e "$insdir_tools/latest")"
	# better stay off custom java path to avoid conflicts with conda openjdk and pre-compiled jars requesting specific versions (IncompatibleClassChangeError)
	# [[ $tp && -e $tp/java ]] && export JAVA_HOME=$(dirname $(readlink -e $tp/java))
	[[ $tp ]] && export PATH="$(readlink -e "$tp/!(java)" | xargs -echo | sed 's/ /:/g'):$PATH"
	export PATH="$(readlink -e "$BASHBONE_DIR/scripts" | xargs -echo | sed 's/ /:/g'):$PATH"

	# make use of local scope during trace to trigger tmp file deletion etc.
	trap 'declare -F _cleanup::${FUNCNAME[0]} &> /dev/null && _cleanup::${FUNCNAME[0]}' RETURN
	if [[ $- =~ i ]]; then
		# do not split in muliple lines. lineno will be wrong
		# since trap needs to persist in shell, make sure return is triggerd only from sourced bashbone functions. otherwise there will be issues with bash completion, vte etc.
		trap 'e=$?; if [[ ${BASH_SOURCE[0]} && "$(readlink -e "${BASH_SOURCE[0]}")" =~ "$BASHBONE_DIR" ]]; then configure::err -x $e -e "$BASHBONE_ERROR" -l $LINENO -f ${FUNCNAME[0]}; return $e; fi' ERR
	else
		# dont call exit directly. allow for back trace through all functions. local scopes are available
		trap 'e=$?; if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then configure::err -x $e -e "$BASHBONE_ERROR" -l $LINENO -s "$0"; exit $e; else configure::err -x $e -e "$BASHBONE_ERROR" -l $LINENO -f ${FUNCNAME[0]}; return $e; fi' ERR
	fi

	return 0
}

configure::instances_by_threads(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
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
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
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
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
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
