#!/usr/bin/env bash
# (c) Konstantin Riege

configure::exit(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-p <pid>      | process id
			-f <function> | to call - ALWAYS LAST OPTION
		EOF
		return 1
	}

	local OPTIND arg mandatory pid exitfun
	while getopts ':p:f:' arg; do
		case $arg in
			p)	((++mandatory)); pid=$OPTARG;;
			f)	exitfun=$OPTARG;;
			:) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	[[ $exitfun ]] && {
		shift $((OPTIND-(OPTIND-4))) # -p P -f F
		$exitfun "$@"
	}

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
			-w <workdir>  | path to initially used
		EOF
		return 1
	}

	local OPTIND arg mandatory wdir="$PWD" fun src lineno error cmd ex
	while getopts 'w:f:s:l:e:x:' arg; do
		case $arg in
			w)	wdir="$OPTARG";;
			f)	fun="$OPTARG";;
			s)	src="$OPTARG";;
			l)	((++mandatory)); lineno=$OPTARG;;
			e)	error="$OPTARG";;
			x)	((++mandatory)); ex=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	[[ $fun ]] && {
		declare -F "$fun"
		local line
		read -r fun line src < <(declare -F "$fun") # requires shopt -s extdebug
		[[ $- =~ i ]] && ((lineno+=line))
	}
	local cmd=$(cd "$wdir"; awk -v l=$lineno '{ if(NR>=l){if($0~/\s\\\s*$/){o=o$0}else{print o$0; exit}}else{if($0~/\s\\\s*$/){o=o$0}else{o=""}}}' $src | sed -E -e 's/\s+/ /g' -e 's/(^\s+|\s+$)//g')
	[[ $fun ]] && src="$src ($fun)"

	sleep 1 # to be very last entry of a logifle
	commander::printerr "${error:-"..an unexpected one"} (exit $ex) @ $src @ line $lineno @ $cmd"

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
