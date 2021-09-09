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
			:)	_usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	[[ $exitfun ]] && {
		shift $((OPTIND-(OPTIND-4))) # -p P -f F
		$exitfun "$@"
	}

	declare -a pids=($(pstree -p $pid | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+" | tail -n +2))
	{ kill -INT "${pids[@]}" && wait "${pids[@]}"; } &> /dev/null || true # includes pids of pstree parser pipeline above, thus throws errors necessary to be catched
	# KILL is harsh and does not work for commander::runcmd abortion of parallel jobs
	# INT does not send NOHUB signal
	# TERM is graceful signal and likewise to INT is often cached by applications

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
		local line
		read -r fun line src < <(declare -F "$fun") # requires shopt -s extdebug
		[[ $- =~ i ]] && ((lineno+=line)) # do not!! use [[ ${BASH_EXECUTION_STRING} ]] || ((lineno+=line))
	}
	if [[ -e "$src" ]]; then # self sourcing by commander::runcmd and commander::qsubcmd with cleanup function may causes src to be removed before reaching this point
		local cmd=$(cd "$wdir"; awk -v l=$lineno '{ if(NR>=l){if($0~/\s\\\s*$/){o=o$0}else{print o$0; exit}}else{if($0~/\s\\\s*$/){o=o$0}else{o=""}}}' $src | sed -E -e 's/\s+/ /g' -e 's/(^\s+|\s+$)//g')
		[[ $fun ]] && src="$src ($fun)"
		commander::printerr "${error:-"..an unexpected one"} (exit $ex) @ $src @ line $lineno @ $cmd"
		commander::printerr "$(eval "echo -e \"$cmd\"")"
	else
		[[ $fun ]] && src="$src ($fun)"
		commander::printerr "${error:-"..an unexpected one"} (exit $ex) @ $src @ line $lineno"
	fi

	return 0
}

configure::instances_by_threads(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-i <instances> | number of targeted
			-t <threads>   | per instance targeted
			-T <threads>   | available
		EOF
		return 1
	}

	local OPTIND arg instances ithreads=1 maxthreads
	while getopts 'i:t:T:' arg; do
		case $arg in
			i) instances=$((OPTARG==0?1:OPTARG));;
			t) ithreads=$OPTARG;;
			T) maxthreads=$OPTARG;;
			*) _usage;;
		esac
	done

	local t=$(grep -cF processor /proc/cpuinfo)
	[[ ! $maxthreads ]] && maxthreads=$t
	[[ $maxthreads -gt $t ]] && maxthreads=$t
	[[ ! $instances ]] && instances=$maxthreads

	local maxinstances=$maxthreads
	[[ $maxinstances -gt $(( (maxthreads+10)/ithreads==0?1:(maxthreads+10)/ithreads )) ]] && maxinstances=$(( (maxthreads+10)/ithreads==0?1:(maxthreads+10)/ithreads )) #+10 for better approximation
	[[ $instances -gt $maxthreads ]] && instances=$maxthreads
	[[ $instances -gt $maxinstances ]] && instances=$maxinstances
	ithreads=$((maxthreads/instances))

	echo "$instances $ithreads"
	return 0
}

configure::memory_by_instances(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-i <instances> | number of targeted
			-T <threads>   | available
			-M <memory>    | available
		EOF
		return 1
	}

	local OPTIND arg mandatory instances maxthreads maxmemory
	while getopts 'i:T:M:' arg; do
		case $arg in
			i) ((++mandatory)); instances=$((OPTARG==0?1:OPTARG));;
			T) maxthreads=$OPTARG;;
			M) maxmemory=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage


	local m
	if [[ ! $maxmemory ]]; then
		m=$(grep -F -i memavailable /proc/meminfo | awk '{printf("%d",$2*0.95/1024)}')
		maxmemory=$m
	fi
	m=$(grep -F -i memtotal /proc/meminfo | awk '{printf("%d",$2*0.95/1024)}')
	[[ $maxmemory -gt $m ]] && maxmemory=$m
	local t=$(grep -cF processor /proc/cpuinfo)
	[[ ! $maxthreads ]] && maxthreads=$t
	[[ $maxthreads -gt $t ]] && maxthreads=$t

	[[ $instances -gt $maxthreads ]] && instances=$maxthreads
	local imemory=$((maxmemory/instances))

	echo "$instances $imemory"
	return 0
}

configure::instances_by_memory(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-T <threads> | available
			-m <memory>  | per instance targeted
			-M <memory>  | available
		EOF
		return 1
	}

	local OPTIND arg mandatory maxthreads memory maxmemory
	while getopts 'T:m:M:' arg; do
		case $arg in
			T) maxthreads=$OPTARG;;
			m) ((++mandatory)); memory=$OPTARG;;
			M) maxmemory=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	local m
	if [[ ! $maxmemory ]]; then
		m=$(grep -F -i memavailable /proc/meminfo | awk '{printf("%d",$2*0.95/1024)}')
		maxmemory=$m
	fi
	m=$(grep -F -i memtotal /proc/meminfo | awk '{printf("%d",$2*0.95/1024)}')
	[[ $maxmemory -gt $m ]] && maxmemory=$m
	local t=$(grep -cF processor /proc/cpuinfo)
	[[ ! $maxthreads ]] && maxthreads=$t
	[[ $maxthreads -gt $t ]] && maxthreads=$t

	[[ $memory -gt $maxmemory ]] && memory=$maxmemory
	local instances=$((maxmemory/memory))
	[[ $instances -gt $maxthreads ]] && instances=$maxthreads
	local ithreads=$((maxthreads/instances))

	echo "$instances $ithreads"
	return 0
}

configure::jvm(){
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-i <instances> | number of targeted
			-t <threads>   | per instance targeted
			-T <threads>   | available
			-m <memory>    | per instance targeted
			-M <memory>    | available
		EOF
		return 1
	}

	local OPTIND arg instances ithreads=1 maxthreads memory=1 maxmemory
	while getopts 'i:t:T:m:M:' arg; do
		case $arg in
			i)	instances=$((OPTARG==0?1:OPTARG));;
			t)	ithreads=$OPTARG;;
			T)	maxthreads=$OPTARG;;
			m)	memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			*)	_usage
		esac
	done

	local m
	if [[ ! $maxmemory ]]; then
		m=$(grep -F -i memavailable /proc/meminfo | awk '{printf("%d",$2*0.95/1024)}')
		maxmemory=$m
	fi
	m=$(grep -F -i memtotal /proc/meminfo | awk '{printf("%d",$2*0.95/1024)}')
	[[ $maxmemory -gt $m ]] && maxmemory=$m
	local t=$(grep -cF processor /proc/cpuinfo)
	[[ ! $maxthreads ]] && maxthreads=$t
	[[ $maxthreads -gt $t ]] && maxthreads=$t
	[[ ! $instances ]] && instances=$maxthreads

	local jmem jgct jcgct
	local maxinstances=$((maxmemory/memory))
	[[ $maxinstances -gt $(( (maxthreads+10)/ithreads==0?1:(maxthreads+10)/ithreads )) ]] && maxinstances=$(( (maxthreads+10)/ithreads==0?1:(maxthreads+10)/ithreads )) #+10 for better approximation
	[[ ! $instances ]] && instances=$maxinstances
	[[ $instances -gt $maxthreads ]] && instances=$maxthreads
	[[ $instances -gt $maxinstances ]] && instances=$maxinstances

	ithreads=$((maxthreads/instances))
	jmem=$((maxmemory/instances))
	jgct=$(((3+5*ithreads/8)/instances)) # For N <= 8 parallel GC will use just as many, i.e., N GC threads. For N > 8 available processors, the number of GC threads will be computed as 3+5N/8
	[[ $jgct -eq 0 ]] && jgct=1
	jcgct=$((jgct/4)) # Sets n to approximately 1/4 of the number of parallel garbage collection threads (ParallelGCThreads).
	[[ $jcgct -eq 0 ]] && jcgct=1

	echo "$instances $ithreads $jmem $jgct $jcgct"
	return 0
}
