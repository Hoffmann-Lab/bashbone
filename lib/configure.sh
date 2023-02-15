#!/usr/bin/env bash
# (c) Konstantin Riege

function configure::instances_by_threads(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-i <instances> | number of targeted
			-t <threads>   | per instance targeted
			-T <threads>   | available
		EOF
		return 1
	}

	local OPTIND arg instances ithreads=1 maxthreads dryrun=false
	while getopts 'i:t:T:d' arg; do
		case $arg in
			i) instances=$((OPTARG==0?1:OPTARG));;
			t) ithreads=$OPTARG;;
			T) maxthreads=$OPTARG;;
			d) dryrun=true;;
			*) _usage;;
		esac
	done

	local t=$($dryrun && [[ $maxthreads ]] && echo $maxthreads || grep -cF processor /proc/cpuinfo)
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

function configure::memory_by_instances(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-i <instances> | number of targeted
			-T <threads>   | available
			-M <memory>    | available
		EOF
		return 1
	}

	local OPTIND arg mandatory instances maxthreads maxmemory dryrun=false
	while getopts 'i:T:M:d' arg; do
		case $arg in
			i) ((++mandatory)); instances=$((OPTARG==0?1:OPTARG));;
			T) maxthreads=$OPTARG;;
			M) maxmemory=$OPTARG;;
			d) dryrun=true;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	local m=$($dryrun && [[ $maxmemory ]] && echo $maxmemory || grep -F MemAvailable /proc/meminfo | awk '{printf("%d",$2*0.9/1024)}')
	[[ ! $maxmemory ]] && maxmemory=$m
	[[ $maxmemory -gt $m ]] && maxmemory=$m
	local t=$($dryrun && [[ $maxthreads ]] && echo $maxthreads || grep -cF processor /proc/cpuinfo)
	[[ ! $maxthreads ]] && maxthreads=$t
	[[ $maxthreads -gt $t ]] && maxthreads=$t

	[[ $instances -gt $maxthreads ]] && instances=$maxthreads
	local imemory=$((maxmemory/instances))

	echo "$instances $imemory"
	return 0
}

function configure::instances_by_memory(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-T <threads> | available
			-m <memory>  | per instance targeted
			-M <memory>  | available
		EOF
		return 1
	}

	local OPTIND arg mandatory maxthreads memory maxmemory dryrun=false
	while getopts 'T:m:M:d' arg; do
		case $arg in
			T) maxthreads=$OPTARG;;
			m) ((++mandatory)); memory=$OPTARG;;
			M) maxmemory=$OPTARG;;
			d) dryrun=true;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 1 ]] && _usage

	local m=$($dryrun && [[ $maxmemory ]] && echo $maxmemory || grep -F MemAvailable /proc/meminfo | awk '{printf("%d",$2*0.9/1024)}')
	[[ ! $maxmemory ]] && maxmemory=$m
	[[ $maxmemory -gt $m ]] && maxmemory=$m
	local t=$($dryrun && [[ $maxthreads ]] && echo $maxthreads || grep -cF processor /proc/cpuinfo)
	[[ ! $maxthreads ]] && maxthreads=$t
	[[ $maxthreads -gt $t ]] && maxthreads=$t

	[[ $memory -gt $maxmemory ]] && memory=$maxmemory
	local instances=$((maxmemory/memory))
	[[ $instances -gt $maxthreads ]] && instances=$maxthreads
	local ithreads=$((maxthreads/instances))

	echo "$instances $ithreads"
	return 0
}

function configure::jvm(){
	function _usage(){
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

	local OPTIND arg instances ithreads=1 maxthreads memory=1 maxmemory dryrun=false
	while getopts 'i:t:T:m:M:d' arg; do
		case $arg in
			i)	instances=$((OPTARG==0?1:OPTARG));;
			t)	ithreads=$OPTARG;;
			T)	maxthreads=$OPTARG;;
			m)	memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			d)	dryrun=true;;
			*)	_usage
		esac
	done

	local m=$($dryrun && [[ $maxmemory ]] && echo $maxmemory || grep -F MemAvailable /proc/meminfo | awk '{printf("%d",$2*0.9/1024)}')
	[[ ! $maxmemory ]] && maxmemory=$m
	[[ $maxmemory -gt $m ]] && maxmemory=$m
	local t=$($dryrun && [[ $maxthreads ]] && echo $maxthreads || grep -cF processor /proc/cpuinfo)
	[[ ! $maxthreads ]] && maxthreads=$t
	[[ $maxthreads -gt $t ]] && maxthreads=$t
	[[ ! $instances ]] && instances=$maxthreads

	local jmem jgct jcgct
	[[ $memory -gt $maxmemory ]] && memory=$maxmemory
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
