#!/usr/bin/env bash
# (c) Konstantin Riege

progress::_bar() {
	local mod=0
	while true; do
		((++mod))
		case $mod in
			[15]) echo -en "\r|";;
			[26]) echo -en "\r/";;
			[37]) echo -en "\r-";;
			4) echo -en "\r\\";;
			8) echo -en "\r\\"; mod=0;;
		esac
		sleep 0.2
	done

	return 0
}

progress::log() {
	local funcname=${FUNCNAME[0]}
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v [0|1|2] | verbosity level
			-o <file>  | path to
		EOF
		return 0
	}
	local OPTIND arg mandatory log verbosity
	while getopts 'v:o:' arg; do
		case $arg in
			v)	((++mandatory)); verbosity=$OPTARG;;
			o)	((++mandatory)); log="$OPTARG";;
			*)	_usage;	return 1;;
		esac
	done
	[[ $mandatory -lt 2 ]] && { _usage; return 1; }

	case $verbosity in
		0)	progress::_bar &
			{ tail -f $log 2>&1 | grep -E --line-buffered '^\s*(:INFO:|:ERROR:|:BENCHMARK:|:WARNING:)'; } &
			;;
		1)	progress::_bar &
			{ tail -f $log 2>&1 | grep -E --line-buffered '^\s*(:INFO:|:CMD:|:ERROR:|:BENCHMARK:|:WARNING:)'; } &
			;;
		2)	{ tail -f $log 2>&1; } &;;
		*)	_usage;	return 1;;
	esac

	return 0
}
