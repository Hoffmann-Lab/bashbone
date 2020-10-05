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
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v [0|1|2]    | verbosity level
			-o <logfile>  | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory log verbosity
	while getopts 'v:o:' arg; do
		case $arg in
			v)	((++mandatory)); verbosity=$OPTARG;;
			o)	((++mandatory)); log="$OPTARG";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	commander::printinfo "$(date)" > "$log"
	# do not grep :ERROR: since this goes to stderr and will be printed anyways
	case $verbosity in
		0)	progress::_bar &
			{ tail -f "$log" | grep -E --line-buffered '^\s*(:INFO:|:BENCHMARK:|:WARNING:)'; } &
			;;
		1)	progress::_bar &
			{ tail -f "$log" | grep -E --line-buffered '^\s*(:INFO:|:CMD:|:BENCHMARK:|:WARNING:)'; } &
			;;
	esac

	return 0
}

progress::observe() {
	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v [0|1|2]    | verbosity level
			-o <logfile>  | path to
			-f <function> | and parameters
		EOF
		return 1
	}

	local OPTIND arg mandatory log verbosity fun
	while getopts 'v:o:f:' arg; do
		case $arg in
			v)	((++mandatory)); verbosity=$OPTARG;;
			o)	((++mandatory)); log="$OPTARG";;
			f)	fun=$OPTARG
				[[ $mandatory -lt 2 ]] && _usage
				shift $((OPTIND-1))
				case $verbosity in
					0)	$fun "$@" 2> >(tee -ai "$log" | grep -E --line-buffered '^\s*:ERROR:' >&2 || true) >> "$log";;
					1)	$fun "$@" 2> >(tee -ai "$log" | grep -E --line-buffered '^\s*:ERROR:' >&2 || true) >> "$log";;
					2)	$fun "$@" 2> >(tee -ai "$log" >&2) | tee -ai "$log";;
					*)	_usage;;
				esac
				return 0
			;;
			*)	_usage;;
		esac
	done

	_usage
}
