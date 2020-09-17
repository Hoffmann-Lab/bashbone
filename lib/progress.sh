#!/usr/bin/env bash
# (c) Konstantin Riege

progress::_bar() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

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
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	_usage(){
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-v [0|1|2] | verbosity level
			-o <file>  | path to
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

	# do not grep :ERROR: since this goes to stderr and will be printed anyways
	case $verbosity in
		0)	progress::_bar &
			{ tail -f $log | grep -E --line-buffered '^\s*(:INFO:|:BENCHMARK:|:WARNING:)'; } &
			;;
		1)	progress::_bar &
			{ tail -f $log | grep -E --line-buffered '^\s*(:INFO:|:CMD:|:BENCHMARK:|:WARNING:)'; } &
			;;
		2)	{ tail -f $log | grep -v -E --line-buffered '^\s*:ERROR:'; } &;;
		*)	_usage;	return 1;;
	esac

	return 0
}
