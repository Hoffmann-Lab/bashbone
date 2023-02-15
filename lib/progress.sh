#!/usr/bin/env bash
# (c) Konstantin Riege

function progress::_bar(){
	trap 'return 0' INT

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

function progress::log(){
	declare -a pids
	local tmpdir
	function _cleanup::progress::log(){
		rm -rf "$tmpdir"
	}

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v [0|1|2]    | verbosity level
			-o <logfile>  | path to
			-f <function> | and parameters - needs to be last!
		EOF
		return 1
	}

	local OPTIND arg mandatory log verbosity fun
	while getopts 'v:o:f:' arg; do
		case $arg in
			v)	((++mandatory)); verbosity=$OPTARG;;
			o)	((++mandatory)); log="$OPTARG"; mkdir -p "$(dirname "$log")";;
			f)	((++mandatory)); fun=$OPTARG; shift $((OPTIND-1)); break;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	tmpdir="$(mktemp -d -p "${TMPDIR:-/tmp}" fifo.XXXXXXXXXX)"
	mkfifo "$tmpdir/stderr" "$tmpdir/stdout"
	case $verbosity in
		0)	{ progress::_bar & } 2>/dev/null
			{ tee -ia "$log" < "$tmpdir/stdout" | { trap '(exit 130)' INT; grep --color=never -E --line-buffered '^\s*(:INFO:|:BENCHMARK:|:WARNING:)' || true; } & } 2> /dev/null
			{ tee -ia "$log" < "$tmpdir/stderr" | { trap '(exit 130)' INT; grep --color=never -E --line-buffered '^\s*:ERROR:' || true; } & } >&2 2> /dev/null
		;;
		1)	{ { trap '(exit 130)' INT; progress::_bar; } & } 2>/dev/null
			{ tee -ia "$log" < "$tmpdir/stdout" | { trap '(exit 130)' INT; grep --color=never -E --line-buffered '^\s*(:INFO:|:CMD:|:BENCHMARK:|:WARNING:)' || true; } & } 2> /dev/null
			{ tee -ia "$log" < "$tmpdir/stderr" | { trap '(exit 130)' INT; grep --color=never -E --line-buffered '^\s*:ERROR:' || true; } & } >&2 2> /dev/null
		;;
		2)	{ tee -ia "$log" < "$tmpdir/stdout" & } 2> /dev/null
			{ tee -ia "$log" < "$tmpdir/stderr" & } >&2 2> /dev/null
		;;
		*)	_usage;;
	esac

	$fun "$@" 2> "$tmpdir/stderr" > "$tmpdir/stdout"

	return 0
}
