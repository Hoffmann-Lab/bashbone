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
	local tmpdir
	_cleanup::progress::log(){
		rm -rf "$tmpdir"
	}

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
			o)	((++mandatory)); log="$OPTARG"; mkdir -p "$(dirname "$log")";;
			f)	((++mandatory)); fun=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	tmpdir="$(mktemp -d -p "/dev/shm" fifo.XXXXXXXXXX)"
	mkfifo "$tmpdir/stderr" "$tmpdir/stdout"
	case $verbosity in
		0)	jobs | grep -F 'progress::_bar' || progress::_bar &
			tee -ia "$log" < "$tmpdir/stdout" | grep -E --line-buffered '^\s*(:INFO:|:BENCHMARK:|:WARNING:)' &
			tee -ia "$log" < "$tmpdir/stderr" | grep -E --line-buffered '^\s*:ERROR:' >&2 &
			;;
		1)	jobs | grep -F 'progress::_bar' || progress::_bar &
			tee -ia "$log" < "$tmpdir/stdout" | grep -E --line-buffered '^\s*(:INFO:|:CMD:|:BENCHMARK:|:WARNING:)' &
			tee -ia "$log" < "$tmpdir/stderr" | grep -E --line-buffered '^\s*:ERROR:' >&2 &
		;;
		2)	tee -ia "$log" < "$tmpdir/stdout" &
			tee -ia "$log" < "$tmpdir/stderr" >&2 &
		;;
		*)	_usage;;
	esac

	shift $((OPTIND-1))
	$fun "$@" 2> "$tmpdir/stderr" > "$tmpdir/stdout"

	return 0
}
