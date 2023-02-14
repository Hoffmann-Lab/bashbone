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
	declare -a pids
	local tmpdir
	_cleanup::progress::log(){
		rm -rf "$tmpdir"
		[[ $pids ]] && { env kill -PIPE ${pids[@]} && wait ${pids[@]}; } &> /dev/null || true
	}

	_usage(){
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
			f)	((++mandatory)); fun=$OPTARG; shift $((OPTIND-1));;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	tmpdir="$(mktemp -d -p "/tmp" fifo.XXXXXXXXXX)"
	mkfifo "$tmpdir/stderr" "$tmpdir/stdout"
	case $verbosity in
		0)	{ progress::_bar & } 2>/dev/null # do not use subshell here. will not terminated
			pids+=($!)
			# use subshells to avoid job control messages
			( tee -ia "$log" < "$tmpdir/stdout" | { grep -E --line-buffered '^\s*(:INFO:|:BENCHMARK:|:WARNING:)' || true; } & )
			pids+=($!)
			( tee -ia "$log" < "$tmpdir/stderr" | { grep -E --line-buffered '^\s*:ERROR:' >&2 || true; } & )
			pids+=($!)
		;;
		1)	{ progress::_bar & } 2>/dev/null
			pids+=($!)
			( tee -ia "$log" < "$tmpdir/stdout" | { grep -E --line-buffered '^\s*(:INFO:|:CMD:|:BENCHMARK:|:WARNING:)' || true; } & )
			pids+=($!)
			( tee -ia "$log" < "$tmpdir/stderr" | { grep -E --line-buffered '^\s*:ERROR:' >&2 || true; } & )
			pids+=($!)
		;;
		2)	( tee -ia "$log" < "$tmpdir/stdout" & )
			pids+=($!)
			( tee -ia "$log" < "$tmpdir/stderr" >&2 & )
			pids+=($!)
		;;
		*)	_usage;;
	esac

	$fun "$@" 2> "$tmpdir/stderr" > "$tmpdir/stdout"

	return 0
}
