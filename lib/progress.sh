#!/usr/bin/env bash
# (c) Konstantin Riege

function progress::_wheel(){
	trap 'return 0' TERM

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

function progress::_log(){
	trap - RETURN

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v [0..2]     | verbosity level
			-o <logfile>  | path to
			-f <fifo>     | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory log verbosity fifo
	while getopts 'v:o:f:' arg; do
		case $arg in
			v)	((++mandatory)); verbosity=$OPTARG;;
			o)	((++mandatory)); log="$OPTARG";;
			f)	((++mandatory)); fifo=$OPTARG;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	case $verbosity in
		0)	tee -ia "$log" < "$fifo" || true;;
		1)	tee -ia "$log" < "$fifo" | grep --color=never -E --line-buffered '^\s*(:INFO:|:CMD:|:BENCHMARK:|:WARNING:)' || true;;
		2)	tee -ia "$log" < "$fifo" | grep --color=never -E --line-buffered '^\s*:ERROR:' || true;;
		*)	_usage;;
	esac

	return 0
}

function progress::log(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-v [0..2]     | verbosity level
			                0 : bashbone :CMD:|:INFO:|:BENCHMARK:|:WARNING:|:ERROR:
			                1 : stderr + bashbone :CMD:|:INFO:|:BENCHMARK:|:WARNING:|:ERROR:
			                2 : full, i.e. stderr + stdout + bashbone :CMD:|:INFO:|:BENCHMARK:|:WARNING:|:ERROR:
			-o <logfile>  | path to
			-r            | override existing logs
			-f <function> | and parameters - needs to be last!
		EOF
		return 1
	}

	local OPTIND arg mandatory log verbosity=2 fun override=false
	while getopts 'v:o:f:r' arg; do
		case $arg in
			v)	verbosity=$OPTARG;;
			o)	((++mandatory)); log="$OPTARG"; mkdir -p "$(dirname "$log")";;
			r)	override=true;;
			f)	((++mandatory)); fun=$OPTARG; shift $((OPTIND-1)); break;;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	local tmpdir="$(mktemp -d -p "${TMPDIR:-/tmp}" fifo.XXXXXXXXXX)"
	mkfifo "$tmpdir/stderr" "$tmpdir/stdout"
	$override && rm -f "$log"

	case $verbosity in
		0)
			{ progress::_wheel & } 2>/dev/null
			echo "kill -TERM $! &> /dev/null" >> "$BASHBONE_CLEANUP"
			set -m
			{ progress::_log -v 1 -o "$log" -f "$tmpdir/stdout" & } 2>/dev/null
			{ progress::_log -v 2 -o "$log" -f "$tmpdir/stderr" & } >&2 2>/dev/null
			set +m
		;;
		1)	{ progress::_wheel & } 2>/dev/null
			echo "kill -TERM $! &> /dev/null" >> "$BASHBONE_CLEANUP"
			set -m
			{ progress::_log -v 1 -o "$log" -f "$tmpdir/stdout" & } 2>/dev/null
			{ progress::_log -v 0 -o "$log" -f "$tmpdir/stderr" & } >&2 2>/dev/null
			set +m
		;;
		*)	set -m
			{ progress::_log -v 0 -o "$log" -f "$tmpdir/stdout" & } 2>/dev/null
			{ progress::_log -v 0 -o "$log" -f "$tmpdir/stderr" & } >&2 2>/dev/null
			set +m
		;;
		*)	_usage;;
	esac

	# manually wrap. aliases are not expanded on runtime
	_bashbone_wrapper $fun "$@" 2> "$tmpdir/stderr" > "$tmpdir/stdout"

	return 0
}
