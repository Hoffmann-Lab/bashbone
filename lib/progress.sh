#!/usr/bin/env bash
# (c) Konstantin Riege

progress::bar() {
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

progress::log0() {
	tail -f $1 2>&1 | grep -E --line-buffered '^\s*(:INFO:|:ERROR:|:BENCHMARK:|:WARNING:)'

	return 0
}

progress::log1() {
	tail -f $1 2>&1 | grep -E --line-buffered '^\s*(:INFO:|:CMD:|:ERROR:|:BENCHMARK:|:WARNING:)'

	return 0
}
