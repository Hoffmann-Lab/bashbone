#!/usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$0")/../bashbone_lite.sh" -a "$@" || exit 1

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename $0) shuffle and split single or paired fastq(.gz|.bz2) files into two pseudo-replicate fastq(.gz) files

		VERSION
		0.1.0

		SYNOPSIS
		$(basename $0) -1 <fastq> [-2 <fastq>] -o <prefix> [-p <tmpdir> -t <threads> -z]

		OPTIONS
		-h          | this help
		-1 <path>   | input SE or first mate fastq(.gz|.bz2)
		-2 <path>   | input mate pair fastq(.gz|.bz2)
		-z          | compress output using pigz with [-t] threads (fallback: gzip)
		-t <value>  | compression threads (default: $t)
		-o <path>   | prefix of output fastq, extended by [1|2].(R1.|R2.)fastq(.gz)
		-p <path>   | temp directory (default: /tmp)

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	return 1
}

t=$(grep -cF processor /proc/cpuinfo 2> /dev/null || echo 1)
while getopts 1:2:o:t:p:n:zh ARG; do
	case $ARG in
		1) i="$OPTARG";;
		2) j="$OPTARG";;
		o) o="$OPTARG";;
		p) tmpdir="$OPTARG";;
		t) t=$OPTARG;;
		z) z=true;;
		h) { usage || exit 0; };;
		*) usage;
	esac
done
[[ $# -eq 0 ]] && { usage || exit 0; }
if [[ ! $i || ! $o ]]; then
	usage
fi

mkdir -p "$(dirname "$o")"

[[ $tmpdir ]] && {
	mkdir -p "$tmpdir"
	params="-p $tmpdir"
} || params=""
tmp="$(mktemp $params shufnsplit.XXXXXXXXXX)"
#trap 'rm -f "$tmp"*' INT TERM EXIT

${z:-false} && {
	if [[ $j ]]; then
		pigz -h 2> /dev/null && z="pigz -k -c -p $(((t+3)/4))" || z="gzip -k -c"
	else
		pigz -h 2> /dev/null && z="pigz -k -c -p $(((t+1)/2))" || z="gzip -k -c"
	fi
	e=".fastq.gz"
} || {
	z="cat"
	e=".fastq"
}

open=$({ readlink -e "$i" | file -b --mime-type -f - | grep -oF -e 'gzip' -e 'bzip2' || echo cat; } | sed '/cat/!{s/gzip/pigz -p 1/; s/$/ -cd/}')

if [[ $j ]]; then
	paste <($open "$i" | sed -E '/^\s*$/d' | paste - - - -) <($open "$j" | sed -E '/^\s*$/d' | paste - - - -) | shuf > "$tmp"
	split -a 1 --numeric-suffixes=1 --additional-suffix="$e" -n l/2 --filter="bash -c \" awk -F '\\\t' -v OFS='\\\n' '{print \\\$1,\\\$2,\\\$3,\\\$4; print \\\$5,\\\$6,\\\$7,\\\$8 > \\\"/dev/fd/2\\\"}' > >($z > \$FILE.R1) 2> >($z > \$FILE.R2) \"" "$tmp" "$o"
else
	$open "$i" | sed -E '/^\s*$/d' | paste - - - - | shuf > "$tmp"
	split -a 1 --numeric-suffixes=1 --additional-suffix="$e" -n l/2 --filter="tr '\t' '\n' | $z > \$FILE" "$tmp" "$o"
fi
