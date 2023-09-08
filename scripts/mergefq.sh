#!/usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$0")/../bashbone_lite.sh" -a "$@" || exit 1

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename $0) merges and sorts by read name or unmerges paired fastq(.gz|.bz2) files

		VERSION
		0.1.1

		SYNOPSIS MERGE
		$(basename $0) -[tmd] <value> -i <fastq1> -j <fastq2> -o <fastq>

		SYNOPSIS UNMERGE
		$(basename $0) -u 1 -i <fastq> -o <fastq1>
		$(basename $0) -u 2 -i <fastq> -o <fastq2>

		OPTIONS
		-t [value]  | threads ($t)
		-m [string] | amount of memory to use eg. 64000M (available: $m)
		-d [path]   | tmpdir with size eq -m ($PWD)
		-i [path]   | input fastq_R1 (see -j) or merged fastq (see -u)
		-j [path]   | input fastq_R2 - triggers merge mode
		-o [path]   | out fastq
		-z          | compress output using pigz with $t threads (fallback: gzip)
		-u [1|2]    | unmerge fastq (see -i) to fastq_R1 or fastq_R2 (see -o)
		-h          | this help

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	return 1
}

i=''
j=''
o=''
z=''
h=''
u=0
d="${TMPDIR:-/tmp}"
t=$(grep -cF processor /proc/cpuinfo)
m=$(grep -F MemTotal /proc/meminfo | awk '{printf("%dM",$2/1024*0.95)}')
while getopts i:j:o:u:d:t:m:zh ARG; do
	case $ARG in
		i) i="$OPTARG";;
		j) j="$OPTARG";;
		o) o="$OPTARG"; mkdir -p "$(dirname "$o")";;
		u) u=$OPTARG;;
		d) d="$OPTARG";;
		t) t=$OPTARG;;
		m) m=$OPTARG;;
		z) z=1;;
		h) { usage || exit 0; };;
		*) usage;
	esac
done
[[ $# -eq 0 ]] && { usage || exit 0; }
if [[ ! $i || ! $o || $u -gt 2 ]]; then
	usage
fi

if [[ $z ]]; then
	pigz -h &> /dev/null && z="pigz -k -c -p $t" || z="gzip -k -c"
else
	z='cat'
fi

open=$({ readlink -e "$i" | file -b --mime-type -f - | grep -oF -e 'gzip' -e 'bzip2' || echo cat; } | sed '/cat/!{s/gzip/pigz -p 1/; s/$/ -cd/}')

if [[ $u -gt 0 ]]; then
	if [[ $u -eq 2 ]]; then
		$open "$i" | sed -E '/^\s*$/d' | paste - - - - | sed -n '2~2p' | tr '\t' '\n' | $z > "$o"
	else
		$open "$i" | sed -E '/^\s*$/d' | paste - - - - | sed -n '1~2p' | tr '\t' '\n' | $z > "$o"
	fi
else
	#paste <($open "$i") <($open "$j") | paste - - - - - - - - | awk -F '\t' -v OFS='\n' '{print $1,$3,$5,$7; print $2,$4,$6,$8}' | $z > "$o"
	$open "$i" "$j" | sed -E '/^\s*$/d' | paste - - - - | LC_ALL=C sort -k1,1 -S $m -T "$d" --parallel=$t | tr '\t' '\n' | $z > "$o"
fi
