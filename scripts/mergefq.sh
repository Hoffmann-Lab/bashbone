#!/usr/bin/env bash
# (c) Konstantin Riege

usage(){
cat <<- EOF
	DESCRIPTION
	$(basename $0) merges or unmerges paired fastq(.gz|.bz2) files

	VERSION
	0.1.0

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
	exit 0
}

i=''
j=''
o=''
z=''
h=''
u=0
d=$PWD
t=$(cat /proc/cpuinfo | grep -cF processor 2> /dev/null || echo 1)
m=$(grep -i memtotal /proc/meminfo | awk '{printf("%d",$2*0.8/1024)}' 2> /dev/null || echo "2G")
while getopts i:j:o:u:d:t:m:zh ARG; do
	case $ARG in
		i) i="$OPTARG";;
		j) j="$OPTARG";;
		o) o="$OPTARG";;
		u) u=$OPTARG;;
		d) d="$OPTARG";;
		t) t=$OPTARG;;
		m) m=$OPTARG;;
		z) z=1;;
		h) h=1;;
	esac
done

if [[ $# -eq 0 ]] || [[ $h ]] || [[ ! $i ]] || [[ ! $o ]] || [[ $u -gt 2 ]]; then
	usage
fi

if [[ $z ]]; then
	[[ $(which pigz 2> /dev/null) ]] && z="pigz -k -c -p $t" || z="gzip -k -c"
else
	z='cat'
fi

open=$(readlink -e "$i" | file -f - | grep -Eo '(gzip|bzip)' && echo -cd || echo cat)

if [[ $u -gt 0 ]]; then
	if [[ $u -eq 2 ]]; then
		exec $open "$i" | sed -r '/^\s*$/d' | paste - - - - | sed -n '2~2p' | tr '\t' '\n' | $z > "$o"
	else
		exec $open "$i" | sed -r '/^\s*$/d' | paste - - - - | sed -n '1~2p' | tr '\t' '\n' | $z > "$o"
	fi
else
	exec $open "$i" "$j" | sed -r '/^\s*$/d' | paste - - - - | LC_ALL=C sort -k1,1 -S $m -T "$d" --parallel=$t | tr '\t' '\n' | $z > "$o"
fi

exit
