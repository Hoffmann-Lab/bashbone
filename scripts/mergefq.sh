#!/usr/bin/env bash
# (c) Konstantin Riege
usage() {	
	echo "merge paired fastq(.gz|.bz2) files:"
	echo "$(basename $0) -[tmd] <value> -i <fastq1> -j <fastq2> -o <fastq>"
	echo "unmerge:"
	echo "$(basename $0) -u 1 -i <fastq> -o <fastq1>"
	echo "$(basename $0) -u 2 -i <fastq> -o <fastq2>"
	echo "options:"
	echo "-t [value]  | threads ($t)"
	echo "-m [string] | amount of memory to use eg. 64000M (available: $m)"
	echo "-d [path]   | tmpdir with size eq -m ($PWD)"
	echo "-i [path]   | input fastq_R1 (see -j) or merged fastq (see -u)"
	echo "-j [path]   | input fastq_R2 - triggers merge mode"
	echo "-o [path]   | out fastq"
	echo "-z          | compress output using pigz with $t threads (fallback: gzip)"
	echo "-u [1|2]    | unmerge fastq (see -i) to fastq_R1 or fastq_R2 (see -o)"
	echo "-h          | this help"
	exit 1
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
		exec $open "$i" | paste - - - - | sed -n '2~2p' | tr '\t' '\n' | $z > "$o"
	else 
		exec $open "$i" | paste - - - - | sed -n '1~2p' | tr '\t' '\n' | $z > "$o"
	fi
else
	exec $open "$i" "$j" | paste - - - - | LC_ALL=C sort -k1,1 -S $m -T "$d" --parallel=$t | tr '\t' '\n' | $z > "$o"
fi

exit
