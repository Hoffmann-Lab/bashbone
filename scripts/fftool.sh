#!/usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$0")/../bashbone_lite.sh" -a "$@" || exit 1

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename $0) flat file index or random access

		VERSION
		0.2.0

		SYNOPSIS
		$(basename $0) [-i] [-f <file>] [-o <file>]

		OPTIONS
		-h           | this help
		-i           | index input file (see -f) or from stdin (requires -o)
		-s <size>    | index record size. default 1000
		-f <infile>  | path to file to be indexed (see -i) or accessed (requires -f). default: stdin (requires -o)
		-o <outfile> | if input comes from stdin, path to output file (mutual exclusive to -f)
		-r <range>   | for random data access. format: <#|inf>L@<#>L i.e. get a certain number of lines, all until EOF by inf keyword respectively, after skipping a certain number of lines.
		             | example: 5L@20L means get 5 lines after skipping 20 lines

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	return 1
}

while getopts f:o:r:s:ih ARG; do
	case $ARG in
		f)	f="$OPTARG";;
		r)	range="$OPTARG";;
		i)	i=true;;
		o)	o="$OPTARG"; mkdir -p "$(dirname "$o")";;
		s)	recordsize=$OPTARG;;
		h)	{ usage || exit 0; };;
		*)	usage;
	esac
done
[[ $# -eq 0 ]] && { usage || exit 0; }
if [[ ! $f && ! $o ]]; then
	usage
fi
recordsize=${recordsize:-1000}

if ${i:-false}; then
	if [[ $f ]]; then
		# xsv index $f
		awk -v l=$recordsize '{c=c+length($0)+1}NR%l==0{print NR"\t"c}END{print NR"\t"c}' "$f" > "${f%.*}.ffi"
	else
		tee -i >(awk -v l=$recordsize '{c=c+length($0)+1}NR%l==0{print NR"\t"c}END{print NR"\t"c}' > "${o%.*}.ffi") > "$o" | cat
	fi
else
	if [[ ! $r ]]; then
		usage
	fi
	recordsize=$(head -1 "${f%.*}.ffi" | cut -f 1)
	read -r n skip < <(sed -E 's/L@?/ /g' <<< "$range")
	# if [[ -e "${f%.*}.ffi" ]]; then
	exec {FD}<>"$f"
	perl -slane 'exit if $F[0]>$skip; $offset=$F[1]; END{open(my $fh, "<&=", $fd); seek($fh,$offset,0); if($skip%$rs!=0){ while(<$fh>){last if ++$i >= $skip%$rs}}}' -- -fd=$FD -skip=$skip -rs=$recordsize "${f%.*}.ffi"
	if [[ "$n" == "inf" ]]; then
		# xsv slice -n -s $skip $f
		cat <&$FD
	else
		# xsv slice -n -s $skip -l $n $f
		head -n $n <&$FD
	fi
	exec {FD}>&-
	# else
	# 	if [[ "$n" == "inf" ]]; then
	# 		tail -n +$((skip+1)) "$f"
	# 	else
	# 		# xsv slice -n -s $skip -l $n $f
	# 		tail -n +$((skip+1)) "$f" | head -n $n
	# 	fi
	# fi
fi
