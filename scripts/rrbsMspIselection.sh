#!/usr/bin/env bash
# (c) Konstantin Riege

usage(){
cat <<- EOF
	DESCRIPTION
	$(basename $0) select single or pairs of reads with proper MspI cutting site from fastq(.gz|.bz2) files
	Due to end-repair incorporation of meth. or unmeth. Cs, selected R2 reads will be 5' clipped by two nucleotides.
	To remove SE or R1 3' end-repaired cutting site, clip off YG prefixed adapter.

	This script is also able to detect and remove diversity adapter sequences as introduced by some sample prep protocols like Methyl-Seq from NuGEN Technologies.

	MspI
	CGG...........C
	  C...........CCG

	End-repair and BS
	RETPADA-CGG...........TCG-ADAPTER
	RETPADA-CGG...........TTG-ADAPTER
	RETPADA-TGG...........TCG-ADAPTER
	RETPADA-TGG...........TTG-ADAPTER
	--->R1                     R2<---

	R1 starts with CGG or TGG
	R2 starts with CGA or CAA

	VERSION
	0.1.0

	SYNOPSIS MERGE
	$(basename $0) -i <fastq1> -j <fastq2> -o <fastq1> -p <fastq2>


	OPTIONS
	-t [value]  | threads ($t) used for output compression using pigz (fallback: gzip)
	-d [value]  | optional. maximum length of diversity adapter
	-i [path]   | input SE or first mate fastq(.gz|bz2)
	-j [path]   | input second mate fastq(.gz|bz2) or unmapped R2 to be SE processed
	-o [path]   | output path for gzip compressed SE or first mate fastq.gz file
	-p [path]   | output path for gzip compressed second mate fastq.gz file
	-h          | this help

	REFERENCES
	(c) Konstantin Riege
	konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 1
}

i=''
j=''
o=''
p=''
h=''
d=0
t=$(cat /proc/cpuinfo | grep -cF processor 2> /dev/null || echo 1)
while getopts i:j:o:p:d:t:h ARG; do
	case $ARG in
		i) i="$OPTARG"; readlink -e "$i" &> /dev/null || { echo "error. can not find $i" >&2; exit 1; };;
		j) j="$OPTARG"; readlink -e "$j" &> /dev/null || { echo "error. can not find $j" >&2; exit 1; };;
		o) o="$OPTARG"; mkdir -p "$(dirname "$o")" || { echo "error. can write to $o" >&2; exit 1; };;
		p) p="$OPTARG"; mkdir -p "$(dirname "$p")" || { echo "error. can write to $p" >&2; exit 1; };;
		d) d=$OPTARG;;
		t) t=$OPTARG;;
		h) (usage); exit 0;;
		*) usage;
	esac
done
[[ $# -eq 0 ]] && usage
[[ $i && ! $o ]] && usage
[[ $j && ! $p ]] && usage


[[ $i ]] && open=$(readlink -e "$i" | file -f - | grep -Eo '(gzip|bzip)' && echo -cd || echo cat) || open=$(readlink -e "$j" | file -f - | grep -Eo '(gzip|bzip)' && echo -cd || echo cat)

if [[ $i ]] && [[ $j ]]; then
	pigz -h 2> /dev/null && z="pigz -k -c -p $(((t+1)/2))" || z="gzip -k -c"

	if [[ $d -gt 0 ]]; then
		# 3rd base of cutting site should be high quality
		paste <($open "$i" | paste - - - -) <($open "$j" | paste - - - -) | perl -F'\t' -slane '$r1=substr($F[1],0,$d+3); next unless $r1=~/(CGG|TGG)/; $r1p=$-[0]; $r2=substr($F[5],0,$d+3); next unless $r2=~/(CAA|CGA)/; $r2p=$-[0]; $F[1]=substr($F[1],$r1p); $F[3]=substr($F[3],$r1p); $F[5]=substr($F[5],$r2p+2); $F[7]=substr($F[7],$r2p+2); print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7])' -- -d=$d > >($z > "$o") 2> >($z > "$p")
	else
		# less stringent R1 motif if R2 motif ok or vice versa
		paste <($open "$i" | paste - - - -) <($open "$j" | paste - - - -) | perl -F'\t' -lane 'next unless $F[1]=~/^(CG\w|TG\w|NGG|CNG|TNG)/; if($F[1]=~/^(C|T)GG/){next unless $F[5]=~/^(CA\w|CG\w|NAA|NGA|CNA)/}else{next unless $F[5]=~/^(CAA|CGA)/} $F[5]=substr($F[5],2); $F[7]=substr($F[7],2); print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7])' > >($z > "$o") 2> >($z > "$p")
	fi
else
	pigz -h &> /dev/null && z="pigz -k -c -p $t" || z="gzip -k -c"

	# for SE or unmapped R1 input
	[[ $i ]] && $open "$i" | paste - - - - | perl -F'\t' -slane '$r1=substr($F[1],0,$d+3); next unless $r1=~/(CGG|TGG)/; $r1p=$-[0]; $F[1]=substr($F[1],$r1p); $F[3]=substr($F[3],$r1p); print join("\n",@F)' -- -d=$d > >($z > "$o")
	# for unmapped R2 input

	[[ $j ]] && $open "$j" | paste - - - - | perl -F'\t' -slane '$r2=substr($F[1],0,$d+3); next unless $r2=~/(CAA|CGA)/; $r2p=$-[0]; $F[1]=substr($F[1],$r2p+2); $F[3]=substr($F[3],$r2p+2); print join("\n",@F)' -- -d=$d > >($z > "$p")
fi

exit $((${PIPESTATUS[@]/%/+}0))