#!/usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$0")/../bashbone_lite.sh" -a "$@" || exit 1

# Taq1
# CGA...........A
#   T...........TCG

# End-repair and BS
# RETPADA-CGA...........TCG-ADAPTER
# RETPADA-CGA...........TTG-ADAPTER
# RETPADA-TGA...........TCG-ADAPTER
# RETPADA-TGA...........TTG-ADAPTER
# --->R1                     R2<---

# R1 starts with CGA or TGA
# R2 starts with CGA or CAA

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename $0) select single or pairs of reads with proper MspI cutting site (single-digest) from fastq(.gz|.bz2) files
		Due to end-repair incorporation of meth. or unmeth. Cs, selected R2 reads will be 5' clipped by two nucleotides.
		This script is also able to detect and remove diversity adapter sequences as introduced by some sample prep protocols like Methyl-Seq from NuGEN Technologies.

		Required post-processing:
		To remove SE or R1 3' end-repaired cutting site, clip off YG or NN prefixed adapter afterwards.
		R2 clipping could also be done afterwards by e.g. utilizing Cutadapt with parameter -U 2.
		-> both tasks, as well as adapter sequence inference, can be solved by e.g. utilizing TrimGalore! with --rrbs option.

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
		0.1.2

		SYNOPSIS MERGE
		$(basename $0) -i <fastq1> -j <fastq2> -o <fastq1> -p <fastq2>


		OPTIONS
		-t [value]  | threads ($t) used for output compression using pigz (fallback: gzip)
		-d [value]  | optional. maximum length of diversity adapter (implies strict motif selection)
		-i [path]   | input SE or first mate fastq(.gz|bz2)
		-j [path]   | input mate pair fastq(.gz|bz2) or without -i to process R2 in SE manner
		-o [path]   | output path for gzip compressed SE or first mate fastq.gz file
		-p [path]   | output path for gzip compressed second mate fastq.gz file
		-c [value]  | number bases to be clipped at R2 5' (default: 2)
		-s          | strict motif selection for PE data i.e. (C|T)GG and C(G|A)A instead of wildcard (C|T)G* and C(G|A)A or (C|T)GG and C(G|A)*
		-v          | invert selection
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
p=''
h=''
d=0
v=false
s=false
t=$(grep -cF processor /proc/cpuinfo)
c=2
while getopts i:j:o:p:d:t:c:hvs ARG; do
	case $ARG in
		i) i="$OPTARG";;
		j) j="$OPTARG";;
		o) o="$OPTARG"; mkdir -p "$(dirname "$o")";;
		p) p="$OPTARG"; mkdir -p "$(dirname "$p")";;
		d) d=$OPTARG;;
		t) t=$OPTARG;;
		c) c=$OPTARG;;
		h) { usage || exit 0; };;
		v) v=true;;
		s) s=true;;
		*) usage;
	esac
done
[[ $# -eq 0 ]] && { usage || exit 0; }
[[ $i && ! $o ]] && usage
[[ $j && ! $p ]] && usage


[[ $i ]] && open=$({ readlink -e "$i" | file -b --mime-type -f - | grep -oF -e 'gzip' -e 'bzip2' || echo cat; } | sed '/cat/!{s/gzip/pigz -p 1/; s/$/ -cd/}') || open=$({ readlink -e "$j" | file -b --mime-type -f - | grep -oF -e 'gzip' -e 'bzip2' || echo cat; } | sed '/cat/!{s/gzip/pigz -p 1/; s/$/ -cd/}')

if [[ $i ]] && [[ $j ]]; then
	pigz -h 2> /dev/null && z="pigz -k -c -p $(((t+1)/2))" || z="gzip -k -c"

	if [[ $d -gt 0 ]]; then
		if $v; then
			paste <($open "$i" | paste - - - -) <($open "$j" | paste - - - -) | perl -F'\t' -slane '$r1=substr($F[1],0,$d+3); unless($r1=~/(CGG|TGG)/){print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7]); next} $r1p=$-[0]; $r2=substr($F[5],0,$d+3); unless($r2=~/(CAA|CGA)/){print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7])}' -- -d=$d > >($z > "$o") 2> >($z > "$p") | cat
		else
			# 3rd base of cutting site should be high quality
			paste <($open "$i" | paste - - - -) <($open "$j" | paste - - - -) | perl -F'\t' -slane '$r1=substr($F[1],0,$d+3); next unless $r1=~/(CGG|TGG)/; $r1p=$-[0]; $r2=substr($F[5],0,$d+3); next unless $r2=~/(CAA|CGA)/; $r2p=$-[0]; $F[1]=substr($F[1],$r1p); $F[3]=substr($F[3],$r1p); $F[5]=substr($F[5],$r2p+$c); $F[7]=substr($F[7],$r2p+$c); print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7])' -- -d=$d -c=$c > >($z > "$o") 2> >($z > "$p") | cat
		fi
	else
		if $s; then
			if $v; then
				paste <($open "$i" | paste - - - -) <($open "$j" | paste - - - -) | perl -F'\t' -lane 'unless($F[1]=~/^(C|T)GG/){print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7]); next} unless($F[5]=~/^C(G|A)A/){print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7])}' > >($z > "$o") 2> >($z > "$p") | cat
			else
				paste <($open "$i" | paste - - - -) <($open "$j" | paste - - - -) | perl -F'\t' -slane 'next unless $F[1]=~/^(C|T)GG/; next unless $F[5]=~/^C(GA)A/; $F[5]=substr($F[5],$c); $F[7]=substr($F[7],$c); print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7])' -- -c=$c > >($z > "$o") 2> >($z > "$p") | cat
			fi
		else
			if $v; then
				paste <($open "$i" | paste - - - -) <($open "$j" | paste - - - -) | perl -F'\t' -lane 'unless($F[1]=~/^(CG\w|TG\w|NGG|CNG|TNG)/){print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7]); next} if($F[1]=~/^(C|T)GG/){unless($F[5]=~/^(CA\w|CG\w|NAA|NGA|CNA)/){print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7]); next}}else{unless($F[5]=~/^(CAA|CGA)/){print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7])}}' > >($z > "$o") 2> >($z > "$p") | cat
			else
				# less stringent R1 motif if R2 motif ok or vice versa
				paste <($open "$i" | paste - - - -) <($open "$j" | paste - - - -) | perl -F'\t' -slane 'next unless $F[1]=~/^(CG\w|TG\w|NGG|CNG|TNG)/; if($F[1]=~/^(C|T)GG/){next unless $F[5]=~/^(CA\w|CG\w|NAA|NGA|CNA)/}else{next unless $F[5]=~/^(CAA|CGA)/} $F[5]=substr($F[5],$c); $F[7]=substr($F[7],$c); print join("\n",@F[0..3]); print STDERR join("\n",@F[4..7])' -- -c=$c > >($z > "$o") 2> >($z > "$p") | cat
			fi
		fi
	fi
else
	pigz -h &> /dev/null && z="pigz -k -c -p $t" || z="gzip -k -c"

	if $v; then
		if [[ $i ]]; then
			$open "$i" | paste - - - - | perl -F'\t' -slane '$r1=substr($F[1],0,$d+3); unless($r1=~/(CGG|TGG)/){print join("\n",@F)}' -- -d=$d > >($z > "$o") | cat
		else
			$open "$j" | paste - - - - | perl -F'\t' -slane '$r2=substr($F[1],0,$d+3); unless($r2=~/(CAA|CGA)/){print join("\n",@F)}' -- -d=$d > >($z > "$p") | cat
		fi
	else
		# for SE or unmapped R1 input
		if [[ $i ]]; then
			$open "$i" | paste - - - - | perl -F'\t' -slane '$r1=substr($F[1],0,$d+3); next unless $r1=~/(CGG|TGG)/; $r1p=$-[0]; $F[1]=substr($F[1],$r1p); $F[3]=substr($F[3],$r1p); print join("\n",@F)' -- -d=$d > >($z > "$o") | cat
		else
			# for unmapped R2 input
			$open "$j" | paste - - - - | perl -F'\t' -slane '$r2=substr($F[1],0,$d+3); next unless $r2=~/(CAA|CGA)/; $r2p=$-[0]; $F[1]=substr($F[1],$r2p+$c); $F[3]=substr($F[3],$r2p+$c); print join("\n",@F)' -- -d=$d -c=$c > >($z > "$p") | cat
		fi
	fi
fi
