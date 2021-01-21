#!/usr/bin/env bash
# (c) Konstantin Riege

[[ $# -lt 2 ]] && {
cat << EOF
usage:
$(basename $0) <gtf> <fasta|fai|bam|sam> [<exon|gene> [+|-]<offset>]

version:
0.3.0

description:
- given a gtf file, this script merges exons on transcript level and prints exons, introns, genes, and intergenic regions in gtf format to stdout
- a fasta or fasta index (fai) or any matching bam/sam file is used to infer chromosome sizes
- optionally an offset first elongates gene or exon features up (-) and/or downstream (+) by optional [+|-] offset prefix

output:
- merged exons of multiple transcripts per gene
- introns according to merged exons
- genes according to merged exons
- intergenic regions according to merged genes

requirements:
- input file paths must not include white spaces or other special characters!
- needs mergexons.pl script next to this script
- needs bedtools to be in PATH (check: $([[ $(samtools view &> /dev/null; echo $?) -eq 0 ]] && echo ok || echo failed))
- needs samtools to be in PATH in case of sam/bam input (check: $([[ $(bedtools &> /dev/null; echo $?) -eq 0 ]] && echo ok || echo failed))

EOF
    exit 0
}

error(){
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && echo "$1" 1>&2 && exit 1
}

[[ $(bedtools &> /dev/null; echo $?) -gt 0 ]] && echo "bedtools requiered to be in PATH" 1>&2 && exit 1
gtf=${1:?"error: gtf necessary as first argument"} || exit 1
fabam=${2:?"error: fasta, fai or bam or sam file necessary as second argument"} || exit 1
feature=${3:-gene}
offset=${4:-0}

tmp=$(mktemp --suffix=.gtf2features)
error "cannot create tmp files"
trap 'rm -f $tmp.*' INT TERM EXIT
echo "temp file prefix: $tmp" 1>&2

#create 0-based-start-pos bed
echo "read chromosomes" 1>&2
if [[ $(head -1 $fabam) =~ ^\> ]]; then
	awk '{if($0~/^>/){if(c>0){print o"\t0\t"c"\t.\t.\t+"; print o"\t0\t"c"\t.\t.\t-"} o=substr($1,2,length($1)); c=0}else{c=c+length($0)}}END{print o"\t0\t"c"\t.\t.\t+"; print o"\t0\t"c"\t.\t.\t-"}' $fabam > $tmp.chromosomes
	error "chromosome extraction failed"
elif [[ $(head -1 $fabam | awk '{print NF}') -eq 5 ]]; then
	awk -v OFS='\t' '{print $1,"0",$2,".",".","+"; print $1,"0",$2,".",".","-"}' $fabam > $tmp.chromosomes
	error "chromosome extraction failed"
else
	[[ $(samtools view &> /dev/null; echo $?) -gt 0 ]] && echo "samtools requiered to be in PATH" 1>&2 && exit 1
	samtools view -H $fabam | awk '/^@SQ/{$1=""; print}' | cut -d ':' -f 2,3 | sed -r 's/(\S+)\s+LN:(.+)/\1\t0\t\2\t.\t.\t+\n\1\t0\t\2\t.\t.\t-/' > $tmp.chromosomes
	error "chromosome extraction failed"
fi
cut -f 1,3 $tmp.chromosomes | uniq > $tmp.chromosomesinfo

if [[ $offset != "0" ]]; then
	echo "enlarge ${feature}s" 1>&2
	params="-b $offset"
	[[ "$offset" =~ ^+ ]] && params="-s -l 0 -r ${offset/#+/}"
	[[ "$offset" =~ ^- ]] && params="-s -l ${offset/#-/} -r 0"
	awk -v f=$feature '$3==f' $gtf | sort -k1,1 -k4,4n -k5,5n | bedtools slop -i - -g $tmp.chromosomesinfo $params > $tmp.gtf
	error "bedtools failed on enlargement"
	awk -v f=$feature '$3!=f' $gtf >> $tmp.gtf
	gtf=$tmp.gtf
fi

echo "merge exons" 1>&2
# perl script takes care of sloped genes to adapt first and/or last exon, gene according to sloped exons respectively
# does not require sorted input
$(dirname $(readlink -e $0))/mergexons.pl $gtf > $tmp.final
error "mergexons.pl failed"


echo "extract intergenic regions" 1>&2
#returns 0-based-start-pos bed
awk '$3=="gene"' $tmp.final | sort -k1,1 -k4,4n -k5,5n | bedtools merge -s -i - -c 6,6,7 -o distinct -delim ';' > $tmp.genes
error "bedtools failed on genes"

#take care of returned 0-based-start-pos -> $F[1]++
bedtools subtract -s -a $tmp.chromosomes -b $tmp.genes | perl -lane '$i{$F[0]}++; $F[1]++; print join"\t",($F[0],"merger","intergenic",@F[1..3],$F[5],".","feature_id \"intergenic:".$i{$F[0]}."@".$F[0]."\";")' >> $tmp.final
error "bedtools failed on intergenic regions"

echo "sort results" 1>&2
sort -k1,1 -k4,4n -k5,5n $tmp.final
error "final sort failed"

echo "success" 1>&2
exit 0