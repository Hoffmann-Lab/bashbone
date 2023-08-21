#! /usr/bin/env bash
# (c) Konstantin Riege

function quantify::featurecounts(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-x <strandness> | hash per bam of
			-g <gtf>      | path to
			-l <level>    | feature (default: exon)
			-f <feature>  | feature (default: gene, needs <feature>_id tag)
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}" gtf level="exon" feature="gene"
	declare -n _mapper_featurecounts _strandness_featurecounts
	while getopts 'S:s:t:r:x:g:l:f:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			r) ((++mandatory)); _mapper_featurecounts=$OPTARG;;
			x) ((++mandatory)); _strandness_featurecounts=$OPTARG;;
			g) ((++mandatory)); gtf="$OPTARG";;
			l) level=$OPTARG;;
			f) feature=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "quantifying reads"

	# featurecounts cannot handle more than 64 threads
	declare -n _bams_featurecounts="${_mapper_featurecounts[0]}"
	local instances ithreads i=$((${#_mapper_featurecounts[@]}*${#_bams_featurecounts[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -i $i -t 64 -T $threads)
	[[ $instances -lt $i && $ithreads -gt 64 ]] && ((++instances)) && ithreads=$((threads/instances))
	[[ $ithreads -gt 64 ]] && ithreads=64

	declare -a cmdchk=("featureCounts -v |& grep -oE 'v[.0-9]+'")
	local version=$(commander::runcmd -c subread -a cmdchk)
	[[ "$(echo -e "v2.0.1\n$version" | sort -Vr | head -1)" == "v2.0.1" ]] && version=" " || version="--countReadPairs "

	declare -a tdirs cmd1
	local m f o params x
	for m in "${_mapper_featurecounts[@]}"; do
		declare -n _bams_featurecounts=$m
		mkdir -p "$outdir/$m"
		for f in "${_bams_featurecounts[@]}"; do
			o="$outdir/$m/$(basename "$f" .bam).${feature}counts"
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.featurecounts)")

			# infer SE or PE
			params=''
			x=$(samtools view -F 4 "$f" | head -10000 | cat <(samtools view -H "$f") - | samtools view -c -f 1)
			[[ $x -gt 0 ]] && params+="-p $version"
			[[ "$feature" == "$level" ]] && params+='-f -O '

			# use -M and --ignoreDup to count also multi-mappings and dulpicates ie. everythin given from bam
			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				featureCounts
					$params
					-Q 0
					--minOverlap 10
					--maxMOp 999
					-s ${_strandness_featurecounts["$f"]}
					-T $ithreads
					-t $level
					-g ${feature}_id
					-M
					--ignoreDup
					--tmpDir "${tdirs[-1]}"
					-a "$gtf"
					-o "$o"
					"$f"
			CMD
				awk 'NR>2{
					print \$1"\t"\$NF
				}' $o > $o.htsc
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c subread -v -b -i $instances -a cmd1
	fi

	return 0
}

function quantify::normalize(){
	quantify::tpm "$@"
}

function quantify::tpm(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-g <gtf>      | path to
			-i <countsdir>| path to
			-l <level>    | feature (default: exon)
			-f <feature>  | feature (default: gene, needs <feature>_id tag)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads countsdir gtf level="exon" feature="gene" tmpdir="${TMPDIR:-/tmp}"
	declare -n _mapper_tpm
	while getopts 'S:s:t:r:g:i:l:f:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			r) ((++mandatory)); _mapper_tpm=$OPTARG;;
			g) ((++mandatory)); gtf="$OPTARG";;
			i) ((++mandatory)); countsdir="$OPTARG";;
			l) level="$OPTARG";;
			f) feature="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	declare -n _bams_tpm=${_mapper_tpm[0]}
	if [[ ${#_bams_tpm[@]} -gt 1 ]]; then
		commander::printinfo "calculating transcripts per million and variance stabilized counts"
	else
		commander::printinfo "calculating transcripts per million"
	fi

	local m f countfile header sample csv
	declare -a cmd1 cmd2
	for m in "${_mapper_tpm[@]}"; do
		declare -n _bams_tpm=$m
		csv="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.csv)"
		header="id"
		for f in "${_bams_tpm[@]}"; do
			sample="$(basename "${f%.*}")"
			header+="\t$sample"
			countfile="$(find -L "$countsdir/$m" -maxdepth 1 -name "$sample*.${feature}counts.htsc" -print -quit | grep .)"
			echo "$countfile,$countfile" >> "$csv"

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
				tpm.pl "$gtf" $level ${feature}_id "$countfile" > "$countfile.tpm"
			CMD
		done

		if [[ ${#_bams_tpm[@]} -gt 1 ]]; then
			commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				Rscript - <<< '
					options(warn=-1);
					suppressMessages(library("DESeq2"));
					suppressMessages(library("argparser"));
					args = arg_parser("Calculate variance stabilized counts", hide.opts=T);
					args = add_argument(args, "--csv", short="-c", help="path to input csv", flag=F);
					args = parse_args(args);
					dds = DESeqDataSetFromHTSeqCount(sampleTable=read.table(args$csv, header=F, sep=",", stringsAsFactors=F, check.names=F, quote=""), directory="", design=~1);
					dds = estimateSizeFactors(dds);
					vsd = rbind(colnames(dds),assay(varianceStabilizingTransformation(dds, blind=FALSE)));
					ret = apply(vsd, 2, function(x){ write.table(data.frame(rownames(dds), x[2:length(x)], check.names=F), file=paste0(x[1],".vsc"), col.names=F, row.names=F, quote=F, sep='\t') });
				'
			CMD
				-c "$csv"
			CMD
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $threads -a cmd1
	fi

	return 0
}

function quantify::bamcoverage(){
	# bedg != samtools bedcov <- which is sum(histogram = samtools depth), use bedtools multicov
	# use
	# bedtools multicov -D -F -bams <bam>[ <bam>..] -bed <bins> > <bedg> (slow)
	# or
	# bamCoverage --bam $f -o <*.bw|*.bedg> -p 24 -bs <binsize> (merges adjacent bins of same coverage into up to 5000000bp bins, ignores mates with an insert >1000bp)
	# or
	# featureCounts from SAF formatted bins, which allows also fractional counts

	# bamCoverage                      featureCounts
	# ...                              ...
	# chrY  26672600  26672800  3      chrY  26672600  26672800  3
	# chrY  26672800  30000000  0      chrY  26672800  26673000  0
	# chrY  30000000  35000000  0      chrY  26673000  26673200  0
	# chrY  35000000  40000000  0      ...
	# chrY  40000000  45000000  0
	# chrY  45000000  46674400  0      chrY  46674200  46674400  0
	# chrY  46674400  46674600  16     chrY  46674400  46674600  16
	# ...

	# bug in v<=2.0.3 with fractional counts on largest overlaps - read spanning two features get 0 and 0.5 instead of 0 and 1

	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			description:
			assign reads to bins or bed/bed-like file given features and write counts plus TPM/BPM normalized to bedGraph and bigWig files
			- all reads are counted! including:
			- singletons, duplicated, multi-mapped, bad-quality, cross-chromosome, circular and highly mismatched ones

			${FUNCNAME[1]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-r <mapper>     | array of bams within array of
			-x <strandness> | hash per bam of. default: 0
			                  0 - unstranded
			                  1 - stranded
			                  2 - reversely stranded
			-o <overlap>    | of a read, that spans multiple features, which assigns it to one or multiple feature. default: 1
			                  0   - largest
			                  0-1 - fraction
			                  1-N - basepairs
			-w <windowsize> | for bins to be quantified, unless given a bed/bed-like file by -b option. default: 1
			-b <bed>        | or bed-like file with 3 or 6 columns. latter for strand-specific quantification. see also -w
			-i              | index of -b option given file is 1 and not 0 (e.g. narrowPeak, gtf, gff)
			-f              | count reads fractional if, according to -o option, they can be assigned to multiple features (+1/n)
			                  NOTE: not supported for version < 2.0.4, when combined with -o 0
			-u              | only unambiguously assigend reads will be counted
			-p              | count read-pairs/fragments instead of reads (experimental!)
			                  NOTE1: +1 instead of +2 is counted only if mates are unambiguously assigned to same feature, else +1,+1
			                  NOTE2: fractional quantification is not available
			                  NOTE3: unless -o option corresponds to basepairs, assignment is derived from the sum of both mate lengths
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false tmpdir="${TMPDIR:-/tmp}" threads overlap=1 windowsize=1 bed index=0 fractional=false unambiguous=false pairs=false
	declare -n _mapper_bam2bedg _strandness_bam2bedg
	while getopts 'S:s:t:r:x:o:w:b:ifup' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			r) ((++mandatory)); _mapper_bam2bedg=$OPTARG;;
			x) _strandness_bam2bedg=$OPTARG;;
			o) overlap=$OPTARG;;
			w) windowsize=$OPTARG;;
			b) bed="$OPTARG";;
			i) index=1;;
			f) fractional=true;;
			u) unambiguous=true;;
			p) pairs=true;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	declare -n _bams_bam2bedg="${_mapper_bam2bedg[0]}"
	local instances ithreads i=$((${#_mapper_bam2bedg[@]}*${#_bams_bam2bedg[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -i $i -t 64 -T $threads)
	[[ $instances -lt $i && $ithreads -gt 64 ]] && ((++instances)) && ithreads=$((threads/instances))
	[[ $ithreads -gt 64 ]] && ithreads=64

	declare -a cmdchk=("featureCounts -v |& grep -oE 'v[.0-9]+'")
	local version=$(commander::runcmd -c subread -a cmdchk)
	$fractional && [[ $overlap -eq 0 ]] && [[ "$(echo -e "v2.0.3\n$version" | sort -Vr | head -1)" == "v2.0.3" ]] && _usage
	[[ "$(echo -e "v2.0.1\n$version" | sort -Vr | head -1)" == "v2.0.1" ]] && version="old" || version="new"

	commander::printinfo "convertig bam to bedg/bw"

	local m f params x o
	declare -a cmd1 cmd2 cmd3 tdirs
	for m in "${_mapper_bam2bedg[@]}"; do
		declare -n _bams_bam2bedg=$m
		for f in "${_bams_bam2bedg[@]}"; do
			o="${f%.*}"
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bam2bedg)")

			if [[ $bed ]]; then
				commander::makecmd -a cmd1 -s ' ' -o "${tdirs[-1]}/regions.saf" -c {COMMANDER[0]}<<- CMD {COMMANDER[2]}<<- 'CMD'
					samtools view -H "${_bams_bam2bedg[0]}" | sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}' | helper::sort -t $threads -k1,1 > "${tdirs[-1]}/chr.sizes";
					helper::sort -f "$bed" -t $threads -k1,1 -k2,2n -k3,3n | awk -v i=$index -v OFS='\t'
				CMD
					'BEGIN{i=i==1?0:1}{id=$1":"$2"-"$3; if($4){id=$4}; s="."; if($6){s=$6}; print id,$1,$2+i,$3,".",".",s}'
				CMD
			else
				commander::makecmd -a cmd1 -s '|' -o "${tdirs[-1]}/regions.saf" -c {COMMANDER[0]}<<- CMD {COMMANDER[2]}<<- 'CMD'
					samtools view -H "${_bams_bam2bedg[0]}" | sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}' | helper::sort -t $threads -k1,1 > "${tdirs[-1]}/chr.sizes";
					bedtools makewindows -w $windowsize -g "${tdirs[-1]}/chr.sizes"
				CMD
					awk -v OFS='\t' '{id=$1":"$2"-"$3; print id,$1,$2+1,$3,".",".","."}'
				CMD
			fi

			params=''
			$unambiguous || params+=" -O"
			$fractional && params+=" --fraction"
			params+=" $(echo $overlap | awk '{if($1==0){print "--largestOverlap"}else{print $1<1? "--fracOverlap "$1 : "--minOverlap "$1}}')"
			if [[ $version == "old" ]]; then
				$pairs && params+=" -p"
			else
				x=$(samtools view -F 4 "$f" | head -10000 | cat <(samtools view -H "$f") - | samtools view -c -f 1)
				if [[ $x -gt 0 ]]; then
					params+=" -p"
					$pairs && params+=" --countReadPairs"
				fi
			fi

			commander::makecmd -a cmd2 -s '|' -o "$o.coverage.tpm.bed" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD'
				featureCounts
					$params
					-Q 0
					--maxMOp 999
					-s $(declare -p ${!_strandness_bam2bedg} &> /dev/null && echo ${_strandness_bam2bedg["$f"]} || echo 0)
					-T $ithreads
					-f
					-M
					--ignoreDup
					--tmpDir "${tdirs[-1]}"
					-F SAF
					-a "${tdirs[-1]}/regions.saf"
					-o /dev/stdout
					"$f"
			CMD
				tee -ia >(awk -v OFS='\t' 'NR>2{print \$2,\$3-1,\$4,sprintf("%d",\$7==int(\$7)?\$7:\$7+1)}' > "$o.coverage.counts.bed")
			CMD
				awk '
					NR>2{
						o[NR-3]=$2"\t"$3-1"\t"$4;
						v[NR-3]=sprintf("%d",$7==int($7)?$7:$7+1)/$6;
						sum=sum+v[NR-3];
					}
					END{
						sum=sum/1000000;
						for(i=0;i<=NR-3;i++){
							print o[i]"\t"v[i]/sum;
						}
					}
				'
			CMD

			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
				bedGraphToBigWig "$o.coverage.counts.bed" "${tdirs[-1]}/chr.sizes" "$o.coverage.counts.bw"
			CMD
			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
				bedGraphToBigWig "$o.coverage.tpm.bed" "${tdirs[-1]}/chr.sizes" "$o.coverage.tpm.bw"
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -i 1 -a cmd1
		commander::runcmd -c subread -v -b -i $instances -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
	fi

	return 0
}
