#! /usr/bin/env bash
# (c) Konstantin Riege

quantify::featurecounts() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-5 <skip>     | true/false md5sums, gtf prep respectively
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-x <strandness> | hash per bam of
			-g <gtf>      | path to
			-l <level>    | feature (default: exon)
			-f <tag>      | feature (default: gene_id)
			-p <tmpdir>   | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads outdir tmpdir gtf level="exon" featuretag="gene_id"
	declare -n _mapper_featurecounts _strandness_featurecounts
	while getopts 'S:s:t:r:x:g:l:f:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			r) ((++mandatory)); _mapper_featurecounts=$OPTARG;;
			x) ((++mandatory)); _strandness_featurecounts=$OPTARG;;
			g) ((++mandatory)); gtf="$OPTARG";;
			l) level=$OPTARG;;
			f) featuretag=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage && return 1

	# featurecounts cannot handle more than 64 threads
	local instances ithreads m f
	for m in "${_mapper_featurecounts[@]}"; do
		declare -n _bams_featurecounts=$m
		((instances+=${#_bams_featurecounts[@]}))
	done
	instances=$((instances*2))
	read -r instances ithreads < <(configure::instances_by_threads -i $((instances==0?1:instances)) -t 64 -T $threads)

	declare -a cmd1 tdirs
	local mf f o params
	for m in "${_mapper_featurecounts[@]}"; do
		declare -n _bams_featurecounts=$m
		mkdir -p "$outdir/$m"
		for f in "${_bams_featurecounts[@]}"; do
			o="$outdir/$m/$(basename "$f" .bam).${featuretag/_id/}counts"
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.featurecounts)")

			# infer SE or PE
			params=''
			[[ $(samtools view -F 4 "$bam" | head -10000 | cat <(samtools view -H "$bam") - | samtools view -c -f 1) -gt 0 ]] && params+='-p '
			[[ "$featuretag" != "gene_id" ]] && params+='-f -O '

			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				featureCounts
					$params
					-Q 0
					--minOverlap 10
					-s ${_strandness_featurecounts["$f"]}
					-T $ithreads
					-t $level
					-g $featuretag
					--tmpDir "${tdirs[-1]}"
					-a "$gtf"
					-o "$o"
					$f
			CMD
				awk 'NR>2{
					print \$1"\t"\$NF
				}' $o > $o.htsc
			CMD
		done
	done

	$skip && {
		commander::printcmd -a cmd1
	} || {
		{	commander::runcmd -c subread -v -b -t $instances -a cmd1
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}

quantify::tpm() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-g <gtf>      | path to
			-i <countsdir>| path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads countsdir gtf
	declare -n _mapper_tpm
	while getopts 'S:s:t:r:g:i:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			r) ((++mandatory)); _mapper_tpm=$OPTARG;;
			g) ((++mandatory)); gtf="$OPTARG";;
			i) ((++mandatory)); countsdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage && return 1

	commander::printinfo "calculating transcripts per million"

	local m f countfile
	declare -a cmd1
	for m in "${_mapper_tpm[@]}"; do
		declare -n _bams_tpm=$m
		for f in "${_bams_tpm[@]}"; do
			countfile="$countsdir/$m/$(basename $f)"
			countfile=$(readlink -e "${countfile%.*}"*.+(genecounts|counts).+(reduced|htsc) | head -1)
			commander::makecmd -a _cmds1_tpm -s '|' -c {COMMANDER[0]}<<- CMD
				tpm.pl "$gtf" "$countfile" > "$countfile.tpm"
			CMD
		done
	done

	$skip && {
		commander::printcmd -a cmd1
	} || {
		{	commander::runcmd -v -b -t $threads -a cmd1
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
