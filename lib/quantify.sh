#! /usr/bin/env bash
# (c) Konstantin Riege

quantify::featurecounts() {
	declare -a tdirs
	_cleanup::quantify::featurecounts(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
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
			-p <tmpdir>   | path to
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir gtf level="exon" feature="gene"
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
			f) feature=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "quantifying reads"

	# featurecounts cannot handle more than 64 threads
	local instances ithreads m f
	for m in "${_mapper_featurecounts[@]}"; do
		declare -n _bams_featurecounts=$m
		((instances+=${#_bams_featurecounts[@]}))
	done
	instances=$((instances*2))
	read -r instances ithreads < <(configure::instances_by_threads -i $((instances==0?1:instances)) -t 64 -T $threads)

	declare -a cmdchk=("featureCounts -v |& grep -oE 'v[.0-9]+'")
	local version=$(commander::runcmd -c subread -a cmdchk)
	[[ "$(echo -e "v2.0.1\n$version" | sort -Vr | head -1)" == "v2.0.1" ]] && version=" " || version="--countReadPairs "

	declare -a cmd1
	local mf f o params x
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

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				featureCounts
					$params
					-Q 0
					--minOverlap 10
					-s ${_strandness_featurecounts["$f"]}
					-T $ithreads
					-t $level
					-g ${feature}_id
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

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -c subread -v -b -i $instances -a cmd1
	fi

	return 0
}

quantify::tpm() {
	_usage() {
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

	local OPTIND arg mandatory skip=false threads countsdir gtf level="exon" feature="gene"
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

	commander::printinfo "calculating transcripts per million"

	local m f countfile header sample
	declare -a cmd1 cmd2
	for m in "${_mapper_tpm[@]}"; do
		declare -n _bams_tpm=$m
		header="id"
		for f in "${_bams_tpm[@]}"; do
			sample="$(basename "${f%.*}")"
			header+="\t$sample"
			countfile="$(find -L "$countsdir/$m" -maxdepth 1 -name "$sample*.${feature}counts.htsc" -print -quit | grep .)"

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
				tpm.pl "$gtf" $level ${feature}_id "$countfile" > "$countfile.tpm"
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -i $threads -a cmd1
	fi
	#commander::runcmd -v -b -i $threads -a cmd2

	return 0
}
