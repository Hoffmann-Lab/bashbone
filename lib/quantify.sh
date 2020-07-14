#! /usr/bin/env bash
# (c) Konstantin Riege

quantify::featurecounts() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-5 <skip>     | true/false md5sums, gtf prep respectively
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-g <gtf>      | path to
			-l <level>    | feature (default: exon)
			-f <tag>      | feature (default: gene_id)
			-p <tmpdir>   | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads outdir tmpdir gtf level featuretag
	declare -n _mapper_featurecounts
	while getopts 'S:s:5:t:r:g:l:f:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			5) $OPTARG && skipmd5=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			r) ((mandatory++)); _mapper_featurecounts=$OPTARG;;
			g) ((mandatory++)); gtf="$OPTARG";;
			l) level=$OPTARG;;
			f) featuretag=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done

	[[ $mandatory -lt 5 ]] && _usage && return 1

	commander::print "inferring library preparation method and quantifying reads"

	$skipmd5 && {
		commander::warn "skip checking md5 sums and thus annotation preparation"
	} || {
		commander::print "checking md5 sums"
		local thismd5gtf
		[[ -s $gtf ]] && thismd5gtf=$(md5sum "$gtf" | cut -d ' ' -f 1)
		if [[ ! -s ${gtf%.*}.transcripts.gtf ]] || [[ "$thismd5gtf" && "$thismd5gtf" != "$md5gtf" ]]; then
			commander::print "preparing annotation for exon level quantification"

			declare -a cmdprep
			commander::makecmd -a cmdprep -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				dexseq_prepare_annotation2.py
					-r no
					-f "$tmpdir/tmp.gtf"
					$gtf /dev/null
			CMD
				sed -r 's/(.+gene_id\s+")([^"]+)(.+exon_number\s+")([^"]+)(.+)/\1\2\3\4\5; exon_id "\2:\4"/' "$tmpdir/tmp.gtf" > "${gtf%.*}.transcripts.gtf"
			CMD
				rm -f "$tmpdir/tmp.gtf"
			CMD

			commander::runcmd -v -b -t $threads -a cmdprep || {
				commander::printerr "$funcname failed at annotation preparation"
				return 1
			}
		fi
	}

	declare -a cmd1 cmd2
	local m f
	for m in "${_mapper_featurecounts[@]}"; do
		declare -n _bams_featurecounts=$m
		mkdir -p "$tmpdir/$m"
		for f in "${_bams_featurecounts[@]}"; do
			alignment::_inferexperiment \
				-1 cmd1 \
				-2 cmd2 \
				-p "$tmpdir/$m" \
				-i "$f" \
				-g "$gtf"
		done
	done

	declare -A strandness

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	local l
			declare -a a mapdata
			commander::runcmd -v -b -t $threads -a cmd1 && \
			echo ":INFO: running commands of array cmd2" && \
			commander::printcmd -a cmd2 && \
			conda activate py3 && \
			mapfile -t mapdata < <(commander::runcmd -t $threads -a cmd2)
			for l in "${mapdata[@]}"; do
				a=($l)
				strandness["${a[@]:1}"]="${a[0]}"
			done
			conda activate py2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	# featurecounts cannot handle more than 64 threads
	local instances ithreads m f
	for m in "${_mapper_featurecounts[@]}"; do
		declare -n _bams_featurecounts=$m
		((instances+=${#_bams_featurecounts[@]}))
	done
	instances=$((instances*2))
	# $((instances==0?1:instances)) in case of misuse of this function for annotation preparation without bam files
	read -r instances ithreads < <(configure::instances_by_threads -i $((instances==0?1:instances)) -t 64 -T $threads)

	declare -a cmd3
	for m in "${_mapper_featurecounts[@]}"; do
		declare -n _bams_featurecounts=$m
		mkdir -p "$outdir/$m"
		for f in "${_bams_featurecounts[@]}"; do
			o=$outdir/$m/$(basename $f)
			quantify::_featurecounts \
				-1 cmd3 \
				-t $ithreads \
				-i "$f" \
				-s ${strandness["$f"]:-?} \
				-g "$gtf" \
				-l ${level:-exon} \
				-f ${featuretag:-gene_id} \
				-p "$tmpdir" \
				-o "${o%.*}.counts"
			quantify::_featurecounts \
				-1 cmd3 \
				-t $ithreads \
				-i "$f" \
				-s ${strandness["$f"]:-?} \
				-g "${gtf%.*}.transcripts.gtf" \
				-l ${level:-exon} \
				-f exon_id \
				-p "$tmpdir" \
				-o "${o%.*}.exoncounts"
		done
	done

	$skip && {
		commander::printcmd -a cmd3
	} || {
		{	conda activate py2r && \
			commander::runcmd -v -b -t $instances -a cmd3 && \
			conda activate py2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

quantify::_featurecounts() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-1 <cmds1>      | array of
			-t <threads>    | number of
			-s <strandness> | [0|1|2]
			                  (0) unstranded, 
			                  (1) FR stranded, 
			                  (2) FR reversely stranded
			-i <bam>        | path to
			-g <gtf>        | path to
			-l <level>      | feature (default: exon)
			-f <tag>        | feature (default: gene_id)
			-p <tmpdir>     | path to
			-o <outfile>    | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory threads bam strandness gtf outfile level featuretag
	declare -n _cmds1_featurecounts
	while getopts '1:t:s:i:g:l:f:o:p:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_featurecounts=$OPTARG;;
			t) ((mandatory++)); threads=$OPTARG; [[ $threads -gt 64 ]] && threads=64;;
			s) ((mandatory++)); strandness=$OPTARG;;
			i) ((mandatory++)); bam="$OPTARG";;
			g) ((mandatory++)); gtf="$OPTARG";;
			l) level=$OPTARG;;
			f) featuretag=$OPTARG;;
			o) ((mandatory++)); outfile="$OPTARG";;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1

	# infer SE or PE
	local params=''
	[[ $(samtools view -F 4 -h "$bam" | head -10000 | samtools view -c -f 1) -gt 0 ]] && params+='-p '
	[[ "$featuretag" != "gene_id" ]] && params+='-f -O '

	commander::makecmd -a _cmds1_featurecounts -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
		featureCounts
			$params
			-Q 0
			--minOverlap 10
			-s $strandness
			-T $threads
			-t ${level:-exon}
			-g ${featuretag:-gene_id}
			--tmpDir "$tmpdir"
			-a "$gtf"
			-o "$outfile"
			$bam
	CMD
		awk 'NR>2{
			print \$1"\t"\$NF
		}' $outfile > $outfile.htsc
	CMD

	return 0
}

quantify::tpm() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
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
			t) ((mandatory++)); threads=$OPTARG;;
			r) ((mandatory++)); _mapper_tpm=$OPTARG;;
			g) ((mandatory++)); gtf="$OPTARG";;
			i) ((mandatory++)); countsdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage && return 1

	commander::print "calculating transcripts per million"

	local m f countfile
	declare -a cmd1
	for m in "${_mapper_tpm[@]}"; do
		declare -n _bams_tpm=$m
		for f in "${_bams_tpm[@]}"; do
			countfile="$countsdir/$m/$(basename $f)"
			countfile=$(readlink -e "${countfile%.*}"*.counts.+(reduced|htsc) | head -1)
			quantify::_tpm \
				-1 cmd1 \
				-g $gtf \
				-i $countfile
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

quantify::_tpm() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-1 <cmds1>     | array of
			-g <gtf>       | path to
			-i <countfile> | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory countfile gtf
	declare -n _cmds1_tpm
	while getopts '1:i:g:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_tpm=$OPTARG;;
			g) ((mandatory++)); gtf="$OPTARG";;
			i) ((mandatory++)); countfile="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage && return 1

	commander::makecmd -a _cmds1_tpm -s '|' -c {COMMANDER[0]}<<- CMD
		tpm.pl $gtf $countfile > $countfile.tpm
	CMD

	return 0
}
