#! /usr/bin/env bash
# (c) Konstantin Riege

fusions::starfusion(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip>     | true/false return
			-s <softskip>     | true/false only print commands
			-t <threads>      | number of
			-g <genome>       | path to
			-o <outdir>       | path to
			-1 <fastq1>       | array of
			-2 <fastq2>       | array of
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads genome outdir
	declare -n _fq1_starfusion _fq2_starfusion
	while getopts 'S:s:t:g:o:1:2:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((mandatory++)); threads=$OPTARG;;
			g)	((mandatory++)); genome="$OPTARG";;
			o)	((mandatory++)); outdir="$OPTARG/starfusion"; mkdir -p "$outdir" || return 1;;
			1)	((mandatory++)); _fq1_starfusion=$OPTARG;;
			2)	_fq2_starfusion=$OPTARG;;
			*)	_usage; return 1;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage && return 1
	commander::printinfo "detecting gene fusions star-fusion"

	# as of v1.8 and still true for v1.9, STAR 2.7.2b is requiered
	# NEEDS CTAT_genome_lib
	# STAR-Fusion can use multimapped fusions and uses the following STAR paremeters
	# --chimSegmentMin 12 \  # ** essential to invoke chimeric read detection & reporting **
	# --chimJunctionOverhangMin 8 \
	# --chimOutJunctionFormat 1 \   # **essential** includes required metadata in Chimeric.junction.out file.
	# --alignSJDBoverhangMin 10 \
	# --alignMatesGapMax 100000 \   # avoid readthru fusions within 100k
	# --alignIntronMax 100000 \
	# --alignSJstitchMismatchNmax 5 -1 5 5 \   # settings improved certain chimera detections
	# --outSAMattrRGline ID:GRPundef \
	# --chimMultimapScoreRange 3 \
	# --chimScoreJunctionNonGTAG -4 \
	# --chimMultimapNmax 20 \
	# --chimNonchimScoreDropMin 10 \
	# --peOverlapNbasesMin 12 \
	# --peOverlapMMp 0.1 \
	# --alignInsertionFlush Right \
	# --alignSplicedMateMapLminOverLmate 0 \
	# --alignSplicedMateMapLmin 30 <- STAR default is 0.66 as fraction of read and not a fixed length
	# -> incompatible with arriba (--chimOutType WithinBAM SoftClip)
	# -> no kickstart for fucioncatcher
	# => no kickstart possible for all three tools
	# since we do not limit alignMatesGapMax/alignIntronMax, output may contain commonly observed read-through fusions due to polymerase misses STOP
	# -> apply arriba blacklist filter

	local i odir b e
	declare -a cmd1=()
	for i in "${!_fq1_starfusion[@]}"; do
		helper::basename -f "${_fq1_starfusion[$i]}" -o b -e e
		odir="$outdir/$b"
		mkdir -p $odir

		if [[ ${_fq2_starfusion[$i]} ]]; then
			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cd "$odir"
			CMD
				STAR-Fusion
				--genome_lib_dir "$(dirname "$genome")"
				--CPU $threads
				--left_fq "${_fq1_starfusion[$i]}"
				--right_fq "${_fq2_starfusion[$i]}"
				--output_dir "$odir"
				--FusionInspector validate
				--examine_coding_effect
			CMD
		else
			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cd "$odir"
			CMD
				STAR-Fusion
				--genome_lib_dir "$(dirname "$genome")"
				--CPU $threads
				--left_fq "${_fq1_starfusion[$i]}"
				--output_dir "$odir"
				--FusionInspector validate
				--examine_coding_effect
			CMD
		fi
	done

	$skip && {
		commander::printcmd -a cmd1
	} || {
		{	conda activate starfusion && \
			commander::runcmd -v -b -t 1 -a cmd1 && \
			conda activate py2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

fusions::arriba(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip>     | true/false return
			-s <softskip>     | true/false only print commands
			-5 <skip>         | true/false md5sums, indexing respectively
			-t <threads>      | number of
			-g <genome>       | path to
			-a <gtf>          | path to
			-o <outdir>       | path to
			-p <tmpdir>       | path to
			-f <size>         | assumed mean fragment
			-1 <fastq1>       | array of
			-2 <fastq2>       | array of
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome gtf outdir tmpdir fragmentsize
	declare -n _fq1_arriba _fq2_arriba
	while getopts 'S:s:5:t:g:a:o:p:f:1:2:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((mandatory++)); threads=$OPTARG;;
			g)	((mandatory++)); genome="$OPTARG";;
			a)	((mandatory++)); gtf="$OPTARG";;
			o)	((mandatory++)); outdir="$OPTARG/arriba"; mkdir -p "$outdir" || return 1;;
			p)	((mandatory++)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			f)	((mandatory++)); fragmentsize=$OPTARG;;
			1)	((mandatory++)); _fq1_arriba=$OPTARG;;
			2)	_fq2_arriba=$OPTARG;;
			*)	_usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1
	commander::printinfo "detecting gene fusions arriba"

	# compatible with CTAT_genome_lib
	# Arriba STAR kickoff needs WithinBAM SoftClip or SeparateSAMold instead of Junctions (STAR-Fusion)
	# returning both in one STAR command impossible: related question (Junctions WithinBAM) @ https://github.com/alexdobin/STAR/issues/685
	# --chimMultimapNmax and --peOverlap* options incompatible with WithinBAM or SeparateSAMold
	# Most STAR-based fusion detection tools only consider chimeric alignments as evidence for gene fusions and are blind to focal deletions,
	# hence. As a workaround, these tools recommend reducing the value of the parameter --alignIntronMax. But this impairs the quality of alignment,
	# because it reduces the scope that STAR searches to find a spliced alignment.

	# star array is implicitly defined globally in alignment::star function
	# to be able to use a temporary star array here, declare it explicitly i.e. locally
	declare -a mapper star
	declare -A strandness
	local params=""
	params+=" --outFilterMultimapNmax 1"
	params+=" --outFilterMismatchNmax 3"
	params+=" --chimOutJunctionFormat 1"
	params+=" --chimScoreMin 1"
	params+=" --chimScoreSeparation 1"
	params+=" --chimSegmentMin 10"
	params+=" --chimJunctionOverhangMin 10"
	params+=" --chimNonchimScoreDropMin 10"
	params+=" --alignSJstitchMismatchNmax 5 -1 5 5"
	params+=" --chimScoreDropMax 30"
	params+=" --chimOutType WithinBAM SoftClip"
	params+=" --chimScoreJunctionNonGTAG 0"
	params+=" --chimSegmentReadGapMax 3"

	{	alignment::star \
			-S false \
			-s $skip \
			-5 $skipmd5 \
			-1 _fq1_arriba \
			-2 _fq2_arriba \
			-o "$outdir" \
			-p "$tmpdir" \
			-t $threads \
			-g "$genome" \
			-x $genome.star.idx \
			-r mapper \
			-c "$params" && \

		alignment::postprocess \
			-S false \
			-s $skip \
			-j sort \
			-t $threads \
			-p "$tmpdir" \
			-o "$outdir" \
			-r mapper && \

		alignment::postprocess \
			-S false \
			-s $skip \
			-j index \
			-t $threads \
			-p "$tmpdir" \
			-o "$outdir" \
			-r mapper && \

		alignment::inferstrandness \
			-S false \
			-s $skip \
			-t $threads \
			-r mapper \
			-x strandness \
			-g "$gtf" \
			-p "$tmpd"ir
	} || return 1

	local m f o
	declare -a cmd3
	for f in "${star[@]}"; do
		o="$outdir/$(basename "$f" .bam)"
		commander::makecmd -a cmd3 -s '&&' -c {COMMANDER[0]}<<- CMD
			arriba
				-a "$genome"
				-g "$gtf"
				-b "\$CONDA_PREFIX/var/lib/arriba/blacklist_hg38_GRCh38_2018-11-04.tsv.gz"
				-x "$f"
				-o "$o.fusions.tsv"
				-O "$o.fusions.discarded.tsv"
				-s $(case ${strandness["$f"]} in 0) echo "no";; 1) echo "yes";; 2) echo "reverse";; *) echo "?";; esac)
				-F $fragmentsize
				-T
				-P
		CMD
	done

	$skip && {
		commander::printcmd -a cmd3
	} || {
		{	conda activate arriba && \
			commander::runcmd -v -b -t $threads -a cmd3 && \
			conda activate py2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
