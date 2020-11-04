#! /usr/bin/env bash
# (c) Konstantin Riege

alignment::slice(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <slicesinfo>| hash per bam of
			-p <tmpdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory tmpdir
	declare -n _mapper_slice _bamslices_slice
	while getopts 'S:s:t:m:r:c:p:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			r)	((++mandatory)); _mapper_slice=$OPTARG;;
			c)	((++mandatory)); _bamslices_slice=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "slicing alignments"

	local minstances instances mthreads ithreads m
	read -r minstances mthreads < <(configure::instances_by_memory -t $threads -m $memory)
	for m in "${_mapper_slice[@]}"; do
		declare -n _bams_slice=$m
		((instances+=${#_bams_slice[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -n _bams_slice=${_mapper_slice[0]}
	local tdir f i bed o
	declare -a mapdata cmd1 cmd2 cmd3 cmd4

	mkdir -p "$tmpdir/genome"
	samtools view -H "${_bams_slice[0]}" | sed -En '/^@SQ/{s/.+\tSN:(\S+)\s+LN:(\S+).*/\1\t0\t\2/p}' > $tmpdir/genome/chr.bed
	# do not use process substitution as Rscript argument - sometimes R swallows parts
	commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			suppressMessages(library("knapsack"));
			args <- commandArgs(TRUE);
			slices <- as.numeric(args[1]);
			odir <- args[2];
			df <- read.table(args[3], header=F, sep="\t", stringsAsFactors=F);

			len <- as.integer(df[,3]/10+1);
			srt <- sort(len,decreasing=T,index.return=1);
			len <- srt$x;
			df <- df[srt$ix,];

			lmax <- max(len);
			n <- slices+1;
			bins <- c();
			i <- 1;
			while(n > slices){
				i <- i+1;
				b <- binpacking(len, cap=as.integer(lmax/2*i+1));
				n <- b$nbins;
				bins <- b$xbins;
			};
			for (j in 1:n){
				write.table(df[bins %in% j,], row.names = F, col.names = F, file=file.path(odir,paste("slice",j,"bed",sep=".")), quote=F, sep="\t");
			};
		'
	CMD
		$minstances "$tmpdir/genome" "$tmpdir/genome/chr.bed"
	CMD

	rm -f "$tmpdir/genome/slice".*.bed
	commander::runcmd -c r -v -b -t $threads -a cmd1

	for m in "${_mapper_slice[@]}"; do
		declare -n _bams_slice=$m
		tdir="$tmpdir/$m"
		mkdir -p "$tdir" # do not use mktemp which triggers cleanup, since slices might be reused later
		for f in "${_bams_slice[@]}"; do
			o="$tdir"/$(basename "$f")
			o="${o%.*}"

			_bamslices_slice["$f"]="$o.slices.info"
			ls -v "$tmpdir/genome/slice".*.bed | sed -E "s@.+\.([0-9]+)\.bed@$o.slice.\1.bam@" > "$o.slices.info"

			alignment::_index -1 cmd2 -t $ithreads -i "$f"

			for bed in $(ls -v "$tmpdir/genome/slice".*.bed); do
				i=$(basename "$bed" .bed | rev | cut -d '.' -f 1 | rev)
				commander::makecmd -a cmd3 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					samtools view
						-@ $ithreads
						-L "$bed"
						-M
						-b
						"$f"
						> "$o.slice.$i.bam"
				CMD
					samtools index -@ $ithreads "$o.slice.$i.bam" "$o.slice.$i.bai"
				CMD
			done
		done
	done

	if $skip; then
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -t $instances -a cmd2
		commander::runcmd -v -b -t $instances -a cmd3
	fi

	return 0
}

alignment::rmduplicates(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory tmpdir outdir regex='\S+:(\d+):(\d+):(\d+)\s*.*'
	declare -n _mapper_rmduplicates _bamslices_rmduplicates
	while getopts 'S:s:t:m:x:r:c:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			x)	regex="$OPTARG";;
			r)	((++mandatory)); _mapper_rmduplicates=$OPTARG;;
			c)	((++mandatory)); _bamslices_rmduplicates=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "removing duplicates"

	local minstances mthreads ithreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)
	local nfh=$(($(ulimit -n)/minstances))
	[[ ! $nfh ]] || [[ $nfh -le 1 ]] && nfh=$((1024/minstances))

	local m i o slice instances ithreads odir
	for m in "${_mapper_rmduplicates[@]}"; do
		declare -n _bams_rmduplicates=$m
		((instances+=${#_bams_rmduplicates[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2
	for m in "${_mapper_rmduplicates[@]}"; do
		declare -n _bams_rmduplicates=$m
		odir="$outdir/$m"
		mkdir -p "$odir"
		for i in "${!_bams_rmduplicates[@]}"; do
			tomerge=()
			while read -r slice; do
				commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					picard
						-Xmx${jmem}m
						-XX:ParallelGCThreads=$jgct
						-XX:ConcGCThreads=$jcgct
						-Djava.io.tmpdir="$tmpdir"
						MarkDuplicates
						I="$slice"
						O="$slice.rmdup"
						M="$slice.metrics"
						READ_NAME_REGEX='$regex'
						REMOVE_DUPLICATES=true
						ASSUME_SORT_ORDER=coordinate
						VALIDATION_STRINGENCY=SILENT
						VERBOSITY=WARNING
						MAX_FILE_HANDLES=$nfh
				CMD
					mv "$slice.rmdup" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_rmduplicates[${_bams_rmduplicates[$i]}]}"

			o="$odir"/$(basename "${_bams_rmduplicates[$i]}")
			o="${o%.*}.rmdup.bam"

			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_rmduplicates["$o"]="${_bamslices_rmduplicates[${_bams_rmduplicates[$i]}]}"
			_bams_rmduplicates[$i]="$o"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c picard -v -b -t $minstances -a cmd1
		commander::runcmd -v -b -t $instances -a cmd2
	fi

	return 0
}

alignment::clipmateoverlaps() {
	# must not handle supplementary i.e. circular/chimeric transcripts due to potentially uneven lists
	# if not flagged as supplementary and mate2 is mapped totally before mate1, both will be flagged as unmapped
	# fully covered by mate will be flagged as unmapped
	# adjusted poolsize will consume ~15gb memory
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-o <outbase>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory outdir
	declare -n _mapper_clipmateoverlaps _bamslices_clipmateoverlaps
	declare -A nidx tidx
	while getopts 'S:s:t:m:r:c:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			r)	((++mandatory)); _mapper_clipmateoverlaps=$OPTARG;;
			c)	((++mandatory)); _bamslices_clipmateoverlaps=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "clipping ends of overlapping mate pairs"

	local m i o slice odir instances ithreads minstances mthreads
	read -r minstances mthreads < <(configure::instances_by_memory -t $threads -m $memory)

	for m in "${_mapper_clipmateoverlaps[@]}"; do
		declare -n _bams_clipmateoverlaps=$m
		((instances+=${#_bams_clipmateoverlaps[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3
	for m in "${_mapper_clipmateoverlaps[@]}"; do
		declare -n _bams_clipmateoverlaps=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_bams_clipmateoverlaps[@]}"; do
			tomerge=()

			while read -r slice; do
				# use stdout with bam extension to enfoce bam output
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
					bam clipOverlap
						--in "$slice"
						--out -.bam
						--unmapped
						--excludeFlags 0x80C
						--poolSize 10000000
						--stats
					> "$slice.mateclipped"
				CMD
					[[ $? -le 2 ]] && true || false
				CMD

				commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					mv "$slice.mateclipped" "$slice"
				CMD
					samtools index -@ $ithreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_clipmateoverlaps[${_bams_clipmateoverlaps[$i]}]}"

			o="$odir/$(basename "${_bams_clipmateoverlaps[$i]}")"
			o="${o%.*}.mateclipped.bam"

			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_clipmateoverlaps["$o"]="${_bamslices_clipmateoverlaps[${_bams_clipmateoverlaps[$i]}]}"
			_bams_clipmateoverlaps[$i]="$o"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -c bamutil -v -b -t $minstances -a cmd1
		commander::runcmd -v -b -t $instances -a cmd2
		commander::runcmd -v -b -t $instances -a cmd3
	fi

	return 0
}

alignment::reorder() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory genome tmpdir outdir i
	declare -n _mapper_reorder _bamslices_reorder
	declare -A nidx tidx
	while getopts 'S:s:t:g:m:r:c:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			r)	((++mandatory)); _mapper_reorder=$OPTARG;;
			c)	((++mandatory)); _bamslices_reorder=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage

	commander::printinfo "reordering alignments"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o slice odir instances ithreads
	for m in "${_mapper_reorder[@]}"; do
		declare -n _bams_reorder=$m
		((instances+=${#_bams_reorder[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2
	for m in "${_mapper_reorder[@]}"; do
		declare -n _bams_reorder=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_bams_reorder[@]}"; do
			tomerge=()

			while read -r slice; do
				commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					picard
						-Xmx${jmem}m
						-XX:ParallelGCThreads=$jgct
						-XX:ConcGCThreads=$jcgct
						-Djava.io.tmpdir="$tmpdir"
						ReorderSam
						I="$slice"
						O="$slice.ordered"
						R="$genome"
						VALIDATION_STRINGENCY=SILENT
						VERBOSITY=WARNING
				CMD
					mv "$slice.ordered" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_reorder[${_bams_reorder[$i]}]}"

			o="$odir/$(basename "${_bams_reorder[$i]}")"
			o="${o%.*}.ordered.bam"

			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_reorder["$o"]="${_bamslices_reorder[${_bams_reorder[$i]}]}"
			_bams_reorder[$i]="$o"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c picard -v -b -t $minstances -a cmd1
		commander::runcmd -v -b -t $instances -a cmd2
	fi

	return 0
}

alignment::addreadgroup() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (requieres -2)
			-2 <tumoridx>  | array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory tmpdir outdir i rgprefix=''
	declare -n _mapper_addreadgroup _bamslices_addreadgroup _nidx_addreadgroup _tidx_addreadgroup
	declare -A nidx tidx
	while getopts 'S:s:t:m:n:r:1:2:c:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			n)	rgprefix=$OPTARG;;
			r)	((++mandatory)); _mapper_addreadgroup=$OPTARG;;
			1)	_nidx_addreadgroup=$OPTARG;;
			2)	_tidx_addreadgroup=$OPTARG;;
			c)	((++mandatory)); _bamslices_addreadgroup=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage

	if [[ ! $_nidx_addreadgroup ]]; then
		declare -n _bams_addreadgroup="${_mapper_addreadgroup[0]}"
		for i in "${!_bams_addreadgroup[@]}"; do
			nidx[$i]=1
		done
	else
		[[ ! $_tidx_addreadgroup ]] && _usage
		for i in "${!_nidx_addreadgroup[@]}"; do
			nidx[$i]=1
		done
		for i in "${!_tidx_addreadgroup[@]}"; do
			tidx[$i]=1
		done
	fi

	commander::printinfo "replacing read group tags"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o rgprefix slice instances ithreads odir
	for m in "${_mapper_addreadgroup[@]}"; do
		declare -n _bams_addreadgroup=$m
		((instances+=${#_bams_addreadgroup[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2
	for m in "${_mapper_addreadgroup[@]}"; do
		declare -n _bams_addreadgroup=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_bams_addreadgroup[@]}"; do
			if [[ ! $rgprefix ]]; then
				[[ ${nidx[$i]} ]] && rgprefix=NORMAL || rgprefix=TUMOR
			fi

			tomerge=()
			o="$(basename "${_bams_addreadgroup[$i]}")"
			o="${o%.*}"
			while read -r slice; do
				commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					picard
						-Xmx${jmem}m
						-XX:ParallelGCThreads=$jgct
						-XX:ConcGCThreads=$jcgct
						-Djava.io.tmpdir="$tmpdir"
						AddOrReplaceReadGroups
						I="$slice"
						O="$slice.rg"
						RGID=${rgprefix}id
						RGLB=${rgprefix}lib
						RGPL=illumina
						RGPU=${rgprefix}unit
						RGSM=${rgprefix}
						VALIDATION_STRINGENCY=SILENT
						VERBOSITY=WARNING
				CMD
					mv "$slice.rg" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_addreadgroup[${_bams_addreadgroup[$i]}]}"

			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$odir/$o.rg.bam"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD
			_bamslices_addreadgroup["$odir/$o.rg.bam"]="${_bamslices_addreadgroup[${_bams_addreadgroup[$i]}]}"
			_bams_addreadgroup[$i]="$odir/$o.rg.bam"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c picard -v -b -t $minstances -a cmd1
		commander::runcmd -v -b -t $instances -a cmd2
	fi

	return 0
}

alignment::splitncigar() {
	declare -a tdirs
	_cleanup::alignment::splitncigar(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory genome tmpdir outdir i
	declare -n _mapper_splitncigar _bamslices_splitncigar
	declare -A nidx tidx
	while getopts 'S:s:t:g:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			r)	((++mandatory)); _mapper_splitncigar=$OPTARG;;
			c)	((++mandatory)); _bamslices_splitncigar=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage

	commander::printinfo "splitting N-cigar alignments"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o slice instances ithreads odir
	for m in "${_mapper_splitncigar[@]}"; do
		declare -n _bams_splitncigar=$m
		((instances+=${#_bams_splitncigar[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2
	for m in "${_mapper_splitncigar[@]}"; do
		declare -n _bams_splitncigar=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_bams_splitncigar[@]}"; do
			tomerge=()

			# v3.X ReassignOneMappingQuality: e.g. misused 255 as unique flag to 60
			# -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60
			# v4.X does it automatically for 255, but better do it via STAR parameter outSAMmapqUnique
			# ALLOW_N_CIGAR_READS removed
			# maybe test later: --read-filter NonChimericOriginalAlignmentReadFilter
			# 					--read-filter MappingQualityReadFilter --minimum-mapping-quality 0
			# only for FR-PE data else produces empty file!	--read-filter MateDifferentStrandReadFilter
			while read -r slice; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
				commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					gatk
						--java-options '
							-Xmx${jmem}m
							-XX:ParallelGCThreads=$jgct
							-XX:ConcGCThreads=$jcgct
							-Djava.io.tmpdir="$tmpdir"
						'
						SplitNCigarReads
						-I "$slice"
						-O "$slice.nsplit"
						-R "$genome"
						-verbosity ERROR
						--tmp-dir "${tdirs[-1]}"
				CMD
					mv "$slice.nsplit" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_splitncigar[${_bams_splitncigar[$i]}]}"

			o="$odir/$(basename "${_bams_splitncigar[$i]}")"
			o="${o%.*}.nsplit.bam"

			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_splitncigar["$o"]="${_bamslices_splitncigar[${_bams_splitncigar[$i]}]}"
			_bams_splitncigar[$i]="$o"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c gatk -v -b -t $minstances -a cmd1
		commander::runcmd -v -b -t $instances -a cmd2
	fi

	return 0
}

alignment::leftalign() {
	declare -a tdirs
	_cleanup::alignment::leftalign(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory genome tmpdir outdir i
	declare -n _mapper_leftalign _bamslices_leftalign
	declare -A nidx tidx
	while getopts 'S:s:t:g:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			r)	((++mandatory)); _mapper_leftalign=$OPTARG;;
			c)	((++mandatory)); _bamslices_leftalign=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage

	commander::printinfo "leftaligning alignments"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o slice odir instances ithreads
	for m in "${_mapper_leftalign[@]}"; do
		declare -n _bams_leftalign=$m
		((instances+=${#_bams_leftalign[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2
	for m in "${_mapper_leftalign[@]}"; do
		declare -n _bams_leftalign=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_bams_leftalign[@]}"; do
			tomerge=()

			while read -r slice; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
				commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
							'
						LeftAlignIndels
						-I "$slice"
						-O "$slice.leftaln"
						-R "$genome"
						-verbosity ERROR
						--tmp-dir "${tdirs[-1]}"
				CMD
					mv "$slice.leftaln" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_leftalign[${_bams_leftalign[$i]}]}"

			o="$odir/$(basename "${_bams_leftalign[$i]}")"
			o="${o%.*}.leftaln.bam"

			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_leftalign["$o"]="${_bamslices_leftalign[${_bams_leftalign[$i]}]}"
			_bams_leftalign[$i]="$o"
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c gatk -v -b -t $minstances -a cmd1
		commander::runcmd -v -b -t $instances -a cmd2
	fi

	return 0
}

alignment::bqsr() {
	local tmpfile
	declare -a tdirs
	_cleanup::alignment::bqsr(){
		rm -f "$tmpfile"
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-d <dbsnp>     | path to
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory genome dbsnp tmpdir outdir i
	declare -n _mapper_bqsr _bamslices_bqsr
	declare -A nidx tidx
	while getopts 'S:s:t:g:d:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			d)	dbsnp="$OPTARG";;
			r)	((++mandatory)); _mapper_bqsr=$OPTARG;;
			c)	((++mandatory)); _bamslices_bqsr=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage

	if [[ ! $dbsnp ]]; then
		local tmpfile="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.vcf)"
		dbsnp="$tmpfile"
		echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > "$dbsnp"
		bgzip -f -@ $threads < "$dbsnp" > "$dbsnp.gz"
		tabix -f -p vcf "$dbsnp.gz"
		dbsnp="$dbsnp.gz"
	fi

	commander::printinfo "base quality score recalibration"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o slice odir instances ithreads
	for m in "${_mapper_bqsr[@]}"; do
		declare -n _bams_bqsr=$m
		((instances+=${#_bams_bqsr[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3
	for m in "${_mapper_bqsr[@]}"; do
		declare -n _bams_bqsr=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		for i in "${!_bams_bqsr[@]}"; do
			tomerge=()

			while read -r slice; do
				# https://gatkforums.broadinstitute.org/gatk/discussion/7131/is-indel-realignment-removed-from-gatk4
				# https://software.broadinstitute.org/gatk/blog?id=7847
				#
				# known polymorphic sites used to exclude regions around
				# https://gatkforums.broadinstitute.org/gatk/discussion/1247/what-should-i-use-as-known-variants-sites-for-running-tool-x
				# gatk known-sites bundle: https://software.broadinstitute.org/gatk/documentation/article.php?id=1213
				#
				# BaseRecalibratorSpark with --spark-master local[$mthreads] is BETA!!
				# --bqsr-baq-gap-open-penalty (Phred Scaled) Default value is 40. 30 is perhaps better for whole genome call sets
				# ApplyBQSRSpark fails as of v4.1.2.0
				# -> i.e. also true for BQSRPipelineSpark which does both steps in one
				# GatherBQSRReports - Gathers scattered BQSR recalibration reports into a single file
				#
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
							'
						BaseRecalibrator
						-R "$genome"
						--known-sites "$dbsnp"
						-I "$slice"
						-O "$slice.bqsreport"
						-verbosity ERROR
						--tmp-dir "${tdirs[-1]}"
				CMD

				commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					rm -rf "$slice.bqsr.parts"
				CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
							'
						ApplyBQSR
						-bqsr "$slice.bqsreport"
						-I "$slice"
						-O "$slice.bqsr"
						-verbosity ERROR
						--tmp-dir "${tdirs[-1]}"
				CMD
					mv "$slice.bqsr" "$slice"
				CMD
					samtools index -@ $mthreads "$slice" "${slice%.*}.bai"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_bqsr[${_bams_bqsr[$i]}]}"

			o="$odir/$(basename "${_bams_bqsr[$i]}")"
			o="${o%.*}.bqsr.bam"
			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
				samtools merge
					-@ $ithreads
					-f
					-c
					-p
					"$o"
					$(printf '"%s" ' "${tomerge[@]}")
			CMD

			_bamslices_bqsr["$o"]="${_bamslices_bqsr[${_bams_bqsr[$i]}]}"
			_bams_bqsr[$i]="$o"
		done
	done

	if $skip; then
        commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -c gatk -v -b -t $minstances -a cmd1
		commander::runcmd -c gatk -v -b -t $minstances -a cmd2
		commander::runcmd -v -b -t $instances -a cmd3
	fi

	return 0
}
