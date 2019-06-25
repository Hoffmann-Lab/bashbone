#! /usr/bin/env bash
# (c) Konstantin Riege

alignment::segemehl() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-5 <skip>       | true/false md5sums, indexing respectively
			-t <threads>    | number of
			-a <accuracy>   | optional: 80 to 100
			-i <insertsize> | optional: 50 to 200000+
			-r <mapper>     | array of bams within array of
			-g <genome>     | path to
			-x <genomeidx>  | path to
			-p <nosplit>    | optional: true/false
			-o <outdir>     | path to
			-1 <fastq1>     | array of
			-2 <fastq2>     | array of
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false skipmd5=false threads genome genomeidx outdir accuracy insertsize nosplitaln=false
	declare -n _fq1_segemehl _fq2_segemehl _segemehl
	while getopts 'S:s:5:t:g:x:a:p:i:r:o:1:2:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((mandatory++)); threads=$OPTARG;;
			g)	((mandatory++)); genome="$OPTARG";;
			x)	((mandatory++)); genomeidx="$OPTARG";;
			o)	((mandatory++)); outdir="$OPTARG/segemehl";;
			a)	accuracy=$OPTARG;;
			p)	nosplitaln=$OPTARG;;
			i)	insertsize=$OPTARG;;
			r)	((mandatory++))
				declare -n _mapper_segemehl=$OPTARG
				_mapper_segemehl+=(segemehl)
				_segemehl=segemehl
			;;
			1) 	((mandatory++)); _fq1_segemehl=$OPTARG;;
			2) 	((mandatory++)); _fq2_segemehl=$OPTARG;;
			*) 	_usage; return 1;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage && return 1

	commander::print "mapping segemehl"

	$skipmd5 && {
		commander::warn "skip checking md5 sums and genome indexing respectively"
	} || {
		commander::print "checking md5 sums"
		thismd5genome=$(md5sum "$genome" 2> /dev/null | cut -d ' ' -f 1)
		[[ -s "$genomeidx" ]] && thismd5segemehl=$(md5sum "$genomeidx" 2> /dev/null | cut -d ' ' -f 1)
		if [[ "$thismd5genome" != "$md5genome" || "$thismd5segemehl" != "$md5segemehl" ]]; then
			commander::print "indexing genome for segemehl"
			declare -a cmdidx
			commander::makecmd -a cmdidx -s '|' -c {COMMANDER[0]}<<- CMD
				segemehl -x "$genomeidx" -d "$genome"
			CMD
			commander::runcmd -v -b -t $threads -a cmdidx || commander::printerr "$funcname failed at indexing"
			thismd5segemehl=$(md5sum "$genomeidx" | cut -d ' ' -f 1)
			sed -i "s/md5segemehl=.*/md5segemehl=$thismd5segemehl/" $genome.md5.sh
		fi
	}

	# read not properly paired - additional tag:
	# YI:i:0 (orientation)
	# YI:i:1 (insertsize)
	# YI:i:2 (orientation + insertsize)
	# YI:i:3 (chimeric)
	# for splice site detection use $o.sngl.bed
	declare -a cmd1
	local o e params
	for i in "${!_fq1_segemehl[@]}"; do
		helper::basename -f "${_fq1_segemehl[$i]}" -o o -e e
		o="$outdir/$o.bam"
        $nosplitaln && params='' || params="-S $o " #segemehl trims suffix
		[[ $accuracy ]] && params+="-A $accuracy "
		if [[ ${_fq2_segemehl[$i]} ]]; then
			[[ $insertsize ]] && params+="-I $insertsize "
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				segemehl
				$params
				-i "$genomeidx"
				-d "$genome"
				-q "${_fq1_segemehl[$i]}"
				-p "${_fq2_segemehl[$i]}"
				-t $threads
				-b
				-o "$o"
			CMD
		else
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				segemehl
				$params
				-i "$genomeidx"
				-d "$genome"
				-q "${_fq1_segemehl[$i]}"
				-t $threads
				-b
				-o "$o"
			CMD
		fi
		_segemehl+=("$o")
	done

	$skip && {
		commander::printcmd -a cmd1 
	} || {
		{	mkdir -p "$outdir" && \
			commander::runcmd -v -b -t 1 -a cmd1 
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

alignment::star() {
	return 1
	# DO USE outSAMmapqUnique 60 instead of 255 - necessary for gatk implemented MappingQualityAvailableReadFilter
	# thuerefore used als q filter during unification - see alignment::uniqify
}

alignment::postprocess() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-j <job>      | [uniqify|sort|index]
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-p <tmpdir>   | path to
			-o <outdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir job
	declare -n _mapper_process
	while getopts 'S:s:t:j:r:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			j) ((mandatory++)); job=${OPTARG,,*};;
			r) ((mandatory++)); _mapper_process=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	local instances ithreads m i outbase tmpbase newbam
	for m in "${_mapper_process[@]}"; do
		declare -n _bams_process=$m
		((instances+=${#_bams_process[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	commander::print "$job alignments"

	declare -a cmd1 cmd2
	for m in "${_mapper_process[@]}"; do
		declare -n _bams_process=$m
		mkdir -p "$outdir/$m" "$tmpdir/$m"
		for i in "${!_bams_process[@]}"; do
			outbase=$(basename "${_bams_process[$i]}")
			outbase="${outbase%.*}"
			tmpbase=$tmpdir/$m/$outbase
			outbase=$outdir/$m/$outbase
			
			case $job in 
				uniqify)
					alignment::_uniqify \
						-1 cmd1 \
						-2 cmd2 \
						-t $ithreads \
						-m $m \
						-i "${_bams_process[$i]}" \
						-o $outbase \
						-r newbam
					_bams_process[$i]="$newbam"
				;;
				sort)
					instances=1
					alignment::_sort \
						-1 cmd1 \
						-t $threads \
						-i "${_bams_process[$i]}" \
						-o $outbase \
						-p $tmpbase \
						-r newbam
					_bams_process[$i]="$newbam"
				;;
				index)
					alignment::_index \
						-1 cmd1 \
						-t $ithreads \
						-i "${_bams_process[$i]}" \
				;;
				*) _usage; return 1;;
			esac
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	mkdir -p "$outdir" && \
			commander::runcmd -v -b -t $instances -a cmd1 && \
			commander::runcmd -v -b -t $instances -a cmd2
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

alignment::_uniqify() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-1 <cmds1>    | array of
			-2 <cmds2>    | array of
			-t <threads>  | number of
			-m <mapper>   | name of
			-i <sam|bam>  | alignment file
			-o <outbase>  | path to
			-r <var>      | returned alignment file
		EOF
		return 0
	}

	local OPTIND arg mandatory threads sambam outbase m
	declare -n _cmds1_uniqify _cmds2_uniqify _returnfile_uniqify
	while getopts '1:2:t:i:o:r:m:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_uniqify=$OPTARG;;
			2) ((mandatory++)); _cmds2_uniqify=$OPTARG;;
			t) ((mandatory++)); threads=$OPTARG;;
			i) ((mandatory++)); sambam="$OPTARG";;
			m) m=${OPTARG,,*};;
			o) ((mandatory++)); outbase="$OPTARG";;
			r) _returnfile_uniqify=$OPTARG; ;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	_returnfile_uniqify="$outbase.unique.bam"

	readlink -e "$sambam" | file -f - | grep -qF compressed || {
		commander::makecmd -a _cmds1_uniqify -s '|' -c {COMMANDER[0]}<<- CMD
			samtools view
				-@ $threads
				-b
				-F 4
				"$sambam"
				> "$outbase.bam"
		CMD
	}

	commander::makecmd -a _cmds2_uniqify -s '|' -c {COMMANDER[0]}<<- CMD
		samtools view
			-@ $threads
			-b
			-f 4
			"$sambam"
			> "$outbase.unmapped.bam"
	CMD

	# infer SE or PE filter
	local params=''
	[[ $(samtools view -F 4 -h "$sambam" | head -10000 | samtools view -c -f 1) -gt 0 ]] && params+='-f 2 '

	if [[ "$m" =~ bwa || "$m" =~ bowtie || "$m" =~ tophat || "$m" =~ hisat || "$m" =~ star || \
		$(samtools view -F 4 -h "$sambam" | head -10000 | grep -cE '\s+NH:i:[0-9]+\s+' ) -eq 0 ]]; then

		#extract uniques just by MAPQ
		[[ "$m" == "star" ]] && params+='-q 60' || params+='-q 1'
		commander::makecmd -a _cmds2_uniqify -s '&&' -c {COMMANDER[0]}<<- CMD
			samtools view
				$params
				-@ $ithreads
				-F 4
				-F 256
				-b
				"$sambam"
				> "$_returnfile_uniqify"
		CMD
	else
		# sed is faster than grep here
		commander::makecmd -a _cmds2_uniqify -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
			LC_ALL=C samtools view
				$params
				-h
				-@ $ithreads
				-F 4 
				-F 256
				"$sambam"
		CMD
			sed -n '/^@/p; /\tNH:i:1\t/p'
		CMD
			samtools view
				-@ $ithreads
				-b
				> "$_returnfile_uniqify"
		CMD
	fi

	return 0
}

alignment::_sort() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-1 <cmds1>    | array of
			-t <threads>  | number of
			-i <bam>      | binary alignment file
			-o <outbase>  | path to
			-p <tmpbase>  | path to
			-r <var>      | returned alignment file
		EOF
		return 0
	}

	local OPTIND arg mandatory threads bam outbase tmpbase
	declare -n _cmds1_sort _returnfile_sort
	while getopts '1:t:i:o:p:r:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_sort=$OPTARG;;
			t) ((mandatory++)); threads=$OPTARG;;
			i) ((mandatory++)); bam="$OPTARG";;
			o) ((mandatory++)); outbase="$OPTARG";;
			p) ((mandatory++)); tmpbase="$OPTARG";;
			r) _returnfile_sort=$OPTARG; ;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	_returnfile_sort="$outbase.sorted.bam"

	commander::makecmd -a _cmds1_sort -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
		rm -f "$tmpbase.*"
	CMD
		samtools sort
			-@ $threads
			-T "$tmpbase"
			"$bam"
			> "$_returnfile_sort"
	CMD

	return 0
}

alignment::_index() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-1 <cmds1>    | array of
			-t <threads>  | number of
			-i <bam>      | sorted alignment file
		EOF
		return 0
	}

	local OPTIND arg mandatory threads bam
	declare -n _cmds1_index
	while getopts '1:t:i:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_index=$OPTARG;;
			t) ((mandatory++)); threads=$OPTARG;;
			i) ((mandatory++)); bam="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage && return 1

	commander::makecmd -a _cmds1_index -s '|' -c {COMMANDER[0]}<<- CMD
		samtools index
			-@ $threads
			"$bam"
			"${bam%.*}.bai"			
	CMD

	return 0
}

alignment::_inferexperiment() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-1 <cmds1>    | array of
			-2 <cmds1>    | array of
			-i <bam>      | alignment file
			-g <gtf>      | annotation file
			-p <tmpdir>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory bam gtf tmpdir
	declare -n _cmds1_inferexperiment _cmds2_inferexperiment
	while getopts '1:2:t:i:g:p:' arg; do
		case $arg in
			1) ((mandatory++)); _cmds1_inferexperiment=$OPTARG;;
			2) ((mandatory++)); _cmds2_inferexperiment=$OPTARG;;
			i) ((mandatory++)); bam="$OPTARG";;
			g) ((mandatory++)); gtf="$OPTARG";;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	commander::makecmd -a _cmds1_inferexperiment -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		perl -lane '
			next if $F[0] =~ /^(MT|chrM)/i;
			next unless $F[2] eq "exon";
			if($F[6] eq "+"){
				$plus++;
				print join"\t",($F[0],$F[3],$F[4],$F[1],0,$F[6]);
			} else {
				$minus++;
				print join"\t",($F[0],$F[3],$F[4],$F[1],0,$F[6]);
			}
			exit if $plus>2000 && $minus>2000;
		'
	CMD
		"$gtf" > "$tmpdir/$(basename $bam).bed"
	CMD

	# 0 - unstranded
	# 1 - dUTP ++,-+ (FR stranded)
	# 2 - dUTP +-,++ (FR, reversely stranded)
	commander::makecmd -a _cmds2_inferexperiment -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
		{ echo "$bam" &&
			infer_experiment.py
			-q 0 
			-i "$bam"
			-r "$tmpdir/$(basename $bam).bed"; 
		}
	CMD
		perl -lane '
			BEGIN{
				$p=1;
			}
			$f=$_ if $.==1;
			$p=0 if /SingleEnd/i;
			$s1=$F[-1] if $F[-2]=~/^\"\d?\+\+/;
			$s2=$F[-1] if $F[-2]=~/^\"\d?\+-/;
			END{
				print $s1 > 0.7 ? "1 $f" : $s2 > 0.7 ? "2 $f" : "0 $f";
			}
		'
	CMD

	return 0
}

alignment::slice(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>  | true/false return
			-s <softskip>  | truefalse only print commands
			-t <threads>   | number of
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory tmpdir
	declare -n _mapper_slice _bamslices_slice
	while getopts 'S:s:t:m:r:c:p:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			r) ((mandatory++)); _mapper_slice=$OPTARG;;
			c) ((mandatory++)); _bamslices_slice=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	commander::print "slicing alignments"

	local minstances instances mthreads ithreads m
	read -r minstances mthreads < <(configure::instances_by_memory -t $threads -m $memory)
	for m in "${_mapper_slice[@]}"; do
		declare -n _bams_slice=$m
		((instances+=${#_bams_slice[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	local tdir f i chrs o
	declare -a cmd1 cmd2 cmd3
	for m in "${_mapper_slice[@]}"; do
		declare -n _bams_slice=$m
		tdir="$tmpdir/$m"
		mkdir -p "$tdir"

		for f in "${_bams_slice[@]}"; do
			cmd1=()
			commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
				Rscript - <<< '
					suppressMessages(library("knapsack"));
					args <- commandArgs(TRUE);
					slices <- as.numeric(args[1]);
					args <- args[-1];
					chr <- args[rep(c(T,F),length(args)/2)];
					len <- as.numeric(args[rep(c(F,T),length(args)/2)]);
					len <- as.integer(len/10);

					srt <- sort(len,decreasing=T,index.return=1);
					len <- srt$x;
					chr <- chr[srt$ix];

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
						cat(paste(chr[bins %in% j],sep="",collapse=" "), sep="\n");
					};
				'
			CMD
				$minstances $(samtools view -H "$f" | sed -rn '/^@SQ/{s/.+\tSN:(\S+)\s+LN:(\S+).*/\1\t\2/p}')
			CMD

			o="$tdir"/$(basename "$f")
			o="${o%.*}"
			rm -f "$o.slices.info"
			_bamslices_slice["$f"]="$o.slices.info"

			alignment::_index -1 cmd2 -t $ithreads -i "$f"

			i=0
			mapfile -t < <(conda activate py2r &>/dev/null && commander::runcmd -a cmd1)
			for chrs in "${MAPFILE[@]}"; do
				((++i))
				echo "$o.slice$i.bam" >> "$o.slices.info"
				# -M for multi reagion iterator is faster and keeps sort order
				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
					samtools view 
						-@ $ithreads
						-M
						-b
						"$f"
						$chrs > "$o.slice$i.bam"
				CMD
			done
		done
	done

	$skip && {
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	} || {
		{	commander::runcmd -v -b -t $instances -a cmd2 && \
			commander::runcmd -v -b -t $instances -a cmd3
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

alignment::rmduplicates(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>  | true/false return
			-s <softskip>  | truefalse only print commands
			-t <threads>   | number of
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory tmpdir outdir regex='\S+:(\d+):(\d+):(\d+)\s*.*'
	declare -n _mapper_rmduplicates _bamslices_rmduplicates
	while getopts 'S:s:t:m:x:r:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			x) regex="$OPTARG";;
			r) ((mandatory++)); _mapper_rmduplicates=$OPTARG;;
			c) ((mandatory++)); _bamslices_rmduplicates=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage && return 1

	commander::print "removing duplicates"

	local minstances mthreads ithreads jmem jgct jcgct 
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)
	local nfh=$(($(ulimit -n)/minstances))
	[[ ! $nfh ]] || [[ $nfh -le 1 ]] && nfh=$((1024/minstances))

	local m i slice instances ithreads odir
	for m in "${_mapper_rmduplicates[@]}"; do
		declare -n _bams_rmduplicates=$m
		((instances+=${#_bams_rmduplicates[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3
	for m in "${_mapper_rmduplicates[@]}"; do
		declare -n _bams_rmduplicates=$m
		odir="$outdir/$m"
		mkdir -p "$odir"
		for i in "${!_bams_rmduplicates[@]}"; do
			tomerge=()
			while read -r slice; do
				alignment::_index -1 cmd1 -t $ithreads -i "$slice"
				
				commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
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

				tomerge+=("$slice")
			done < "${_bamslices_rmduplicates[${_bams_rmduplicates[$i]}]}"

			o="$odir"/$(basename "${_bams_rmduplicates[$i]}")
			o="${o%.*}.rmdup.bam"

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

			_bamslices_rmduplicates["$o"]="${_bamslices_rmduplicates[${_bams_rmduplicates[$i]}]}"
			_bams_rmduplicates[$i]="$o"
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	} || {
		{	commander::runcmd -v -b -t $instances -a cmd1 && \
			commander::runcmd -v -b -t $minstances -a cmd2 && \
			commander::runcmd -v -b -t $instances -a cmd3
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

alignment::reorder() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::printerr {COMMANDER[0]}<<- EOF
			$funcname usage: 
			-S <hardskip>  | true/false return
			-s <softskip>  | truefalse only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory genome tmpdir outdir i
	declare -n _mapper_reorder _bamslices_reorder
	declare -A nidx tidx
	while getopts 'S:s:t:g:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			r) ((mandatory++)); _mapper_reorder=$OPTARG;;
			c) ((mandatory++)); _bamslices_reorder=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG";;
			o) ((mandatory++)); outdir="$OPTARG";;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1

	commander::print "reordering alignments"

	local minstances mthreads jmem jgct jcgct 
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i slice odir tdir instances ithreads dinstances djmem djgct djcgct
	for m in "${_mapper_reorder[@]}"; do
		declare -n _bams_reorder=$m
		((instances+=${#_bams_reorder[@]}))
	done
	read -r dinstances ithreads djmem djgct djcgct < <(configure::jvm -i $instances -t 1 -T $threads)
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3 cmd4
	for m in "${_mapper_reorder[@]}"; do
		declare -n _bams_reorder=$m
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		for i in "${!_bams_reorder[@]}"; do
			tomerge=()
			o="$(basename "${_bams_reorder[$i]}")"
			o=${o%.*}

			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
				ln -sfn "$genome" "$tdir/$o.fa"
			CMD
				samtools faidx "$tdir/$o.fa"
			CMD
				rm -f "$tdir/$o.dict"
			CMD
				picard
					-Xmx${djmem}m
					-XX:ParallelGCThreads=$djgct
					-XX:ConcGCThreads=$djcgct
					-Djava.io.tmpdir="$tmpdir"
					CreateSequenceDictionary
					R="$tdir/$o.fa"
					VERBOSITY=WARNING
			CMD

			while read -r slice; do
				alignment::_index -1 cmd2 -t $ithreads -i "$slice"

				commander::makecmd -a cmd3 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				picard
					-Xmx${jmem}m
					-XX:ParallelGCThreads=$jgct
					-XX:ConcGCThreads=$jcgct
					-Djava.io.tmpdir="$tmpdir"
					ReorderSam
					I="$slice"
					O="$slice.ordered"
					R="$tdir/$o.fa"
					VALIDATION_STRINGENCY=SILENT
					VERBOSITY=WARNING
				CMD
					mv "$slice.ordered" "$slice"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_reorder[${_bams_reorder[$i]}]}"

			o="$odir/$o.ordered.bam"
			# slices have full sam header info used by merge to maintain the global sort order
			commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD
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

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	} || {
		{	commander::runcmd -v -b -t $dinstances -a cmd1 && \
			commander::runcmd -v -b -t $instances -a cmd2 && \
			commander::runcmd -v -b -t $minstances -a cmd3 && \
			commander::runcmd -v -b -t $instances -a cmd4
		} || { 
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}
