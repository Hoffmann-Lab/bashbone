#! /usr/bin/env bash
# (c) Konstantin Riege

function bigwig::_apply(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-c <cmpfiles> | array of
			-i <bwdir>    | path to pileup bigwig files
			-o <outdir>   | path to
			-j <job>      | to be applied [mean|median|sum|stddev|min|max]
			-e            | expect missing values due to differing regions quantified
		EOF
		return 1
	}

	local OPTIND arg mandatory threads tmpdir="${TMPDIR:-/tmp}" skip=false skipmd5=false bwdir outdir missing=false job
	declare -n _mapper_bwapply _cmpfiles_bwapply
	while getopts 'S:s:t:r:c:i:o:j:e' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_bwapply=$OPTARG;;
			c)	((++mandatory)); _cmpfiles_bwapply=$OPTARG;;
			i)	((++mandatory)); bwdir="$OPTARG";;
			j)	((++mandatory)); job="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			e)	missing=true;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "creating per condition mean bigwig tracks"

	declare -a cmd1 tomerge tdirs
	declare -A visited
	local m c odir sample condition library replicate factors instances ithreads
	for m in "${_mapper_bwapply[@]}"; do
		odir="$outdir/$m"
		mkdir -p "$odir"

		for c in $(cat "${_cmpfiles_bwapply[@]}" | perl -F'\t' -lane 'next if exists $m{$F[1]}; $m{$F[1]}=1; print $F[1]'); do
			visited=()
			tomerge=()
			header="id"
			meanheader="id"
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bwapply)")

			while read -r sample condition library replicate factors; do
				[[ ${visited["$sample"]} ]] && continue || visited["$sample"]=1
				tomerge+=("$(find -L "$bwdir/$m" -maxdepth 1 -name "$sample*.pileup.tpm.bw" -print -quit | grep .)")
			done < <(cat "${_cmpfiles_bwapply[@]}" | awk -v c=$c '$2==c' | sort -k4,4V)
			[[ ${#tomerge[@]} -gt $ithreads ]] && ithreads=${#tomerge[@]}

			bigwig::apply \
				-1 cmd1 \
				-t $threads \
				-p ${tdirs[-1]} \
				-f tomerge \
				-o "$odir/$c.$job.pileup.tpm.bw" \
				-j $job \
				-e $missing
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
	else
		$missing && instances=1 || read -r instances ithreads < <(configure::instances_by_threads -t $ithreads -T $threads)
		commander::runcmd -v -b -i $instances -a cmd1
	fi

	return 0
}

function bigwig::apply(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			description
			bigwig::apply returns a bigwig file after applying [mean|median|sum|stddev|min|max] on all intervals of given bigwig or bedgraph files

			synopsis
			bigwig::apply <OPTIONS> [<files>]

			options
			-t <threads>  | optional. number of threads to use (see -e). default: 1
			-f <files>    | optional. array of bigwig files. default: positional arguments
			-o <outfile>  | required. path to output bigwig file
			-j <job>      | required. job to be applied. one of [mean|median|sum|stddev|min|max]
			-e <expect>   | optional. expect and handle missing values due to differing regions quantified. one of [true|false]. default: false
			              | NOTE: if set to false, input files must have suffix .bw
			              | NOTE: if set to false, number of threads == number of input files
			-g <genome>   | optional. supply indexed fasta as switch to enable and operate on bedgraph input

			developer options
			-1 <cmds>     | array to append crafted command to instead of execution
			-p <tmpdir>   | path to unique temporary directory
		EOF
		return 1
	}

	local OPTIND arg mandatory threads=1 tmpdir outfile missing=false job execute=true genome
	declare -n _cmds1_bwapply _files_bwapply
	while getopts '1:t:p:f:j:o:e:g:' arg; do
		case $arg in
			1)	execute=false; _cmds1_bwapply=$OPTARG;;
			t)	threads=$OPTARG;;
			p)	tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			j)	((++mandatory)); job="$OPTARG";;
			f)	_files_bwapply="$OPTARG";;
			o)	((++mandatory)); outfile="$OPTARG"; mkdir -p "$(dirname "$outfile")";;
			e)	missing="$OPTARG";;
			g)	genome="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 2 ]] && _usage

	if $execute; then
		unset _cmds1_bwapply
		declare -a _cmds1_bwapply
		tmpdir=$(mktemp -d -p "${TMPDIR:-/tmp}" cleanup.XXXXXXXXXX.bwapply)
	else
		[[ $tmpdir ]] || _usage
	fi

	shift $((OPTIND-1))
	if [[ $1 ]]; then
		unset _files_bwapply
		declare -a _files_bwapply=("$@")
	else
		[[ $_files_bwapply ]] || _usage
	fi

	if [[ $genome ]]; then
			[[ $job == "stddev" ]] && job="sttdev"

			commander::makecmd -a _cmds1_bwapply -m -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				apply(){
					local JOB_NR=\$(wc -l < "\$FILE")
					paste <(cut -f 1-3 "\$FILE") <(cut -f 4- "\$FILE" | datamash transpose | datamash --narm $job 1-\$JOB_NR | datamash transpose | sed 's/\<nan\>/NA/g')
				}
				export -f apply
			CMD
				bedtools unionbedg -filler NA -i $(printf '"%s" ' "${_files_bwapply[@]}") | helper::lapply -d 20000 -f -t $threads -k -c apply | bg2bw -i /dev/stdin -c "$genome.fai" -o "$outfile"
			CMD
	else
		if $missing; then
			[[ $job == "stddev" ]] && job="sttdev"

			commander::makecmd -a _cmds1_bwapply -s ' ' -m -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
				apply(){
					local JOB_NR=\$(wc -l < "\$FILE")
					paste <(cut -f 1-3 "\$FILE") <(cut -f 4- "\$FILE" | datamash transpose | datamash --narm $job 1-\$JOB_NR | datamash transpose | sed 's/\<nan\>/NA/g')
				}
				export -f apply
				bigWigInfo -chroms "${_files_bwapply[0]}"
			CMD
					| sed -n '/chromCount/,/basesCovered/p'	| sed -nE '2,${$d;s/^\s*(\S+).*\s([0-9]+)$/\1\t\2/p}'
			CMD
					> "$tmpdir/sizes"
				bedtools unionbedg -filler NA -i $(printf '<(bwcat -i "%s") ' "${_files_bwapply[@]}")
			CMD
					| helper::lapply -d 20000 -f -t $threads -k -c apply | bg2bw -i /dev/stdin -c "$tmpdir/sizes" -o "$outfile"
			CMD
		else
			# attention. default values of missing data is 0 and cannot be set to nan
			# wiggletools write_bg - mean default 0 file1.bw default 0 file2.bw
			commander::makecmd -a _cmds1_bwapply -s ' ' -m -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD
				bigWigInfo -chroms "${_files_bwapply[0]}"
			CMD
					| sed -n '/chromCount/,/basesCovered/p' | sed -nE '2,${$d;s/^\s*(\S+).*\s([0-9]+)$/\1\t\2/p}'
			CMD
					> "$tmpdir/sizes"
				while read -r c l;
			CMD
				do
					wiggletools write_bg - seek \$c 1 \$l $job $(printf '"%s" ' "${_files_bwapply[@]}")
				done < "$tmpdir/sizes"
			CMD
					| bg2bw -i /dev/stdin -c "$tmpdir/sizes" -o "$outfile"
			CMD
		fi
	fi

	if $execute; then
		commander::runcmd -i 1 -a _cmds1_bwapply
	fi

	return 0
}

function bigwig::profiles(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-i <bwdir>    | path to pileup bigwig files
			-n <names>    | array of
			-g <gtf>      | path to with transcript feature (and transcript_id feature tag) for TSS and transcript profiling. mutually exclusive to -b
			-b <bedfiles> | array of paths to bed or bed-like files e.g. narrowPeak. mutually exclusive to -g
			-o <outdir>   | path to
			-p            | pearson correlation plot of coverage files
			-a            | restrict gtf based profiling to canonical transcripts (default: all protein-coding and lncRNA transcripts)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}" gtf bed bwdir pearson=false maxmemory canonicals=false
	declare -n _mapper_profiles _strandness_profiles _bedfiles_profiles _names_profiles
	while getopts 'S:s:t:r:n:g:b:i:o:M:pa' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_profiles=$OPTARG;;
			g)	gtf="$OPTARG";;
			b)	_bedfiles_profiles=$OPTARG;;
			n)	_names_profiles=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			i)	((++mandatory)); bwdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p)	pearson=true;;
			a)	canonicals=true;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage
	# [[ ! $_bedfiles_profiles && ! $gtf ]] && _usage

	$pearson && commander::printinfo "correlating bigwigs and creating profiles" || commander::printinfo "creating bigwig profiles"

	declare -p _bedfiles_profiles | grep -q '=' || {
		unset _bedfiles_profiles
		declare -a _bedfiles_profiles
	}

	declare -a tdirs cmd1 cmd2 cmd3 cmd4 coverages pileups toprofile tomerge
	local m f b e odir imemory instances="${#_mapper_profiles[@]}"
	read -r instances imemory < <(configure::memory_by_instances -i $instances -T $threads -M "$maxmemory")
	local ithreads=$((threads/instances))

	for m in "${_mapper_profiles[@]}"; do
		declare -n _bams_profiles=$m
		odir="$outdir/$m"
		mkdir -p "$odir"
		coverages=()
		pileups=()
		toprofile=("${_bedfiles_profiles[@]}")
		tomerge=()

		for f in "${_bams_profiles[@]}"; do
			$pearson && coverages+=("$(find -L "$bwdir/$m" -maxdepth 1 -name "$(basename ${f%.*}).coverage.tpm.bw" -print -quit | grep .)")
			[[ $toprofile || $gtf ]] && pileups+=("$(find -L "$bwdir/$m" -maxdepth 1 -name "$(basename ${f%.*}).pileup.tpm.bw" -print -quit | grep .)")
		done

		#[[ ${#pileups[@]} -gt 1 ]] && e="$(echo -e "${pileups[0]}\t${pileups[1]}" | sed -E 's/(.+)\t.*\1/\t\1/' | cut -f 2)" || e=".pileup.tpm.bw"
		e=".pileup.tpm.bw"

		if [[ $gtf ]]; then
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				{ $canonicals && grep -Fw -f <(canonicals.pl "$gtf" exon gene_id | cut -f 2) $gtf || cat $gtf; }
			CMD
				perl -F'\t' -lane '
					next if $F[0] eq "chrM" || $F[0] eq "MT" || $F[0]=~/^\s*#/;
					next unless $F[2] eq "transcript";
					if ($F[-1]=~/gene_biotype "([^"]+)/){
						$b=$1;
						next if $b=~/(ribozyme|overlapping|non_coding)/;
						if ($b=~/RNA/){next unless $b=~/lncRNA|lincNRA/};
					}
					$t=$F[0].":".($F[3]-1)."-".$F[4];
					next if exists $m{$t};
					$m{$t}=1;
					if ($F[-1]=~/transcript_id "([^"]+)/){
						$t=$1;
					}
					print join"\t",($F[0],$F[3]-1,$F[4],$t,0,$F[6])
				'
			CMD
				helper::sort -t $ithreads -M "$imemory" -k1,1 -k2,2n -k3,3n > "$odir/transcripts.bed"
			CMD

			toprofile=("$odir/transcripts.bed")
		fi

		for f in "${toprofile[@]}"; do
			b="$(basename "${f%.*}")"

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.computematrix)")
			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				cut -f 1-6 $f | helper::lapply
					-t 2
					-o "${tdirs[-1]}"
					-f
					-c '
					TMPDIR="${tdirs[-1]}" computeMatrix reference-point
						-p $(((threads+1)/2))
						-S $(printf '"%s" ' "${pileups[@]}")
						-R "\$FILE"
						-bs 10
						--beforeRegionStartLength 3000
						--afterRegionStartLength 3000
						--missingDataAsZero
						-o "\$FILE.matrix.gz"
				' 2> >(sed -un '/^Skipping/!p' >&2) | cat;
			CMD
				bgzip -@ 1 -kcd "${tdirs[-1]}"/*.matrix.gz | grep -v '^@' | helper::pgzip -t $threads -o "${tdirs[-1]}/data.gz";
				n=\$(gztool -l "${tdirs[-1]}/data.gz" |& sed -nE '/Number of lines/{s/.*:\s+([0-9]+).*/\1/p}');
				bgzip -@ 1 -kcd "${tdirs[-1]}"/*.matrix.gz | head -1 | sed -E 's/"group_boundaries":\[[0-9,]+\]/"group_boundaries":[0,'\$n']/' | bgzip -@ 1 -kc > "$odir/$b.tss.matrix.gz";
				cat "${tdirs[-1]}/data.gz" >> "$odir/$b.tss.matrix.gz";
			CMD
			# either do not use --skipZeros or add || true to avoid error "The region group genes had no matching entries!" i.e. all regions in bed chunk return a zero coverage
			# extremely slow and memory hungry: computeMatrixOperations rbind -m "${tdirs[-1]}"/*.matrix.gz -o "$odir/tss.matrix.gz"
			# not necessary when using deeptools plotHeatmap: computeMatrixOperations sort -m "${tdirs[-1]}/tss.matrix.gz" -R "$odir/transcripts.bed" -o "$odir/tss.matrix.gz"
			# use unbuffered sed to avoid "skipping <id> due to being absent in the computeMatrix output" warnings

			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				TMPDIR="${tdirs[-1]}" plotHeatmap
					-m "$odir/$b.tss.matrix.gz"
					-o "$odir/$b.tss.heatmap.pdf"
					--plotFileFormat pdf
					--colorMap RdBu
					--refPointLabel TSS
					--samplesLabel $([[ $_names_profiles ]] && printf '"%s" ' "${_names_profiles[@]}" || { basename -a "${pileups[@]%$e}" | xargs -I {} printf "'%s' " {}; })
			CMD

			# multithreading capacities limited. parallel instances are much faster and according to
			# https://janbio.home.blog/2021/03/19/speed-up-deeptools-computematrix/
			# runtime sweet spot is at 5000 chunks
			# -> split into threads-fold chunks

			# commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
			# 	computeMatrix scale-regions
			# 		-p $threads
			# 		-S $(printf '"%s" ' "${pileups[@]}")
			# 		-R "$odir/transcripts.bed"
			# 		-bs 1
			# 		--beforeRegionStartLength 3000
			# 		--regionBodyLength 8000
			# 		--afterRegionStartLength 3000
			# 		--skipZeros
			# 		-o "$odir/transcripts.matrix.gz"
			# CMD
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.computematrix)")
			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				awk -F '\t' -v OFS='\t' '\$4="region_"NR' "$f" | cut -f 1-6 $f | helper::lapply
					-t 2
					-o "${tdirs[-1]}"
					-f
					-c '
						TMPDIR="${tdirs[-1]}" computeMatrix scale-regions
						-p $(((threads+1)/2))
						-S $(printf '"%s" ' "${pileups[@]}")
						-R "\$FILE"
						-bs 10
						--beforeRegionStartLength 3000
						--regionBodyLength 8000
						--afterRegionStartLength 3000
						--missingDataAsZero
						-o "\$FILE.matrix.gz"
				' 2> >(sed -un '/^Skipping/!p' >&2) | cat;
			CMD
				bgzip -@ 1 -kcd "${tdirs[-1]}"/*.matrix.gz | grep -v '^@' | helper::pgzip -t $threads -o "${tdirs[-1]}/data.gz";
				n=\$(gztool -l "${tdirs[-1]}/data.gz" |& sed -nE '/Number of lines/{s/.*:\s+([0-9]+).*/\1/p}');
				bgzip -@ 1 -kcd "${tdirs[-1]}"/*.matrix.gz | head -1 | sed -E 's/"group_boundaries":\[[0-9,]+\]/"group_boundaries":[0,'\$n']/' | bgzip -@ 1 -kc > "$odir/$b.matrix.gz";
				cat "${tdirs[-1]}/data.gz" >> "$odir/$b.matrix.gz";
			CMD

			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				TMPDIR="${tdirs[-1]}" plotProfile
					-m "$odir/$b.matrix.gz"
					--outFileNameData "$odir/$b.profile.tsv"
					-o "$odir/$b.profile.pdf"
					--plotFileFormat pdf
					--numPlotsPerRow 2
					--samplesLabel $([[ $_names_profiles ]] && printf '"%s" ' "${_names_profiles[@]}" || { basename -a "${pileups[@]%$e}" | xargs -I {} printf "'%s' " {}; })
			CMD
		done

		if [[ ${#toprofile[@]} -gt 1 ]]; then
			commander::makecmd -a cmd3 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				{	for f in $(printf '"%s" ' "${toprofile[@]}"); do
						b="\$(basename "\${f%.*}")";
						bgzip -@ 1 -cd "$odir/\$b.tss.matrix.gz" | head -1 | paste <(echo \$b) -;
			CMD
					done | perl -F'\t' -lane '
						BEGIN{
							@b=(0);
						}
						push @n,"\"$F[0]\"";
						/group_boundaries":\[\d+,(\d+)\]/;
						push @b,$b[-1]+$1;
						$l=$F[1];
						END{
							$b=join(",",@b);
							$l=~s/"group_boundaries":\[[0-9,]+\]/"group_boundaries":[$b]/;
							$n=join(",",@n);
							$l=~s/"group_labels":\[[^]]+\]/"group_labels":[$n]/;
							print $l;
						}
					';
			CMD
					for f in $(printf '"%s" ' "${toprofile[@]}"); do
						b="\$(basename "\${f%.*}")";
						bgzip -@ 1 -kcd "$odir/\$b.tss.matrix.gz" | grep -v ^@;
					done;
				} | helper::pgzip -t $threads -o "$odir/tss.matrix.gz"
			CMD

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.plot)")
			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				TMPDIR="${tdirs[-1]}" plotHeatmap
					-m "$odir/tss.matrix.gz"
					-o "$odir/tss.heatmap.pdf"
					--plotFileFormat pdf
					--colorMap RdBu
					--refPointLabel TSS
					--samplesLabel $([[ $_names_profiles ]] && printf '"%s" ' "${_names_profiles[@]}" || { basename -a "${pileups[@]%$e}" | xargs -I {} printf "'%s' " {}; })
					--regionsLabel $(basename -a "${toprofile[@]%.*}" | xargs -I {} printf "'%s' " {})
			CMD

			commander::makecmd -a cmd3 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				{	for f in $(printf '"%s" ' "${toprofile[@]}"); do
						b="\$(basename "\${f%.*}")";
						bgzip -@ 1 -cd "$odir/\$b.matrix.gz" | head -1 | paste <(echo \$b) -;
			CMD
					done | perl -F'\t' -lane '
						BEGIN{
							@b=(0);
						}
						push @n,"\"$F[0]\"";
						/group_boundaries":\[\d+,(\d+)\]/;
						push @b,$b[-1]+$1;
						$l=$F[1];
						END{
							$b=join(",",@b);
							$l=~s/"group_boundaries":\[[0-9,]+\]/"group_boundaries":[$b]/;
							$n=join(",",@n);
							$l=~s/"group_labels":\[[^]]+\]/"group_labels":[$n]/;
							print $l;
						}
					';
			CMD
					for f in $(printf '"%s" ' "${toprofile[@]}"); do
						b="\$(basename "\${f%.*}")";
						bgzip -@ 1 -kcd "$odir/\$b.matrix.gz" | grep -v ^@;
					done;
				} | helper::pgzip -t $threads -o "$odir/matrix.gz"
			CMD

			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				TMPDIR="${tdirs[-1]}" plotProfile
					-m "$odir/matrix.gz"
					--outFileNameData "$odir/profile.tsv"
					-o "$odir/profile.pdf"
					--plotFileFormat pdf
					--numPlotsPerRow 2
					--samplesLabel $([[ $_names_profiles ]] && printf '"%s" ' "${_names_profiles[@]}" || { basename -a "${pileups[@]%$e}" | xargs -I {} printf "'%s' " {}; })
					--regionsLabel $(basename -a "${toprofile[@]%.*}" | xargs -I {} printf "'%s' " {})
			CMD
		fi

		if [[ ${#coverages[@]} -gt 1 ]]; then
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.correlation)")

			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bwcat -i "${coverages[0]}" | cut -f 1-3 > "${tdirs[-1]}/regions.bed"
			CMD
				TMPDIR="${tdirs[-1]}" multiBigwigSummary BED-file
					-p $threads
					--bwfiles $(printf '"%s" ' "${coverages[@]}")
					--BED "${tdirs[-1]}/regions.bed"
					-out "$odir/coverage.npz"
			CMD

			#e="$(echo -e "${coverages[0]}\t${coverages[1]}" | sed -E 's/(.+)\t.*\1/\t\1/' | cut -f 2)"
			e=".coverage.tpm.bw"
			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				TMPDIR="${tdirs[-1]}" plotCorrelation
					-in "$odir/coverage.npz"
					--corMethod pearson
					--skipZeros
					--plotTitle "Pearson Correlation of TPM normalized counts"
					--whatToPlot heatmap
					--colorMap RdBu_r
					-o "$odir/coverage.correlation.pearson.pdf"
					--outFileCorMatrix "$odir/coverage.correlation.pearson.tsv"
					--plotFileFormat pdf
					--labels $([[ $_names_profiles ]] && printf '"%s" ' "${_names_profiles[@]}" || { basename -a "${coverages[@]%$e}" | xargs -I {} printf "'%s' " {}; })
			CMD
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		commander::runcmd -v -b -i $instances -a cmd1
		commander::runcmd -c deeptools -v -b -i 1 -a cmd2
		commander::runcmd -v -b -i 1 -a cmd3
		commander::runcmd -c deeptools -v -b -i 1 -a cmd4
		# plotting heavily consumes memory!!
	fi

	return 0
}
