#! /usr/bin/env bash
# (c) Konstantin Riege

function quantify::featurecounts(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-r <mapper>     | array of bams within array of
			-x <strandness> | hash per bam of
			-g <gtf>        | path to
			-l <level>      | feature (default: exon)
			-f <feature>    | feature (default: gene, needs <feature>_id tag)
			-o <outdir>     | path to
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
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
			${FUNCNAME[-2]} usage:
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
	[[ $# -eq 0 ]] && { _usage || return 0; }
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
					vsd = rbind(colnames(dds),assay(varianceStabilizingTransformation(dds, blind=TRUE)));
					ret = apply(vsd, 2, function(x){ write.table(data.frame(rownames(dds), x[2:length(x)], check.names=F), file=paste0(x[1],".vsc"), col.names=F, row.names=F, quote=F, sep="\t") });
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
	# bedg != samtools bedcov <- which is sum(histogram = samtools depth)
	# use
	# bedtools multicov -D -F -bams <bam>[ <bam>..] -bed <bins> > <bedg> (slow)
	# or
	# bamCoverage --bam $f -o <*.bw|*.bedg> -p 24 -bs <binsize> (merges adjacent bins of same coverage into up to 5000000bp bins, ignores mates with an insert >1000bp)
	# or
	# featureCounts from SAF formatted bins, which allows also fractional counts (slow for very small windows)

	# bamCoverage merges bins of equal counts up to a length of 5000000 likewise kent tools/ ucsc bigWigToBedGraph does
	# -> use cgpbigwig: bwcat tool!
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
			assign reads to bins or bed/bed-like file given features and write counts plus TPM/BPM normalized to bigWig files
			- all reads are counted! including:
			- singletons, duplicated, multi-mapped, bad-quality, cross-chromosome, circular and highly mismatched ones
			in case of ATAC data with the intention to quantify mono-nucleosome positions
			- reads will be centered
			- when give, nucleosome free regions will be subtracted
			to maintain orignal/unmerged bins/windows when convertig bigWig to BedGraph, utilize cgpBigWig bwcat

			${FUNCNAME[-2]} usage:
			-S <hardskip>     | true/false return
			-s <softskip>     | true/false only print commands
			-t <threads>      | number of
			-M <maxmemory>    | amount of
			-r <mapper>       | array of bams within array of
			-x <strandness>   | hash per bam of. default: 0
			                    0 - unstranded
			                    1 - stranded
			                    2 - reversely stranded
			-o <outdir>       | path to
			-O <overlap>      | of a read, that spans multiple features, which assigns it to one or multiple feature. default: 1
			                    0   - largest
			                    0-1 - fraction
			                    1-N - basepairs
			-w <windowsize>   | for bins to be quantified, unless given a bed/bed-like file by -b option. default: 100
			                    NOTE: when using -w 1, deepTools bamCoverage will be utilized to create a pileup (see also -e, ignores -x -O -f -u -p)
			-e <fragmentsize> | when using -w 1 or -m, extend paired-end reads to fragment size and single-end or singleton reads to given fragmentsize
			-b <bed>          | or bed-like file with 3 or 6 columns. latter for strand-specific quantification (see also -x, mutually exclusive to -w)
			-i                | first reference base index of -b option given bed-like file is 1 and not 0 (e.g. gtf, gff)
			-f                | count reads fractional if, according to -o option, they can be assigned to multiple features (+1/n)
			                    NOTE: not supported for version < 2.0.4, when combined with -o 0
			-u                | only unambiguously assigend reads will be counted
			-p                | count read-pairs/fragments instead of reads (experimental!)
			                    NOTE1: +1 instead of +2 is counted only if mates are unambiguously assigned to same feature, else +1,+1
			                    NOTE2: fractional quantification is not available
			                    NOTE3: unless -o option corresponds to basepairs, assignment is derived from the sum of both mate lengths
			-m <monoidx>      | array of MNase or mono-nucleosome filtered bam idices within -r
			                    NOTE: deepTools bamCoverage will be utilized and reads will be centered (ignores -x -O -w -b -1 -f -u -p -d)
			-n <nfridx>       | array of nucleosome free region filtered or any other background bam idices within -r (requires also -m)
			                    NOTE: DANPOS will be utilized and reads will be centered (ignores -x -O -w -b -1 -f -u -p -d)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false tmpdir="${TMPDIR:-/tmp}" threads overlap=1 windowsize=100 bed index=0 fractional=false unambiguous=false pairs=false outdir fragmentsize maxmemory
	declare -n _mapper_bam2bedg _strandness_bam2bedg _nidx_bam2bedg _midx_bam2bedg
	while getopts 'S:s:t:M:r:x:o:O:w:b:e:n:m:ifup' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_bam2bedg=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			x)	_strandness_bam2bedg=$OPTARG;;
			o)	outdir="$OPTARG"; mkdir -p "$outdir";;
			O)	overlap=$OPTARG;;
			w)	windowsize=$OPTARG;;
			b)	bed="$OPTARG";;
			i)	index=1;;
			f)	fractional=true;;
			u)	unambiguous=true;;
			p)	pairs=true;;
			e)	fragmentsize=$OPTARG;;
			n)	_nidx_bam2bedg=$OPTARG;;
			m)	_midx_bam2bedg=$OPTARG;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 2 ]] && _usage
	[[ $_nidx_bam2bedg ]] && [[ ! $_midx_bam2bedg ]] && _usage

	if [[ $_midx_bam2bedg ]]; then
		local m i f n params x o odir e="pileup"

		declare -n _bams_bam2bedg=${_mapper_bam2bedg[0]}
		x=$(samtools view -F 4 "${_bams_bam2bedg[0]}" | head -10000 | cat <(samtools view -H "${_bams_bam2bedg[0]}") - | samtools view -c -f 1)
		[[ $x -gt 0 ]] && params="--paired 1" || params="--frsz ${fragmentsize:-200}"
		local minstances mthreads memory=$(samtools view -H "${_bams_bam2bedg[0]}" | sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}' | awk '{i+=$2}END{i=i/100000*3; if(i<1000){print 1000}else{printf "%0.f\n",i}}')
		read -r minstances mthreads < <(configure::instances_by_memory -T $threads -m $memory -M "$maxmemory")

		declare -a tdirs cmd1 cmd2 cmd3 cmd4
		for m in "${_mapper_bam2bedg[@]}"; do
			declare -n _bams_bam2bedg=$m
			for i in "${!_midx_bam2bedg[@]}"; do
				f="${_bams_bam2bedg[${_midx_bam2bedg[$i]}]}"
				o="$(basename "${f%.*}")"
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bam2bedg)")
				odir="$outdir/$m"
				[[ $outdir ]] || odir="$(dirname "$f")"
				mkdir -p "$odir"

				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
					samtools view -H "$f"
				CMD
					sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}'
				CMD
					helper::sort -t $threads -k1,1 > "${tdirs[-1]}/chr.sizes"
				CMD

				if [[ ${_nidx_bam2bedg[$i]} ]]; then
					n="${_bams_bam2bedg[${_nidx_bam2bedg[$i]}]}"
					# region and peak calling after nucleosome position estimation from subtracted data not necessary here
					# for position calling (danpos.py dpos) robert was using -jd 20 i.e. minimum distance of two nucleosome centers before merging into one histone signal
					# should be 147/+147/2. assuming that tn5 can bind nucleosome surface 120/2+120/2 or if conservative 100/2+100/2 which is danpos default..
					#
					# --extend $fragmentsize in case of TF or chip possible. default is 80 (publication mentions 74 as 147/2)
					# whereas bamCoverage --MNase uses only 3! and --center only read length
					#
					# howto save memory or speedup: Split!the!data!by!chromosomes!and!run!a!subset!of!chromosomes!each!time

					mkdir -p "${tdirs[-1]}/counts" "${tdirs[-1]}/tpm"
					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						ln -sn "$(realpath -se "$f")" "${tdirs[-1]}/counts/nucleosome.bam"; ln -sn "$(realpath -se "${f%.*}.bai")" "${tdirs[-1]}/counts/nucleosome.bai"
					CMD
						ln -sn "$(realpath -se "$n")" "${tdirs[-1]}/counts/background.bam"; ln -sn "$(realpath -se "${n%.*}.bai")" "${tdirs[-1]}/counts/background.bai"
					CMD
						danpos.py dtriple
							"${tdirs[-1]}/counts/nucleosome.bam"
							-b "${tdirs[-1]}/counts/background.bam"
							-a 1
							-o "${tdirs[-1]}/counts"
							-n N
							--extend 74
							-r 0
							-k 0
							-j 0
							$params
					CMD
						wigToBigWig "${tdirs[-1]}/counts/pooled/nucleosome.bgsub.smooth.wig" "${tdirs[-1]}/chr.sizes" "$odir/$o.$e.counts.bw"
					CMD

					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						ln -sn "$(realpath -se "$f")" "${tdirs[-1]}/tpm/nucleosome.bam"; ln -sn "$(realpath -se "${f%.*}.bai")" "${tdirs[-1]}/tpm/nucleosome.bai"
					CMD
						ln -sn "$(realpath -se "$n")" "${tdirs[-1]}/tpm/background.bam"; ln -sn "$(realpath -se "${n%.*}.bai")" "${tdirs[-1]}/tpm/background.bai"
					CMD
						danpos.py dtriple
							"${tdirs[-1]}/tpm/nucleosome.bam"
							-b "${tdirs[-1]}/tpm/background.bam"
							-jd 20
							-a 1
							-o "${tdirs[-1]}/tpm"
							-n F
							--extend 74
							-r 0
							-k 0
							-j 0
							$params
					CMD
						wigToBigWig "${tdirs[-1]}/tpm/pooled/nucleosome.bgsub.smooth.wig" "${tdirs[-1]}/chr.sizes" "$odir/$o.$e.tpm.bw"
					CMD

					# commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					# 	wigToBigWig "${tdirs[-1]}/counts/pooled/nucleosome.bgsub.smooth.wig" "${tdirs[-1]}/chr.sizes" "$odir/$o.$e.counts.bw"
					# CMD
					# 	bigWigToBedGraph "$odir/$o.$e.counts.bw" "$odir/$o.$e.counts.bed"
					# CMD
					# commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					# 	wigToBigWig "${tdirs[-1]}/tpm/pooled/nucleosome.bgsub.smooth.wig" "${tdirs[-1]}/chr.sizes" "$odir/$o.$e.tpm.bw"
					# CMD
					# 	bigWigToBedGraph "$odir/$o.$e.tpm.bw" "$odir/$o.$e.tpm.bed"
					# CMD
				else
					# to ensure compatibility with ucsc tools, that require chromosome order to be sorted lexicographical
					# -> better use cgpbigwig and direct bigwig output
					# commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					# 	rm -f "$odir/$o.$e.counts.bed"
					# CMD
					# 	mapfile -t chr < <(cut -f 1 "${tdirs[-1]}/chr.sizes")
					# CMD
					# 	for r in "\${chr[@]}"; do
					# 		bamCoverage
					# 			-r \$r
					# 			-b "$f"
					# 			-bs 1
					# 			-p $threads
					# 			--maxFragmentLength 200000
					# 			--normalizeUsing None
					# 			--outFileFormat bedgraph
					# 			--extendReads ${fragmentsize:-200}
					# 			--centerReads
					# 			-o "${tdirs[-1]}/$e.counts.bed";
					# 		cat "${tdirs[-1]}/$e.counts.bed" >> "$odir/$o.$e.counts.bed";
					# 	done
					# CMD
					commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
						bamCoverage
							-b "$f"
							-bs 1
							-p $threads
							--maxFragmentLength 200000
							--normalizeUsing None
							--extendReads ${fragmentsize:-200}
							--centerReads
							-o "$odir/$o.$e.counts.bw"
					CMD

					commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
						bamCoverage
							-b "$f"
							-bs 1
							-p $threads
							--maxFragmentLength 200000
							--exactScaling
							--normalizeUsing CPM
							--outFileFormat bedgraph
							--extendReads ${fragmentsize:-200}
							--centerReads
							-o "$odir/$o.$e.tpm.bw"
					CMD

					# commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
					# 	bedGraphToBigWig "$odir/$o.$e.counts.bed" "${tdirs[-1]}/chr.sizes" "$odir/$o.$e.counts.bw"
					# CMD
					# commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
					# 	bedGraphToBigWig "$odir/$o.$e.tpm.bed" "${tdirs[-1]}/chr.sizes" "$odir/$o.$e.tpm.bw"
					# CMD
				fi
			done
		done

		if $skip; then
			commander::printcmd -a cmd1
			commander::printcmd -a cmd2
			commander::printcmd -a cmd3
			# commander::printcmd -a cmd4
		else
			commander::runcmd -v -b -i 1 -a cmd1
			commander::runcmd -c danpos -v -b -i $minstances -a cmd2
			commander::runcmd -c deeptools -v -b -i 1 -a cmd3
			# commander::runcmd -c danpos -v -b -i $threads -a cmd4
		fi

		return 0
	fi

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

	local m f params x o e odir
	[[ $bed || $windowsize -gt 1 ]] && e="coverage" || e="pileup"

	declare -a cmd1 cmd2 cmd3 cmd4 tdirs
	for m in "${_mapper_bam2bedg[@]}"; do
		declare -n _bams_bam2bedg=$m
		for f in "${_bams_bam2bedg[@]}"; do
			o="$(basename "${f%.*}")"
			# tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bam2bedg)")
			odir="$outdir/$m"
			[[ $outdir ]] || odir="$(dirname "$f")"
			mkdir -p "$odir"
			if [[ $bed || $windowsize -gt 1 ]]; then
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bam2bedg)")

				if [[ $bed ]]; then
					commander::makecmd -a cmd1 -s ' ' -o "${tdirs[-1]}/regions.saf" -c {COMMANDER[0]}<<- CMD {COMMANDER[2]}<<- 'CMD'
						samtools view -H "$f" | sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}' | helper::sort -t $threads -k1,1 > "${tdirs[-1]}/chr.sizes";
						helper::sort -f "$bed" -t $threads -k1,1 -k2,2n -k3,3n | awk -v i=$index -v OFS='\t'
					CMD
						'BEGIN{i=i==1?0:1}{id=$1":"$2"-"$3; if($4){id=$4}; s="."; if($6){s=$6}; print id,$1,$2+i,$3,".",".",s}'
					CMD
				else
					commander::makecmd -a cmd1 -s '|' -o "${tdirs[-1]}/regions.saf" -c {COMMANDER[0]}<<- CMD {COMMANDER[2]}<<- 'CMD'
						samtools view -H "$f" | sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}' | helper::sort -t $threads -k1,1 > "${tdirs[-1]}/chr.sizes";
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

				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD
					featureCounts
						$params
						-Q 0
						--maxMOp 999
						-s $(declare -p _strandness_bam2bedg | grep -q '=' && echo ${_strandness_bam2bedg["$f"]} || echo 0)
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
					tee -ia >(awk -v OFS='\t' 'NR>2{print \$2,\$3-1,\$4,sprintf("%d",\$7==int(\$7)?\$7:\$7+1)}' | tee -i "$odir/$o.$e.counts.bed" | bg2bw -i /dev/stdin -o "$odir/$o.$e.counts.bw" -c "${tdirs[-1]}/chr.sizes")
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
					tee -i "$odir/$o.$e.tpm.bed" | bg2bw -i /dev/stdin -o "$odir/$o.$e.tpm.bw" -c "${tdirs[-1]}/chr.sizes"
				CMD
			else
				# commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				# 	samtools view -H "$f"
				# CMD
				# 	sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}'
				# CMD
				# 	helper::sort -t $threads -k1,1 > "${tdirs[-1]}/chr.sizes"
				# CMD

				[[ $fragmentsize ]] && params="--extendReads $fragmentsize" || params=''

				# commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				# 	rm -f "$odir/$o.$e.counts.bed"
				# CMD
				# 	mapfile -t chr < <(cut -f 1 "${tdirs[-1]}/chr.sizes")
				# CMD
				# 	for r in "\${chr[@]}"; do
				# 		bamCoverage
				# 			$params
				# 			-r \$r
				# 			-b "$f"
				# 			-bs $windowsize
				# 			-p $threads
				# 			--maxFragmentLength 200000
				# 			--normalizeUsing None
				# 			--outFileFormat bedgraph
				# 			-o "${tdirs[-1]}/$e.counts.bed";
				# 		cat "${tdirs[-1]}/$e.counts.bed" >> "$odir/$o.$e.counts.bed";
				# 	done
				# CMD

				commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
					bamCoverage
						$params
						-b "$f"
						-bs $windowsize
						-p $threads
						--maxFragmentLength 200000
						--normalizeUsing None
						-o "$odir/$o.$e.counts.bw"
				CMD

				commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
					bamCoverage
						$params
						-b "$f"
						-bs $windowsize
						-p $threads
						--maxFragmentLength 200000
						--exactScaling
						--normalizeUsing CPM
						-o  "$odir/$o.$e.tpm.bw"
				CMD
			fi

			# commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
			# 	bedGraphToBigWig "$odir/$o.$e.counts.bed" "${tdirs[-1]}/chr.sizes" "$odir/$o.$e.counts.bw"
			# CMD
			# commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
			# 	bedGraphToBigWig "$odir/$o.$e.tpm.bed" "${tdirs[-1]}/chr.sizes" "$odir/$o.$e.tpm.bw"
			# CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		# commander::printcmd -a cmd4
	else
		commander::runcmd -v -b -i 1 -a cmd1
		commander::runcmd -c subread -v -b -i $instances -a cmd2
		commander::runcmd -c deeptools -v -b -i 1 -a cmd3
		# commander::runcmd -v -b -i $threads -a cmd4
	fi

	return 0
}

function quantify::profiles(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-r <mapper>   | array of bams within array of
			-c <bwdir>    | path to pileup bigwig files
			-g <gtf>      | path to with transcript feature (and transcript_id feature tag) for TSS and transcript profiling
			-b <bed>      | path to
			-o <outdir>   | path to
			-p <pearson>  | correlation plot of coverage files
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}" gtf bed bwdir pearson=false
	declare -n _mapper_profiles _strandness_profiles
	while getopts 'S:s:t:r:g:b:c:o:p' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			r) ((++mandatory)); _mapper_profiles=$OPTARG;;
			g) gtf="$OPTARG";;
			b) bed="$OPTARG";;
			c) ((++mandatory)); bwdir="$OPTARG";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p) pearson=true;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage

	declare -a tdirs cmd1 cmd2 cmd3 coverages pileups toprofile
	local m f b e odir

	for m in "${_mapper_profiles[@]}"; do
		declare -n _bams_profiles=$m
		odir="$outdir/$m"
		mkdir -p "$odir"
		coverages=()
		pileups=()
		toprofile=()
		for f in "${_bams_profiles[@]}"; do
			$pearson && coverages+=("$(find -L "$bwdir/$m" -maxdepth 1 -name "$(basename ${f%.*}).coverage.tpm.bw" -print -quit | grep .)")
			pileups+=("$(find -L "$bwdir/$m" -maxdepth 1 -name "$(basename ${f%.*}).pileup.tpm.bw" -print -quit | grep .)")
		done

		if [[ $gtf ]]; then
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				cat "$gtf"
			CMD
				perl -F'\t' -lane '
					next if $F[0] eq "chrM" || $F[0]=~/^\s*#/;
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
				helper::sort -t $threads -k1,1 -k2,2n -k3,3n > "$odir/transcripts.bed"
			CMD

			toprofile+=("$odir/transcripts.bed")

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.computematrix)")
			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
				n=\$(wc -l < "$odir/transcripts.bed");
				t=$threads;
			CMD
				i=1;
				while [[ $((n/(t*i))) -gt 5000 ]]; do ((i++)); done;
				n=$(perl -M'POSIX qw(ceil)' -se 'print ceil($n/($t*$i))' -- -n=$n -t=$t -i=$i);
			CMD
				cat "$odir/transcripts.bed" | parallel --tmpdir "${tdirs[-1]}" -P $threads -I {} --pipe --cat -N \$n computeMatrix reference-point
					-p 1
					-S $(printf '"%s" ' "${pileups[@]}")
					-R "{}"
					-bs 1
					--beforeRegionStartLength 1000
					--afterRegionStartLength 1000
					--skipZeros
					--missingDataAsZero
					-o "{}.matrix.gz" 2> >(sed -un '/^Skipping/!p' >&2) | cat;
			CMD
				pigz -p 1 -cd "${tdirs[-1]}"/*.matrix.gz | grep -v '^@' | helper::pgzip -t $threads -o "${tdirs[-1]}/data.gz";
				n=\$(gztool -l "${tdirs[-1]}/data.gz" |& sed -nE '/Number of lines/{s/.*:\s+([0-9]+).*/\1/p}');
				pigz -p 1 -cd "${tdirs[-1]}"/*.matrix.gz | head -1 | sed -E 's/"group_boundaries":\[[0-9,]+\]/"group_boundaries":[0,'\$n']/' | pigz -p 1 -kc > "$odir/tss.matrix.gz";
				cat "${tdirs[-1]}/data.gz" >> "$odir/tss.matrix.gz";
			CMD
			# extremely slow and memory hungry: computeMatrixOperations rbind -m "${tdirs[-1]}"/*.matrix.gz -o "$odir/tss.matrix.gz"
			# not necessary when using deeptools plotHeatmap: computeMatrixOperations sort -m "${tdirs[-1]}/tss.matrix.gz" -R "$odir/transcripts.bed" -o "$odir/tss.matrix.gz"
			# use unbuffered sed to avoid "skipping <id> due to being absent in the computeMatrix output" warnings

			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
				plotHeatmap
					-m "$odir/tss.matrix.gz"
					-o "$odir/tss.heatmap.pdf"
					--plotFileFormat pdf
					--colorMap RdBu
					--refPointLabel TSS
					--samplesLabel $(basename -a "${pileups[@]/%.$e.*/}" | xargs -I {} printf "'%s' " {})
			CMD
			# --outFileNameMatrix "$odir/tss.heatmap.matrix.gz"
		fi

		[[ $bed ]] && toprofile+=("$bed")

		for f in "${toprofile[@]}"; do
			b="$(basename "${f%.*}")"

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
			commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
				n=\$(wc -l < "$f");
				t=$threads;
			CMD
				i=1;
				while [[ $((n/(t*i))) -gt 5000 ]]; do ((i++)); done;
				n=$(perl -M'POSIX qw(ceil)' -se 'print ceil($n/($t*$i))' -- -n=$n -t=$t -i=$i);
			CMD
				awk -F '\t' -v OFS='\t' '\$4="region_"NR' "$f" | parallel --tmpdir "${tdirs[-1]}" -P $threads -I {} --pipe --cat -N \$n computeMatrix scale-regions
					-p 1
					-S $(printf '"%s" ' "${pileups[@]}")
					-R "{}"
					-bs 1
					--beforeRegionStartLength 3000
					--regionBodyLength 8000
					--afterRegionStartLength 3000
					--skipZeros
					--missingDataAsZero
					-o "{}.matrix.gz" 2> >(sed -un '/^Skipping/!p' >&2) | cat;
			CMD
				pigz -p 1 -cd "${tdirs[-1]}"/*.matrix.gz | grep -v '^@' | helper::pgzip -t $threads -o "${tdirs[-1]}/data.gz";
				n=\$(gztool -l "${tdirs[-1]}/data.gz" |& sed -nE '/Number of lines/{s/.*:\s+([0-9]+).*/\1/p}');
				pigz -p 1 -cd "${tdirs[-1]}"/*.matrix.gz | head -1 | sed -E 's/"group_boundaries":\[[0-9,]+\]/"group_boundaries":[0,'\$n']/' | pigz -p 1 -kc > "$odir/transcripts.matrix.gz";
				cat "${tdirs[-1]}/data.gz" >> "$odir/transcripts.matrix.gz";
			CMD

			[[ ${#pileups[@]} -gt 1 ]] && e="$(echo -e "${pileups[0]}\t${pileups[1]}" | sed -E 's/(.+)\t.*\1/\t\1/' | cut -f 2)" || e=".pileup.tpm.bw"
			# e=${pileups[0]##*.}

			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
				plotProfile
					-m "$odir/$b.matrix.gz"
					--outFileNameData "$odir/$b.profile.tsv"
					-o "$odir/$b.profile.pdf"
					--plotFileFormat pdf
					--samplesLabel $(basename -a "${pileups[@]/%.$e.*/}" | xargs -I {} printf "'%s' " {})
			CMD
		done

		if [[ ${#coverages[@]} -gt 1 ]]; then
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
				multiBigwigSummary BED-file
					-p $threads
					--bwfiles $(printf '"%s" ' "${coverages[@]}")
					--BED "${coverages[0]%.*}.bed"
					-out "$odir/coverage.npz"
			CMD

			e="$(echo -e "${coverages[0]}\t${coverages[1]}" | sed -E 's/(.+)\t.*\1/\t\1/' | cut -f 2)"
			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
				plotCorrelation
					-in "$odir/coverage.npz"
					--corMethod pearson
					--skipZeros
					--plotTitle "Pearson Correlation of TPM normalized counts"
					--whatToPlot heatmap
					--colorMap RdBu_r
					-o "$odir/coverage.correlation.pearson.pdf"
					--outFileCorMatrix "$odir/coverage.correlation.pearson.tsv"
					--plotFileFormat pdf
					--labels $(basename -a "${coverages[@]/%.$e.*/}" | xargs -I {} printf "'%s' " {})
			CMD
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -c deeptools -v -b -i 1 -a cmd2
		commander::runcmd -c deeptools -v -b -i $threads -a cmd3
		# plotting may heavily consume memory!!
	fi

	return 0
}
