#! /usr/bin/env bash
# (c) Konstantin Riege

function quantify::salmon(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} quantifies either transcriptomic alignments or pre-processed read data utilizing the quasi-mapping software salmon

			-S <hardskip>     | optional. default: false
			                  | [true|false] do nothing and return
			-s <softskip>     | optional. default: false
			                  | [true|false] do nothing but check for files and print commands
			-5 <skip>         | optional, default: false
			                  | [true|false] skip md5sum check, indexing respectively
			-t <threads>      | mandatory
			                  | number of threads
			-M <maxmemory>    | optional. default: all available
			                  | amount of memory to allocate
			-r <mapper>       | mandatory
			                  | array of array names which contain alignment paths. salmon will be added unless -1 option
			                  | mapper+=(salmon); salmon=(/outdir/salmon/1.counts /outdir/salmon/2.counts ..)
			-g <reference>    | mandatory
			                  | path to reference in fasta format
			-a <gtf>          | mandatory unless reference is transcriptomic
			                  | path to annotation in gtf format
			-f <feature>      | optional. default: gene
			                  | feature for which counts should be summed up. gtf needs <feature>_id tag
			-i <transcripts>  | optional. default: false
			                  | [true|false] reference is transcriptomic
			                  |   true:  transcript extraction prior to indexing not necessary
			                  |   false: requires gtf
			-x <genomeidxdir> | mandatory
			                  | path to salmon or salmonTE genome index
			-o <outdir>       | mandatory
			                  | path to output directory. subdirectory salmon will be created
			-1 <fastq1>       | optional.
			                  | array which contains single or first mate fastq(.gz) paths
			                  | unless given, count name sorted, transcriptomic alignments from array of array names which contain alignment paths (see -r)
			                  | alignments must have been processed by alignment::postprocess -j namesort
			-2 <fastq2>       | optional
			                  | array which contains mate pair fastq(.gz) paths
			-d <strandness>   | optional. default when used with -s: ?
			                  | [0|1|2] to define library strandness method to skip auto inference. 0 = unstranded, 1 = stranded/fr second strand or 2 = reversely stranded /fr first strand
			-F                | optional
			                  | force indexing even if md5sums match. ignored upon -5
			-P <parameter>    | optional
			                  | additional salmonTE parameter
		EOF
		return 1
	}

	# default regex for ensembl trascript fasta
	local OPTIND arg mandatory skip=false skipmd5=false threads genome gtf genomeidxdir outdir forceidx=false countaln=false inparams feature=gene tmpdir="${TMPDIR:-/tmp}" transcriptome=false maxmemory default
	declare -n _fq1_salmon _fq2_salmon _mapper_salmon
	declare -g -a salmon=()
	while getopts 'S:s:5:t:g:a:x:d:o:1:2:f:r:M:P:i:Fh' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			5)	$OPTARG && skipmd5=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_salmon=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			a)	gtf="$OPTARG";;
			x)	((++mandatory)); genomeidxdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			1)	_fq1_salmon=$OPTARG;;
			2)	_fq2_salmon=$OPTARG;;
			f)	feature="$OPTARG";;
			i)	transcriptome="$OPTARG";;
			M)	maxmemory=$OPTARG;;
			P)	inparams="$OPTARG";;
			F)	forceidx=true;;
			d)	default=$OPTARG;;
			h)	{ _usage || return 0; };;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage
	$transcriptome || [[ $gtf ]] || _usage

	commander::printinfo "quantifying reads by salmon"

	if $skipmd5; then
		commander::warn "skip checking md5 sums and genome indexing respectively"
	else
		commander::printinfo "checking md5 sums"
		[[ ! -s "$genome.md5.sh" ]] && cp "$BASHBONE_DIR/lib/md5.sh" "$genome.md5.sh"
		source "$genome.md5.sh"

		local doindex=$forceidx idxparams thismd5genome thismd5salmon thismd5gtf
		thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
		[[ -s "$genomeidxdir/seq.bin" ]] && thismd5salmon=$(md5sum "$genomeidxdir/seq.bin" | cut -d ' ' -f 1)
		if ! $transcriptome; then
			thismd5gtf=$(md5sum "$gtf" | cut -d ' ' -f 1) && [[ "$thismd5gtf" != "$md5gtf" ]] && doindex=true
		fi
		if [[ "$thismd5genome" != "$md5genome" || ! "$thismd5salmon" || "$thismd5salmon" != "$md5salmon" ]]; then
			doindex=true
		fi

		if $doindex; then
			mkdir -p "$genomeidxdir"
			declare -a cmdidx

			if $transcriptome; then
				commander::printinfo "indexing reference for salmon"

				commander::makecmd -a cmdidx -s ';' -c {COMMANDER[0]}<<- CMD
					salmon index -t "$genome" -i "$genomeidxdir" -p $threads
				CMD
			else
				commander::printinfo "extracting and indexing transcriptome for salmon"

				# use much faster decoy index
				# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
				# and create clades to satisfy salmonTE
				commander::makecmd -a cmdidx -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					genome2transcriptome.pl -v -f "$genome" -g "$gtf" -o "$genomeidxdir";
				CMD
					perl -F'\t' -slane '
						BEGIN{print "name,class,clade"}
						next unless $F[2] eq "gene";
						$F[-1]=~/${f}_id "([^"]+)/;
						print "$F[0],$1,$f"
					'
				CMD
					-- -f=$feature "$genomeidxdir/transcriptome.gtf" > "$genomeidxdir/clades.csv";
				CMD
					samtools faidx --fai-idx /dev/stdout "$genome" | cut -f 1 > "$genomeidxdir/decoys.txt";
					cat "$genomeidxdir/transcriptome.gtf" "$gtf" > "$genomeidxdir/decoygenome.gtf";
					cat "$genomeidxdir/transcriptome.fa" "$genome" > "$genomeidxdir/decoygenome.fa";
					salmon index -t "$genomeidxdir/decoygenome.fa" -d "$genomeidxdir/decoys.txt" -i "$genomeidxdir" -p $threads
				CMD
			fi

			commander::runcmd -c salmon -v -b -i 1 -a cmdidx
			commander::printinfo "updating md5 sums"
			thismd5salmon=$(md5sum "$genomeidxdir/seq.bin" | cut -d ' ' -f 1)
			sed -i "s/md5salmon=.*/md5salmon=$thismd5salmon/" "$genome.md5.sh"
			sed -i "s/md5genome=.*/md5genome=$thismd5genome/" "$genome.md5.sh"
			$transcriptome || sed -i "s/md5gtf=.*/md5gtf=$thismd5gtf/" "$genome.md5.sh"
		fi
	fi

	declare -a tdirs cmd1 cmd2
	local o b e i m f x imemory instances ithreads params

	if [[ $_fq1_salmon ]]; then
		outdir="$outdir/salmon"
		mkdir -p "$outdir"
		_mapper_salmon+=(salmon)
		cmd1=()
		cmd2=()

		instances=${#_fq1_salmon[@]}
		read -r instances imemory < <(configure::memory_by_instances -i $instances -T $threads -M "$maxmemory")
		ithreads=$((threads/instances))

		params="$inparams"
		if [[ -e "$(dirname "$(which SalmonTE.py)")/reference/$(basename "$genomeidxdir")" ]]; then
			params+=" --reference='$(basename "$genomeidxdir")'"
		else
			tdirs+=("$(mktemp -d -p "$(dirname "$(which SalmonTE.py)")/reference/" "$(basename "$genomeidxdir").XXXXXXXXXX")")
			ln -sfn "$(realpath -s "$genomeidxdir")/"* "${tdirs[0]}"
			params+=" --reference='$(basename "${tdirs[0]}")'"
		fi

		for i in "${!_fq1_salmon[@]}"; do
			helper::basename -f "${_fq1_salmon[$i]}" -o b -e e
			o="$(realpath -s "$outdir/$b.${feature}counts")"
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.salmon)")

			if [[ ${_fq2_salmon[$i]} ]]; then
				# --useVBOpt is default in recent version, installed by bashbone compile
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					cd "${tdirs[-1]}";
					ln -s "$(realpath -se "${_fq1_salmon[$i]}")" sample_R1.fastq.gz;
					ln -s "$(realpath -se "${_fq2_salmon[$i]}")" sample_R2.fastq.gz
				CMD
					SalmonTE.py quant
						$params
						--outpath="${tdirs[-1]}"
						--num_threads=$threads
						--exprtype=count
						sample_R1.fastq.gz
						sample_R2.fastq.gz
				CMD
					sed "s/sample/$b/" MAPPING_INFO.csv
				CMD
					cat "sample/quant.sf" | tee >(awk -F '\\t' 'NR>1{printf "%s\\t%0.f\n",\$1,\$NF}' | helper::sort -t $threads -M "$maxmemory" -k1,1 > "$o.htsc") > "$o"
				CMD
				# sample/quant.sf is default salmon output. ./EXPR.csv is a joined csv count matrix from salmonTE in case of multiple inputs...
				# if salmonTE reference is used as input, salmonTE can be further used to test for differential TE expression on consensus family level
				# (..which is an outdated method. see: 10.1093/bib/bbab417)
				# therefore, join all jounts (similar to experiments.htsc) into a file called EXPR.csv with header TE,sample1,sample2,..
			else
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					cd "${tdirs[-1]}";
					ln -s "$(realpath -se "${_fq1_salmon[$i]}")" sample.fastq.gz
				CMD
					SalmonTE.py quant
						$params
						--outpath="${tdirs[-1]}"
						--num_threads=$threads
						--exprtype=count
						sample.fastq.gz
				CMD
					sed "s/sample/$b/" MAPPING_INFO.csv
				CMD
					cat "sample/quant.sf" | tee >(awk -F '\\t' 'NR>1{printf "%s\\\t%0.f\n",\$1,\$NF}' | helper::sort -t $threads -M "$maxmemory" -k1,1 > "$o.htsc") > "$o"
				CMD
			fi

			salmon+=("$o")

			if [[ $gtf ]]; then
				commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
					perl -F'\t' -slane '
						BEGIN{
							open F,"<$g" or die $!;
							while(<F>){
								chomp;
								@F=split/\t/;
								next if exists $m{$F[0]};
								if ($F[-1]=~/${f}_id "([^"]+)/){
									$m{$F[0]}=$1;
									$c{$1}=0;
								}
							}
							close F;
						}
						next if $.==1;
						$c{$m{$F[0]}}+=$F[-1];
						END{
							printf "%s\t%0.f\n",$_,$c{$_} for keys %c;
						}
					'
				CMD
					-- -g="$($transcriptome && echo "$gtf" || echo "$genomeidxdir/transcriptome.gtf")" -f="$feature" "$o" | helper::sort -t $ithreads -M $imemory -k1,1 > "$o.htsc"
				CMD
			fi
		done

		if $skip; then
			commander::printcmd -a cmd1
			commander::printcmd -a cmd2
		else
			commander::runcmd -c salmon -v -b -i 1 -a cmd1
			commander::runcmd -v -b -i $threads -a cmd2
		fi
	elif [[ ${_mapper_salmon[0]} ]]; then

		declare -A strandness

		if [[ $gtf ]]; then
			alignment::inferstrandness \
				-S false \
				-s $skip \
				-d "$default" \
				-t $threads \
				-r _mapper_salmon \
				-x strandness \
				-g "$($transcriptome && echo "$gtf" || echo "$genomeidxdir/transcriptgenome.gtf")"
		fi

		declare -n _bams_salmon="${_mapper_salmon[0]}"
		local instances=$((${#_mapper_salmon[@]}*${#_bams_salmon[@]}))
		read -r instances imemory < <(configure::memory_by_instances -i $instances -T $threads -M "$maxmemory")
		ithreads=$((threads/instances))

		for m in "${_mapper_salmon[@]}"; do
			declare -n _bams_salmon=$m
			mkdir -p "$outdir/$m"
			for f in "${_bams_salmon[@]}"; do
				o="$outdir/$m/$(basename "$f" .bam).${feature}counts"
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.salmon)")

				x=$(samtools view -F 4 "$f" | head -10000 | cat <(samtools view -H "$f") - | samtools view -c -f 1)
				[[ $x -gt 0 ]] && params="$inparams -l I" || params="$inparams -l "

				case ${strandness["$f"]} in
					0) params+="U";;
					1) params+="SF";;
					2) params+="SR";;
					*) params="$inparams -l A";;
				esac

				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					salmon quant
						$params
						--discardOrphans
						--minAssignedFrags 1
						-a "$f"
						-t "$($transcriptome && echo "$genome" || echo "$genomeidxdir/transcriptome.fa")"
						--threads $threads
						-o "${tdirs[-1]}"
				CMD
					cat "${tdirs[-1]}/quant.sf" | tee >(awk -F '\\t' 'NR>1{printf "%s\\\t%0.f\n",\$1,\$NF}' | helper::sort -t $threads -M "$maxmemory" -k1,1 > "$o.htsc") > "$o"
				CMD

				if [[ $gtf ]]; then
					commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
						perl -F'\t' -slane '
							BEGIN{
								open F,"<$g" or die $!;
								while(<F>){
									chomp;
									@F=split/\t/;
									next if exists $m{$F[0]};
									if ($F[-1]=~/${f}_id "([^"]+)/){
										$m{$F[0]}=$1;
										$c{$1}=0;
									}
								}
								close F;
							}
							next if $.==1;
							$c{$m{$F[0]}}+=$F[-1];
							END{
								printf "%s\t%0.f\n",$_,$c{$_} for keys %c;
							}
						'
					CMD
						-- -g="$($transcriptome && echo "$gtf" || echo "$genomeidxdir/transcriptome.gtf")" -f="$feature" "$o" | helper::sort -t $ithreads -M $imemory -k1,1 > "$o.htsc"
					CMD
				fi
			done
		done

		if $skip; then
			commander::printcmd -a cmd1
			commander::printcmd -a cmd2
		else
			commander::runcmd -c salmon -v -b -i 1 -a cmd1
			commander::runcmd -v -b -i $threads -a cmd2
		fi
	fi

	return 0
}

function quantify::featurecounts(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>   | true/false return
			-s <softskip>   | true/false only print commands
			-t <threads>    | number of
			-M <maxmemory>  | amount of
			-r <mapper>     | array of bams within array of
			-x <strandness> | hash per bam of
			-g <gtf>        | path to
			-l <level>      | feature (default: exon)
			-f <feature>    | feature (default: gene, needs <feature>_id tag)
			-i <transcripts>| true/false input is transcriptomic i.e. count fractional by transcript_id tag in gtf and sum up per feature
			-o <outdir>     | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}" gtf level="exon" feature="gene" transcriptome=false maxmemory
	declare -n _mapper_featurecounts _strandness_featurecounts
	while getopts 'S:s:t:r:x:g:l:f:o:i:M:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_featurecounts=$OPTARG;;
			x)	((++mandatory)); _strandness_featurecounts=$OPTARG;;
			g)	((++mandatory)); gtf="$OPTARG";;
			l)	level=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			f)	feature=$OPTARG;;
			i)	transcriptome=$OPTARG;;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "quantifying reads by featureCounts"

	# featurecounts cannot handle more than 64 threads
	declare -n _bams_featurecounts="${_mapper_featurecounts[0]}"
	local instances ithreads imemory i=$((${#_mapper_featurecounts[@]}*${#_bams_featurecounts[@]}))
	read -r instances imemory < <(configure::memory_by_instances -i $i -T $threads -M "$maxmemory")
	local ithreads2=$((threads/instances))
	read -r instances ithreads < <(configure::instances_by_threads -i $i -t 64 -T $threads)
	[[ $instances -lt $i && $ithreads -gt 64 ]] && ((++instances)) && ithreads=$((threads/instances))
	[[ $ithreads -gt 64 ]] && ithreads=64

	declare -a cmdchk=("featureCounts -v |& grep -oE 'v[.0-9]+'")
	local version=$(commander::runcmd -c subread -a cmdchk)
	[[ "$(echo -e "v2.0.1\n$version" | sort -Vr | head -1)" == "v2.0.1" ]] && version=" " || version="--countReadPairs "

	declare -a tdirs cmd1 cmd2
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
			[[ $x -gt 0 ]] && params+="-p $version "
			if $transcriptome; then
				# mimics salmon fraction count
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					featureCounts
						$params
						-Q 0
						--minOverlap 10
						--maxMOp 999
						-s ${_strandness_featurecounts["$f"]}
						-T $ithreads
						-t $level
						-g transcript_id
						--fraction
						-M
						--ignoreDup
						--tmpDir "${tdirs[-1]}"
						-a "$gtf"
						-o "$o"
						"$f"
				CMD

				commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
					perl -F'\t' -slane '
						BEGIN{
							open F,"<$g" or die $!;
							while(<F>){
								chomp;
								@F=split/\t/;
								next if exists $m{$F[0]};
								if ($F[-1]=~/${f}_id "([^"]+)/){
									$m{$F[0]}=$1;
								}
							}
							close F;
						}
						next if $.<3;
						$c{$m{$F[0]}}+=$F[-1];
						END{
							printf "%s\t%0.f\n",$_,$c{$_} for keys %c;
						}
					'
				CMD
					-- -g="$gtf" -f="$feature" "$o" | helper::sort -t $ithreads2 -M $imemory -k1,1 > "$o.htsc";
				CMD
			else
				[[ "$feature" == "$level" ]] && params+='-f -O '

				# use -M and --ignoreDup to count also multi-mappings and dulpicates ie. everythin given from bam
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
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

				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					awk 'NR>2{
						print \$1"\t"\$NF
					}' "$o" | helper::sort -t $ithreads2 -M $imemory -k1,1 > "$o.htsc"
				CMD
			fi
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	else
		commander::runcmd -c subread -v -b -i $instances -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
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
			-c <feature>      | use with -w oder -b to output HTSeq-count (htsc) files using feature name instead of bigWig files
			-e <fragmentsize> | when using -w 1 or -m, extend paired-end reads to fragment size and single-end or singleton reads to given fragmentsize
			-b <bed>          | or bed-like file with 3 or 6 columns. latter for strand-specific quantification (see also -x, mutually exclusive to -w)
			-i                | first reference base index of -b option given bed-like file is 1 and not 0 (e.g. gtf, gff)
			-f                | count reads fractional if, according to -o option, they can be assigned to multiple features (+1/n)
			                    NOTE: not supported for version < 2.0.4, when combined with -O 0
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

	local OPTIND arg mandatory skip=false tmpdir="${TMPDIR:-/tmp}" threads overlap=1 windowsize=100 bed index=0 fractional=false unambiguous=false pairs=false outdir fragmentsize maxmemory feature
	declare -n _mapper_bam2bedg _strandness_bam2bedg _nidx_bam2bedg _midx_bam2bedg
	while getopts 'S:s:t:M:r:x:o:O:w:b:e:n:m:c:ifup' arg; do
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
			c)	feature=$OPTARG;;
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
					helper::sort -t $threads -M "$maxmemory" -k1,1 > "${tdirs[-1]}/chr.sizes"
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

	commander::printinfo "convertig bam to bw"

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
						samtools view -H "$f" | sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}' | helper::sort -t $threads -M "$maxmemory" -k1,1 > "${tdirs[-1]}/chr.sizes";
						helper::sort -f "$bed" -t $threads -M "$maxmemory" -k1,1 -k2,2n -k3,3n | awk -v i=$index -v OFS='\t'
					CMD
						'BEGIN{i=i==1?0:1}{id=$1":"$2"-"$3; if($4){id=$4}; s="."; if($6){s=$6}; print id,$1,$2+i,$3,".",".",s}'
					CMD
				else
					commander::makecmd -a cmd1 -s '|' -o "${tdirs[-1]}/regions.saf" -c {COMMANDER[0]}<<- CMD {COMMANDER[2]}<<- 'CMD'
						samtools view -H "$f" | sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}' | helper::sort -t $threads -M "$maxmemory" -k1,1 > "${tdirs[-1]}/chr.sizes";
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

				if [[ $feature ]]; then
					commander::makecmd -a cmd2 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD
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
						| tee -i >(cat > "$odir/$o.${feature}counts") >(awk -v OFS='\t' 'NR>2{print \$1,sprintf("%d",\$7==int(\$7)?\$7:\$7+1)}' > "$odir/$o.${feature}counts.htsc")
					CMD
						| awk '
							NR>2{
								o[NR-3]=$1;
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
						> "$odir/$o.${feature}counts.htsc.tpm"
					CMD
				else
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
						tee -i >(awk -v OFS='\t' 'NR>2{print \$2,\$3-1,\$4,sprintf("%d",\$7==int(\$7)?\$7:\$7+1)}' | bg2bw -i /dev/stdin -o "$odir/$o.$e.counts.bw" -c "${tdirs[-1]}/chr.sizes")
					CMD
						awk '
							NR>2{
								o[NR-3]=$1;
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
						bg2bw -i /dev/stdin -o "$odir/$o.$e.tpm.bw" -c "${tdirs[-1]}/chr.sizes"
					CMD
				fi
			else
				# commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				# 	samtools view -H "$f"
				# CMD
				# 	sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}'
				# CMD
				# 	helper::sort -t $threads -M "$maxmemory" -k1,1 > "${tdirs[-1]}/chr.sizes"
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
			-g <gtf>      | path to with transcript feature (and transcript_id feature tag) for TSS and transcript profiling. mutually exclusive to -b
			-b <bedfiles> | array of paths to bed or bed-like files e.g. narrowPeak. mutually exclusive to -g
			-o <outdir>   | path to
			-p            | pearson correlation plot of coverage files
			-a            | restrict gtf based profiling to canonical transcripts (default: all protein-coding and lncRNA transcripts)
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads outdir tmpdir="${TMPDIR:-/tmp}" gtf bed bwdir pearson=false maxmemory canonicals=false
	declare -n _mapper_profiles _strandness_profiles _bedfiles_profiles
	while getopts 'S:s:t:r:g:b:c:o:M:pa' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			r)	((++mandatory)); _mapper_profiles=$OPTARG;;
			g)	gtf="$OPTARG";;
			b)	_bedfiles_profiles=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			c)	((++mandatory)); bwdir="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			p)	pearson=true;;
			a)	canonicals=true;;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage
	[[ ! $_bedfiles_profiles && ! $gtf ]] && _usage

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
			pileups+=("$(find -L "$bwdir/$m" -maxdepth 1 -name "$(basename ${f%.*}).pileup.tpm.bw" -print -quit | grep .)")
		done

		[[ ${#pileups[@]} -gt 1 ]] && e="$(echo -e "${pileups[0]}\t${pileups[1]}" | sed -E 's/(.+)\t.*\1/\t\1/' | cut -f 2)" || e=".pileup.tpm.bw"

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
					computeMatrix reference-point
						-p $(((threads+1)/2))
						-S $(printf '"%s" ' "${pileups[@]}")
						-R "\$FILE"
						-bs 1
						--beforeRegionStartLength 3000
						--afterRegionStartLength 3000
						--skipZeros
						--missingDataAsZero
						-o "\$FILE.matrix.gz"
				' 2> >(sed -un '/^Skipping/!p' >&2) | cat;
			CMD
				pigz -p 1 -cd "${tdirs[-1]}"/*.matrix.gz | grep -v '^@' | helper::pgzip -t $threads -o "${tdirs[-1]}/data.gz";
				n=\$(gztool -l "${tdirs[-1]}/data.gz" |& sed -nE '/Number of lines/{s/.*:\s+([0-9]+).*/\1/p}');
				pigz -p 1 -cd "${tdirs[-1]}"/*.matrix.gz | head -1 | sed -E 's/"group_boundaries":\[[0-9,]+\]/"group_boundaries":[0,'\$n']/' | pigz -p 1 -kc > "$odir/$b.tss.matrix.gz";
				cat "${tdirs[-1]}/data.gz" >> "$odir/$b.tss.matrix.gz";
			CMD
			# extremely slow and memory hungry: computeMatrixOperations rbind -m "${tdirs[-1]}"/*.matrix.gz -o "$odir/tss.matrix.gz"
			# not necessary when using deeptools plotHeatmap: computeMatrixOperations sort -m "${tdirs[-1]}/tss.matrix.gz" -R "$odir/transcripts.bed" -o "$odir/tss.matrix.gz"
			# use unbuffered sed to avoid "skipping <id> due to being absent in the computeMatrix output" warnings

			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				plotHeatmap
					-m "$odir/$b.tss.matrix.gz"
					-o "$odir/$b.tss.heatmap.pdf"
					--plotFileFormat pdf
					--colorMap RdBu
					--refPointLabel TSS
					--samplesLabel $(basename -a "${pileups[@]/%$e/}" | xargs -I {} printf "'%s' " {})
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
						computeMatrix scale-regions
						-p $(((threads+1)/2))
						-S $(printf '"%s" ' "${pileups[@]}")
						-R "\$FILE"
						-bs 1
						--beforeRegionStartLength 3000
						--regionBodyLength 8000
						--afterRegionStartLength 3000
						--skipZeros
						--missingDataAsZero
						-o "\$FILE.matrix.gz"
				' 2> >(sed -un '/^Skipping/!p' >&2) | cat;
			CMD
				pigz -p 1 -cd "${tdirs[-1]}"/*.matrix.gz | grep -v '^@' | helper::pgzip -t $threads -o "${tdirs[-1]}/data.gz";
				n=\$(gztool -l "${tdirs[-1]}/data.gz" |& sed -nE '/Number of lines/{s/.*:\s+([0-9]+).*/\1/p}');
				pigz -p 1 -cd "${tdirs[-1]}"/*.matrix.gz | head -1 | sed -E 's/"group_boundaries":\[[0-9,]+\]/"group_boundaries":[0,'\$n']/' | pigz -p 1 -kc > "$odir/$b.matrix.gz";
				cat "${tdirs[-1]}/data.gz" >> "$odir/$b.matrix.gz";
			CMD

			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				plotProfile
					-m "$odir/$b.matrix.gz"
					--outFileNameData "$odir/$b.profile.tsv"
					-o "$odir/$b.profile.pdf"
					--plotFileFormat pdf
					--numPlotsPerRow 2
					--samplesLabel $(basename -a "${pileups[@]/%$e/}" | xargs -I {} printf "'%s' " {})
			CMD
		done

		if [[ ${#toprofile[@]} -gt 1 ]]; then
			commander::makecmd -a cmd3 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				{	for f in $(printf '"%s" ' "${toprofile[@]}"); do
						b="\$(basename "\${f%.*}")";
						pigz -p 1 -cd "$odir/\$b.tss.matrix.gz" | head -1 | paste <(echo \$b) -;
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
						pigz -p 1 -kcd "$odir/\$b.tss.matrix.gz" | grep -v ^@;
					done;
				} | helper::pgzip -t $threads -o "$odir/tss.matrix.gz"
			CMD

			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				plotHeatmap
					-m "$odir/tss.matrix.gz"
					-o "$odir/tss.heatmap.pdf"
					--plotFileFormat pdf
					--colorMap RdBu
					--refPointLabel TSS
					--samplesLabel $(basename -a "${pileups[@]/%$e/}" | xargs -I {} printf "'%s' " {})
					--regionsLabel $(basename -a "${toprofile[@]%.*}" | xargs -I {} printf "'%s' " {})
			CMD

			commander::makecmd -a cmd3 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				{	for f in $(printf '"%s" ' "${toprofile[@]}"); do
						b="\$(basename "\${f%.*}")";
						pigz -p 1 -cd "$odir/\$b.matrix.gz" | head -1 | paste <(echo \$b) -;
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
						pigz -p 1 -kcd "$odir/\$b.matrix.gz" | grep -v ^@;
					done;
				} | helper::pgzip -t $threads -o "$odir/matrix.gz"
			CMD

			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
				plotProfile
					-m "$odir/matrix.gz"
					--outFileNameData "$odir/profile.tsv"
					-o "$odir/profile.pdf"
					--plotFileFormat pdf
					--numPlotsPerRow 2
					--samplesLabel $(basename -a "${pileups[@]/%$e/}" | xargs -I {} printf "'%s' " {})
					--regionsLabel $(basename -a "${toprofile[@]%.*}" | xargs -I {} printf "'%s' " {})
			CMD
		fi

		if [[ ${#coverages[@]} -gt 1 ]]; then
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.correlation)")

			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bwcat -i "${coverages[0]}" | cut -f 1-3 > "${tdirs[-1]}/regions.bed"
			CMD
				multiBigwigSummary BED-file
					-p $threads
					--bwfiles $(printf '"%s" ' "${coverages[@]}")
					--BED "${tdirs[-1]}/regions.bed"
					-out "$odir/coverage.npz"
			CMD

			e="$(echo -e "${coverages[0]}\t${coverages[1]}" | sed -E 's/(.+)\t.*\1/\t\1/' | cut -f 2)"
			commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
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
					--labels $(basename -a "${coverages[@]/%$e/}" | xargs -I {} printf "'%s' " {})
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
