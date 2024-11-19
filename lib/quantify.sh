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
						next unless $F[2] eq $f;
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
			quantify::bamcoverage assigns reads to bins, bed files, or piles them up to report TPM/BPM (and input-) normalized bigWig files
			in case of ATAC data with the intention to quantify mono-nucleosome positions or TF ChIP-seq experiments
			- reads can be centered
			- nucleosome free regions or any other input background can be subtracted

			-S <hardskip>     | true/false return
			-s <softskip>     | true/false only print commands
			-t <threads>      | number of cpus. default: 1
			-r <mapper>       | array of bams within array of
			-o <outdir>       | path to

			WINDOW/PILEUP MODE
			-w <windowsize>   | pileup a bam file or count reads within equal sized bins. default: 1 (pileup)
			-e <fragmentsize> | extends reads to fragment size
			-c                | for paired-end data, centers reads using inferred fragment size, 150nt otherwise (adjust via -e option)
			-i <tidx>         | array of IP* bam idices within -r or any other treatment like MNase or mono-nucleosome filtered regions
			-a <nidx>         | array of normal bam idices within -r or any other background to be subtracted like nucleosome free regions
			-m <memory>       | estimate of memory usage per chunk when using -a option. determines parallel jobs (see -M). default: 10000
			-M <maxmemory>    | amount of total memory allowed to use when using -a option. default: all available

			WINDOW MODE
			-w <windowsize>   | count reads within equal sized bins
			-e <fragmentsize> | extends reads to fragment size
			-O <overlap>      | overlap strategy for how to assign reads spanning multiple features. default: 1
			                    0   - largest overlap
			                    0-1 - fraction overlap
			                    1-N - basepair overlap
			-f                | count reads fractional if, according to -O option, they can be assigned to multiple features (+1/n)
			                    NOTE: not supported for featureCounts version < 2.0.4, when combined with -O 0
			-u                | count only unambiguously assigend reads
			-n <feature>      | triggers switch from bigWig to HTSeq-count (htsc) format with given name as prefix for output files

			BED MODE
			-b <bed>          | bed file with 3 or 6 columns, latter for strand-specific quantification
			-x <strandness>   | associative array of strandness information per bam file to do strand-specific quantification. default: 0
			                  | 0 = unstranded, 1 = stranded/fr second strand or 2 = reversely stranded /fr first strand
			-O <overlap>      | overlap strategy for how to assign reads spanning multiple features. default: 1
			                    0   - largest overlap
			                    0-1 - fraction overlap
			                    1-N - basepair overlap
			-f                | count reads fractional if, according to -O option, they can be assigned to multiple features (+1/n)
			                    NOTE: not supported for featureCounts version < 2.0.4, when combined with -O 0
			-u                | count only unambiguously assigend reads
			-n <feature>      | triggers switch from bigWig to HTSeq-count (htsc) format with given feature name as prefix for output files
		EOF
		return 1
	}
	# -p | count read-pairs/fragments instead of reads (experimental!)
	# NOTE1: +1 instead of +2 is counted only if mates are unambiguously assigned to same feature, else +1,+1
	# NOTE2: fractional quantification is not available
	# NOTE3: unless -o option corresponds to basepairs, assignment is derived from the sum of both mate lengths

	local OPTIND arg mandatory skip=false tmpdir="${TMPDIR:-/tmp}" threads=1 overlap=1 windowsize=1 bed index=0 fractional=false unambiguous=false pairs=false outdir fragmentsize maxmemory feature center=false smooth=0 memory
	declare -n _mapper_bamcoverage _strandness_bamcoverage _nidx_bamcoverage _tidx_bamcoverage
	while getopts 'S:s:t:m:M:r:x:o:O:w:b:e:n:a:i:cfup' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	threads=$OPTARG;;
			r)	((++mandatory)); _mapper_bamcoverage=$OPTARG;;
			m)	memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			x)	_strandness_bamcoverage=$OPTARG;;
			o)	outdir="$OPTARG"; mkdir -p "$outdir";;
			O)	overlap=$OPTARG;;
			w)	windowsize=$OPTARG;;
			b)	bed="$OPTARG";;
			f)	fractional=true;;
			u)	unambiguous=true;;
			p)	pairs=true;;
			e)	fragmentsize=$OPTARG;;
			a)	_nidx_bamcoverage=$OPTARG;;
			i)	_tidx_bamcoverage=$OPTARG;;
			n)	feature=$OPTARG;;
			c)	center=true;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 1 ]] && _usage
	[[ $_nidx_bamcoverage ]] && [[ ! $_tidx_bamcoverage ]] && _usage

	commander::printinfo "converting bam to bw"

	declare -n _bams_bamcoverage=${_mapper_bamcoverage[0]}
	if [[ ! $_nidx_bamcoverage && ! $_tidx_bamcoverage ]]; then
		declare -a tidx_bamcoverage=("${!_bams_bamcoverage[@]}") # use all bams as unpaired input unless -a and -i
		_tidx_bamcoverage=tidx_bamcoverage
	fi

	local chrsizes="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.sizes)" saf="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.regions)"
	local x readlength e m f c params o odir
	declare -a cmd1 cmd2 cmd3 cmd4 cmd5 cmd6 tdirs chrs tomerge

	x=$(samtools view -F 4 "${_bams_bamcoverage[0]}" | head -10000 | cat <(samtools view -H "${_bams_bamcoverage[0]}") - | samtools view -c -f 1)
	readlength=$(samtools view "${_bams_bamcoverage[0]}" | head -10000 | awk -F '\t' '{l=length($10);ml=ml<l?l:ml}END{print ml}')
	samtools view -H "${_bams_bamcoverage[0]}" | sed -nE '/^@SQ/{s/.*SN:(\S+).*LN:(\S+)\s*.*/\1\t\2/p}' > "$chrsizes"
	mapfile -t chrs < <(cut -f 1 "$chrsizes")

	if [[ $bed || $windowsize -gt 1 ]]; then
		memory=${memory:-5000}
		e="coverage"
		if [[ $bed ]]; then
			commander::makecmd -a cmd5 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD
				helper::sort -f "$bed" -t $threads -M "$maxmemory" -k1,1 -k2,2n -k3,3n
			CMD
				| awk -F '\t' '{f+=FNR==1?1:0} f==1{id=$1":"$2"-"$3; if($4){id=$4}; s="."; if($6){s=$6}; o[$1]=o[$1]id"\t"$1"\t"$2+1"\t"$3"\t.\t.\t"s"\n"} f==2{printf o[$1]}'
			CMD
				- "$chrsizes" > "$saf"
			CMD
			# commander::makecmd -a cmd5 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[2]}<<- CMD
			# 	awk -v OFS='\t' '{id=$1":"$2"-"$3; if($4){id=$4}; s="."; if($6){s=$6}; print id,$1,$2+1,$3,".",".",s}'
			# CMD
			# 	"$bed" > "$saf"
			# CMD
		else
			commander::makecmd -a cmd5 -s '|' -o "$saf" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
				bedtools makewindows -w $windowsize -g "$chrsizes"
			CMD
				awk -v OFS='\t' '{id=$1":"$2"-"$3; print id,$1,$2+1,$3,".",".","."}'
			CMD
		fi
	else
		memory=${memory:-10000}
		e="pileup"
	fi

	local instances ithreads instances2 ithreads2 i=$((${#_mapper_bamcoverage[@]}*${#_bams_bamcoverage[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -t 10 -T $threads)
	read -r instances2 ithreads2 < <(configure::instances_by_memory -m "$memory" -M "$maxmemory" -T $threads)
	read -r instances3 ithreads3 < <(configure::instances_by_threads -i $i -t 64 -T $threads)
	[[ $instances3 -lt $i && $ithreads3 -gt 64 ]] && ((++instances3)) && ithreads3=$((threads/instances3))
	[[ $ithreads3 -gt 64 ]] && ithreads3=64

	for m in "${_mapper_bamcoverage[@]}"; do
		declare -n _bams_bamcoverage=$m
		# odir="$outdir/$m"
		# mkdir -p "$odir"
		for i in "${!_tidx_bamcoverage[@]}"; do
			f="${_bams_bamcoverage[${_tidx_bamcoverage[$i]}]}"
			[[ $outdir ]] && odir="$outdir/$m" && mkdir -p "$odir" || odir="$(dirname "$f")"
			o="$(basename "${f%.*}")"
			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bamcoverage)")

			if $center; then
				# MNase-seq or mono-nucleosome filtered ATAC-seq or TF ChIP-seq
				# for the latter, extend size should actually be fragmentsize/2 but number bp wrapping a histone (~147) is also fine
				# info: not really recommended for SE data

				if [[ $_nidx_bamcoverage ]]; then
					n="${_bams_bamcoverage[${_nidx_bamcoverage[$i]}]}"
					[[ $x -gt 0 ]] && params="--paired 1" || params="--frsz ${fragmentsize:-150}"

					# danpos finds center using PE tlen infered frsz or given frsz (for SE) and extends center in both directions by extend/2 (default: --extend 80)
					# without given -b for subtraction, outfile will be "${tdirs[-1]}/pooled/$(basename "${f%.*}").smooth.wig"
					# without smooth, outfile will be "${tdirs[-1]}/pooled/$(basename "${f%.*}").bgsub.wig"
					# trailing zeroes are not reported! -> apply wiggletools
					tomerge=()
					cmd5=()
					for c in "${chrs[@]}"; do
						commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
							{	samtools view -@ $ithreads -H "$f" | grep -E "^@SQ\s+SN:$c\s+";
								samtools view -@ $ithreads "$f" $c;
							} | samtools view -@ $ithreads --write-index -o "${tdirs[-1]}/$c.bam##idx##${tdirs[-1]}/$c.bai"
						CMD
							{	samtools view -@ $ithreads -H "$n" | grep -E "^@SQ\s+SN:$c\s+";
								samtools view -@ $ithreads "$n" $c;
							} | samtools view -@ $ithreads --write-index -o "${tdirs[-1]}/$c.bg.bam##idx##${tdirs[-1]}/$c.bg.bai"
						CMD

						commander::makecmd -a cmd2 -s ' ' -c <<- CMD
							cd "${tdirs[-1]}";
							read -r size reads < <(samtools idxstats "$c.bam" | head -1 | cut -f 2-3);
							if [[ \$reads -eq 0 ]]; then
								echo -e "$c\t0\t\$size\t0" > "$c.bedg";
							else
								mkdir -p pooled;
								exec 11>&1;
								ln -sn /dev/fd/11 "pooled/$c.bgsub.wig";
								danpos.py dtriple
									"$c.bam"
									-b "$c.bg.bam"
									--extend $((${fragmentsize:-150}/2))
									-a $windowsize
									-o "${tdirs[-1]}"
									-n F
									-jd 20
									-z $smooth
									-r 0
									-k 0
									-j 0
									$params
								11> >(wiggletools cat - | sed -E "\\\${s/\S+\t([0\.]+)\\\$/\$size\t\1/;t;p;s/\S+\t(\S+)\t\S+\\\$/\1\t\$size\t0/}" > "$c.bedg") | cat;
								exec 11>&-;
								rm "$c.bam" "$c.bg.bam";
							fi
						CMD

						tomerge+=("${tdirs[-1]}/$c.bedg")
					done
					# bigWigMerge much faster than bwjoin -f "${tdirs[-1]}/chr.sizes" -p "${tdirs[-1]}" -o "$odir/$o.$e.tpm.bw"
					# and slightly faster than wiggletools sum
					# in any case this requires 11> >(wigToBigWig /dev/stdin "$c.size" "$c.bw") | cat;
					# -> drawback: reordering of chromosomes by bigWigMerge -threshold=-1 $(printf '"%s" ' "${tomerge[@]}") /dev/stdout
					# -> use bedgraph files
					commander::makecmd -a cmd3 -c <<- CMD
						cat $(printf '"%s" ' "${tomerge[@]}") | bg2bw -i /dev/stdin -c "$chrsizes" -o "$odir/$o.$e.tpm.bw"
					CMD

					# alternative, but due to low value range, coverage gets saw like shapes, which looks very unnatural, even with smoothing
					# plus log2 ratio produces negative values
					# likewise to bamCoverage, bamCompare centers reads by tlen inferred frzs and does not extend them. in case of singletons or SE fallback to given frzs
					# commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					# 	TMPDIR="${tdirs[-1]}" bamCompare
					# 		-b1 "$f"
					# 		-b2 "$n"
					# 		-bs $windowsize
					# 		-p $threads
					# 		--normalizeUsing BPM
					# 		--exactScaling
					# 		--extendReads ${fragmentsize:-150}
					# 		--centerReads
					# 		--maxFragmentLength 1000
					# 		--scaleFactorsMethod None
					# 		--smoothLength 10
					# 		--outFileFormat bigwig
					# 		-o "${tdirs[-1]}/$o.$e.tpm.bw"
					# CMD
					# commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					# 	bwcat -i "${tdirs[-1]}/$o.pileup.tpm.bw"
					# CMD
					# 	awk -v OFS='\t' '$NF<0{$NF=0}{print}'
					# CMD
					# 	bg2bw -i /dev/stdin -c "$chrsizes" -o "$odir/$o.$e.tpm.bw"
					# CMD
				else
					# --exactScaling ??
					# --smoothLength $smooth
					commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
						TMPDIR="${tdirs[-1]}" bamCoverage
							-b "$f"
							-bs $windowsize
							-p $threads
							--normalizeUsing BPM
							--extendReads ${fragmentsize:-150}
							--centerReads
							--maxFragmentLength 1000
							--outFileFormat bigwig
							-o "$odir/$o.$e.tpm.bw"
					CMD
				fi
			else
				# no centering
				if [[ $_nidx_bamcoverage ]]; then
					n="${_bams_bamcoverage[${_nidx_bamcoverage[$i]}]}"

					tomerge=()
					cmd5=()
					for c in "${chrs[@]}"; do
						commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
							{	samtools view -@ $ithreads -H "$f" | grep -E "^@SQ\s+SN:$c\s+";
								samtools view -@ $ithreads "$f" $c;
							} | samtools view -@ $ithreads --write-index -o "${tdirs[-1]}/$c.bam##idx##${tdirs[-1]}/$c.bai"
						CMD
							{	samtools view -@ $ithreads -H "$n" | grep -E "^@SQ\s+SN:$c\s+";
								samtools view -@ $ithreads "$n" $c;
							} | samtools view -@ $ithreads --write-index -o "${tdirs[-1]}/$c.bg.bam##idx##${tdirs[-1]}/$c.bg.bai"
						CMD

						# in case of given fragment length, i.e. to pileup broad signals, extend to frsz
						commander::makecmd -a cmd2 -s ' ' -c <<- CMD
							cd "${tdirs[-1]}";
							read -r size reads < <(samtools idxstats "$c.bam" | head -1 | cut -f 2-3);
							if [[ \$reads -eq 0 ]]; then
								echo -e "$c\t0\t\$size\t0" > "$c.bedg";
							else
								mkdir -p pooled;
								exec 11>&1;
								ln -sn /dev/fd/11 "pooled/$c.bgsub.wig";
								danpos.py dtriple
									"$c.bam"
									-b "$c.bg.bam"
									--frsz $readlength
									--extend ${fragmentsize:-$readlength}
									-a $windowsize
									-o "${tdirs[-1]}"
									-n F
									-jd 20
									-z $smooth
									-r 0
									-k 0
									-j 0
								11> >(wiggletools cat - | sed -E "\\\${s/\S+\t([0\.]+)\\\$/\$size\t\1/;t;p;s/\S+\t(\S+)\t\S+\\\$/\1\t\$size\t0/}" > "$c.bedg") | cat;
								exec 11>&-;
								rm "$c.bam" "$c.bg.bam";
							fi
						CMD

						tomerge+=("${tdirs[-1]}/$c.bedg")
					done

					commander::makecmd -a cmd3 -c <<- CMD
						cat $(printf '"%s" ' "${tomerge[@]}") | bg2bw -i /dev/stdin -c "$chrsizes" -o "$odir/$o.$e.tpm.bw"
					CMD
				else
					if [[ ! $bed && $windowsize -eq 1 ]]; then
						# in case of given fragment length, i.e. to pileup broad signals, extend to frsz
						[[ $fragmentsize ]] && params=" --extendReads $fragmentsize" || params=""
						# --exactScaling ??
						# --smoothLength $smooth
						commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
							TMPDIR="${tdirs[-1]}" bamCoverage
								$params
								-b "$f"
								-bs $windowsize
								-p $threads
								--normalizeUsing BPM
								--maxFragmentLength 1000
								--outFileFormat bigwig
								-o "$odir/$o.$e.tpm.bw"
						CMD
					else
						local version
						if [[ ! $version ]]; then
							declare -a cmdchk=("featureCounts -v |& grep -oE 'v[.0-9]+'")
							version=$(commander::runcmd -c subread -a cmdchk)
							# $fractional && [[ $overlap -eq 0 ]] && [[ "$(echo -e "v2.0.3\n$version" | sort -Vr | head -1)" == "v2.0.3" ]] && _usage
							[[ "$(echo -e "v2.0.1\n$version" | sort -Vr | head -1)" == "v2.0.1" ]] && version="old" || version="new"
						fi

						# in case of given fragment length, i.e. to pileup broad signals, extend to frsz
						[[ $fragmentsize ]] && params="--readExtension3 $((fragmentsize-readlength))" || params=""
						$unambiguous || params+=" -O"
						$fractional && params+=" --fraction"
						params+=" $(echo $overlap | awk '{if($1==0){print "--largestOverlap"}else{print $1<1? "--fracOverlap "$1 : "--minOverlap "$1}}')"
						if [[ $version == "old" ]]; then
							$pairs && params+=" -p"
						else
							# ATTENTION: from v2.0.3 on, no SE quantification possible if data is PE due to internal PE checker that leads to termination
							x=$(samtools view -F 4 "$f" | head -10000 | cat <(samtools view -H "$f") - | samtools view -c -f 1)
							if [[ $x -gt 0 ]]; then
								params+=" -p"
								$pairs && params+=" --countReadPairs"
							fi
						fi

						if [[ $feature ]]; then
							commander::makecmd -a cmd6 -s ' ' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD
								featureCounts
									$params
									-Q 0
									--maxMOp 999
									-s $(declare -p _strandness_bamcoverage | grep -q '=' && echo ${_strandness_bamcoverage["$f"]} || echo 0)
									-T $ithreads3
									-f
									-M
									--ignoreDup
									--tmpDir "${tdirs[-1]}"
									-F SAF
									-a "$saf"
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
							commander::makecmd -a cmd6 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD' {COMMANDER[3]}<<- CMD
								featureCounts
									$params
									-Q 0
									--maxMOp 999
									-s $(declare -p _strandness_bamcoverage | grep -q '=' && echo ${_strandness_bamcoverage["$f"]} || echo 0)
									-T $ithreads3
									-f
									-M
									--ignoreDup
									--tmpDir "${tdirs[-1]}"
									-F SAF
									-a "$saf"
									-o /dev/stdout
									"$f"
							CMD
								tee -i >(awk -v OFS='\t' 'NR>2{print \$2,\$3-1,\$4,sprintf("%d",\$7==int(\$7)?\$7:\$7+1)}' | bg2bw -i /dev/stdin -o "$odir/$o.$e.counts.bw" -c "$chrsizes")
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
								bg2bw -i /dev/stdin -o "$odir/$o.$e.tpm.bw" -c "$chrsizes"
							CMD
						fi
					fi
				fi
			fi
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
		commander::printcmd -a cmd6
	else
		commander::runcmd -v -b -i $instances -a cmd1
		commander::runcmd -c danpos -v -b -i $instances2 -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
		commander::runcmd -c deeptools -v -b -i 1 -a cmd4
		commander::runcmd -v -b -i 1 -a cmd5
		commander::runcmd -c subread -v -b -i $instances3 -a cmd6
	fi

	return 0
}