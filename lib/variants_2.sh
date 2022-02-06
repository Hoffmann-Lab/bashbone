#! /usr/bin/env bash
# (c) Konstantin Riege

variants::haplotypecaller() {
	local tmpfile
	declare -a tdirs
	_cleanup::variants::haplotypecaller(){
		rm -f $tmpfile
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
			-e <isdna>     | true for dna, false for rna (default: true)
			-m <memory>    | amount of
			-M <maxmemory> | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory maxmemory genome dbsnp tmpdir outdir isdna=true
	declare -n _mapper_haplotypecaller _bamslices_haplotypecaller
	while getopts 'S:s:t:g:d:e:m:M:r:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			m) ((++mandatory)); memory=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			e) isdna="$OPTARG";;
			r) ((++mandatory)); _mapper_haplotypecaller=$OPTARG;;
			c) ((++mandatory)); _bamslices_haplotypecaller=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage

	if [[ ! $dbsnp ]]; then
		tmpfile="$(mktemp -p "$tmpdir" cleanup.XXXXXXXXXX.vcf.gz)"
		dbsnp="$tmpfile"
		echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" | bgzip -f -@ $threads > "$tmpfile"
		tabix -f -p vcf "$tmpfile"
		dbsnp="$tmpfile"
	fi

	commander::printinfo "calling variants haplotypecaller"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory -M "$maxmemory")

	local m i o t e slice instances ithreads odir tdir
	for m in "${_mapper_haplotypecaller[@]}"; do
		declare -n _bams_haplotypecaller=$m
		((instances+=${#_bams_haplotypecaller[@]}))
	done
	read -r i memory < <(configure::memory_by_instances -i $((instances*4)) -T $threads -M "$maxmemory") # for final bcftools sort of vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6
	for m in "${_mapper_haplotypecaller[@]}"; do
		declare -n _bams_haplotypecaller=$m
		tdir="$tmpdir/$m"
		odir="$outdir/$m/gatk"
		mkdir -p "$odir" "$tdir"

		for i in "${!_bams_haplotypecaller[@]}"; do
			tomerge=()

			while read -r slice; do
				# all alt Phredscaled Likelihoods ordering:
				# reference homozygous (0/0)
				# ref and alt 1 heterzygous (0/1)
				# alt 1 homozygous (1/1)
				# ref and alt 2 heterzygous (0/2)
				# alt 1 and alt 2 heterozygous (1/2)
				# alt 2 homozygous (2/2)
				# gatk bug as of v4.1.2.0 --max-reads-per-alignment-start 0 not a valid option
				# HaplotypeCallerSpark with --spark-master local[$mthreads] is BETA and differs in results!!
				# Spark does not yet support -D "$dbsnp"
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
								-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
							'
						HaplotypeCaller
						-I "$slice"
						-O "$slice.unfiltered.vcf"
						-R "$genome"
						-D "$dbsnp"
						-A Coverage
						-A DepthPerAlleleBySample
						-A StrandBiasBySample
						--dont-use-soft-clipped-bases true
						--smith-waterman FASTEST_AVAILABLE
						--max-reads-per-alignment-start 0
						--min-base-quality-score 20
						--minimum-mapping-quality 0
						--standard-min-confidence-threshold-for-calling 20
						--native-pair-hmm-threads $mthreads
						--max-alternate-alleles 3
						--all-site-pls true
						-verbosity INFO
						--tmp-dir "${tdirs[-1]}"
				CMD

				if $isdna; then
					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
						gatk
							--java-options '
									-Xmx${jmem}m
									-XX:ParallelGCThreads=$jgct
									-XX:ConcGCThreads=$jcgct
									-Djava.io.tmpdir="$tmpdir"
									-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
								'
							VariantFiltration
							-V "$slice.unfiltered.vcf"
							-O "$slice.vcf"
							-filter "QD < 2.0"
							--filter-name "QD2"
							-filter "FS > 30.0"
							--filter-name "FS30"
							-filter "SOR > 3.0"
							--filter-name "SOR3"
							-filter "MQ < 40.0"
							--filter-name "MQ40"
							-filter "MQRankSum < -3.0"
							--filter-name "MQRankSum-3"
							-filter "ReadPosRankSum < -3.0"
							--filter-name "ReadPosRankSum-3"
							-verbosity INFO
							--tmp-dir "${tdirs[-1]}"
					CMD
				else
					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
						gatk
							--java-options '
									-Xmx${jmem}m
									-XX:ParallelGCThreads=$jgct
									-XX:ConcGCThreads=$jcgct
									-Djava.io.tmpdir="$tmpdir"
									-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
								'
							VariantFiltration
							-V "$slice.unfiltered.vcf"
							-O "$slice.vcf"
							--cluster-window-size 35
							--cluster-size 3
							-filter "QD < 2.0"
							--filter-name "QD2"
							-filter "FS > 30.0"
							--filter-name "FS30"
							-verbosity INFO
							--tmp-dir "${tdirs[-1]}"
					CMD
				fi

				commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				if [[ $dbsnp ]]; then
					commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl
					CMD
						vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$slice.fixed.nomulti.normed.vcf"
					CMD
				else
					commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
					CMD
				fi

				tomerge+=("$slice")
			done < "${_bamslices_haplotypecaller[${_bams_haplotypecaller[$i]}]}"

			o="$(basename "${_bams_haplotypecaller[$i]}")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
				#DO NOT PIPE - DATALOSS!
				commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$t.$e" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T "${tdirs[-1]}" -m ${memory}M -o "$o.$e" "$t.$e"
				CMD

				commander::makecmd -a cmd7 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
		commander::printcmd -a cmd6
		commander::printcmd -a cmd7
	else
		commander::runcmd -c gatk -v -b -t $minstances -a cmd1
		commander::runcmd -c gatk -v -b -t $minstances -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
		commander::runcmd -v -b -t $threads -a cmd4
		commander::runcmd -v -b -t $threads -a cmd5
		commander::runcmd -v -b -t $threads -a cmd6
		commander::runcmd -v -b -t $instances -a cmd7
	fi

	return 0
}

variants::mutect() {
	declare -a tdirs
	_cleanup::variants::mutect(){
		rm -rf "${tdirs[@]}"
	}

	# You do not need to make you own panel of normals (unless you have a huge number of samples,
	# it may even be counterproductive than our generic public panel).
	# Instead you may use gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf.
	# For the -germline-resource you should use gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.
	# -> this seems to be replacing old dbSNP input, which does not have AF info field
	# For the -V and -L arguments to GetPileupSummaries you may use gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf.
	# -> can this be used? ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/GATK/00-common_all.vcf.gz
	# -> cam this be used? ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/homo_sapiens_somatic.vcf.gz
	# BUT b37 => HG19 ...
	# https://software.broadinstitute.org/gatk/download/bundle
	# https://console.cloud.google.com/storage/browser/gatk-best-practices/somatic-hg38/
	# BUT this is b37 (similar to hg19) liftovered data (picard LiftoverVCF ? ENSEMBL/NCBI remap?)
	# from helsinki workshop to create gatk tutotial webpage : https://software.broadinstitute.org/gatk/documentation/article?id=11136
	# data downloaded to /misc/paras/data/genomes/GRCh38.p12/GATK-resources

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-m <memory>    | amount of
			-M <maxmemory> | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (requieres -2)
			-2 <tumoridx>  | array of
			-c <sliceinfo> | array of
			-d <pondb>     | path to
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory maxmemory genome tmpdir outdir pondb dbsnp
	declare -n _mapper_mutect _bamslices_mutect _nidx_mutect _tidx_mutect
	while getopts 'S:s:t:g:d:n:m:M:r:1:2:c:p:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			M)	maxmemory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			d)	dbsnp="$OPTARG";;
			n)	pondb="$OPTARG";;
			r)	((++mandatory)); _mapper_mutect=$OPTARG;;
			c)	((++mandatory)); _bamslices_mutect=$OPTARG;;
			1)	((++mandatory)); _nidx_mutect=$OPTARG;;
			2)	((++mandatory)); _tidx_mutect=$OPTARG;;
			p)	((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 9 ]] && _usage

	commander::printinfo "calling variants mutect"

	local minstances mthreads jmem jgct jcgct minstances2 mthreads2 jmem2 jgct2 jcgct2
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory -M "$maxmemory")
	local params params2 m i o t nslice slice odir tdir ithreads instances=$((${#_mapper_mutect[@]}*${#_tidx_mutect[@]}))
	read -r i memory < <(configure::memory_by_instances -i $((instances*4)) -T $threads -M "$maxmemory") # for final bcftools sort of vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -T $threads)
	read -r minstances2 mthreads2 jmem2 jgct2 jcgct2 < <(configure::jvm -i $instances -T $threads -M "$maxmemory")

	[[ $pondb ]] && params="-pon '$pondb'"
	[[ -s "$genome.af_only_gnomad.vcf.gz" ]] && params+=" --germline-resource '$genome.af_only_gnomad.vcf.gz'"

	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6 cmd7 cmd8 cmd9 cmd10
	for m in "${_mapper_mutect[@]}"; do
		declare -n _bams_mutect=$m
		odir="$outdir/$m/gatk"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		for i in "${!_tidx_mutect[@]}"; do
			tomerge=()

			while read -r nslice slice; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
				# normal name defined for RGSM sam header entry by alignment::addreadgroup
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
								-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
							'
						Mutect2
						$params
						-I "$nslice"
						-I "$slice"
						-normal NORMAL
						-tumor TUMOR
						-O "$slice.unfiltered.vcf"
						-R "$genome"
						-A Coverage
						-A DepthPerAlleleBySample
						-A StrandBiasBySample
						--f1r2-tar-gz "$slice.f1r2.tar.gz"
						--smith-waterman FASTEST_AVAILABLE
						--f1r2-median-mq 0
						--f1r2-min-bq 20
						--callable-depth 10
						--min-base-quality-score 20
						--minimum-mapping-quality 0
						--minimum-allele-fraction 0.1
						--native-pair-hmm-threads $mthreads
						--verbosity INFO
						--tmp-dir "${tdirs[-1]}"
						--max-mnp-distance 0
						--ignore-itr-artifacts true
						--dont-use-soft-clipped-bases true
				CMD
				# vs gatk3: removed -cosmic and -dbsnp. The best practice now is to pass in gnomAD as the germline-resource

				# --max-mnp-distance 0
				# default:1 - set to 0 to list adjacent SNPs as single entries rather than an MNP
				# useful for my SNP based genotype tree inference and similar to vardict -X 0

				# --ignore-itr-artifacts to avoid "Fasta index file could not be opened" error after running for ~2h due to hardlimit of parallel filehandlers (ulimit -Hn)
				# This flag disables a read transformer that repairs a very rare artifact in which transposons and other repeated elements self-hybridize and then T4 polymerase incorrectly fills in overhanging ends using the repeat DNA as the template.
				# An other solution would be to raise ulimit to >40k (!) which requires root. fd will not be nuked by garbace collection even if gc threads increased.

				# to apply LearnReadOrientationMode use --f1r2-tar-gz option

				# --dont-use-soft-clipped-bases recommended for RNA based calling due to splitNcigar which elongates both splits towards full length read alignments and mask them\n\t\t\t\t# since we are doing clipmateoverlap, which introduces soft-clips, always enable this option

				if [[ -s "$genome.small_common.vcf.gz" ]]; then
					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
						gatk
							--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
								-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
							'
							GetPileupSummaries
							-I "$slice"
							-V "$genome.small_common.vcf.gz"
							-L "$genome.small_common.vcf.gz"
							-O "$slice.pileups.table"
							--verbosity INFO
							--tmp-dir "${tdirs[-1]}"
					CMD


					##########
					#TODO awk 'NR>2' *slice?.bam.pileups.table | sort -k1,1V -k2,2n  BUT -k1,1 must be in order of $genome.list
					# then predict contamination as whole and filter on full set
					##########
					commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
						gatk
							--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
								-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
							'
							CalculateContamination
							-I "$slice.pileups.table"
							-O "$slice.contamination.table"
							--tumor-segmentation "$slice.tumorsegments.table"
							--verbosity INFO
							--tmp-dir "${tdirs[-1]}"
					CMD
					params2=" --contamination-table '$slice.contamination.table' --tumor-segmentation '$slice.tumorsegments.table'"
				fi

				commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
					gatk
						--java-options '
							-Xmx${jmem}m
							-XX:ParallelGCThreads=$jgct
							-XX:ConcGCThreads=$jcgct
							-Djava.io.tmpdir="$tmpdir"
							-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
						'
						LearnReadOrientationModel
						-I "$slice.f1r2.tar.gz"
						-O "$slice.f1r2_model.tar.gz"
						--tmp-dir "${tdirs[-1]}"
				CMD

				commander::makecmd -a cmd5 -s ';' -c {COMMANDER[0]}<<- CMD
					gatk
						--java-options '
							-Xmx${jmem}m
							-XX:ParallelGCThreads=$jgct
							-XX:ConcGCThreads=$jcgct
							-Djava.io.tmpdir="$tmpdir"
							-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
						'
						FilterMutectCalls
						$params2
						--ob-priors "$slice.f1r2_model.tar.gz"
						-R "$genome"
						-V "$slice.unfiltered.vcf"
						-stats "$slice.unfiltered.vcf.stats"
						--min-median-base-quality 20
						--min-median-mapping-quality 0
						-O "$slice.vcf"
						--verbosity INFO
						--tmp-dir "${tdirs[-1]}"
				CMD
				#In order to tweak results in favor of more sensitivity users may set -f-score-beta to a value greater than its default of 1 (beta is the relative weight of sensitivity versus precision in the harmonic mean). Setting it lower biases results toward greater precision.
				# optional: --tumor-segmentation segments.table

				commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd7 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				if [[ $dbsnp ]]; then
					commander::makecmd -a cmd8 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl
					CMD
						vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$slice.fixed.nomulti.normed.vcf"
					CMD
				else
					commander::makecmd -a cmd8 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
					CMD
				fi

				tomerge+=("$slice")
			done < <(paste "${_bamslices_mutect[${_bams_mutect[${_nidx_mutect[$i]}]}]}" "${_bamslices_mutect[${_bams_mutect[${_tidx_mutect[$i]}]}]}")

			o="$(basename "${_bams_mutect[${_tidx_mutect[$i]}]}")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
			commander::makecmd -a cmd9 -s ';' -c {COMMANDER[0]}<<- CMD
				gatk
					--java-options '
						-Xmx${jmem2}m
						-XX:ParallelGCThreads=$jgct2
						-XX:ConcGCThreads=$jcgct2
						-Djava.io.tmpdir="$tmpdir"
						-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
					'
					MergeMutectStats
					$(printf ' -stats "%s" ' "${tomerge[@]/%/.unfiltered.vcf.stats}")
					-O "$o.vcf.stats"
					--verbosity INFO
					--tmp-dir "${tdirs[-1]}"
			CMD

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				tdirs+=("$(mktemp -d -p "$tdir" cleanup.XXXXXXXXXX.bcftools)")
				#DO NOT PIPE - DATALOSS!
				commander::makecmd -a cmd10 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$t.$e" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T "${tdirs[-1]}" -m ${memory}M -o "$o.$e" "$t.$e"
				CMD

				commander::makecmd -a cmd11 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done

		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
		commander::printcmd -a cmd6
		commander::printcmd -a cmd7
		commander::printcmd -a cmd8
		commander::printcmd -a cmd9
		commander::printcmd -a cmd10
		commander::printcmd -a cmd11
	else
		commander::runcmd -c gatk -v -b -t $minstances -a cmd1
		commander::runcmd -c gatk -v -b -t $minstances -a cmd2
		commander::runcmd -c gatk -v -b -t $minstances -a cmd3
		commander::runcmd -c gatk -v -b -t $minstances -a cmd4
		commander::runcmd -c gatk -v -b -t $minstances -a cmd5
		commander::runcmd -v -b -t $threads -a cmd6
		commander::runcmd -v -b -t $threads -a cmd7
		commander::runcmd -v -b -t $threads -a cmd8
		commander::runcmd -c gatk -v -b -t $minstances2 -a cmd9
		commander::runcmd -v -b -t $threads -a cmd10
		commander::runcmd -v -b -t $instances -a cmd11
	fi

	return 0
}

variants::bcftools() {
	declare -a tdirs
	_cleanup::variants::bcftools(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-M <maxmemory> | amount of
			-g <genome>    | path to
			-d <dbsnp>     | path to
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (for somatic calling. requieres -2)
			-2 <tumoridx>  | array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory genome dbsnp tmpdir outdir
	declare -n _mapper_bcftools _nidx_bcftools _tidx_bcftools
	while getopts 'S:s:t:M:g:d:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((++mandatory)); _mapper_bcftools=$OPTARG;;
			1) _nidx_bcftools=$OPTARG;;
			2) _tidx_bcftools=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "calling variants bcftools"

	local m=${_mapper_bcftools[0]}
	declare -n _bams_bcftools=$m
	[[ ! $_tidx_bcftools ]] && unset _tidx_bcftools && declare -a _tidx_bcftools=(${!_bams_bcftools[@]}) # use all indices for germline calling

	local ithreads i memory instances=$((${#_mapper_bcftools[@]}*${#_tidx_bcftools[@]}))
	read -r i memory < <(configure::memory_by_instances -i $((instances*4)) -T $threads -M "$maxmemory") # for final bcftools sort of vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads) # for final bgzip

	tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed)")
	declare -a cmd1slice cmd2slice
	commander::makecmd -a cmd1slice -s ' ' -o "${tdirs[0]}/tmp" -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		perl -lane '
			$i=0;
			for (map {$_*100000} 1..($F[1]-1)/100000){
				$i=$_;
				print join"\t",($F[0],$_-100000,$_-1);
			}
			print join"\t",($F[0],$i,$F[1]-1) if $i<$F[1]-1;
		'
	CMD
		"$genome.fai"
	CMD

	commander::makecmd -a cmd2slice -s ' ' -c {COMMANDER[0]}<<- CMD
		split
		--numeric-suffixes=1
		-a ${#threads}
		-n l/$threads
		--additional-suffix=.bed
		--filter='bedtools merge -d 1 -i - > \$(dirname \$FILE)/\$(basename \$FILE | sed "s/^x0*/slice./")'
		"${tdirs[0]}/tmp" "${tdirs[0]}/x"
	CMD

	if $skip; then
		commander::printcmd -a cmd1slice
		commander::printcmd -a cmd2slice
	else
		commander::runcmd -v -b -t $threads -a cmd1slice
		commander::runcmd -v -b -t $threads -a cmd2slice
	fi


	local m i j f nf o t e slice odir tdir
	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6
	for m in "${_mapper_bcftools[@]}"; do
		declare -n _bams_bcftools=$m
		tdir="$tmpdir/$m"
		odir="$outdir/$m/bcftools"
		mkdir -p "$odir" "$tdir"

		for i in "${!_tidx_bcftools[@]}"; do
			f="${_bams_bcftools[${_tidx_bcftools[$i]}]}"
			o="$(basename "$f")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"
			tomerge=()

			for j in $(seq 1 $threads); do
				slice="$t.slice.$j"

				# may add --ploidy-file  CHROM, FROM, TO, SEX, PLOIDY  (unlisted regions have ploidy 2)
				# should be used with --samples-file <(echo "SAMPLE M"), when used without, chrY vcf files are empty.
				# workaround: set sex to A and X,Y,M to 1. warning: GT will be 1 only (not 1/1 or 0/1)
				# chrX 1 999999999 M 1
				# X 1 999999999 M 1
				# chrY 1 999999999 M 1
				# Y 1 999999999 M 1
				# chrY 1 999999999 F 0
				# Y 1 999999999 F 0
				# chrMT 1 999999999 M 1
				# chrM 1 999999999 M 1
				# MT 1 999999999 M 1
				# M 1 999999999 M 1
				# chrM 1 999999999 F 1
				# MT 1 999999999 F 1
				# M 1 999999999 F 1

				if [[ $_nidx_bcftools ]]; then # somatic
					nf="${_bams_bcftools[${_nidx_bcftools[$i]}]}"

					# print orphans A and overlaps x, since this can be adressed by the user via post-aln functions
					commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						bcftools mpileup
							-A
							-x
							-B
							-Q 20
							-m 10
							-d 100000
							-L 100000
							-a 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR'
							-f "$genome"
							-R "${tdirs[0]}/slice.$j.bed"
							-O u
							"$nf" "$f"
					CMD
						bcftools call
							-m
							-f GQ,GP
							-v
					CMD
						vcffilter -v -f "DP < 10"
					CMD
						vcfsamplediff VCFSAMPLEDIFF NORMAL TUMOR - > "$slice.vcf"
					CMD
					# use to filter for germline risk
					# do not use --strict, which requires that no observation in the germline support the somatic alternate i.e. no 0/1 0/1
					# grep -E '(^#|VCFSAMPLEDIFF=somatic)' > "$slice.vcf"
					#CMD # do not hard filter here. vcfix.pl adds vcfix_somatic FILTER tag if somatic or loh and GQ MAF thresholds passed
				else
					commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
						bcftools mpileup
							-A
							-x
							-B
							-Q 20
							-m 10
							-d 100000
							-L 100000
							-a 'DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR'
							-f "$genome"
							-R "${tdirs[0]}/slice.$j.bed"
							-O u
							"$f"
					CMD
						bcftools call
							-m
							-f GQ,GP
							-v
					CMD
						vcffilter -v -f "DP < 10" > "$slice.vcf"
					CMD
				fi
				# attention: conda bcftools may spawn too many file descriptors accoring to ulimit

				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				if [[ $dbsnp ]]; then
					commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl
					CMD
						vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$slice.fixed.nomulti.normed.vcf"
					CMD
				else
					commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
					CMD
				fi

				tomerge+=("$slice")
			done

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
				#DO NOT PIPE - DATALOSS!
				commander::makecmd -a cmd5 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$t.$e" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T "${tdirs[-1]}" -m ${memory}M -o "$o.$e" "$t.$e"
				CMD

				commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done
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
		commander::runcmd -v -b -t $threads -a cmd1
		commander::runcmd -v -b -t $threads -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
		commander::runcmd -v -b -t $threads -a cmd4
		commander::runcmd -v -b -t $threads -a cmd5
		commander::runcmd -v -b -t $instances -a cmd6
	fi

	return 0
}

variants::freebayes() {
	declare -a tdirs
	_cleanup::variants::freebayes(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-M <maxmemory> | amount of
			-g <genome>    | path to
			-d <dbsnp>     | path to
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (for somatic calling. requieres -2)
			-2 <tumoridx>  | array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory genome dbsnp tmpdir outdir
	declare -n _mapper_freebayes _nidx_freebayes _tidx_freebayes
	while getopts 'S:s:t:M:g:d:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((++mandatory)); _mapper_freebayes=$OPTARG;;
			1) _nidx_freebayes=$OPTARG;;
			2) _tidx_freebayes=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "calling variants freebayes"

	local m=${_mapper_freebayes[0]}
	declare -n _bams_freebayes=$m
	[[ ! $_tidx_freebayes ]] && unset _tidx_freebayes && declare -a _tidx_freebayes=(${!_bams_freebayes[@]}) # use all indices for germline calling

	local ithreads i memory instances=$((${#_mapper_freebayes[@]}*${#_tidx_freebayes[@]}))
	read -r i memory < <(configure::memory_by_instances -i $((instances*4)) -T $threads -M "$maxmemory") # for final bcftools sort of vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads) # for final bgzip

	tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed)")
	declare -a cmd1slice cmd2slice
	commander::makecmd -a cmd1slice -s ' ' -o "${tdirs[0]}/tmp" -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		perl -lane '
			$i=0;
			for (map {$_*100000} 1..($F[1]-1)/100000){
				$i=$_;
				print join"\t",($F[0],$_-100000,$_-1);
			}
			print join"\t",($F[0],$i,$F[1]-1) if $i<$F[1]-1;
		'
	CMD
		"$genome.fai"
	CMD

	commander::makecmd -a cmd2slice -s ' ' -c {COMMANDER[0]}<<- CMD
		split
		--numeric-suffixes=1
		-a ${#threads}
		-n l/$threads
		--additional-suffix=.bed
		--filter='bedtools merge -d 1 -i - > \$(dirname \$FILE)/\$(basename \$FILE | sed "s/^x0*/slice./")'
		"${tdirs[0]}/tmp" "${tdirs[0]}/x"
	CMD

	if $skip; then
		commander::printcmd -a cmd1slice
		commander::printcmd -a cmd2slice
	else
		commander::runcmd -v -b -t $threads -a cmd1slice
		commander::runcmd -v -b -t $threads -a cmd2slice
	fi


	local m i j f nf o t e slice odir tdir
	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6
	for m in "${_mapper_freebayes[@]}"; do
		declare -n _bams_freebayes=$m
		tdir="$tmpdir/$m"
		odir="$outdir/$m/freebayes"
		mkdir -p "$odir" "$tdir"

		for i in "${!_tidx_freebayes[@]}"; do
			f="${_bams_freebayes[${_tidx_freebayes[$i]}]}"
			o="$(basename "$f")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"
			tomerge=()

			for j in $(seq 1 $threads); do
				slice="$t.slice.$j"

				if [[ $_nidx_freebayes ]]; then # somatic
					nf="${_bams_freebayes[${_nidx_freebayes[$i]}]}"
					commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
						freebayes
							-f "$genome"
							-t "${tdirs[0]}/slice.$j.bed"
							--genotype-qualities
							--pooled-discrete
							--pooled-continuous
							--report-genotype-likelihood-max
							--allele-balance-priors-off
							--use-best-n-alleles 3
							--min-base-quality 20
							--min-coverage 10
							--min-alternate-fraction 0.1
							--harmonic-indel-quality
							--min-repeat-entropy 1
							--no-partial-observations
							-P 0.0001
							"$f" "$nf"
					CMD
						vcfsamplediff VCFSAMPLEDIFF NORMAL TUMOR - > "$slice.vcf"
					CMD
					# use to filter for germline risk
					# do not use --strict, which requires that no observation in the germline support the somatic alternate i.e. no 0/1 0/1
					# grep -E '(^#|VCFSAMPLEDIFF=somatic)' > "$slice.vcf"
					# do not hard filter here. vcfix.pl adds vcfix_somatic FILTER tag if somatic or loh and GQ MAF thresholds passed
				else
					commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
						freebayes
							-f "$genome"
							-t "${tdirs[0]}/slice.$j.bed"
							--genotype-qualities
							--use-best-n-alleles 3
							--min-base-quality 20
							--min-coverage 10
							--min-alternate-fraction 0.1
							--harmonic-indel-quality
							--min-repeat-entropy 1
							--no-partial-observations
							-P 0.0001
							"$f"
						> "$slice.vcf"
					CMD
				fi

				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				if [[ $dbsnp ]]; then
					commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl
					CMD
						vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$slice.fixed.nomulti.normed.vcf"
					CMD
				else
					commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
					CMD
				fi

				tomerge+=("$slice")
			done

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
				#DO NOT PIPE - DATALOSS!
				commander::makecmd -a cmd5 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$t.$e" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T "${tdirs[-1]}" -m ${memory}M -o "$o.$e" "$t.$e"
				CMD

				commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done
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
		commander::runcmd -c freebayes -v -b -t $threads -a cmd1
		commander::runcmd -v -b -t $threads -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
		commander::runcmd -v -b -t $threads -a cmd4
		commander::runcmd -v -b -t $threads -a cmd5
		commander::runcmd -v -b -t $instances -a cmd6
	fi

	return 0
}

variants::varscan() {
	declare -a tdirs
	_cleanup::variants::varscan(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-m <memory>    | per threads amount of (default: 2048)
			-M <maxmemory> | amount of
			-g <genome>    | path to
			-d <dbsnp>     | path to
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (for somatic calling. requieres -2)
			-2 <tumoridx>  | array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory maxmemory genome dbsnp tmpdir outdir
	declare -n _mapper_varscan _nidx_varscan _tidx_varscan
	while getopts 'S:s:t:M:g:d:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			m) memory=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((++mandatory)); _mapper_varscan=$OPTARG;;
			1) _nidx_varscan=$OPTARG;;
			2) _tidx_varscan=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "calling variants varscan"

	local m=${_mapper_varscan[0]}
	declare -n _bams_varscan=$m
	[[ ! $_tidx_varscan ]] && unset _tidx_varscan && declare -a _tidx_varscan=(${!_bams_varscan[@]}) # use all indices for germline calling

	local minstances mthreads jmem jgct jcgct ithreads i memory1 memory2 instances=$((${#_mapper_varscan[@]}*${#_tidx_varscan[@]}))
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m ${memory:-2048} -M "$maxmemory") # good estimate for varscan
	read -r i memory1 < <(configure::memory_by_instances -i $((minstances*instances)) -T $threads -M "$maxmemory") # for bcftools sort of initial vcf slices
	read -r i memory2 < <(configure::memory_by_instances -i $((instances*4)) -T $threads -M "$maxmemory") # for final bcftools sort of vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads) # for final bgzip

	tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed)")
	declare -a cmd1slice cmd2slice
	commander::makecmd -a cmd1slice -s ' ' -o "${tdirs[0]}/tmp" -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		perl -lane '
			$i=0;
			for (map {$_*100000} 1..($F[1]-1)/100000){
				$i=$_;
				print join"\t",($F[0],$_-100000,$_-1);
			}
			print join"\t",($F[0],$i,$F[1]-1) if $i<$F[1]-1;
		'
	CMD
		"$genome.fai"
	CMD

	commander::makecmd -a cmd2slice -s ' ' -c {COMMANDER[0]}<<- CMD
		split
		--numeric-suffixes=1
		-a ${#minstances}
		-n l/$minstances
		--additional-suffix=.bed
		--filter='bedtools merge -d 1 -i - > \$(dirname \$FILE)/\$(basename \$FILE | sed "s/^x0*/slice./")'
		"${tdirs[0]}/tmp" "${tdirs[0]}/x"
	CMD

	if $skip; then
		commander::printcmd -a cmd1slice
		commander::printcmd -a cmd2slice
	else
		commander::runcmd -v -b -t $threads -a cmd1slice
		commander::runcmd -v -b -t $threads -a cmd2slice
	fi


	local m i j f nf o t e slice odir tdir
	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6 cmd7 cmd8
	for m in "${_mapper_varscan[@]}"; do
		declare -n _bams_varscan=$m
		tdir="$tmpdir/$m"
		odir="$outdir/$m/varscan"
		mkdir -p "$odir" "$tdir"

		for i in "${!_tidx_varscan[@]}"; do
			f="${_bams_varscan[${_tidx_varscan[$i]}]}"
			o="$(basename "$f")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"
			tomerge=()

			for j in $(seq 1 $minstances); do
				slice="$t.slice.$j"

				if [[ $_nidx_varscan ]]; then # somatic
					nf="${_bams_varscan[${_nidx_varscan[$i]}]}"

					# print orphans A and overlaps x, since this can be adressed by the user via post-aln functions
					commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
						samtools mpileup
							-A
							-x
							-B
							-Q 20
							-d 100000
							--reference "$genome"
							-l "${tdirs[0]}/slice.$j.bed"
							"$nf" "$f"
						> "$slice.mpileup"
					CMD

					# pileup can be piped into varscan. but if disk is too busy, varscan simply stops and does not wait for input stream
					commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
						varscan
							-Xmx${jmem}m
							-XX:ParallelGCThreads=$jgct
							-XX:ConcGCThreads=$jcgct
							-Djava.io.tmpdir="$tmpdir"
						somatic "$slice.mpileup"
							--mpileup 1
							--min-coverage-normal 10
							--min-coverage-tumor 10
							--min-var-freq 0.1
							--min-freq-for-hom 0.7
							--strand-filter 1
							--p-value 1
							--somatic-p-value 0.05
							--output-snp "$slice.snp"
							--output-indel "$slice.indel"
							--output-vcf 1
					CMD

					tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
					# do not use bcftools concat -o "$slice.concat" "$slice.snp.reheader" "$slice.indel.reheader", since positions can overlap which is not handled by bcftools unless input is bgzip
					commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						bcftools reheader -f "$genome.fai" -o "$slice.snp.reheader" "$slice.snp.vcf"
					CMD
						bcftools reheader -f "$genome.fai" -o "$slice.indel.reheader" "$slice.indel.vcf"
					CMD
						cat <(bcftools view -h "$slice.snp.reheader") <(bcftools view -H "$slice.snp.reheader" | awk -v OFS='\t' '{\$7="."; print}') <(bcftools view -H "$slice.indel.reheader" | awk -v OFS='\t' '{\$7="."; print}') > "$slice.concat"
					CMD
						bcftools sort -T "${tdirs[-1]}" -m ${memory1}M -o "$slice.vcf" "$slice.concat"
					CMD
					#	grep -E '(^#|SS=2)' "$slice.unfiltered" > "$slice.vcf"
					#CMD # do not hard filter here. vcfix.pl adds vcfix_somatic FILTER tag if SS=2
				else
					commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
						samtools mpileup
							-A
							-x
							-B
							-Q 20
							-d 100000
							--reference "$genome"
							-l "${tdirs[0]}/slice.$j.bed"
							"$f"
						> "$slice.mpileup"
					CMD
					# pileup can be piped into varscan. but if disk is too busy, varscan simply stops and does not wait for input stream
					commander::makecmd -a cmd2 -s '|' -o "$slice.toreheader" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
						varscan
							-Xmx${jmem}m
							-XX:ParallelGCThreads=$jgct
							-XX:ConcGCThreads=$jcgct
							-Djava.io.tmpdir="$tmpdir"
						mpileup2cns "$slice.mpileup"
							--min-coverage 10
							--min-avg-qual 20
							--min-var-freq 0.1
							--min-freq-for-hom 0.7
							--p-value 1
							--strand-filter 1
							--output-vcf 1
							--variants 1
					CMD
						awk -v OFS='\t' '{if($7=="PASS"){$7="."}; print}'
					CMD
					# all sites are labeled with PASS

					# do not pipe into bcftools : --fai cannot be used when reading from stdin
					commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
						bcftools reheader -f "$genome.fai" -o "$slice.vcf" "$slice.toreheader"
					CMD
				fi

				commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				if [[ $dbsnp ]]; then
					commander::makecmd -a cmd6 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl
					CMD
						vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$slice.fixed.nomulti.normed.vcf"
					CMD
				else
					commander::makecmd -a cmd6 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
					CMD
				fi

				tomerge+=("$slice")
			done

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
				#DO NOT PIPE - DATALOSS!
				commander::makecmd -a cmd7 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$t.$e" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T "${tdirs[-1]}" -m ${memory2}M -o "$o.$e" "$t.$e"
				CMD

				commander::makecmd -a cmd8 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
		commander::printcmd -a cmd6
		commander::printcmd -a cmd7
		commander::printcmd -a cmd8
	else
		commander::runcmd -v -b -t $threads -a cmd1
		commander::runcmd -c varscan -v -b -t $minstances -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
		commander::runcmd -v -b -t $threads -a cmd4
		commander::runcmd -v -b -t $threads -a cmd5
		commander::runcmd -v -b -t $threads -a cmd6
		commander::runcmd -v -b -t $threads -a cmd7
		commander::runcmd -v -b -t $instances -a cmd8
	fi

	return 0
}

variants::vardict() {
	declare -a tdirs
	_cleanup::variants::vardict(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-m <memory>    | per thread amount of (default: 12288)
			-M <maxmemory> | amount of
			-g <genome>    | path to
			-d <dbsnp>     | path to
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (for somatic calling. requieres -2)
			-2 <tumoridx>  | array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory maxmemory genome dbsnp tmpdir outdir
	declare -n _mapper_vardict _nidx_vardict _tidx_vardict
	while getopts 'S:s:t:M:g:d:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			m) memory=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((++mandatory)); _mapper_vardict=$OPTARG;;
			1) _nidx_vardict=$OPTARG;;
			2) _tidx_vardict=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "calling variants vardict"

	local m=${_mapper_vardict[0]}
	declare -n _bams_vardict=$m
	[[ ! $_tidx_vardict ]] && unset _tidx_vardict && declare -a _tidx_vardict=(${!_bams_vardict[@]}) # use all indices for germline calling

	local minstances mthreads jmem jgct jcgct ithreads i memory1 memory2 instances=$((${#_mapper_vardict[@]}*${#_tidx_vardict[@]}))
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m ${memory:=12288} -M "$maxmemory")
	# I often estimate needed memory by about 1Mb for 1 bp of region (here are included possible intermediate structures for variations, reference and softclips).
	# Also you have to take account of amount of threads that you run, because memory will be used by all of them and will be split between threads
	# https://github.com/AstraZeneca-NGS/VarDictJava/issues/216
	# -> 10kb chunks -> 10G per thread + JVM overhead by GC
	read -r i memory1 < <(configure::memory_by_instances -i $((minstances*instances)) -T $threads -M "$maxmemory") # for bcftools sort of initial vcf slices
	read -r i memory2 < <(configure::memory_by_instances -i $((instances*4)) -T $threads -M "$maxmemory") # for final bcftools sort of vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads) # for final bgzip

	tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed)")
	declare -a cmd1slice cmd2slice
	commander::makecmd -a cmd1slice -s ' ' -o "${tdirs[0]}/tmp" -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		perl -M'List::Util qw(max)' -slane '
			$l = $m-2048 >= 9500 ? 10000 : int(($m-1024)/1000)*1000;
			$i=0;
			for (map {$_*$l} 1..($F[1]-1)/$l){
				$i=$_;
				print join"\t",($F[0],max(0,$_-$l-100),$_-1);
			}
			print join"\t",($F[0],$i-100,$F[1]-1) if $i<$F[1]-1;
		'
	CMD
		-- -m=$jmem "$genome.fai"
	CMD
	# despite of mutlithreading support, execute in parallel on multiple chunks, to limit memory footprint which heavily depends on region size and coverage.
	# small chunks within one bed file wille be individually processed. thus a small overlap is recommended to no miss any variant.
	# overlaps will by internally handled -> https://github.com/AstraZeneca-NGS/VarDictJava/issues/64

	commander::makecmd -a cmd2slice -s ' ' -c {COMMANDER[0]}<<- CMD
		split
		--numeric-suffixes=1
		-a ${#minstances}
		-n l/$minstances
		--additional-suffix=.bed
		--filter='cat > \$(dirname \$FILE)/\$(basename \$FILE | sed "s/^x0*/slice./")'
		"${tdirs[0]}/tmp" "${tdirs[0]}/x"
	CMD

	if $skip; then
		commander::printcmd -a cmd1slice
		commander::printcmd -a cmd2slice
	else
		commander::runcmd -v -b -t $threads -a cmd1slice
		commander::runcmd -v -b -t $threads -a cmd2slice
	fi


	local m i j f nf o t e slice odir tdir
	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6 cmd7
	for m in "${_mapper_vardict[@]}"; do
		declare -n _bams_vardict=$m
		tdir="$tmpdir/$m"
		odir="$outdir/$m/vardict"
		mkdir -p "$odir" "$tdir"

		for i in "${!_tidx_vardict[@]}"; do
			f="${_bams_vardict[${_tidx_vardict[$i]}]}"
			o="$(basename "$f")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"
			tomerge=()

			for j in $(seq 1 $minstances); do
				slice="$t.slice.$j"

				if [[ $_nidx_vardict ]]; then # somatic
					nf="${_bams_vardict[${_nidx_vardict[$i]}]}"
					# prefer -fisher option (vardict skips corrupt sites) over vardict | testsomatic.R | var2vcf_paired.pl , due to R: Error in fisher.test
					# use -X 0 to not combine snps within 3 bp window into indel or complex, which in addiation reduces the risk of StringIndexOutOfBoundsException
					# vardict will fail upon > 10 errors per region
					commander::makecmd -a cmd1 -s '|' -o "$slice.vcf" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD'
						JAVA_OPTS="-Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.tmpdir='$tmpdir' -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
						vardict-java
							-G "$genome"
							-f 0.1
							-V 0.05
							-F 0x200
							-N TUMOR
							-th 1
							-q 20
							-U
							-X 0
							-fisher
							-v
							-c 1
							-S 2
							-E 3
							-b "$f|$nf"
							"${tdirs[0]}/slice.$j.bed"
					CMD
						var2vcf_paired.pl
							-q 20
							-Q 0
							-d 10
							-v 3
							-F 0.3
							-f 0.1
							-M
							-A
							-N "TUMOR|NORMAL"
					CMD
						awk -v OFS='\t' '{if(/^#CHROM/){s=1}; if(s){t=$11; $11=$10; $10=t}; print}'
					CMD
					#	grep -E '(^#|STATUS=[^;]*[sS]omatic)' > "$slice.vcf"
					#CMD # do not hard filter here. vcfix.pl adds vcfix_somatic FILTER tag if STATUS=[^;]*[sS]omatic
				else
					# prefer -fisher option (vardict skips corrupt sites) over vardict | teststrandbias.R | var2vcf_valid.pl , due to R: Error in fisher.test
					commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
						JAVA_OPTS="-Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.tmpdir='$tmpdir' -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
						vardict-java
						-G "$genome"
						-f 0.1
						-F 0x200
						-N SAMPLE
						-th 1
						-q 20
						-U
						-X 0
						-fisher
						-v
						-c 1
						-S 2
						-E 3
						-b "$f"
						"${tdirs[0]}/slice.$j.bed"
					CMD
						var2vcf_valid.pl
							-q 20
							-Q 0
							-d 10
							-v 2
							-F 0.3
							-f 0.1
							-A
							-E
							-N SAMPLE
						> "$slice.vcf"
					CMD
				fi
				# -A print all variants at the same site. otherwise reported one is chosen by MAF.
				# note: this does not always mean multiple alleles per se (like bcftools norm), but sometimes duplicons with slightly differnt format and info fields

				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools reheader -f "$genome.fai" -o "$slice.reheader" "$slice.vcf"
				CMD
					bcftools sort -T "${tdirs[-1]}" -m ${memory1}M -o "$slice.vcf" "$slice.reheader"
				CMD

				commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				if [[ $dbsnp ]]; then
					commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl
					CMD
						vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$slice.fixed.nomulti.normed.vcf"
					CMD
				else
					commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
						vcffixup "$slice.fixed.nomulti.vcf"
					CMD
						vt normalize -q -n -r "$genome" -
					CMD
						vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
					CMD
				fi

				tomerge+=("$slice")
			done

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
				#DO NOT PIPE - DATALOSS!
				commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$t.$e" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T "${tdirs[-1]}" -m ${memory2}M -o "$o.$e" "$t.$e"
				CMD

				commander::makecmd -a cmd7 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
		commander::printcmd -a cmd6
		commander::printcmd -a cmd7
	else
		commander::runcmd -c vardict -v -b -t $minstances -a cmd1
		commander::runcmd -v -b -t $threads -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
		commander::runcmd -v -b -t $threads -a cmd4
		commander::runcmd -v -b -t $threads -a cmd5
		commander::runcmd -v -b -t $threads -a cmd6
		commander::runcmd -v -b -t $instances -a cmd7
	fi

	return 0
}

variants::vardict_threads() {
	# may be used as fallback to reduce disk i/o
	# use one region file and multiple threads on it
	# -> something limits parallelization to ~14 threads. chunks may not be the reason: tested with up to 100k
	# memory footprint raises with multiple threads
	declare -a tdirs
	_cleanup::variants::vardict(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-M <maxmemory> | amount of
			-g <genome>    | path to
			-d <dbsnp>     | path to
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (for somatic calling. requieres -2)
			-2 <tumoridx>  | array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory genome dbsnp tmpdir outdir
	declare -n _mapper_vardict _nidx_vardict _tidx_vardict
	while getopts 'S:s:t:M:g:d:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((++mandatory)); _mapper_vardict=$OPTARG;;
			1) _nidx_vardict=$OPTARG;;
			2) _tidx_vardict=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "calling variants vardict"

	local m=${_mapper_vardict[0]}
	declare -n _bams_vardict=$m
	[[ ! $_tidx_vardict ]] && unset _tidx_vardict && declare -a _tidx_vardict=(${!_bams_vardict[@]}) # use all indices for germline calling

	local minstances mthreads jmem jgct jcgct ithreads i memory instances=$((${#_mapper_vardict[@]}*${#_tidx_vardict[@]}))
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i 1 -T $threads -M "$maxmemory")
	read -r i memory < <(configure::memory_by_instances -i $instances -T $threads -M "$maxmemory") # for final bcftools sort of vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads) # for final bgzip

	tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bed)")
	declare -a cmd1slice
	commander::makecmd -a cmd1slice -s ' ' -o "${tdirs[0]}/regions.bed" -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		perl -M'List::Util qw(max)' -lane '
			$i=0;
			for (map {$_*10000} 1..($F[1]-1)/10000){
				$i=$_;
				print join"\t",($F[0],max(0,$_-10000-100),$_-1);
			}
			print join"\t",($F[0],$i-100,$F[1]-1) if $i<$F[1]-1;
		'
	CMD
		"$genome.fai"
	CMD
	# TODO either do chunks depending on jmem/threads or restrict -th parameter to maxmem/jmem

	if $skip; then
		commander::printcmd -a cmd1slice
	else
		commander::runcmd -v -b -t $threads -a cmd1slice
	fi

	local m i f nf o t e odir tdir
	declare -a cmd1 cmd2 cmd3 cmd4 cmd5 cmd6
	for m in "${_mapper_vardict[@]}"; do
		declare -n _bams_vardict=$m
		tdir="$tmpdir/$m"
		odir="$outdir/$m/vardict"
		mkdir -p "$odir" "$tdir"

		for i in "${!_tidx_vardict[@]}"; do
			f="${_bams_vardict[${_tidx_vardict[$i]}]}"
			o="$(basename "$f")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"
			tomerge=()

			if [[ $_nidx_vardict ]]; then # somatic
				nf="${_bams_vardict[${_nidx_vardict[$i]}]}"
				# prefer -fisher option (vardict skips corrupt sites) over vardict | testsomatic.R | var2vcf_paired.pl , due to R: Error in fisher.test
				# use -X 0 to not combine snps within 3 bp window into indel or complex, which in addiation reduces the risk of StringIndexOutOfBoundsException
				# vardict will fail upon > 10 errors per region
				commander::makecmd -a cmd1 -s '|' -o "$t.toreheader" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD'
					JAVA_OPTS="-Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.tmpdir='$tmpdir' -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
					vardict-java
						-G $genome
						-f 0.1
						-V 0.05
						-F 0x200
						-N TUMOR
						-th $threads
						-q 20
						-U
						-X 0
						-fisher
						-v
						-c 1
						-S 2
						-E 3
						-b "$f|$nf"
						"${tdirs[0]}/regions.bed"
				CMD
					var2vcf_paired.pl
						-q 20
						-Q 0
						-d 10
						-v 3
						-F 0.3
						-f 0.1
						-M
						-A
						-N "TUMOR|NORMAL"
				CMD
					awk -v OFS='\t' '{if(/^#CHROM/){s=1}; if(s){t=$11; $11=$10; $10=t}; print}'
				CMD
				#	grep -E '(^#|STATUS=[^;]*[sS]omatic)' > "$slice.vcf"
				#CMD # do not hard filter here. vcfix.pl adds vcfix_somatic FILTER tag if STATUS=[^;]*[sS]omatic
			else
				# prefer -fisher option (vardict skips corrupt sites) over vardict | teststrandbias.R | var2vcf_valid.pl , due to R: Error in fisher.test
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					JAVA_OPTS="-Xmx${jmem}m -XX:ParallelGCThreads=$jgct -XX:ConcGCThreads=$jcgct -Djava.io.tmpdir='$tmpdir' -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
					vardict-java
					-G $genome
					-f 0.1
					-F 0x200
					-N SAMPLE
					-th $threads
					-q 20
					-U
					-X 0
					-fisher
					-v
					-c 1
					-S 2
					-E 3
					-b "$f"
					"${tdirs[0]}/regions.bed"
				CMD
					var2vcf_valid.pl
						-q 20
						-Q 0
						-d 10
						-v 2
						-F 0.3
						-f 0.1
						-A
						-E
						-N SAMPLE
					> "$t.toreheader"
				CMD
			fi
			# -A print all variants at the same site. otherwise reported one is chosen by MAF.
			# note: this does not always mean multiple alleles per se (like bcftools norm), but sometimes duplicons with slightly differnt format and info fields

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bcftools reheader -f "$genome.fai" -o "$t.vcf" "$t.toreheader"
			CMD
				bcftools sort -T "${tdirs[-1]}" -m ${memory}M -o "$o.vcf" "$t.vcf"
			CMD

			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
				vcfix.pl -i "$o.vcf" > "$o.fixed.vcf"
			CMD

			commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bcftools norm -f "$genome" -c s -m-both "$o.fixed.vcf"
			CMD
				vcfix.pl -i - > "$o.fixed.nomulti.vcf"
			CMD

			if [[ $dbsnp ]]; then
				commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					vcffixup "$o.fixed.nomulti.vcf"
				CMD
					vt normalize -q -n -r "$genome" -
				CMD
					vcfixuniq.pl
				CMD
					vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$o.fixed.nomulti.normed.vcf"
				CMD
			else
				commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					vcffixup "$o.fixed.nomulti.vcf"
				CMD
					vt normalize -q -n -r "$genome" -
				CMD
					vcfixuniq.pl > "$o.fixed.nomulti.normed.vcf"
				CMD
			fi

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done
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
		commander::runcmd -c vardict -v -b -t 1 -a cmd1
		commander::runcmd -v -b -t $threads -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
		commander::runcmd -v -b -t $threads -a cmd4
		commander::runcmd -v -b -t $threads -a cmd5
		commander::runcmd -v -b -t $instances -a cmd6
	fi

	return 0
}

variants::platypus() {
	declare -a tdirs
	_cleanup::variants::platypus(){
		rm -rf "${tdirs[@]}"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-M <maxmemory> | amount of
			-g <genome>    | path to
			-d <dbsnp>     | path to
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (for somatic calling. requieres -2)
			-2 <tumoridx>  | array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory genome dbsnp tmpdir outdir
	declare -n _mapper_platypus _nidx_platypus _tidx_platypus
	while getopts 'S:s:t:M:g:d:m:r:1:2:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			M) maxmemory="$OPTARG";;
			g) ((++mandatory)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((++mandatory)); _mapper_platypus=$OPTARG;;
			1) _nidx_platypus=$OPTARG;;
			2) _tidx_platypus=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "calling variants platypus"

	local m=${_mapper_platypus[0]}
	declare -n _bams_platypus=$m
	[[ ! $_tidx_platypus ]] && unset _tidx_platypus && declare -a _tidx_platypus=(${!_bams_platypus[@]}) # use all indices for germline calling

	local ithreads i memory instances=$((${#_mapper_platypus[@]}*${#_tidx_platypus[@]}))
	read -r i memory < <(configure::memory_by_instances -i $((instances*4)) -T $threads -M "$maxmemory") # for final bcftools sort of vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 10 -T $threads) # for final bgzip

	local m i f nf o t e odir tdir
	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6
	for m in "${_mapper_platypus[@]}"; do
		declare -n _bams_platypus=$m
		tdir="$tmpdir/$m"
		odir="$outdir/$m/platypus"
		mkdir -p "$odir" "$tdir"

		for i in "${!_tidx_platypus[@]}"; do
			f="${_bams_platypus[${_tidx_platypus[$i]}]}"
			o="$(basename "$f")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"

			# can make use of --regions=$(awk '{print $1":"$2"-"$3}' $bed | xargs -echo | sed 's/ /,/g')
			if [[ $_nidx_platypus ]]; then # somatic
				nf="${_bams_platypus[${_nidx_platypus[$i]}]}"
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					platypus callVariants
						--logFileName /dev/null
						--output "$t.toreheader"
						--refFile "$genome"
						--bamFiles "$nf","$f"
						--minGoodQualBases 0
						--minBaseQual 20
						--minMapQual 0
						--minReads 10
						--nCPU $threads
						--filterVarsByCoverage 1
						--minVarFreq 0.1
						--filterDuplicates 0
						--useEMLikelihoods 1
						--minPosterior 3
						--verbosity 1
						--filterDuplicates 0
				CMD
				# cannot handle --output >(bcftools > out.vcf)

				# use to filter for germline risk
				# do not use --strict, which requires that no observation in the germline support the somatic alternate i.e. no 0/1 0/1
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools reheader -f "$genome.fai" "$t.toreheader"
				CMD
					vcfsamplediff VCFSAMPLEDIFF NORMAL TUMOR - > "$o.vcf"
				CMD
				#	grep -E '(^#|VCFSAMPLEDIFF=somatic)' > "$o.vcf"
				#CMD # do not hard filter here. vcfix.pl adds vcfix_somatic FILTER tag if somatic or loh and GQ MAF thresholds passed
			else
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					platypus callVariants
						--logFileName /dev/null
						--output "$t.toreheader"
						--refFile "$genome"
						--bamFiles "$f"
						--minGoodQualBases 0
						--minBaseQual 20
						--minMapQual 0
						--minReads 10
						--nCPU $threads
						--filterVarsByCoverage 1
						--minVarFreq 0.1
						--useEMLikelihoods 1
						--minPosterior 3
						--verbosity 1
						--filterDuplicates 0
				CMD

				commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
					bcftools reheader -f "$genome.fai" -o "$o.vcf" "$t.toreheader"
				CMD
			fi

			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
				vcfix.pl -i "$o.vcf" > "$o.fixed.vcf"
			CMD

			commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bcftools norm -f "$genome" -c s -m-both "$o.fixed.vcf"
			CMD
				vcfix.pl -i - > "$o.fixed.nomulti.vcf"
			CMD

			if [[ $dbsnp ]]; then
				commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					vcffixup "$o.fixed.nomulti.vcf"
				CMD
					vt normalize -q -n -r "$genome" -
				CMD
					vcfixuniq.pl
				CMD
					vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$o.fixed.nomulti.normed.vcf"
				CMD
			else
				commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					vcffixup "$o.fixed.nomulti.vcf"
				CMD
					vt normalize -q -n -r "$genome" -
				CMD
					vcfixuniq.pl > "$o.fixed.nomulti.normed.vcf"
				CMD
			fi

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				commander::makecmd -a cmd6 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done
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
		commander::runcmd -c platypus -v -b -t 1 -a cmd1
		commander::runcmd -v -b -t $threads -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
		commander::runcmd -v -b -t $threads -a cmd4
		commander::runcmd -v -b -t $threads -a cmd5
		commander::runcmd -v -b -t $instances -a cmd6
	fi

	return 0
}
