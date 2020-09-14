#! /usr/bin/env bash
# (c) Konstantin Riege

variants::vcfzip() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-z <var>      | of path to file
			example:
			$funcname -t 4 -v f1 -v f2
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads
	declare -a tozip_vcfzip
	while getopts 'S:s:t:z:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			z) ((mandatory++)); tozip_vcfzip+=("$OPTARG");;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage && return 1

	commander::printinfo "compressing and indexing vcf"

	declare -a cmd1 cmd2
	local f
	for f in "${tozip_vcfzip[@]}"; do
		declare -n _f_vcfzip=$f
		readlink -e "$_f_vcfzip" | file -f - | grep -qF 'compressed' || {
			commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
				bgzip -f -@ $threads < "$_f_vcfzip" > "$_f_vcfzip.gz"
			CMD
			_f_vcfzip="$_f_vcfzip.gz"
		}
		commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
			tabix -f -p vcf "$_f_vcfzip"
		CMD
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	commander::runcmd -v -b -t 1 -a cmd1 && \
			commander::runcmd -v -b -t $threads -a cmd2
		} || {
			commander::printerr "$funcname failed"
			return 1
		}
	}

	return 0
}

variants::haplotypecaller() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
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
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory genome dbsnp tmpdir outdir i
	declare -n _mapper_haplotypecaller _bamslices_haplotypecaller
	declare -A nidx tidx
	while getopts 'S:s:t:g:d:m:r:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((mandatory++)); _mapper_haplotypecaller=$OPTARG;;
			c) ((mandatory++)); _bamslices_haplotypecaller=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			o) ((mandatory++)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1
	if [[ ! $dbsnp ]]; then
		dbsnp="$tmpdir/$(basename "$genome").vcf"
		echo -e "##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > "$dbsnp"
		bgzip -f -@ $threads < "$dbsnp" > "$dbsnp.gz"
		tabix -f -p vcf "$dbsnp.gz"
		dbsnp="$dbsnp.gz"
	fi

	commander::printinfo "calling variants haplotypecaller"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o t e slice odir instances ithreads odir tdir
	for m in "${_mapper_haplotypecaller[@]}"; do
		declare -n _bams_haplotypecaller=$m
		((instances+=${#_bams_haplotypecaller[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 1 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6 tdirs
	for m in "${_mapper_haplotypecaller[@]}"; do
		declare -n _bams_haplotypecaller=$m
		tdir="$tmpdir/$m"
		odir="$outdir/$m"
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
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
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
						-O "$slice.vcf"
						-R "$genome"
						-D "$dbsnp"
						-A Coverage
						-A DepthPerAlleleBySample
						-A StrandBiasBySample
						--min-base-quality-score 20
						--native-pair-hmm-threads $mthreads
						--max-alternate-alleles 3
						--all-site-pls true
						-verbosity INFO
						--tmp-dir "${tdirs[-1]}"
				CMD

				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					vcffixup "$slice.fixed.nomulti.vcf"
				CMD
					vt normalize -q -n -r "$genome" -
				CMD
					vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
				CMD

				tomerge+=("$slice")
			done < "${_bamslices_haplotypecaller[${_bams_haplotypecaller[$i]}]}"

			o=$(basename "${_bams_haplotypecaller[$i]}")
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"

			for e in vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
				#DO NOT PIPE - DATALOSS!
				commander::makecmd -a cmd5 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$t.$e" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T "${tdirs[-1]}" -m ${memory}M -o "$o.$e" "$t.$e"
				CMD

				commander::makecmd -a cmd6 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
		commander::printcmd -a cmd5
		commander::printcmd -a cmd6
	} || {
		{	commander::runcmd -c gatk -v -b -t $minstances -a cmd1 && \
			commander::runcmd -v -b -t $threads -a cmd2 && \
			commander::runcmd -v -b -t $threads -a cmd3 && \
			commander::runcmd -v -b -t $threads -a cmd4 && \
			commander::runcmd -v -b -t $minstances -a cmd5 && \
			commander::runcmd -v -b -t $instances -a cmd6
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}

variants::panelofnormals() {
# The panel of normals not only represents common germline variant sites,
# it presents commonly noisy sites in sequencing data, e.g. mapping artifacts or
# other somewhat random but systematic artifacts of sequencing.
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
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
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory genome dbsnp tmpdir outdir i
	declare -n _mapper_panelofnormals _bamslices_panelofnormals
	while getopts 'S:s:t:g:m:r:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			r) ((mandatory++)); _mapper_panelofnormals=$OPTARG;;
			c) ((mandatory++)); _bamslices_panelofnormals=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			o) ((mandatory++)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage && return 1

	commander::printinfo "calling panel of normals"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)

	local m i o t e slice odir tdir
	declare -a tomerge cmd1 cmd2 cmd3 tdirs
	for m in "${_mapper_panelofnormals[@]}"; do
		declare -n _bams_panelofnormals=$m
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		for i in "${!_bams_panelofnormals[@]}"; do
			tomerge=()

			while read -r slice; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
					gatk
						--java-options '
								-Xmx${jmem}m
								-XX:ParallelGCThreads=$jgct
								-XX:ConcGCThreads=$jcgct
								-Djava.io.tmpdir="$tmpdir"
								-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
							'
						Mutect2
						-I "$slice"
						-O "$slice.pon.vcf"
						-R "$genome"
						-A Coverage
						-A DepthPerAlleleBySample
						-A StrandBiasBySample
						--min-base-quality-score 20
						--native-pair-hmm-threads $mthreads
						-verbosity INFO
						--tmp-dir "${tdirs[-1]}"
						--max-mnp-distance 0
				CMD

				tomerge+=("$slice.pon.vcf")
			done < "${_bamslices_panelofnormals[${_bams_panelofnormals[$i]}]}"

			o="$(basename "${_bams_panelofnormals[$i]}")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
			#DO NOT PIPE - DATALOSS!
			commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bcftools concat -o "$t.pon.vcf" $(printf '"%s" ' "${tomerge[@]}")
			CMD
				bcftools sort -T "${tdirs[-1]}" -m ${memory}M -o "$o.pon.vcf" "$t.pon.vcf"
			CMD
		done
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
	} || {
		{	commander::runcmd -c gatk -v -b -t $minstances -a cmd1 && \
			commander::runcmd -v -b -t $threads -a cmd2
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}

variants::makepondb() {
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-r <mapper>    | array of sorted, indexed bams within array of
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory genome dbsnp tmpdir outdir i
	declare -n _mapper_makepondb
	while getopts 'S:s:t:g:r:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			r) ((mandatory++)); _mapper_makepondb=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			o) ((mandatory++)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage && return 1

	commander::printinfo "calling panel of normals"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m 4000)

	local m i o odir params instances ithreads
	for m in "${_mapper_makepondb[@]}"; do
		declare -n _bams_makepondb=$m
		((instances+=${#_bams_makepondb[@]}))
	done
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 1 -T $threads)

	local m i o t odir params
	declare -a tomerge cmd1 cmd2 cmd3 tdirs
	for m in "${_mapper_makepondb[@]}"; do
		declare -n _bams_makepondb=$m
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		tomerge=()
		for i in "${!_bams_makepondb[@]}"; do
			o="$(basename "${_bams_makepondb[$i]}")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"
			[[ ! -s "$o.pon.vcf" ]] && commander::printerr "file does not exists $o.pon.vcf" && return 1

			commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[4]}<<- CMD
				bcftools reheader -s <(echo "NORMAL NORMAL$i") -o "$t.pon.vcf" "$o.pon.vcf"
			CMD
				mv "$t.pon.vcf" "$o.pon.vcf"
			CMD
				bgzip -f -@ $ithreads < "$o.pon.vcf" > "$o.pon.vcf.gz"
			CMD
				tabix -f -p vcf "$o.pon.vcf.gz"
			CMD

			tomerge+=("$o.pon.vcf.gz")
		done

		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
		commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			rm -rf "$odir/pondb"
		CMD
			gatk
				--java-options '
					-Xmx${jmem}m
					-XX:ParallelGCThreads=$jgct
					-XX:ConcGCThreads=$jcgct
					-Djava.io.tmpdir="$tmpdir"
					-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
				'
				GenomicsDBImport
				-L "$genome.list"
				-R "$genome"
				$(printf ' -V "%s" ' "${tomerge[@]}")
				--genomicsdb-workspace-path "$odir/pondb"
				--verbosity INFO
				--tmp-dir "${tdirs[-1]}"
		CMD
		#rm pondb becase this does not work yet: --overwrite-existing-genomicsdb-workspace true

		[[ -s "$genome.af_only_gnomad.vcf.gz" ]] && params=" --germline-resource '$genome.af_only_gnomad.vcf.gz'" || params=''
		commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
			gatk
				--java-options '
					-Xmx${jmem}m
					-XX:ParallelGCThreads=$jgct
					-XX:ConcGCThreads=$jcgct
					-Djava.io.tmpdir="$tmpdir"
					-DGATK_STACKTRACE_ON_USER_EXCEPTION=true
				'
				CreateSomaticPanelOfNormals
				$params
				-V "gendb://$odir/pondb"
				-R "$genome"
				-O "$odir/pon.vcf.gz"
				--verbosity INFO
				--tmp-dir "${tdirs[-1]}"
		CMD
	done

	$skip && {
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd2
	} || {
		{	commander::runcmd -v -b -t $instances -a cmd1 && \
			commander::runcmd -c gatk -v -b -t $minstances -a cmd2 && \
			commander::runcmd -c gatk -v -b -t $minstances -a cmd3
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}

variants::mutect() {
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

	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			$funcname usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-m <memory>    | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-1 <normalidx> | array of (requieres -2)
			-2 <tumoridx>  | array of
			-c <sliceinfo> | array of
			-d <mypon>     | true/false
			-p <tmpdir>    | path to
			-o <outbase>   | path to
		EOF
		return 0
	}

	local OPTIND arg mandatory skip=false threads memory genome tmpdir outdir i mypon=false
	declare -n _mapper_mutect _bamslices_mutect _nidx_mutect _tidx_mutect
	while getopts 'S:s:t:g:d:m:r:1:2:c:d:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((mandatory++)); threads=$OPTARG;;
			m) ((mandatory++)); memory=$OPTARG;;
			g) ((mandatory++)); genome="$OPTARG";;
			r) ((mandatory++)); _mapper_mutect=$OPTARG;;
			c) ((mandatory++)); _bamslices_mutect=$OPTARG;;
			d) mypon=$OPTARG;;
			1) ((mandatory++)); _nidx_mutect=$OPTARG;;
			2) ((mandatory++)); _tidx_mutect=$OPTARG;;
			p) ((mandatory++)); tmpdir="$OPTARG"; mkdir -p "$tmpdir" || return 1;;
			o) ((mandatory++)); outdir="$OPTARG"; mkdir -p "$outdir" || return 1;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 9 ]] && _usage && return 1

	commander::printinfo "calling variants mutect"

	local minstances mthreads jmem jgct jcgct minstances2 mthreads2 jmem2 jgct2 jcgct2
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory)
	read -r minstances2 mthreads2 jmem2 jgct2 jcgct2 < <(configure::jvm -T $threads -m 4000)

	local params params2 m i o t slice odir tdir ithreads instances=$((${#_mapper_mutect[@]}*${#_tidx_mutect[@]}))
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -t 1 -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3 cmd4 cmd5 cmd6 cmd7 cmd8 cmd9 cmd10 tdirs
	for m in "${_mapper_mutect[@]}"; do
		declare -n _bams_mutect=$m
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		$mypon && {
			# TODO later on the fly ? needs to be completed by pon calling,
			# but sice pon creation is germline calling i'd like to keep things seperated
			# {	variants::makepondb \
			# 		-s $skip \
			# 		-t $threads \
			# 		-g "$genome" \
			# 		-r _mapper_mutect \
			# 		-1 _nidx_mutect \ <- requiered
			# 		-p "$tmpdir" \
			# 		-o "$outdir"
			# } || return 1
			[[ ! -s "$odir/pon.vcf.gz" ]] && commander::printerr "cannot find panel of normals $odir/pon.vcf.gz" && return 1
			params=" -pon '$odir/pon.vcf.gz'"
		} || {
			if [[ -s "$genome.pon.vcf.gz" ]]; then
				params=" -pon '$genome.pon.vcf.gz'"
			else
				params=''
				commander::warn "proceeding without panel of normals - expect file $genome.pon.vcf.gz"
			fi
		}

		for i in "${!_tidx_mutect[@]}"; do
			tomerge=()

			while read -r nslice slice; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
				# normal name defined for RGSM sam header entry by alignment::addreadgroup
				commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD
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
						--min-base-quality-score 20
						--native-pair-hmm-threads $mthreads
						--verbosity INFO
						--tmp-dir "${tdirs[-1]}"
						--max-mnp-distance 0
				CMD
				# --max-mnp-distance 0
				# default:1 - set to 0 to list adjacent SNPs as single entries rather than an MNP
				# useful for my SNP based genotype tree inference

				if [[ -s "$genome.somatic_common.vcf.gz" ]]; then
					commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD
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
							-V "$genome.somatic_common.vcf.gz"
							-L "$genome.somatic_common.vcf.gz"
							-O "$slice.pileups.table"
							--verbosity INFO
							--tmp-dir "${tdirs[-1]}"
					CMD


					##########
					#TODO awk 'NR>2' *slice?.bam.pileups.table | sort -k1,1V -k2,2n  BUT -k1,1 must be in order of $genome.list
					# then predict contamination as whole and filter on full set
					##########
					commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD
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
				else
					params2=''
				fi

				commander::makecmd -a cmd4 -s '|' -c {COMMANDER[0]}<<- CMD
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
						-R "$genome"
						-V "$slice.unfiltered.vcf"
						-stats "$slice.unfiltered.vcf.stats"
						-O "$slice.vcf"
						--verbosity INFO
						--tmp-dir "${tdirs[-1]}"
				CMD
				#In order to tweak results in favor of more sensitivity users may set -f-score-beta to a value greater than its default of 1 (beta is the relative weight of sensitivity versus precision in the harmonic mean). Setting it lower biases results toward greater precision.
				# optional: --tumor-segmentation segments.table

				commander::makecmd -a cmd5 -s '|' -c {COMMANDER[0]}<<- CMD
					vcfix.pl -i "$slice.vcf" > "$slice.fixed.vcf"
				CMD

				commander::makecmd -a cmd6 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools norm -f "$genome" -c s -m-both "$slice.fixed.vcf"
				CMD
					vcfix.pl -i - > "$slice.fixed.nomulti.vcf"
				CMD

				commander::makecmd -a cmd7 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					vcffixup "$slice.fixed.nomulti.vcf"
				CMD
					vt normalize -q -n -r "$genome" -
				CMD
					vcfixuniq.pl > "$slice.fixed.nomulti.normed.vcf"
				CMD

				tomerge+=("$slice")
			done < <(paste "${_bamslices_mutect[${_bams_mutect[${_nidx_mutect[$i]}]}]}" "${_bamslices_mutect[${_bams_mutect[${_tidx_mutect[$i]}]}]}")

			o="$(basename "${_bams_mutect[${_tidx_mutect[$i]}]}")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
			commander::makecmd -a cmd8 -s '|' -c {COMMANDER[0]}<<- CMD
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
				commander::makecmd -a cmd9 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools concat -o "$t.$e" $(printf '"%s" ' "${tomerge[@]/%/.$e}")
				CMD
					bcftools sort -T "${tdirs[-1]}" -m ${memory}M -o "$o.$e" "$t.$e"
				CMD

				commander::makecmd -a cmd10 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
			done

		done
	done

	$skip && {
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
	} || {
		{	commander::runcmd -c gatk -v -b -t $minstances -a cmd1 && \
			commander::runcmd -c gatk -v -b -t $minstances -a cmd2 && \
			commander::runcmd -c gatk -v -b -t $minstances -a cmd3 && \
			commander::runcmd -c gatk -v -b -t $minstances -a cmd4 && \
			commander::runcmd -v -b -t $threads -a cmd5 && \
			commander::runcmd -v -b -t $threads -a cmd6 && \
			commander::runcmd -v -b -t $threads -a cmd7 && \
			commander::runcmd -c gatk -v -b -t $minstances2 -a cmd8 && \
			commander::runcmd -v -b -t $instances -a cmd9 && \
			commander::runcmd -v -b -t $instances -a cmd10
		} || {
			rm -rf "${tdirs[@]}"
			commander::printerr "$funcname failed"
			return 1
		}
	}

	rm -rf "${tdirs[@]}"
	return 0
}
