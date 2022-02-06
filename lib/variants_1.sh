#! /usr/bin/env bash
# (c) Konstantin Riege

variants::vcfzip() {
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-i <vcf>       | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads vcf
	while getopts 'S:s:t:i:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); vcf="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	commander::printinfo "compressing vcf file"

	declare -a cmd1
	helper::makevcfzipcmd -a cmd1 -t $threads -z vcf

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -t $threads -a cmd1
	fi

	return 0
}

variants::vcfnorm() {
	declare -a tdirs
	_cleanup::variants::vcfnorm(){
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
			-r <vcfs>      | array of
			-z <zip>       | true/false compress output
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads dbsnp genome zip=false
	declare -n _vcfs_vcfnorm
	while getopts 'S:s:t:g:d:r:z:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((++mandatory)); _vcfs_vcfnorm="$OPTARG";;
			z) zip="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "normalizing vcf files"

	local instances ithreads
	read -r instances ithreads < <(configure::instances_by_threads -i $((${#_vcfs_vcfnorm[@]}*3)) -t 10 -T $threads)

	declare -a cmd1 cmd2 cmd3 cmd4
	local vcf o e
	for vcf in "${_vcfs_vcfnorm[@]}"; do
		helper::basename -f "$vcf" -o o -e e
		o="$(dirname "$vcf")/$o"

		commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			bcftools view $vcf
		CMD
			vcfix.pl -i - > "$o.fixed.vcf"
		CMD

		commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			bcftools norm -f "$genome" -c s -m-both "$o.fixed.vcf"
		CMD
			vcfix.pl -i - > "$o.fixed.nomulti.vcf"
		CMD

		if [[ $dbsnp ]]; then
			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
				vcffixup "$o.fixed.nomulti.vcf"
			CMD
				vt normalize -q -n -r "$genome" -
			CMD
				vcfixuniq.pl
			CMD
				vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$o.fixed.nomulti.normed.vcf"
			CMD
		else
			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
				vcffixup "$o.fixed.nomulti.vcf"
			CMD
				vt normalize -q -n -r "$genome" -
			CMD
				vcfixuniq.pl > "$o.fixed.nomulti.normed.vcf"
			CMD
		fi

		if $zip; then
			for e in fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
				commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
					bgzip -f -@ $ithreads < "$o.$e" > "$o.$e.gz"
				CMD
					tabix -f -p vcf "$o.$e.gz"
				CMD
					rm "$o.$e"
				CMD
			done
		fi
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
		commander::printcmd -a cmd4
	else
		commander::runcmd -v -b -t $threads -a cmd1
		commander::runcmd -v -b -t $threads -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
		commander::runcmd -v -b -t $instances -a cmd4
	fi

	return 0
}

variants::panelofnormals() {
	declare -a tdirs
	_cleanup::variants::panelofnormals(){
		rm -rf "${tdirs[@]}"
	}

	# The panel of normals not only represents common germline variant sites,
	# it presents commonly noisy sites in sequencing data, e.g. mapping artifacts or
	# other somewhat random but systematic artifacts of sequencing.
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
			-c <sliceinfo> | array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory maxmemory genome tmpdir outdir
	declare -n _mapper_panelofnormals _bamslices_panelofnormals
	while getopts 'S:s:t:g:m:M:r:c:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			m) ((++mandatory)); memory=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			r) ((++mandatory)); _mapper_panelofnormals=$OPTARG;;
			c) ((++mandatory)); _bamslices_panelofnormals=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 7 ]] && _usage

	commander::printinfo "calling panel of normals"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory -M "$maxmemory")

	local m i o t e slice odir instances ithreads odir tdir
	for m in "${_mapper_panelofnormals[@]}"; do
		declare -n _bams_panelofnormals=$m
		((instances+=${#_bams_panelofnormals[@]}))
	done
	read -r i memory < <(configure::memory_by_instances -i $((instances*4)) -T $threads -M "$maxmemory") # for final bcftools sort of vcf fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf
	read -r instances ithreads < <(configure::instances_by_threads -i $instances -T $threads)

	declare -a tomerge cmd1 cmd2 cmd3
	for m in "${_mapper_panelofnormals[@]}"; do
		declare -n _bams_panelofnormals=$m
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		for i in "${!_bams_panelofnormals[@]}"; do
			tomerge=()

			while read -r slice; do
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
						Mutect2
						-I "$slice"
						-O "$slice.vcf"
						-R "$genome"
						-A Coverage
						-A DepthPerAlleleBySample
						-A StrandBiasBySample
						--smith-waterman FASTEST_AVAILABLE
						--f1r2-median-mq 0
						--f1r2-min-bq 20
						--callable-depth 10
						--min-base-quality-score 20
						--minimum-mapping-quality 0
						--native-pair-hmm-threads $mthreads
						-verbosity INFO
						--tmp-dir "${tdirs[-1]}"
						--max-mnp-distance 0
						--ignore-itr-artifacts
				CMD

				tomerge+=("$slice.vcf")
			done < "${_bamslices_panelofnormals[${_bams_panelofnormals[$i]}]}"

			o="$(basename "${_bams_panelofnormals[$i]}")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.bcftools)")
			#DO NOT PIPE - DATALOSS!
			commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bcftools concat -o "$t.vcf" $(printf '"%s" ' "${tomerge[@]}")
			CMD
				bcftools sort -T "${tdirs[-1]}" -m ${memory}M -o "$o.vcf" "$t.vcf"
			CMD

			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bgzip -f -@ $ithreads < "$o.vcf" > "$o.vcf.gz"
			CMD
				tabix -f -p vcf "$o.vcf.gz"
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -c gatk -v -b -t $minstances -a cmd1
		commander::runcmd -v -b -t $threads -a cmd2
		commander::runcmd -v -b -t $instances -a cmd3
	fi

	return 0
}

variants::makepondb() {
	declare -a tdirs
	_cleanup::variants::makepondb(){
		rm -rf "${tdirs[@]}"
		rm -f "$outdir"/*/blocked
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-M <maxmemory> | amount of
			-g <genome>    | path to
			-r <mapper>    | array of sorted, indexed bams within array of
			-p <tmpdir>    | path to
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory genome tmpdir outdir
	declare -n _mapper_makepondb
	while getopts 'S:s:t:M:g:r:p:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			r) ((++mandatory)); _mapper_makepondb=$OPTARG;;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 5 ]] && _usage

	commander::printinfo "creating panel of normals database"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i ${#_mapper_makepondb[@]} -T $threads -M "$maxmemory")

	local m i o t odir params
	declare -a cmd1 cmd2 cmd3
	for m in "${_mapper_makepondb[@]}"; do
		declare -n _bams_makepondb=$m
		odir="$outdir/$m"
		mkdir -p "$odir"

		tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
		for i in "${!_bams_makepondb[@]}"; do
			o="$(basename "${_bams_makepondb[$i]}")"
			echo -e "$m.${o%.*}\t$odir/${o%.*}.vcf.gz" >> ${tdirs[-1]}/import.list
		done

		while [[ -e "$odir/blocked" ]]; do sleep 2; done # do not update unless initialized

		if [[ ! -e "$odir/pondb" && ! -e "$odir/blocked" ]]; then
			touch $odir/blocked # claimed by first parallel instance reaching this line

			bcftools view -h "$odir/${o%.*}.vcf.gz" | head -n -1 > "${tdirs[-1]}/vcf"
			echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tINITIALIZE" >> "${tdirs[-1]}/vcf"
			bgzip -f -@ $threads < "${tdirs[-1]}/vcf" > "${tdirs[-1]}/vcf.gz"
			tabix -f -p vcf "${tdirs[-1]}/vcf.gz"

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
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
					-V "${tdirs[-1]}/vcf.gz"
					--reader-threads $mthreads
					--genomicsdb-workspace-path "$odir/pondb"
					--verbosity INFO
					--tmp-dir "${tdirs[-1]}"
			CMD
				rm -f "$odir/blocked"
			CMD
		fi

		# alternative rm "$odir/pondb" && --genomicsdb-workspace-path "$odir/pondb" or --overwrite-existing-genomicsdb-workspace true
		commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
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
				--sample-name-map "${tdirs[-1]}/import.list"
				--reader-threads $mthreads
				--genomicsdb-update-workspace-path "$odir/pondb"
				--verbosity INFO
				--tmp-dir "${tdirs[-1]}"
		CMD

		# there is no conflict if multiple instances are updating pondb in parallel (non-blocking). however do not write pon.vcf in parallel

		[[ -s "$genome.af_only_gnomad.vcf.gz" ]] && params=" --germline-resource '$genome.af_only_gnomad.vcf.gz'" || params=''
		commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
			while [[ -e "$odir/blocked" ]]; do sleep 2; done
		CMD
			touch "$odir/blocked"
		CMD
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
			rm -f "$odir/blocked"
		CMD
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd2
	else
		commander::runcmd -c gatk -v -b -t $minstances -a cmd1
		commander::runcmd -c gatk -v -b -t $minstances -a cmd2
		commander::runcmd -c gatk -v -b -t $minstances -a cmd3
	fi

	return 0
}

variants::tree(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-m <memory>   | amount of
			-r <mapper>   | array of bams within array of
			-i <vcfdir>   | path to
			-j <caller>   | name of
			-o <outdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory vcfdir caller outdir
    declare -n _mapper_tree _cmpfiles_tree _idfiles_tree
	while getopts 'S:s:t:m:r:i:j:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			r)	((++mandatory)); _mapper_tree=$OPTARG;;
			i)	((++mandatory)); vcfdir="$OPTARG";;
			j)	((++mandatory)); caller="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/genotree"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "reconstructing genotype tree"

	declare -a cmd1 cmd2 cmd3 cmd4 cmd5 cmd6 cmd7 cmd8 cmd9
	declare -a tomerge names
	local m f o e tmp vcffile odir

	for m in "${_mapper_tree[@]}"; do
		declare -n _bams_tree=$m
		odir=$outdir/$m
		mkdir -p $odir
		e=$(echo -e "$(basename ${_bams_tree[0]})\t$(basename ${_bams_tree[1]})" | rev | sed -E 's/(.+).+\t\1.+/\1/' | rev)
		tomerge=()
		names=()
		for f in "${_bams_tree[@]}"; do
			tmp="$vcfdir/$m/$caller/$(basename "$f")"
			vcffile="$(readlink -e "${tmp%.*}".fixed.nomulti.normed.vcf)" || vcffile="$(readlink -e "${tmp%.*}".fixed.nomulti.normed.vcf.gz)"

			names+=("$(basename "$f" $e)")
			o="$odir/${names[-1]}"
			tomerge+=("$o")

			commander::makecmd -a cmd1 -s '|' -o "$o.snv" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
				bcftools view "$vcffile"
			CMD
				perl -F"\t" -M"List::MoreUtils qw(indexes)" -lane '
					if (/^#/){
						($i) = indexes {/^FORMAT$/} @F;
						next;
					}
					@t=split/:/,$F[$i];
					@v=split/:/,$F[-1];
					($cov) = indexes {/^COV$/} @t;
					$cov=$v[$cov];
					($ad) = indexes {/^AD$/} @t;
					$ad=(split/,/,$v[$ad])[1];
					($asf) = indexes {/^ASF$/} @t;
					$asf=$v[$asf];
					print "$F[0]\t".($F[1]-1)."\t$F[1]" if $F[-1]=~/^1\/1/ && length($F[3])==1 && length($F[4])==1 && $cov>=10 && $ad>=3 && $asf >= 0.2;
				'
			CMD

			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
				samtools bedcov -Q 0 "$odir/SNV.bed" "$f" > "$o.snvcov"
			CMD

			commander::makecmd -a cmd6 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- 'CMD' {COMMANDER[4]}<<- CMD
				bcftools view "$vcffile"
			CMD
				awk '/^#/ || (length($4)==1 && length($5)==1 && $4~/[ACGT]/ && $5~/[ACGT]/)'
			CMD
				bedtools intersect -wa -a - -b "$odir/SNVCOV.filtered.bed"
			CMD
				awk '{print $1":"$2":"$4"\t"$5}'
			CMD
				sort -k1,1 > "$o.snvcov.filtered.nt"
			CMD
		done

		commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			sort --parallel=$threads -S ${memory}M -T /dev/shm/ -k1,1 -k2,2n -k3,3n $(printf '"%s" ' "${tomerge[@]/%/.snv}")
		CMD
			uniq > $odir/SNV.bed
		CMD

		commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD
			bedtools unionbedg -i $(printf '"%s" ' "${tomerge[@]/%/.snvcov}") > "$odir/SNVCOV.full.bed"
		CMD

		commander::makecmd -a cmd5 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
			perl -lane '
				print join"\t",@F[0..2] unless grep {$_<10} @F[3..$#F]
			'
		CMD
			"$odir/SNVCOV.full.bed" > "$odir/SNVCOV.filtered.bed"
		CMD

		commander::makecmd -a cmd7 -s '|' -o "$odir/REF.nt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD'
			cut -f 1 $(printf '"%s" ' "${tomerge[@]/%/.snvcov.filtered.nt}")
		CMD
			sort --parallel=$threads -S ${memory}M -T /dev/shm/ -u
		CMD
			sed -E 's/(.+):(\w)$/\1:\2\t\2/'
		CMD

		# replace all 0 by ref nt and remove invariable sites ie. same nt in all samples
		commander::makecmd -a cmd8 -s '|' -o "$odir/ALN.fa" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- 'CMD'
			helper::multijoin -h "$(echo "POS REFERENCE ${names[*]}" | sed 's/ /\t/g')" -f "$odir/REF.nt" $(printf '"%s" ' "${tomerge[@]/%/.snvcov.filtered.nt}")
		CMD
			sed -E -e ':goto; s/^(\S+\s+)(\w)(.*\s+)0/\1\2\3\2/; t goto;' -e '/^\S+(\s+\w)\1*\1$/d'
		CMD
			cut -f 2-
		CMD
			datamash transpose
		CMD
			perl -lane 'print ">$F[0]"; shift @F; print join"",@F'
		CMD

		commander::makecmd -a cmd9 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
			rm -f "$odir/"*.ASCGTRCAT.*
		CMD
			raxmlHPC-PTHREADS-AVX2 -T $threads -f a -N 1000 -x 1234 -p 1234 -s "$odir/ALN.fa" -o REFERENCE -w "$odir" -n ASCGTRCAT.tree -m ASC_GTRCAT -V --asc-corr=lewis
		CMD
			sed -E -e 's/:[0-9\.]+//g' -e 's/\[([0-9]+)\]/\1/g' "$odir/RAxML_bipartitionsBranchLabels.ASCGTRCAT.tree" > "$odir/TREE.newick"
		CMD
			newicktopdf -pc 1 -boot -notitle "$odir/TREE.newick"
		CMD
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
	else
		commander::runcmd -v -b -t $threads -a cmd1
		commander::runcmd -v -b -t 1 -a cmd2
		commander::runcmd -v -b -t $threads -a cmd3
		commander::runcmd -v -b -t $threads -a cmd4
		commander::runcmd -v -b -t $threads -a cmd5
		commander::runcmd -v -b -t $threads -a cmd6
		commander::runcmd -v -b -t $threads -a cmd7
		commander::runcmd -v -b -t $threads -a cmd8
		commander::runcmd -c raxml -v -b -t 1 -a cmd9
	fi

	return 0
}
