#! /usr/bin/env bash
# (c) Konstantin Riege

function variants::vcfnorm(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-d <dbsnp>     | path to
			-r <vcfs>      | array of
			-z             | compress output
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads dbsnp genome zip=false tmpdir="${TMPDIR:-/tmp}" maxmemory
	declare -n _vcfs_vcfnorm
	while getopts 'S:s:t:g:d:r:M:z' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			M) maxmemory=$OPTARG;;
			t) ((++mandatory)); threads=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			d) dbsnp="$OPTARG";;
			r) ((++mandatory)); _vcfs_vcfnorm="$OPTARG";;
			z) zip=true;;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "normalizing vcf files"

	local instances ithreads ignore n=3
	[[ $dbsnp ]] && n=4
	read -r instances ithreads imemory ignore < <(configure::jvm -i $((${#_vcfs_vcfnorm[@]}*n)) -t 10 -T $threads -M "$maxmemory")

	declare -a tdirs cmd1 cmd2 cmd3 cmd4
	local f o e
	for f in "${_vcfs_vcfnorm[@]}"; do
		helper::basename -f "$f" -o o -e e
		o="$(dirname "$f")/$o"

		commander::makecmd -a cmd1 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			bcftools view "$f"
		CMD
			vcfix.pl -i - > "$o.fixed.vcf"
		CMD

		commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			bcftools norm -f "$genome" -c s -m-both "$o.fixed.vcf"
		CMD
			vcfix.pl -i - > "$o.fixed.nomulti.vcf"
		CMD

		if [[ $dbsnp ]]; then
			commander::makecmd -a cmd3 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- CMD
				vcffixup "$o.fixed.nomulti.vcf"
			CMD
				vt normalize -q -n -r "$genome" -
			CMD
				vcfixuniq.pl
			CMD
				tee -i "$o.fixed.nomulti.normed.vcf"
			CMD
				vcftools --vcf - --exclude-positions "$dbsnp" --recode --recode-INFO-all --stdout > "$o.fixed.nomulti.normed.nodbsnp.vcf"
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
			for e in fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf $([[ $dbsnp ]] && echo fixed.nomulti.normed.nodbsnp.vcf); do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.vcfnorm)")
				commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools sort -T "${tdirs[-1]}" -m ${imemory}M "$o.$e" | bgzip -k -c -@ $ithreads | tee "$o.$e.gz" | bcftools index --threads $ithreads -f -t -o "$o.$e.gz.tbi"
				CMD
					rm "$o.$e"
				CMD
			done
		else
			for e in fixed.vcf fixed.nomulti.vcf fixed.nomulti.normed.vcf $([[ $dbsnp ]] && echo fixed.nomulti.normed.nodbsnp.vcf); do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.vcfnorm)")
				commander::makecmd -a cmd4 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
					bcftools sort -T "${tdirs[-1]}" -m ${imemory}M "$o.$e" > "${tdirs[-1]}/$(basename "$o.$e")"
				CMD
					mv "${tdirs[-1]}/$(basename "$o.$e")" "$o.$e"
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
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $instances -a cmd4
	fi

	return 0
}

function variants::panelofnormals(){
	# The panel of normals not only represents common germline variant sites,
	# it presents commonly noisy sites in sequencing data, e.g. mapping artifacts or
	# other somewhat random but systematic artifacts of sequencing.
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-g <genome>    | path to
			-m <memory>    | amount of
			-M <maxmemory> | amount of
			-r <mapper>    | array of sorted, indexed bams within array of
			-c <sliceinfo> | array of
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads memory maxmemory genome tmpdir="${TMPDIR:-/tmp}" outdir
	declare -n _mapper_panelofnormals _bamslices_panelofnormals
	while getopts 'S:s:t:g:m:M:r:c:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			m) ((++mandatory)); memory=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			r) ((++mandatory)); _mapper_panelofnormals=$OPTARG;;
			c) ((++mandatory)); _bamslices_panelofnormals=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "calling panel of normals"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -T $threads -m $memory -M "$maxmemory")
	declare -n _bams_panelofnormals="${_mapper_panelofnormals[0]}"
	local ignore instances1 ithreads1 imemory1 instances2 ithreads2 imemory2
	local instances=$((${#_mapper_panelofnormals[@]}*${#_bams_panelofnormals[@]}))
	read -r instances1 ithreads1 imemory1 ignore < <(configure::jvm -i $((instances*minstances)) -t 10 -T $threads) # for slice wise bcftools sort + bgzip
	read -r instances2 ithreads2 imemory2 ignore < <(configure::jvm -i $instances -t 10 -T $threads) # for final bcftools concat of all 4 vcf files

	local m i o t slice odir tdir
	declare -a tdirs tomerge cmd1 cmd2 cmd3
	for m in "${_mapper_panelofnormals[@]}"; do
		declare -n _bams_panelofnormals=$m
		odir="$outdir/$m"
		tdir="$tmpdir/$m"
		mkdir -p "$odir" "$tdir"

		for i in "${!_bams_panelofnormals[@]}"; do
			o="$(basename "${_bams_panelofnormals[$i]}")"
			t="$tdir/${o%.*}"
			o="$odir/${o%.*}"
			tomerge=()

			while read -r slice; do
				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
				commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
					MALLOC_ARENA_MAX=4 gatk
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

				tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
				commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD
					bcftools sort -T "${tdirs[-1]}" -m ${imemory1}M "$slice.$e"
				CMD
					bgzip -k -c -@ $ithreads1
				CMD
					tee "$slice.$e.gz"
				CMD
					bcftools index --threads $ithreads1 -f -t -o "$slice.$e.gz.tbi"
				CMD

				tomerge+=("$slice.vcf.gz")
			done < "${_bamslices_panelofnormals[${_bams_panelofnormals[$i]}]}"

			tdirs+=("$(mktemp -d -p "$tmpdir" cleanup.XXXXXXXXXX.gatk)")
			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				bcftools concat --threads $ithreads2 -O z -a -d all -o "$t.vcf.gz" $(printf '"%s" ' "${tomerge[@]}")
			CMD
				bcftools sort -T "${tdirs[-1]}" -m ${imemory2}M "$t.vcf.gz" | bgzip -k -c -@ $ithreads2 | tee "$o.vcf.gz" | bcftools index --threads $ithreads2 -f -t -o "$o.vcf.gz.tbi"
			CMD
		done
	done

	if $skip; then
		commander::printcmd -a cmd1
		commander::printcmd -a cmd2
		commander::printcmd -a cmd3
	else
		commander::runcmd -c gatk -v -b -i $minstances -a cmd1
		commander::runcmd -v -b -i $threads -a cmd2
		commander::runcmd -v -b -i $instances -a cmd3
	fi

	return 0
}

function variants::makepondb(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
			-S <hardskip>  | true/false return
			-s <softskip>  | true/false only print commands
			-t <threads>   | number of
			-M <maxmemory> | amount of
			-g <genome>    | path to
			-r <mapper>    | array of sorted, indexed bams within array of
			-o <outdir>    | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory skip=false threads maxmemory genome tmpdir="${TMPDIR:-/tmp}" outdir
	declare -n _mapper_makepondb
	while getopts 'S:s:t:M:g:r:o:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			M) maxmemory=$OPTARG;;
			g) ((++mandatory)); genome="$OPTARG";;
			r) ((++mandatory)); _mapper_makepondb=$OPTARG;;
			o) ((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*) _usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 4 ]] && _usage

	commander::printinfo "creating panel of normals database"

	local minstances mthreads jmem jgct jcgct
	read -r minstances mthreads jmem jgct jcgct < <(configure::jvm -i ${#_mapper_makepondb[@]} -T $threads -M "$maxmemory")

	local m i o t odir params
	declare -a tdirs cmd1 cmd2 cmd3
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
			touch "$odir/blocked" # claimed by first parallel instance reaching this line
			echo "rm -f '$odir/blocked'" >> "$BASHBONE_CLEANUP"

			bcftools view -h "$odir/${o%.*}.vcf.gz" | head -n -1 > "${tdirs[-1]}/vcf"
			echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tINITIALIZE" >> "${tdirs[-1]}/vcf"
			bgzip -k -c -@ $threads "${tdirs[-1]}/vcf" | tee "${tdirs[-1]}/vcf.gz" | bcftools index -f -t -o "${tdirs[-1]}/vcf.gz.tbi"

			commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
				MALLOC_ARENA_MAX=4 gatk
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
			MALLOC_ARENA_MAX=4 gatk
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
			MALLOC_ARENA_MAX=4 gatk
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
		commander::runcmd -c gatk -v -b -i $minstances -a cmd1
		commander::runcmd -c gatk -v -b -i $minstances -a cmd2
		commander::runcmd -c gatk -v -b -i $minstances -a cmd3
	fi

	return 0
}

function variants::tree(){
	function _usage(){
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[-2]} usage:
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

	local OPTIND arg mandatory skip=false threads memory vcfdir caller outdir tmpdir="${TMPDIR:-/tmp}"
    declare -n _mapper_tree
	while getopts 'S:s:t:m:r:i:j:o:' arg; do
		case $arg in
			S)	$OPTARG && return 0;;
			s)	$OPTARG && skip=true;;
			t)	((++mandatory)); threads=$OPTARG;;
			m)	((++mandatory)); memory=$OPTARG;;
			r)	((++mandatory)); _mapper_tree=$OPTARG;;
			i)	((++mandatory)); vcfdir="$OPTARG";;
			j)	((++mandatory)); caller="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			*)	_usage;;
		esac
	done
	[[ $# -eq 0 ]] && { _usage || return 0; }
	[[ $mandatory -lt 6 ]] && _usage

	commander::printinfo "reconstructing genotype tree"

	declare -a cmd1 cmd2 cmd3 cmd4 cmd5 cmd6 cmd7 cmd8 cmd9
	declare -a tomerge names
	local m f o e vcffile odir

	for m in "${_mapper_tree[@]}"; do
		declare -n _bams_tree=$m
		odir="$outdir/$m/$caller"
		mkdir -p "$odir"
		e="$(echo -e "$(basename "${_bams_tree[0]}")\t$(basename "${_bams_tree[1]}")" | rev | sed -E 's/(.+).+\t\1.+/\1/' | rev)"
		tomerge=()
		names=()
		for f in "${_bams_tree[@]}"; do
			vcffile=$(find -L "$vcfdir/$m/$caller" -name "$(basename "$f" "$e")*.vcf" -or -name "$(basename "$f" "$e")*.vcf.gz" | grep . | perl -e 'print "".(sort { length($b) <=> length($a) } <>)[0]')
			names+=("$(basename "$f" "$e")")
			o="$odir/${names[-1]}"
			tomerge+=("$o")

			# commander::makecmd -a cmd1 -s '|' -o "$o.snv" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
			# 	bcftools view "$vcffile"
			# CMD
			# 	perl -F"\t" -M"List::MoreUtils qw(indexes)" -lane '
			# 		if (/^#/){
			# 			($i) = indexes {/^FORMAT$/} @F;
			# 			next;
			# 		}
			# 		@t=split/:/,$F[$i];
			# 		@v=split/:/,$F[-1];
			# 		($cov) = indexes {/^COV$/} @t;
			# 		$cov=$v[$cov];
			# 		($ad) = indexes {/^AD$/} @t;
			# 		$ad=(split/,/,$v[$ad])[1];
			# 		($asf) = indexes {/^ASF$/} @t;
			# 		$asf=$v[$asf];
			# 		print "$F[0]\t".($F[1]-1)."\t$F[1]" if $F[-1]=~/^1(\/|\|)1/ && length($F[3])==1 && length($F[4])==1 && $cov>=10 && $ad>=3 && $asf >= 0.2;
			# 	'
			# CMD
			commander::makecmd -a cmd1 -s '|' -o "$o.snv" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD'
				bcftools view -U -v snps -g hom -f PASS,vcfix_ok,. -H "$vcffile"
			CMD
				awk -v OFS='\t' '{print $1,($2-1),$2}'
			CMD

			commander::makecmd -a cmd3 -s ';' -c {COMMANDER[0]}<<- CMD
				samtools bedcov -Q 0 -g 1796 -j "$odir/SNV.bed" "$f" > "$o.snvcov"
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

		# sort --parallel=$threads -S ${memory}M -T "$tmpdir" -k1,1 -k2,2n -k3,3n $(printf '"%s" ' "${tomerge[@]/%/.snv}")
		commander::makecmd -a cmd2 -s '|' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
			cat $(printf '"%s" ' "${tomerge[@]/%/.snv}")
		CMD
			helper::sort -t $threads -M $memory -k1,1 -k2,2n -k3,3n
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

		# sort --parallel=$threads -S ${memory}M -T "$tmpdir" -u
		commander::makecmd -a cmd7 -s '|' -o "$odir/REF.nt" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- 'CMD'
			cut -f 1 $(printf '"%s" ' "${tomerge[@]/%/.snvcov.filtered.nt}")
		CMD
			helper::sort -t $threads -M $memory -u
		CMD
			sed -E 's/(.+):(\w)$/\1:\2\t\2/'
		CMD

		# replace all 0 by ref nt and remove invariable sites ie. same nt in all samples
		commander::makecmd -a cmd8 -s '|' -o "$odir/ALN.fa" -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- 'CMD' {COMMANDER[2]}<<- CMD {COMMANDER[3]}<<- CMD {COMMANDER[4]}<<- 'CMD'
			helper::multijoin -t $((threads/${#_mapper_tree[@]}+1)) -e 0 -h "$(echo "POS REFERENCE ${names[*]}" | sed 's/ /\t/g')" "$odir/REF.nt" $(printf '"%s" ' "${tomerge[@]/%/.snvcov.filtered.nt}")
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
		commander::runcmd -v -b -i $threads -a cmd1
		commander::runcmd -v -b -i 1 -a cmd2
		commander::runcmd -v -b -i $threads -a cmd3
		commander::runcmd -v -b -i $threads -a cmd4
		commander::runcmd -v -b -i $threads -a cmd5
		commander::runcmd -v -b -i $threads -a cmd6
		commander::runcmd -v -b -i $threads -a cmd7
		commander::runcmd -v -b -i ${#_mapper_tree[@]} -a cmd8
		commander::runcmd -c raxml -v -b -i 1 -a cmd9
	fi

	return 0
}
