#! /usr/bin/env bash
# (c) Konstantin Riege

genome::mkdict() {
	local dict
	_cleanup::genome::mkdict(){
		rm -f "$dict"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-5 <skip>     | true/false md5sums, indexing respectively
			-t <threads>  | number of
			-i <genome>   | path to
			-p <tmpdir>   | path to
		EOF
		return 1
	}

	local OPTIND arg mandatory threads genome tmpdir skip=false skipmd5=false
	while getopts 'S:s:5:t:i:p:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			5) $OPTARG && skipmd5=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); genome="$OPTARG";;
			p) ((++mandatory)); tmpdir="$OPTARG";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "creating genome dictionary"

	if $skipmd5; then
		commander::warn "skip checking md5 sums and genome dictionary creation respectively"
	else
		commander::printinfo "checking md5 sums"

		local instances ithreads jmem jgct jcgct
		read -r instances ithreads jmem jgct jcgct < <(configure::jvm -T $threads)
		declare -a cmd1 cmd2

		dict="$(mktemp -u -p "$tmpdir" cleanup.XXXXXXXXXX.dict)"
		commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
			picard
				-Xmx${jmem}m
				-XX:ParallelGCThreads=$jgct
				-XX:ConcGCThreads=$jcgct
				-Djava.io.tmpdir="$tmpdir"
				CreateSequenceDictionary
				R="$genome"
				O="$dict"
				VERBOSITY=WARNING
		CMD

		commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			grep -Eo 'SN:\S+' "$dict" | cut -d ':' -f 2- > "$genome.list"
		CMD
			mv "$dict" "${genome%.*}.dict"
		CMD

		commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
			samtools faidx "$genome"
		CMD

		if $skip; then
			commander::printcmd -a cmd1
			commander::printcmd -a cmd2
		else
			commander::runcmd -c picard -v -b -t $threads -a cmd1
			local md5dict thismd5genome thismd5dict
			md5dict=$(md5sum "$dict" | cut -d ' ' -f 1)
			thismd5genome=$(md5sum "$genome" | cut -d ' ' -f 1)
			[[ -s "${genome%.*}.dict" ]] && thismd5dict=$(md5sum "${genome%.*}.dict" | cut -d ' ' -f 1)
			if [[ "$thismd5genome" != "$md5genome" || ! "$thismd5dict" || "$thismd5dict" != "$md5dict" ]]; then
				commander::runcmd -v -b -t $threads -a cmd2
			fi
		fi
	fi

	return 0
}

genome::view(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-m <memory>   | amount of
			-g <genome>   | path to indexed fasta or igv .genome file
			-o <outdir>   | path to
			-l <files>    | array of gtf/bed/bam/narrowPeak/... file paths to load
			-i <ids>      | array of gene ids to goto and make snapshots
			-p <pos>      | array of positions to goto and make snapshots (chrom:start-stop)
			-d <number>   | delay seconds between positions
			-e <exit>     | true/false automatically after last position
			-s <snapshot> | true/false per position/id
		EOF
		return 1
	}

	local OPTIND arg mandatory genome gtf delay=0 outdir autoexit=false snapshots=false memory
	declare -n _ids_view _pos_view _files_view
	while getopts 'm:g:o:l:i:p:d:es' arg; do
		case $arg in
			m)	memory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			l)	_files_view=$OPTARG;;
			i)	_ids_view=$OPTARG;;
			p)	_pos_view=$OPTARG;;
			d)	delay=$((OPTARG*1000));;
			s)	snapshots=true;;
			e)	autoexit=true;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage
	[[ ! $gtf && ${#_ids_view[@]} -gt 0 ]] && _usage

	[[ $memory ]] || {
		local instances memory
		read -r instances memory < <(configure::memory_by_instances -i 1 -T 1)
	}

	# hack to prevent loading last genome before loading user defined stuff - problem:
	declare -a cmd=("unset DISPLAY; igv")
	local pref="$(commander::runcmd -c igv -a cmd |& grep 'IGV Directory:' | awk -F': ' '{print $NF}')"/prefs.properties
	[[ ! -s "$pref" ]] && echo -e "##RNA\n##THIRD_GEN" > "$pref"
	#[[ -e $pref ]] && sed "s@(DEFAULT_GENOME_KEY=).+@\1$genome@" "$pref" # loads silently without progress bar and does not show igv window unless done
	# deletion defaults to hg19 and downloads it unless available
	# sed -i -E -e '/DEFAULT_GENOME_KEY/d' "$pref"

	# see https://github.com/igvteam/igv/blob/master/src/main/resources/org/broad/igv/prefs/preferences.tab
	# slurp (:a;ba) section from ##RNA to e.g. ##THIRD_GEN or EOF into holdspace (H;n), then switch to holdspace (x) and do replacements/appends (s).
	# finally print (p), switch back to pattern space (x) and handle EOF
	sed -i -nE -e '/DEFAULT_GENOME_KEY/d;' \
	-e '1{s/(DETAILS_BEHAVIOR=.*|^)/DETAILS_BEHAVIOR=CLICK\n/p};' -e '/DETAILS_BEHAVIOR/d;' \
	-e '/##RNA/!p; /##RNA/{:slurp; $!{H;n;/##/!{b slurp}}; ${/##/!{H}};x; s/^\n//;' \
	-e 's/(\s*)(SAM.SHOW_MISMATCHES)=[^\n]*/\1\2=true/g; /SAM.SHOW_MISMATCHES=/!{s/\s*$/\nSAM.SHOW_MISMATCHES=true/};' \
	-e 's/(\s*)(DETAILS_BEHAVIOR)=[^\n]*/\1\2=CLICK/g; /DETAILS_BEHAVIOR=/!{s/\s*$/\nDETAILS_BEHAVIOR=CLICK/};' \
	-e 's/(\s*)(SAM.MAX_VISIBLE_RANGE)=[^\n]*/\1\2=1000/g; /SAM.MAX_VISIBLE_RANGE=/!{s/\s*$/\nSAM.MAX_VISIBLE_RANGE=1000/};' \
	-e 's/(\s*)(SAM.BASE_QUALITY_MIN)=[^\n]*/\1\2=5/g; /SAM.BASE_QUALITY_MIN=/!{s/\s*$/\nSAM.BASE_QUALITY_MIN=5/};' \
	-e 's/(\s*)(SAM.SORT_OPTION)=[^\n]*/\1\2=FIRST_OF_PAIR_STRAND/g; /SAM.SORT_OPTION=/!{s/\s*$/\nSAM.SORT_OPTION=FIRST_OF_PAIR_STRAND/};' \
	-e 's/(\s*)(SAM.DOWNSAMPLE_READS)=[^\n]*/\1\2=true/g; /SAM.DOWNSAMPLE_READS=/!{s/\s*$/\nSAM.DOWNSAMPLE_READS=true/};' \
	-e 's/(\s*)(SAM.COLOR_BY)=[^\n]*/\1\2=FIRST_OF_PAIR_STRAND/g; /SAM.COLOR_BY=/!{s/\s*$/\nSAM.COLOR_BY=FIRST_OF_PAIR_STRAND/};' \
	-e 'p;x;$!{p};${/##/p}}' \
	"$pref"

	f="$outdir/run.batch"

	# the genome path will always be interpreted as an id and thus must not be quoted!
	cat <<- EOF > "$f"
		new
		snapshotDirectory "$outdir"
		setSleepInterval 0

		genome $genome
		load $(printf '"%s",' "${_files_view[@]}" | sed 's/,$//')
		sort FIRSTOFPAIRSTRAND
	EOF

	local i
	for i in "${_ids_view[@]}"; do
		echo "goto $(grep -E -m 1 $'\t'gene$'\t'".+gene_id \"$i\"" $gtf | awk '{print $1":"$4"-"$5}')" >> $f
		if [[ $delay -gt 0 ]]; then
			cat <<- EOF >> $f
				setSleepInterval $delay
				echo
				setSleepInterval 0
			EOF
		fi
		$snapshots && echo "snapshot $i.jpg" >> $f
	done
	for i in "${_pos_view[@]}"; do
		echo "goto $i" >> $f
		if [[ $delay -gt 0 ]]; then
			cat <<- EOF >> $f
				setSleepInterval $delay
				echo
				setSleepInterval 0
			EOF
		fi
		$snapshots && echo "snapshot $i.jpg" >> $f
	done
	$autoexit && echo "exit" >> $f

	declare -a cmd1
	commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
		java
			--module-path="\$CONDA_PREFIX/lib/igv"
			-Xmx${memory}m
			@"\$CONDA_PREFIX/lib/igv/igv.args"
			-Dapple.laf.useScreenMenuBar=true
			-Djava.net.preferIPv4Stack=true
			--module=org.igv/org.broad.igv.ui.Main
			--batch "$f"
	CMD

	commander::runcmd -c igv -v -a cmd1
}
