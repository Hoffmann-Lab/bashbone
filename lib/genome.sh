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
			-F            | force
		EOF
		return 1
	}

	local OPTIND arg mandatory threads genome tmpdir skip=false skipmd5=false force=false
	while getopts 'S:s:5:t:i:p:F' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			5) $OPTARG && skipmd5=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); genome="$OPTARG";;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			F) force=true;;
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
		read -r instances ithreads jmem jgct jcgct < <(configure::jvm -i 1 -T $threads)
		declare -a cmd1 cmd2

		dict="$(mktemp -u -p "$tmpdir" cleanup.XXXXXXXXXX.dict)"
		commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
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

		commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD
			grep -Eo 'SN:\S+' "$dict" | cut -d ':' -f 2- > "$genome.list"
		CMD
			mv "$dict" "${genome%.*}.dict"
		CMD

		commander::makecmd -a cmd2 -s ';' -c {COMMANDER[0]}<<- CMD
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
			if $force || [[ "$thismd5genome" != "$md5genome" || ! "$thismd5dict" || "$thismd5dict" != "$md5dict" ]]; then
				commander::runcmd -v -b -t $threads -a cmd2
			fi
		fi
	fi

	return 0
}

genome::view(){
	local dict
	_cleanup::genome::view(){
		[[ -e "$pref.bak" ]] && mv "$pref.bak" "$pref"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-m <memory>   | amount of
			-g <genome>   | path to indexed fasta or igv .genome file
			-o <outdir>   | path to batch skript and snapshot files
			-l <files>    | array of gtf/bed/bam/narrowPeak/... file paths to load
			-i <ids>      | array of gene ids to goto and make snapshots
			-j <label>    | array of label for snapshots of gene ids to goto
			-p <pos>      | array of positions to goto and make snapshots (chrom:start-stop)
			-q <pos>      | array of label for snapshots of positions to goto
			-d <number>   | delay seconds between positions
			-r <range>    | of visibility in kb (default: 1000)
			-v <visable>  | pixels per panel (default: 1000)
			-e            | automatically exit after last position
			-s            | do snapshots per position/id
		EOF
		return 1
	}

	local OPTIND arg mandatory genome gtf delay=0 outdir autoexit=false snapshots=false memory range hight=1000
	declare -n _ids_view _pos_view _ids_label_view _pos_label_view _files_view
	while getopts 'm:g:o:l:i:j:p:q:d:r:v:es' arg; do
		case $arg in
			m)	memory=$OPTARG;;
			g)	((++mandatory)); genome="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir";;
			l)	_files_view=$OPTARG;;
			i)	_ids_view=$OPTARG;;
			j)	_ids_label_view=$OPTARG;;
			p)	_pos_view=$OPTARG;;
			q)	_pos_label_view=$OPTARG;;
			d)	delay=$((OPTARG*1000));;
			r)	range=$OPTARG;;
			v)	hight=$OPTARG;;
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

	# sed "s@(DEFAULT_GENOME_KEY=).+@\1$genome@" "$pref" # loads silently without progress bar and does not show igv window unless done
	# deletion defaults to hg19 and downloads it unless available
	# sed -i -E -e '/DEFAULT_GENOME_KEY/d' "$pref"
	# see also https://github.com/igvteam/igv/blob/master/src/main/resources/org/broad/igv/prefs/preferences.tab

	[[ -e "$pref" ]] && mv "$pref" "$pref.bak"
	cat <<- EOF > "$pref"
		SAM.SHOW_MISMATCHES=true
		SAM.MAX_VISIBLE_RANGE=${range:-1000}
		SAM.SHOW_SOFT_CLIPPED=false
		SAM.SHOW_ALL_BASES=false
		SAM.SHOW_CENTER_LINE=true
		SAM.DOWNSAMPLE_READS=false
		SAM.FILTER_DUPLICATES=false
		SAM.FILTER_FAILED_READS=false
		SAM.ALIGNMENT_SCORE_THRESHOLD=0
		SAM.QUALITY_THRESHOLD=0
		SAM.SORT_OPTION=START
		SAM.COLOR_BY=FIRST_OF_PAIR_STRAND
		DETAILS_BEHAVIOR=CLICK

		##RNA
		SAM.SHOW_JUNCTION_TRACK=true
		SAM.SHOW_COV_TRACK=true
		SAM.SHOW_ALIGNMENT_TRACK=true
		SAM.MAX_VISIBLE_RANGE=${range:-1000}

		##THIRD_GEN
		SAM.DOWNSAMPLE_READS=false
		SAM.MAX_VISIBLE_RANGE=${range:-1000}
	EOF

	local i x l f="$outdir/run.batch"

	# the genome path will always be interpreted as an id and thus must not be quoted!
	cat <<- EOF > "$f"
		new
		snapshotDirectory "$outdir"
		setSleepInterval 0
		maxPanelHeight $hight

		genome $genome
	EOF
	for i in "${_files_view[@]}"; do
		echo "load \"$i\"" >> "$f"
	done
	echo "sort FIRSTOFPAIRSTRAND" >> "$f"

	for x in "${!_ids_view[@]}"; do
		i="${_ids_view[$x]}"
		l="${_ids_label_view[$x]}"
		echo "goto $(grep -E -m 1 $'\t'gene$'\t'".+gene_id \"$i\"" $gtf | awk '{print $1":"$4"-"$5}')" >> $f
		if [[ $delay -gt 0 ]]; then
			cat <<- EOF >> $f
				setSleepInterval $delay
				echo
				setSleepInterval 0
			EOF
		fi
		if $snapshots; then
			cat <<- EOF >> $f
				snapshot ${l:+$l.}$i.png
				snapshot ${l:+$l.}.$i.svg
			EOF
		fi
	done
	for x in "${!_pos_view[@]}"; do
		i="${_pos_view[$x]}"
		l="${_pos_label_view[$x]}"
		echo "goto $i" >> $f
		if [[ $delay -gt 0 ]]; then
			cat <<- EOF >> $f
				setSleepInterval $delay
				echo
				setSleepInterval 0
			EOF
		fi
		if $snapshots; then
			cat <<- EOF >> $f
				collapse
				snapshot ${l:+$l.}$i.png
				snapshot ${l:+$l.}$i.svg
			EOF
		fi
	done
	$autoexit && echo "exit" >> $f

	declare -a cmd1
	commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
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
