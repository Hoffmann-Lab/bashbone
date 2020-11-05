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
			-g <genome>     | path to
			-f <annotation> | path to gtf or igv .genome file
			-o <outdir>     | path to
			-l <files>      | array of bed/bam/narrowPeak/... file paths to load
			-i <ids>        | array of gene ids to goto and make snapshots
			-p <pos>        | array of positions to goto and make snapshots (chrom:start-stop)
			-d <number>     | delay milliseconds between positions
			-e <exit>       | automatically after last position
			-s <snapshot>   | true/false per position/id
		EOF
		return 1
	}

	local OPTIND arg mandatory genome gtf delay=0 outdir autoexit=false snapshots=false
	declare -n _ids_view _pos_view _files_view
	while getopts 'g:f:o:l:i:p:d:es' arg; do
		case $arg in
			g)	((++mandatory)); genome="$OPTARG";;
			f)	gtf="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG/igv"; mkdir -p "$outdir";;
			l)	_files_view=$OPTARG;;
			i)	_ids_view=$OPTARG;;
			p)	_pos_view=$OPTARG;;
			d)	delay="$OPTARG";;
			s)	snapshots=true;;
			e)	autoexit=true;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage
	[[ ! $gtf && ${#_ids_view[@]} -gt 0 ]] && _usage

	f=$outdir/run.batch

	cat <<- EOF > $f
		new
		snapshotDirectory $outdir
		setSleepInterval 0

		genome $genome
		load $gtf
		load $(printf '"%s",' "${_files_view[@]}")
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

	declare -a cmd1=('igv --batch "$f"')
	commander::runcmd -c igv -v -a cmd1
}
