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

genome::indexgtf(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-5 <skip>     | true/false md5sums, indexing respectively
			-t <threads>  | number of
			-i <gtf>      | path to
			-F            | force
		EOF
		return 1
	}

	local OPTIND arg mandatory threads gtf skip=false skipmd5=false force=false
	while getopts 'S:s:5:t:i:F' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			5) $OPTARG && skipmd5=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			i) ((++mandatory)); gtf="$OPTARG";;
			F) force=true;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage

	commander::printinfo "indexing gtf"

	if $skipmd5; then
		commander::warn "skip checking md5 sums and gtf index creation respectively"
	else
		commander::printinfo "checking md5 sums"

		declare -a cmd1

		dict="$(mktemp -u -p "$tmpdir" cleanup.XXXXXXXXXX.dict)"
		commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD {COMMANDER[1]}<<- CMD {COMMANDER[2]}<<- CMD
			bgzip -k -c -@ $threads "$gtf" > "$gtf.gz"
		CMD
			tabix -f "$gtf.gz"
		CMD
			mv "$gtf.gz.tbi" "$gtf.gz.igv"
		CMD

		if $skip; then
			commander::printcmd -a cmd1
		else
			local thismd5gtf=$(md5sum "$gtf" | cut -d ' ' -f 1)
			if $force || [[ "$thismd5gtf" && "$thismd5gtf" != "$md5gtf" ]]; then
				commander::runcmd -v -b -t 1 -a cmd1
			fi
		fi
	fi

	return 0
}

genome::mkgodb(){
	local tmpdir
	_cleanup::genome::mkdict(){
		rm -f "$tmpdir/org.My.eg.db"
	}

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-S <hardskip> | true/false return
			-s <softskip> | true/false only print commands
			-t <threads>  | number of
			-p <tmpdir>   | path to
			-g <gofile>   | path to 4-column tab separated file with BP, MF and CC or see column 3 below
			                ..
			                ENSG00000199065 GO:0005615 cellular_component extracellular space
			                ENSG00000199065 GO:1903231 molecular_function mRNA binding involved in posttranscriptional gene silencing
			                ENSG00000199065 GO:0035195 biological_process gene silencing by miRNA
			                ..
		EOF
		return 1
	}

	local OPTIND arg mandatory threads gofile skip=false
	while getopts 'S:s:t:p:g:' arg; do
		case $arg in
			S) $OPTARG && return 0;;
			s) $OPTARG && skip=true;;
			t) ((++mandatory)); threads=$OPTARG;;
			g) ((++mandatory)); gofile="$OPTARG";;
			p) ((++mandatory)); tmpdir="$OPTARG"; mkdir -p "$tmpdir";;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 3 ]] && _usage

	commander::printinfo "creating genome go orgdb"

	declare -a cmd1
	commander::makecmd -a cmd1 -s ' ' -c {COMMANDER[0]}<<- 'CMD' {COMMANDER[1]}<<- CMD
		Rscript - <<< '
			args <- commandArgs(TRUE);
			threads <- args[1];
			tmpdir <- args[2];
			libdir <- args[3];
			gofile <- args[4];

			unlink(libdir, recursive = T);
			dir.create(libdir, recursive = T, mode="0750");
			BiocManager::install(c("GO.db","AnnotationForge"), dependencies=T, ask=F, force=T, Ncpus=threads, clean=T, destdir=tmpdir, lib=libdir);

			.libPaths(c(libdir,.libPaths()));
			suppressMessages(library(GO.db));
			suppressMessages(library(AnnotationForge));

			df <- read.table(gofile, sep="\t", stringsAsFactors=F, check.names=F, quote="", header=F, na.strings="", col.names=c("GID","GO","ONTOLOGY","DESCRIPTION"));
			df <- unique(df[!is.na(df$GO) & !is.na(df$ONTOLOGY),c(1,2)]);
			df$EVIDENCE <- "IEA";

			makeOrgPackage(
			  go=df,
			  version="0.1",
			  maintainer="Anony Mous <anon@anony.mous>",
			  author="Anony Mous <anon@anony.mous>",
			  outputDir = tmpdir,
			  tax_id="0",
			  genus="M",
			  species="y",
			  goTable="go"
			);

			unlink(file.path(libdir,"org.My.eg.db"), recursive = T);
			install.packages(file.path(tmpdir,"org.My.eg.db"), repos=NULL, Ncpus=1, clean=T, destdir=tmpdir, lib=libdir);
			unlink(file.path(tmpdir,"org.My.eg.db"), recursive = T);
		'
	CMD
		$threads "$tmpdir" "$(dirname "$gofile")/oRgdb" "$gofile"
	CMD

	if $skip; then
		commander::printcmd -a cmd1
	else
		commander::runcmd -v -b -t 1 -a cmd1
	fi

	return 0
}


genome::view(){
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
			${FUNCNAME[1]} usage:
			-m <memory>   | amount of
			-i <genome>   | path to indexed fasta or igv .genome file
			-g <gtf>      | path to fasta matching, indexed gtf file
			-o <outdir>   | path to batch skript, current igv config and snapshot files
			-f <files>    | array of gtf/bed/bam/narrowPeak/... file paths to load
			-x <ids>      | array of gene ids to goto and make snapshots
			-y <label>    | array of label for snapshots of gene ids to goto
			-p <pos>      | array of positions to goto and optional make snapshots (chrom:start-stop)
			-q <label>    | array of label for snapshots of positions to goto
			-d <number>   | delay seconds between positions
			-r <range>    | of visibility in kb (default: 1000)
			-v <visable>  | pixels per panel (default: 1000)
			-e            | automatically exit after last position
			-s            | do snapshots per position/id
			-n            | enable searching for gene names (during startup, igv seems to be frozen)
		EOF
		return 1
	}

	local OPTIND arg mandatory genome gtf delay=0 outdir autoexit=false snapshots=false memory range hight=1000 searchable=false
	declare -n _ids_view _pos_view _ids_label_view _pos_label_view _files_view
	while getopts 'm:i:g:o:f:x:y:p:q:d:r:v:esn' arg; do
		case $arg in
			m)	memory=$OPTARG;;
			i)	((++mandatory)); genome="$OPTARG";;
			g)	gtf="$OPTARG";;
			o)	((++mandatory)); outdir="$OPTARG"; mkdir -p "$outdir/snapshots";;
			f)	_files_view=$OPTARG;;
			x)	_ids_view=$OPTARG;;
			y)	_ids_label_view=$OPTARG;;
			p)	_pos_view=$OPTARG;;
			q)	_pos_label_view=$OPTARG;;
			d)	delay=$((OPTARG*1000));;
			r)	range=$OPTARG;;
			v)	hight=$OPTARG;;
			s)	snapshots=true;;
			e)	autoexit=true;;
			n)	searchable=true;;
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 2 ]] && _usage
	[[ ! $gtf && ${#_ids_view[@]} -gt 0 ]] && _usage

	[[ $memory ]] || {
		local instances memory
		read -r instances memory < <(configure::memory_by_instances -i 1 -T 1)
	}

	rm -rf "$outdir/igv"
	mkdir -p "$outdir/igv/genomes"
	echo -e "$(basename "$genome")\t$outdir/igv/genomes/current.json\tcurrent" > "$outdir/igv/genomes/user-defined-genomes.txt"

	cat <<- EOF > "$outdir/igv/prefs.properties"
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
		DEFAULT_GENOME_KEY=current

		##RNA
		SAM.SHOW_JUNCTION_TRACK=true
		SAM.SHOW_COV_TRACK=true
		SAM.SHOW_ALIGNMENT_TRACK=true
		SAM.MAX_VISIBLE_RANGE=${range:-1000}

		##THIRD_GEN
		SAM.DOWNSAMPLE_READS=false
		SAM.MAX_VISIBLE_RANGE=${range:-1000}
	EOF

	local i x l f="$outdir/igv/igv.batch"

	cat <<- EOF > "$outdir/igv/genomes/current.json"
		{
		  id: "current",
		  name: "$(basename "$genome")",
		  fastaURL: "$(realpath -s "$genome")",
		  indexURL: "$(realpath -s "$genome.fai")"$([[ $gtf ]] && echo ',')
	EOF
	if [[ $gtf ]]; then
		cat <<- EOF >> "$outdir/igv/genomes/current.json"
			  tracks: [
			    {
			      name: "$(basename "$gtf")",
			      url: "$(realpath -s "$gtf.gz")",
		EOF
		$searchable || echo "      indexURL: \"$(realpath -s "$gtf.gz.igv")\"," >> "$outdir/igv/genomes/current.json"
		cat <<- EOF >> "$outdir/igv/genomes/current.json"
			      type: "annotation",
			      format: "gtf",
			      order: Number.MAX_VALUE
			    }
			  ]
		EOF
	fi
	echo '}' >> "$outdir/igv/genomes/current.json"

	# do not use "new"
	# no requirement for "genome current", because this is default
	cat <<- EOF > "$f"
		snapshotDirectory "$outdir/snapshots"
		setSleepInterval 0
		maxPanelHeight $hight
	EOF
	for i in "${_files_view[@]}"; do
		echo "load \"$(realpath -s "$i")\"" >> "$f"
	done
	echo "sort FIRSTOFPAIRSTRAND" >> "$f"

	for x in "${!_ids_view[@]}"; do
		i="${_ids_view[$x]}"
		l="${_ids_label_view[$x]}"
		echo "goto $(grep -E -m 1 $'\t'gene$'\t'".+gene_id \"$i\"" "$gtf" | awk '{print $1":"$4"-"$5}')" >> "$f"
		if [[ $delay -gt 0 ]]; then
			cat <<- EOF >> $f
				setSleepInterval $delay
				echo
				setSleepInterval 0
			EOF
		fi
		if $snapshots; then
			cat <<- EOF >> "$f"
				snapshot "${l:+$l.}$i.png"
				snapshot "${l:+$l.}.$i.svg"
			EOF
		fi
	done
	for x in "${!_pos_view[@]}"; do
		i="${_pos_view[$x]}"
		l="${_pos_label_view[$x]}"
		echo "goto $i" >> "$f"
		if [[ $delay -gt 0 ]]; then
			cat <<- EOF >> "$f"
				setSleepInterval $delay
				echo
				setSleepInterval 0
			EOF
		fi
		if $snapshots; then
			cat <<- EOF >> "$f"
				collapse
				snapshot "${l:+$l.}$i.png"
				snapshot "${l:+$l.}$i.svg"
			EOF
		fi
	done
	$autoexit && echo "exit" >> "$f"

	declare -a cmd1
	commander::makecmd -a cmd1 -s ';' -c {COMMANDER[0]}<<- CMD
		java
			--module-path="\$CONDA_PREFIX/lib/igv"
			-Xmx${memory}m
			@"\$CONDA_PREFIX/lib/igv/igv.args"
			-Dapple.laf.useScreenMenuBar=true
			-Djava.net.preferIPv4Stack=true
			--module=org.igv/org.broad.igv.ui.Main
			--igvDirectory "$outdir/igv"
			--batch "$outdir/igv/igv.batch"
	CMD

	commander::runcmd -c igv -v -a cmd1

	return 0
}
