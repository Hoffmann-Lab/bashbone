#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$0")/../bashbone_lite.sh" -x cleanup -a "$@" || exit 1

############################################

cleanup(){
	rm -rf "$tmp"
	[[ $mkfg ]] && printf '%s' "${mkfg[@]}" > "$HOME/.ncbi/user-settings.mkfg"
	return 0
}

############################################

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename "$0") retrieves fastq data from ncbi sra based on GSE, GSM, SRR or SRX accession numbers
		  - support for parallel download instances
		  - if installed, sra-toolkit via ncbi gov resource is priorized over retrieval via ebi uk mirror
		  - ebi uk mirror can be defined as fallback upon fastq-dump errors

		VERSION
		0.6.0

		REQUIREMENTS
		Depends on chosen options
		  - esearch (from eutilities https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/)
		  - pigz (for compression threads > 1)
		  - wget or curl or fastq-dump (from stra-toolkit https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/)

		SYNOPSIS
		$(basename "$0") [OPTIONS] [GSE|GSM|SRR|SRX|SRP] [..]
		$(basename "$0") -l [OPTIONS] [SRR.sra] [..]

		OPTIONS
		-d [path] : download into this directory (default: "$PWD")
		-p [num]  : number of maximum parallel download instances (default: 2)
		-t [num]  : pigz compression threads per fastq-dump download instance (default: no pigz)
		            HINT1: single threaded gzip compression may reduce download speed
		-m [path] : path to temporary directory (default: "$PWD/tmp.XXXXXXXXXX.sradump")
		-l        : convert local sra files to compressed fastq files
		            HINT1: filename must start with SRR accession number followed by an extension
		            HINT2: misuse this parameter iff all accession numbers are SRR ids and eutilities fail due to server errors
		-s        : show received meta information for given accession numbers and exit
		-r        : do not fetch meta information for given SRR* accession numbers
		-o [file] : additionally write meta information for given accession numbers to file
		-a [list] : fetch all, biological and technical/UMI/barcode reads, and reassign read numbers according to given comma separated list
		            HINT1: can create up to 6 files depending on single-end/paired-end with/without UMIs/barcode and single-/dual-index
		            HINT2: example for single-index paired-end data without UMIs assuming 1=R1, 2=index, 3=R2 use: -a 1,3,2
		            HINT3: if usure, simply provide 1,2,3,4,5,6
		-w        : unless -a, download from ebi uk mirror utilizing wget
		            HINT1: outperformes fastq-dump but uses less stable connections which may requires additional runs
		            HINT2: use -p 1 in a final run to ensure all files were downloaded correctly
		-c        : unless -a, download from ebi uk mirror utilizing curl
		            HINT1: experimental!
		-f        : unless -a, use ebi uk mirror as fallback upon fastq-dump failures. requires -w or -c
		-z        : unless -a, download sra file from amazon aws utilizing awscli
		            HINT1: uses -p 1 due to download of multiple chunks in parallel
		            HINT2: use -l in a second run to convert the local sra files

		EXAMPLES
		$(basename "$0") -p 2 -t 4 -m /tmp GSM1446883 SRR1528586 SRX663213
		$(basename "$0") -l -p 2 -t 4 -m /tmp /path/to/SRR1528586.sra /path/to/SRR1528588.sra

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	return 1
}

while getopts o:p:m:t:d:a:srwcfhlz ARG; do
	case $ARG in
		o) outfile="$OPTARG";;
		p) instances=$OPTARG;;
		t) threads=$OPTARG;;
		m) t=$OPTARG;;
		a) mapfile -t -d ',' mateid < <(printf "%s" "$OPTARG");;
		d) outdir="$OPTARG";;
		l) files=true;;
		w) ebi=true; method=wget;;
		c) ebi=true; method=curl;;
		z) aws=true; method=awscli;;
		s) nodownload=true;;
		r) outfile=/dev/null; nofetch=true;;
		f) fallback=true;;
		h) { usage || exit 0; };;
		*) usage;;
	esac
done
shift $((OPTIND-1))
[[ $# -eq 0 ]] && { usage || exit 0; }

declare -a srr title
instances=${instances:-2}
threads=${threads:-1}
tmp="$(mktemp -d -p "${t:-$PWD}" tmp.XXXXXXXXXX.sradump)"
files=${files:-false}
$files && ebi=false && aws=false && fallback=false
ebi=${ebi:-false}
aws=${aws:-false}
[[ $mateid ]] && ebi=false && aws=false && fallback=false
nodownload=${nodownload:-false}
nofetch=${nofetch:-false}
fallback=${fallback:-false}
BASHBONE_ERROR="requires -w or -c"
$fallback && $ebi
outfile="$(realpath -s "${outfile:-/dev/null}")"
outdir="$(realpath -s "${outdir:-$PWD}")"
BASHBONE_ERROR="cannot create output directory"
mkdir -p "$outdir"
BASHBONE_ERROR="cannot create directory for output file"
mkdir -p "$(dirname "$outfile")"
[[ "$outfile" != "/dev/null" ]] && rm -f "$outfile"
BASHBONE_ERROR="cannot create output file"
touch "$outfile"
resume=true
unset BASHBONE_ERROR

if $files; then
	for i in $(seq 1 $#); do
		srr+=("${!i}")
	done
else
	for i in $(seq 1 $#); do
		id="${!i}"
		$nofetch && [[ "$id" =~ ^SRR ]] && {
			echo -e "$id\t$id" >&2
			srr+=("$id")
		} || {
			if esearch -db sra -query "$id" | efetch -format xml | grep -q .; then
				ids=($id)
			else
				# alternative for GEO accessions: Rscript -e 'library(GEOquery); geo_data <- getGEO(geo_accession); info <- geo_data@header[["characteristics_ch1"]];'
				ids=($(esearch -db gds -query "$id" | efetch | grep -oE 'Sample\s*Accession:\s*\S*' 1 | awk '{print $NF}'))
			fi
			for id in "${ids[@]}"; do
				mapfile -t mapdata < <(esearch -db sra -query "$id" | efetch -format xml | grep .)
				mapfile -t mapdata < <(join -t $'\t' <(printf '%s\n' "${mapdata[@]}" | xtract -pattern EXPERIMENT_PACKAGE -element RUN@accession,EXPERIMENT/TITLE,Member@sample_name,Member@sample_title | sed -E ':a;s/([^\t]+)\t\1/\1/;ta' | awk -F '\t' '{if($4){$3=$4}; print $1"\tsample_name=\""$3"\";\tsample_title=\""$2"\";"}') <(printf '%s\n' "${mapdata[@]}" | xtract -pattern EXPERIMENT_PACKAGE -element RUN@accession,SAMPLE_ATTRIBUTE/TAG,SAMPLE_ATTRIBUTE/VALUE | perl -F'\t' -lane 'for (1..$#F/2){$F[$_]=~s/\s+/_/g; $F[$_].="=\"$F[$_+$#F/2]\";"} print join"\t",@F[0..$#F/2]'))
				srr+=($(printf '%s\n' "${mapdata[@]}" | cut -f 1))
				printf '%s\n' "${mapdata[@]}" | sed "s/^/$id\t/" | tee -a "$outfile" >&2
			done
		}
	done
	$nodownload && exit 0
fi

! $ebi || $fallback && {
	mkdir -p $HOME/.ncbi
	touch "$HOME/.ncbi/user-settings.mkfg"
	mapfile mkfg < "$HOME/.ncbi/user-settings.mkfg"
	rm "$HOME/.ncbi/user-settings.mkfg"

	#vdb-config -i & sleep 2; kill $! # does not work when piped
	vdb-config --restore-defaults &> /dev/null
	# from conda ossuuid (dependency of sra-toolkit) or use systems uuidgen if installed
	vdb-config -s /LIBS/GUID=$(uuid)
	vdb-config -s /http/timeout/read=20000
	vdb-config -s /repository/user/main/public/root="$tmp" # change cache temp directory
	vdb-config -s /repository/user/cache-disabled=true # comment to enable caching
}

fastqdump_sngl(){
	# if local file should be converted:
	# filename must be SRRxxxxx.sra because read ids become filename SRRxxxxx and otherwise --defline-seq fails to request origial read ids via web service
	local id cmd="fastq-dump --split-3"
	[[ $($cmd 2>&1 | grep -c unrecognized) -gt 0 ]] && cmd="fastq-dump --split-e"
	for id in "${srr[@]}"; do
		# may add Not-filtered and 0-control-bits and barcode $ri[:N:0:$sg]
		# with accession, id and if available cellid:lane:tile:x:y
		echo "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --gzip -O '$outdir' '$id'" >&2
		echo -ne "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --gzip -O '$outdir' '$id'\0"
	done | xargs -0 -P $instances -I {} bash -c {}
	return 0
}

fastqdump_mult(){
	local id cmd="fastq-dump --split-3"
	[[ $($cmd 2>&1 | grep -c unrecognized) -gt 0 ]] && cmd="fastq-dump --split-e"

	# works, but in case of SE data creates SRR_1.fastq.gz file and runs data stream always through sed
	#fastqdump_all "$cmd"
	#return 0

	for id in "${srr[@]}"; do
		if [[ $($cmd --defline-seq '@$ac.$si[.$sn] $ri' --defline-qual '+' -X 1 --stdout "$id" 2>/dev/null | wc -l || return 1) -eq 8 ]]; then
			# or: fastq-dump | paste - - - - - - - - | tee >(cut -f 1-4 | tr '\t' '\n' | pigz -c > R1) >(cut -f 5-8 | tr '\t' '\n' | pigz > R2)
			echo "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout '$id' | paste - - - - | tee >(sed -n '1~2{s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/$(basename "$id" .sra)_1.fastq.gz') >(sed -n '2~2{s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/$(basename "$id" .sra)_2.fastq.gz') > /dev/null | cat" | tr -s '\\' '\\' >&2
			echo -ne "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout '$id' | paste - - - - | tee >(sed -n '1~2{s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/$(basename "$id" .sra)_1.fastq.gz') >(sed -n '2~2{s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/$(basename "$id" .sra)_2.fastq.gz') > /dev/null | cat\0"
		else
			echo "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout '$id' > >(pigz -p $threads -c > '$outdir/$(basename "$id" .sra).fastq.gz') | cat" | tr -s '\\' '\\' >&2
			echo -ne "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout '$id' > >(pigz -p $threads -c > '$outdir/$(basename "$id" .sra).fastq.gz') | cat\0"
		fi
	done | xargs -0 -P $instances -I {} bash -c {}
	return 0
}

fastqdump_all(){
	# example for single-end with UMIs SRR20073591
	local id i j n ithreads cmd="${1:-"fastq-dump --split-files"}"
	[[ $mateid ]] || mateid=(1 2 3 4 5 6)
	local -a teecmds
	for id in "${srr[@]}"; do
		n=$($cmd --defline-seq '@$ac.$si[.$sn] $ri' --defline-qual '+' -X 1 --stdout "$id" 2>/dev/null | wc -l || return 1)
		ithreads=$((threads/(n/4)))
		[[ $ithreads -eq 0 ]] && ithreads=1;
		teecmds=()
		j=0
		for i in $(seq 1 4 $n); do
			teecmds+=("cut -f $i-$((i+3)) | sed 's/[0-9]\\\t/${mateid[$j]}\\\t/' | tr '\\\t' '\\\n' | pigz -p $ithreads -c > '$outdir/$(basename "$id" .sra)_${mateid[$j]}.fastq.gz'")
			((++j))
		done
		echo "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout '$id' 2>/dev/null | paste $(yes - | head -n $n | xargs echo) | parallel --tmpdir '$tmp' --termseq INT,1000,TERM,0 --halt now,fail=1 --pipe --tee ::: $(printf '"%s" ' "${teecmds[@]}")" | tr -s '\\' '\\' >&2
		echo -ne "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout '$id' 2>/dev/null | paste $(yes - | head -n $n | xargs echo) | parallel --tmpdir '$tmp' --termseq INT,1000,TERM,0 --halt now,fail=1 --pipe --tee ::: $(printf '"%s" ' "${teecmds[@]}")\0"
	done | xargs -0 -P $instances -I {} bash -c {}

	return 0
}

fasterqdump(){
	# from docs: To approach feature parity with fastq-dump, fasterq-dump supports a flexible defline. There are 2 new commandline-parameters for that:
	# --seq-defline FORMAT and --qual-defline FORMAT
	# -> these options don't exist! but if original read ids will be automatically added (space seperated) if available

	# fasterq-dump -t [TMP] --threads 1 --split-spot --print-read-nr --skip-technical --stdout [ID|FILE]
	# -P|--print-read-nr adds /1 /2 ... /5 to read id. better would be @ID 1 and @ID 2 instead of @ID/1 and @ID/2
	# -> to spare parsing dont use!
	# to include technical: --include-technical (example 1-> R1 , 2-> tech , 3-> R2)

	# if local file should be converted:
	# filename must be SRRxxxxx.sra. file must be in working directory
	# -> read ids become filename SRRxxxxx.sra -> create SRRxxxxx directory and symlink file as SRRxxxxx.sra into it

	# runtime benchmark for local SRR24441955 conversion
	# fastq-dump: 1m48s
	# fasterq-dump -t /dev/shm --threads 1: 2m24s
	# fasterq-dump -t /dev/shm --threads 30: 1m53s
	# -> due to tmp usage fasterq-dump is slow and has a high disk footprint

	# runtime benchmark for download and conversion SRR24441955
	# fastq-dump: 13m13s
	# fasterq-dump -t /dev/shm --threads 30: 13m6s
	# -> parity

	# summary:
	# less flexible in defline
	# adds length info to read
	# does not compress
	# does not clean up tmp upon abort
	# tmp directory contains uncompressed fastq data before merged or send to stdout
	# -> aside being able to download multiple files parallel, buth then for sure not to stdout, fasterq-dump sucks

	return 0
}

s3dump_awscli(){
	# wget https://sra-pub-run-odp.s3.amazonaws.com/sra/$id/$id
	local id
	for id in "${srr[@]}"; do
		echo "aws s3 sync 's3://sra-pub-run-odp/sra/$id' '$outdir' --no-sign-request && mv '$outdir/$id' '$outdir/$id.sra'" >&2
		echo -ne "aws s3 sync 's3://sra-pub-run-odp/sra/$id' '$outdir' --no-sign-request && mv '$outdir/$id' '$outdir/$id.sra'\0"
	done | xargs -0 -P $instances -I {} bash -c {}
	return 0
}

ftpdump_wget(){
	local params id url i=-1 tdir
	$resume && params="-c" || params=""
	for id in ${srr[@]}; do # do not quote. in case srr=("$(fastqdump ...)") terminates succesfully, srr==("") -> id==""
		url=$([[ $(echo -n $id | wc -c) -lt 10 ]] && echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$id || echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$(printf '%03i' $(echo ${id:9} | sed 's/^0*//'))/$id)
		# attemp to tackle colliding .listing files at $outdir by sleep
		echo "wget $params -q -P '$outdir' --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --recursive --no-directories --no-parent --level=1 --reject 'index.htm*' --accept-regex '$id.*\.fastq\.gz' '$url/'" >&2
		echo -ne "sleep $((++i%instances*2)); wget $params -q -P '$outdir' --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --recursive --no-directories --no-parent --level=1 --reject 'index.htm*' --accept-regex '$id.*\.fastq\.gz' '$url/'\0"
		# --glob has the same listing problem and additionally adds risk for 404 due to wildcards not supported in HTTP (despite of ftp url, where glob should work)
		# echo "wget $params -q -P '$outdir' --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on '$url/$id*.fastq.gz'" >&2
		# echo -ne "wget $params -q -P '$outdir' --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping '$url/$id*.fastq.gz'\0"
	done | xargs -0 -P $instances -I {} bash -c {}
	return 0
}

ftpdump_curl(){
	local id params url
	$resume && params="-C -" || params=""
	for id in ${srr[@]}; do # do not quote. in case srr=("$(fastqdump ...)") terminates succesfully, srr==("") -> id==""
		url=$([[ $(echo -n $id | wc -c) -lt 10 ]] && echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$id || echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$(printf '%03i' $(echo ${id:9} | sed 's/^0*//'))/$id)
		echo "curl $params --progress-bar --retry-connrefused --connect-timeout 60 --retry-delay 10 --retry 10 --create-dirs -o '$outdir/#1#2#3' '$url/{$id}{,_1,_2}{.fastq.gz}'" >&2
		echo -ne "curl $params --progress-bar --retry-connrefused --connect-timeout 60 --retry-delay 10 --retry 10 --create-dirs -o '$outdir/#1#2#3' '$url/{$id}{,_1,_2}{.fastq.gz}' || true\0"
	done | xargs -0 -P $instances -I {} bash -c {}
	# no colliding .listing files, but since at least one of the three files of {$id}{,_1,_2} will be always missed, curl throws errors (in case the latter one was missing).
	# thus use "|| true" or use --parallel (--parallel-max $instances) which always succeeds. when used, one may completely replace xargs by curl --parallel [..] -o $o1 $url1 -o $o2 $url2
	# though, a log parser cannot classify error 78 aka server 550 aka file not into a true server error or an accepted misbehaviour by {$id}{,_1,_2}
	return 0
}

compression=$([[ $threads -eq 1 ]] && echo sngl || echo mult)
if $fallback; then
	# redirect stderr to subshell via process substitution to be used as stdout
	# use tee to directly print stdout again as stderr and further filter stdout for failed SRR ids via awk
	# print message in awk to stderr and check if accession number is truely SRR by second awk command utilizing regex
	# second awk print SSR ids to stdout or other ID to stderr
	# finally, srr array holds SRR numbers to be re-downloaded via ftpdump
	srr=("$({ fastqdump_$compression || true; } 2> >(tee /dev/fd/2 | awk -v x="$(basename "$0")" '/failed (S|E|D)RR[0-9]+$/{print "\nDONT WORRY! "x" WILL RETRY TO DOWNLOAD "$NF" FROM A DIFFERNT SOURCE" > "/dev/fd/2"; print $NF}') | awk '{if(/^(S|E|D)RR[0-9]+$/){print}else{print > "/dev/fd/2"}}')")
	resume=false
	BASHBONE_ERROR="@ ftpdump_$method during rescuing of fastqdump_$compression"
	ftpdump_$method
	exit 0
else
	if [[ $mateid ]]; then
		BASHBONE_ERROR="@ fastqdump_all"
		fastqdump_all
	elif $ebi; then
		BASHBONE_ERROR="@ ftpdump_$method"
		ftpdump_$method
	elif $aws; then
		BASHBONE_ERROR="@ s3dump_$method"
		s3dump_$method
	else
		BASHBONE_ERROR="@ fastqdump_$compression"
		fastqdump_$compression
	fi

fi

exit 0
