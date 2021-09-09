#! /usr/bin/env bash
# (c) Konstantin Riege
shopt -s extglob

die(){
	#in contrast to exit, kill works also if triggered within a subshell - $(sratoolkit [..])
	echo -e "$*" >&2
	kill -USR1 $$
}

cleanup(){
	[[ -e $tmp ]] && rm -rf "$tmp"
	[[ $mkfg ]] && printf '%s' "${mkfg[@]}" > "$HOME/.ncbi/user-settings.mkfg"
}

trap 'exit $?' INT TERM
trap 'exit 1' USR1

trap '
	cleanup
	sleep 1
	declare -a pids=($(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+" | tail -n +2))
	{ kill -KILL "${pids[@]}" && wait "${pids[@]}"; } &> /dev/null
	printf "\r"
' EXIT

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename "$0") retrieves fastq data from ncbi sra based on GSM, SRR or SRX accession numbers
		  - initially, information for given accession numbers will be printed to stderr
		  - support for parallel download instances
		  - if installed, sra-toolkit via ncbi gov resource is priorized over retrieval via ebi uk mirror
		  - ebi uk mirror is used as fallback upon fastq-dump errors

		VERSION
		0.1.6

		REQUIREMENTS
		Depends on chosen options
		  - esearch (from eutilities https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/)
		  - pigz (for threads > 1) and fastq-dump/fasterq-dump (from stra-toolkit https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/)
		  - wget or curl

		SYNOPSIS
		$(basename "$0") [OPTIONS] [GSM|SRR|SRX|SRP] [GSM|SRR|SRX|SRP] [..]

		OPTIONS
		-s        : show received information for given accession numbers and exit
		-o [file] : additionally print information for given accession numbers to file
		-d [path] : download into this directory
		-p [num]  : number of maximum parallel download instances (default: 2)
		-n        : no ebi mirror fallback upon fastq-dump failure
		-t [num]  : pigz compression threads for fastq-dump (default: no pigz)
		-m [path] : path to temporary directory (default: "$PWD")
		-e        : priorize ebi uk mirror utilizing wget
		            HINT: use -p 1 in a second run to ensure all files were downloaded correctly
		-c        : experimental! priorize ebi uk mirror utilizing curl
		-f        : experimental! switch to fasterq-dump
		            NOTE: not yet recommended!
		              - uncompressed output
		              - misconfigured read deflines

		EXAMPLE
		  $(basename "$0") -p 4 -t 4 -m /dev/shm GSM1446883 SRR1528586 SRX663213

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 1
}

unset OPTIND
while getopts o:p:m:t:d:sfecnh ARG; do
	case $ARG in
		o) outfile="$OPTARG";;
		p) instances="$OPTARG";;
		t) threads="$OPTARG";;
		m) t="$OPTARG";;
		d) outdir="$OPTARG";;
		f) faster=true;;
		e) ebi=true; method=wget;;
		c) ebi=true; method=curl;;
		s) nodownload=true;;
		n) nofallback=true;;
		h) (usage); exit 0;;
		*) usage;;
	esac
done
shift $((OPTIND-1))

[[ $# -eq 0 || -z $1 || $(command -v esearch > /dev/null; echo $?) -gt 0 || ( $(command -v fastq-dump > /dev/null; echo $?) -gt 0 && $(command -v wget > /dev/null; echo $?) -gt 0 ) ]] && {
	usage
}

declare -a srr title
instances=${instances:-2}
threads=${threads:-1}
tmp="$(mktemp -d -p "${t:-$PWD}" tmp.XXXXXXXXXX.sradump)"
faster=${faster:-false}
ebi=${ebi:-false}
method=${method:-wget}
nodownload=${nodownload:-false}
nofallback=${nofallback:-false}
outfile="${outfile:-/dev/null}"
outdir="${outdir:-$PWD}"
mkdir -p "$outdir" || die "ERROR cannot create output directory"
mkdir -p $(dirname "$outfile") || die "ERROR cannot create directory for output file"
[[ "$outfile" != "/dev/null" ]] && rm -f $outfile
touch $outfile || die "ERROR cannot create output file"
resume=true

for i in $(seq 1 $#); do
	id="${!i}"
	# [[ "$id" =~ ^SRR ]] && {
	# 	echo -e "$id\t$id" >&2
	# 	srr+=("$id")
	# } || {
		i="${#srr[@]}"
		srr+=($(esearch -db sra -query "$id" | efetch --format docsum | grep -oE '(S|E|D)RR[^"]+'))
		n=$(esearch -db sra -query "$id" | efetch --format info | grep -oE 'sample_title="[^"]+' | cut -d '"' -f 2 | sort -Vu | xargs echo)
		printf "$id\t%s\t$n\n" "${srr[@]:$i}" | tee -a "$outfile" >&2
	# }
done

$nodownload && exit 0

! $ebi || ! $nofallback && {
	mkdir -p $HOME/.ncbi
	touch "$HOME/.ncbi/user-settings.mkfg"
	mapfile mkfg < "$HOME/.ncbi/user-settings.mkfg"
	rm "$HOME/.ncbi/user-settings.mkfg"


	#vdb-config -i & sleep 2; kill $! # does not work when piped
	vdb-config --restore-defaults &> /dev/null
	vdb-config -s /LIBS/GUID=$(uuid) # from conda ossuuid (dependency of sra-toolkit) or use systems uuidgen if installed
	vdb-config -s /http/timeout/read=20000
	vdb-config -s /repository/user/main/public/root="$tmp" # change cache temp directory
	vdb-config -s /repository/user/cache-disabled=true # comment to enable caching
}

fastqdump_sngl(){
	local id cmd="fastq-dump --split-3"
	[[ $($cmd 2>&1 | grep -c unrecognized) -gt 0 ]] && cmd="fastq-dump --split-e"
	for id in ${srr[@]}; do
		# may add Not-filtered and 0-control-bits and barcode $ri[:N:0:$sg]
		# with accession, id and if available cellid:lane:tile:x:y
		echo "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --gzip -O '$outdir' $id" >&2
		echo -ne "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --gzip -O '$outdir' $id\0"
	done | xargs -0 -P $instances -I {} bash -c {} || return 1
	return 0
}

fastqdump_mult(){
	local id cmd="fastq-dump --split-3"
	[[ $($cmd 2>&1 | grep -c unrecognized) -gt 0 ]] && cmd="fastq-dump --split-e"
	for id in ${srr[@]}; do
		if [[ $($cmd --defline-seq '@$ac.$si[.$sn] \$ri' --defline-qual '+' -X 1 --stdout "$id" 2>/dev/null | wc -l || return 1) -gt 4 ]]; then
			echo "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout '$id' | paste - - - - | tee >(sed -n '1~2{s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/${id}_1.fastq.gz') >(sed -n '2~2{s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/${id}_2.fastq.gz') > /dev/null | cat" >&2
			echo -ne "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout $id | paste - - - - | tee >(sed -n '1~2{s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/${id}_1.fastq.gz') >(sed -n '2~2{s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/${id}_2.fastq.gz') > /dev/null | cat\0"
		else
			echo "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout $id > >(pigz -p $threads -c > '$outdir/$id.fastq.gz') | cat" >&2
			echo -ne "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout $id > >(pigz -p $threads -c > '$outdir/$id.fastq.gz') | cat\0"
		fi
	done | xargs -0 -P $instances -I {} bash -c {} || return 1
	return 0
}

fasterqdump(){
	local id
	for id in ${srr[@]}; do
		echo "fasterq-dump -t $tmp -p -f -P -O '$outdir' '$id'" >&2
		echo -ne "fasterq-dump -t $tmp -p -f -P -O '$outdir' '$id'\0"
	done | xargs -0 -P $instances -I {} bash -c {} || return 1
	return 0
}

ftpdump_wget(){
	local params id url i=-1
	$resume && params="-c" || params=""
	for id in ${srr[@]}; do # do not quote. in case srr=("$(fastqdump ...)") terminates succesfully, srr==("") -> id==""
		url=$([[ $(echo -n $id | wc -c) -lt 10 ]] && echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$id || echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$(printf '%03i' $(echo ${id:9} | sed 's/^0*//'))/$id)
		# attemp to tackle colliding .listing files by sleep
		echo "wget $params -q -P '$outdir' --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --recursive --no-directories --no-parent --level=1 --reject 'index.htm*' --accept-regex '$id.*\.fastq\.gz' '$url/'" >&2
		echo -ne "sleep $((++i%instances*2)); wget $params -q -P '$outdir' --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --recursive --no-directories --no-parent --level=1 --reject 'index.htm*' --accept-regex '$id.*\.fastq\.gz' '$url/'\0"
		# --glob has the same listing problem and additionally adds risk for 404 due to wildcards not supported in HTTP (despite of ftp url, where glob should work)
		# echo "wget $params -q -P '$outdir' --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on '$url/$id*.fastq.gz'" >&2
		# echo -ne "wget $params -q -P '$outdir' --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping '$url/$id*.fastq.gz'\0"
	done | xargs -0 -P $instances -I {} bash -c {} || return 1
	return 0
}

ftpdump_curl(){
	local id params url
	$resume && params="-C -" || params=""
	for id in ${srr[@]}; do # do not quote. in case srr=("$(fastqdump ...)") terminates succesfully, srr==("") -> id==""
		url=$([[ $(echo -n $id | wc -c) -lt 10 ]] && echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$id || echo ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/$(printf '%03i' $(echo ${id:9} | sed 's/^0*//'))/$id)
		echo "curl $params --progress-bar --retry-connrefused --connect-timeout 60 --retry-delay 10 --retry 10 --create-dirs -o '$outdir/#1#2#3' '$url/{$id}{,_1,_2}{.fastq.gz}'" >&2
		echo -ne "curl $params --progress-bar --retry-connrefused --connect-timeout 60 --retry-delay 10 --retry 10 --create-dirs -o '$outdir/#1#2#3' '$url/{$id}{,_1,_2}{.fastq.gz}'\0"
	done | xargs -0 -P $instances -I {} bash -c {}
	# no colliding .listing files, but since at least one of the three files of {$id}{,_1,_2} will be always missed, curl throws errors (in case the latter one was missing).
	# thus do not "|| return 1" unless logged and parsed or use --parallel (--parallel-max $instances) which always succeeds. when used, one may completely replace xargs by curl --parallel [..] -o $o1 $url1 -o $o2 $url2
	# though, a log parser cannot classify error 78 aka server 550 aka file not into a true server error or an accepted misbehaviour by {$id}{,_1,_2}
	return 0
}

$faster && [[ $(command -v fasterq-dump > /dev/null; echo $?) -eq 0 ]] && {
	fasterqdump && exit 0 || die "UNFORSEEN ERROR"
}

$ebi || [[ $(command -v fastq-dump > /dev/null; echo $?) -gt 0 ]] && {
	ftpdump_$method && exit 0 || die "UNFORSEEN ERROR"
}

# redirect stderr to subshell via process substitution to be used as stdout
# use tee to directly print stdout again as stderr and further filter stdout for failed SRR ids via awk
# print message in awk to stderr and check if accession number is truely SRR by second awk command utilizing regex
# second awk print SSR ids to stdout or other ID to stderr
# finally, srr array holds SRR numbers to be re-downloaded via ftpdump
compression=$([[ $threads -eq 1 ]] && echo sngl || echo mult)
$nofallback && {
	fastqdump_$compression && exit 0 || die "UNFORSEEN ERROR"
} || {
	srr=("$(fastqdump_$compression 2> >(tee /dev/fd/2 | awk -v x="$(basename "$0")" '/failed (S|E|D)RR[0-9]+$/{print "\nDONT WORRY! "x" WILL RETRY TO DOWNLOAD "$NF" FROM A DIFFERNT SOURCE" > "/dev/fd/2"; print $NF}') | awk '{if(/^(S|E|D)RR[0-9]+$/){print}else{print > "/dev/fd/2"}}')")
	resume=false
	ftpdump_$method && exit 0 || die "UNFORSEEN ERROR"
}
