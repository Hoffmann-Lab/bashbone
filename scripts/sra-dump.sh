#! /usr/bin/env bash
# (c) Konstantin Riege

set -o pipefail -o errtrace # for ERR trap inheritance
# -o functrace for DEBUG and RETURN trap inheritance in case of err backtrace implementation

# shopt -s extdebug # for advanced LINENO handling via read -r fun line src < <(declare -F ${FUNCNAME[0]})
shopt -s extglob # e.g. +(a|b|c)*
# shopt -s expand_aliases # to allow aliases in scripts

cleanup(){
	[[ $tmp && -e "$tmp" ]] && rm -rf "$tmp"
	[[ $mkfg ]] && printf '%s' "${mkfg[@]}" > "$HOME/.ncbi/user-settings.mkfg"
	return 0
}

trap '
	cleanup
	declare -a pids=($(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+" | tail -n +2))
	{ kill -INT "${pids[@]}" && wait "${pids[@]}"; } &> /dev/null
	printf "\r"
' EXIT

trap 'exit $?' INT TERM
trap 'exit 1' USR1
trap 'e=$?; [[ $ERROR ]] && echo "ERROR $ERROR"; if [[ $e -ne 141 ]]; then [[ $BASHPID -eq $$ ]] && exit $e || kill -USR1 $$; fi' ERR
# ignore SIGPIPE caused by e.g. 'samtools view bam | head' and take care of parent exit upon subshell ERR
# to trigger ERR during xargs runtime let jobs exit with 255 e.g. echo 'CMD || exit 255' | xargs -I {} bash -c {}

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename "$0") retrieves fastq data from ncbi sra based on GSM, SRR or SRX accession numbers
		  - initially, information for given accession numbers will be printed to stderr
		  - support for parallel download instances
		  - if installed, sra-toolkit via ncbi gov resource is priorized over retrieval via ebi uk mirror
		  - ebi uk mirror is used as fallback upon fastq-dump errors

		VERSION
		0.2.0

		REQUIREMENTS
		Depends on chosen options
		  - esearch (from eutilities https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/)
		  - pigz (for threads > 1) and fastq-dump/fasterq-dump (from stra-toolkit https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/)
		  - wget or curl

		SYNOPSIS
		$(basename "$0") [OPTIONS] [GSM|SRR|SRX|SRP] [..]
		$(basename "$0") -f [OPTIONS] [SRA] [..]

		OPTIONS
		-d [path] : download into this directory (default: "$PWD")
		-p [num]  : number of maximum parallel download instances (default: 2)
		-t [num]  : pigz compression threads per fastq-dump download instance (default: no pigz)
		-m [path] : path to temporary directory (default: "$PWD/tmp.XXXXXXXXXX.sradump")
		-f        : convert local sra files to compressed fastq files
		-s        : show received meta information for given accession numbers and exit
		-r        : do not fetch meta information for given SRR* accession numbers
		-o [file] : additionally print information for given accession numbers to file
		-a [list] : fetch all, biological and technical reads, and reassign mate ids according to a comma separated list e.g. 1,3,2
		-b        : use ebi uk mirror fallback utilizing wget (unless -c) upon fastq-dump failures
		-e        : unless -a, priorize ebi uk mirror utilizing wget
		            HINT1: outperformes fastq-dump but uses less stable connections which may requires additional runs
		            HINT2: use -p 1 in a final run to ensure all files were downloaded correctly
		-c        : experimental! unless -a, priorize ebi uk mirror utilizing curl

		EXAMPLES
		$(basename "$0") -p 2 -t 4 -m /dev/shm GSM1446883 SRR1528586 SRX663213
		$(basename "$0") -f -p 2 -t 4 -m /dev/shm /path/to/SRR1528586.sra /path/to/SRR1528588.sra

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 1
}

unset OPTIND
while getopts o:p:m:t:d:a:srecbhf ARG; do
	case $ARG in
		o) outfile="$OPTARG";;
		p) instances=$OPTARG;;
		t) threads=$OPTARG;;
		m) t=$OPTARG;;
		a) mapfile -t -d ',' mateid < <(printf "%s" "$OPTARG");;
		d) outdir="$OPTARG";;
		f) files=true;;
		e) ebi=true; method=wget;;
		c) ebi=true; method=curl;;
		s) nodownload=true;;
		r) outfile=/dev/null; nofetch=true;;
		b) ebi=true; fallback=true;;
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
files=${files:-false}
$files && ebi=false && fallback=false
ebi=${ebi:-false}
[[ $mateid ]] && ebi=false && fallback=false
method=${method:-wget}
nodownload=${nodownload:-false}
nofetch=${nofetch:-false}
fallback=${fallback:-false}
outfile="${outfile:-/dev/null}"
outdir="${outdir:-$PWD}"
ERROR="cannot create output directory"
mkdir -p "$outdir"
ERROR="cannot create directory for output file"
mkdir -p "$(dirname "$outfile")"
[[ "$outfile" != "/dev/null" ]] && rm -f "$outfile"
ERROR="cannot create output file"
touch "$outfile"
resume=true
unset ERROR

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
			i="${#srr[@]}"
			srr+=($(esearch -db sra -query "$id" | efetch --format docsum | grep -oE '(S|E|D)RR[^"]+'))
			n=$(esearch -db sra -query "$id" | efetch --format info | grep -oE 'sample_title="[^"]+' | cut -d '"' -f 2 | sort -Vu | xargs echo)
			printf "$id\t%s\t$n\n" "${srr[@]:$i}" | tee -a "$outfile" >&2
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
	vdb-config -s /LIBS/GUID=$(uuid) # from conda ossuuid (dependency of sra-toolkit) or use systems uuidgen if installed
	vdb-config -s /http/timeout/read=20000
	vdb-config -s /repository/user/main/public/root="$tmp" # change cache temp directory
	vdb-config -s /repository/user/cache-disabled=true # comment to enable caching
}

fastqdump_sngl(){
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
	for id in "${srr[@]}"; do
		if [[ $($cmd --defline-seq '@$ac.$si[.$sn] \$ri' --defline-qual '+' -X 1 --stdout "$id" 2>/dev/null | wc -l || return 1) -gt 4 ]]; then
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
	local id cmd="fastq-dump --split-files"
	for id in "${srr[@]}"; do
		if [[ $($cmd --defline-seq '@$ac.$si[.$sn] \$ri' --defline-qual '+' -X 1 --stdout "$id" 2>/dev/null | wc -l || return 1) -gt 6 ]]; then
			echo "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout '$id' | paste - - - - | tee >(sed -n '1~3{s/[0-9]\\\t/${mateid[0]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+2)/3)) -c > '$outdir/$(basename "$id" .sra)_${mateid[0]}.fastq.gz') >(sed -n '2~3{s/[0-9]\\\t/${mateid[1]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+2)/3)) -c > '$outdir/$(basename "$id" .sra)_${mateid[1]}.fastq.gz') >(sed -n '3~3{s/[0-9]\\\t/${mateid[2]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+2)/3)) -c > '$outdir/$(basename "$id" .sra)_${mateid[2]}.fastq.gz') > /dev/null | cat" | tr -s '\\' '\\' >&2
			echo -ne "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout '$id' | paste - - - - | tee >(sed -n '1~3{s/[0-9]\\\t/${mateid[0]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+2)/3)) -c > '$outdir/$(basename "$id" .sra)_${mateid[0]}.fastq.gz') >(sed -n '2~3{s/[0-9]\\\t/${mateid[0]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+2)/3)) -c > '$outdir/$(basename "$id" .sra)_${mateid[1]}.fastq.gz') >(sed -n '3~3{s/[0-9]\\\t/${mateid[2]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+2)/3)) -c > '$outdir/$(basename "$id" .sra)_${mateid[2]}.fastq.gz') > /dev/null | cat\0"
		else
			echo "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout $id | paste - - - - | tee >(sed -n '1~2{s/[0-9]\\\t/${mateid[0]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/$(basename "$id" .sra)_${mateid[0]}.fastq.gz') >(sed -n '2~2{s/[0-9]\\\t/${mateid[1]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/$(basename "$id" .sra)_${mateid[1]}.fastq.gz') > /dev/null | cat" | tr -s '\\' '\\' >&2
			echo -ne "$cmd --defline-seq '@\$ac.\$si[.\$sn] \$ri' --defline-qual '+' --stdout $id | paste - - - - | tee >(sed -n '1~2{s/[0-9]\\\t/${mateid[0]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/$(basename "$id" .sra)_${mateid[0]}.fastq.gz') >(sed -n '2~2{s/[0-9]\\\t/${mateid[1]}\\\t/;s/\\\t/\\\n/gp}' | pigz -p $(((threads+1)/2)) -c > '$outdir/$(basename "$id" .sra)_${mateid[1]}.fastq.gz') > /dev/null | cat\0"
		fi
	done | xargs -0 -P $instances -I {} bash -c {}
	return 0
}

fasterqdump(){
	local id
	for id in "${srr[@]}"; do
		echo "fasterq-dump -t '$tmp' -p -f -P -O '$outdir' '$id'" >&2
		echo -ne "fasterq-dump -t '$tmp' -p -f -P -O '$outdir' '$id'\0"
	done | xargs -0 -P $instances -I {} bash -c {}
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

$faster && [[ $(command -v fasterq-dump > /dev/null; echo $?) -eq 0 ]] && {
	ERROR="@ fasterqdump"
	fasterqdump
	exit 0
}

$ebi || [[ $(command -v fastq-dump > /dev/null; echo $?) -gt 0 ]] && {
	ERROR="@ ftpdump_$method"
	ftpdump_$method
	exit 0
}

compression=$([[ $threads -eq 1 ]] && echo sngl || echo mult)
$fallback && {
	# redirect stderr to subshell via process substitution to be used as stdout
	# use tee to directly print stdout again as stderr and further filter stdout for failed SRR ids via awk
	# print message in awk to stderr and check if accession number is truely SRR by second awk command utilizing regex
	# second awk print SSR ids to stdout or other ID to stderr
	# finally, srr array holds SRR numbers to be re-downloaded via ftpdump
	srr=("$({ fastqdump_$compression || true; } 2> >(tee /dev/fd/2 | awk -v x="$(basename "$0")" '/failed (S|E|D)RR[0-9]+$/{print "\nDONT WORRY! "x" WILL RETRY TO DOWNLOAD "$NF" FROM A DIFFERNT SOURCE" > "/dev/fd/2"; print $NF}') | awk '{if(/^(S|E|D)RR[0-9]+$/){print}else{print > "/dev/fd/2"}}')")
	resume=false
	ERROR="@ ftpdump_$method during rescuing of fastqdump_$compression"
	ftpdump_$method
	exit 0
} || {
	[[ $mateid ]] && {
		ERROR="@ fastqdump_all"
		fastqdump_all
		exit 0
	} || {
		ERROR="@ fastqdump_$compression"
		fastqdump_$compression
		exit 0
	}
}

exit 0
