#! /usr/bin/env bash
# (c) Konstantin Riege
trap 'cleanup; kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT
trap 'exit $?' INT TERM
trap 'exit 1' USR1

cleanup(){
	[[ "$tmp" ]] && {
		rm -rf "$tmp"
		[[ "$cfg" ]] && {
			sed -i -r "s@/repository/user/main/public/root\s*=.+@$cfg@" "$HOME/.ncbi/user-settings.mkfg"
		} || {
			sed -i '$ d' "$HOME/.ncbi/user-settings.mkfg"
		}
	}
}

die(){
	#in contrast to exit, kill works also if triggered within a subshell - $(sratoolkit [..])
	echo -e "UNFORSEEN ERROR" >&2
	kill -USR1 $$
}

usage(){
	cat <<- EOF
		SYNOPSIS
		  $(basename "$0") retrieves fastq data from ncbi sra based on GSM, SRR or SRX accession numbers
		  - initially, information on given accession numbers will be shown
		  - downloads will be saved to current working directory ("$PWD")
		  - support for parallel download instances
		  - if installed, sra-toolkit via ncbi gov resource is priorized over retrieval via ebi uk mirror
		  - ebi uk mirror is used as fallback upon fastq-dump errors (not true for fasterq-dump)
		REQUIREMENTS
		  - esearch (from eutilities https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/)
		  - fastq-dump/fasterq-dump (optional, from stra-toolkit https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/)
		  - wget
		USAGE
		  $(basename "$0") [OPTIONS] [GSM|SRR|SRX|SRP] [GSM|SRR|SRX|SRP] [..]
		OPTIONS
		  -p [num]  : optional, number of maximum parallel download instances (2)
		  -t [path] : optional, path to temporary directory ("$PWD")
		  -e        : optional, switch to ebi uk mirror utilizing wget
		  -f        : optional, switch to fasterq-dump 
		              NOTE: not yet recommended!
		                - uncompressed output
		                - misconfigured read deflines
		EXAMPLE
		  $(basename "$0") -p 4 -t /dev/shm GSM1446883 SRR1528586 SRX663213
		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 0
}

unset OPTIND
while getopts p:t:feh ARG; do
	case $ARG in
		p) instances="$OPTARG";;
		t) tmp="$OPTARG";;
		f) faster=true;;
		e) ebi=true;;
		h|*) usage;;
	esac
done
shift $((OPTIND-1))

[[ $# -eq 0 || -z $1 || $(command -v esearch > /dev/null; echo $?) -gt 0 || ( $(command -v fastq-dump > /dev/null; echo $?) -gt 0 && $(command -v wget > /dev/null; echo $?) -gt 0 ) ]] && {
	usage
}

declare -a srr title
instances=${instances:-2}
threads=${threads:-2}
tmp="$(mktemp -d -p "${tmp:-$PWD}")"
faster=${faster:-false}
ebi=${ebi:-false}

[[ -s "$HOME/.ncbi/user-settings.mkfg" ]] && {
	cfg="$(grep -E '/repository/user/main/public/root\s*=.+' "$HOME/.ncbi/user-settings.mkfg")"
	[[ "$cfg" ]] && {
		sed -i -r "s@(/repository/user/main/public/root\s*=\s*)(.+)@\1'$tmp'@" "$HOME/.ncbi/user-settings.mkfg"
	} || echo "/repository/user/main/public/root = $tmp" >> "$HOME/.ncbi/user-settings.mkfg"
} || {
	mkdir -p $HOME/.ncbi
	echo "/repository/user/main/public/root = '$tmp'" > "$HOME/.ncbi/user-settings.mkfg"
}

for i in $(seq 1 $#); do
	id="${!i}"
	[[ "$id" =~ ^SRR ]] && {
		echo -e "$id\t$id" >&2
		srr+=("$id")
	} || {
		i="${#srr[@]}"
		srr+=("$(esearch -db sra -query "$id" | efetch --format docsum | grep -oE 'SRR[^"]+')")
		n=$(esearch -db sra -query "$id" | efetch --format info | grep -oE 'sample_title="[^"]+' | cut -d '"' -f 2 | sort -Vu | xargs echo)
		printf "$id\t%s\t$n\n" "${srr[@]:$i}" >&2
	}
done

fastqdump(){
	cmd="fastq-dump --split-3"
	[[ $($cmd 2>&1 | grep -c unrecognized) -gt 0 ]] && cmd="fastq-dump --split-e"
	for id in "${srr[@]}"; do
		echo "$cmd --defline-seq '@\$ac.\$si[ \$sg ][ \$sn]/\$ri' --defline-qual '+' --gzip -O '$PWD' '$id'" >&2
		echo -ne "$cmd --defline-seq '@\$ac.\$si[ \$sg ][ \$sn]/\$ri' --defline-qual '+' --gzip -O '$PWD' '$id'\0"
	done | xargs -0 -P $instances -I {} bash -c {} || return 1
}

fasterqdump(){
	for id in "${srr[@]}"; do
		echo "fasterq-dump -t $tmp -p -f -P -O '$PWD' '$id'" >&2
		echo -ne "fasterq-dump -t $tmp -p -f -P -O '$PWD' '$id'\0"
	done | xargs -0 -P $instances -I {} bash -c {} || return 1
}

ftpdump(){
	for id in "${srr[@]}"; do
		#echo "wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N --glob on ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${srr:0:6}/00${id:9}/$id/$id*.fastq.gz"
		#echo -ne "wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N --glob on ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${srr:0:6}/00${id:9}/$id/$id*.fastq.gz\0"
		echo "wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex '$id.*\.fastq\.gz' ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/00${id:9}/$id/" >&2
		echo -ne "wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex '$id.*\.fastq\.gz' ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/00${id:9}/$id/\0"
	done | xargs -0 -P $instances -I {} bash -c {} || return 1
}

#redirect stderr to subshell via process substitution and filter for failed SRR ids to print them to stdout instead
# | print SSR ids to stdout and rest to stderr -> array of to redo ids
$faster && [[ $(command -v fasterq-dump > /dev/null; echo $?) -eq 0 ]] && {
	fasterqdump && exit 0 || die
}
$ebi || [[ $(command -v fastq-dump > /dev/null; echo $?) -gt 0 ]] && {
	ftpdump && exit 0 || die
}
srr=("$(fastqdump 2> >(tee /dev/fd/2 | awk -v x="$(basename "$0")" '/failed SRR[0-9]+$/{print "\nDONT WORRY! "x" WILL RETRY TO DOWNLOAD "$NF" FROM A DIFFERNT SOURCE" > "/dev/fd/2"; print $NF}') | awk '{if(/^SRR[0-9]+$/){print}else{print > "/dev/fd/2"}}')")
ftpdump || die

exit 0
