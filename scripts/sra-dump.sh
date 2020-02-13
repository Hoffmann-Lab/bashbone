#! /usr/bin/env bash
# (c) Konstantin Riege
trap 'kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT
trap 'exit $?' INT TERM
trap 'exit 1' USR1

die(){
	#in contrast to exti, kill works also if triggered within subshell
	echo "UNFORSEEN ERROR" >&2
	kill -USR1 $$
}

[[ ! $1 || $1 =~ ^- || $(which esearch >/dev/null; echo $?) -gt 0 || $(which fastq-dump >/dev/null; echo $?) -gt 0 || $(which curl >/dev/null; echo $?) -gt 0 || $(which wget >/dev/null; echo $?) -gt 0 ]] && {
	cat <<- EOF
		SYNOPSIS
		  $(basename $0) retrieves fastq data from ncbi sra based on GSM, SRR or SRX accession numbers
		  - initially, information on given accession numbers will be shown
		  - downloads will be saved to current working directory ($PWD)
		  - up to four files will be downloaded in parallel
		  - if installed, sra-toolkit via ncbi gov resource is priorized over wget from ebi uk mirror
		REQUIREMENTS
		  - esearch (from eutilities)
		  - fastq-dump (optional, from stra-toolkit)
		  - curl
		  - wget
		USAGE
		  $(basename $0) [GSM|SRR|SRX] [..]
		EXAMPLE
		  $(basename $0) GSM1446883 SRR1528586 SRX663213
		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 1
}

declare -a srr title

for i in $(seq 1 $#); do
	id=${!i}
	[[ $id =~ ^SRR ]] && {
		echo -e "$id\t$id" >&2
		srr+=($id)
	} || {
		i=${#srr[@]}
		srr+=($(esearch -db sra -query $id | efetch --format docsum | grep -oE 'SRR[^"]+'))
		n=$(esearch -db sra -query $id | efetch --format info | grep -oE 'sample_title="[^"]+' | cut -d '"' -f 2 | sort -Vu | xargs -echo)
		printf "$id\t%s\t$n\n" ${srr[@]:$i} >&2
	}
done

sratoolkit(){
	for id in ${srr[@]}; do
		echo "fastq-dump --split-3 --defline-seq '@\$ac.\$si[ \$sg ][ \$sn]/\$ri' --defline-qual '+' --gzip -O $PWD -A $id" >&2
		echo -ne "fastq-dump --split-3 --defline-seq '@\$ac.\$si[ \$sg ][ \$sn]/\$ri' --defline-qual '+' --gzip -O $PWD -A $id\0"
	done | xargs -0 -P 4 -I {} bash -c {} || return 1
}

ftpdownload(){
	for id in ${srr[@]}; do
		#echo "wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N --glob on ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${srr:0:6}/00${id:9}/$id/$id*.fastq.gz"
		#echo -ne "wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N --glob on ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${srr:0:6}/00${id:9}/$id/$id*.fastq.gz\0"
		echo "wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex '$id.*\.fastq\.gz' ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/00${id:9}/$id/" >&2
		echo -ne "wget -q --show-progress --progress=bar:force --waitretry 1 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex '$id.*\.fastq\.gz' ftp://ftp.sra.ebi.ac.uk/vol1/fastq/${id:0:6}/00${id:9}/$id/\0"
	done | xargs -0 -P 4 -I {} bash -c {} || return 1
}

#redirect stderr to subshell via process substitution and filter for failed SRR ids to print them to stdout instead
# | print SSR ids to stdout and rest to stderr -> array of to redo ids
[[ $(fastq-dump &> /dev/null; echo $?) -eq 0 ]] && {
	srr=($(sratoolkit 2> >(tee /dev/fd/2 | awk -v x=$(basename $0) '/failed SRR[0-9]+$/{print "\nDONT WORRY! "x" WILL RETRY TO DOWNLOAD "$NF" FROM A DIFFERNT SOURCE" > "/dev/fd/2"; print $NF}') | awk '{if(/^SRR[0-9]+$/){print}else{print > "/dev/fd/2"}}'))
} || {
	ftpdownload || kill -USR1 $$
}

exit 0
