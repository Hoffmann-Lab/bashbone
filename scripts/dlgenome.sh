#! /usr/bin/env bash
# (c) Konstantin Riege
shopt -s extglob
trap 'die' INT TERM
#trap 'kill -PIPE 0' EXIT # kills parental processes as well - shlog conflict
#trap 'kill -PIPE -- -$$' EXIT # kill all childs - works only if $$ is process group leader
trap 'kill -PIPE $(pstree -p $$ | grep -Eo "\([0-9]+\)" | grep -Eo "[0-9]+") &> /dev/null' EXIT # parse pstree
# AVOID DOUBLE FORKS -> run(){cmd &}; run & -> i.e. cmd gets new process group and cannot be killed

die(){
	rm -rf $outdir/tmp
	echo -ne "\e[0;31m"
	echo ":ERROR: $*" | tee -a $log >&2
	echo -ne "\e[m"
	exit 1
}

usage(){
cat <<- EOF
	DESCRIPTION
	$(basename $0) downloads most recent human or mouse genome and annotation including gene ontology and optionally dbSNP

	VERSION
	0.3.1

	SYNOPSIS
	$(basename $0) -g [hg19|hg38|mm10] -s -d

	INPUT OPTIONS
	-o | --out [path]              : output directory
	-g | --genome [hg19|hg38|mm10] : choose GRCh37/hg19 or GRCh38/hg38 or GRCm38/mm10
	-s | --dbsnp                   : additionally download Ensembl dbSNP (very large files)
	-n | --ncbi                    : download NCBI dbSNP instead of Ensembl dbSNP
	-d | --descriptions            : additionally download gene description and ontology information
	-t | --threads [value]         : threads - predicted default: $threads
	-v | --verbose                 : enable verbose mode
	-h | --help                    : prints this message

	REFERENCES
	(c) Konstantin Riege
	konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 0
}

progressbar() {
	mod=0
	while true; do
		((++mod))
		case $mod in
			[15]) echo -en "\r|";;
			[26]) echo -en "\r/";;
			[37]) echo -en "\r-";;
			4) echo -en "\r\\";;
			8) echo -en "\r\\"; mod=0;;
		esac
		sleep 0.2
	done

	return 0
}

progresslog() {
	tail -f $1 2>&1 | grep -E --line-buffered '^\s*(:INFO|:ERROR|Elapsed \(wall|Maximum resident)'

	return 0
}

######
# go #
######

dlgenome::go(){
	[[ ! $go ]] && return 0

	out=$outdir/$1
	export R_LIBS=$outdir/tmp
	mkdir -p $R_LIBS
	dataset="hsapiens_gene_ensembl"
	[[ $g == "mm10" ]] && dataset="mmusculus_gene_ensembl"

	cat <<- EOF > $outdir/tmp/download.R || return 1
		source("https://bioconductor.org/biocLite.R")
		biocLite("biomaRt", suppressUpdates=TRUE)
		library("biomaRt")
		ensembldataset <- "$dataset"
		ensemblversion <- as.numeric(tail(unlist(strsplit(as.character(listMarts()\$version[1]),"\\\s+")),1))
		ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
		ensemblgenes <- listDatasets(ensembl)
		ensemblgenomeversion <- gsub("(^\\\(|\\\)\$)","",tail(unlist(strsplit(as.character(ensemblgenes\$description[ensemblgenes\$dataset==ensembldataset]),"\\\s+")),1))
		ensembl <- useDataset(ensembldataset,mart=ensembl)
		print("downloading datasets, please wait...")
		goids <- getBM(mart=ensembl,attributes=c("ensembl_gene_id","go_id","namespace_1003"))
		descriptions <- getBM(mart=ensembl,attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"))
		write.table(goids,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.go")
		write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")
	EOF
	echo ":INFO: downloading gene ontology information"
	Rscript $outdir/tmp/download.R || return 1

	return 0
}

##########
#  hg38  #
##########

dlgenome::hg38() {
{	cd $outdir && \
	genome=$(curl -s https://www.ensembl.org/info/website/archives/assembly.html | grep -oP 'GRCh38.p[0-9]+' | sort -V | tail -1) && \
	echo ":INFO: downloading $genome genome" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex 'Homo_sapiens.GRCh38.dna.chromosome.(MT|X|Y|[0-9]+).fa.gz$' ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/ && \
	echo ":INFO: extracting chrM" && \
	gzip -dc Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1

	for i in {1..22} X Y; do
		echo ":INFO: extracting chr$i"
		gzip -dc Homo_sapiens.GRCh38.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
		[[ ${PIPESTATUS[0]} -eq 0 ]] || return 1
	done

{	echo ":INFO: downloading annotation" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex 'Homo_sapiens.GRCh38.[0-9]+.chr.gtf.gz$' ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/ && \
	ensembl=$(ls -v Homo_sapiens.GRCh38.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+') && \
	echo ":INFO: extracting annotation" && \
	gzip -dc Homo_sapiens.GRCh38.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};' && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1
	
	rm -f $genome.fa.gtf
	for i in M {1..22} X Y; do
		echo ":INFO: sorting chr$i"
		LC_ALL=C sort -S 1000M --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf || return 1
	done

	cat <<- EOF > $genome.README || return 1
		$(date)
		$USER
		$genome
		Ensembl v$ensembl
		ftp://ftp.ensembl.org/pub/current_fasta/fasta/homo_sapiens/dna/
		ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
	EOF
	
	rm -f Homo_sapiens.GRCh38.*.chr.gtf.gz
	rm -f Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz
	rm -f chr*.gtf

	dlgenome::go $genome.fa.gtf || return 1

	return 0
}

dlgenome::hg38ensembl() {
{	dlgenome::hg38 && \
	cd $outdir && \
	echo ":INFO: downloading dbSNP" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex 'homo_sapiens-chr(MT|X|Y|[0-9]+).vcf.gz$' ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/
} || return 1
	
	echo ":INFO: extracting chrM"
	gzip -dc homo_sapiens-chrMT.vcf.gz | awk '$1~/^#/ || $NF~/^dbSNP/ {OFS="\t"; if($1!~/^#/){if($1=="MT"){$1="chrM"}else{$1="chr"$1}} print}' > $genome.fa.vcf
	for i in {1..22} X Y; do
		echo ":INFO: extracting chr$i"
		gzip -dc homo_sapiens-chr$i.vcf.gz | awk '$1!~/^#/ && $NF~/^dbSNP/ {OFS="\t"; if($1!~/^#/){if($1=="MT"){$1="chrM"}else{$1="chr"$1}} print}' >> $genome.fa.vcf
		[[ ${PIPESTATUS[0]} -eq 0 ]] || return 1
	done

cat << EOF >> $genome.README || return 1
$(date)
$(grep -m 1 -Eo 'dbSNP_[0-9]+' $genome.fa.vcf | sed 's/_/ v\./')
ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/
EOF

	rm -f homo_sapiens-chr*.vcf.gz

	return 0
}

dlgenome::hg38ensembl_v92() {
{	dlgenome::hg38 && \
	cd $outdir && \
	echo ":INFO: downloading dbSNP" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/homo_sapiens.vcf.gz && \
	echo ":INFO: extracting dbSNP" && \
	gzip -dc homo_sapiens.vcf.gz | awk '$1~/^#/ || $NF~/^dbSNP/ {OFS="\t"; if($1!~/^#/){if($1=="MT"){$1="chrM"}else{$1="chr"$1}} print}' > $genome.fa.vcf && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1

	cat <<- EOF >> $genome.README || return 1
		$(date)
		$(grep -m 1 -Eo 'dbSNP_[0-9]+' $genome.fa.vcf | sed 's/_/ v\./')
		ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/
	EOF

	rm -f homo_sapiens.vcf.gz

	return 0
}

dlgenome::hg38ncbi() {
{	dlgenome::hg38 && \
	cd $outdir && \
	echo ":INFO: downloading dbSNP" && \
	url=$(curl -s ftp://ftp.ncbi.nih.gov/snp/organisms/ | grep -oE 'human_[0-9]+_b[0-9]+_GRCh38p7' | sort -V | tail -1) && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/00-common_all.vcf.gz && \
	echo ":INFO: extracting dbSNP" && \
	gzip -dc 00-common_all.vcf.gz | awk '{OFS="\t"; if($1!~/^#/){if($1=="MT"){$1="chrM"}else{$1="chr"$1}} print}' > $genome.fa.vcf && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1

	cat <<- EOF >> $genome.README || return 1
		$(date)
		$(grep -m 1 -F dbSNP_BUILD_ID $genome.fa.vcf | sed 's/_BUILD_ID=/ \.v/')
		ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/
	EOF

	rm -f 00-common_all.vcf.gz

	return 0
}

##########
#  hg19  #
##########

dlgenome::hg19() {
{	cd $outdir && \
	genome='GRCh37.p13' && \
	echo ":INFO: downloading $genome genome" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex 'Homo_sapiens.GRCh37.dna.chromosome.(MT|X|Y|[0-9]+).fa.gz$' ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/ && \
	echo ":INFO: extracting chrM" && \
	gzip -dc Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1

	for i in {1..22} X Y; do
		echo ":INFO: extracting chr$i"
		gzip -dc Homo_sapiens.GRCh37.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
		[[ ${PIPESTATUS[0]} -eq 0 ]] || return 1
	done

{	echo ":INFO: downloading annotation" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex 'Homo_sapiens.GRCh37.[0-9]+.chr.gtf.gz$' ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/ && \
	ensembl=$(ls -v Homo_sapiens.GRCh37.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+') && \
	echo ":INFO: extracting annotation" && \
	gzip -dc Homo_sapiens.GRCh37.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};' && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1
	
	rm -f $genome.fa.gtf
	for i in M {1..22} X Y; do
		echo ":INFO: sorting chr$i"
		LC_ALL=C sort -S 1000M --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf || return 1
	done

	cat <<- EOF > $genome.README || return 1
		$(date)
		$USER
		$genome
		Ensembl v$ensembl
		ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/
		ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/
	EOF
	
	rm -f Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz
	rm -f Homo_sapiens.GRCh37.*.chr.gtf.gz
	rm -f chr*.gtf

	dlgenome::go $genome.fa.gtf || return 1

	return 0
}

dlgenome::hg19ensembl() {
{	dlgenome::hg19 && \
	cd $outdir && \
	echo ":INFO: downloading dbSNP" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/homo_sapiens.vcf.gz && \
	echo ":INFO: extracting dbSNP" && \
	gzip -dc homo_sapiens.vcf.gz | awk '$1~/^#/ || $NF~/^dbSNP/ {OFS="\t"; if($1!~/^#/){if($1=="MT"){$1="chrM"}else{$1="chr"$1}} print}' > $genome.fa.vcf && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1

	cat <<- EOF >> $genome.README || return 1
		$(date)
		$(grep -m 1 -Eo 'dbSNP_[0-9]+' $genome.fa.vcf | sed 's/_/ v\./')
		ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/
	EOF

	rm -f homo_sapiens.vcf.gz

	return 0
}

dlgenome::hg19ncbi() {
{	dlgenome::hg19 && \
	cd $outdir && \
	echo ":INFO: downloading dbSNP" && \
	url=$(curl -s ftp://ftp.ncbi.nih.gov/snp/organisms/ | grep -oE 'human_[0-9]+_b[0-9]+_GRCh37p13' | sort -V | tail -1) && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/00-common_all.vcf.gz && \
	echo ":INFO: extracting dbSNP" && \
	gzip -dc 00-common_all.vcf.gz | awk '{OFS="\t"; if($1!~/^#/){if($1=="MT"){$1="chrM"}else{$1="chr"$1}} print}' > $genome.fa.vcf && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1

	cat <<- EOF > $genome.README || return 1
		$(date)
		$(grep -m 1 -F dbSNP_BUILD_ID $genome.fa.vcf | sed -r 's/.+ID=(.+)/dbSNP v\.\1/')
		ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/
	EOF

	rm -f 00-common_all.vcf.gz

	return 0
}

##########
#  mm10  #
##########

dlgenome::mm10() {
{	cd $outdir && \
	genome=$(curl -s https://www.ensembl.org/info/website/archives/assembly.html | grep -oP 'GRCm38.p[0-9]+' | sort -V | tail -1) && \
	echo ":INFO: downloading $genome genome" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex 'Mus_musculus.GRCm38.dna.chromosome.(MT|X|Y|[0-9]+).fa.gz$' ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/ && \
	echo ":INFO: extracting chrM" && \
	gzip -dc Mus_musculus.GRCm38.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1

	for i in {1..19} X Y; do
		echo ":INFO: extracting chr$i"
		gzip -dc Mus_musculus.GRCm38.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
		[[ ${PIPESTATUS[0]} -eq 0 ]] || return 1
	done

{	echo ":INFO: downloading annotation" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N --recursive --no-directories --no-parent --level 1 --reject 'index.htm*' --accept-regex 'Mus_musculus.GRCm38.[0-9]+.chr.gtf.gz$' ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/ && \
	ensembl=$(ls -v Mus_musculus.GRCm38.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+') && \
	echo ":INFO: extracting annotation" && \
	gzip -dc Mus_musculus.GRCm38.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};' && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1
	
	rm -f $genome.fa.gtf
	for i in M {1..19} X Y; do
		echo ":INFO: sorting chr$i"
		LC_ALL=C sort -S 1000M --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf || return 1
	done

	cat <<- EOF > $genome.README || return 1
		$(date)
		$USER
		$genome
		Ensembl v$ensembl
		ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/
		ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/
	EOF
	
	rm -f Mus_musculus.GRCm38.dna.chromosome.*.fa.gz
	rm -f Mus_musculus.GRCm38.*.chr.gtf.gz
	rm -f chr*.gtf

	dlgenome::go $genome.fa.gtf || return 1

	return 0
}

dlgenome::mm10ensembl() {
{	dlgenome::mm10 && \
	cd $outdir && \
	echo ":INFO: downloading dbSNP" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N ftp://ftp.ensembl.org/pub/current_variation/vcf/mus_musculus/mus_musculus.vcf.gz && \
	echo ":INFO: extracting dbSNP" && \
	gzip -dc mus_musculus.vcf.gz | awk '$1~/^#/ || $NF~/^dbSNP/ {OFS="\t"; if($1!~/^#/){if($1=="MT"){$1="chrM"}else{$1="chr"$1}} print}' > $genome.fa.vcf && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1

	cat <<- EOF >> $genome.README || return 1
		$(date)
		$(grep -m 1 -Eo 'dbSNP_[0-9]+' $genome.fa.vcf | sed 's/_/ v\./')
		ftp://ftp.ensembl.org/pub/current_variation/vcf/mus_musculus/
	EOF
	
	rm -f mus_musculus.vcf.gz

	return 0
}

dlgenome::mm10ncbi() {
{	dlgenome::mm10 && \
	cd $outdir && \
	echo ":INFO: downloading dbSNP" && \
	wget --progress=bar:force --timeout 1 --waitretry 0 --tries 5 --retry-connrefused -N ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz && \
	echo ":INFO: extracting dbSNP" && \
	gzip -dc 00-All.vcf.gz | awk '{OFS="\t"; if($1!~/^#/){if($1=="MT"){$1="chrM"}else{$1="chr"$1}} print}' > $genome.fa.vcf && \
	[[ ${PIPESTATUS[0]} -eq 0 ]] || false
} || return 1

	cat <<- EOF > $genome.README || return 1
		$(date)
		$(grep -m 1 -F dbSNP_BUILD_ID $genome.fa.vcf | sed -r 's/.+ID=(.+)/dbSNP v\.\1/')
		ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/
	EOF
	
	rm -f 00-All.vcf.gz

	return 0
}

############# MAIN #############

checkopt (){
	local arg=false
	case $1 in
		-h | --h | -help | --help) usage;;
		-s | --s | -dbsnp | --dbsnp) release='ensembl';;
		-n | --n | -ncbi | --ncbi) release='ncbi';;
		-d | --d | -descriptions | --descriptions) go=1;;
		-v | --v | -verbose | --verbose) v=1;;

		-o | --o | -out | --out) arg=true; outdir=$2;;
		-g | --g | -genome | --genome) arg=true; g=$2;;
		-t | --t | -threads | --threads) arg=true; threads=$2;;
		-*) echo ":ERROR: illegal option $1"; return 1;; 
		*) echo ":ERROR: illegal option $2"; return 1;;
	esac
	$arg && {
		[[ ! $2 ]] && echo ":ERROR: argument missing for option $1" && return 1
		[[ "$2" =~ ^- ]] && echo ":ERROR: illegal argument $2 for option $1" && return 1
		return 0
	} || {
		[[ $2 ]] && [[ ! "$2" =~ ^- ]] && echo ":ERROR: illegal argument $2 for option $1" && return 1
		return 0
	}

	return 0
}

[[ ! $OSTYPE =~ linux ]] && die "unsupported operating system"
threads=$(cat /proc/cpuinfo | grep -cF processor)
release=''

[[ $# -eq 0 ]] && usage
[[ $# -eq 1 ]] && [[ ! $1 =~ ^- ]] && die "illegal option $1"
for i in $(seq 1 $#); do
	if [[ ${!i} =~ ^- ]]; then
		j=$((i+1))
		checkopt "${!i}" "${!j}" || exit 1
	else 
		((++i))
	fi
done

[[ ! $g ]] && echo ":ERROR: mandatory parameter -g missing" && exit 1
[[ ! $outdir ]] && echo ":ERROR: mandatory parameter -o missing" && exit 1
mkdir -p $outdir || {
	echo ":ERROR: cannot create $outdir"
	exit 1
}
outdir=$(readlink -e $outdir)
log=$outdir/dlgenome.log
touch $log || {
	echo ":ERROR: cannot create $log"
	exit 1
}

if [[ $v ]]; then
	dlgenome::$g$release 2>&1 | tee -a $log
	[[ ${PIPESTATUS[0]} -gt 0 ]] && die
else
	echo ":INFO: check log by executing: tail -f $log"
	progressbar &
	progresslog $log &
	dlgenome::$g$release &> $log
	[[ $? -gt 0 ]] && die
fi

echo ":INFO: success" | tee -a $log
exit 0
