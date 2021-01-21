#! /usr/bin/env bash
# (c) Konstantin Riege
shopt -s extglob

die(){
	echo ":ERROR: $*"
	exit 1
}

cleanup(){
	[[ $outdir ]] && rm -rf "$outdir/tmp"
}

trap 'die "killed"' INT TERM

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
	$(basename $0) downloads most recent human or mouse genome and annotation including gene ontology and optionally dbSNP

	VERSION
	0.4.0

	SYNOPSIS
	$(basename $0) -v -r [hg19|hg38|mm10] -g -a -d -s -n

	INPUT OPTIONS
	-h | --help                       : prints this message
	-t | --threads [value]            : threads - predicted default: $threads
	-o | --out [path]                 : output directory - default: $PWD
	-r | --reference [hg19|hg38|mm10] : choose GRCh37/hg19 or GRCh38/hg38 or GRCm38/mm10
	-g | --genome                     : download Ensembl genome
	-c | --ctat                       : switch to CTAT genome and indices (~30GB)
	-a | --annotation                 : download Ensembl gtf
	-d | --descriptions               : download Ensembl gene description and ontology information (requires R in PATH)
	-s | --dbsnp                      : download Ensembl dbSNP
	-n | --ncbi                       : switch to NCBI dbSNP

	REFERENCES
	(c) Konstantin Riege
	konstantin.riege{a}leibniz-fli{.}de
	EOF
	exit 1
}

######
# go #
######

dlgenome::_go(){
	out=$outdir/$genome.fa.gtf

	export R_LIBS=$outdir/tmp
	mkdir -p $R_LIBS

	[[ $(R --version | head -1 | awk '$3<3.5{print 1}') ]] && {
		cat <<- EOF > $outdir/tmp/download.R || return 1
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("biocLite", quietly = TRUE)) source("https://bioconductor.org/biocLite.R")
				biocLite("biomaRt", suppressUpdates=TRUE)
			}
			library("biomaRt")
			ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
			ensembl <- useDataset("$dataset",mart=ensembl)
			cat("downloading datasets, please wait...\n")
			goids <- getBM(mart=ensembl,attributes=c("ensembl_gene_id","go_id","namespace_1003","name_1006"))
			descriptions <- getBM(mart=ensembl,attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"))
			write.table(goids,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.go")
			write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")
		EOF
	} || {
		cat <<- EOF > $outdir/tmp/download.R || return 1
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cloud.r-project.org", Ncpus=$threads, clean=T)
				library("BiocManager")
				BiocManager::install(c("biomaRt"), Ncpus=$threads, clean=T)
			}
			library("biomaRt")
			ensembl <- useMart("ENSEMBL_MART_ENSEMBL")
			ensembl <- useDataset("$dataset",mart=ensembl)
			cat("downloading datasets, please wait...\n")
			goids <- getBM(mart=ensembl,attributes=c("ensembl_gene_id","go_id","namespace_1003","name_1006"))
			descriptions <- getBM(mart=ensembl,attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"))
			write.table(goids,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.go")
			write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")
		EOF
	}

	echo ":INFO: downloading gene ontology and descriptions"
	Rscript $outdir/tmp/download.R || return 1

	return 0
}

dlgenome::hg38.go(){
	cd $outdir
	genome=GRCh38
	dataset="hsapiens_gene_ensembl"
	dlgenome::_go
}

dlgenome::hg19.go(){
	cd $outdir
	genome=GRCh37
	dataset="hsapiens_gene_ensembl"
	dlgenome::_go
}

dlgenome::mm10.go(){
	cd $outdir
	genome=GRCm38
	dataset="mmusculus_gene_ensembl"
	dlgenome::_go
}

##########
#  hg38  #
##########

dlgenome::hg38.genome() {
	cd $outdir
	genome=GRCh38

	echo ":INFO: downloading $genome genome"
	wget -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz' || return 1

	echo ":INFO: extracting genome"
	gzip -dc Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	for i in {1..22} X Y; do
		gzip -dc Homo_sapiens.GRCh38.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
		[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	done

	ensembl=$(ls -v Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.README || return 1
		$(date)
		$USER
		Ensembl $genome genome
		ftp://ftp.ensembl.org/pub/current_fasta/fasta/homo_sapiens/dna/
	EOF

	rm -f Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz
	return 0
}

dlgenome::hg38.gtf() {
	cd $outdir
	genome=GRCh38

	echo ":INFO: downloading annotation"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.*.chr.gtf.gz' || return 1

	echo ":INFO: extracting annotation"
	gzip -dc Homo_sapiens.GRCh38.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};'
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1

	echo ":INFO: sorting annotation"
	rm -f $genome.fa.gtf
	for i in M {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf || return 1
	done

	ensembl=$(ls -v Homo_sapiens.GRCh38.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.gtf.README || return 1
		$(date)
		$USER
		Ensembl $genome gtf v$ensembl
		ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
	EOF

	rm -f Homo_sapiens.GRCh38.*.chr.gtf.gz
	rm -f chr*.gtf
	return 0
}

dlgenome::hg38.dbsnp.ensembl() {
	cd $outdir
	genome=GRCh38

	echo ":INFO: downloading dbSNP"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/homo_sapiens-chr*.vcf.gz'

	echo ":INFO: extracting dbSNP"
	gzip -dc homo_sapiens-chrMT.vcf.gz | awk -F '\t' -v OFS='\t' '$1!~/^#/ && $NF~/^dbSNP/ {if($1=="MT"){$1="chrM"}else{$1="chr"$1} print}' > chrM.vcf
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	for i in {1..22} X Y; do
		gzip -dc homo_sapiens-chr$i.vcf.gz | awk -F '\t' -v OFS='\t' '$1!~/^#/ && $NF~/^dbSNP/ {if($1=="MT"){$1="chrM"}else{$1="chr"$1} print}' > chr$i.vcf
		[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	done

	echo ":INFO: sorting dbSNP"
	gzip -dc homo_sapiens-chrMT.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	for i in M {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf || return 1
	done

	cat <<- EOF >> $genome.fa.vcf.README || return 1
		$(date)
		$USER
		Ensembl $genome $(grep -m 1 -Eo 'dbSNP_[0-9]+' $genome.fa.vcf | sed 's/_/ v\./')
		ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/
	EOF

	rm -f homo_sapiens-chr*.vcf.gz
	rm -f chr*.vcf
	return 0
}

dlgenome::hg38.dbsnp.ncbi() {
	cd $outdir
	genome=GRCh38

	echo ":INFO: downloading dbSNP"
	url=$(curl -sl ftp://ftp.ncbi.nih.gov/snp/organisms/ | grep -E 'human_[0-9]+_b[0-9]+_GRCh38' | sort -V | grep -v b151 | tail -1) # 151 gz is corrupt
	wget -c -q --show-progress  --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/00-common_all.vcf.gz

	echo ":INFO: extracting dbSNP"
	gzip -dc 00-common_all.vcf.gz | perl -F'\t' -lane 'next if /^#/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.vcf";} print O join("\t",@F); END{close O};'
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1

	echo ":INFO: sorting dbSNP"
	gzip -dc 00-common_all.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	[[ -e chrM.vcf ]] && LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chrM.vcf >> $genome.fa.vcf
	for i in {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf || return 1
	done

	cat <<- EOF >> $genome.fa.vcf.README || return 1
		$(date)
		$USER
		NCBI $genome $(grep -m 1 -F dbSNP_BUILD_ID $genome.fa.vcf | sed -E 's/.*(dbSNP)_BUILD_ID=(.+)/\1 v.\2/')
		ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/
	EOF

	rm -f 00-common_all.vcf.gz
	rm -f chr*.vcf
	return 0
}

dlgenome::hg38.ctat() {
	cd $outdir
	genome=GRCh38

	echo ":INFO: downloading CTAT"
	url='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/'
	file=$(curl -sl $url | grep -E 'GRCh38[^\"]+\.plug-n-play\.tar\.gz' | sort -Vr | head -1)
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused "$url$file"

	echo ":INFO: extracting CTAT"
	tar -xzf $file || return 1
	rm -f $file
	mv $(basename $file .tar.gz) ${genome}_CTAT_genome_lib

	cd ${genome}_CTAT_genome_lib
	cat <<- EOF >> $genome.README || return 1
		$(date)
		$USER
		CTAT $genome $(basename $file .plug-n-play.tar.gz)
		$url
	EOF

	rm -f $file

	return 0
}

##########
#  hg19  #
##########

dlgenome::hg19.genome() {
	cd $outdir
	genome=GRCh37

	echo ":INFO: downloading $genome genome"
	wget -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz' || return 1

	echo ":INFO: extracting genome"
	gzip -dc Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	for i in {1..22} X Y; do
		gzip -dc Homo_sapiens.GRCh37.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
		[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	done

	ensembl=$(ls -v Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.README || return 1
		$(date)
		$USER
		Ensembl $genome genome
		ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/
	EOF

	rm -f Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz
	return 0
}

dlgenome::hg19.gtf() {
	cd $outdir
	genome=GRCh37

	echo ":INFO: downloading annotation"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.*.chr.gtf.gz' || return 1

	echo ":INFO: extracting annotation"
	gzip -dc Homo_sapiens.GRCh37.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};'
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1

	echo ":INFO: sorting annotation"
	rm -f $genome.fa.gtf
	for i in M {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf || return 1
	done

	ensembl=$(ls -v Homo_sapiens.GRCh37.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.gtf.README || return 1
		$(date)
		$USER
		Ensembl $genome gtf v$ensembl
		ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/
	EOF

	rm -f Homo_sapiens.GRCh37.*.chr.gtf.gz
	rm -f chr*.gtf
	return 0
}

dlgenome::hg19.dbsnp.ensembl() {
	cd $outdir
	genome=GRCh37

	echo ":INFO: downloading dbSNP"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/homo_sapiens-chr*.vcf.gz'

	echo ":INFO: extracting dbSNP"
	gzip -dc homo_sapiens-chrMT.vcf.gz | awk -F '\t' -v OFS='\t' '$1!~/^#/ && $NF~/^dbSNP/ {if($1=="MT"){$1="chrM"}else{$1="chr"$1} print}' > chrM.vcf
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	for i in {1..22} X Y; do
		gzip -dc homo_sapiens-chr$i.vcf.gz | awk -F '\t' -v OFS='\t' '$1!~/^#/ && $NF~/^dbSNP/ {if($1=="MT"){$1="chrM"}else{$1="chr"$1} print}' > chr$i.vcf
		[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	done

	echo ":INFO: sorting dbSNP"
	gzip -dc homo_sapiens-chrMT.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	for i in M {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf || return 1
	done

	cat <<- EOF >> $genome.fa.vcf.README || return 1
		$(date)
		$USER
		Ensembl $genome $(grep -m 1 -Eo 'dbSNP_[0-9]+' $genome.fa.vcf | sed 's/_/ v\./')
		ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/
	EOF

	rm -f homo_sapiens-chr*.vcf.gz
	rm -f chr*.vcf
	return 0
}

dlgenome::hg19.dbsnp.ncbi() {
	cd $outdir
	genome=GRCh37

	echo ":INFO: downloading dbSNP"
	url=$(curl -sl ftp://ftp.ncbi.nih.gov/snp/organisms/ | grep -E 'human_[0-9]+_b[0-9]+_GRCh37' | sort -V | tail -1)
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/00-common_all.vcf.gz"

	echo ":INFO: extracting dbSNP"
	gzip -dc 00-common_all.vcf.gz | perl -F'\t' -lane 'next if /^#/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.vcf";} print O join("\t",@F); END{close O};'
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1

	echo ":INFO: sorting dbSNP"
	gzip -dc 00-common_all.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	[[ -e chrM.vcf ]] && LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chrM.vcf >> $genome.fa.vcf
	for i in {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf || return 1
	done

	cat <<- EOF >> $genome.fa.vcf.README || return 1
		$(date)
		$USER
		NCBI $genome $(grep -m 1 -F dbSNP_BUILD_ID $genome.fa.vcf | sed -E 's/.*(dbSNP)_BUILD_ID=(.+)/\1 v.\2/')
		ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/
	EOF

	rm -f 00-common_all.vcf.gz
	rm -f chr*.vcf
	return 0
}

dlgenome::hg19.ctat() {
	cd $outdir
	genome=GRCh37

	echo ":INFO: downloading CTAT"
	url='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/'
	file=$(curl -sl $url | grep -E 'GRCh37[^\"]+\.plug-n-play\.tar\.gz' | sort -Vr | head -1)
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused "$url$file"

	echo ":INFO: extracting CTAT"
	tar -xzf $file || return 1
	rm -f $file
	mv $(basename $file .tar.gz) ${genome}_CTAT_genome_lib

	cd ${genome}_CTAT_genome_lib
	cat <<- EOF >> $genome.README || return 1
		$(date)
		$USER
		CTAT $genome $(basename $file .plug-n-play.tar.gz)
		$url
	EOF

	rm -f $file

	return 0
}

##########
#  mm10  #
##########

dlgenome::mm10.genome() {
	cd $outdir
	genome=GRCm38

	echo ":INFO: downloading $genome genome"
	wget -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.*.fa.gz' || return 1

	echo ":INFO: extracting genome"
	gzip -dc Mus_musculus.GRCm38.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	for i in {1..19} X Y; do
		gzip -dc Mus_musculus.GRCm38.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
		[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1
	done

	ensembl=$(ls -v Mus_musculus.GRCm38.dna.chromosome.*.fa.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.README || return 1
		$(date)
		$USER
		Ensembl $genome genome
		ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/
	EOF

	rm -f Mus_musculus.GRCm38.dna.chromosome.*.fa.gz
	return 0
}

dlgenome::mm10.gtf() {
	cd $outdir
	genome=GRCm38

	echo ":INFO: downloading annotation"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm38.*.chr.gtf.gz' || return 1

	echo ":INFO: extracting annotation"
	gzip -dc Mus_musculus.GRCm38.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};'
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1

	echo ":INFO: sorting annotation"
	rm -f $genome.fa.gtf
	for i in M {1..19} X Y; do
		LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf || return 1
	done

	ensembl=$(ls -v Mus_musculus.GRCm38.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.gtf.README || return 1
		$(date)
		$USER
		Ensembl $genome gtf v$ensembl
		ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/
	EOF

	rm -f Mus_musculus.GRCm38.*.chr.gtf.gz
	rm -f chr*.gtf
	return 0
}

dlgenome::mm10.dbsnp.ensembl() {
	cd $outdir
	genome=GRCm38

	echo ":INFO: downloading dbSNP"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_variation/vcf/mus_musculus/mus_musculus.vcf.gz'

	echo ":INFO: extracting dbSNP"
	gzip -dc mus_musculus.vcf.gz | perl -F'\t' -lane 'next if /^#/ || $F[-1]!~/^dbSNP/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.vcf";} print O join("\t",@F); END{close O};'
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1

	echo ":INFO: sorting dbSNP"
	gzip -dc mus_musculus.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	[[ -e chrM.vcf ]] && LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chrM.vcf >> $genome.fa.vcf
	for i in {1..19} X Y; do
		LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf || return 1
	done

	cat <<- EOF >> $genome.fa.vcf.README || return 1
		$(date)
		$USER
		Ensembl $genome $(grep -m 1 -Eo 'dbSNP_[0-9]+' $genome.fa.vcf | sed 's/_/ v\./')
		ftp://ftp.ensembl.org/pub/current_variation/vcf/mus_musculus/
	EOF

	rm -f mus_musculus.vcf.gz
	rm -f chr*.vcf
	return 0
}

dlgenome::mm10.dbsnp.ncbi() {
	cd $outdir
	genome=GRCm38

	echo ":INFO: downloading dbSNP"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping 'ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz'

	echo ":INFO: extracting dbSNP"
	gzip -dc 00-All.vcf.gz | perl -F'\t' -lane 'next if /^#/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.vcf";} print O join("\t",@F); END{close O};'
	[[ $((${PIPESTATUS[@]/%/+}0)) -gt 0 ]] && return 1

	echo ":INFO: sorting dbSNP"
	gzip -dc 00-All.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	[[ -e chrM.vcf ]] && LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chrM.vcf >> $genome.fa.vcf
	for i in {1..19} X Y; do
		LC_ALL=C sort -S 1000M -T /dev/shm --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf || return 1
	done

	cat <<- EOF >> $genome.fa.vcf.README || return 1
		$(date)
		$USER
		NCBI $genome $(grep -m 1 -F dbSNP_BUILD_ID $genome.fa.vcf | sed -E 's/.*(dbSNP)_BUILD_ID=(.+)/\1 v.\2/')
		ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/
	EOF

	rm -f 00-All.vcf.gz
	rm -f chr*.vcf
	return 0
}

dlgenome::mm10.ctat() {
	cd $outdir
	genome=GRCm38

	echo ":INFO: downloading CTAT"
	url='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/'
	file=$(curl -sl $url | grep -E 'Mouse[^\"]+\.plug-n-play\.tar\.gz' | sort -Vr | head -1)
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused "$url$file"

	echo ":INFO: extracting CTAT"
	tar -xzf $file || return 1
	rm -f $file
	mv $(basename $file .tar.gz) ${genome}_CTAT_genome_lib

	cd ${genome}_CTAT_genome_lib
	cat <<- EOF >> $genome.README || return 1
		$(date)
		$USER
		CTAT $genome $(basename $file .plug-n-play.tar.gz)
		$url
	EOF

	rm -f $file

	return 0
}

############# MAIN #############

checkopt() {
	local arg=false
	case $1 in
		-h | --h | -help | --help) (usage); exit 0;;
		-t | --t | -threads | --threads) arg=true; threads=$2;;
		-o | --o | -out | --out) arg=true; outdir=$2;;
		-r | --r | -reference | --reference) arg=true; ref=$2;;
		-g | --g | -genome | --genome) fun+=("genome");;
		-c | --c | -ctat | --ctat) fun+=("ctat");;
		-a | --a | -annotation | --annotation) fun+=("gtf");;
		-s | --s | -dbsnp | --dbsnp) db="ensembl";;
		-n | --n | -ncbi | --ncbi) db='ncbi';;
		-d | --d | -descriptions | --descriptions) fun+=("go");;
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
fun=()
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
[[ $db ]] && fun+=("dbsnp.$db")
[[ ! $ref ]] && echo ":ERROR: mandatory parameter -r missing" && exit 1
outdir=${outdir:-$PWD}
mkdir -p $outdir || die ":ERROR: cannot create $outdir"
outdir=$(readlink -e $outdir)
log=$outdir/dlgenome.log
rm -f $log
touch $log || die ":ERROR: cannot create $log"


for f in "${fun[@]}"; do
	dlgenome::$ref.$f 2>&1 | tee -ai $log
	[[ ${PIPESTATUS[0]} -gt 0 ]] && die "fetching $f failed"
done
echo ":INFO: success" | tee -ai $log
exit 0
