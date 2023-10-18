#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$0")/../bashbone_lite.sh" -x cleanup -a "$@" || exit 1

############################################

cleanup(){
	[[ $outdir ]] && rm -rf "$outdir/tmp"
	return 0
}

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename $0) downloads most recent human or mouse genome and annotation including gene ontology and optionally dbSNP

		VERSION
		0.6.0

		SYNOPSIS
		$(basename $0) -v -r [hg19|hg38|mm10] -g -a -d -s -n

		INPUT OPTIONS
		-h | --help               : prints this message
		-t | --threads [value]    : threads - predicted default: $threads
		-o | --out [path]         : output directory - default: $PWD
		-r | --reference [string] : choose hg19 hg38 mm10 mm11
		-g | --genome             : download Ensembl genome
		-c | --ctat               : switch to CTAT genome and indices (~30GB)
		-a | --annotation         : download Ensembl gtf
		-d | --descriptions       : download Ensembl gene description and Ensembl gene ontology information (requires R in PATH)
		-m | --msigdb             : download Ensembl gene description and MSigDB gene ontology information (requires R in PATH)
		-e | --enrichr            : download Ensembl gene description and Enrichr human gene ontology information (requires R in PATH)
		-s | --dbsnp              : download Ensembl dbSNP
		-n | --ncbi               : switch to NCBI dbSNP

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF
	return 1
}

######
# go #
######

dlgenome::_go.ensembl(){
	out="$outdir/$genome.fa.gtf"

	export R_LIBS="$outdir/tmp"
	mkdir -p "$R_LIBS"
	# listDatasets(useEnsembl(biomart = "genes"))
	# listDatasets(useEnsembl(biomart = "ensembl"))

	[[ $(R --version | head -1 | awk '$3<3.5{print 1}') ]] && {
		cat <<- EOF > $outdir/tmp/download.R
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("biocLite", quietly = TRUE)) source("https://bioconductor.org/biocLite.R")
				biocLite("biomaRt", suppressUpdates=TRUE)
			}
			library("biomaRt")
			v <- "$version"
			if ("$version" == "latest") v <- listEnsemblArchives()\$version[2]
			ensembl <- useEnsembl(biomart="genes", dataset="$dataset", version=v)
			goids <- data.frame()
			descriptions <- data.frame()
			for (chr in grep("^(MT|X|Y|\\\d+)$",listFilterOptions(mart = ensembl, filter = "chromosome_name"), value=T, perl=T)){
				cat(paste0("downloading datasets of chr",chr,", please wait...\n"))
				goids <- rbind(goids, getBM(mart=ensembl, attributes=c("ensembl_gene_id","go_id","namespace_1003","name_1006"), filters=c("chromosome_name"), values=list(chr)))
				descriptions <- rbind(descriptions, getBM(mart=ensembl, attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters=c("chromosome_name"), values=list(chr)))
			}
			descriptions[,2][descriptions[2]==""] <- descriptions[,1][descriptions[2]==""]
			write.table(goids,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.go")
			write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")
			sink("$out.go.README")
			cat("$(date)\n")
			cat("$USER\n")
			cat(paste0("Ensembl go v",v,"\n")))
			cat("via biomart\n")
			sink()
		EOF
	} || {
		cat <<- EOF > $outdir/tmp/download.R
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cloud.r-project.org", Ncpus=$threads, clean=T)
				library("BiocManager")
				BiocManager::install(c("biomaRt"), Ncpus=$threads, clean=T)
			}
			library("biomaRt")
			v <- "$version"
			if ("$version" == "latest") v <- listEnsemblArchives()\$version[2]
			ensembl <- useEnsembl(biomart="genes", dataset="$dataset", version=v)
			goids <- data.frame()
			descriptions <- data.frame()
			for (chr in grep("^(MT|X|Y|\\\d+)$",listFilterOptions(mart = ensembl, filter = "chromosome_name"), value=T, perl=T)){
				cat(paste0("downloading datasets of chr",chr,", please wait...\n"))
				goids <- rbind(goids, getBM(mart=ensembl, attributes=c("ensembl_gene_id","go_id","namespace_1003","name_1006"), filters=c("chromosome_name"), values=list(chr)))
				descriptions <- rbind(descriptions, getBM(mart=ensembl, attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters=c("chromosome_name"), values=list(chr)))
			}
			descriptions[,2][descriptions[2]==""] <- descriptions[,1][descriptions[2]==""]
			write.table(goids,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.go")
			write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")
			sink("$out.go.README")
			cat("$(date)\n")
			cat("$USER\n")
			cat(paste0("Ensembl go v",v,"\n"))
			cat("via biomart\n")
			sink()
		EOF
	}

	echo ":INFO: downloading gene ontology and descriptions"
	Rscript "$outdir/tmp/download.R"
	return 0
}

dlgenome::_go.enrichr(){
	out="$outdir/$genome.fa.gtf"

	export R_LIBS="$outdir/tmp"
	mkdir -p "$R_LIBS"

	[[ $(R --version | head -1 | awk '$3<3.5{print 1}') ]] && {
		cat <<- EOF > $outdir/tmp/download.R
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("biocLite", quietly = TRUE)) source("https://bioconductor.org/biocLite.R")
				biocLite("biomaRt", suppressUpdates=TRUE)
			}
			library("biomaRt")
			v <- "$version"
			if ("$version" == "latest") v <- listEnsemblArchives()\$version[2]
			ensembl <- useEnsembl(biomart="genes", dataset="$dataset", version=v)
			descriptions <- data.frame()
			for (chr in grep("^(MT|X|Y|\\\d+)$",listFilterOptions(mart = ensembl, filter = "chromosome_name"), value=T, perl=T)){
				cat(paste0("downloading datasets of chr",chr,", please wait...\n"))
				descriptions <- rbind(descriptions, getBM(mart=ensembl, attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters=c("chromosome_name"), values=list(chr)))
			}
			descriptions[,2][descriptions[2]==""] <- descriptions[,1][descriptions[2]==""]
			write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")
		EOF
	} || {
		cat <<- EOF > $outdir/tmp/download.R
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cloud.r-project.org", Ncpus=$threads, clean=T)
				library("BiocManager")
				BiocManager::install(c("biomaRt"), Ncpus=$threads, clean=T)
			}
			library("biomaRt")
			v <- "$version"
			if ("$version" == "latest") v <- listEnsemblArchives()\$version[2]
			ensembl <- useEnsembl(biomart="genes", dataset="$dataset", version=v)
			descriptions <- data.frame()
			for (chr in grep("^(MT|X|Y|\\\d+)$",listFilterOptions(mart = ensembl, filter = "chromosome_name"), value=T, perl=T)){
				cat(paste0("downloading datasets of chr",chr,", please wait...\n"))
				descriptions <- rbind(descriptions, getBM(mart=ensembl, attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters=c("chromosome_name"), values=list(chr)))
			}
			descriptions[,2][descriptions[2]==""] <- descriptions[,1][descriptions[2]==""]
			write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")
		EOF
	}

	echo ":INFO: downloading gene ontology and descriptions"
	Rscript "$outdir/tmp/download.R"

	wget -O "$outdir/tmp/ncbi.info.gz" -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
	gzip -dc "$outdir/tmp/ncbi.info.gz" | perl -lanE 'next unless $F[5]=~/Ensembl:(ENSG\d+)/; say "$F[2]\t$1"' > "$outdir/tmp/ncbi.info"

	rm -f "$out.go"
	for domain in Biological_Process Molecular_Function Cellular_Component; do
		version=$(date +%Y)
		while :; do	wget -O "$outdir/tmp/$domain.gmt" -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_${domain}_$version" && break || ((version--)); done
		perl -F'\t' -slanE '
			BEGIN{
				open F,"<$info" or die $!;
				while(<F>){
					chomp;
					@F=split/\t/;
					$g2i{$F[0]}=$F[1];
				}
				close F;
			}
			$F[0]=~s/\s+\((GO:\d+)\)$//;
			for (@F[2..$#F]){
				$i=$g2i{$_};
				next unless $i;
				say join"\t",($i,$1,$domain,$F[0]);
			}
		' -- -info="$outdir/tmp/ncbi.info" -domain="${domain,,}" "$outdir/tmp/$domain.gmt" >> "$out.go"
	done

	cat <<-EOF > "$out.go.README"
		$(date)
		$USER
		Enrichr go v$version
	EOF
	return 0
}

dlgenome::_go.msigdb(){
	out="$outdir/$genome.fa.gtf"

	export R_LIBS="$outdir/tmp"
	mkdir -p "$R_LIBS"

	[[ $(R --version | head -1 | awk '$3<3.5{print 1}') ]] && {
		cat <<- EOF > $outdir/tmp/download.R
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("biocLite", quietly = TRUE)) source("https://bioconductor.org/biocLite.R")
				biocLite("biomaRt", suppressUpdates=TRUE)
			}
			library("biomaRt")
			v <- "$version"
			if ("$version" == "latest") v <- listEnsemblArchives()\$version[2]
			ensembl <- useEnsembl(biomart="genes", dataset="$dataset", version=v)
			descriptions <- data.frame()
			for (chr in grep("^(MT|X|Y|\\\d+)$",listFilterOptions(mart = ensembl, filter = "chromosome_name"), value=T, perl=T)){
				cat(paste0("downloading datasets of chr",chr,", please wait...\n"))
				descriptions <- rbind(descriptions, getBM(mart=ensembl, attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters=c("chromosome_name"), values=list(chr)))
			}
			descriptions[,2][descriptions[2]==""] <- descriptions[,1][descriptions[2]==""]
			write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")
		EOF
	} || {
		cat <<- EOF > $outdir/tmp/download.R
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cloud.r-project.org", Ncpus=$threads, clean=T)
				library("BiocManager")
				BiocManager::install(c("biomaRt"), Ncpus=$threads, clean=T)
			}
			library("biomaRt")
			v <- "$version"
			if ("$version" == "latest") v <- listEnsemblArchives()\$version[2]
			ensembl <- useEnsembl(biomart="genes", dataset="$dataset", version=v)
			descriptions <- data.frame()
			for (chr in grep("^(MT|X|Y|\\\d+)$",listFilterOptions(mart = ensembl, filter = "chromosome_name"), value=T, perl=T)){
				cat(paste0("downloading datasets of chr",chr,", please wait...\n"))
				descriptions <- rbind(descriptions, getBM(mart=ensembl, attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters=c("chromosome_name"), values=list(chr)))
			}
			descriptions[,2][descriptions[2]==""] <- descriptions[,1][descriptions[2]==""]
			write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")
		EOF
	}

	echo ":INFO: downloading gene ontology and descriptions"
	Rscript "$outdir/tmp/download.R"

	version=$(curl -s https://data.broadinstitute.org/gsea-msigdb/msigdb/release/ | grep -oE ">[0-9][^<]+($msig)" | sed 's/^>//' | sort -Vr | head -1)
	[[ "$msig" == "Hs" ]] && msig="human" || msig="mouse"
	wget -P "$outdir/tmp/" -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --recursive --no-directories --no-parent --level=1 --reject 'index.htm*' --accept-regex ".*\.go\.(bp|mf|cc)\.v$version.symbols.gmt" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/$version/"
	wget -O "$outdir/tmp/go.obo" -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping http://purl.obolibrary.org/obo/go.obo

	grep -Fx '[Term]' -A 3 "$outdir/tmp/go.obo" | grep -vFx -- '--' | paste - - - - | cut -f 2- | sed -E 's/(id:|name:|namespace:)\s//g' | perl -F'\t' -lane '$F[1]=~s/#//; $F[1]=~s/^OBSOLETE\s*//i; print join"\t",@F' > "$outdir/tmp/go.terms.orig"
	grep -Fx '[Term]' -A 3 "$outdir/tmp/go.obo" | grep -vFx -- '--' | paste - - - - | cut -f 2- | sed -E 's/(id:|name:|namespace:)\s//g' | perl -F'\t' -lane '$F[1]=uc($F[1]=~s/\W+/_/gr); $F[1]=~s/^_+//; $F[1]=~s/_+$//; $F[1]=~s/_+/_/g; $F[1]=~s/^OBSOLETE_//; print join"\t",@F' > "$outdir/tmp/go.terms"

	rm -f "$out.go"
	for gmt in "$outdir/tmp/"*.v$version.symbols.gmt; do
		cut -f 1 "$gmt" | sed -E 's/^GO.._//' > "$outdir/tmp/gmt.terms"
		grep -Fw -f "$outdir/tmp/gmt.terms" "$outdir/tmp/go.terms" > "$outdir/tmp/gmt.terms.parsed"
		grep -vFx -f <(cut -f 2 "$outdir/tmp/gmt.terms.parsed") "$outdir/tmp/gmt.terms" > "$outdir/tmp/gmt.terms.missing"
		c=$(cut -d . -f 3 <<< "$gmt" | perl -lane 'print "GO".uc($_)')
		x=$(wc -l < "$outdir/tmp/gmt.terms.missing")
		i=0
		while read -r n; do
			echo "requesting $((++i)) of $x missing $c terms" >&2
			g=$(curl -s "https://www.gsea-msigdb.org/gsea/msigdb/$msig/geneset/${c}_$n" | grep -oP 'GO:\d+</td>' | cut -d '<' -f 1 | sort -u)
			[[ $g ]] || echo ":WARNING: no term found for $n"
			grep -m 1 -Fw $g "$outdir/tmp/go.terms" | awk -F '\t' -v n=$n -v OFS='\t' '{print $1,n,$3}' >> "$outdir/tmp/gmt.terms.parsed"
		done < "$outdir/tmp/gmt.terms.missing"
		perl -F'\t' -slanE '
			BEGIN{
				open F,"<$terms" or die $!;
				while(<F>){
					chomp;
					@F=split/\t/;
					$t2g{$F[1]}=$F[0];
				}
				close F;
				open F,"<$go" or die $!;
				while(<F>){
					chomp;
					@F=split/\t/;
					$g2x{$F[0]}=join("\t",($F[0],$F[2],$F[1]));
				}
				close F;
				open F,"<$info" or die $!;
				while(<F>){
					chomp;
					@F=split/\t/;
					$F[1]=$F[0] unless $F[1];
					$g2i{$F[1]}=$F[0];
				}
				close F;
			}
			$x=$g2x{$t2g{$F[0]=~s/^GO.._//r}};
			next unless $x;
			for (@F[2..$#F]){
				$i=$g2i{$_};
				next unless $i;
				say $i."\t".$x;
			}
		' -- -terms="$outdir/tmp/gmt.terms.parsed" -go="$outdir/tmp/go.terms.orig" -info="$out.info" "$gmt" >> "$out.go"
	done

	cat <<-EOF > "$out.go.README"
		$(date)
		$USER
		MSigDB go v$version
	EOF
	return 0
}

dlgenome::hg38.go(){
	cd "$outdir"
	genome=GRCh38
	dataset="hsapiens_gene_ensembl"
	version="latest"
	msig=Hs
	dlgenome::_go.$godb
	return 0
}

dlgenome::hg19.go(){
	cd "$outdir"
	genome=GRCh37
	dataset="hsapiens_gene_ensembl"
	version="GRCh37"
	msig=Hs
	dlgenome::_go.$godb
	return 0
}

dlgenome::mm10.go(){
	cd "$outdir"
	genome=GRCm38
	dataset="mmusculus_gene_ensembl"
	version="102"
	msig=Mm
	dlgenome::_go.$godb
	return 0
}

dlgenome::mm11.go(){
	cd "$outdir"
	genome=GRCm39
	dataset="mmusculus_gene_ensembl"
	version="latest"
	msig=Mm
	dlgenome::_go.$godb
	return 0
}

##########
#  hg38  #
##########

dlgenome::hg38.genome() {
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading $genome genome"
	wget -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz'

	echo ":INFO: extracting genome"
	gzip -dc Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa
	for i in {1..22} X Y; do
		gzip -dc Homo_sapiens.GRCh38.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
	done

	ensembl=$(ls -v Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.README
		$(date)
		$USER
		Ensembl $genome genome
		ftp://ftp.ensembl.org/pub/current_fasta/fasta/homo_sapiens/dna/
	EOF

	rm -f Homo_sapiens.GRCh38.dna.chromosome.*.fa.gz
	return 0
}

dlgenome::hg38.gtf() {
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading annotation"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/Homo_sapiens.GRCh38.*.chr.gtf.gz'

	echo ":INFO: extracting annotation"
	gzip -dc Homo_sapiens.GRCh38.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};'

	echo ":INFO: sorting annotation"
	rm -f $genome.fa.gtf
	for i in M {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf
	done

	ensembl=$(ls -v Homo_sapiens.GRCh38.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.gtf.README
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
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading dbSNP"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/homo_sapiens-chr*.vcf.gz'

	echo ":INFO: extracting dbSNP"
	gzip -dc homo_sapiens-chrMT.vcf.gz | awk -F '\t' -v OFS='\t' '$1!~/^#/ && $NF~/^dbSNP/ {if($1=="MT"){$1="chrM"}else{$1="chr"$1} print}' > chrM.vcf
	for i in {1..22} X Y; do
		gzip -dc homo_sapiens-chr$i.vcf.gz | awk -F '\t' -v OFS='\t' '$1!~/^#/ && $NF~/^dbSNP/ {if($1=="MT"){$1="chrM"}else{$1="chr"$1} print}' > chr$i.vcf
	done

	echo ":INFO: sorting dbSNP"
	gzip -dc homo_sapiens-chrMT.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	for i in M {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf
	done

	cat <<- EOF >> $genome.fa.vcf.README
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
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading dbSNP"
	url=$(curl -sl ftp://ftp.ncbi.nih.gov/snp/organisms/ | grep -E 'human_[0-9]+_b[0-9]+_GRCh38' | sort -V | tail -1) # 151 gz is corrupt

	wget -c -q --show-progress  --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/00-common_all.vcf.gz

	echo ":INFO: extracting dbSNP"
	gzip -dc 00-common_all.vcf.gz | perl -F'\t' -lane 'next if /^#/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.vcf";} print O join("\t",@F); END{close O};'

	echo ":INFO: sorting dbSNP"
	gzip -dc 00-common_all.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	[[ -e chrM.vcf ]] && LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chrM.vcf >> $genome.fa.vcf
	for i in {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf
	done

	cat <<- EOF >> $genome.fa.vcf.README
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
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading CTAT"
	url='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/'
	file=$(curl -sl $url | grep -E 'GRCh38[^\"]+\.plug-n-play\.tar\.gz' | sort -Vr | head -1)
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused "$url$file"

	echo ":INFO: extracting CTAT"
	tar -xzf $file
	rm -f $file
	mv $(basename $file .tar.gz) ${genome}_CTAT_genome_lib

	cd ${genome}_CTAT_genome_lib
	cat <<- EOF >> $genome.README
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
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading $genome genome"
	wget -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz'

	echo ":INFO: extracting genome"
	gzip -dc Homo_sapiens.GRCh37.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa
	for i in {1..22} X Y; do
		gzip -dc Homo_sapiens.GRCh37.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
	done

	ensembl=$(ls -v Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.README
		$(date)
		$USER
		Ensembl $genome genome
		ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/
	EOF

	rm -f Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz
	return 0
}

dlgenome::hg19.gtf() {
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading annotation"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/Homo_sapiens.GRCh37.*.chr.gtf.gz'

	echo ":INFO: extracting annotation"
	gzip -dc Homo_sapiens.GRCh37.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};'

	echo ":INFO: sorting annotation"
	rm -f $genome.fa.gtf
	for i in M {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf
	done

	ensembl=$(ls -v Homo_sapiens.GRCh37.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.gtf.README
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
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading dbSNP"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/grch37/current/variation/vcf/homo_sapiens/homo_sapiens-chr*.vcf.gz'

	echo ":INFO: extracting dbSNP"
	gzip -dc homo_sapiens-chrMT.vcf.gz | awk -F '\t' -v OFS='\t' '$1!~/^#/ && $NF~/^dbSNP/ {if($1=="MT"){$1="chrM"}else{$1="chr"$1} print}' > chrM.vcf
	for i in {1..22} X Y; do
		gzip -dc homo_sapiens-chr$i.vcf.gz | awk -F '\t' -v OFS='\t' '$1!~/^#/ && $NF~/^dbSNP/ {if($1=="MT"){$1="chrM"}else{$1="chr"$1} print}' > chr$i.vcf
	done

	echo ":INFO: sorting dbSNP"
	gzip -dc homo_sapiens-chrMT.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	for i in M {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf
	done

	cat <<- EOF >> $genome.fa.vcf.README
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
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading dbSNP"
	url=$(curl -sl ftp://ftp.ncbi.nih.gov/snp/organisms/ | grep -E 'human_[0-9]+_b[0-9]+_GRCh37' | sort -V | tail -1)
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "ftp://ftp.ncbi.nih.gov/snp/organisms/$url/VCF/00-common_all.vcf.gz"

	echo ":INFO: extracting dbSNP"
	gzip -dc 00-common_all.vcf.gz | perl -F'\t' -lane 'next if /^#/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.vcf";} print O join("\t",@F); END{close O};'

	echo ":INFO: sorting dbSNP"
	gzip -dc 00-common_all.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	[[ -e chrM.vcf ]] && LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chrM.vcf >> $genome.fa.vcf
	for i in {1..22} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf
	done

	cat <<- EOF >> $genome.fa.vcf.README
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
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading CTAT"
	url='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/'
	file=$(curl -sl $url | grep -E 'GRCh37[^\"]+\.plug-n-play\.tar\.gz' | sort -Vr | head -1)
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused "$url$file"

	echo ":INFO: extracting CTAT"
	tar -xzf $file
	rm -f $file
	mv $(basename $file .tar.gz) ${genome}_CTAT_genome_lib

	cd ${genome}_CTAT_genome_lib
	cat <<- EOF >> $genome.README
		$(date)
		$USER
		CTAT $genome $(basename $file .plug-n-play.tar.gz)
		$url
	EOF

	rm -f $file
	return 0
}

##########
#  mm11  #
##########

dlgenome::mm11.genome() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: downloading $genome genome"
	wget -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.*.fa.gz'

	echo ":INFO: extracting genome"
	gzip -dc Mus_musculus.GRCm39.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa
	for i in {1..19} X Y; do
		gzip -dc Mus_musculus.GRCm39.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
	done

	ensembl=$(ls -v Mus_musculus.GRCm39.dna.chromosome.*.fa.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.README
		$(date)
		$USER
		Ensembl $genome genome
		ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/
	EOF

	rm -f Mus_musculus.GRCm39.dna.chromosome.*.fa.gz
	return 0
}

dlgenome::mm11.gtf() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: downloading annotation"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_gtf/mus_musculus/Mus_musculus.GRCm39.*.chr.gtf.gz'

	echo ":INFO: extracting annotation"
	gzip -dc Mus_musculus.GRCm39.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};'

	echo ":INFO: sorting annotation"
	rm -f $genome.fa.gtf
	for i in M {1..19} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf
	done

	ensembl=$(ls -v Mus_musculus.GRCm39.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.gtf.README
		$(date)
		$USER
		Ensembl $genome gtf v$ensembl
		ftp://ftp.ensembl.org/pub/grch37/current/gtf/homo_sapiens/
	EOF

	rm -f Mus_musculus.GRCm39.*.chr.gtf.gz
	rm -f chr*.gtf
	return 0
}

dlgenome::mm11.dbsnp.ensembl() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: downloading dbSNP"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_variation/vcf/mus_musculus/mus_musculus.vcf.gz'

	echo ":INFO: extracting dbSNP"
	# as of v109 no dbsnp annotation in vcf
	gzip -dc mus_musculus.vcf.gz | perl -F'\t' -lane 'next if /^#/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.vcf";} print O join("\t",@F); END{close O};'

	echo ":INFO: sorting dbSNP"
	gzip -dc mus_musculus.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	[[ -e chrM.vcf ]] && LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chrM.vcf >> $genome.fa.vcf
	for i in {1..19} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf
	done

	cat <<- EOF >> $genome.fa.vcf.README
		$(date)
		$USER
		Ensembl $genome $(grep -m 1 -Eo 'dbSNP_[0-9]+' $genome.fa.vcf | sed 's/_/ v\./')
		ftp://ftp.ensembl.org/pub/current_variation/vcf/mus_musculus/
	EOF

	rm -f mus_musculus.vcf.gz
	rm -f chr*.vcf
	return 0
}

dlgenome::mm11.dbsnp.ncbi() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: ncbi dbsnp not (yet?) available"
	return 0
}

dlgenome::mm11.ctat() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: downloading CTAT"
	url='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/'
	file=$(curl -sl $url | grep -E 'Mouse_GRCm39[^\"]+\.plug-n-play\.tar\.gz' | sort -Vr | head -1)
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused "$url$file"

	echo ":INFO: extracting CTAT"
	tar -xzf $file
	rm -f $file
	mv $(basename $file .tar.gz) ${genome}_CTAT_genome_lib

	cd ${genome}_CTAT_genome_lib
	cat <<- EOF >> $genome.README
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
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: downloading $genome genome"
	wget -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.*.fa.gz'

	echo ":INFO: extracting genome"
	gzip -dc Mus_musculus.GRCm38.dna.chromosome.MT.fa.gz | sed "s/>.*/>chrM/" > $genome.fa
	for i in {1..19} X Y; do
		gzip -dc Mus_musculus.GRCm38.dna.chromosome.$i.fa.gz | sed "s/>.*/>chr$i/" >> $genome.fa
	done

	ensembl=$(ls -v Mus_musculus.GRCm38.dna.chromosome.*.fa.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.README
		$(date)
		$USER
		Ensembl $genome genome
		ftp://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/
	EOF

	rm -f Mus_musculus.GRCm38.dna.chromosome.*.fa.gz
	return 0
}

dlgenome::mm10.gtf() {
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: downloading annotation"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.*.chr.gtf.gz'

	echo ":INFO: extracting annotation"
	gzip -dc Mus_musculus.GRCm38.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};'

	echo ":INFO: sorting annotation"
	rm -f $genome.fa.gtf
	for i in M {1..19} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf
	done

	ensembl=$(ls -v Mus_musculus.GRCm38.*.chr.gtf.gz | tail -1 | grep -Eo '\.[0-9]+')
	cat <<- EOF > $genome.fa.gtf.README
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
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: downloading dbSNP"
	# wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/current_variation/vcf/mus_musculus/mus_musculus.vcf.gz'
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/release-102/variation/vcf/mus_musculus/mus_musculus.vcf.gz'

	echo ":INFO: extracting dbSNP"
	gzip -dc mus_musculus.vcf.gz | perl -F'\t' -lane 'next if /^#/ || $F[-1]!~/^dbSNP/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.vcf";} print O join("\t",@F); END{close O};'

	echo ":INFO: sorting dbSNP"
	gzip -dc mus_musculus.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	[[ -e chrM.vcf ]] && LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chrM.vcf >> $genome.fa.vcf
	for i in {1..19} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf
	done

	cat <<- EOF >> $genome.fa.vcf.README
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
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: downloading dbSNP"
	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping 'ftp://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz'

	echo ":INFO: extracting dbSNP"
	gzip -dc 00-All.vcf.gz | perl -F'\t' -lane 'next if /^#/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.vcf";} print O join("\t",@F); END{close O};'

	echo ":INFO: sorting dbSNP"
	gzip -dc 00-All.vcf.gz | head -1000 | grep '^#' > $genome.fa.vcf
	[[ -e chrM.vcf ]] && LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chrM.vcf >> $genome.fa.vcf
	for i in {1..19} X Y; do
		LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k2,2n chr$i.vcf >> $genome.fa.vcf
	done

	cat <<- EOF >> $genome.fa.vcf.README
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
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: ctat not available anymore"
	return 0
}

############# MAIN #############

checkopt() {
	local arg=false
	case $1 in
		-h | --h | -help | --help) { usage || exit 0; };;
		-t | --t | -threads | --threads) arg=true; threads=$2;;
		-o | --o | -out | --out) arg=true; outdir=$2;;
		-r | --r | -reference | --reference) arg=true; ref=$2;;
		-g | --g | -genome | --genome) fun+=("genome");;
		-c | --c | -ctat | --ctat) fun+=("ctat");;
		-a | --a | -annotation | --annotation) fun+=("gtf");;
		-s | --s | -dbsnp | --dbsnp) db="ensembl";;
		-n | --n | -ncbi | --ncbi) db='ncbi';;
		-d | --d | -descriptions | --descriptions) godb="ensembl";;
		-m | --m | -msigdb | --msigdb) godb="msigdb";;
		-e | --e | -enrichr | --enrichr) godb="enrichr";;
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

threads=$(grep -cF processor /proc/cpuinfo)
fun=()
[[ $# -eq 0 ]] && { usage || exit 0; }
[[ $# -eq 1 ]] && [[ ! $1 =~ ^- ]] && false
for i in $(seq 1 $#); do
	if [[ ${!i} =~ ^- ]]; then
		j=$((i+1))
		checkopt "${!i}" "${!j}" || exit 1
	else
		((++i))
	fi
done
[[ $db ]] && fun+=("dbsnp.$db")
[[ $godb ]] && fun+=("go")
BASHBONE_ERROR="mandatory parameter -r missing"
[[ $ref ]]
outdir="${outdir:-$PWD}"
BASHBONE_ERROR="cannot create $outdir"
mkdir -p "$outdir"
outdir="$(readlink -e "$outdir")"
log="$outdir/dlgenome.log"
rm -f "$log"
BASHBONE_ERROR="cannot create $log"
touch "$log"

for f in "${fun[@]}"; do
	BASHBONE_ERROR="dlgenome::$ref.$f failed"
	dlgenome::$ref.$f 2>&1 | tee -ai "$log"
done
echo ":INFO: success" | tee -ai "$log"

exit 0
