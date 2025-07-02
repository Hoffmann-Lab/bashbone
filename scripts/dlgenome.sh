#! /usr/bin/env bash
# (c) Konstantin Riege

source "$(dirname "$0")/../bashbone_lite.sh" -a "$@" || exit 1

############################################

usage(){
	cat <<- EOF
		DESCRIPTION
		$(basename $0) downloads most recent human or mouse genome, annotation, gene ontology, orthologs, dbSNP

		VERSION
		1.0.1

		SYNOPSIS
		$(basename $0) -r [hg19|hg38|mm10|mm11] -o <PATH> -t <THREADS> [-g|-c] [-a] [-d|-m|-e] [-s|-n|-u|-k]

		INPUT OPTIONS
		-h | --help               : prints this message
		-t | --threads [value]    : threads - predicted default: $threads
		-o | --out [path]         : output directory - default: $PWD
		-r | --reference [string] : choose one of hg19 hg38 mm10 mm11

		(the following options require bgzip in PATH)
		-g | --genome             : download Ensembl genome
		-c | --ctat               : switch to CTAT genome and indices (~30GB)

		(the following options require bgzip and optional tabix in PATH)
		-a | --annotation         : download Ensembl gtf

		(the following options require R in PATH)
		-d | --descriptions       : download Ensembl gene description and Ensembl ontology information
		-m | --msigdb             : switch to MSigDB gene ontology information and other collections (recommended since Ensembl has too many and incomplete sets)
		-e | --enrichr            : switch to Enrichr gene ontology and other collections
		                          : NOTE: for mouse data, genes will be substituted by Ensembl orthologs if possible

		(the following options require bcftools, bgzip, tabix in PATH)
		-s | --dbsnp              : download Ensembl common dbSNP (hg), EVA validated SNPs (mm) respectively (mm10 only)
		                          : NOTE: for hg, each orignal entry is listed in 1000Genomes and contain its african, american, european, asian population MAF
		                            -> filtered for max(MAF)>0.01
		-n | --ncbi               : switch to NCBI common dbSNP (hg), validated SNPs (mm) respectively (mm10 only)
		                          : NOTE: for hg, an original entry may be listed in 1000Genomes. tagged as common if representative major project MAF>0.01
		                            -> filtered for COMMON and 1000Genomes
		-u | --ucsc               : switch to UCSC common dbSNP v142 (mm10 only)
		                          : NOTE: each entry is listed in 1000Genomes and is common due to a representative major project MAF>0.01
		                           -> although UCSC uses NCBI sources, selection of variants slightly differs
		-k | --gatk               : switch to GATK SNP bundle (hg only)

		REFERENCES
		(c) Konstantin Riege
		konstantin.riege{a}leibniz-fli{.}de
	EOF

	return 1
}

################################################################################
# GO
################################################################################

alias dlgenome::_go.ensembl="_bashbone_wrapper dlgenome::_go.ensembl"
function dlgenome::_go.ensembl(){
	out="$outdir/$genome.fa.gtf"

	tdir="$(mktemp -d -p "${TMPDIR:-/tmp}" XXXXXXXXXX.rlibs)"
	export R_LIBS="$tdir"
	mkdir -p "$R_LIBS"
	# listDatasets(useEnsembl(biomart = "genes"))
	# listDatasets(useEnsembl(biomart = "ensembl"))

	[[ $(R --version | head -1 | awk '$3<3.5{print 1}') ]] && {
		cat <<- EOF > "$tdir/download.R"
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("biocLite", quietly = TRUE)) source("https://bioconductor.org/biocLite.R")
				biocLite("biomaRt", suppressUpdates=TRUE)
			}
		EOF
	} || {
		cat <<- EOF > "$tdir/download.R"
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cloud.r-project.org", Ncpus=$threads, clean=T)
				library("BiocManager")
				BiocManager::install(c("biomaRt"), Ncpus=$threads, clean=T)
			}
		EOF
	}
	cat <<- EOF >> "$tdir/download.R"
		library("biomaRt")
		v <- "$version"
		for (try in 1:10){
			tryCatch({
				if ("$version" == "latest"){
					ensembl <- useEnsembl(biomart="genes", dataset="$dataset")
					v <- listEnsemblArchives()\$version[2]
				} else {
					ensembl <- useEnsembl(biomart="genes", dataset="$dataset", version=v)
				}
				try <- 0
				break
			}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 10"))})
		}
		if(try!=0) quit("no",1)
		cat(paste0(":INFO: v",v," from Ensembl\n"))

		goids <- data.frame()
		descriptions <- data.frame()
		for (chr in grep("^(MT|X|Y|\\\d+)$",listFilterOptions(mart = ensembl, filter = "chromosome_name"), value=T, perl=T)){
			cat(paste0(":INFO: working on chr",chr,"\n"))
			df <- data.frame()
			for (try in 1:10){
				tryCatch({
					df <- getBM(mart=ensembl, attributes=c("ensembl_gene_id","go_id","namespace_1003","name_1006"), filters=c("chromosome_name"), values=list(chr))
					try <- 0
					break
				}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 10"))})
			}
			if(try!=0) quit("no",1)
			goids <- rbind(goids, df)
			df <- data.frame()
			for (try in 1:10){
				tryCatch({
					df <- getBM(mart=ensembl, attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters=c("chromosome_name"), values=list(chr))
					try <- 0
					break
				}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 10"))})
			}
			if(try!=0) quit("no",1)
			descriptions <- rbind(descriptions, df)
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

		orthos <- data.frame()
		for (v in as.vector(na.omit(as.integer(listEnsemblArchives()\$version)))){
			cat(paste0(":INFO: working on orthologs from Ensembl v",v,"\n"))
			for (try in 1:2){
				tryCatch({
					ensembl_Hs <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=v)
					ensembl_Mm <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl", version=v)
					orthos <- getLDS(mart = ensembl_Hs, martL = ensembl_Mm, attributes = c("ensembl_gene_id","external_gene_name"), attributesL = c("ensembl_gene_id","external_gene_name"))
					try <- 0
					break
				}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 2"))})
			}
			if(try==0) break
		}
		if(try!=0) quit("no",1)
		write.table(orthos,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.orthologs")
	EOF

	echo ":INFO: downloading gene ontology and descriptions"
	Rscript "$tdir/download.R"

	awk -F '\t' -v OFS='\t' '$2~/\S/ && $4~/\S/ && tolower($2)==tolower($4) {print $1,$3}' "$out.orthologs" > "$tdir/orthologs"
	cut -f 1 "$tdir/orthologs" | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" | cut -f 2 | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" > "$out.orthologs.unique"
	cut -f 1 "$out.orthologs.unique" | grep -vFw -f - "$out.orthologs" | grep -vFw -f <(cut -f 2 "$out.orthologs.unique") | cut -f 1,3 > "$tdir/orthologs"
	cut -f 1 "$tdir/orthologs" | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" | cut -f 2 | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" >> "$out.orthologs.unique"

	return 0
}

alias dlgenome::_go.enrichr="_bashbone_wrapper dlgenome::_go.enrichr"
function dlgenome::_go.enrichr(){
	out="$outdir/$genome.fa.gtf"

	tdir="$(mktemp -d -p "${TMPDIR:-/tmp}" XXXXXXXXXX.rlibs)"
	export R_LIBS="$tdir"
	mkdir -p "$R_LIBS"

	[[ $(R --version | head -1 | awk '$3<3.5{print 1}') ]] && {
		cat <<- EOF > "$tdir/download.R"
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("biocLite", quietly = TRUE)) source("https://bioconductor.org/biocLite.R")
				biocLite("biomaRt", suppressUpdates=TRUE)
			}
		EOF
	} || {
		cat <<- EOF > "$tdir/download.R"
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cloud.r-project.org", Ncpus=$threads, clean=T)
				library("BiocManager")
				BiocManager::install(c("biomaRt"), Ncpus=$threads, clean=T)
			}
		EOF
	}
	cat <<- EOF >> "$tdir/download.R"
		library("biomaRt")
		v <- "$version"
		for (try in 1:10){
			tryCatch({
				if ("$version" == "latest"){
					ensembl <- useEnsembl(biomart="genes", dataset="$dataset")
					v <- listEnsemblArchives()\$version[2]
				} else {
					ensembl <- useEnsembl(biomart="genes", dataset="$dataset", version=v)
				}
				try <- 0
				break
			}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 10"))})
		}
		if(try!=0) quit("no",1)
		cat(paste0(":INFO: v",v," from Ensembl\n"))

		descriptions <- data.frame()
		for (chr in grep("^(MT|X|Y|\\\d+)$",listFilterOptions(mart = ensembl, filter = "chromosome_name"), value=T, perl=T)){
			cat(paste0(":INFO: working on chr",chr,"\n"))
			df <- data.frame()
			for (try in 1:10){
				tryCatch({
					df <- getBM(mart=ensembl, attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters=c("chromosome_name"), values=list(chr))
					try <- 0
					break
				}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 10"))})
			}
			if(try!=0) quit("no",1)
			descriptions <- rbind(descriptions, df)
		}
		descriptions[,2][descriptions[2]==""] <- descriptions[,1][descriptions[2]==""]
		write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")

		orthos <- data.frame()
		for (v in as.vector(na.omit(as.integer(listEnsemblArchives()\$version)))){
			cat(paste0(":INFO: working on orthologs from Ensembl v",v,"\n"))
			for (try in 1:2){
				tryCatch({
					ensembl_Hs <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=v)
					ensembl_Mm <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl", version=v)
					orthos <- getLDS(mart = ensembl_Hs, martL = ensembl_Mm, attributes = c("ensembl_gene_id","external_gene_name"), attributesL = c("ensembl_gene_id","external_gene_name"))
					try <- 0
					break
				}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 2"))})
			}
			if(try==0) break
		}
		if(try!=0) quit("no",1)
		write.table(orthos,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.orthologs")
	EOF

	echo ":INFO: downloading gene ontology and descriptions"
	Rscript "$tdir/download.R"

	awk -F '\t' -v OFS='\t' '$2~/\S/ && $4~/\S/ && tolower($2)==tolower($4) {print $1,$3}' "$out.orthologs" > "$tdir/orthologs"
	cut -f 1 "$tdir/orthologs" | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" | cut -f 2 | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" > "$out.orthologs.unique"
	cut -f 1 "$out.orthologs.unique" | grep -vFw -f - "$out.orthologs" | grep -vFw -f <(cut -f 2 "$out.orthologs.unique") | cut -f 1,3 > "$tdir/orthologs"
	cut -f 1 "$tdir/orthologs" | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" | cut -f 2 | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" >> "$out.orthologs.unique"

	# wget -O "$tdir/ncbi.info.gz" -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
	# gzip -dc "$tdir/ncbi.info.gz" | perl -lanE 'next unless $F[5]=~/Ensembl:(ENSG\d+)/; say "$F[2]\t$1"' > "$tdir/ncbi.info"
	url="https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" | bgzip -kcd -@ 1 | perl -lanE 'next unless $F[5]=~/Ensembl:(ENSG\d+)/; say "$F[2]\t$1"' > "$tdir/ncbi2ensembl"

	cat <<-EOF > "$out.go.README"
		$(date)
		$USER
	EOF

	[[ "$msig" == "Hs" ]] && msig="Human" || msig="Mouse"
	rm -f "$out.go"
	for collection in GO_Biological_Process GO_Molecular_Function GO_Cellular_Component Reactome GWAS_Catalog ChEA CellMarker WikiPathways KEGG; do
		version=$(date +%Y)
		while :; do
			url="https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=${collection}_$version"
			[[ "$collection" =~ ^(WikiPathways|KEGG) ]] && url+="_$msig"
			wget -O "$tdir/$collection.gmt" -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" && echo ":INFO: working on $collection $version" && break
			((version--))
		done
		perl -F'\t' -slanE '
			BEGIN{
				open F,"<$ortho" or die $!;
				while(<F>){
					chomp;
					@F=split/\t/;
					$h2m{$F[0]}=$F[1];
				}
				close F;
				open F,"<$info" or die $!;
				while(<F>){
					chomp;
					@F=split/\t/;
					$g2i{$F[0]}=$F[1];
				}
				close F;
			}
			$n=$F[0];
			if($F[0]=~s/\s+\((GO:\d+)\)$//){
				$n=$1;
				$collection=lc($collection=~s/^GO_//r);
			}
			for (@F[2..$#F]){
				$i=$g2i{$_};
				next unless $i;
				$i=$h2m{$i} if $species eq "Mouse";
				next unless $i;
				say join"\t",($i,$n=~s/\s+/_/gr,$collection,$F[0]);
			}
		' -- -info="$tdir/ncbi2ensembl" -species="$msig" -ortho="$out.orthologs.unique" -collection="$collection" "$tdir/$collection.gmt" >> "$out.go"
		echo "Enrichr $collection v$version" >> "$out.go.README"
	done

	return 0
}

alias dlgenome::_go.msigdb="_bashbone_wrapper dlgenome::_go.msigdb"
function dlgenome::_go.msigdb(){
	out="$outdir/$genome.fa.gtf"

	tdir="$(mktemp -d -p "${TMPDIR:-/tmp}" XXXXXXXXXX.rlibs)"
	export R_LIBS="$tdir"
	mkdir -p "$R_LIBS"

	[[ $(R --version | head -1 | awk '$3<3.5{print 1}') ]] && {
		cat <<- EOF > "$tdir/download.R"
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("biocLite", quietly = TRUE)) source("https://bioconductor.org/biocLite.R")
				biocLite("biomaRt", suppressUpdates=TRUE)
			}
		EOF
	} || {
		# does not work anymore. use default: NULL. if ("$version" == "latest") v <- NULL
		cat <<- EOF > "$tdir/download.R"
			if (!requireNamespace("biomaRt", quietly = TRUE)) {
				if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos="http://cloud.r-project.org", Ncpus=$threads, clean=T)
				library("BiocManager")
				BiocManager::install(c("biomaRt"), Ncpus=$threads, clean=T)
			}
		EOF
	}
	cat <<- EOF >> "$tdir/download.R"
		library("biomaRt")
		v <- "$version"
		for (try in 1:10){
			tryCatch({
				if ("$version" == "latest"){
					ensembl <- useEnsembl(biomart="genes", dataset="$dataset")
					v <- listEnsemblArchives()\$version[2]
				} else {
					ensembl <- useEnsembl(biomart="genes", dataset="$dataset", version=v)
				}
			}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 10"))})
			try <- 0
			break
		}
		if(try!=0) quit("no",1)
		cat(paste0(":INFO: v",v," from Ensembl\n"))

		descriptions <- data.frame()
		for (chr in grep("^(MT|X|Y|\\\d+)$",listFilterOptions(mart = ensembl, filter = "chromosome_name"), value=T, perl=T)){
			cat(paste0(":INFO: working on chr",chr,"\n"))
			df <- data.frame()
			for (try in 1:10){
				tryCatch({
					df <- getBM(mart=ensembl, attributes=c("ensembl_gene_id","external_gene_name","gene_biotype","description"), filters=c("chromosome_name"), values=list(chr))
					try <- 0
					break
				}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 10"))})
			}
			if(try!=0) quit("no",1)
			descriptions <- rbind(descriptions, df)
		}
		descriptions[,2][descriptions[2]==""] <- descriptions[,1][descriptions[2]==""]
		write.table(descriptions,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.info")

		orthos <- data.frame()
		for (v in as.vector(na.omit(as.integer(listEnsemblArchives()\$version)))){
			cat(paste0(":INFO: working on orthologs from Ensembl v",v,"\n"))
			for (try in 1:2){
				tryCatch({
					ensembl_Hs <- useEnsembl(biomart="genes", dataset="hsapiens_gene_ensembl", version=v)
					ensembl_Mm <- useEnsembl(biomart="genes", dataset="mmusculus_gene_ensembl", version=v)
					orthos <- getLDS(mart = ensembl_Hs, martL = ensembl_Mm, attributes = c("ensembl_gene_id","external_gene_name"), attributesL = c("ensembl_gene_id","external_gene_name"))
					try <- 0
					break
				}, error = function(e){message(paste0(":WARNING: connection failed at try ",try," out of 2"))})
			}
			if(try==0) break
		}
		if(try!=0) quit("no",1)
		write.table(orthos,quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t",file="$out.orthologs")
	EOF

	echo ":INFO: downloading gene ontology and descriptions"
	# Rscript "$tdir/download.R"

	awk -F '\t' -v OFS='\t' '$2~/\S/ && $4~/\S/ && tolower($2)==tolower($4) {print $1,$3}' "$out.orthologs" > "$tdir/orthologs"
	cut -f 1 "$tdir/orthologs" | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" | cut -f 2 | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" > "$out.orthologs.unique"
	cut -f 1 "$out.orthologs.unique" | grep -vFw -f - "$out.orthologs" | grep -vFw -f <(cut -f 2 "$out.orthologs.unique") | cut -f 1,3 > "$tdir/orthologs"
	cut -f 1 "$tdir/orthologs" | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" | cut -f 2 | sort | uniq -c | awk '$1==1{print $2}' | grep -Fw -f - "$tdir/orthologs" >> "$out.orthologs.unique"

	version=$(curl -s https://data.broadinstitute.org/gsea-msigdb/msigdb/release/ | grep -oE ">[0-9][^<]+($msig)" | sed 's/^>//' | sort -Vr | head -1)
	[[ "$msig" == "Hs" ]] && msig="human" || msig="mouse"
	wget -P "$tdir/" -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping --recursive --no-directories --no-parent --level=1 --reject 'index.htm*' --accept-regex ".*\.go\.(bp|mf|cc)\.v$version.symbols.gmt" "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/$version/"
	wget -O "$tdir/go.obo" -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "http://purl.obolibrary.org/obo/go.obo"

	grep -Fx '[Term]' -A 3 "$tdir/go.obo" | grep -vFx -- '--' | paste - - - - | cut -f 2- | sed -E 's/(id:|name:|namespace:)\s//g' | perl -F'\t' -lane '$F[1]=~s/#//; $F[1]=~s/^OBSOLETE\s*//i; print join"\t",@F' > "$tdir/go.terms.orig"
	grep -Fx '[Term]' -A 3 "$tdir/go.obo" | grep -vFx -- '--' | paste - - - - | cut -f 2- | sed -E 's/(id:|name:|namespace:)\s//g' | perl -F'\t' -lane '$F[1]=uc($F[1]=~s/\W+/_/gr); $F[1]=~s/^_+//; $F[1]=~s/_+$//; $F[1]=~s/_+/_/g; $F[1]=~s/^OBSOLETE_//; print join"\t",@F' > "$tdir/go.terms"

	rm -f "$out.go"
	for gmt in "$tdir/"*.v$version.symbols.gmt; do
		cut -f 1 "$gmt" | sed -E 's/^GO.._//' > "$tdir/gmt.terms"
		grep -Fw -f "$tdir/gmt.terms" "$tdir/go.terms" > "$tdir/gmt.terms.parsed"
		grep -vFx -f <(cut -f 2 "$tdir/gmt.terms.parsed") "$tdir/gmt.terms" > "$tdir/gmt.terms.missing"
		c=$(basename "$gmt" | cut -d . -f 3 | perl -lane 'print "GO".uc($_)')
		x=$(wc -l < "$tdir/gmt.terms.missing")
		i=0
		while read -r n; do
			echo "requesting $((++i)) of $x missing $c terms" >&2
			g=$(curl -s "https://www.gsea-msigdb.org/gsea/msigdb/$msig/geneset/${c}_$n" | grep -oP 'GO:\d+</td>' | cut -d '<' -f 1 | sort -u)
			[[ $g ]] || echo ":WARNING: no term found for $n"
			grep -m 1 -Fw $g "$tdir/go.terms" | awk -F '\t' -v n=$n -v OFS='\t' '{print $1,n,$3}' >> "$tdir/gmt.terms.parsed"
		done < "$tdir/gmt.terms.missing"
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
		' -- -terms="$tdir/gmt.terms.parsed" -go="$tdir/go.terms.orig" -info="$out.info" "$gmt" >> "$out.go"
	done

	while read -r collection gmt; do
		echo ":INFO: working on $gmt"
		wget -O - -c -q --show-progress --progress=bar:force --timeout=60 --waitretry=10 --tries=10 --retry-connrefused --timestamping "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/$version/$gmt" | perl -F'\t' -slanE '
			BEGIN{
				open F,"<$info" or die $!;
				while(<F>){
					chomp;
					@F=split/\t/;
					$F[1]=$F[0] unless $F[1];
					$g2i{$F[1]}=$F[0];
				}
				close F;
			}
			for (@F[2..$#F]){
				$i=$g2i{$_};
				next unless $i;
				say join"\t",($i,$F[0],$collection,lc($F[0]=~s/_+/ /gr));
			}
		' -- -info="$out.info" -collection="$collection" >> "$out.go"
	done < <(curl -s "https://www.gsea-msigdb.org/gsea/msigdb/$msig/collections.jsp" | grep -e "name=" -e symbols.gmt | grep -B 1 symbols.gmt | grep -v -Fx -- '--' | sed -E 's@.*/([^/]+symbols\.gmt).*@\1@;s@.*>([^>]+)<[^<]+$@\1@' | paste - - | grep -v -E -e "^(a|A)ll gene sets" -e "(o|O)ntology gene sets" -e "GO" | perl -F'\t' -lane '$F[0]=~s/\s+subset\s+of\s+\S+//; $F[0]=~s/^[^:]:\s+//; $F[0]=~s/\W+/_/g; $F[0]="MSigDB_".$F[0]; print join"\t",@F')

	cat <<-EOF > "$out.go.README"
		$(date)
		$USER
		MSigDB go v$version
	EOF
	return 0
}

alias dlgenome::hg19.go="_bashbone_wrapper dlgenome::hg19.go"
function dlgenome::hg19.go(){
	cd "$outdir"
	genome=GRCh37
	dataset="hsapiens_gene_ensembl"
	version="GRCh37"
	msig=Hs
	dlgenome::_go.$godb
	return 0
}

alias dlgenome::hg38.go="_bashbone_wrapper dlgenome::hg38.go"
function dlgenome::hg38.go(){
	cd "$outdir"
	genome=GRCh38
	dataset="hsapiens_gene_ensembl"
	version="latest"
	msig=Hs
	dlgenome::_go.$godb
	return 0
}

alias dlgenome::mm10.go="_bashbone_wrapper dlgenome::mm10.go"
function dlgenome::mm10.go(){
	cd "$outdir"
	genome=GRCm38
	dataset="mmusculus_gene_ensembl"
	version="102"
	msig=Mm
	dlgenome::_go.$godb
	return 0
}

alias dlgenome::mm11.go="_bashbone_wrapper dlgenome::mm11.go"
function dlgenome::mm11.go(){
	cd "$outdir"
	genome=GRCm39
	dataset="mmusculus_gene_ensembl"
	version="latest"
	msig=Mm
	dlgenome::_go.$godb
	return 0
}

################################################################################
#  GENOME
################################################################################

alias ddlgenome::hg19.genome="_bashbone_wrapper dlgenome::hg19.genome"
function dlgenome::hg19.genome() {
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading $genome genome"
	url="https://ftp.ensembl.org/pub/grch37/"
	# url+="$(curl -s "$url" | grep -oE 'release-[0-9]+' | sort -Vr | head -1)/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome"
	url+="$(curl -s "${url/https/ftp}" | grep -E "current\s*->\s*" | awk '{print $NF}')/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome"
	version=$(grep -oE 'release-[0-9]+' <<< "$url" | cut -d '-' -f 2-)
	echo ":INFO: v$version from Ensembl"

	rm -f "$genome.fa"
	for i in MT {1..22} X Y; do
		echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i)
		wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url.$i.fa.gz" | bgzip -kcd -@ 1  | sed "s/^>MT.*/>chrM/;t;s/^>.*/>chr$i/" >> "$genome.fa"
	done


	cat <<- EOF > "$genome.fa.README"
		$(date)
		$USER
		Ensembl v$version $genome
		$(dirname "$url")
	EOF

	return 0
}

alias dlgenome::hg38.genome="_bashbone_wrapper dlgenome::hg38.genome"
function dlgenome::hg38.genome() {
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading $genome genome"
	url="https://ftp.ensembl.org/pub/"
	# url+="$(curl -s "$url" | grep -oE 'release-[0-9]+' | sort -Vr | head -1)/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome"
	url+="$(curl -s "${url/https/ftp}" | grep -E "current\s*->\s*" | awk '{print $NF}')/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome"
	version=$(grep -oE 'release-[0-9]+' <<< "$url" | cut -d '-' -f 2-)
	echo ":INFO: v$version from Ensembl"

	rm -f "$genome.fa"
	for i in MT {1..22} X Y; do
		echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i)
		wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url.$i.fa.gz" | bgzip -kcd -@ 1 | sed "s/^>MT.*/>chrM/;t;s/^>.*/>chr$i/" >> "$genome.fa"
	done

	cat <<- EOF > "$genome.fa.README"
		$(date)
		$USER
		Ensembl v$version $genome
		$(dirname "$url")
	EOF

	return 0
}

alias dlgenome::mm10.genome="_bashbone_wrapper dlgenome::mm10.genome"
function dlgenome::mm10.genome() {
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: downloading $genome genome"
	url="https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome"
	version=$(grep -oE 'release-[0-9]+' <<< "$url" | cut -d '-' -f 2-)
	echo ":INFO: v$version from Ensembl"

	rm -f "$genome.fa"
	for i in MT {1..19} X Y; do
		echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i)
		wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url.$i.fa.gz" | bgzip -kcd -@ 1 | sed "s/^>MT.*/>chrM/;t;s/^>.*/>chr$i/" >> "$genome.fa"
	done

	cat <<- EOF > "$genome.fa.README"
		$(date)
		$USER
		Ensembl v$version $genome
		$(dirname "$url")
	EOF

	return 0
}

alias dlgenome::mm11.genome="_bashbone_wrapper dlgenome::mm11.genome"
function dlgenome::mm11.genome() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: downloading $genome genome"
	url="https://ftp.ensembl.org/pub/"
	# url+="$(curl -s "$url" | grep -oE 'release-[0-9]+' | sort -Vr | head -1)/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome"
	url+="$(curl -s "${url/https/ftp}" | grep -E "current\s*->\s*" | awk '{print $NF}')/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome"
	version=$(grep -oE 'release-[0-9]+' <<< "$url" | cut -d '-' -f 2-)
	echo ":INFO: v$version from Ensembl"

	rm -f "$genome.fa"
	for i in MT {1..19} X Y; do
		echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i)
		wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url.$i.fa.gz" | bgzip -kcd -@ 1 | sed "s/^>MT.*/>chrM/;t;s/^>.*/>chr$i/" >> "$genome.fa"
	done

	cat <<- EOF > "$genome.fa.README"
		$(date)
		$USER
		Ensembl v$version $genome
		$(dirname "$url")
	EOF

	return 0
}

################################################################################
# GTF
################################################################################

alias dlgenome::hg19.gtf="_bashbone_wrapper dlgenome::hg19.gtf"
function dlgenome::hg19.gtf() {
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading annotation"
	url="https://ftp.ensembl.org/pub/grch37/"
	# url+="$(curl -s "$url" | grep -oE 'release-[0-9]+' | sort -Vr | head -1)/gtf/homo_sapiens/"
	url+="$(curl -s "${url/https/ftp}" | grep -E "current\s*->\s*" | awk '{print $NF}')/gtf/homo_sapiens/"
	url+="$(curl -s "$url" | grep -oE 'Homo_sapiens.GRCh37.[0-9]+.chr.gtf.gz' | sort -Vr | head -1)"
	version="$(basename "$url" | cut -d . -f 3)"
	echo ":INFO: v$version from Ensembl"

	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" \
		| bgzip -kcd -@ 1 \
		| sed -n '/^#/!{s/^MT/chrM/p;t;s/^/chr/p}' \
		| LC_ALL=C sort -S 4096M -T "${TMPDIR:-/tmp}" --parallel=$threads -k1,1 -k4,4n -k5,5n \
		| perl -e 'while(<STDIN>){@F=split; push $m{$F[0]}->@*,$_} print $m{$_}->@* for @ARGV' chrM chr{1..22} chrX chrY \
		| tee >(bgzip -kc -@ $threads > "$genome.fa.gtf.gz") > "$genome.fa.gtf" | cat
	[[ -x "$(command -v tabix)" ]] && tabix -f "$genome.fa.gtf.gz"

	cat <<- EOF > "$genome.fa.gtf.README"
		$(date)
		$USER
		Ensembl v$version $genome
		$url
	EOF

	return 0
}

alias dlgenome::hg38.gtf="_bashbone_wrapper dlgenome::hg38.gtf"
function dlgenome::hg38.gtf() {
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading annotation"

	url="https://ftp.ensembl.org/pub/"
	# url+="$(curl -s "$url" | grep -oE 'release-[0-9]+' | sort -Vr | head -1)/gtf/homo_sapiens/"
	url+="$(curl -s "${url/https/ftp}" | grep -E "current\s*->\s*" | awk '{print $NF}')/gtf/homo_sapiens/"
	url+="$(curl -s "$url" | grep -oE 'Homo_sapiens.GRCh38.[0-9]+.chr.gtf.gz' | sort -Vr | head -1)"
	version="$(basename "$url" | cut -d . -f 3)"
	echo ":INFO: v$version from Ensembl"

	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" \
		| bgzip -kcd -@ 1 \
		| sed -n '/^#/!{s/^MT/chrM/p;t;s/^/chr/p}' \
		| LC_ALL=C sort -S 4096M -T "${TMPDIR:-/tmp}" --parallel=$threads -k1,1 -k4,4n -k5,5n \
		| perl -e 'while(<STDIN>){@F=split; push $m{$F[0]}->@*,$_} print $m{$_}->@* for @ARGV' chrM chr{1..22} chrX chrY \
		| tee >(bgzip -kc -@ $threads > "$genome.fa.gtf.gz") > "$genome.fa.gtf" | cat
	[[ -x "$(command -v tabix)" ]] && tabix -f "$genome.fa.gtf.gz"

	cat <<- EOF > "$genome.fa.gtf.README"
		$(date)
		$USER
		Ensembl v$version $genome
		$url
	EOF

	return 0
}

alias dlgenome::mm10.gtf="_bashbone_wrapper dlgenome::mm10.gtf"
function dlgenome::mm10.gtf() {
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: downloading annotation"
	# wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping --glob=on 'ftp://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.*.chr.gtf.gz'

	# echo ":INFO: extracting annotation"
	# gzip -dc Mus_musculus.GRCm38.*.chr.gtf.gz | perl -F'\t' -lane 'next unless $F[0]=~/^(\d+|X|Y|MT)$/; $F[0]="M" if $F[0] eq "MT"; $F[0]="chr$F[0]"; unless($f eq $F[0]){close O; $f=$F[0]; open O,">$f.gtf";} print O join("\t",@F); END{close O};'

	# echo ":INFO: sorting annotation"
	# rm -f $genome.fa.gtf
	# for i in M {1..19} X Y; do
	# 	LC_ALL=C sort -S 1000M -T "${TMPDIR:-/tmp}" --parallel=$threads -k4,4n -k5,5n chr$i.gtf >> $genome.fa.gtf
	# done

	url="https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/"
	url+=$(curl -s "$url" | grep -oE 'Mus_musculus.GRCm38.[0-9]+.chr.gtf.gz' | sort -Vr | head -1)
	version="$(basename "$url" | cut -d . -f 3)"
	echo ":INFO: v$version from Ensembl"

	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" \
		| bgzip -kcd -@ 1 \
		| sed -n '/^#/!{s/^MT/chrM/p;t;s/^/chr/p}' \
		| LC_ALL=C sort -S 4096M -T "${TMPDIR:-/tmp}" --parallel=$threads -k1,1 -k4,4n -k5,5n \
		| perl -e 'while(<STDIN>){@F=split; push $m{$F[0]}->@*,$_} print $m{$_}->@* for @ARGV' chrM chr{1..19} chrX chrY \
		| tee >(bgzip -kc -@ $threads > "$genome.fa.gtf.gz") > "$genome.fa.gtf" | cat
	[[ -x "$(command -v tabix)" ]] && tabix -f "$genome.fa.gtf.gz"

	cat <<- EOF > "$genome.fa.gtf.README"
		$(date)
		$USER
		Ensembl v$version $genome
		$url
	EOF

	return 0
}

alias dlgenome::mm11.gtf="_bashbone_wrapper dlgenome::mm11.gtf"
function dlgenome::mm11.gtf() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: downloading annotation"
	url="https://ftp.ensembl.org/pub/"
	# url+="$(curl -s "$url" | grep -oE 'release-[0-9]+' | sort -Vr | head -1)/gtf/mus_musculus/"
	url+="$(curl -s "${url/https/ftp}" | grep -E "current\s*->\s*" | awk '{print $NF}')/gtf/mus_musculus/"
	url+="$(curl -s "$url" | grep -oE 'Mus_musculus.GRCm39.[0-9]+.chr.gtf.gz' | sort -Vr | head -1)"
	version="$(basename "$url" | cut -d . -f 3)"
	echo ":INFO: v$version from Ensembl"

	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" \
		| bgzip -kcd -@ 1 \
		| sed -n '/^#/!{s/^MT/chrM/p;t;s/^/chr/p}' \
		| LC_ALL=C sort -S 4096M -T "${TMPDIR:-/tmp}" --parallel=$threads -k1,1 -k4,4n -k5,5n \
		| perl -e 'while(<STDIN>){@F=split; push $m{$F[0]}->@*,$_} print $m{$_}->@* for @ARGV' chrM chr{1..19} chrX chrY \
		| tee >(bgzip -kc -@ $threads > "$genome.fa.gtf.gz") > "$genome.fa.gtf" | cat
	[[ -x "$(command -v tabix)" ]] && tabix -f "$genome.fa.gtf.gz"

	cat <<- EOF > "$genome.fa.gtf.README"
		$(date)
		$USER
		Ensembl v$version $genome
		$url
	EOF

	return 0
}

################################################################################
# DBSNP

##### old sources. don't use v150/151 anymore and further, don't store as uncompressed vcf
##### 151 gz is corrupt
# url=$(curl -sl ftp://ftp.ncbi.nih.gov/snp/organisms/ | grep -E 'human_[0-9]+_b150_GRCh38' | sort -V | tail -1)

# according to ncbi dbsnp 150 CAF info holds 1000genomes MAF only for REF and ALTs. from v151 on, common tagged from additional filtered by MAF>0.01
# from v152 on there is no CAF anymore, but single project MAFs (see json files on ftp server)
# ncbi 150: n=37463647 - ucsc says ncbi holds 38M common and 23M from 1000genomes. now apply MAF i.e. CAF filter.
#           n=14851093   via bcftools view -H -i 'CAF[0]!~"^0\.99" && OTH=0 && CFL=0' GRCh38.fa.dbSNP150common.vcf.gz
# ucsc 150: n=14375783   (https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp150Common.txt.gz | awk '$2~/^chr[^_]+$/')

##### new sources
#####
# ncbi 155: 21155945 - versus 14133245 common coming from ucsc always containing 1000Genomes (https://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb)
#                      according to ucsc: ncbi has 15M variants with MAF>0.01 in the 1000Genomes Phase 3 dataset
#           14203886   via bcftools view -H https://ftp.ncbi.nlm.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz | grep -F 1000Genomes | perl -lane '$F[-1]=~/FREQ=([^;]+)/; @m=$1=~/:(0\.\d+)/g; print $_ if (sort {$a <=> $b} @m)[0]<=0.99'
# ncbi 156: 21340526 -
#           15194872   via .. analogous ..
# embl 156: 82663547 - within 1000GENOMES-phase_3.vcf.gz -> some have MAF others ethnicity wise AF listed (africa, america, asia, europa, south-asia)
#           15009538   via bcftools view -i '(MAF!="." && MAX(MAF)>0.01) || (MAF="." && (MAX(AFR)>0.01 || MAX(AMR)>0.01 || MAX(EAS)>0.01 || MAX(EUR)>0.01 || MAX(SAS)>0.01))'

# ATTENTION ncbi and ensembl/ucsc dbsnp slightly differ at indel definitions
# chr1 12345 id12345 AT ATT
# chr1 12344 id12345  T  TT

##### gatk bundles
##### https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle
# HG38
# https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
# https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi
# https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
# https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi
## mutect germline-resource
# https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz
# https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi
# ===> check out also https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/gnomad/v4.1/genomes/
## to check for contamination
# https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz
# https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi
# ===> full exac here https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_genotype/
## pon
# https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz
# https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi
# HG19
# https://storage.googleapis.com/gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
# https://storage.googleapis.com/gatk-legacy-bundles/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz
## mutect germline-resource
# https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf
## to check for contamination
# https://storage.googleapis.com/gatk-best-practices/somatic-b37/small_exac_common_3.vcf
## pon
# https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf

################################################################################

alias dlgenome::hg19.dbsnp.gatk="_bashbone_wrapper dlgenome::hg19.dbsnp.gatk"
function dlgenome::hg19.dbsnp.gatk() {
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading GATK SNP bundle"

	tfile="$(mktemp -p "${TMPDIR:-/tmp}" XXXXXXXXXX.vcf.gz)"

	url="https://storage.googleapis.com/gatk-legacy-bundles/b37/1000G_phase1.snps.high_confidence.b37.vcf.gz"
	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" | bgzip -kcd -@ 1 | sed -E 's/\s+$//' | bgzip -kc -@ $threads > "$tfile"
	tabix -f -p vcf "$tfile"
	{	bcftools view -h "$tfile" 2> >(grep -v "PL should be declared" >&2 || true) | grep -vF -e '##contig='
		for i in MT {1..22} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H "$tfile" $i 2> >(grep -v "PL should be declared" >&2 || true) | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$(basename "$url")"
	tabix -f -p vcf "$(basename "$url")"
	ln -sfnr "$(basename "$url")" SNPs_gold.vcf.gz
	ln -sfnr "$(basename "$url").tbi" SNPs_gold.vcf.gz.tbi


	url="https://storage.googleapis.com/gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"
	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" | bgzip -kcd -@ 1 | sed -E 's/\s+$//' | bgzip -kc -@ $threads > "$tfile"
	tabix -f -p vcf "$tfile"
	{	bcftools view -h "$tfile" | grep -vF -e '##contig='
		for i in MT {1..22} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H "$tfile" $i | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$(basename "$url")"
	tabix -f -p vcf "$(basename "$url")"
	ln -sfnr "$(basename "$url")" InDel_gold.vcf.gz
	ln -sfnr "$(basename "$url").tbi" InDel_gold.vcf.gz.tbi


	url="https://storage.googleapis.com/gatk-best-practices/somatic-b37/small_exac_common_3.vcf"
	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" | bgzip -kc -@ $threads > "$tfile"
	tabix -f -p vcf "$tfile"
	{	bcftools view -h "$tfile" | grep -vF -e '##contig='
		for i in MT {1..22} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H "$tfile" $i | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$(basename "$url").gz"
	tabix -f -p vcf "$(basename "$url").gz"
	ln -sfnr "$(basename "$url").gz" ExAC_common.vcf.gz
	ln -sfnr "$(basename "$url").gz.tbi" ExAC_common.vcf.gz.tbi


	url="https://storage.googleapis.com/gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf"
	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" | bgzip -kc -@ $threads > "$tfile"
	tabix -f -p vcf "$tfile"
	{	bcftools view -h "$tfile" | grep -vF -e '##contig='
		for i in MT {1..22} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H "$tfile" $i | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$(basename "$url").gz"
	tabix -f -p vcf "$(basename "$url").gz"
	ln -sfnr "$(basename "$url").gz" gnomAD_with_AF.vcf.gz
	ln -sfnr "$(basename "$url").gz.tbi" gnomAD_with_AF.vcf.gz.tbi


	url="https://storage.googleapis.com/gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf"
	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" | bgzip -kc -@ $threads > "$tfile"
	tabix -f -p vcf "$tfile"
	{	bcftools view -h "$tfile" | grep -vF -e '##contig='
		for i in MT {1..22} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H "$tfile" $i | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$(basename "$url").gz"
	tabix -f -p vcf "$(basename "$url").gz"
	ln -sfnr "$(basename "$url").gz" PON.vcf.gz
	ln -sfnr "$(basename "$url").gz.tbi" PON.vcf.gz.tbi
}

alias dlgenome::hg19.dbsnp.ncbi="_bashbone_wrapper dlgenome::hg19.dbsnp.ncbi"
function dlgenome::hg19.dbsnp.ncbi() {
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading dbSNP"

	url="https://ftp.ncbi.nlm.nih.gov/snp/archive/"
	url+="$(curl -s "$url" | grep -oE 'b[0-9]+' | sort -Vr | head -1)/VCF/"
	url+="$(curl -s "$url" | grep -oE '[^">]+\.gz' | sort -V | head -1)"
	# sort -V -> .25.vcf.gz -> hg19
	# sort -Vr -> .[>30].vcf.gz -> hg38
	version="$(bcftools view -h "$url" | grep -m 1 -F dbSNP_BUILD_ID | sed -E 's/.*(dbSNP)_BUILD_ID=(.+)/\1 v.\2/')"
	echo ":INFO: v$version from NCBI"

	echo "rm -f '$(basename "$url")'*" >> "$BASHBONE_CLEANUP"
	{	bcftools view -h "$url" | grep -vF '##contig='
		for i in NC_012920.1 $(bcftools view -h "$url" | grep -F '##contig=' | grep -oE "NC_[^>]+" | grep -vF NC_012920.1 | sort -V); do
			echo ":INFO: working on "$(sed -E 's/^NC_012920\S+/chrM/;t;s/^NC_000023\S+/chrX/;t;s/^NC_000024\S+/chrY/;t;s/^NC_000+([0-9]+)\.[0-9]+/chr\1/' <<< $i) >&2
			bcftools view --threads $threads -H "$url" $i 2> >(grep -v "Extreme INFO" >&2 || true) | grep -F COMMON | grep -F 1000Genomes | sed -E 's/^NC_012920\S+/chrM/;t;s/^NC_000023\S+/chrX/;t;s/^NC_000024\S+/chrY/;t;s/^NC_000+([0-9]+)\.[0-9]+/chr\1/'
		done
	} | bgzip -kc -@ $threads > "$genome.fa.vcf.gz"
	tabix -f -p vcf "$genome.fa.vcf.gz"

	cat <<- EOF >> "$genome.fa.vcf.README"
		$(date)
		$USER
		NCBI v$version $genome
		COMMON only i.e. MAF > 0.01 and listed in 1000Genomes
		$url
	EOF

	return 0
}

alias dlgenome::hg19.dbsnp.ensembl="_bashbone_wrapper dlgenome::hg19.dbsnp.ensembl"
function dlgenome::hg19.dbsnp.ensembl() {
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading dbSNP"
	# last time listed under releases was 2021 -> v110 (https://ftp.ensembl.org/pub/grch37/release-110/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz)
	url="https://ftp.ensembl.org/pub/grch37/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"
	version=$(bcftools view -h "$url" | grep -m 1 -oE dbSNP_[0-9]+ | cut -d _ -f 2)
	echo ":INFO: v$version from Ensembl"

	echo "rm '$(basename "$url")'*" >> "$BASHBONE_CLEANUP" # bcftools fetches csi or tbi from server
	{	bcftools view -h "$url" 2> >(grep -v "Invalid tag name" >&2 || true) | grep -vF -e '##contig=' -e "##INFO=<ID=HGMD-PUBLIC"
		for i in MT {1..22} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H -i '(MAF!="." && MAX(MAF)>0.01) || (MAF="." && (MAX(AFR)>0.01 || MAX(AMR)>0.01 || MAX(EAS)>0.01 || MAX(EUR)>0.01 || MAX(SAS)>0.01))' "$url" $i 2> >(grep -v "Invalid tag name" >&2 || true) | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$genome.fa.vcf.gz"
	tabix -f -p vcf "$genome.fa.vcf.gz"

	cat <<- EOF >> "$genome.fa.vcf.README"
		$(date)
		$USER
		Ensembl v$version $genome
		$url
		common only i.e. MAF > 0.01 or if MAF absent one of populations > 0.01
	EOF

	return 0
}

alias dlgenome::hg19.dbsnp.ucsc="_bashbone_wrapper dlgenome::hg19.dbsnp.ucsc"
function dlgenome::hg19.dbsnp.ucsc() {
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: UCSC dbSNP not supported"
	# https://groups.google.com/a/soe.ucsc.edu/g/genome-announce/c/40coV6ZVtFY
	# https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=dbSnp155Composite
	# http://hgdownload.soe.ucsc.edu/gbdb/hg19/snp/dbSnp155Common.bb

	return 0
}

alias dlgenome::hg38.dbsnp.gatk="_bashbone_wrapper dlgenome::hg38.dbsnp.gatk"
function dlgenome::hg38.dbsnp.gatk() {
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading GATK SNP bundle"

	url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
	{	bcftools view -h "$url" | grep -vF -e '##contig='
		for i in chrM chr{1..22} chrX chrY; do
			echo ":INFO: working on $i" >&2
			bcftools view --threads $threads -H "$url" $i
		done
	} | bgzip -kc -@ $threads > "$(basename "$url")"
	tabix -f -p vcf "$(basename "$url")"
	ln -sfnr "$(basename "$url")" SNPs_gold.vcf.gz
	ln -sfnr "$(basename "$url").tbi" SNPs_gold.vcf.gz.tbi


	# for BQSR also this can be supplied https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
	url="https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
	{	bcftools view -h "$url" | grep -vF -e '##contig='
		for i in chrM chr{1..22} chrX chrY; do
			echo ":INFO: working on $i" >&2
			bcftools view --threads $threads -H "$url" $i
		done
	} | bgzip -kc -@ $threads > "$(basename "$url")"
	tabix -f -p vcf "$(basename "$url")"
	ln -sfnr "$(basename "$url")" InDel_gold.vcf.gz
	ln -sfnr "$(basename "$url").tbi" InDel_gold.vcf.gz.tbi


	url="https://storage.googleapis.com/gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz"
	{	bcftools view -h "$url" | grep -vF -e '##contig='
		for i in chrM chr{1..22} chrX chrY; do
			echo ":INFO: working on $i" >&2
			bcftools view --threads $threads -H "$url" $i
		done
	} | bgzip -kc -@ $threads > "$(basename "$url")"
	tabix -f -p vcf "$(basename "$url")"
	ln -sfnr "$(basename "$url")" ExAC_common.vcf.gz
	ln -sfnr "$(basename "$url").tbi" ExAC_common.vcf.gz.tbi


	url="https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz"
	{	bcftools view -h "$url" | grep -vF -e '##contig='
		for i in chrM chr{1..22} chrX chrY; do
			echo ":INFO: working on $i" >&2
			bcftools view --threads $threads -H "$url" $i
		done
	} | bgzip -kc -@ $threads > "$(basename "$url")"
	tabix -f -p vcf "$(basename "$url")"
	ln -sfnr "$(basename "$url")" gnomAD_with_AF.vcf.gz
	ln -sfnr "$(basename "$url").tbi" gnomAD_with_AF.vcf.gz.tbi


	url="https://storage.googleapis.com/gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"
	{	bcftools view -h "$url" | grep -vF -e '##contig='
		for i in chrM chr{1..22} chrX chrY; do
			echo ":INFO: working on $i" >&2
			bcftools view --threads $threads -H "$url" $i
		done
	} | bgzip -kc -@ $threads > "$(basename "$url")"
	tabix -f -p vcf "$(basename "$url")"
	ln -sfnr "$(basename "$url")" PON.vcf.gz
	ln -sfnr "$(basename "$url").tbi" PON.vcf.gz.tbi
}

alias dlgenome::hg38.dbsnp.ncbi="_bashbone_wrapper dlgenome::hg38.dbsnp.ncbi"
function dlgenome::hg38.dbsnp.ncbi() {
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading dbSNP"

	url="https://ftp.ncbi.nlm.nih.gov/snp/archive/"
	url+="$(curl -s "$url" | grep -oE 'b[0-9]+' | sort -Vr | head -1)/VCF/"
	url+="$(curl -s "$url" | grep -oE '[^">]+\.gz' | sort -Vr | head -1)"
	version="$(bcftools view -h "$url" | grep -m 1 -F dbSNP_BUILD_ID | sed -E 's/.*(dbSNP)_BUILD_ID=(.+)/\1 v.\2/')"
	echo ":INFO: v$version from NCBI"
	# sort -V -> .25.vcf.gz -> hg19
	# sort -Vr -> .[>30].vcf.gz -> hg38

	echo "rm '$(basename "$url")'*" >> "$BASHBONE_CLEANUP" # bcftools fetches csi or tbi from server
	{	bcftools view -h "$url" | grep -vF '##contig='
		for i in NC_012920.1 $(bcftools view -h "$url" | grep -F '##contig=' | grep -oE "NC_[^>]+" | grep -vF NC_012920.1 | sort -V); do
			echo ":INFO: working on "$(sed -E 's/^NC_012920\S+/chrM/;t;s/^NC_000023\S+/chrX/;t;s/^NC_000024\S+/chrY/;t;s/^NC_000+([0-9]+)\.[0-9]+/chr\1/' <<< $i) >&2
			bcftools view --threads $threads -H "$url" $i 2> >(grep -v "Extreme INFO" >&2 || true) | grep -F COMMON | grep -F 1000Genomes | sed -E 's/^NC_012920\S+/chrM/;t;s/^NC_000023\S+/chrX/;t;s/^NC_000024\S+/chrY/;t;s/^NC_000+([0-9]+)\.[0-9]+/chr\1/'
		done
	} | bgzip -kc -@ $threads > "$genome.fa.vcf.gz"
	tabix -f -p vcf "$genome.fa.vcf.gz"

	cat <<- EOF >> "$genome.fa.vcf.README"
		$(date)
		$USER
		NCBI v$version $genome
		COMMON only i.e. MAF > 0.01 and listed in 1000Genomes
		$url
	EOF

	return 0
}

alias dlgenome::hg38.dbsnp.ensembl="_bashbone_wrapper dlgenome::hg38.dbsnp.ensembl"
function dlgenome::hg38.dbsnp.ensembl() {
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading dbSNP"
	# url="https://ftp.ensembl.org/pub/current_variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"
	url="https://ftp.ensembl.org/pub/"
	# url+="$(curl -s "$url" | grep -oE 'release-[0-9]+' | sort -Vr | head -1)/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"
	url+="$(curl -s "${url/https/ftp}" | grep -E "current\s*->\s*" | awk '{print $NF}')/variation/vcf/homo_sapiens/1000GENOMES-phase_3.vcf.gz"
	version="$(bcftools view -h "$url" | grep -m 1 -oE dbSNP_[0-9]+ | cut -d _ -f 2)"
	echo ":INFO: v$version from Ensembl"

	echo "rm '$(basename "$url")'*" >> "$BASHBONE_CLEANUP" # bcftools fetches csi or tbi from server
	{	bcftools view -h "$url" 2> >(grep -v "Invalid tag name" >&2 || true) | grep -vF -e '##contig=' -e "##INFO=<ID=HGMD-PUBLIC"
		for i in MT {1..22} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H -i '(MAF!="." && MAX(MAF)>0.01) || (MAF="." && (MAX(AFR)>0.01 || MAX(AMR)>0.01 || MAX(EAS)>0.01 || MAX(EUR)>0.01 || MAX(SAS)>0.01))' "$url" $i 2> >(grep -v "Invalid tag name" >&2 || true) | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$genome.fa.vcf.gz"
	tabix -f -p vcf "$genome.fa.vcf.gz"

	cat <<- EOF >> "$genome.fa.vcf.README"
		$(date)
		$USER
		Ensembl v$version $genome
		common only i.e. MAF > 0.01 or if MAF absent one of populations > 0.01
		$url
	EOF

	return 0
}

alias dlgenome::hg38.dbsnp.ucsc="_bashbone_wrapper dlgenome::hg38.dbsnp.ucsc"
function dlgenome::hg38.dbsnp.ucsc() {
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: UCSC dbSNP not supported"
	return 0

	# TODO needs more effort due to missing start position ref. see also
	# https://groups.google.com/a/soe.ucsc.edu/g/genome/c/DIZVXERQnnk
	# https://genome.ucsc.edu/cgi-bin/hgTrackUi?db=hg19&g=dbSnp155Composite
	# e.g.
	# vcf: chr1 1234 A AT
	# bb:  chr1 1234 1234 "" "T"

	# http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/dbSnp155Common.bb
	url=$(curl -s "http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/" | grep -oE '[^">]+Common.bb' | sort -Vr | head -1)
	url="http://hgdownload.soe.ucsc.edu/gbdb/hg38/snp/$url"
	# wget -c -q --show-progress  --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping -O $genome.fa.bb "$url"
	# needs dl first random access via url not possible
	{	bcftools view -h "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/00-common_all.vcf.gz" | grep -vF '##contig='
		for i in chrM chr{1..22} chrX chrY; do
			echo "working on $i" >&2
			bigBedToBed -chrom=$i "$genome.fa.bb" stdout | perl -lane '
				print
				@maf=split(/,/,$F[9]);
				@alt_observed=split(/,/,$F[6]);
				if ($#alt_observed == -1){
					$x="deletion"
					@alt_project=split(/,/,$F[10]);
				} else {
					@alt_project=split(/,/,$F[11]);
				}
			' | head
		done
	} #| bgzip -kc -@ $threads > "$genome.fa.vcf.gz"
	return
	tabix -f -p vcf "$genome.fa.vcf.gz"

	# for project order and description do
	# curl -s "https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=varRep&hgta_track=dbSnp155Composite&hgta_table=dbSnp155&hgta_doSchema=describe+table+schema" | sed -n '/<ol>/,/<\/ol>/{p}' | grep 'href' | sed -E 's/^\s*//;:a;s/<[^>]+>//;ta;' | tr -s ' ' '_'
	# remember: VCF position is bed start+1 coordiante
	# chr1	11008	rs575272151	C	G	.	.	RS=575272151;RSPOS=11008;dbSNPBuildID=142;SSR=0;SAO=0;VP=0x050000020005150024000100;GENEINFO=DDX11L1:100287102;WGT=1;VC=SNV;R5;ASP;VLD;G5;KGPhase3;CAF=0.9119,0.08806;COMMON=1
	# chr1	11007	11008	rs575272151	C	2	G,T,	0	31	0.0880591,8.43028e-05,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,-inf,	C,C,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,	G,T,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,	1986	snv	commonSome,rareSome,overlapDiffClass,	151778775275	194

	# according to ucsc, MAF threshold is applied on any project reporting frequenies not just 1000genomes i.e. even some CAF[0]>0.99 are included
	# ucsc bigbed has always 1000genomes MAF (first entry in list) and
	# 1) sometimes more ALT than MAF listed
	# 2) different MAF listed for same ALT and sometimes even REF
	# 3) sometimes REF is minor and has to be excluded from ALT and CAF[0]=1000genomesMAF
	# -> solve by extract from list, sort and report
}

alias dlgenome::mm10.dbsnp.gatk="_bashbone_wrapper dlgenome::mm10.dbsnp.gatk"
function dlgenome::mm10.dbsnp.gatk() {
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: GATK bundle not available"
	return 0
}

alias dlgenome::mm10.dbsnp.ncbi="_bashbone_wrapper dlgenome::mm10.dbsnp.ncbi"
function dlgenome::mm10.dbsnp.ncbi() {
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: downloading dbSNP"
	url="https://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz"
	version="$(bcftools view -h "$url" | grep -m 1 -F dbSNP_BUILD_ID | sed -E 's/.*(dbSNP)_BUILD_ID=(.+)/\1 v.\2/')"
	echo ":INFO: v$version from NCBI"

	# 81432271 -> 5570707 VLD validation tag
	# likewise ucsc does for mm11, ncbi refers to EVA current_ids which needs to be filtered by validation tag RS_VALIDATED
	# https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/mus_musculus/

	echo "rm '$(basename "$url")'*" >> "$BASHBONE_CLEANUP" # bcftools fetches csi or tbi from server
	{	bcftools view -h "$url" | grep -v '##contig='
		for i in MT {1..19} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H -i 'VLD=1' "$url" $i | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$genome.fa.vcf.gz"
	tabix -f -p vcf "$genome.fa.vcf.gz"

	cat <<- EOF >> "$genome.fa.vcf.README"
		$(date)
		$USER
		NCBI v$version $genome
		filtered for VLD
		$url
	EOF

	return 0
}

alias dlgenome::mm10.dbsnp.ensembl="_bashbone_wrapper dlgenome::mm10.dbsnp.ensembl"
function dlgenome::mm10.dbsnp.ensembl() {
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: Ensembl dbSNP not supported"
	return 0

	url="https://ftp.ensembl.org/pub/release-102/variation/vcf/mus_musculus/mus_musculus.vcf.gz"
	# contains 80M variants without any INFO on how to filter them
}

alias dlgenome::mm10.dbsnp.eva="_bashbone_wrapper dlgenome::mm10.dbsnp.eva"
function dlgenome::mm10.dbsnp.eva() {
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: downloading dbSNP"
	url="https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/"
	version="$(curl -s "$url" | grep -oE 'release_[0-9]+' | sort -Vr | head -1)"
	url+="$version/by_species/mus_musculus/"
	url+="$(curl -s "$url" | grep -oE 'GRCm38.p[0-9]+' | sort -Vr | head -1)/"
	url+="$(curl -s "$url" | grep -oE '[^"]+current_ids.vcf.gz' | head -1)"
	echo ":INFO: v$(echo "$version" | cut -d _ -f 2 -) from EVA"

	# likewise ucsc does for mm11, ncbi refers to EVA current_ids which needs to be filtered by validation tag RS_VALIDATED
	# https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/mus_musculus/
	# 82691010 -> 5679799

	echo "rm '$(basename "$url")'*" >> "$BASHBONE_CLEANUP" # bcftools fetches csi or tbi from server
	{	bcftools view -h "$url" | grep -v '##contig='
		for i in MT {1..19} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H -i 'RS_VALIDATED=1' "$url" $i | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$genome.fa.vcf.gz"
	tabix -f -p vcf "$genome.fa.vcf.gz"

	cat <<- EOF >> "$genome.fa.vcf.README"
		$(date)
		$USER
		EVA v$(echo "$version" | cut -d _ -f 2 -) $genome
		filtered for RS_VALIDATED
		$url
	EOF

	return 0
}

alias dlgenome::mm10.dbsnp.ucsc="_bashbone_wrapper dlgenome::mm10.dbsnp.ucsc"
function dlgenome::mm10.dbsnp.ucsc() {
	cd "$outdir"
	genome=GRCm38

	tfile="$(mktemp -p "${TMPDIR:-/tmp}" snp142Common.XXXXXXXXXX.bed)"
	echo ":INFO: downloading dbSNP"
	url="https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/snp142Common.txt.gz"
	# echo "from NCBI v142 for filtering"
	wget -O - -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused --timestamping "$url" | bgzip -kcd -@ 1 > "$tfile"

	# url="https://ftp.ncbi.nih.gov/snp/organisms/archive/mouse_10090/VCF/00-All.vcf.gz"
	# better use https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/mus_musculus/GRCm38.p4/10090_GCA_000001635.6_current_ids.vcf.gz
	# alternative: https://hgdownload.soe.ucsc.edu/gbdb/mm10/bbi/evaSnp6.bb
	# filter by UCSC common
	# 81432271 -> 8404784 VLD using NCBI
	# 82691010 -> 8592575 RS_VALIDATED using EVA

	url="https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/mus_musculus/"
	url+="$(curl -s "$url" | grep -oE 'GRCm38.p[0-9]+' | sort -Vr | head -1)/"
	url+="$(curl -s "$url/" | grep -oE '[^"]+current_ids.vcf.gz' | head -1)"
	version="$(bcftools view -h "$url" | grep -m 1 -F dbSNP_BUILD_ID | sed -E 's/.*(dbSNP)_BUILD_ID=(.+)/\1 v.\2/')"
	echo ":INFO: v$version from UCSC"

	echo "rm '$(basename "$url")'*" >> "$BASHBONE_CLEANUP" # bcftools fetches csi or tbi from server
	{	bcftools view -h "$url" | grep -v '##contig='
		for i in MT {1..19} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H -T <(awk -v i=$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) '$1==i' "$tfile" | sed 's/^chrM/MT/;t;s/^chr//') "$url" $i | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$genome.fa.vcf.gz"
	tabix -f -p vcf "$genome.fa.vcf.gz"

	cat <<- EOF >> "$genome.fa.vcf.README"
		$(date)
		$USER
		UCSC v$version $genome
		filtered for https://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/snp142Common.txt.gz
		$url
	EOF

	return 0
}

alias dlgenome::mm11.dbsnp.gatk="_bashbone_wrapper dlgenome::mm11.dbsnp.gatk"
function dlgenome::mm11.dbsnp.gatk() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: GATK bundle not available"
	return 0
}

alias dlgenome::mm11.dbsnp.ncbi="_bashbone_wrapper dlgenome::mm11.dbsnp.ncbi"
function dlgenome::mm11.dbsnp.ncbi() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: NCBI dbSNP not available"
	return 0
}

alias dlgenome::mm11.dbsnp.ensembl="_bashbone_wrapper dlgenome::mm11.dbsnp.ensembl"
function dlgenome::mm11.dbsnp.ensembl() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: Ensembl dbSNP not supported"
	return 0

	url="https://ftp.ensembl.org/pub/"
	# url+="$(curl -s "$url" | grep -oE 'release-[0-9]+' | sort -Vr | head -1)/variation/vcf/mus_musculus/mus_musculus.vcf.gz"
	url+="$(curl -s "${url/https/ftp}" | grep -E "current\s*->\s*" | awk '{print $NF}')/variation/vcf/mus_musculus/mus_musculus.vcf.gz"
	# contains 80M variants without any INFO on how to filter them
}

alias dlgenome::mm11.dbsnp.eva="_bashbone_wrapper dlgenome::mm11.dbsnp.eva"
function dlgenome::mm11.dbsnp.eva() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: eva dbsnp not supported"
	return 0
	# in contrast to mm10, likewise to ncbi, the vcf file has no validation info included

	echo ":INFO: downloading dbSNP"

	url="https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/"
	version="$(curl -s "$url" | grep -oE 'release_[0-9]+' | sort -Vr | head -1)"
	url+="$version/by_species/mus_musculus/GRCm39/"
	url+="$(curl -s "$url" | grep -oE '[^"]+current_ids.vcf.gz' | head -1)"
	echo ":INFO: v$(echo "$version" | cut -d _ -f 2-) from EVA"

	echo "rm '$(basename "$url")'*" >> "$BASHBONE_CLEANUP" # bcftools fetches csi or tbi from server
	{	bcftools view -h "$url" | grep -v '##contig='
		for i in MT {1..19} X Y; do
			echo ":INFO: working on "$(sed 's/^MT/chrM/;t;s/^/chr/' <<< $i) >&2
			bcftools view --threads $threads -H -i 'RS_VALIDATED=1' "$url" $i | sed 's/^MT/chrM/;t;s/^/chr/'
		done
	} | bgzip -kc -@ $threads > "$genome.fa.vcf.gz"
	tabix -f -p vcf "$genome.fa.vcf.gz"

	cat <<- EOF >> "$genome.fa.vcf.README"
		$(date)
		$USER
		EVA v$(echo "$version" | cut -d _ -f 2-) $genome
		$url
		filtered for RS_VALIDATED
	EOF
}

alias dlgenome::mm11.dbsnp.ucsc="_bashbone_wrapper dlgenome::mm11.dbsnp.ucsc"
function dlgenome::mm11.dbsnp.ucsc() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: UCSC dbSNP not supported"
	return 0
	# bigBedToBed https://hgdownload.soe.ucsc.edu/gbdb/mm39/bbi/evaSnp6.bb
	# ucsc data source is eva
	# https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/release_6/by_species/mus_musculus/GRCm39/10090_GCA_000001635.9_current_ids.vcf.gz
}

################################################################################
# CTAT
################################################################################

alias dlgenome::hg19.ctat="_bashbone_wrapper dlgenome::hg19.ctat"
function dlgenome::hg19.ctat() {
	cd "$outdir"
	genome=GRCh37

	echo ":INFO: downloading CTAT"
	url='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/'
	url+="$(curl -sl "$url" | grep -oE 'GRCh37[^\"]+\.plug-n-play\.tar\.gz' | grep -v STAR | sort -Vr | head -1)"
	version="$(basename "$url" .tar.gz | rev | cut -d _ -f 1 | rev)"
	echo ":INFO: v$version from Broad Institute"

	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused "$url" | bgzip -kcd -@ 1 | tar -x
	mv "$(basename "$url" .tar.gz)" "${genome}_CTAT_genome_lib"

	cd "${genome}_CTAT_genome_lib"
	cat <<- EOF >> "$genome.README"
		$(date)
		$USER
		CTAT v$version $genome
		$url
	EOF

	return 0
}

alias dlgenome::hg38.ctat="_bashbone_wrapper dlgenome::hg38.ctat"
function dlgenome::hg38.ctat() {
	cd "$outdir"
	genome=GRCh38

	echo ":INFO: downloading CTAT"

	url='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/'
	url+="$(curl -sl "$url" | grep -oE 'GRCh38[^\"]+\.plug-n-play\.tar\.gz' | grep -v STAR | sort -Vr | head -1)"
	version="$(basename "$url" .tar.gz | rev | cut -d _ -f 1 | rev)"
	echo ":INFO: v$version from Broad Institute"

	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused "$url" | bgzip -kcd -@ 1 | tar -x
	mv "$(basename "$url" .tar.gz)" "${genome}_CTAT_genome_lib"

	cd "${genome}_CTAT_genome_lib"
	cat <<- EOF >> "$genome.README"
		$(date)
		$USER
		CTAT v$version $genome
		$url
	EOF

	return 0
}

alias dlgenome::mm10.ctat="_bashbone_wrapper dlgenome::mm10.ctat"
function dlgenome::mm10.ctat() {
	cd "$outdir"
	genome=GRCm38

	echo ":INFO: CTAT not available anymore"
	return 0
}

alias dlgenome::mm11.ctat="_bashbone_wrapper dlgenome::mm11.ctat"
function dlgenome::mm11.ctat() {
	cd "$outdir"
	genome=GRCm39

	echo ":INFO: downloading CTAT"
	url='https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/'
	url+="$(curl -s "$url" | grep -oE 'Mouse_GRCm39[^\"]+\.plug-n-play\.tar\.gz' | grep -v STAR | sort -Vr | head -1)"
	version="$(basename "$url" .tar.gz | rev | cut -d _ -f 1 | rev)"
	echo ":INFO: v$version from Broad Institute"

	wget -c -q --show-progress --progress=bar:force --waitretry=10 --tries=10 --retry-connrefused "$url" | bgzip -kcd -@ 1 | tar -x
	mv "$(basename "$url" .tar.gz)" "${genome}_CTAT_genome_lib"

	cd "${genome}_CTAT_genome_lib"
	cat <<- EOF >> "$genome.README"
		$(date)
		$USER
		CTAT v$version $genome
		$url
	EOF

	return 0
}

################################################################################
# MAIN
################################################################################

checkopt() {
	local arg=false
	case $1 in
		-h | --h | -help | --help) { usage || exit 0; };;
		-t | --t | -threads | --threads) arg=true; threads=$2;;
		-o | --o | -out | --out) arg=true; outdir="$2";;
		-r | --r | -reference | --reference) arg=true; ref=$2;;
		-g | --g | -genome | --genome) fun+=("genome");;
		-c | --c | -ctat | --ctat) fun+=("ctat");;
		-a | --a | -annotation | --annotation) fun+=("gtf");;
		-s | --s | -dbsnp | --dbsnp) db="ensembl";;
		-n | --n | -ncbi | --ncbi) db='ncbi';;
		-u | --u | -ucsc | --ucsc) db='ucsc';;
		-k | --k | -gatk | --gatk) db='gatk';;
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
BASHBONE_ERROR="mandatory parameter -r missing"
[[ $ref ]]
if [[ $db ]]; then
	[[ "$db" == "ensembl" && $ref == mm* ]] && db="eva"
	fun+=("dbsnp.$db")
fi
[[ $godb ]] && fun+=("go")
outdir="${outdir:-$PWD}"
BASHBONE_ERROR="cannot create $outdir"
mkdir -p "$outdir"
outdir="$(realpath -se "$outdir")"
log="$outdir/dlgenome.log"
rm -f "$log"
BASHBONE_ERROR="cannot create $log"
touch "$log"

for f in "${fun[@]}"; do
	BASHBONE_ERROR="dlgenome::$ref.$f failed"
	eval dlgenome::$ref.$f 2>&1 | tee -ai "$log"
done
echo ":INFO: success" | tee -ai "$log"

exit 0
