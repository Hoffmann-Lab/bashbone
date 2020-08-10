#! /usr/bin/env bash
# (c) Konstantin Riege

compile::_usage(){
	cat <<- EOF 
		usage:
			-i <path>    | installation base
			-t <threads> | number of
	EOF
	return 0
}

compile::_parse(){
	local OPTIND arg mandatory
	declare -n _insdir_parse _threads_parse
	while getopts 'r:s:i:t:' arg; do
		case $arg in
			r) ((++mandatory)); _insdir_parse="$OPTARG";;
			s) ((++mandatory)); _threads_parse=$OPTARG;;
			i) ((++mandatory)); _insdir_parse="$OPTARG";;
			t) ((++mandatory)); _threads_parse=$OPTARG;;
			*) compile::_usage; return 1;;
		esac
	done
	[[ $mandatory -lt 4 ]] && { compile::_usage; return 1; }

	return 0
}

compile::all(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	{	compile::bashbone -i "$insdir" -t $threads && \
		compile::conda -i "$insdir" -t $threads && \
		compile::java -i "$insdir" -t $threads && \
		compile::trimmomatic -i "$insdir" -t $threads && \
		compile::sortmerna -i "$insdir" -t $threads && \
		compile::segemehl -i "$insdir" -t $threads && \
		compile::knapsack -i "$insdir" -t $threads && \
		compile::dexseq -i "$insdir" -t $threads && \
		compile::wgcna -i "$insdir" -t $threads && \
		compile::dgca -i "$insdir" -t $threads && \
		compile::revigo -i "$insdir" -t $threads && \
		compile::gem -i "$insdir" -t $threads && \
		compile::idr -i "$insdir" -t $threads
	} || return 1

	return 0
}

compile::bashbone() {
	local insdir threads src=$(dirname $(readlink -e $0))
	compile::_parse -r insdir -s threads "$@"

	local version 
	source $src/lib/version.sh

	commander::printinfo "installing bashbone"
	{	rm -rf "$insdir/bashbone-$version" && \
		mkdir -p "$insdir/bashbone-$version" && \
		cp -r "$src"/* "$insdir/bashbone-$version" && \
		mkdir -p "$insdir/latest" && \
		ln -sfn "$insdir/bashbone-$version" "$insdir/latest/bashbone"
	} || return 1
	return 0
}

compile::upgrade(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads || return 1
	return 0
}

compile::conda() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing conda and tools"
	{	url='https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && \
		wget -q -O $insdir/miniconda.sh $url && \
		version=$(bash $insdir/miniconda.sh -h | grep -F Installs | cut -d ' ' -f 3) && \
		rm -rf $insdir/conda && \
		mkdir -p $insdir/conda && \
		bash $insdir/miniconda.sh -b -f -p $insdir/conda && \
		rm $insdir/miniconda.sh && \
		source $insdir/conda/bin/activate && \

		conda env remove -y -n py2 && \
		conda env remove -y -n py3 && \
		conda create -y -n py2 python=2 && \
		conda create -y -n py2r python=2 && \
		conda create -y -n py3 python=3 && \
		
		# tophat2/hisat2 and some R stuff needs python2 whereas cutadapt,idr,rseqc need python3 env
		# star-fusion needs perl-set-intervaltree perl-db-file perl-set-intervaltree perl-uri perl-io-gzip
		#   installation might be fixed manually via perl-app-cpanminus and execution of cpanm Set::IntervalTree URI ...
		conda install -n py2 -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			perl perl-threaded perl-db-file perl-dbi perl-app-cpanminus perl-list-moreutils perl-try-tiny perl-set-intervaltree perl-uri \
			numpy scipy pysam cython matplotlib \
			datamash \
			fastqc rcorrector \
			star star-fusion bwa hisat2 macs2 \
			samtools picard bamutil bedtools \
			ucsc-facount khmer \
			bcftools gatk4 freebayes varscan platypus-variant vardict vardict-java \
			vcflib vt snpeff snpsift
		chmod 755 $insdir/conda/envs/py2/bin/run_rcorrector.pl && \
		conda list -n py2 -f "fastqc|rcorrector|star|star-fusion|bwa|hisat2|macs2|samtools|picard|bamutil|bedtools|ucsc-facount|khmer" | grep -v '^#' > $insdir/condatools.txt && \
		conda list -n py2 -f "bcftools|gatk4|freebayes|varscan|platypus-variant|vardict|vardict-java|vcflib|vt|snpeff|snpsift" | grep -v '^#' >> $insdir/condatools.txt && \

		conda install -n py3 -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			numpy scipy pysam cython matplotlib \
			cutadapt rseqc htseq diego bedtools && \
		conda list -n py3 -f "cutadapt|rseqc|diego|bedtools|htseq" | grep -v '^#' >> $insdir/condatools.txt && \

		conda install -n py2r -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 readline make automake xz zlib bzip2 pigz pbzip2 ncurses htslib ghostscript \
			r-devtools bioconductor-biocinstaller bioconductor-biocparallel \
			bioconductor-genomicfeatures bioconductor-genefilter \
			subread r-wgcna bioconductor-deseq2 bioconductor-dexseq bioconductor-gseabase bioconductor-clusterprofiler \
			r-dplyr r-ggplot2 r-gplots r-rcolorbrewer r-svglite r-pheatmap r-ggpubr r-treemap r-rngtools && \
		conda list -n py2r -f "subread|r-wgcna|bioconductor-deseq2|bioconductor-dexseq" | grep -v '^#' >> $insdir/condatools.txt && \

		conda clean -y -a
	} || return 1

	return 0
}

## HELP for manual R package installation
# from github:
# sometimes INSTALL_opts = '--no-lock' is required, too
# Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('cran/WebGestaltR', threads=$threads)"
# from src as linked at bioconductor.org
# wget -O ~/src.tar.gz http://master.bioconductor.org/packages/release/bioc/src/contrib/EnrichmentBrowser_2.14.0.tar.gz
# Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); install.packages('~/src.tar.gz', repos = NULL, type = 'source')"

compile::java() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing java"
	{	url="https://download.oracle.com/otn-pub/java/jdk/13.0.2+8/d4173c853231432d94f001e99d882ca7/jdk-13.0.2_linux-x64_bin.tar.gz" && \
		wget -q --no-cookies --no-check-certificate --header "Cookie: oraclelicense=accept-securebackup-cookie" -O $insdir/java.tar.gz $url && \
		version=$(echo $url | perl -lane '$_=~/jdk-([^-_]+)/; print $1') && \
		tar -xzf $insdir/java.tar.gz -C $insdir && \
		rm $insdir/java.tar.gz && \
		mkdir -p $insdir/latest && \
		ln -sfn $(ls -vd $insdir/jdk-*/bin | tail -1) $insdir/latest/java
	} || return 1

	return 0
}

compile::_javawrapper() {
	local java=java
	[[ $3 ]] && java="$3"
	cat <<- EOF > "$1" || return 1
		#!/usr/bin/env bash
		java=$java
		[[ \$JAVA_HOME && -e "\$JAVA_HOME/bin/java" ]] && java="\$JAVA_HOME/bin/java"
		declare -a jvm_mem_args jvm_prop_args pass_args
		for arg in \$@; do
			case \$arg in
				-D*) jvm_prop_args+=("\$arg");;
				-XX*) jvm_prop_args+=("\$arg");;
				-Xm*) jvm_mem_args+=("\$arg");;
				*) pass_args+=("\$arg");;
			esac
		done
		[[ ! \$jvm_mem_args ]] && jvm_mem_args+=("-Xms1024m") && jvm_mem_args+=("-Xmx4g")
		exec "\$java" "\${jvm_mem_args[@]}" "\${jvm_prop_args[@]}" -jar "$2" "\${pass_args[@]}"
	EOF
	chmod 755 "$1" || return 1
	
	return 0
}

compile::trimmomatic(){
	# conda trimmomatic wrapper is written in python and thus cannot handle process substitutions
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing trimmomatic"
	{	source $insdir/conda/bin/activate py2 && \
		url='http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip' && \
		wget -q $url -O $insdir/trimmomatic.zip && \
		unzip -o -d $insdir $insdir/trimmomatic.zip && \
		rm $insdir/trimmomatic.zip && \
		cd $(ls -dv $insdir/Trimmomatic-*/ | tail -1) && \
		mkdir -p $insdir/latest bin && \
		compile::_javawrapper bin/trimmomatic $(readlink -e trimmomatic-*.jar) $insdir/latest/java/java && \
		ln -sfn $PWD/bin $insdir/latest/trimmomatic
	}
}

compile::sortmerna() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing sortmerna"
	{	source $insdir/conda/bin/activate py2 && \
		url='https://github.com/biocore/sortmerna/archive/2.1.tar.gz' && \
		wget -q $url -O $insdir/sortmerna.tar.gz && \
		tar -xzf $insdir/sortmerna.tar.gz -C $insdir && \
		rm $insdir/sortmerna.tar.gz && \
		cd $(ls -dv $insdir/sortmerna-*/ | tail -1) && \
		make clean; true && \
		./configure --prefix=$PWD && \
		make -j $threads && \
		make install -i && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/sortmerna
	} || return 1

	for i in rRNA_databases/*.fasta; do
		o=index/$(basename $i .fasta)-L18
		echo -ne "bin/indexdb_rna --ref $i,$o -m 4096 -L 18\0"
	done | xargs -0 -P $threads -I {} bash -c {} || return 1

	return 0
}

compile::segemehl() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing segemehl"
	{	source $insdir/conda/bin/activate py2 && \
		url='http://www.bioinf.uni-leipzig.de/Software/segemehl/downloads/segemehl-0.3.4.tar.gz' && \
		wget -q $url -O $insdir/segemehl.tar.gz && \
		tar -xzf $insdir/segemehl.tar.gz -C $insdir && \
		rm $insdir/segemehl.tar.gz && \
        cd $insdir/segemehl-0.3.4 && \
        export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig && \
		make clean; true && \
		make -j $threads all && \
		mkdir -p bin && \
		mv *.x bin && \
		touch bin/segemehl bin/haarz && \
		chmod 755 bin/* && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/segemehl
	} || return 1

	cat <<- 'EOF' > $insdir/latest/segemehl/segemehl || return 1
		#!/usr/bin/env bash
		[[ $CONDA_PREFIX ]] && export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
		l=$(pkg-config --variable=libdir htslib)
		[[ $l ]] && export LD_LIBRARY_PATH=$l
		$(cd $(dirname $0) && echo $PWD)/segemehl.x $*
	EOF
	cat <<- 'EOF' > $insdir/latest/segemehl/haarz || return 1
		#!/usr/bin/env bash
		[[ $CONDA_PREFIX ]] && export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
		l=$(pkg-config --variable=libdir htslib)
		[[ $l ]] && export LD_LIBRARY_PATH=$l
		$(cd $(dirname $0) && echo $PWD)/haarz.x $*
	EOF

	return 0
}

compile::dexseq() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing dexseq"
	{	source $insdir/conda/bin/activate py2 && \
		#cat <(echo '#!/usr/bin/env python') $insdir/conda/envs/py2r/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py > $insdir/conda/envs/py2/bin/dexseq_prepare_annotation.py && \
		rm -f $insdir/conda/envs/py2r/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py && \
		chmod 755 $insdir/conda/envs/py2/bin/dexseq_prepare_annotation.py && \
		cd $insdir && \
        rm -rf Subread_to_DEXSeq && \
		git clone https://github.com/vivekbhr/Subread_to_DEXSeq && \
		cd Subread_to_DEXSeq && \
		mkdir -p bin && \
		mv *.py bin && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/preparedexseq
	} || return 1

	return 0
}

compile::wgcna() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing wgcna"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('cran/WGCNA', threads=$threads, force=T)"
	} || return 1

	return 0
}

compile::dgca() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing dgca"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('andymckenzie/DGCA', threads=$threads, force=T)"
	} || return 1

	return 0
}

compile::revigo() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing revigo"
	{	source $insdir/conda/bin/activate py2 && \
		cd $insdir && \
        rm -rf revigo && \
		git clone https://gitlab.leibniz-fli.de/kriege/revigo.git && \
		mkdir -p $insdir/latest && \
		ln -sfn $insdir/revigo $insdir/latest/revigo
	} || return 1

	return 0
}

compile::knapsack(){
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing knapsack"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); install.packages('knapsack', repos='http://R-Forge.R-project.org')"
	} || return 1

	return 0
}

compile::gem() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing gem"
	{	source $insdir/conda/bin/activate py2 && \
		url='https://groups.csail.mit.edu/cgs/gem/download/gem.v3.4.tar.gz' && \
		version=$(basename $url | sed -E 's/.+v([0-9]+.+).tar.gz/\1/') && \
		wget -q $url -O $insdir/gem.tar.gz && \
		tar -xzf $insdir/gem.tar.gz -C $insdir && \
		mv $insdir/gem $insdir/gem-$version
		rm $insdir/gem.tar.gz && \
		cd $insdir/gem-$version && \
		mkdir -p bin && \
		wget -q -O bin/Read_Distribution_default.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_default.txt && \
		wget -q -O bin/Read_Distribution_CLIP.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_CLIP.txt && \
		compile::_javawrapper bin/gem $(readlink -e gem.jar) $insdir/latest/java/java && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/gem
	} || return 1

	return 0
}

compile::idr() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing idr"
	{	source $insdir/conda/bin/activate py3 && \
		url='https://github.com/kundajelab/idr/archive/2.0.4.2.tar.gz' && \
		wget -q $url -O $insdir/idr.tar.gz && \
		tar -xzf $insdir/idr.tar.gz -C $insdir && \
		rm $insdir/idr.tar.gz && \
		cd $(ls -vd $insdir/idr*/ | tail -1) && \
		pip install numpy matplotlib && \
		python setup.py install && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/idr
	} || return 1

	return 0
}

### OLD STUFF

compile::annovar() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	commander::printinfo "installing annovar"
	{	url="http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz" && \
		wget -q $url -O $insdir/annovar.tar.gz && \
		tar -xzf $insdir/annovar.tar.gz -C $insdir && \
		rm $insdir/annovar.tar.gz && \
		cd $insdir/annovar && \
		url='http://www.openbioinformatics.org/annovar/download/table_annovar.pl' && \
		wget -q $url -O table_annovar.pl && \
		chmod 755 table_annovar.pl && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/annovar
	} || return 1

	return 0
}

compile::_setup_annovar() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	commander::printinfo "configuring annovar databases"
	{	source $insdir/conda/bin/activate py2 && \
		cd -P $insdir/latest/annovar && \
		#refSeq
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/ && \
		# ./annotate_variation.pl -buildver hg19 -downdb ensGene humandb/ #UCSC - needs seq to get mRNA
		# ./annotate_variation.pl -buildver hg19 -downdb knownGene humandb/ #UCSC - needs seq to get mRNA
		# ./annotate_variation.pl -buildver hg19 -downdb seq humandb/hg19_seq
		# ./retrieve_seq_from_fasta.pl humandb/hg19_ensGene.txt -seqdir humandb/hg19_seq -format ensGene -outfile humandb/hg19_ensGeneMrna.fa
		# ./retrieve_seq_from_fasta.pl humandb/hg19_knownGene.txt -seqdir humandb/hg19_seq -format knownGene -outfile humandb/hg19_knownGeneMrna.fa
		url='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/ensemblToGeneName.txt.gz' && \
		wget -q $url -O humandb/ensemblToGeneName.txt.gz && \
		gzip -d humandb/ensemblToGeneName.txt.gz  && \
		# ./annotate_variation.pl -buildver hg19 -downdb cytoBand humandb/
		# ./annotate_variation.pl -buildver hg19 -downdb genomicSuperDups humandb/ 
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar avsnp147 humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp33a humandb/ && \
		# ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp142 humandb
		./annotate_variation.pl -buildver hg19 -downdb tfbsConsSites humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb targetScanS humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb wgRna humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb gwasCatalog humandb/ && \
		./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20170130 humandb/
	} || return 1

	return 0
}

compile::_setup_snpeff() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@" || return 1

	commander::printinfo "configuring snpeff databases"
	{	source $insdir/conda/bin/activate py2 && \
		java -jar snpEff.jar download -v GRCh37.75 && \
		#java -jar snpEff.jar download -v hg19 #hg19: UCSC, hg19kg: UCSC knownGenes, GRCh37.75: Ensembl 
		url='http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/promoter_predictions/master_known.bed' && \
		wget -q $url -O data/promoter.bed && \
		url='http://ftp.ebi.ac.uk/pub/databases/ensembl/encode/integration_data_jan2011/byDataType/openchrom/jan2011/mirna_tss/miRNA_promoters_hg19_edited_data.bed' && \
		wget -q $url -O data/miRNApromoter.bed && \
		url='ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz' && \
		wget -q $url -O data/clinvar.vcf.gz && \
		url='ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz.tbi' && \
		wget -q $url -O data/clinvar.vcf.gz.tbi && \
		url='ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz' && \
		wget -q $url -O data/exac.vcf.gz && \
		url='ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/ExAC.r1.sites.vep.vcf.gz.tbi' && \
		wget -q $url -O data/exac.vcf.gz.tbi && \
		# url='https://drive.google.com/open?id=0B60wROKy6OqceTNZRkZnaERWREk'
		#(see https://sites.google.com/site/jpopgen/dbNSFP)
		url='ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFPv2.9.3.zip' && \
		wget -q $url -O data/dbnsfp.zip && \
		unzip -q -o -d data dbnsfp.zip && \
		rm data/dbnsfp.zip && \
		head -n 1 data/dbNSFP*_variant.chr1 > data/dbNSFP.txt && \
		cat dbNSFP*_variant.chr* | grep -v "^#" >> data/dbNSFP.txt && \
		rm dbNSFP*_variant.chr* && \
		$MUVAC/bin/samtools/bgzip -f -@ $threads < data/dbNSFP.txt > data/dbNSFP.txt.gz && \
		$MUVAC/bin/samtools/tabix -f -s 1 -b 2 -e 2 data/dbNSFP.txt.gz && \
		# url='http://www.genome.gov/admin/gwascatalog.txt'
		url='ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations.tsv' && \
		wget -q $url -O data/gwas.txt
	} || return 1

	return 0
}

compile::m6aviewer() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing m6aviewer"
	{	source $insdir/conda/bin/activate py2 && \
		url='http://dna2.leeds.ac.uk/m6a/m6aViewer_1_6_1.jar' && \
		mkdir -p $insdir/m6aViewer/bin && \
		wget -q $url -O $insdir/m6aViewer/m6aViewer_1_6_1.jar && \
		cd $insdir/m6aViewer && \
		compile::_javawrapper $PWD/bin/m6aViewer $(readlink -e m6aViewer_1_6_1.jar) $insdir/latest/java/java && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/m6aViewer
	} || return 1

	return 0
}

compile::metpeak() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing metpeak"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('compgenomics/MeTPeak', build_opts = c('--no-resave-data', '--no-manual'), threads=$threads, force=T)"
	} || return 1

	return 0
}

compile::zerone() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing zerone"
	{	source $insdir/conda/bin/activate py2 && \
		cd $insdir && \
		rm -rf zerone && \
		git clone https://github.com/nanakiksc/zerone.git && \
		cd zerone && \
		make clean; true && \
		make -j $threads && \
		mkdir bin && \
		mv zerone bin && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/zerone
	} || return 1

	return 0
}

compile::dpgpc() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing dp_gp_cluster"
	{	source $insdir/conda/bin/activate py2 && \
		cd $insdir && \
		rm -rf DP_GP_cluster && \
		git clone https://github.com/PrincetonUniversity/DP_GP_cluster.git && \
		cd DP_GP_cluster && \
		sed -i -r '18,19{s/^#\s*//}' bin/DP_GP_cluster.py && \
		pip install GPy pandas numpy scipy matplotlib cython sklearn && \
		python setup.py install && \
		touch bin/DP_GP_cluster && \
		chmod 755 bin/* && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/DP_GP_cluster
	} || return 1

	cat <<- 'EOF' > $insdir/latest/DP_GP_cluster/DP_GP_cluster || return 1
		#!/usr/bin/env bash
		export PYTHONPATH=$CONDA_PREFIX/lib/python2.7/site-packages/:$PYTHONPATH
		$(cd $(dirname \$0) && echo $PWD)/DP_GP_cluster.py $*
	EOF
	return 0
}

compile::webgestalt() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing webgestalt"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('cran/WebGestaltR', threads=$threads, force=T)"
	} || return 1

	return 0
}
