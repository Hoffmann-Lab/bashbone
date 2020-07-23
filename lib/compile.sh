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
			r)	((++mandatory)); _insdir_parse="$OPTARG";;
			s)	((++mandatory)); _threads_parse=$OPTARG;;
			i)	((++mandatory)); _insdir_parse="$OPTARG";;
			t)	((++mandatory)); _threads_parse=$OPTARG;;
			*)	compile::_usage; return 1;;
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
		compile::dexseq -i "$insdir" -t $threads && \
		compile::wgcna -i "$insdir" -t $threads && \
		compile::dgca -i "$insdir" -t $threads && \
		compile::revigo -i "$insdir" -t $threads && \
		compile::knapsack -i "$insdir" -t $threads
	} || return 1

	return 0
}

compile::bashbone() {
	local insdir threads src=$(dirname $(readlink -e $0))
	compile::_parse -r insdir -s threads "$@"

	local version 
	source $src/lib/version.sh

	commander::print "installing bashbone"
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

	commander::print "installing conda and tools"
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
			fastqc trimmomatic rcorrector \
			star star-fusion bwa hisat2 \
			samtools picard bamutil \
		chmod 755 $insdir/conda/envs/py2/bin/run_rcorrector.pl && \
		conda list -n py2 -f "fastqc|trimmomatic|rcorrector|star|star-fusion|bwa|hisat2|samtols|picard" | grep -v '^#' > $insdir/condatools.txt && \

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
# Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('cran/WebGestaltR', threads=$threads)"
# from src as linked at bioconductor.org
# wget -O ~/src.tar.gz http://master.bioconductor.org/packages/release/bioc/src/contrib/EnrichmentBrowser_2.14.0.tar.gz
# Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); install.packages('~/src.tar.gz', repos = NULL, type = 'source')"

compile::java() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::print "installing java"
	{	url="https://download.oracle.com/otn-pub/java/jdk/13.0.2+8/d4173c853231432d94f001e99d882ca7/jdk-13.0.2_linux-x64_bin.tar.gz" && \
		wget -q --no-cookies --no-check-certificate --header "Cookie: oraclelicense=accept-securebackup-cookie" -O $insdir/java.tar.gz $url && \
		version=$(echo $url | perl -lane '$_=~/jdk-([^-_]+)/; print $1') && \
		tar -xzf $insdir/java.tar.gz -C $insdir && \
		rm $insdir/java.tar.gz && \
		mkdir -p $insdir/latest && \
		ln -sfn $(ls -vd $insdir/jdk*/bin/ | tail -1) $insdir/latest/java
	} || return 1

	return 0
}

compile::_javawrapper() {
	cat <<- EOF > "$1" || return 1
		#!/usr/bin/env bash
		java=java
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
	# conda trimmomatic wrapper is written in python and thus dont like process substitutions as files
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::print "installing trimmomatic"
	{	source $insdir/conda/bin/activate py2 && \
		url='http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip' && \
		wget -q $url -O $insdir/trimmomatic.zip && \
		unzip -d $insdir $insdir/trimmomatic.zip && \
		rm $insdir/trimmomatic.zip && \
		cd $(ls -dv $insdir/Trimmomatic-*/ | tail -1) && \
		mkdir -p bin && \
		commander::_javawrapper bin/trimmomatic $(readlink -e trimmomatic-*.jar) && \
		ln -sfn $PWD/bin $insdir/latest/trimmomatic
	}
}

compile::sortmerna() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::print "installing sortmerna"
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

compile::sortmerna_new_buggy() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::print "installing sortmerna"
	{	source $insdir/conda/bin/activate py2 && \
		url='https://github.com/biocore/sortmerna/archive/v3.0.3.tar.gz' && \
		wget -q $url -O $insdir/sortmerna.tar.gz && \
		tar -xzf $insdir/sortmerna.tar.gz -C $insdir && \
		rm $insdir/sortmerna.tar.gz && \
		cd $(ls -dv $insdir/sortmerna-*/ | tail -1) && \
		url='https://github.com/biocore/sortmerna/releases/download/v3.0.3/sortmerna-3.0.3-Linux_U16.sh' && \
		wget -q $url -O install.sh && \
		bash install.sh --prefix=$PWD --skip-license && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/sortmerna
	} || return 1

	for i in rRNA_databases/*.fasta; do
		o=index/$(basename $i .fasta)-L18
		echo -ne "bin/indexdb --ref $i,$o -m $memory -L 18\0"
	done | xargs -0 -P $threads -I {} bash -c {} || return 1

	return 0
}

compile::segemehl() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::print "installing segemehl"
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

	commander::print "installing dexseq"
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

	commander::print "installing wgcna"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('cran/WGCNA', threads=$threads, force=T)"
	} || return 1

	return 0
}

compile::dgca() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::print "installing dgca"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('andymckenzie/DGCA', threads=$threads, force=T)"
	} || return 1

	return 0
}

compile::revigo() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::print "installing revigo"
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

	commander::print "installing knapsack"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); install.packages('knapsack', repos='http://R-Forge.R-project.org')"
	} || return 1

	return 0
}
