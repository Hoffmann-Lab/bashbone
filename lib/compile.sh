#! /usr/bin/env bash
# (c) Konstantin Riege

compile::_usage(){
	cat <<- EOF 
		usage:
			-i <path> | installation base
	EOF
	return 0
}

compile::_parse(){
	local OPTIND arg mandatory 
	declare -n _insdir_parse
	while getopts 'v:r:i:' arg; do
		case $arg in
			r)	((++mandatory)); _insdir_parse=$OPTARG;;
			i)	((++mandatory)); _insdir_parse=$OPTARG;;
			*)	compile::_usage; return 1;;
		esac
	done
	[[ $mandatory -lt 2 ]] && { compile::_usage; return 1; }
}

compile::all(){
	local insdir
	compile::_parse -r insdir "$@"

	{	compile::bashbone -i "$insdir" && \
		compile::conda -i "$insdir" && \
		compile::java -i "$insdir" && \
		compile::perlmodules -i "$insdir" && \
		compile::sortmerna -i "$insdir" && \
		compile::segemehl -i "$insdir" && \
		compile::dexseq -i "$insdir" && \
		compile::wgcna -i "$insdir" && \
		compile::dgca -i "$insdir" && \
		compile::revigo -i "$insdir" && \
		compile::gem -i "$insdir" && \
		compile::idr -i "$insdir" && \
		compile::knapsack -i "$insdir"
	} || return 1

	return 0
}

compile::bashbone() {
	local insdir
	compile::_parse -r insdir "$@"

	local version 
	source "$(readlink -e $(dirname $0))"/lib/version.sh

	commander::print "installing bashbone"
	{	mkdir -p "$insdir/bashbone-$version" && \
		cp -r "$(readlink -e $(dirname $0))"/* "$insdir/bashbone-$version" && \
		mkdir -p "$insdir/latest" && \
		ln -sfn "$insdir/bashbone-$version" "$insdir/latest/bashbone"
	} || return 1
	return 0
}

compile::upgrade(){
	local insdir
	compile::_parse -r insdir "$@"

	compile::bashbone -i "$insdir" || return 1
	return 0
}

compile::conda() {
	local insdir
	compile::_parse -r insdir "$@"
	shift $# # necessary for conda activate

	commander::print "installing conda and tools"
	{	mkdir -p $insdir/conda && \
		url='https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh' && \
		wget -q -O $insdir/miniconda.sh $url && \
		bash $insdir/miniconda.sh -b -f -p $insdir/conda && \
		rm $insdir/miniconda.sh && \
		source $insdir/conda/bin/activate && \
		conda create -y -n py2 python=2 && \
		conda create -y -n py2r python=2 && \
		conda create -y -n py3 python=3 && \
		
		# readline 7 causes library version number to be lower on the shared object warnings
		# use quantstack gcc for r and perl module installation - defaults gcc has a weird usage
		# under perl > 5.22 List::MoreUtils installation fails
		# macs2, tophat2/hisat2 and R stuff needs python2 whereas cutadapt,idr,rseqc need python3 env
		conda activate py2 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			gcc-7 libgcc-7 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			make automake zlib ncurses xz bzip2 pigz pbzip2 ghostscript htslib readline=6 perl=5.22 perl-threaded=5.22 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			fastqc trimmomatic rcorrector \
			star bwa hisat2 macs2 \
			samtools picard bedtools && \
		chmod 755 $insdir/conda/envs/py2/bin/run_rcorrector.pl && \
		conda clean -y -a && \

		conda activate py3 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			gcc-7 libgcc-7 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			make automake zlib ncurses xz bzip2 pigz pbzip2 ghostscript readline=6 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			cutadapt rseqc && \
		conda clean -y -a && \

		conda activate py2r && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			gcc-7 libgcc-7 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			make automake zlib ncurses xz bzip2 pigz pbzip2 ghostscript readline=6 && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			r-devtools bioconductor-biocinstaller bioconductor-biocparallel \
			bioconductor-genomicfeatures bioconductor-genefilter \
			subread r-wgcna bioconductor-deseq2 bioconductor-dexseq bioconductor-clusterprofiler && \
		conda install -y --override-channels -c iuc -c bioconda -c main -c conda-forge -c defaults -c quantstack \
			r-dplyr r-ggplot2 r-gplots r-rcolorbrewer r-svglite r-pheatmap r-ggpubr r-treemap r-rngtools && \
		conda clean -y -a && \

		conda deactivate
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
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing java"
	{	url="http://download.oracle.com/otn-pub/java/jdk/11.0.2+9/f51449fcd52f4d52b93a989c5c56ed3c/jdk-11.0.2_linux-x64_bin.tar.gz" && \
		wget -q --no-cookies --no-check-certificate --header "Cookie: oraclelicense=accept-securebackup-cookie" -O $insdir/java.tar.gz $url && \
		tar -xzf $insdir/java.tar.gz -C $insdir && \
		rm $insdir/java.tar.gz && \
		mkdir -p $insdir/latest && \
		ln -sfn $(ls -vd $insdir/jdk*/bin/ | tail -1) $insdir/latest/java
	} || return 1

	return 0
}

compile::_javawrapper() {
	cat <<- EOF > $1 || return 1
		#!/usr/bin/env bash
		set -eu -o pipefail
		export LC_ALL=en_US.UTF-8

		java=java
		if [[ -n \$JAVA_HOME ]]; then
			if [[ -e \$JAVA_HOME/bin/java ]]; then
				java=\$JAVA_HOME/bin/java
			fi
		fi

		jvm_mem_opts=""
		jvm_prop_opts=""
		pass_args=""
		for arg in \$@; do
			case \$arg in
				'-D'*) jvm_prop_opts="\$jvm_prop_opts \$arg";;
				'-XX'*) jvm_prop_opts="\$jvm_prop_opts \$arg";;
				'-Xm'*) jvm_mem_opts="\$jvm_mem_opts \$arg";;
				*) pass_args="\$pass_args \$arg";;
			esac
		done
		[[ ! \$jvm_mem_opts ]] && jvm_mem_opts="-Xms512m -Xmx1g"

		pass_arr=(\$pass_args)
		if [[ \${pass_arr[0]:=} == org* ]]; then
			eval \$java \$jvm_mem_opts \$jvm_prop_opts -cp $2 \$pass_args
		else
			eval \$java \$jvm_mem_opts \$jvm_prop_opts -jar $2 \$pass_args
		fi
		exit
	EOF
	chmod 755 $1 || return 1
	
	return 0
}

compile::perlmodules() {
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing perl modules"
	{	source $insdir/conda/bin/activate py2 && \
		url='cpanmin.us' && \
		mkdir -p $insdir/cpanm && \
		wget -q $url -O $insdir/cpanm/cpanm && \
		chmod 755 $insdir/cpanm/cpanm && \
		$insdir/cpanm/cpanm --reinstall List::MoreUtils Exporter::Tiny Try::Tiny
	} || return 1

	return 0
}

compile::sortmerna() {
	local insdir
	compile::_parse -r insdir "$@"

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
		echo -ne "bin/indexdb_rna --ref $i,$o -m $memory -L 18\0"
	done | xargs -0 -P $threads -I {} bash -c {} || return 1

	return 0
}

compile::sortmerna_new_buggy() {
	local insdir
	compile::_parse -r insdir "$@"

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

compile::segemehl () {
	local insdir
	compile::_parse -r insdir "$@"

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
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing dexseq"
	{	source $insdir/conda/bin/activate py2 && \
		cat <(echo '#!/usr/bin/env python') $insdir/conda/envs/py2r/lib/R/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py > $insdir/conda/envs/py2/bin/dexseq_prepare_annotation.py && \
		chmod 755 $insdir/conda/envs/py2/bin/dexseq_prepare_annotation.py && \
		cd $insdir && \
        rm -rf Subread_to_DEXSeq && \
		git clone https://github.com/vivekbhr/Subread_to_DEXSeq && \
		cd Subread_to_DEXSeq && \
		mkdir -p bin && \
		mv *.py bin && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/preparedexseq && \
		pip install htseq
	} || return 1

	return 0
}

compile::wgcna() {
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing wgcna"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('cran/WGCNA', threads=$threads, force=T)"
	} || return 1

	return 0
}

compile::dgca() {
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing dgca"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('andymckenzie/DGCA', threads=$threads, force=T)"
	} || return 1

	return 0
}

compile::revigo() {
	local insdir
	compile::_parse -r insdir "$@"

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

compile::gem() {
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing gem"
	{	source $insdir/conda/bin/activate py2 && \
		url='https://groups.csail.mit.edu/cgs/gem/download/gem.v3.2.tar.gz' && \
		wget -q $url -O $insdir/gem.tar.gz && \
		tar -xzf $insdir/gem.tar.gz -C $insdir && \
		rm $insdir/gem.tar.gz && \
		cd $insdir/gem && \
		mkdir -p bin && \
		cp Read_Distribution_default.txt bin && \
		compile::_javawrapper $PWD/bin/gem $PWD/gem.jar && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/gem
	} || return 1

	return 0
}

compile::idr() {
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing idr"
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

compile::knapsack(){
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing knapsack"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); install.packages('knapsack', repos='http://R-Forge.R-project.org')"
	} || return 1

	return 0
}

### OLD STUFF

compile::m6aviewer() {
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing m6aviewer"
	{	source $insdir/conda/bin/activate py2 && \
		url='http://dna2.leeds.ac.uk/m6a/m6aViewer_1_6_1.jar' && \
		mkdir -p $insdir/m6aViewer/bin && \
		wget -q $url -O $insdir/m6aViewer/m6aViewer_1_6_1.jar && \
		cd $insdir/m6aViewer && \
		compile::_javawrapper $PWD/bin/m6aViewer $PWD/m6aViewer_1_6_1.jar && \
		mkdir -p $insdir/latest && \
		ln -sfn $PWD/bin $insdir/latest/m6aViewer
	} || return 1

	return 0
}

compile::metpeak() {
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing metpeak"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('compgenomics/MeTPeak', build_opts = c('--no-resave-data', '--no-manual'), threads=$threads, force=T)"
	} || return 1

	return 0
}

compile::zerone() {
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing zerone"
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
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing dp_gp_cluster"
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
	local insdir
	compile::_parse -r insdir "$@"

	commander::print "installing webgestalt"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); library('devtools'); install_github('cran/WebGestaltR', threads=$threads, force=T)"
	} || return 1

	return 0
}
