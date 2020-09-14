#! /usr/bin/env bash
# (c) Konstantin Riege

compile::_parse(){
	local funcname=${FUNCNAME[0]}
	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
		usage:
			-i <path>    | installation base
			-t <threads> | number of
			-u <upgrade> | true/false conda envs
		EOF
		return 0
	}
	local OPTIND arg mandatory
	declare -n _insdir_parse _threads_parse _upgrade_parse
	while getopts 'r:s:c:i:t:u:' arg; do
		case $arg in
			r) ((++mandatory)); _insdir_parse="$OPTARG";;
			s) ((++mandatory)); _threads_parse=$OPTARG;;
			c) _upgrade_parse=$OPTARG;;
			i) ((++mandatory)); _insdir_parse="$OPTARG";;
			t) ((++mandatory)); _threads_parse=$OPTARG;;
			u) _upgrade_parse=$OPTARG;;
			*) _usage; return 1;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage && return 1

	return 0
}

compile::all(){
	local insdir threads
	(	trap 'exit $?' ERR INT TERM
		set -e
		compile::_parse -r insdir -s threads "$@"
		compile::bashbone -i "$insdir" -t $threads
		compile::conda -i "$insdir" -t $threads
		compile::conda_tools -i "$insdir" -t $threads
		compile::java -i "$insdir" -t $threads
		compile::trimmomatic -i "$insdir" -t $threads
		compile::sortmerna -i "$insdir" -t $threads
		compile::segemehl -i "$insdir" -t $threads
		compile::preparedexseq -i "$insdir" -t $threads
		compile::revigo -i "$insdir" -t $threads
		compile::gem -i "$insdir" -t $threads
		compile::idr -i "$insdir" -t $threads
	)
	return $?
}

compile::bashbone() {
	local insdir threads version src=$(dirname $(readlink -e $0))
	(	trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing bashbone"
		compile::_parse -r insdir -s threads "$@"
		source $src/lib/version.sh
		rm -rf "$insdir/bashbone-$version"
		mkdir -p "$insdir/bashbone-$version"
		cp -r "$src"/* "$insdir/bashbone-$version"
		mkdir -p "$insdir/latest"
		ln -sfn "$insdir/bashbone-$version" "$insdir/latest/bashbone"
	)
	return $?
}

compile::upgrade(){
	local insdir threads
	(	trap 'exit $?' ERR INT TERM
		set -e
		compile::_parse -r insdir -s threads "$@"
		compile::bashbone -i "$insdir" -t $threads
		compile::conda_tools -i "$insdir" -t $threads -u true
	)
	return $?
}

compile::conda() {
	local insdir threads url version tmpdir n bin
	(	trap 'rm -rf "$tmpdir"' EXIT
		trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing conda"
		compile::_parse -r insdir -s threads "$@"
		tmpdir="$insdir/tmp"
		mkdir -p "$tmpdir"

		url="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
		wget -q -O "$insdir/miniconda.sh" "$url"
		mkdir -p "$insdir/conda"
		bash "$insdir/miniconda.sh" -b -u -f -p "$insdir/conda" 2>&1
		rm -f "$insdir/miniconda.sh"

		source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
		conda update -y conda

		# base env: install tools used on the fly and compilers for perl modules
		commander::printinfo "setup conda base env"
		conda install -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 gxx_linux-64 gfortran_linux-64 \
			glib pkg-config make automake cmake \
			bzip2 pigz pbzip2 \
			perl-threaded perl-app-cpanminus perl-list-moreutils perl-try-tiny \
			curl ghostscript dos2unix \
			ucsc-facount khmer \
			datamash samtools bedtools \
			htslib bcftools vcflib vt
		cpanm Switch

		# setup r env with compilers for r packages
		# to avoid r downgrades due to modules r built versions, compile them manually (ggpubr requires nlopt)
		commander::printinfo "setup conda r env for deseq, dexseq, survival, wgcna, ggpubr, pheatmap, dplyr, knapsack"
		n=r
		conda create -y -n $n python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda \
			gcc_linux-64 gxx_linux-64 gfortran_linux-64 \
			glib pkg-config make automake cmake \
			bzip2 pigz pbzip2 \
			htslib nlopt r-base
		for bin in perl samtools bedtools; do
			[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
		# basics
		declare -a cmd1
		commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
			Rscript - <<< '
				options(unzip="$(command -v unzip)");
				Sys.setenv(TAR="$(command -v tar)");
				install.packages(c("BiocManager","devtools","codetools"),
					repos="http://cloud.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
			' 2>&1
		CMD
		# bioconductor
		declare -a cmd2
		commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
			Rscript - <<< '
				options(unzip="$(command -v unzip)");
				Sys.setenv(TAR="$(command -v tar)");
				BiocManager::install(c("BiocParallel","genefilter","DESeq2","DEXSeq","TCGAutils","TCGAbiolinks","impute","preprocessCore","GO.db","AnnotationDbi"),
					ask=F, Ncpus=$threads, clean=T, destdir="$tmpdir");
			' 2>&1
		CMD
		# cran - needs to be last since WGCNA depends on bioconductor packages impute,...
		# R-Forge.r does not complaine about knapsack not being compatible with R>=4
		declare -a cmd3
		commander::makecmd -a cmd3 -s '&&' -c {COMMANDER[0]}<<- CMD
			Rscript - <<< '
				options(unzip="$(command -v unzip)");
				Sys.setenv(TAR="$(command -v tar)");
				install.packages(c("knapsack","WGCNA","dplyr","tidyverse","ggpubr","ggplot2","gplots","RColorBrewer","svglite","pheatmap","data.table"),
					repos="http://cloud.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
				install.packages(c("knapsack"), repos="http://R-Forge.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
			' 2>&1
		CMD
		# github
		declare -a cmd4
		commander::makecmd -a cmd4 -s '&&' -c {COMMANDER[0]}<<- CMD
			Rscript - <<< '
				options(unzip="$(command -v unzip)");
				Sys.setenv(TAR="$(command -v tar)");
				devtools::install_github("andymckenzie/DGCA", upgrade="never", force=T, clean=T, destdir="$tmpdir")"
			' 2>&1
		CMD
		commander::runcmd -c $n -t 1 -a cmd1
		commander::runcmd -c $n -t 1 -a cmd2
		commander::runcmd -c $n -t 1 -a cmd3
		commander::runcmd -c $n -t $threads -a cmd4

		commander::printinfo "conda clean up"
		conda clean -y -a
		conda deactivate
	)
	return $?
}

compile::conda_tools() {
	local insdir threads upgrade=false url version tool n bin
	declare -A envs
	(	trap 'exit $?' ERR INT TERM
		set -e

		compile::_parse -r insdir -s threads -c upgrade "$@"
		source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
		while read -r tool; do
			envs[$tool]=true
		done < <(conda info -e | awk -v prefix="^"$insdir '$NF ~ prefix {print $1}')

		# python 3 envs
		for tool in fastqc cutadapt rcorrector star bwa rseqc subread arriba star-fusion picard bamutil macs2 diego gatk4 freebayes varscan; do
			n=${tool//[^[:alpha:]]/}
			$upgrade && ${envs[$n]:=false} && continue

			commander::printinfo "setup conda $tool env"
			conda create -y -n $n python=3
			conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
			# link commonly used base binaries into env
			for bin in perl samtools bedtools; do
				[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
			done
		done
		chmod 755 "$insdir/conda/envs/rcorrector/bin/run_rcorrector.pl" # necessary fix

		tool=vardict
		n=${tool//[^[:alpha:]]/}
		$upgrade && ${envs[$n]:=false} || {
			commander::printinfo "setup conda $tool env"
			conda create -y -n $n python=3
			conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool vardict-java
			for bin in perl samtools bedtools; do
				[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
			done
		}

		tool=snpeff
		n=${tool//[^[:alpha:]]/}
		$upgrade && ${envs[$n]:=false} || {
			commander::printinfo "setup conda $tool env"
			conda create -y -n $n python=3
			conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool snpsift
			for bin in perl samtools bedtools; do
				[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
			done
		}


		# python 2 envs
		tool=platypus-variant
		n=platypus
		$upgrade && ${envs[$n]:=false} || {
			commander::printinfo "setup conda $tool env"
			conda create -y -n $n python=2
			conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
			for bin in perl samtools bedtools; do
				[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
			done
		}

		# this is a pipeline itself and thus will not be part of bashbone
		# commander::printinfo "setup conda fusion-catcher env"
		# tool=fusion-catcher
		# n=${tool//[^[:alpha:]]/}
		# conda create -y -n $n python=2
		# conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
		# for bin in perl samtools bedtools; do
		# 	[[ $(conda list -n $n -f $bin) ]] && ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		# done
		# conda activate $n
		# commander::printinfo "downloading databases"
		# download-human-db.sh
		# rm -f $FC_DB_PATH/*.tar.gz* # env varibale
		# conda deactivate

		commander::printinfo "conda clean up"
		conda clean -y -a
		conda deactivate
	)
	return $?
}

compile::java() {
	local insdir threads url version
	(	trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing java"
		compile::_parse -r insdir -s threads "$@"
		source $insdir/conda/bin/activate base
		url="https://download.oracle.com/otn-pub/java/jdk/14.0.2+12/205943a0976c4ed48cb16f1043c5c647/jdk-14.0.2_linux-x64_bin.tar.gz"
		wget -q --no-cookies --no-check-certificate --header "Cookie: oraclelicense=accept-securebackup-cookie" -O $insdir/java.tar.gz $url
		version=$(echo $url | perl -lane '$_=~/jdk-([^-_]+)/; print $1')
		tar -xzf $insdir/java.tar.gz -C $insdir
		rm $insdir/java.tar.gz
		mkdir -p $insdir/latest
		ln -sfn $(ls -vd $insdir/jdk-*/bin | tail -1) $insdir/latest/java
	)
	return $?
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
	chmod 755 "$1"

	return 0
}

compile::trimmomatic(){
	# conda trimmomatic wrapper is written in python and thus cannot handle process substitutions
	local insdir threads url
	(	trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing trimmomatic"
		compile::_parse -r insdir -s threads "$@"
		source $insdir/conda/bin/activate base
		url='http://www.usadellab.org/cms/?page=trimmomatic'
		url='http://www.usadellab.org/cms/'$(curl -s $url | grep Version | grep -oE '[^"]+Trimmomatic-[0-9]+\.[0-9]+\.zip' | sort -Vr | head -1)
		wget -q $url -O $insdir/trimmomatic.zip
		unzip -o -d $insdir $insdir/trimmomatic.zip
		rm $insdir/trimmomatic.zip
		cd $(ls -dv $insdir/Trimmomatic-*/ | tail -1)
		mkdir -p $insdir/latest bin
		compile::_javawrapper bin/trimmomatic $(readlink -e trimmomatic-*.jar) $insdir/latest/java/java
		ln -sfn $PWD/bin $insdir/latest/trimmomatic
	)
	return $?
}

compile::sortmerna() {
	local insdir threads url
	(	trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing sortmerna"
		echo $insdir
		compile::_parse -r insdir -s threads "$@"
		source $insdir/conda/bin/activate base
		url='https://github.com/biocore/sortmerna/archive/2.1.tar.gz'
		wget -q $url -O $insdir/sortmerna.tar.gz
		tar -xzf $insdir/sortmerna.tar.gz -C $insdir
		rm $insdir/sortmerna.tar.gz
		cd $(ls -dv $insdir/sortmerna-*/ | tail -1)
		make clean || true
		./configure --prefix=$PWD
		make -j $threads
		make install -i
		mkdir -p $insdir/latest
		ln -sfn $PWD/bin $insdir/latest/sortmerna
		commander::printinfo "indexing databases"
		for i in rRNA_databases/*.fasta; do
			o=index/$(basename $i .fasta)-L18
			echo -ne "bin/indexdb_rna --ref $i,$o -m 4096 -L 18\0"
		done | xargs -0 -P $threads -I {} bash -c {}
	)
	return $?
}

compile::segemehl() {
	local insdir threads url
	(	trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing segemehl"
		compile::_parse -r insdir -s threads "$@"
		source $insdir/conda/bin/activate base
		url='http://www.bioinf.uni-leipzig.de/Software/segemehl/downloads/'
		url="$url"$(curl -s $url | grep -oE 'segemehl-[0-9\.]+\.tar\.gz' | sort -Vr | head -1)
		wget -q $url -O $insdir/segemehl.tar.gz
		tar -xzf $insdir/segemehl.tar.gz -C $insdir
		rm $insdir/segemehl.tar.gz
		cd $(ls -dv $insdir/segemehl-*/ | tail -1)
		export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
		make clean || true
		make -j $threads all
		mkdir -p bin
		mv *.x bin
		touch bin/segemehl bin/haarz
		chmod 755 bin/*
		mkdir -p $insdir/latest
		ln -sfn $PWD/bin $insdir/latest/segemehl
		cat <<- 'EOF' > $insdir/latest/segemehl/segemehl
			#!/usr/bin/env bash
			[[ $CONDA_PREFIX ]] && export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
			l=$(pkg-config --variable=libdir htslib)
			[[ $l ]] && export LD_LIBRARY_PATH=$l
			$(cd $(dirname $0) && echo $PWD)/segemehl.x $*
		EOF
		echo hier
		cat <<- 'EOF' > $insdir/latest/segemehl/haarz
			#!/usr/bin/env bash
			[[ $CONDA_PREFIX ]] && export PKG_CONFIG_PATH=$CONDA_PREFIX/lib/pkgconfig
			l=$(pkg-config --variable=libdir htslib)
			[[ $l ]] && export LD_LIBRARY_PATH=$l
			$(cd $(dirname $0) && echo $PWD)/haarz.x $*
		EOF
	)
	return $?
}

compile::preparedexseq() {
	local insdir threads
	(	trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing dexseq"
		compile::_parse -r insdir -s threads "$@"
		source $insdir/conda/bin/activate base
		cd $insdir
        rm -rf Subread_to_DEXSeq
		git clone https://github.com/vivekbhr/Subread_to_DEXSeq
		cd Subread_to_DEXSeq
		mkdir -p bin
		mv *.py bin
		mkdir -p $insdir/latest
		ln -sfn $PWD/bin $insdir/latest/preparedexseq
	)
	return $?
}

compile::revigo() {
	local insdir threads
	(	trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing revigo"
		compile::_parse -r insdir -s threads "$@"
		source $insdir/conda/bin/activate base
		cd $insdir
		rm -rf revigo
		git clone https://gitlab.leibniz-fli.de/kriege/revigo.git
		cd revigo
		mkdir bin
		compile::_javawrapper bin/revigo $(readlink -e RevigoStandalone.jar) $insdir/latest/java/java
		mkdir -p $insdir/latest
		ln -sfn $PWD/bin $insdir/latest/revigo
	)
	return $?
}

compile::gem() {
	local insdir threads url version
	(	trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing gem"
		compile::_parse -r insdir -s threads "$@"
		source $insdir/conda/bin/activate base
		url='https://groups.csail.mit.edu/cgs/gem/download/'
		url="$url"$(curl -s $url | grep -oE "gem.v[0-9\.]+\.tar\.gz" | sort -Vr | head -1)
		version=$(basename $url | sed -E 's/gem.v([0-9\.]+)\.tar\.gz/\1/')
		wget -q $url -O $insdir/gem.tar.gz
		tar -xzf $insdir/gem.tar.gz -C $insdir
		mv $insdir/gem $insdir/gem-$version
		rm $insdir/gem.tar.gz
		cd $insdir/gem-$version
		mkdir -p bin
		wget -q -O bin/Read_Distribution_default.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_default.txt
		wget -q -O bin/Read_Distribution_CLIP.txt https://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_CLIP.txt
		compile::_javawrapper bin/gem $(readlink -e gem.jar) $insdir/latest/java/java
		mkdir -p $insdir/latest
		ln -sfn $PWD/bin $insdir/latest/gem
	)
	return $?
}

compile::idr() {
	local insdir threads url
	(	trap 'exit $?' ERR INT TERM
		set -e
		commander::printinfo "installing idr"
		compile::_parse -r insdir -s threads "$@"
		source $insdir/conda/bin/activate base
		url=https://github.com/kundajelab/idr
		url="$url/"$(curl -s $url/tags | grep -oE "archive\/[0-9\.]+\.tar\.gz" | sort -Vr | head -1)
		wget -q $url -O $insdir/idr.tar.gz
		tar -xzf $insdir/idr.tar.gz -C $insdir
		rm $insdir/idr.tar.gz
		cd $(ls -vd $insdir/idr-*/ | tail -1)
		pip install numpy matplotlib
		python setup.py install
		mkdir -p $insdir/latest
		ln -sfn $PWD/bin $insdir/latest/idr
	)
	return $?
}

### OLD STUFF

compile::conda_old() {
	local insdir threads url version
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

compile::wgcna() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing wgcna"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); devtools::install_github('cran/WGCNA', upgrade='never', force=T, clean=T)"
	} || return 1

	return 0
}

compile::dgca() {
	local insdir threads
	compile::_parse -r insdir -s threads "$@"

	commander::printinfo "installing dgca"
	{	source $insdir/conda/bin/activate py2r && \
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); devtools::install_github('andymckenzie/DGCA', upgrade='never', force=T, clean=T)"
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

compile::_annovar() {
	local insdir threads url
	compile::_parse -r insdir -s threads "$@" || return 1

	commander::printinfo "installing annovar"
	{	url="http://www.openbioinformatics.org/annovar/download/xxxxxxxx/annovar.latest.tar.gz" && \
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
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); devtools::install_github('compgenomics/MeTPeak', build_opts = c('--no-resave-data', '--no-manual'), upgrade='never', force=T, clean=T)"
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
		Rscript -e "options(unzip='$(which unzip)'); Sys.setenv(TAR='$(which tar)'); devtools::install_github('cran/WebGestaltR', upgrade='never', force=T, clean=T)"
	} || return 1

	return 0
}
