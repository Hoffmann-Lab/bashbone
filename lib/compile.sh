#! /usr/bin/env bash
# (c) Konstantin Riege

compile::_parse(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	_usage() {
		commander::print {COMMANDER[0]}<<- EOF
		usage:
			-i <path>    | installation base
			-t <threads> | number of
			-u <upgrade> | true/false conda envs
		EOF
		return 1
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
			*) _usage;;
		esac
	done
	[[ $mandatory -lt 4 ]] && _usage

	return 0
}

compile::all(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads
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

	return 0
}

compile::bashbone() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads version src=$(dirname $(readlink -e $0))
	commander::printinfo "installing bashbone"
	compile::_parse -r insdir -s threads "$@"
	source $src/lib/version.sh
	rm -rf "$insdir/bashbone-$version"
	mkdir -p "$insdir/bashbone-$version"
	cp -r "$src"/* "$insdir/bashbone-$version"
	mkdir -p "$insdir/latest"
	ln -sfn "$insdir/bashbone-$version" "$insdir/latest/bashbone"

	return 0
}

compile::upgrade(){
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads
	compile::_parse -r insdir -s threads "$@"
	compile::bashbone -i "$insdir" -t $threads
	compile::conda_tools -i "$insdir" -t $threads -u true

	return 0
}

compile::conda() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'rm -rf "$tmpdir"; trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads url version tmpdir n bin

	commander::printinfo "installing conda"
	compile::_parse -r insdir -s threads "$@"
	tmpdir="$insdir/tmp"
	mkdir -p "$tmpdir"

	url="https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh"
	wget -q -O "$insdir/miniconda.sh" "$url"
	mkdir -p "$insdir/conda"
	bash "$insdir/miniconda.sh" -b -u -f -p "$insdir/conda"
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
		conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	done
	# basics
	declare -a cmd1
	commander::makecmd -a cmd1 -s '&&' -c {COMMANDER[0]}<<- CMD
		Rscript - <<< '
			options(unzip="$(command -v unzip)");
			Sys.setenv(TAR="$(command -v tar)");
			install.packages(c("BiocManager","devtools","codetools"),
				repos="http://cloud.r-project.org", Ncpus=$threads, clean=T, destdir="$tmpdir");
		'
	CMD
	# bioconductor
	declare -a cmd2
	commander::makecmd -a cmd2 -s '&&' -c {COMMANDER[0]}<<- CMD
		Rscript - <<< '
			options(unzip="$(command -v unzip)");
			Sys.setenv(TAR="$(command -v tar)");
			BiocManager::install(c("BiocParallel","genefilter","DESeq2","DEXSeq","TCGAutils","TCGAbiolinks","impute","preprocessCore","GO.db","AnnotationDbi"),
				ask=F, Ncpus=$threads, clean=T, destdir="$tmpdir");
		'
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
		'
	CMD
	# github
	declare -a cmd4
	commander::makecmd -a cmd4 -s '&&' -c {COMMANDER[0]}<<- CMD
		Rscript - <<< '
			options(unzip="$(command -v unzip)");
			Sys.setenv(TAR="$(command -v tar)");
			devtools::install_github("andymckenzie/DGCA", upgrade="never", force=T, clean=T, destdir="$tmpdir");
		'
	CMD
	commander::runcmd -c $n -t 1 -a cmd1
	commander::runcmd -c $n -t 1 -a cmd2
	commander::runcmd -c $n -t 1 -a cmd3
	commander::runcmd -c $n -t $threads -a cmd4

	commander::printinfo "conda clean up"
	conda clean -y -a
	conda deactivate

	return 0
}

compile::conda_tools() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads upgrade=false url version tool n bin doclean=false
	declare -A envs

	compile::_parse -r insdir -s threads -c upgrade "$@"
	source "$insdir/conda/bin/activate" base # base necessary, otherwise fails due to $@ which contains -i and -t
	while read -r tool; do
		envs[$tool]=true
	done < <(conda info -e | awk -v prefix="^"$insdir '$NF ~ prefix {print $1}')

	# python 3 envs
	for tool in fastqc cutadapt rcorrector star bwa rseqc subread arriba star-fusion picard bamutil macs2 diego gatk4 freebayes varscan; do
		n=${tool//[^[:alpha:]]/}
		$upgrade && ${envs[$n]:=false} && continue
		doclean=true

		commander::printinfo "setup conda $tool env"
		conda create -y -n $n python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
		# link commonly used base binaries into env
		for bin in perl samtools bedtools; do
			conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	done
	chmod 755 "$insdir/conda/envs/rcorrector/bin/run_rcorrector.pl" # necessary fix

	tool=vardict
	n=${tool//[^[:alpha:]]/}
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $tool env"
		conda create -y -n $n python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool vardict-java
		for bin in perl samtools bedtools; do
			conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	tool=snpeff
	n=${tool//[^[:alpha:]]/}
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $tool env"
		conda create -y -n $n python=3
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool snpsift
		for bin in perl samtools bedtools; do
			conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}


	# python 2 envs
	tool=platypus-variant
	n=platypus
	$upgrade && ${envs[$n]:=false} || {
		doclean=true

		commander::printinfo "setup conda $tool env"
		conda create -y -n $n python=2
		conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
		for bin in perl samtools bedtools; do
			conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
		done
	}

	# this is a pipeline itself and thus will not be part of bashbone
	# commander::printinfo "setup conda fusion-catcher env"
	# tool=fusioncatcher
	# n=${tool//[^[:alpha:]]/}
	# conda create -y -n $n python=2
	# conda install -n $n -y --override-channels -c iuc -c conda-forge -c bioconda -c main -c defaults -c r -c anaconda $tool
	# for bin in perl samtools bedtools; do
	# 	conda list -n $n -f $bin | grep -qv '^#' || ln -sfnr "$insdir/conda/bin/$bin" "$insdir/conda/envs/$n/bin/$bin"
	# done
	# conda activate $n
	# commander::printinfo "downloading databases"
	# download-human-db.sh
	# rm -f $FC_DB_PATH/*.tar.gz* # env variable
	# conda deactivate

	$doclean && {
		commander::printinfo "conda clean up"
		conda clean -y -a
	}

	conda deactivate
	return 0
}

compile::java() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads url version

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

	return 0
}

compile::_javawrapper() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

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
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	# conda trimmomatic wrapper is written in python and thus cannot handle process substitutions
	local insdir threads url

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

	return 0
}

compile::sortmerna() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads url

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
	make install -i # ignore errors caused by --prefix=$PWD
	mkdir -p $insdir/latest
	ln -sfn $PWD/bin $insdir/latest/sortmerna
	commander::printinfo "indexing databases"
	for i in rRNA_databases/*.fasta; do
		o=index/$(basename $i .fasta)-L18
		echo -ne "bin/indexdb_rna --ref $i,$o -m 4096 -L 18\0"
	done | xargs -0 -P $threads -I {} bash -c {}

	return 0
}

compile::segemehl() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads url

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

	return 0
}

compile::preparedexseq() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads

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

	return 0
}

compile::revigo() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads

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

	return 0
}

compile::gem() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads url version

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

	return 0
}

compile::idr() {
	set -o pipefail
	local error funcname=${FUNCNAME[0]}
	trap 'trap - ERR; trap - RETURN' RETURN
	trap 'configure::err -x $? -f "$funcname" -l $LINENO -e "$error" -c "$BASH_COMMAND"; return $?' ERR

	local insdir threads url

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

	return 0
}
